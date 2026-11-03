#!/usr/bin/env python3
"""Re-export the regulator-dependent LSME methods into a staging directory.

Requires the saved MatchMaker inputs and the corrected MatchMakerParser checkout.
Never writes into lsme/numeric. See docs/finite-matching-regeneration.md.
"""

import argparse
import ast
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys


BASELINE = "67cebeb779a8a7e83114826e5791a16445665ac9"
PARSER_COMMIT = "19374b57069443c8607b5e35490a398a3b236e83"
MODELS = {
    "GranadaCapitalOmega1": 15,
    "GranadaCapitalOmega2": 16,
    "GranadaCapitalOmega4": 17,
    "GranadaCapitalPhi": 19,
    "GranadaCapitalPi1": 12,
    "GranadaCapitalPi7": 13,
    "GranadaCapitalS1": 2,
    "GranadaCapitalUpsilon": None,
    "GranadaCapitalXi1": 6,
    "Granadaomega1": 9,
    "Granadaomega2": 10,
    "Granadaomega4": 11,
    "Granadavarphi": 4,
    "Granadazeta": 14,
}
VARPHI_OVERLAYS = {
    "alphaOH", "alphaOHBox", "alphaOHD", "alphaOdH", "alphaOdHbar",
    "alphaOeH", "alphaOeHbar", "alphaOuH", "alphaOuHbar",
}
REGULATORS = {"epsilonbar", "invepsilonbar"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def fingerprint(node):
    """Hash an AST without positions, avoiding recursion on long generated sums."""
    digest = hashlib.sha256()
    stack = [node]
    while stack:
        item = stack.pop()
        if isinstance(item, ast.AST):
            fields = list(ast.iter_fields(item))
            digest.update(repr((type(item).__name__, [k for k, _ in fields])).encode())
            stack.extend(value for _, value in reversed(fields))
        elif isinstance(item, list):
            digest.update(repr(("list", len(item))).encode())
            stack.extend(reversed(item))
        else:
            digest.update(repr((type(item).__name__, item)).encode())
    return digest.hexdigest()


def methods(source):
    tree = ast.parse(source)
    classes = [node for node in tree.body if isinstance(node, ast.ClassDef)]
    require(len(classes) == 1, "Expected exactly one matching class")
    return classes[0].name, {
        node.name: node for node in classes[0].body if isinstance(node, ast.FunctionDef)
    }


def attributes(node):
    return {
        item.attr for item in ast.walk(node)
        if isinstance(item, ast.Attribute)
        and isinstance(item.value, ast.Name) and item.value.id == "self"
    }


def terms(node):
    """Yield signed additive terms; leave every coefficient atom untouched."""
    stack = [(node, 1)]
    while stack:
        item, sign = stack.pop()
        if isinstance(item, ast.BinOp) and isinstance(item.op, (ast.Add, ast.Sub)):
            stack.append((item.right, -sign if isinstance(item.op, ast.Sub) else sign))
            stack.append((item.left, sign))
        else:
            yield item, sign


def rational(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return Fraction(node.value)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -rational(node.operand)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        return rational(node.left) / rational(node.right)
    raise ValueError("Not a rational constant")


def laurent_polynomial(node, resolve=True):
    """Exact independent check of expanded sums of products, with rational weights.

    Nonregulator atoms (including flavor indices and logarithms) are hashed
    structurally. No floating evaluation or general symbolic simplifier is used.
    """
    result = defaultdict(Fraction)
    for term, sign in terms(node):
        coefficient = Fraction(sign)
        powers = defaultdict(int)
        stack = [(term, 1)]
        while stack:
            item, exponent = stack.pop()
            if isinstance(item, ast.BinOp) and isinstance(item.op, ast.Mult):
                stack.extend([(item.left, exponent), (item.right, exponent)])
                continue
            if isinstance(item, ast.BinOp) and isinstance(item.op, ast.Pow):
                power = rational(item.right)
                require(power.denominator == 1, "Noninteger power")
                stack.append((item.left, exponent * int(power)))
                continue
            try:
                coefficient *= rational(item) ** exponent
                continue
            except ValueError:
                pass
            if (isinstance(item, ast.Attribute) and isinstance(item.value, ast.Name)
                    and item.value.id == "self" and item.attr in REGULATORS):
                if resolve:
                    powers["regulator"] += exponent * (1 if item.attr == "epsilonbar" else -1)
                else:
                    powers[item.attr] += exponent
            else:
                require(not (attributes(item) & REGULATORS), "Unexpanded regulator factor")
                powers[fingerprint(item)] += exponent
        key = tuple(sorted((atom, power) for atom, power in powers.items() if power))
        result[key] += coefficient
    return {key: value for key, value in result.items() if value}


def expression(method):
    require(len(method.body) == 1 and isinstance(method.body[0], ast.Return),
            f"Expected one return expression in {method.name}")
    return method.body[0].value


def splice_returns(baseline, replacements):
    lines = baseline.splitlines(keepends=True)
    offsets = [0]
    for line in lines:
        offsets.append(offsets[-1] + len(line))
    edits = []
    for old, new, new_source in replacements:
        start = offsets[old.lineno - 1] + old.col_offset
        end = offsets[old.end_lineno - 1] + old.end_col_offset
        new_lines = new_source.splitlines(keepends=True)
        selected = new_lines[new.lineno - 1:new.end_lineno]
        if len(selected) == 1:
            replacement = selected[0][new.col_offset:new.end_col_offset]
        else:
            selected[-1] = selected[-1][:new.end_col_offset]
            selected[0] = selected[0][new.col_offset:]
            replacement = "".join(selected)
        edits.append((start, end, replacement))
    for start, end, replacement in sorted(edits, reverse=True):
        baseline = baseline[:start] + replacement + baseline[end:]
    return baseline


def reconcile(baseline, old_export, new_export, model):
    _, current = methods(baseline)
    _, old = methods(old_export)
    _, new = methods(new_export)
    require(current.keys() <= old.keys() and old.keys() == new.keys(),
            f"Unexpected method coverage for {model}")
    differences = {
        name for name in current if name != "__init__"
        and fingerprint(current[name].body) != fingerprint(old[name].body)
    }
    expected = VARPHI_OVERLAYS if model == "Granadavarphi" else set()
    require(differences == expected, f"Unexpected source differences for {model}: {differences}")
    replacements, changed = [], []
    for name, method in current.items():
        if name == "__init__" or name in differences:
            continue
        require(fingerprint(method.args) == fingerprint(new[name].args),
                f"Changed signature: {model}.{name}")
        if fingerprint(old[name].body) == fingerprint(new[name].body):
            continue
        before, after = expression(method), expression(new[name])
        require(REGULATORS <= attributes(before), f"Unexpected change: {model}.{name}")
        require(laurent_polynomial(before) == laurent_polynomial(after),
                f"Laurent identity failed: {model}.{name}")
        for term, _ in terms(after):
            require(not (REGULATORS <= attributes(term)),
                    f"Unresolved regulator product: {model}.{name}")
        replacements.append((before, after, new_export))
        changed.append(name)
    expected_changed = {
        name for name, method in current.items()
        if REGULATORS <= attributes(method)
    }
    require(set(changed) == expected_changed, f"Missed regulator methods for {model}")
    candidate = splice_returns(baseline, replacements)
    compile(candidate, model + "_matching.py", "exec")
    _, candidate_methods = methods(candidate)
    for name in current.keys() - set(changed):
        require(fingerprint(current[name]) == fingerprint(candidate_methods[name]),
                f"Nonmatching code changed: {model}.{name}")
    return candidate, sorted(changed), sorted(differences)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--parser-root", type=Path, required=True)
    cli.add_argument("--reference-root", type=Path, required=True,
                     help="linear-one-loop-dict-ref checkout (contains results/)")
    cli.add_argument("--output-dir", type=Path, required=True, help="A new staging directory")
    cli.add_argument("--kernel", required=True)
    args = cli.parse_args()
    repo = Path(__file__).resolve().parents[1]
    parser = args.parser_root.resolve()
    reference = args.reference_root.resolve()
    stage = args.output_dir.resolve()
    stage.mkdir(parents=True, exist_ok=False)
    for directory in ("baseline", "old", "new", "candidate"):
        (stage / directory).mkdir()
    parser_head = subprocess.check_output(["git", "-C", str(parser), "rev-parse", "HEAD"], text=True).strip()
    require(parser_head == PARSER_COMMIT, "Use the audited MatchMakerParser commit")
    subprocess.run(["git", "-C", str(parser), "diff", "--exit-code", PARSER_COMMIT,
                    "--", "MatchMakerParser.m", "PythonForm.m", "Warsaw.m"], check=True)
    reference_head = subprocess.check_output(
        ["git", "-C", str(reference), "rev-parse", "HEAD"], text=True).strip()
    manifest = {"baseline_commit": BASELINE, "parser_commit": parser_head,
                "reference_commit": reference_head, "models": []}
    commands = [f"Get[{json.dumps(str(parser / 'MatchMakerParser.m'))}];"]
    for model, number in MODELS.items():
        filename = model + "_matching.py"
        baseline = subprocess.check_output(
            ["git", "-C", str(repo), "show", f"{BASELINE}:lsme/numeric/{filename}"], text=True)
        (stage / "baseline" / filename).write_text(baseline)
        class_name, _ = methods(baseline)
        source = (repo / "tools/matching_inputs/Upsilon" if number is None
                  else reference / "results" / str(number))
        config, data = source / "Config.m", source / "MatchingResult.dat"
        manifest["models"].append({
            "model": model, "file": filename, "class": class_name,
            "source": str(source), "config_sha256": sha256(config.read_bytes()),
            "matching_data_sha256": sha256(data.read_bytes()),
            "baseline_sha256": sha256(baseline.encode()),
        })
        commands += ["$Couplings = <||>; $ExoticParams = {};", f"Get[{json.dumps(str(config))}];"]
        for mode in ("old", "new"):
            output = stage / mode / filename
            action = (f"MatchMakerToPython[{json.dumps(class_name.removesuffix('MatchingResult'))}, "
                      f"{json.dumps(str(data))}, {json.dumps(str(output))}]; "
                      f"Close[{json.dumps(str(data))}]; Close[{json.dumps(str(output))}];")
            if mode == "old":
                action = "Block[{ResolveRegulatorProducts = Identity}, " + action + "];"
            commands.append(f'Print["{mode} {model}"]; Print[AbsoluteTiming[{action}][[1]]];')
    commands.append("Exit[0];")
    script = stage / "export.wls"
    script.write_text("\n".join(commands) + "\n")
    with (stage / "export.log").open("w") as log:
        subprocess.run([args.kernel, "-script", str(script)], stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    for entry in manifest["models"]:
        filename = entry["file"]
        baseline, old, new = [(stage / mode / filename).read_text() for mode in ("baseline", "old", "new")]
        candidate, changed, preserved = reconcile(baseline, old, new, entry["model"])
        (stage / "candidate" / filename).write_text(candidate)
        entry.update(changed_methods=changed, preserved_overlays=preserved,
                     candidate_sha256=sha256(candidate.encode()))
        print(f"{entry['model']}: {len(changed)} corrected methods; {len(preserved)} preserved overlays", flush=True)
    manifest["changed_method_count"] = sum(len(entry["changed_methods"]) for entry in manifest["models"])
    (stage / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Staged {len(manifest['models'])} classes in {stage / 'candidate'}")


if __name__ == "__main__":
    sys.setrecursionlimit(100000)
    main()
