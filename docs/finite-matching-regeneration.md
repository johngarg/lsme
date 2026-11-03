# Finite matching regeneration

The 14 affected scalar classes now retain finite contributions from
`epsilonbar * invepsilonbar`. The corrected MatchMakerParser resolves these
products before rendering Python. The generated expressions still contain
uncancelled positive powers and poles; setting both attributes to zero now
extracts the finite Laurent coefficient.

This is an algebraic export correction. It does not rerun the one-loop
matching, change an evanescent scheme, or implement an additional finite
scheme conversion. Do not add a separate numerical regulator-product
correction to these regenerated coefficients.

## Inputs and preserved corrections

The regeneration starts from LSME commit
`67cebeb779a8a7e83114826e5791a16445665ac9` and uses MatchMakerParser commit
`19374b57069443c8607b5e35490a398a3b236e83`. For 13 models the saved inputs are
from `johngarg/linear-one-loop-dict-ref`, commit
`0732710e43a208ccf38d2a44db7704b84e7eaafc`.

Upsilon requires the different saved result from the local MatchMakerParser
`granada/GranadaCapitalUpsilon` archive. Its old export reproduces **all 105
current coefficient methods exactly**, including the coupling-symmetry
correction. The reference repository's Upsilon result does not reproduce
those corrections and was rejected. The selected data and configuration are
preserved under `tools/matching_inputs/Upsilon/`; their hashes are recorded
in `finite-matching-manifest.json`.

For varphi, nine methods differ from the old saved export because of the later
quartic correction: `alphaOH`, `alphaOHBox`, `alphaOHD`, and the unbarred and
barred `alphaOdH`, `alphaOeH`, and `alphaOuH`. None contains a mixed regulator
product. All nine are preserved unchanged. Its 18 regulator-dependent methods
match the saved input and are regenerated normally.

Only return expressions in 61 existing methods are replaced. Class names,
constructors, defaults, metadata, import paths, and every other method remain
unchanged. This retains the Upsilon mass-helper repair and the recent shared
SM-input and WCxf changes. The general/equal-mass/custodial Matchete classes
are outside this change.

## Reproduce into staging

Use Python 3.9 or newer, NumPy, a Wolfram kernel, the corrected
MatchMakerParser checkout, and the reference repository. The baseline LSME
commit must be available in this checkout's Git history. From the LSME root:

```sh
python tools/regenerate_finite_matching.py \
  --parser-root /path/to/MatchMakerParser \
  --reference-root /path/to/linear-one-loop-dict-ref \
  --output-dir /tmp/lsme-finite-regeneration \
  --kernel /Applications/Wolfram.app/Contents/MacOS/WolframKernel
```

The output directory must not already exist. The script writes old and new
exports, baseline classes, reconciled candidate classes, a Wolfram script,
an export log, and a manifest with source hashes. It never overwrites the
package. Mathematica execution must use the approval/escalation required by
the environment in which it runs.

Before accepting a candidate, the script verifies that the old export
reproduces the baseline methods, apart from the explicit varphi overlays.
It then proves an exact Laurent-polynomial identity for each changed method
using rational coefficients and structurally identified flavor/logarithm
atoms. This independent check substitutes `invepsilonbar = 1/epsilonbar`
in both expressions. It also rejects remaining mixed regulator monomials
and verifies that all other methods and signatures are unchanged.

To validate every flavor component at two complex parameter points:

```sh
python tools/validate_finite_matching.py \
  --candidate-dir /tmp/lsme-finite-regeneration/candidate \
  --report /tmp/lsme-finite-validation.json \
  --all-flavours
```

The numerical reference comes from evaluating the **original** classes in
Git, not the new exporter. For `F(e,p)=A(e)+p B(e)`, with degree at most two
in `e`, six evaluations determine `B'(0)` exactly up to floating arithmetic.
The comparison checks the finite part `A(0)+B'(0)`, the remaining pole
residue, tree contributions, and a probe at `e=0.25, p=0.4`. General complex
matrices exercise flavor ordering; these algebraic test points are not
proposed physical benchmarks.

The compact pytest fixture uses three flavor samples per method at the same
two points. To reconstruct it from the old classes, omit `--all-flavours`
and add `--write-fixtures /tmp/finite-fixtures`. This also updates only the
18 affected entries per point in `varphi_zero_coupling_baseline.json`, using
the historical parameters and the original class with the regulator
identity imposed. The original fixture evaluated the two regulator markers
independently, so retaining those numbers would test the superseded algebra.

After reviewing the staged diff, copy the candidate classes and reviewed
fixtures into the package and run `python -m pytest`. Regeneration is kept
separate from applying changes so later hand corrections cannot be silently
overwritten.

## Validation record

The branch baseline passes 81 package tests. The 28 new parameterized
regression cases fail on the original classes and pass with the regenerated
methods. The complete package suite passes **109 tests**. The full flavor
comparison passes **39,528 checks** at two complex parameter points, with a
maximum absolute discrepancy of approximately **2.78e-16**. Results are saved in
`finite-matching-validation.json`; source provenance, changed method names,
and preserved overlays are in `finite-matching-manifest.json`.

The known unrelated flavor-texture, Dirac metadata, older combined-quadruplet,
and Higgs quadratic-parameter questions remain as described in the Fitmaker
audit. This change does not claim to resolve them.
