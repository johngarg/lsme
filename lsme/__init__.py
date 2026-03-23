from . import numeric

try:
    from . import symbolic
except ModuleNotFoundError as exc:
    if exc.name != f"{__name__}.symbolic":
        raise
    symbolic = None
