#!/usr/bin/env python3
"""Rebuild lsme/numeric/data/sm_table.npz from the installed Wilson.

Run from the repository root after upgrading Wilson or changing the range:

    python scripts/build_sm_table.py [--mu-max-gev 1e6] [--points 64]
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lsme.numeric.sm_table import MU_MAX_GEV, TABLE_PATH, build_sm_table


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mu-max-gev", type=float, default=MU_MAX_GEV)
    parser.add_argument("--points", type=int, default=64)
    parser.add_argument("--validation", type=int, default=20)
    parser.add_argument("--output", type=Path, default=TABLE_PATH)
    args = parser.parse_args()

    table = build_sm_table(args.mu_max_gev, args.points, args.validation)
    table.save(args.output)
    print(f"Wrote {args.output}")
    print(json.dumps(table.metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
