#!/usr/bin/env python3
"""Run the frozen finite checks for the three defect-insertion arguments."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
REPLAYS = (
    "linear8_replay.py",
    "linear9_replay.py",
    "linear10_replay.py",
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--all-mutants",
        action="store_true",
        help="run each replay's base checks and reject every guarded mutant",
    )
    args = parser.parse_args()
    child_args = ["--all-mutants"] if args.all_mutants else []
    failures = []

    for name in REPLAYS:
        print(f"== {name} ==", flush=True)
        result = subprocess.run(
            [sys.executable, str(HERE / name), *child_args],
            cwd=HERE,
            check=False,
        )
        if result.returncode != 0:
            failures.append((name, result.returncode))

    if failures:
        for name, code in failures:
            print(f"FAILED {name}: exit {code}")
        return 1

    if args.all_mutants:
        print("ALL REPLAY MUTANT SUITES REJECTED")
    else:
        print("ALL BASE REPLAY CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
