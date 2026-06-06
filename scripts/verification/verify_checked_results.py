#!/usr/bin/env python3
"""Verify repository-local reproduction artifacts.

This lightweight verifier is intentionally limited to data that is checked in to
this repository. The full Julia data pipeline also needs unregistered Julia
packages plus raw EDF/MAT input data that are not part of the git history.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import h5py

ROOT = Path(__file__).resolve().parents[2]

EXPECTED_CSV = {
    "path": "pat9_snippet_1s.csv",
    "rows": 18,
    "cols": 256,
    "min": 0.0,
    "max": 1.0,
    "sum": 2863.968557173066,
    "sha256": "b4cfdf285cc7524d68a0359f66efabe096a44faee82dd78bf75b7ff4231d39da",
}

EXPECTED_JLD2 = {
    "test_AN_contributions_2022_02_28-131849.jld2": {
        "shape": (10, 14),
        "dtype": "float64",
        "sum": 156.55570543487505,
        "mean": 1.1182550388205361,
        "array_sha256": "aad712280108cceea0a29f6d82a71e9df2f0f1b4113e7c067daef7a67f91db6f",
    },
    "test_contributions_2022_02_28-125224.jld2": {
        "shape": (10, 14),
        "dtype": "float64",
        "sum": 1139097.568704424,
        "mean": 8136.4112050316,
        "array_sha256": "50cd4f0c3db46d4c760f0ed5b8b97c8205c23bcf62ff91371f792ae39a9400d3",
    },
}

README_SCRIPT_REFERENCES = [
    "scripts/contributions_timeseries/contributions_patPAT.jl",
    "scripts/contributions_timeseries/contributions_patPAT_artifacts.jl",
    "scripts/reanalysis/diffs_tricorr.jl",
    "scripts/reanalysis/diffs_tricorr_artifacts.jl",
    "scripts/reanalysis/detecttricorr_seizures.jl",
    "scripts/reanalysis/diffs_aeeg.jl",
    "scripts/reanalysis/diffs_aEEG_artifacts.jl",
    "scripts/reanalysis/detectaeeg_seizures.jl",
]


def assert_close(name: str, actual: float, expected: float, tolerance: float = 1e-9) -> None:
    if abs(actual - expected) > tolerance:
        raise AssertionError(f"{name}: expected {expected}, got {actual}")


def verify_readme_scripts() -> None:
    missing = [script for script in README_SCRIPT_REFERENCES if not (ROOT / script).is_file()]
    if missing:
        raise AssertionError("Missing README-referenced scripts: " + ", ".join(missing))
    print(f"ok: {len(README_SCRIPT_REFERENCES)} README-referenced scripts exist")


def verify_csv() -> None:
    path = ROOT / EXPECTED_CSV["path"]
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != EXPECTED_CSV["sha256"]:
        raise AssertionError(f"{path.name}: expected sha256 {EXPECTED_CSV['sha256']}, got {digest}")

    with path.open(newline="") as handle:
        rows = [[float(value) for value in row] for row in csv.reader(handle)]

    if len(rows) != EXPECTED_CSV["rows"]:
        raise AssertionError(f"{path.name}: expected {EXPECTED_CSV['rows']} rows, got {len(rows)}")
    if any(len(row) != EXPECTED_CSV["cols"] for row in rows):
        raise AssertionError(f"{path.name}: expected every row to have {EXPECTED_CSV['cols']} columns")

    values = [value for row in rows for value in row]
    assert_close(f"{path.name} sum", sum(values), EXPECTED_CSV["sum"])
    assert_close(f"{path.name} min", min(values), EXPECTED_CSV["min"])
    assert_close(f"{path.name} max", max(values), EXPECTED_CSV["max"])
    print(f"ok: {path.name} matches checked-in snippet fingerprint")


def contribution_array(path: Path):
    with h5py.File(path, "r") as handle:
        contribution_ref = handle["contributions"][()]["data"]
        return handle[contribution_ref][()]


def verify_jld2_results() -> None:
    for filename, expected in EXPECTED_JLD2.items():
        path = ROOT / filename
        array = contribution_array(path)
        digest = hashlib.sha256(array.tobytes(order="C")).hexdigest()
        if tuple(array.shape) != expected["shape"]:
            raise AssertionError(f"{filename}: expected shape {expected['shape']}, got {array.shape}")
        if str(array.dtype) != expected["dtype"]:
            raise AssertionError(f"{filename}: expected dtype {expected['dtype']}, got {array.dtype}")
        if digest != expected["array_sha256"]:
            raise AssertionError(f"{filename}: expected array sha256 {expected['array_sha256']}, got {digest}")
        assert_close(f"{filename} sum", float(array.sum()), expected["sum"])
        assert_close(f"{filename} mean", float(array.mean()), expected["mean"])
        print(f"ok: {filename} contributions matrix matches checked-in result fingerprint")


def main() -> None:
    verify_readme_scripts()
    verify_csv()
    verify_jld2_results()
    print("all checked-in verification artifacts reproduced")


if __name__ == "__main__":
    main()
