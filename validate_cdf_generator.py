#!/usr/bin/env python3
"""Vipin Singh Sehrawat
Arithmetic limiting laws over shifted primes: Hausdorff dimension zero, Fourier decay, and absolute continuity
arXiv:2603.11196
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import generate_cdf_overlay as generator
import numpy as np


HERE = Path(__file__).resolve().parent


def check_two_prime_law() -> dict:
    n = 200_000
    samples, _, _, counts = generator.sample_products(
        np.array([3, 5], dtype=np.int64), n, 42, progress=False
    )
    support, observed = np.unique(samples, return_counts=True)
    expected_support = np.array([4 / 15, 1 / 3, 2 / 5, 1 / 2])
    probabilities = np.array([1 / 8, 3 / 8, 1 / 8, 3 / 8])
    assert len(support) == 4
    assert np.allclose(support, expected_support, rtol=0, atol=2e-16)
    z = (observed / n - probabilities) / np.sqrt(probabilities * (1 - probabilities) / n)
    assert np.max(np.abs(z)) < 6
    curve, _ = generator.thinned_ecdf(np.sort(samples))
    for value, cumulative in zip(support, np.cumsum(observed) / n):
        rows = curve[curve[:, 0] == value]
        assert len(rows) >= 2 and np.max(rows[:, 1]) == cumulative
    return {"samples": n, "seed": 42, "support": support.tolist(),
            "expected_probabilities": probabilities.tolist(),
            "observed_counts": observed.tolist(), "count_z_scores": z.tolist(),
            "prime_activation_counts": counts.tolist(), "cdf_at_all_atoms": "passed"}


def check_refill_against_scalar_replay() -> dict:
    n, seed, batch = 1234, 12345, 7
    primes = np.array([3, 5, 7], dtype=np.int64)
    original_batch = generator.GAP_BATCH
    try:
        generator.GAP_BATCH = batch
        actual, _, actual_state, actual_counts = generator.sample_products(
            primes, n, seed, progress=False
        )
    finally:
        generator.GAP_BATCH = original_batch

    rng = np.random.Generator(np.random.PCG64(seed))
    reference = np.full(n, 0.5)
    counts = []
    for prime in primes:
        p = int(prime)
        q, position, active = 1 / (p - 1), -1, []
        finished = False
        while not finished:
            size = min(batch, max(32, math.ceil(1.01 * (n - 1 - position) * q) + 16))


            for gap in rng.geometric(q, size=size):
                position += int(gap)
                if position >= n:
                    finished = True
                    break
                active.append(position)
        reference[active] *= 1 - 1 / p
        counts.append(len(active))
    assert np.array_equal(actual_counts, counts)
    assert actual_state == rng.bit_generator.state
    error = float(np.max(np.abs(actual - reference)))
    assert error <= 2e-16
    return {"samples": n, "seed": seed, "gap_batch": batch,
            "activation_counts": counts, "final_rng_state_matches": True,
            "maximum_value_error": error}


def check_exact_small_ecdf() -> dict:
    curve, _ = generator.thinned_ecdf(np.array([0.2, 0.2, 0.3, 0.5, 0.5]))
    expected = np.array([[0, 0], [.2, 0], [.2, .4], [.3, .4],
                         [.3, .6], [.5, .6], [.5, 1], [.5, 1]])
    assert np.array_equal(curve, expected)
    return {"input": [0.2, 0.2, 0.3, 0.5, 0.5], "exact_coordinates": curve.tolist()}


def check_arithmetic_against_rationals() -> dict:
    limit = 100
    primes, actual = generator.arithmetic_values(limit)
    reference_primes = [p for p in range(3, limit + 1, 2)
                        if all(p % d for d in range(2, math.isqrt(p) + 1))]
    reference = []
    for p in reference_primes:
        totient = sum(math.gcd(k, p - 1) == 1 for k in range(1, p))
        value = Fraction(totient, p - 1)
        for j in range(1, 21):
            value *= 1 - Fraction(1, p**j)
        reference.append(float(value))
    assert np.array_equal(primes, reference_primes)
    error = float(np.max(np.abs(actual - reference)))
    assert error <= 3e-16
    return {"prime_limit": limit, "odd_prime_count": len(reference_primes),
            "reference": "Trial division, gcd totient count, exact rational 20-factor product",
            "maximum_value_error": error}


def check_complete_rerun() -> dict:
    script = HERE / "generate_cdf_overlay.py"
    records = []
    with tempfile.TemporaryDirectory(prefix="cdf-generator-validation-") as directory:
        for name in ("first", "second"):
            destination = Path(directory) / name
            command = [sys.executable, str(script), "--output-dir", str(destination),
                       "--samples", "200000", "--odd-primes", "2",
                       "--prime-limit", "1000", "--save-samples"]
            subprocess.run(command, check=True, capture_output=True, text=True)
            records.append(json.loads((destination / "generation_record.json").read_text()))
    first, second = records
    for key in ("array_sha256", "sampling", "checks", "artifacts"):
        assert first[key] == second[key], f"Nondeterministic {key}"
    return {"samples": 200_000, "odd_primes": 2, "arithmetic_prime_limit": 1000,
            "matched": ["array hashes", "RNG states", "diagnostics", "all artifact hashes"],
            "artifacts": first["artifacts"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=HERE / "validation_report.json")
    args = parser.parse_args()
    report = {
        "status": "passed",
        "validated_utc": datetime.now(timezone.utc).isoformat(),
        "environment": {"python": sys.version, "numpy": np.__version__,
                        "matplotlib": generator.matplotlib.__version__, "platform": platform.platform()},
        "generator_sha256": generator.sha256_file(HERE / "generate_cdf_overlay.py"),
        "validator_sha256": generator.sha256_file(Path(__file__)),
        "checks": {
            "two_prime_law": check_two_prime_law(),
            "geometric_refill_scalar_replay": check_refill_against_scalar_replay(),
            "duplicate_and_endpoint_ecdf": check_exact_small_ecdf(),
            "arithmetic_rational_reference": check_arithmetic_against_rationals(),
            "complete_small_run_reproducibility": check_complete_rerun(),
        },
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"All five validation checks passed; report: {args.report.resolve()}")


if __name__ == "__main__":
    main()
