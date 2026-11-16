#!/usr/bin/env python3
"""Vipin Singh Sehrawat
Arithmetic limiting laws over shifted primes: Hausdorff dimension zero, Fourier decay, and absolute continuity
arXiv:2603.11196
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "monograph-matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
AUTHOR = "Vipin Singh Sehrawat"
TITLE = "Arithmetic limiting laws over shifted primes: dimension zero, Fourier decay, and absolute continuity"
ARXIV_ID = "2603.11196"
DEFAULT_SAMPLES = 10_000_000
DEFAULT_ODD_PRIMES = 10_000
DEFAULT_PRIME_LIMIT = 1_000_000
DEFAULT_SEED = 42
GAP_BATCH = 1_048_576
MATRIX_FACTORS = 20
MAX_DISPLAY_RANKS = 5_001


def sha256_file(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def array_sha256(values: np.ndarray, dtype: str) -> str:

    canonical = np.ascontiguousarray(values, dtype=dtype)
    return hashlib.sha256(memoryview(canonical).cast("B")).hexdigest()


def primes_through(limit: int) -> np.ndarray:
    sieve = np.ones(limit + 1, dtype=np.bool_)
    sieve[:2] = False
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    return np.flatnonzero(sieve).astype(np.int64)


def first_odd_primes(count: int) -> np.ndarray:
    limit = max(32, int((count + 2) * (math.log(count + 2) + math.log(math.log(count + 2)))) + 16)
    while True:
        primes = primes_through(limit)
        odd = primes[primes >= 3]
        if len(odd) >= count:
            return odd[:count]
        limit *= 2


def arithmetic_values(limit: int) -> tuple[np.ndarray, np.ndarray]:

    primes = primes_through(limit)
    phi = np.arange(limit + 1, dtype=np.int64)
    for p in primes:
        phi[p::p] -= phi[p::p] // p
    primes = primes[primes >= 3]
    values = phi[primes - 1].astype(np.float64) / (primes - 1)
    inverse = 1.0 / primes.astype(np.float64)
    power = inverse.copy()
    for _ in range(MATRIX_FACTORS):
        values *= 1.0 - power
        power *= inverse
    return primes, values


def sample_products(
    primes: np.ndarray, sample_count: int, seed: int, *, progress: bool = True
) -> tuple[np.ndarray, dict, dict, np.ndarray]:


    rng = np.random.Generator(np.random.PCG64(seed))
    before = rng.bit_generator.state
    logs = np.zeros(sample_count, dtype=np.float64)
    successes = np.empty(len(primes), dtype=np.int64)
    for prime_index, prime in enumerate(primes):
        p = int(prime)
        q = 1.0 / (p - 1)
        jump = math.log1p(-1.0 / p)
        position = -1
        count = 0
        while True:
            expected = (sample_count - 1 - position) * q
            draw_count = min(GAP_BATCH, max(32, math.ceil(1.01 * expected) + 16))
            gaps = rng.geometric(q, size=draw_count)
            positions = np.cumsum(gaps, dtype=np.int64)
            positions += position
            stop = int(np.searchsorted(positions, sample_count, side="left"))
            valid = positions[:stop]
            logs[valid] += jump
            count += stop
            if stop < draw_count:
                break
            position = int(positions[-1])
        successes[prime_index] = count
        if progress and ((prime_index + 1) % 1000 == 0 or prime_index + 1 == len(primes)):
            print(f"Sampled {prime_index + 1:,}/{len(primes):,} prime coordinates", flush=True)
    np.exp(logs, out=logs)
    logs *= 0.5
    return logs, before, rng.bit_generator.state, successes


def finite_moments(primes: np.ndarray) -> tuple[float, float, float]:
    p = primes.astype(np.float64)
    q = 1.0 / (p - 1.0)
    factor = 1.0 - 1.0 / p
    atom = math.exp(math.fsum(math.log1p(-float(x)) for x in q))
    mean = 0.5 * math.exp(math.fsum(math.log1p(float(x)) for x in q * (factor - 1.0)))
    second = 0.25 * math.exp(math.fsum(math.log1p(float(x)) for x in q * (factor * factor - 1.0)))
    return atom, mean, second


def thinned_ecdf(sorted_values: np.ndarray) -> tuple[np.ndarray, int]:


    n = len(sorted_values)
    stride = max(1, math.ceil(n / (MAX_DISPLAY_RANKS - 1)))
    ranks = np.unique(np.r_[np.arange(0, n, stride, dtype=np.int64), n - 1])
    x = np.unique(sorted_values[ranks])
    left = np.searchsorted(sorted_values, x, side="left") / n
    right = np.searchsorted(sorted_values, x, side="right") / n
    points = np.empty((2 * len(x) + 2, 2), dtype=np.float64)
    points[0] = (0.0, 0.0)
    points[1:-1:2, 0] = x
    points[2:-1:2, 0] = x
    points[1:-1:2, 1] = left
    points[2:-1:2, 1] = right
    points[-1] = (0.5, 1.0)
    if not (np.all(np.diff(points[:, 0]) >= 0) and np.all(np.diff(points[:, 1]) >= 0)):
        raise AssertionError("Thinned ECDF is not monotone")
    return points, stride


def draw_figure(empirical: np.ndarray, monte_carlo: np.ndarray, destination: Path) -> None:
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "cm", "font.size": 10.5,
        "axes.linewidth": 0.7, "pdf.fonttype": 42, "path.simplify": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.major.width": 0.65, "ytick.major.width": 0.65,
    })
    fig = plt.figure(figsize=(343.44 / 72.0, 242.64 / 72.0))
    ax = fig.add_axes([42.014274 / 343.44, 36.396349 / 242.64,
                       (336.24 - 42.014274) / 343.44, (235.44 - 36.396349) / 242.64])
    ax.plot(empirical[:, 0], empirical[:, 1], color="royalblue", lw=1.35,
            label="empirical CDF")
    ax.plot(monte_carlo[:, 0], monte_carlo[:, 1], color="crimson", lw=1.35,
            linestyle=(0, (8, 5)), label="Monte Carlo CDF approximation")
    ax.set(xlim=(0, 0.52), ylim=(0, 1.02))
    ax.set_xticks(np.arange(0.1, 0.51, 0.1))
    ax.set_yticks(np.arange(0.0, 1.01, 0.2))
    ax.set_xlabel(r"$\alpha$", fontsize=14, labelpad=5)
    ax.set_ylabel(r"$G(\alpha)$", fontsize=14, labelpad=4)
    ax.legend(loc="upper left", frameon=False, fontsize=9.3, handlelength=2.7,
              borderaxespad=0.5, labelspacing=0.3)
    fig.savefig(destination / "cdf_overlay.pdf", metadata={
        "Author": AUTHOR, "Title": TITLE, "Subject": f"arXiv:{ARXIV_ID}",
        "Creator": None, "Producer": None,
        "CreationDate": None, "ModDate": None,
    })
    fig.savefig(destination / "cdf_overlay.png", dpi=200, metadata={
        "Author": AUTHOR, "Title": TITLE, "Description": f"arXiv:{ARXIV_ID}",
        "Software": None,
    })
    plt.close(fig)


def run(args: argparse.Namespace) -> None:
    started = datetime.now(timezone.utc)
    clock_start = time.perf_counter()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    head = first_odd_primes(args.odd_primes)
    if args.odd_primes == DEFAULT_ODD_PRIMES and int(head[-1]) != 104743:
        raise AssertionError("Wrong odd-prime cutoff")
    print("Computing arithmetic sample", flush=True)
    empirical_primes, empirical = arithmetic_values(args.prime_limit)
    if args.prime_limit == DEFAULT_PRIME_LIMIT and len(empirical) != 78497:
        raise AssertionError("Wrong arithmetic prime count")
    print(f"Sampling {args.samples:,} products with PCG64({args.seed})", flush=True)
    samples, rng_before, rng_after, successes = sample_products(head, args.samples, args.seed)
    if not np.all(np.isfinite(samples) & (samples > 0) & (samples <= 0.5)):
        raise AssertionError("Invalid simulated support")
    atom, expected_mean, expected_second = finite_moments(head)
    mean = float(np.mean(samples))
    endpoint_count = int(np.count_nonzero(samples == 0.5))
    mean_se = math.sqrt((expected_second - expected_mean**2) / args.samples)
    atom_se = math.sqrt(atom * (1 - atom) / args.samples)
    mean_z = (mean - expected_mean) / mean_se
    atom_z = (endpoint_count / args.samples - atom) / atom_se

    if abs(mean_z) > 6 or abs(atom_z) > 6:
        raise AssertionError(f"Monte Carlo check exceeds six standard errors: mean={mean_z}, atom={atom_z}")
    print(f"Statistical checks: mean z={mean_z:.3f}, endpoint z={atom_z:.3f}", flush=True)
    samples_hash_unsorted = array_sha256(samples, "<f8")
    empirical_hash_unsorted = array_sha256(empirical, "<f8")
    samples.sort()
    empirical.sort()
    empirical_curve, empirical_stride = thinned_ecdf(empirical)
    monte_carlo_curve, simulation_stride = thinned_ecdf(samples)
    np.savetxt(output / "empirical_curve.csv", empirical_curve, delimiter=",",
               header="alpha,ecdf", comments="", fmt="%.17g")
    np.savetxt(output / "monte_carlo_curve.csv", monte_carlo_curve, delimiter=",",
               header="alpha,ecdf", comments="", fmt="%.17g")
    np.savetxt(output / "prime_success_counts.csv", np.column_stack((head, successes)),
               delimiter=",", header="prime,activated_sample_count", comments="", fmt="%d")
    if args.save_samples:
        np.save(output / "sorted_monte_carlo_samples.npy", samples, allow_pickle=False)
        np.save(output / "sorted_arithmetic_values.npy", empirical, allow_pickle=False)
    draw_figure(empirical_curve, monte_carlo_curve, output)
    artifacts = ["cdf_overlay.pdf", "cdf_overlay.png", "empirical_curve.csv",
                 "monte_carlo_curve.csv", "prime_success_counts.csv"]
    if args.save_samples:
        artifacts += ["sorted_monte_carlo_samples.npy", "sorted_arithmetic_values.npy"]
    record = {
        "schema": "monograph-cdf-generation-v1",
        "monograph": {"author": AUTHOR, "title": TITLE, "arxiv_id": ARXIV_ID},
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter() - clock_start,
        "environment": {"python": sys.version, "numpy": np.__version__,
                        "matplotlib": matplotlib.__version__, "platform": platform.platform(),
                        "machine": platform.machine(), "float_dtype": "IEEE-754 binary64"},
        "source": {"generator_sha256": sha256_file(Path(__file__))},
        "parameters": {"samples": args.samples, "odd_prime_count": args.odd_primes,
                       "first_prime": int(head[0]), "last_prime": int(head[-1]),
                       "arithmetic_prime_limit": args.prime_limit,
                       "arithmetic_prime_count": len(empirical_primes),
                       "matrix_product_factors": MATRIX_FACTORS, "seed": args.seed,
                       "generator": "numpy.random.Generator(numpy.random.PCG64(seed))"},
        "sampling": {"method": "Independent geometric waiting times for each ascending prime; q=1/(p-1), support {1,2,...}, position initially -1.",
                     "maximum_gap_batch": GAP_BATCH,
                     "refill_size": "min(1048576, max(32, ceil(1.01*(samples-1-position)*q)+16))",
                     "overshoot": "After first position >= samples, all remaining gaps in that batch are discarded; next prime starts at position -1 using the same RNG stream.",
                     "accumulation": "float64 log products, math.log1p(-1.0/p), ascending prime order; X=0.5*numpy.exp(log_product)",
                     "state_before": rng_before, "state_after": rng_after},
        "thinning": {"method": "Select equally spaced sample ranks plus final rank; retain both exact left and right ECDF counts at every selected x; add (0,0) and (0.5,1). Connect retained points linearly for display.",
                     "maximum_rank_targets": MAX_DISPLAY_RANKS,
                     "arithmetic_rank_stride": empirical_stride,
                     "simulation_rank_stride": simulation_stride,
                     "arithmetic_plot_points": len(empirical_curve),
                     "simulation_plot_points": len(monte_carlo_curve),
                     "pdf_dimensions_pt": [343.44, 242.64]},
        "checks": {"support": "all simulated samples finite and in (0,0.5]",
                   "cdf_monotonicity": "passed",
                   "expected_endpoint_mass": atom, "observed_endpoint_count": endpoint_count,
                   "observed_endpoint_mass": endpoint_count / args.samples,
                   "endpoint_fraction_standard_error": atom_se, "endpoint_z_score": atom_z,
                   "expected_mean": expected_mean, "observed_mean": mean,
                   "expected_second_moment": expected_second,
                   "mean_standard_error": mean_se, "mean_z_score": mean_z,
                   "statistical_check": "Both absolute z scores <=6. Statistical diagnostics, not certified probability enclosures.",
                   "simple_wasserstein_upper_bound": 1 / (2 * int(head[-1])),
                   "finite_head_mean_over_y_log_y": expected_mean / (int(head[-1]) * math.log(int(head[-1])))},
        "array_sha256": {"odd_primes_int64_le": array_sha256(head, "<i8"),
                         "arithmetic_primes_int64_le": array_sha256(empirical_primes, "<i8"),
                         "arithmetic_values_float64_le_in_prime_order": empirical_hash_unsorted,
                         "monte_carlo_float64_le_in_sample_order": samples_hash_unsorted,
                         "sorted_monte_carlo_float64_le": array_sha256(samples, "<f8")},
        "artifacts": {name: {"sha256": sha256_file(output / name), "bytes": (output / name).stat().st_size}
                      for name in artifacts},
    }
    (output / "generation_record.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"Generated figure, data, and record in {output}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "generated")
    parser.add_argument("--samples", type=int, default=DEFAULT_SAMPLES)
    parser.add_argument("--odd-primes", type=int, default=DEFAULT_ODD_PRIMES)
    parser.add_argument("--prime-limit", type=int, default=DEFAULT_PRIME_LIMIT)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--save-samples", action="store_true", help="Also save full sorted arrays (about 81 MB at defaults)")
    args = parser.parse_args()
    if args.samples < 100 or args.odd_primes < 2 or args.prime_limit < 5 or args.seed < 0:
        parser.error("Require samples>=100, odd-primes>=2, prime-limit>=5, and seed>=0")
    run(args)


if __name__ == "__main__":
    main()
