#!/usr/bin/env python3
"""Deterministically verify the four certified CDF enclosures of the paper:
*The shifted-prime Erdős–Wintner law for primitive-root determinant
densities: extremal order, dimension zero, and Fourier decay* (V. S. Sehrawat).

The verification follows the proof of the certified-enclosure proposition:

* split off the Bernoulli variable at p=3;
* convolve the head primes 5 <= p <= 50,000 on a monotone dyadic grid;
* control the non-negative tail by the stated Bennett bound;
* at beta=3.0001, use the stated forced-zero event and Markov bound with
  secondary cutoff 2,000,000.

Logarithms used to choose grid cells are enclosed by exact rational atanh
series to more than 200 bits.  Probability masses are propagated in IEEE-754
binary64; the reported CDF bounds are widened by 1e-9.  That pad is
deliberately generous: the head recurrence is an l^1-contraction whose
standard binary64 forward-error bound is under 128u per prime, with
u = 2^-53, hence below 7.4e-11 over the 5,131 head primes.

Exit codes: 0 verified, 1 containment failure, 2 usage error, 3 output could
not be written.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from pathlib import Path

import numpy as np


HEAD_CUTOFF = 50_000
TAIL_CUTOFF = 80_000
SECONDARY_CUTOFF = 2_000_000
GRID_BITS = 20
GRID_SCALE = 1 << GRID_BITS
LOG_BITS = 220
SHIFT = Fraction(3, 20_000)  # 0.00015 in the Bennett lower bound
FIXED_SCALE = 10**60
FLOAT_CDF_PAD = Decimal("1e-9")
TRANSCENDENTAL_PAD = Decimal("1e-40")

PUBLISHED = {
    "3": (Decimal("0.49146"), Decimal("0.49149")),
    "4": (Decimal("0.95378"), Decimal("0.95385")),
    "5": (Decimal("0.99922"), Decimal("0.99924")),
    "3.0001": (Decimal("0.5559"), Decimal("0.5645")),
}

# beta, beta/2 (the first split threshold), beta/3 (the second threshold)
TARGETS = (
    ("3", Fraction(3, 2), Fraction(1, 1)),
    ("4", Fraction(2, 1), Fraction(4, 3)),
    ("5", Fraction(5, 2), Fraction(5, 3)),
    ("3.0001", Fraction(30_001, 20_000), Fraction(30_001, 30_000)),
)


def sieve_primes(limit: int) -> list[int]:
    """Return all primes <= limit by a deterministic Eratosthenes sieve."""
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            count = (limit - start) // p + 1
            sieve[start : limit + 1 : p] = b"\x00" * count
    return [n for n in range(2, limit + 1) if sieve[n]]


def log_ratio_interval(ratio: Fraction, bits: int = LOG_BITS) -> tuple[Fraction, Fraction]:
    """Enclose log(ratio) by an exact rational atanh series.

    For ratio >= 1, put z=(ratio-1)/(ratio+1).  Then

        log(ratio) = 2 * sum_{k>=0} z^(2k+1)/(2k+1).

    The omitted positive tail is at most the first omitted power divided by
    its denominator and by 1-z^2.  Iteration continues until its width is
    below 2^-bits.
    """
    if ratio < 1:
        lo, hi = log_ratio_interval(1 / ratio, bits)
        return -hi, -lo
    if ratio == 1:
        return Fraction(0), Fraction(0)

    z = (ratio - 1) / (ratio + 1)
    z2 = z * z
    term = z
    partial = Fraction(0)
    target = Fraction(1, 1 << bits)
    k = 0
    while True:
        partial += 2 * term / (2 * k + 1)
        first_omitted = term * z2
        tail = 2 * first_omitted / ((2 * k + 3) * (1 - z2))
        if tail < target:
            return partial, partial + tail
        term = first_omitted
        k += 1


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def certified_grid_floor(
    ratio: Fraction, subtract: Fraction = Fraction(0)
) -> tuple[int, Fraction, Fraction]:
    """Return floor(GRID_SCALE*(log(ratio)-subtract)) with a certificate."""
    bits = LOG_BITS
    while True:
        lo, hi = log_ratio_interval(ratio, bits)
        lo -= subtract
        hi -= subtract
        floor_lo = floor_fraction(lo * GRID_SCALE)
        floor_hi = floor_fraction(hi * GRID_SCALE)
        if floor_lo == floor_hi:
            return floor_lo, lo, hi
        bits += 64


def ceil_div(a: int, b: int) -> int:
    return -(-a // b)


def fixed_upper_sum(primes: list[int], lo: int, hi: int, power: int) -> Fraction:
    """Exact fixed-point upper bound for sum 1/(p-1)^power, lo < p <= hi."""
    numerator = 0
    for p in primes:
        if p <= lo:
            continue
        if p > hi:
            break
        numerator += ceil_div(FIXED_SCALE, (p - 1) ** power)
    return Fraction(numerator, FIXED_SCALE)


def fixed_product_lower(primes: list[int], hi_exclusive: int) -> Fraction:
    """Fixed-point lower bound for product_{5<=p<hi} (p-2)/(p-1)."""
    value = FIXED_SCALE
    for p in primes:
        if p < 5:
            continue
        if p >= hi_exclusive:
            break
        value = value * (p - 2) // (p - 1)
    return Fraction(value, FIXED_SCALE)


def fraction_to_decimal(value: Fraction, precision: int = 80) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = precision
        return Decimal(value.numerator) / Decimal(value.denominator)


def bennett_upper(mean_upper: Fraction, variance_upper: Fraction, b_upper: Fraction) -> Decimal:
    """Bennett upper bound for P(R_T > SHIFT), enlarged by a decimal pad."""
    t_lower = SHIFT - mean_upper
    if t_lower <= 0:
        return Decimal(1)
    with localcontext() as ctx:
        ctx.prec = 80
        t = fraction_to_decimal(t_lower)
        v = fraction_to_decimal(variance_upper)
        b = fraction_to_decimal(b_upper)
        u = b * t / v
        h = (Decimal(1) + u) * (Decimal(1) + u).ln() - u
        bound = (-(v / (b * b)) * h).exp()
        return bound + TRANSCENDENTAL_PAD


@dataclass(frozen=True)
class HeadBounds:
    lower: dict[tuple[str, str], Decimal]
    upper: dict[tuple[str, str], Decimal]
    total_mass_lower_grid: Decimal
    total_mass_upper_grid: Decimal


def head_grid_bounds(head_primes: list[int]) -> HeadBounds:
    """Compute monotone lower/upper head-CDF bounds on the dyadic grid."""
    if np.finfo(np.float64).eps != 2.0**-52:
        raise RuntimeError("This verifier requires IEEE-754 binary64.")

    query_indices: dict[tuple[str, str], int] = {}
    max_threshold = Fraction(0)
    for name, first_ratio, second_ratio in TARGETS:
        for which, ratio in (("first", first_ratio), ("second", second_ratio)):
            index, lo, _ = certified_grid_floor(ratio)
            query_indices[(name, which)] = index
            max_threshold = max(max_threshold, lo)
            if which == "first" or name in {"4", "5"}:
                shifted_index, shifted_lo, _ = certified_grid_floor(ratio, SHIFT)
                query_indices[(name, which + "_shifted")] = shifted_index
                max_threshold = max(max_threshold, shifted_lo)

    # The largest requested threshold is log(5/2)<0.917.  A small exact
    # guard leaves all queried cells inside the retained non-negative grid.
    grid_length = floor_fraction(max_threshold * GRID_SCALE) + 4

    floor_shifts: list[int] = []
    ceil_shifts: list[int] = []
    for p in head_primes:
        shift_floor, _, _ = certified_grid_floor(Fraction(p, p - 1))
        floor_shifts.append(shift_floor)
        ceil_shifts.append(shift_floor + 1)

    def propagate(shifts: list[int]) -> tuple[np.ndarray, Decimal]:
        mass = np.zeros(grid_length, dtype=np.float64)
        mass[0] = 1.0
        for p, shift in zip(head_primes, shifts):
            q = np.float64(1.0 / (p - 1))
            old = mass.copy()
            mass *= np.float64(1.0) - q
            if shift < grid_length:
                mass[shift:] += old[: grid_length - shift] * q
        total = Decimal(repr(math.fsum(float(x) for x in mass)))
        return mass, total

    # R_up >= R_H gives the lower CDF bound; R_down <= R_H gives the upper.
    lower_mass, lower_total = propagate(ceil_shifts)
    upper_mass, upper_total = propagate(floor_shifts)

    def collect(mass: np.ndarray, direction: str) -> dict[tuple[str, str], Decimal]:
        result: dict[tuple[str, str], Decimal] = {}
        for key, index in query_indices.items():
            if index < 0:
                raw = Decimal(0)
            else:
                stop = min(index + 1, len(mass))
                raw = Decimal(repr(math.fsum(float(x) for x in mass[:stop])))
            if direction == "lower":
                result[key] = max(Decimal(0), raw - FLOAT_CDF_PAD)
            else:
                result[key] = min(Decimal(1), raw + FLOAT_CDF_PAD)
        return result

    return HeadBounds(
        lower=collect(lower_mass, "lower"),
        upper=collect(upper_mass, "upper"),
        total_mass_lower_grid=lower_total,
        total_mass_upper_grid=upper_total,
    )


def format_decimal(
    value: Decimal, places: int = 12, *, rounding: str | None = None
) -> str:
    quantum = Decimal(1).scaleb(-places)
    if rounding is None:
        return str(value.quantize(quantum))
    return str(value.quantize(quantum, rounding=rounding))


def verify() -> str:
    primes = sieve_primes(SECONDARY_CUTOFF)
    head_primes = [p for p in primes if 5 <= p <= HEAD_CUTOFF]
    head = head_grid_bounds(head_primes)

    mean_upper = (
        fixed_upper_sum(primes, HEAD_CUTOFF, TAIL_CUTOFF, 2)
        + Fraction(1, TAIL_CUTOFF - 1)
    )
    variance_upper = (
        fixed_upper_sum(primes, HEAD_CUTOFF, TAIL_CUTOFF, 3)
        + Fraction(1, 2 * (TAIL_CUTOFF - 1) ** 2)
    )
    b_upper = Fraction(1, HEAD_CUTOFF)
    tail_probability = bennett_upper(mean_upper, variance_upper, b_upper)

    # At beta=3.0001, theta=log(30001/30000).  Terms with p<30001
    # individually exceed theta and are forced to vanish.  Conditional on
    # that event, Markov is applied only to the remaining non-negative sum.
    small_ratio = Fraction(30_001, 30_000)
    theta_lo, _ = log_ratio_interval(small_ratio)
    forced_zero_lower = fixed_product_lower(primes, 30_001)
    remaining_mean_upper = (
        fixed_upper_sum(primes, 30_000, SECONDARY_CUTOFF, 2)
        + Fraction(1, SECONDARY_CUTOFF - 1)
    )
    markov_factor_lower = max(
        Fraction(0), Fraction(1) - remaining_mean_upper / theta_lo
    )
    small_cdf_lower = fraction_to_decimal(
        forced_zero_lower * markov_factor_lower
    ) - TRANSCENDENTAL_PAD

    computed: dict[str, tuple[Decimal, Decimal]] = {}
    for name, _, _ in TARGETS:
        first_lower = max(
            Decimal(0),
            head.lower[(name, "first_shifted")] - tail_probability,
        )
        first_upper = head.upper[(name, "first")]

        if name == "3":
            # The infinite remaining sum is positive almost surely.
            second_lower = Decimal(0)
            second_upper = Decimal(0)
        elif name == "3.0001":
            second_lower = max(Decimal(0), small_cdf_lower)
            second_upper = head.upper[(name, "second")]
        else:
            second_lower = max(
                Decimal(0),
                head.lower[(name, "second_shifted")] - tail_probability,
            )
            second_upper = head.upper[(name, "second")]

        computed[name] = (
            (first_lower + second_lower) / 2,
            (first_upper + second_upper) / 2,
        )

    lines = [
        "Deterministic verification of the four certified CDF enclosures",
        f"head_cutoff={HEAD_CUTOFF}",
        f"tail_cutoff={TAIL_CUTOFF}",
        f"secondary_cutoff={SECONDARY_CUTOFF}",
        f"grid=2^-{GRID_BITS}",
        f"log_enclosure_bits={LOG_BITS}",
        f"head_prime_count={len(head_primes)}",
        f"float_cdf_pad={FLOAT_CDF_PAD}",
        f"head_retained_mass_up_grid={format_decimal(head.total_mass_lower_grid, 15)}",
        f"head_retained_mass_down_grid={format_decimal(head.total_mass_upper_grid, 15)}",
        f"tail_mean_upper={format_decimal(fraction_to_decimal(mean_upper), 18)}",
        f"tail_variance_upper={format_decimal(fraction_to_decimal(variance_upper), 24)}",
        f"bennett_shift={format_decimal(fraction_to_decimal(SHIFT), 8)}",
        f"bennett_tail_probability_upper={format_decimal(tail_probability, 15)}",
        f"forced_zero_probability_lower={format_decimal(fraction_to_decimal(forced_zero_lower), 15)}",
        f"small_tail_mean_upper={format_decimal(fraction_to_decimal(remaining_mean_upper), 18)}",
        f"small_threshold_cdf_lower={format_decimal(small_cdf_lower, 15)}",
        "",
        "beta | verifier enclosure | manuscript enclosure | result",
    ]

    all_pass = True
    for name, _, _ in TARGETS:
        lower, upper = computed[name]
        pub_lower, pub_upper = PUBLISHED[name]
        passed = lower >= pub_lower and upper <= pub_upper
        all_pass &= passed
        lines.append(
            f"{name} | "
            f"[{format_decimal(lower, rounding=ROUND_FLOOR)}, "
            f"{format_decimal(upper, rounding=ROUND_CEILING)}] "
            f"| [{pub_lower}, {pub_upper}] | {'PASS' if passed else 'FAIL'}"
        )

    lines.extend(
        [
            "",
            "ALL FOUR ENCLOSURES VERIFIED" if all_pass else "VERIFICATION FAILED",
        ]
    )
    if not all_pass:
        raise AssertionError("\n".join(lines))
    return "\n".join(lines) + "\n"


def main() -> int:
    # Validate arguments before the ~10 s verification, so a usage error is
    # reported immediately.
    path: Path | None = None
    if len(sys.argv) == 3 and sys.argv[1] == "--output":
        path = Path(sys.argv[2])
    elif len(sys.argv) != 1:
        sys.stderr.write(f"usage: {Path(sys.argv[0]).name} [--output PATH]\n")
        return 2

    output = verify()

    if path is None:
        sys.stdout.write(output)
        return 0
    try:
        path.write_text(output, encoding="utf-8")
    except OSError as exc:
        # Distinct from exit 1, which signals a containment failure.
        sys.stderr.write(f"error: could not write output: {exc}\n")
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())