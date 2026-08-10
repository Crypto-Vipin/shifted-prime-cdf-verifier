# Certified CDF enclosures — deterministic verifier

Independent verification of the four certified CDF enclosures stated in the
paper *The shifted-prime Erdős–Wintner law for primitive-root determinant
densities: extremal order, dimension zero, and Fourier decay* (V. S. Sehrawat).

The verifier is self-certifying: it recomputes each enclosure from scratch and
checks it against the interval printed in the paper. On success it reports
`ALL FOUR ENCLOSURES VERIFIED` and exits 0; on any containment failure it
raises `AssertionError` and exits 1. It does not modify or compile the paper.

## Usage

Requires Python 3.9 or later and NumPy:

```sh
pip install -r requirements.txt
```

Then:

```sh
python3 verify_cdf.py
```

Runtime is about ten seconds. To write the report to a file instead of
standard output:

```sh
python3 verify_cdf.py --output OUTPUT_PATH
```

Exit codes are 0 (verified), 1 (containment failure), 2 (usage error), and
3 (output could not be written). The verifier asserts that its probability
arrays use IEEE-754 binary64. Output is deterministic and carries no version
strings or timestamps, so any compliant IEEE-754 environment reproduces it
byte for byte.

## Discretization

The verifier realizes the discretization fixed in the proof:

| Parameter | Value |
| --- | --- |
| Head cutoff `L` | 50,000 |
| Tail cutoff `M` | 80,000 |
| Secondary cutoff `M*` | 2,000,000 |
| Dyadic grid | `2^-20` |
| Logarithmic enclosure width | below `2^-220` |
| Bennett shift `s` | 0.00015 |
| binary64 allowance | `1e-9` |

## Mathematical certificate

Write

```
R = sum_{p>=5} B_p log(p/(p-1)),    Pr(B_p = 1) = 1/(p-1).
```

Splitting off `B_3` gives the paper's exact identity

```
F(beta) = (1/2) Pr(R <= log(beta/2)) + (1/2) Pr(R <= log(beta/3)).
```

The head law is convolved on the dyadic grid `2^-20 Z_{>=0}`. Rounding every
logarithmic jump upward gives a stochastic upper coupling and therefore a lower
CDF bound; rounding every jump downward gives a lower coupling and therefore an
upper CDF bound. Every logarithm and every queried threshold is enclosed by the
exact rational identity

```
log r = 2 sum_{k>=0} z^(2k+1)/(2k+1),    z = (r-1)/(r+1),
```

with positive remainder bounded by `2 z^(2n+3) / ((2n+3)(1-z^2))`. The
resulting enclosure width is below `2^-220`.

For the non-negative tail `R_T` the verifier uses the paper's bounds

```
mu <= sum_{50000 < p <= 80000} 1/(p-1)^2 + 1/79999,
v  <= sum_{50000 < p <= 80000} 1/(p-1)^3 + 1/(2 * 79999^2),
b  <= 1/50000,
```

and Bennett's inequality at the fixed shift `s = 0.00015`. All rational sums
are rounded upward on an exact `1e-60` fixed-point lattice.

For the secondary threshold at `beta = 3.0001`, terms with
`log(p/(p-1)) > log(30001/30000)` — equivalently `p < 30001` — are forced to
vanish. Their probability is bounded below by exact fixed-point multiplication.
Markov's inequality is then applied only to the remaining sum, whose mean is
bounded using primes through `M* = 2,000,000` and the integral tail
`1/(M*-1)`.

## Floating-point error allowance

Grid shifts and queried thresholds use exact rational logarithm enclosures; the
tail sums and the forced-zero product use exact fixed-point rational bounds.
The Bennett transcendental expression is evaluated with 80-digit decimal
arithmetic and enlarged by an outward pad of `1e-40`. Only the head probability
recurrence is binary64. One update has the non-negative stochastic form

```
y_j = (1-q) x_j + q x_{j-h},    q = 1/(p-1).
```

The exact map is an `l^1`-contraction. Standard binary64 error bounds for the
two coefficient operations, two products, and one addition give less than
`128u` of accumulated `l^1` error per prime, where `u = 2^-53`. Over the 5,131
head primes this is below `7.4e-11`; underflow contributes less than `1e-290`.
The verifier widens every head CDF by `1e-9`, over thirteen times this bound.
`math.fsum` is used for the final non-negative sums.

The resulting verifier intervals lie strictly inside all four intervals printed
in the paper. Lower endpoints are printed rounded downward and upper endpoints
rounded upward, so the printed intervals are themselves valid enclosures.

## Citation

V. S. Sehrawat, *The shifted-prime Erdős–Wintner law for primitive-root
determinant densities: extremal order, dimension zero, and Fourier decay*.

## License

MIT.