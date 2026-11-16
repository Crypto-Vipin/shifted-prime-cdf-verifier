Vipin Singh Sehrawat

*Arithmetic limiting laws over shifted primes: Hausdorff dimension zero, Fourier decay, and absolute continuity*

[arXiv:2603.11196](https://arxiv.org/abs/2603.11196)

Python 3.11 or later. Run from this directory:

```sh
python3 -m pip install -r requirements.txt
python3 verify_cdf.py
python3 generate_cdf_overlay.py
python3 validate_cdf_generator.py
```

`verify_cdf.py` checks the four certified CDF enclosures.
`generate_cdf_overlay.py` creates the figure, CSV data, and generation record in
`generated/`, using all odd primes up to 1,000,000 and 10,000,000 Monte Carlo
samples over the first 10,000 odd primes, with PCG64 seed 42.
`validate_cdf_generator.py` runs five implementation and reproducibility checks.
Monte Carlo diagnostics are not certified bounds.
