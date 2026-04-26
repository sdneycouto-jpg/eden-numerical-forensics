# Validated factorization results

This document records the empirical results of QS-EDEN in C/GMP on
balanced semiprime targets, validated against the GitHub repository's
implementation in `c_module/qs_scale.c`.

All results are from real factorization runs, each yielding correct
prime factors verified by `p * q == N`.

## Results table

| Bits | Time | Example factored |
|------|------|------------------|
| 50 bits | <50ms | (multiple cases) |
| 60 bits | 275ms | 459359383744961369 = 828853211 x 554210779 |
| 70 bits | 1.9s | 720677329433092577089 = 22840560601 x 31552523689 |
| 80 bits | 10.1s | 991564502565306583519079 = 933373627363 x 1062344674733 |
| 90 bits | 57.6s | 540587929957694873795321321 = 24023117810663 x 22502821416367 |

## Configuration parameters

| Bits | B (prime bound) | M (sieve radius) | Factor base size | Smooths found | Kernel relations |
|------|-----------------|------------------|-------------------|---------------|------------------|
| 60 | 2,000 | 100,000 | 146 | 388 | 243 |
| 70 | 5,000 | 300,000 | 350 | 1,774 | 1,425 |
| 80 | 12,000 | 800,000 | 718 | 1,679 | 968 |
| 90 | 30,000 | 2,000,000 | 1,626 | 9,624 | 7,999 |

## Scaling

Time grows roughly 5x per 10 bits, consistent with the QS theoretical
complexity ~exp(c * sqrt(ln N * ln ln N)).

Estimated extrapolations on the same code, single core:

- 100 bits: ~3-5 minutes (memory-bound at large M)
- 110 bits: ~15-25 minutes
- 128 bits: ~30-60 minutes (would need careful tuning)
- 256+ bits: out of QS regime; would need MPQS + large primes
- 512+ bits: out of practical QS; NFS territory
- 1024+ bits: NFS only, cluster compute

## Pipeline executed

For each N:

1. **Build factor base**: primes p <= B where N is a quadratic residue
   mod p (Legendre symbol = 1).
2. **Direct sieve**: for each x in [-M, M], compute Q(x) = (x + isqrt(N))^2 - N
   and try to fully factor over the factor base.
3. **Reduce mod 2**: each exponent vector becomes a binary vector in GF(2).
4. **Gaussian elimination over GF(2)** with `uint64_t` bitsets for speed.
   Find subsets of vectors whose XOR is zero.
5. **For each kernel vector**: compute X = product of a_i mod N,
   Y = product of p^(exp_i / 2) mod N. Try gcd(X - Y, N) and gcd(X + Y, N).
   At least one is a non-trivial factor with high probability.

## Compile and run

```bash
sudo apt install libgmp-dev
cd c_module
make
./qs_scale 60 2000 100000     # 60 bits in ~275ms
./qs_scale 70 5000 300000     # 70 bits in ~2s
./qs_scale 80 12000 800000    # 80 bits in ~10s
./qs_scale 90 30000 2000000   # 90 bits in ~58s
```

## Where this fits in the auditor

The QS-EDEN module covers the gap between the orquestrador's existing
algorithms (Fermat, Pollard rho, p-1) and the regime where Pari/GP
SQUFOF+ECM become inefficient. Specifically:

- Up to 60 bits balanced: orquestrador (Pari/GP SQUFOF) handles this faster
- 60-90 bits balanced: QS-EDEN is competitive
- 90-120 bits balanced: QS-EDEN with tuned parameters
- 120+ bits balanced: would need MPQS (not in this release)

For Fermat-weak, p-1-smooth, shared-prime, and small-factor defects,
the existing orquestrador continues to handle them in milliseconds
regardless of N's bit length.

## Honest limits

This QS implementation is educational/demonstrative, not a competitor
to msieve, YAFU, or CADO-NFS. Those tools have decades of optimization:

- Self-initializing sieve
- Large prime variation (single + double large primes)
- Polynomial selection (MPQS)
- Block Lanczos / Block Wiedemann for matrix step
- Highly tuned assembly inner loops

For production factoring of cryptographic-scale numbers, use those
tools. This module exists to show the pipeline end-to-end and provide
a baseline for the auditor's QS-tier capability.
