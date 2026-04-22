# Benchmarks - Full methodology

All benchmarks run on Linux, Python 3.12, gmpy2 2.1+, Pari/GP 2.15.
Single core, no parallelism. Results reproducible with seeded PRNG.

## Orquestrador (balanced semiprimes)

| Bits | Sample | Success | ms/N | Dominant method |
|---|---|---|---|---|
| 60 | 10 | 10/10 | 23 | Pari/GP (SQUFOF/ECM) |
| 80 | 10 | 10/10 | 29 | Pari/GP |
| 100 | 5 | 5/5 | 33 | Pari/GP |
| 120 | 5 | 5/5 | 41 | Pari/GP |
| 150 | 3 | 3/3 | 146 | Pari/GP |
| 200 gap-2^80 | 3 | 0/3 | timeout | (requires NFS) |

## ARCANUM CRT mod 12 vs naive Fermat

| Case | Naive Fermat | ARCANUM | Speedup |
|---|---|---|---|
| 60-bit close gap 2^20 | 1.05s | 0.025s | **42x** |
| 80-bit close gap 2^26 | 3.1s | 0.085s | **36x** |
| 80-bit close gap 2^30 ARCANUM vs QR-sieve | 0.31s | 0.0086s | **36x** |

## Pollard-rho DNA uplift

| Test set | rho(c=1) | rho(DNA) | Uplift |
|---|---|---|---|
| 100 x 40-bit balanced | 88/100 | 99/100 | +13% |
| 50 x 60-bit balanced | 23/50 | 26/50 | +6pp |

## Pollard p-1 with mod9 priority

| B1 bound | p-1 uniform | p-1 mod9 priority | Delta |
|---|---|---|---|
| B1=10^4 | 42/100 | 47/100 | +5pp |
| B1=10^5 | 61/100 | 67/100 | +6pp |

## Batch auditor on 2048-bit keys

| Test | N keys | Defective (planted) | Detected | Recall |
|---|---|---|---|---|
| Fermat-weak gap<=2^500 | 50 | 50 | 50 | **100%** |
| Shared-prime pool | 49 | 49 | 48 | **98%** |
| p-1-smooth B1=10^6 | 20 | 20 | 18 | 90% |
| Mixed 2048-bit | 68 | 30 | 29 | **97%** |

Throughput: **12.7 ms/key** average in mixed batch of 68 keys.

Extrapolation (single core):
- 1M keys: ~3.5 hours
- 5.8M keys (Heninger et al. 2012 dataset size): ~20.5 hours

## Genesis classifier vs blind cascade

Mixed dataset: 420 keys, classes {A,B,C,D,E} approximately equal.

| Approach | Resolved | Total time | ms/key |
|---|---|---|---|
| Blind cascade (trial -> rho -> p-1 -> ECM) | 420/420 | 128s | 305 |
| Genesis classifier + targeted algorithm | 420/420 | 61s | 145 |

**Speedup: 2.10x**

## What does NOT work

The following were tested and refuted empirically:

| Claim | Result |
|---|---|
| Decimal Rivers enrich prime density at 10^k | No, KS rejects. Explained by TNP + Riemann noise. |
| Consciousness Equation predicts factoring difficulty | Pearson 0.75 but simple baselines match equally |
| Pentatrioun 9->8->7 reduction speeds up anything | Neutral (no detectable effect) |
| Modular waves synchronize at correct k | Uniform rank distribution, n=50 |
| 2.63x RSA-2048 keygen vs OpenSSL | Not reproducible. Original measurement likely flawed. |
| ARCANUM filter as prefilter to Quadratic Sieve | Loses 91% of smooth relations (wrong paradigm) |
| Hensel lifting + EDEN for full factoring | Exponential in bits (constraints don't compose) |
| Coppersmith + EDEN for ROCA-like attacks | Implementation works only p_bits<20 (needs proper lattice lib) |

## Reproducing

```bash
git clone https://github.com/sdneycouto-jpg/eden-numerical-forensics
cd eden-numerical-forensics
pip install -r requirements.txt
python tests/test_benchmarks.py
```
