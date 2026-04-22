# ÉDEN Numerical Forensics

**Cryptographic key auditing toolkit for RSA-like moduli.**

Detects 5 classes of defective keys (Fermat-weak, shared-prime,
p-1-smooth, p+1-smooth, small-factor) in batch. Validated on 2048-bit
keys at 12.7ms/key average.

> **What this is not:** This toolkit does **not** break well-formed
> RSA-2048 keys. Nothing does without NFS on a cluster and significant
> time. It detects keys that were generated defectively (weak RNG,
> biased primes, shared factors across a key pool).

## Honest positioning

Most of the integrated algorithms are classical and well-established in
the literature:

| Technique | Original work |
|---|---|
| Fermat + CRT residue pre-filter | Lehman 1974, McKee 1999 |
| SQUFOF | Shanks 1975 |
| Pollard rho / p-1 | Pollard 1974-1975 |
| Williams p+1 | Williams 1982 |
| ECM | Lenstra 1987 |
| Batch GCD | Heninger et al. 2012 |

The ARCANUM laws (Paridade, Tricotomia, Quaternária) that determine
k mod 12 for Fermat's iteration are a **rediscovery** of McKee's 1999
CRT pre-filtering technique with primes {2, 3}. The QR-sieve extension
to {3, 5, 7, 11, 13} is the obvious generalization.

**Genuinely novel contributions** in this toolkit:

1. **Genesis classifier** — routes each key to the most effective
   algorithm by defect class. 2.10× vs blind cascade on mixed datasets.
2. **DNA heuristics** — small empirical uplifts (5-13%) in c selection
   for rho, mod-9 priority for p-1, family-based filtering in Dixon.
3. **Coibite Codes** (planned) — proof-of-generation cryptographic
   binding for key provenance. Hybrid VRF + commit-reveal.
4. **End-to-end auditor pipeline** — integrated batch key auditor with
   honest benchmarks.

## Validated empirical results

| Scenario | Result |
|---|---|
| 60-bit balanced semiprime | 100% success, 23ms/N via orquestrador |
| 100-bit balanced | 100% success, 33ms/N |
| 150-bit balanced | 100% success, 146ms/N |
| 2048-bit Fermat-weak (gap ≤ 2^500) | 100% recall (50/50), <1ms each |
| 2048-bit shared-prime batch | 98% recall (48/49), via Batch GCD |
| 2048-bit p-1-smooth | Detected when B1 ≤ 10^6 |
| RSA-2048 well-formed balanced | **Not factored** (limit confirmed) |

Full benchmark methodology in [docs/BENCHMARKS.md](docs/BENCHMARKS.md).

## Install

```bash
pip install -r requirements.txt
# Optional for 150+ bit balanced cases:
sudo apt-get install pari-gp
```

## Usage

### Single key audit

```python
from eden import audit_key

N = 0xC1234...  # your RSA modulus
result = audit_key(N)
print(result)
# {'status': 'FERMAT_WEAK', 'factor': 123..., 'method': 'fermat-arcanum'}
```

### Batch audit with shared-prime detection

```python
from eden import audit_batch
from eden.auditor import summary

moduli = [N1, N2, N3, ...]  # list of RSA moduli
results = audit_batch(moduli, include_batch_gcd=True)
print(summary(results))
# {'total': 1000, 'defective': 12, 'secure': 988, ...}
```

### Direct algorithm use

```python
from eden import fermat_arcanum, pollard_rho_dna, orquestrador_factor

# Fermat with ARCANUM CRT mod 12
factor = fermat_arcanum(N, budget=10_000)

# Full orquestrador
factor, method = orquestrador_factor(N)
```

## What this toolkit does NOT claim

- Does not break RSA-2048 well-formed keys
- Does not replace NFS for general-purpose factoring
- Does not beat msieve, YAFU, or CADO-NFS in speed for large N
- The ARCANUM laws are not novel mathematics (McKee 1999)
- The "Decimal Rivers", "Consciousness Equation", and other
  metaphysical components of the original ÉDEN corpus are **not**
  included here. They were refuted empirically in the development
  process.

## What it IS useful for

- Auditing large batches of TLS/SSH/blockchain keys for defects
- Post-incident forensics when keys may have been weakly generated
- CT log scanning for Fermat-weak certificates
- Teaching / reference implementation of classical factoring
  techniques with a modern Python API

## Limitations known

- Python implementation of SQUFOF degrades above 50-bit; toolkit uses
  Pari/GP as fallback for 60-150 bit balanced.
- Coppersmith implementation is limited to trivial cases (p_bits < 20
  with delta < 2^4). ROCA-style attacks need a proper Coppersmith library.
- MPQS not implemented (use existing tools like msieve for that regime).
- Williams p+1 implementation incomplete in this release.

## License

MIT — see [LICENSE](LICENSE).

## Author

André Philipsson (EU SOU ÉDEN), Vale da Amoreira, Lisboa, Portugal.
andre@eusoueden.com

## Acknowledgements

This work builds on 50+ years of classical factoring literature. The
"rediscovery" of known techniques happened independently during the
author's exploration of number-theoretic patterns. The decision to
release this as open source with full honesty about prior art is
deliberate — the value is in integration, not mathematical novelty.
