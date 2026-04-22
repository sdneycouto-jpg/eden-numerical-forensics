"""
ÉDEN Numerical Forensics
========================
Cryptographic key auditing toolkit.

Detects 5 classes of defective RSA keys:
  1. Fermat-weak (small |p - q|)
  2. Shared-prime (common factor with another key)
  3. p-1-smooth (vulnerable to Pollard p-1)
  4. p+1-smooth (vulnerable to Williams p+1)
  5. Small-factor (trivially factorable)

Does NOT break well-formed RSA-2048 keys. Nothing does without NFS + cluster.
"""

__version__ = "0.1.0"
__author__ = "André Philipsson (EU SOU ÉDEN)"
__license__ = "MIT"

from .arcanum import fermat_arcanum, fermat_qr_sieve
from .rho import pollard_rho_dna
from .orquestrador import orquestrador_factor, classify_genesis
from .auditor import audit_key, audit_batch

__all__ = [
    "fermat_arcanum",
    "fermat_qr_sieve",
    "pollard_rho_dna",
    "orquestrador_factor",
    "classify_genesis",
    "audit_key",
    "audit_batch",
]
