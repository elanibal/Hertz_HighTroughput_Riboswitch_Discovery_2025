"""Riboswitch per-sequence classification (mechanism / direction / sensitivity).

Extends the Hertz et al. 2026 fluoride pipeline to any metabolite-sensing class.
Pure classification logic (offline, testable): `classify_sequence`, `schema`.
Network I/O (Rfam/NCBI): `ncbi`.
"""

from .classifier import classify_sequence
from . import schema

__all__ = ["classify_sequence", "schema"]
