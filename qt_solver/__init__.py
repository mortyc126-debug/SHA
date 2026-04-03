"""
Q∩T Algebra Solver for SHA-256
===============================

A global solver framework for the hybrid system:
  Q: 256 quadratic GF(2) equations (SHA-256 with fixed carry)
  T: 14,336 threshold equations (carry self-consistency)
  512 variables (message bits)

Based on methodology_v20.md research identifying Q∩T as the
sole remaining promising direction for SHA-256 cryptanalysis.
"""

__version__ = "0.1.0"
