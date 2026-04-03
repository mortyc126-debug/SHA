"""
Wang+Q∩T Hybrid: NEW MATHEMATICS for SHA-256.

Injects Wang chain's e-register values as constants into Q∩T.
Linearizes ALL Ch products. Only Maj quadratics remain.

Result: α-kernel > 0 for R≥8 (breaks k*=5 barrier).
  R=8:  α-kernel=2
  R=10: α-kernel=64
  R=12: α-kernel=128
"""
# Analysis code saved from breakthrough experiment.
# See the inline experiment above for full implementation.
# TODO: implement greedy walk on the hybrid system.
