Before you make any comparisons to ASOL:

- Keep in mind that ASOL was developed to solve distillation problems. These
  problems are large and sparse. ASOL will show relatively poor performance on
  tiny problems WITH THE DEFAULT SETTINGS. If you wish to make comparisons,
  solve the extr30 and mss20 benchmarks first with your solver, with its default
  settings. Then, use the same settings when solving the tiny problems.
  
- You should compile ASOL optimized to your platform. The binary distributions
  available from the website of ASOL are not appropriate for comparisons.
  
- ASOL is a proof-of-concept implementation, with a brute-force pruning step. 
  There are several ways to improve its efficiency.
