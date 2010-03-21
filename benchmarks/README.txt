Before you make any comparisons to ASOL:

- Keep in mind that ASOL was developed to solve distillation problems. These
  problems are large and sparse. ASOL will show relatively poor performance on
  smaller problems WITH THE DEFAULT SETTINGS. If you wish to make comparisons,
  solve the extr30 and mss20 benchmarks first with your solver, with its default
  settings. Then, use the same settings when solving the tiny problems.
  
- You should compile ASOL optimized to your platform. The binary distributions
  available from the website of ASOL are not appropriate for comparisons. It is
  optimized for size (small) and not for speed. The HSL subroutines would be
  faster than MUMPS but they cannot be included due to license limiations.
  
- ASOL is a proof-of-concept implementation, with a brute-force pruning step. 
  There are several ways to improve its efficiency.
