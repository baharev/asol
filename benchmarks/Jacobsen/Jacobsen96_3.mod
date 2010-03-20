#################### DESCRIPTION ###############################################
#
# Benchmarks Jacobsen59_2 and Jacobsen96_3 
#
# Last updated: 20 Mar 2010
#
# The Jacobsen96_3.mod has 5 solutions, Jacobsen59_2.mod has 3 solutions. They
# correspond to the steady states of a methanol-propanol column (mass reflux,
# energy-balances included) discussed in:
#
# E. W. Jacobsen, S. Skogestad;
# Multiple Steady States in Ideal Two-Product Distillation;
# AIChE Journal, 1991, 37, 499-511.
# 
# One solution is infeasible in practice (negative flow rates).
# 
# Solving this problem with the interval methods is discussed in:
# 
# A. Baharev, L. Kolev, E. Rév;
# Computing multiple steady states in homogeneous azeotropic and ideal
# two-product distillation
# (preprint available from http://reliablecomputing.eu)
#
# High resolution bifurcation diagrams are presented in the above manuscript.
#
################## PROBLEM STATISTICS ##########################################
#
# Presolve eliminates 0 constraints and 1 variable.
# Substitution eliminates 24 variables.
# Adjusted problem:
# 19 variables:
#         16 nonlinear variables
#         3 linear variables
# 19 constraints, all nonlinear; 90 nonzeros
# 0 objectives.
#
################# PARAMETERS FROM SPECIFICATION ################################

param Lw  := 96.0;
param V_N := 3.0;

param N   := 8;
param N_F := 5;

param  F{1..N};
param  z{1..N};
param qF{1..N};

param alpha := 3.55;
param M1 := 32.04;
param M2 := 60.10;

param V_L{1..N}; param V_U{1..N};

############### VARIABLES ######################################################

var x{1..N} >= 1.0e-4, <= 1.0;
var V{j in 1..N} >= V_L[j], <= V_U[j];

var D >= 0.0, <= 1.12;
var Q >= 0.0, <= 2.30;
var d >= 0.0, <= 0.52;
var q >= 0.0, <= 0.187;

############# DEFINED VARIABLES (THEY ARE ELIMINATED BY SUBTITUTION) ###########

var y{j in 1..N} = (alpha*x[j])/(1.0+(alpha-1.0)*x[j]);

var HL{j in 0..N-1} =
  if j == 0 then
    0.1667*exp(-1.087*y[1])
  else 
    0.1667*exp(-1.087*x[j]);

var HV{j in 1..N} =
    0.1349*exp(-3.98*x[j])+0.4397*exp(-0.088*x[j]);

############## EQUATIONS #######################################################

CS_d:
  d = D*y[1];

CS_Q:
  Q = V[1]*(HV[1]-HL[0]);

CS_q:
  q = D*HL[0];

M_eq{j in 1..N-1}:
  sum{k in 1..j}F[k]*z[k] + V[j+1]*y[j+1] = d + (sum{k in 1..j}F[k]+V[j+1]-D)*x[j];

M_tot:
  F[N_F]*z[N_F] = d + (F[N_F]-D)*x[N];

H_eq{j in 1..N-1}:
  sum{k in 1..j} qF[k] + V[j+1]*HV[j+1] = Q + q + (sum{k in 1..j}F[k]+V[j+1]-D)*HL[j];

spec_Lw:
  Lw = (V[1]-D)*(M1*y[1]+M2*(1.0-y[1]));

################### DATA SECTION ###############################################

for {j in 1..N} {
  let   F[j] := 0.0;
  let   z[j] := 0.0;
  let  qF[j] := 0.0;
  let V_L[j] := 2.0;
  let V_U[j] := 4.0;
}

let V_L[N] := V_N;
let V_U[N] := V_N;

let  F[N_F] := 1.0;
let  z[N_F] := 0.5;
let qF[N_F] := F[N_F]*0.1667*exp(-1.087*z[N_F]);

option show_stats 1;
option presolve 10;
option substout 1;
option var_bounds 2;
option nl_comments 0;
write gJacobsen96_3;

# option solver ipopt;

#solve;

# display x;
# display V;
# display D;
# display Q;
# display d;
# display q;
#
# END OF FILE
