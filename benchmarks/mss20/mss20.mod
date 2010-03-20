#################### DESCRIPTION ###############################################
#
# Benchmark mss20
#
# Last updated: 20 Mar 2010
#
# The model has three solutions corresponding to the lower concentration steady
# state (LSS), upper stable steady state or higher steady state (HSS), and
# unstable steady state (USS) of a distillation column with specified methanol -
# methyl butyrate - toluene feed composition. The distillate molar flow rate,
# denoted by D, can be chosen as bifurcation parameter.
#
# Conventional inside-out and simultaneous correction procedures were reported
# to miss the unstable steady state solution. Previously, all steady state
# branches were computed with an appropriate continuation method. A solution
# path was traced with respect to the bifurcation parameter. The initial
# estimates were carefully chosen, special attention was paid to the turning
# points and branch switching. See:
#
# Vadapalli A, Seader JD;
# A generalized framework for computing bifurcation diagrams using process
# simulation programs
# Computers and Chemical Engineering 2001, 25, 445--464
#
# Kannan A, Joshi MR, Reddy GR, Shah DM;
## Multiple-Steady-States Identification in Homogeneous Azeotropic Distillation
# Using a Process Simulator
# Ind Eng Chem Res. 2005, 44, 4386--4399
#
# Güttinger TE, Dorn C, Morari M;
# Experimental Study of Multiple Steady States in Homogeneous Azeotropic
# Distillation
# Ind Eng Chem Res. 1997, 36, 794--802
#
# Doedel EJ, Wang XJ, Fairgrieve TF;
# AUTO94: Software for Continuation and Bifurcation Problems in Ordinary
# Differential Equations
# Technical Report CRPC#95#1, Center for Research on Parallel Computing,
# California Institute of Technology, Pasadena CA 91125, 1995
#
# Except the number of stages and the feed stage location, the model and its
# parameters correspond to the "Auto model" of Güttinger, Dorn, Morari.
#
# The method proposed in the manuscript cited below seems to be the first
# alternative to continuation methods that is capable to find all solutions to
# the this problem automatically, without any assistance from the user. See:
#
# A. Baharev, L. Kolev, E. Rév;
# Computing multiple steady states in homogeneous azeotropic and ideal
# two-product distillation;
# (preprint available from http://reliablecomputing.eu)
#
################### PROBLEM STATISTICS #########################################
#
# Presolve eliminates 60 constraints and 60 variables.
# Substitution eliminates 380 variables.
# Adjusted problem:
# 140 variables, all nonlinear
# 140 constraints; 594 nonzeros
#	  120 nonlinear constraints
#	  20 linear constraints
# 0 objectives.
#
################# PARAMETERS FROM SPECIFICATION ################################

# NUMBER OF STAGES
param N := 20;

# FEED STAGE LOCATION
param N_F := 10;

# NUMBER OF COMPONENTS
param C := 3;

# DISTILLATE MOLAR FLOW RATE
param D;

let D := 0.455;

# VAPOR FLOW RATE, FROM SPECIFICATION
param V := 1.38;

# FEED LOCATION AND VALUES ARE GIVEN IN THE DATA SECTION

param F{j in 0..N};
param f{i in 1..C, j in 0..N};

# AUXILIARY PARAMETERS

param L{j in 0..N} = V - D + sum{k in 0..j} F[k];
param B := F[N_F] - D;

################## PARAMETERS ##################################################

param a{1..C};
param b{1..C};
param c{1..C};

param r{1..C, 1..C};
param s{1..C, 1..C};

param P := 100000.0;

# LOWER/UPPER BOUNDS ON THE VARIABLES, INITIAL ESTIMATES (IF NEEDED)
# VALUES ARE GIVEN IN THE DATA SECTION

param x_L{1..C, 1..N}; param x_U{1..C, 1..N}; param x_0{1..C, 1..N};
param K_L{1..C, 1..N}; param K_U{1..C, 1..N}; param K_0{1..C, 1..N};

param T_0{1..N};

############### VARIABLES ######################################################

var x{i in 1..C, j in 1..N} >= x_L[i,j], <= x_U[i,j], := x_0[i,j];

var K{i in 1..C, j in 1..N} >= K_L[i,j], <= K_U[i,j], := K_0[i,j];

var T{j in 1..N} >= 336.3, <= 383.4, := T_0[j];

####### DEFINED VARIABLES (THEY ARE ELIMINATED BY PRESOLVE / SUBTITUTION) ######

var p{i in 1..C, j in 1..N} = exp(a[i]+b[i]/(T[j]+c[i]));

var rcp_T{j in 1..N} = 1.0/T[j];

var Lambda{i1 in 1..C, i2 in 1..C, j in 1..N} = exp(r[i1,i2]+s[i1,i2]*rcp_T[j]);

var sum_xLambda{i in 1..C, j in 1..N} = sum{i1 in 1..C} (x[i1,j]*Lambda[i,i1,j]);

var rcp_sum_xLambda{i in 1..C, j in 1..N} = 1.0/sum_xLambda[i,j];

var gamma{i in 1..C, j in 1..N} =
  exp( -log(sum_xLambda[i,j]) + 1.0 - (sum{i2 in 1..C} (x[i2,j]*Lambda[i2,i,j]*rcp_sum_xLambda[i2,j])) );
  
############## EQUATIONS #######################################################

# AUXILIARY EQUATIONS

E_aux_K{j in 1..N, i in 1..C}:
	K[i,j] - gamma[i,j]*(p[i,j]/P) = 0.0;

# MATERIAL BALANCES

M_tot{i in 1..C}:
	D*(K[i,1]*x[i,1]) + B*x[i,N] - f[i,N_F] = 0.0;

# NOTE THE UNUSUAL FORMULATION
M_eq{j in 1..N-1, i in 1..C}:
	L[j]*x[i,j] + sum{i1 in j+1..N} f[i,i1] - B*x[i,N] - V*(K[i,j+1]*x[i,j+1]) = 0.0;

# SUMMATION EQUATIONS

S_x_eq{j in 1..N}:
	sum{i in 1..C} x[i,j] - 1.0 = 0.0;

################### DATA SECTION ###############################################

data;

let a[1] := 23.4832;
let a[2] := 20.5110;
let a[3] := 20.9064;

let b[1] := -3634.01;
let b[2] := -2664.30;
let b[3] := -3096.52;

let c[1] := -33.768;
let c[2] := -79.483;
let c[3] := -53.668;

let r[1,2] :=  0.7411;
let r[1,3] :=  0.9645;
let r[2,3] := -1.4350;

let r[2,1] := -1.0250;
let r[3,1] := -0.9645;
let r[3,2] :=  2.7470;

let r[1,1] := 0.0;
let r[2,2] := 0.0;
let r[3,3] := 0.0;

let s[1,2] := -477.00;
let s[1,3] := -903.1024;
let s[2,3] :=  768.20;

let s[2,1] :=  72.78;
let s[3,1] := -140.9995;
let s[3,2] := -1419.0;

let s[1,1] := 0.0;
let s[2,2] := 0.0;
let s[3,3] := 0.0;

# LOWER AND UPPER BOUNDS ON THE VARIABLES

for {j in 1..N} {
	let x_L[1,j] := 0.0;
	let x_U[1,j] := 0.9998;
	let x_L[2,j] := 1.0e-4;
	let x_U[2,j] := 0.9999;
	let x_L[3,j] := 1.0e-4;
	let x_U[3,j] := 0.9999;
}

# THIS BOUND SEEMS TO BE REASONABLE FROM ENGINEENERING POINT OF VIEW

let x_L[1,1] := 0.83;

# THESE BOUNDS WERE CALCULATED IN ADVANCE

for {j in 1..N} {
	let K_L[1,j] := 0.9784;
	let K_L[2,j] := 0.2445;
	let K_L[3,j] := 0.2745;
	
	let K_U[1,j] := 40.52;
	let K_U[2,j] := 1.317;
	let K_U[3,j] := 1.975;
}

# COMES FROM x_L[1,1] (K=y/x, y<1) 

let K_U[1,1] := 1.21;

# FEED VALUES, LOCATION

for {i in 1..C, j in 0..N}
	let f[i,j] := 0.0;

let f[1, N_F] := 0.4098370;
let f[2, N_F] := 0.01229769;
let f[3, N_F] := 0.06090665;

for {j in 0..N}
	let F[j] := sum{i in 1..C} f[i,j];

################################################################################

# DUMB INITIAL ESTIMATES (IF NEEDED)

for {i in 1..C, j in 1..N}
	let x_0[i,j] := 0.33;

for {j in 1..N}
	let T_0[j] := 337.0;

for {i in 1..C, j in 1..N} # FOR THE LOWER BRANCH
	let K_0[i,j] := 1.0;

# UNCOMMENT ONE OF THESE TO SWITCH TO ANOTHER BRANCH
	
# for {j in 1..N} {  # FOR THE BRANCH IN THE MIDDLE
	# let K_0[1,j] := 40.52;
# }

# for {j in N_F..N} {# FOR THE UPPER BRANCH
	# let K_0[1,j] := 40.52;
# }

################################################################################

option show_stats 1;
option presolve 10;
option substout 1;
option var_bounds 2;
option nl_comments 0;
option nl_permute 0;
#option auxfiles c;
#write gmss20;

# option solver ipopt;
# solve;
#
# END OF FILE
