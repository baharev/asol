g9 0 2 0 19 20090130 0 4 0 241	# problem jacobsen_comments
 19 19 0 0 19	# vars, constraints, objectives, ranges, eqns
 19 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 16 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 90 0	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 2 0 22 0	# common exprs: b,c,o,c1,o1
V19 0 0	#y[1]
o3	#/
o2	#*
n3.55
v0	#x[1]
o0	# + 
n1
o2	#*
n2.55
v0	#x[1]
V20 0 0	#HL[0]
o2	#*
n0.1667
o44	#exp
o2	#*
n-1.087
v19	#y[1]
C0	#CS_d
o16	#-
o2	#*
v15	#D
v19	#y[1]
V21 0 2	#HV[1]
o0	# + 
o2	#*
n0.1349
o44	#exp
o2	#*
n-3.98
v0	#x[1]
o2	#*
n0.4397
o44	#exp
o2	#*
n-0.088
v0	#x[1]
C1	#CS_Q
o16	#-
o2	#*
v8	#V[1]
o1	# - 
v21	#HV[1]
v20	#HL[0]
C2	#CS_q
o16	#-
o2	#*
v15	#D
v20	#HL[0]
V22 0 4	#y[2]
o3	#/
o2	#*
n3.55
v1	#x[2]
o0	# + 
n1
o2	#*
n2.55
v1	#x[2]
C3	#M_eq[1]
o1	# - 
o2	#*
v9	#V[2]
v22	#y[2]
o2	#*
o1	# - 
v9	#V[2]
v15	#D
v0	#x[1]
V23 0 5	#y[3]
o3	#/
o2	#*
n3.55
v2	#x[3]
o0	# + 
n1
o2	#*
n2.55
v2	#x[3]
C4	#M_eq[2]
o1	# - 
o2	#*
v10	#V[3]
v23	#y[3]
o2	#*
o1	# - 
v10	#V[3]
v15	#D
v1	#x[2]
V24 0 6	#y[4]
o3	#/
o2	#*
n3.55
v3	#x[4]
o0	# + 
n1
o2	#*
n2.55
v3	#x[4]
C5	#M_eq[3]
o1	# - 
o2	#*
v11	#V[4]
v24	#y[4]
o2	#*
o1	# - 
v11	#V[4]
v15	#D
v2	#x[3]
V25 0 7	#y[5]
o3	#/
o2	#*
n3.55
v4	#x[5]
o0	# + 
n1
o2	#*
n2.55
v4	#x[5]
C6	#M_eq[4]
o1	# - 
o2	#*
v12	#V[5]
v25	#y[5]
o2	#*
o1	# - 
v12	#V[5]
v15	#D
v3	#x[4]
V26 0 8	#y[6]
o3	#/
o2	#*
n3.55
v5	#x[6]
o0	# + 
n1
o2	#*
n2.55
v5	#x[6]
C7	#M_eq[5]
o1	# - 
o2	#*
v13	#V[6]
v26	#y[6]
o2	#*
o1	# - 
o0	# + 
n1
v13	#V[6]
v15	#D
v4	#x[5]
V27 0 9	#y[7]
o3	#/
o2	#*
n3.55
v6	#x[7]
o0	# + 
n1
o2	#*
n2.55
v6	#x[7]
C8	#M_eq[6]
o1	# - 
o2	#*
v14	#V[7]
v27	#y[7]
o2	#*
o1	# - 
o0	# + 
n1
v14	#V[7]
v15	#D
v5	#x[6]
V28 0 10	#y[8]
o3	#/
o2	#*
n3.55
v7	#x[8]
o0	# + 
n1
o2	#*
n2.55
v7	#x[8]
C9	#M_eq[7]
o1	# - 
o2	#*
n2
v28	#y[8]
o2	#*
o1	# - 
n3
v15	#D
v6	#x[7]
C10	#M_tot
o16	#-
o2	#*
o1	# - 
n1
v15	#D
v7	#x[8]
V29 0 12	#HL[1]
o2	#*
n0.1667
o44	#exp
o2	#*
n-1.087
v0	#x[1]
V30 0 12	#HV[2]
o0	# + 
o2	#*
n0.1349
o44	#exp
o2	#*
n-3.98
v1	#x[2]
o2	#*
n0.4397
o44	#exp
o2	#*
n-0.088
v1	#x[2]
C11	#H_eq[1]
o1	# - 
o2	#*
v9	#V[2]
v30	#HV[2]
o2	#*
o1	# - 
v9	#V[2]
v15	#D
v29	#HL[1]
V31 0 13	#HL[2]
o2	#*
n0.1667
o44	#exp
o2	#*
n-1.087
v1	#x[2]
V32 0 13	#HV[3]
o0	# + 
o2	#*
n0.1349
o44	#exp
o2	#*
n-3.98
v2	#x[3]
o2	#*
n0.4397
o44	#exp
o2	#*
n-0.088
v2	#x[3]
C12	#H_eq[2]
o1	# - 
o2	#*
v10	#V[3]
v32	#HV[3]
o2	#*
o1	# - 
v10	#V[3]
v15	#D
v31	#HL[2]
V33 0 14	#HL[3]
o2	#*
n0.1667
o44	#exp
o2	#*
n-1.087
v2	#x[3]
V34 0 14	#HV[4]
o0	# + 
o2	#*
n0.1349
o44	#exp
o2	#*
n-3.98
v3	#x[4]
o2	#*
n0.4397
o44	#exp
o2	#*
n-0.088
v3	#x[4]
C13	#H_eq[3]
o1	# - 
o2	#*
v11	#V[4]
v34	#HV[4]
o2	#*
o1	# - 
v11	#V[4]
v15	#D
v33	#HL[3]
V35 0 15	#HL[4]
o2	#*
n0.1667
o44	#exp
o2	#*
n-1.087
v3	#x[4]
V36 0 15	#HV[5]
o0	# + 
o2	#*
n0.1349
o44	#exp
o2	#*
n-3.98
v4	#x[5]
o2	#*
n0.4397
o44	#exp
o2	#*
n-0.088
v4	#x[5]
C14	#H_eq[4]
o1	# - 
o2	#*
v12	#V[5]
v36	#HV[5]
o2	#*
o1	# - 
v12	#V[5]
v15	#D
v35	#HL[4]
V37 0 16	#HL[5]
o2	#*
n0.1667
o44	#exp
o2	#*
n-1.087
v4	#x[5]
V38 0 16	#HV[6]
o0	# + 
o2	#*
n0.1349
o44	#exp
o2	#*
n-3.98
v5	#x[6]
o2	#*
n0.4397
o44	#exp
o2	#*
n-0.088
v5	#x[6]
C15	#H_eq[5]
o1	# - 
o2	#*
v13	#V[6]
v38	#HV[6]
o2	#*
o1	# - 
o0	# + 
n1
v13	#V[6]
v15	#D
v37	#HL[5]
V39 0 17	#HL[6]
o2	#*
n0.1667
o44	#exp
o2	#*
n-1.087
v5	#x[6]
V40 0 17	#HV[7]
o0	# + 
o2	#*
n0.1349
o44	#exp
o2	#*
n-3.98
v6	#x[7]
o2	#*
n0.4397
o44	#exp
o2	#*
n-0.088
v6	#x[7]
C16	#H_eq[6]
o1	# - 
o2	#*
v14	#V[7]
v40	#HV[7]
o2	#*
o1	# - 
o0	# + 
n1
v14	#V[7]
v15	#D
v39	#HL[6]
V41 0 18	#HL[7]
o2	#*
n0.1667
o44	#exp
o2	#*
n-1.087
v6	#x[7]
V42 0 18	#HV[8]
o0	# + 
o2	#*
n0.1349
o44	#exp
o2	#*
n-3.98
v7	#x[8]
o2	#*
n0.4397
o44	#exp
o2	#*
n-0.088
v7	#x[8]
C17	#H_eq[7]
o1	# - 
o2	#*
n2
v42	#HV[8]
o2	#*
o1	# - 
n3
v15	#D
v41	#HL[7]
C18	#spec_Lw
o16	#-
o2	#*
o1	# - 
v8	#V[1]
v15	#D
o0	# + 
o2	#*
n32.04
v19	#y[1]
o2	#*
n60.1
o1	# - 
n1
v19	#y[1]
r	#19 ranges (rhs's)
4 0
4 0
4 0
4 0
4 0
4 0
4 0
4 -0.5
4 -0.5
4 -0.5
4 -0.5
4 0
4 0
4 0
4 0
4 -0.09680472351714639
4 -0.09680472351714639
4 -0.09680472351714639
4 -59
b	#19 bounds (on variables)
0 0.0001 1
0 0.0001 1
0 0.0001 1
0 0.0001 1
0 0.0001 1
0 0.0001 1
0 0.0001 1
0 0.0001 1
0 1.5 2.5
0 1.5 2.5
0 1.5 2.5
0 1.5 2.5
0 1.5 2.5
0 1.5 2.5
0 1.5 2.5
0 0 1
0 0 2.3
0 0 0.5
0 0 0.187
k18	#intermediate Jacobian column lengths
6
10
14
18
22
26
30
33
35
37
39
41
43
45
47
65
73
82
J0 3
0 0
15 0
17 1
J1 3
0 0
8 0
16 1
J2 3
0 0
15 0
18 1
J3 5
0 0
1 0
9 0
15 0
17 -1
J4 5
1 0
2 0
10 0
15 0
17 -1
J5 5
2 0
3 0
11 0
15 0
17 -1
J6 5
3 0
4 0
12 0
15 0
17 -1
J7 5
4 0
5 0
13 0
15 0
17 -1
J8 5
5 0
6 0
14 0
15 0
17 -1
J9 4
6 0
7 0
15 0
17 -1
J10 3
7 0
15 0
17 -1
J11 6
0 0
1 0
9 0
15 0
16 -1
18 -1
J12 6
1 0
2 0
10 0
15 0
16 -1
18 -1
J13 6
2 0
3 0
11 0
15 0
16 -1
18 -1
J14 6
3 0
4 0
12 0
15 0
16 -1
18 -1
J15 6
4 0
5 0
13 0
15 0
16 -1
18 -1
J16 6
5 0
6 0
14 0
15 0
16 -1
18 -1
J17 5
6 0
7 0
15 0
16 -1
18 -1
J18 3
0 0
8 0
15 0
