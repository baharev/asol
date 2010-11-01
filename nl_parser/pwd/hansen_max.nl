g9 0 2 0 6 20100708 0 1 0 1	# problem hansen_max
 8 6 0 0 6	# vars, constraints, objectives, ranges, eqns
 4 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 4 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 18 0	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
C0	#eq_xy
o16	#-
o2	#*
v0	#x
v1	#y
C1	#eq_x2
o16	#-
o5	#^
v0	#x
n2
C2	#eq_y2
o16	#-
o5	#^
v1	#y
n2
C3	#eq_num
n0
C4	#eq_den
n0
C5	#eq_z
o3	#/
v6	#num
v7	#den
x8	# initial guess
0 1.33073
1 1
2 5.1892
3 1.33073
4 1.77084
5 1
6 21.2839
7 4.10157
r	#6 ranges (rhs's)
4 0
4 0
4 0
4 0
4 0
4 0
b	#8 bounds (on variables)
0 1 10
0 1 10
0 5.189197 328.695
0 1 100
0 1 100
0 1 100
0 15.567 1050
0 3 202.4
k7	#intermediate Jacobian column lengths
3
6
7
10
12
14
16
J0 3
0 0
1 0
3 1
J1 2
0 0
4 1
J2 2
1 0
5 1
J3 4
0 -5
3 -14
5 4
6 1
J4 4
1 -1
3 -1
4 -1
7 1
J5 3
2 -1
6 0
7 0
