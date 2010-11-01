g9 0 2 0 2 20100708 0 1 0 1	# problem dummy
 2 2 0 0 2	# vars, constraints, objectives, ranges, eqns
 2 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 2 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 4 0	# nonzeros in Jacobian, gradients
 0 1	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
C0	#eq1
o0	# + 
o5	#^
v0	#x
n2
o5	#^
v1	#y
n2
C1	#eq2
o5	#^
v0	#x
n2
x2	# initial guess
0 1
1 1
r	#2 ranges (rhs's)
4 1
4 0
b	#2 bounds (on variables)
0 -1e+08 1e+08
0 -1e+08 1e+08
k1	#intermediate Jacobian column lengths
2
J0 2
0 0
1 0
J1 2
0 0
1 -1
