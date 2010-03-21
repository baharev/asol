//=============================================================================
//
//  Examples showing the modules of the nonlinear system solver ASOL
//
//  Copyright (C) 2010  Ali Baharev
//  All rights reserved. E-mail: <my_first_name.my_last_name@gmail.com>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//=============================================================================
//
// Version: 0.01
// Last updated on: 20 Mar 2010
//

#include <iostream>
#include <assert.h>

#define HAVE_CXSC
#ifdef  HAVE_CXSC
#include "interval_ext.hpp"
using namespace cxsc;
#endif

#include "affine.hpp"
#include "dag.hpp"

#define HAVE_GLPK
#ifdef  HAVE_GLPK
#include "lp_pruning.hpp"
#endif

using namespace std;
using namespace aa_lib;
using namespace dag_builder;

int Main() {

#ifdef HAVE_CXSC
	{
		cout << endl << "===================================================" << endl;
		cout << endl << "EXAMPLE 1" << endl << endl;
		cout << "Given the following equations:" << endl;
		cout << endl;
		cout << "   x^2 + y^2 = 1" << endl;
		cout << "   x^2 = y" << endl;
		cout << endl;
		cout << "we are trying to tighten the bounds on x and y with constraint propagation." << endl;
		cout << "The equations are taken from:" << endl;
		cout << "F. Goualard, Interval Multivalued Inverse Functions, SCAN'08" << endl;

		cout << endl << "Reading input example.nl" << endl;

		// Build the interval DAG
		dag<interval> ia_dag("example.nl");

		const int n = ia_dag.number_of_vars();

		const interval* const var = ia_dag.get_variables();

		cout << endl << "Input intervals (big M-like):" << endl << endl;
		for (int i=0; i<n; ++i) {
			cout << "   " << var[i] << endl;
		}

		cout << endl << "Performing propagation on the interval DAG" << endl;
		ia_dag.propagate();

		cout << endl << "Lower and upper bounds after propagation:" << endl << endl;
		for (int i=0; i<n; ++i) {
			cout << "   " << var[i] << endl;
		}

		cout << endl << "===================================================" << endl;
		cout << endl << "EXAMPLE 2" << endl;
		cout << endl << "Overwriting these variable bounds with:" << endl;
		cout << endl << "-1.0 <= x <= -0.5 and 0.5 <= y <= 1.0" << endl;
		
		interval ia_var[2];

		ia_var[0] = interval(-1.0, -0.5);
		ia_var[1] = interval( 0.5,  1.0);

		ia_dag.set_variables(ia_var);

		cout << endl << "Then we apply propagation 100 times" << endl;
		for (int i=1; i<=100; ++i)
			ia_dag.propagate();

		cout << endl << "Lower and upper bounds after propagation:" << endl << endl;
		for (int i=0; i<n; ++i) {
			cout << "   " << var[i] << endl;
		}
	}
#endif

	cout << endl << "===================================================" << endl;
	cout << endl << "EXAMPLE 3" << endl << endl;
	cout << "We are linearizing the expressions on the left hand side with affine arithmetic" << endl;
	cout << endl;
	cout << "   x^2 + y^2 - 1 = 0" << endl;
	cout << "   x^2 - y       = 0" << endl;
	cout << endl;
	cout << "over the following box:" << endl << endl;
	cout << "-1.0 <= x <= -0.5 and 0.5 <= y <= 1.0" << endl;

	// Build the affine DAG
	dag<affine> aa_dag("example.nl");

	const int n = aa_dag.number_of_vars();

	// Noise variables do not implement reference counting
	// there is no need for that in case of solving equations
	affine::set_max_used_index(0);

	// Overwriting the original [-1.0e8, 1.0e8] variable bounds
	// with -1.0 <= x <= 0.0 and 0.0 <= y <= 1.0
	affine* const aa_var = new affine[n];

	aa_var[0] = affine(-1.0, -0.5);
	aa_var[1] = affine( 0.5,  1.0);

	aa_dag.set_variables(aa_var);

	cout << endl;
	cout << "Affine forms store index-value pairs, IA shows the enclosed range" << endl;
	cout << endl << "Input affine form for x:" << endl;
	
	cout << endl << aa_var[0] << flush;
	cout << "It means that x = -0.75 + 0.25*e1 where -1 <= e1 <= 1" << endl;
	cout << endl;
	cout << "Similarly     y =  0.75 + 0.25*e2 where -1 <= e2 <= 1" << endl;
	cout << endl << "Input affine form for y:" << endl;
	cout << endl << aa_var[1] << flush;

	cout << "The affine linearization yields the following linear constraints" << endl;
	cout << "enclosing all solutions to the original system of nonlinear equations: " << endl;
	cout << endl;
	cout << " 0.18750 - 0.375*e1 + 0.375*e2 + 0.03125*e3 + 0.03125*e4              = 0.0" << endl;
	cout << "-0.15625 - 0.375*e1 - 0.250*e2                           + 0.03125*e5 = 0.0" << endl;
	cout << endl;
	cout << "                  -1 <= e1, e2, e3, e4, e5 <= 1" << endl;

	// Get the linearized constraints
	affine* const linearized_constraint = new affine[n];

	aa_dag.get_constraints(linearized_constraint);

	cout << endl << "The affine linearization using the affine DAG:" << endl << endl;
	for (int i=0; i<n; ++i) {

		cout << "Eq #" << i << endl;
		cout << linearized_constraint[i] << flush;
	}

	cout << "If we solved the following 4 LP problems: min e1, max e1, min e2, max e2" << endl;
	cout << "subject to the linear constraints and variable bounds obtained with" << endl;
	cout << "affine arithmetic above, we would conclude that the variable bounds " << endl;
	cout << "can be tightened:" << endl;
	cout << endl;
	cout << "   -0.791667 <= x <= -0.733333" << endl;
	cout << "    0.575    <= y <=  0.65" << endl;
	cout << endl;
	cout << "without losing any solutions." << endl;

#ifdef HAVE_GLPK
	{
		// Actually, here is how the above LP problems can be solved;
		// this requires that GLPK is installed

		LP_op lp_operator(aa_dag.number_of_vars(), aa_dag.number_of_nzeros());

		bool toDelete = false;

		lp_operator.prune(linearized_constraint, aa_var, toDelete);

		cout << endl << "Variable bounds after LP pruning:" << endl;

		for (int i=0; i<n; ++i)
			cout << "[ " << aa_var[i].inf() << ", " << aa_var[i].sup() << "]" << endl;
	}
#endif

	cout << endl << "===================================================" << endl;
	cout << endl << "EXAMPLE 4" << endl;

	cout << endl << "The affine linearization is now performed in C++," << flush;
	cout << endl << "it has to be the same as above:" << endl << endl;

	// reset the noise variable counter to zero
	affine::set_max_used_index(0);

	affine x = affine(-1.0, -0.5);
	affine y = affine( 0.5,  1.0);

	linearized_constraint[0] = sqr(x) + sqr(y) - affine(1.0);
	linearized_constraint[1] = sqr(x) - y;

	for (int i=0; i<n; ++i) {

		cout << "Eq #" << i << endl;
		cout << linearized_constraint[i] << flush;
	}

	cout << endl << "===================================================" << endl;
	cout << endl << "EXAMPLE 5" << endl;
	cout << endl;
	cout << "The C++ code below imitates what the constructor of the dag performs," << endl;
	cout << "the expression for this simple example is x^2+y^2-1 and" << endl;
	cout << "   -1.0 <= x <= -0.5 and 0.5 <= y <= 1.0" << endl << endl;

	// reset the noise variable counter to zero
	affine::set_max_used_index(0);

	// Once the .nl file is processed, a similar procedure is performed
	// in the constructor of the dag
	//
	// private members of dag
	//
	// var: variables
	// num: numerical constants
	// tmp: temporary values
	//
	// Building the DAG:
	affine var[] = { affine(-1.0, -0.5), affine( 0.5,  1.0) };

	affine num[] = { affine(2.0), affine(-1.0) };
	
	affine tmp[4];

	node<affine>** Node = new node<affine>* [4];

	// var[0]^2 + var[1]^2 - 1
	//
	// tmp[0] = var[0]^2
	Node[0] = new PowInt<affine>(tmp  , var  , num  );

	// tmp[1] = var[1]^2
	Node[1] = new PowInt<affine>(tmp+1, var+1, num  );

	// tmp[2] = tmp[0] + tmp[1]
	Node[2] = new Plus<affine>  (tmp+2, tmp  , tmp+1);

	// tmp[3] = tmp[2] - 1
	Node[3] = new Plus<affine>  (tmp+3, tmp+2, num+1);

	// The trick is to get the offsets right, for example 3, 2, 1 in the call
	// just above. The major part of the .nl file processing is about calculating
	// these offsets.

	// Now the DAG of the equation is built, we are ready to evaluate
	// the equation with respect to the variable bounds

	cout << "It also imitates the evaluate_all() private member function" << endl;
	cout << endl;

	// The dag does not need to know the type of its nodes
	for (int i=0; i<4; ++i)
		Node[i]->evaluate();

	cout << "The result should be the same as Eq #0 before:" << endl << endl;
	cout << Node[3]->value() << endl;

	// As in the destructor
	for (int i=0; i<4; ++i)
		delete Node[i];
	delete[] Node;

	// End of the examples
	//================================================================

	delete[] aa_var;
	delete[] linearized_constraint;

	return 0;
}


int main() {

	cout << "Examples showing the modules of the nonlinear system solver ASOL" << endl;
	cout << endl;
	cout << "Copyright (C) 2010  Ali Baharev" << endl;
	cout << "All rights reserved. E-mail: <my_first_name.my_last_name@gmail.com>" << endl;
	cout << endl;
	cout << "This program has ABSOLUTELY NO WARRANTY." << endl;
	cout << endl;
	cout << "This program is free software; you may re-distribute it under the terms" << endl;
	cout << "of the GNU General Public License version 3 or later." << endl;
	cout << endl;
	cout << "Built on " << __DATE__ << " " << __TIME__ << endl << endl;

	Main();

	cerr << endl << "Please press enter to quit..." << endl;
	cin.get();

	return 0;
}
