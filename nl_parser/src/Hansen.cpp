//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2010 Ali Baharev
// All rights reserved. E-mail: <my_first_name.my_last_name@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
//==============================================================================

#include  <iostream>
#include "interval_ext.hpp"
#include "dag.hpp"
using namespace std;
using namespace cxsc;

void Hansen_HC_numerator(interval& x, interval& y, interval& n) {

	interval x_new = (49*n+25)/(49*(14*y+5))+2*y/7-5/interval(49);

	cout << "x_new: " << x_new << endl;

	x = Intersect(x, x_new);

	interval d = 49*x*(x+20/49)-4*n;

	cout << "d: " << d << endl;

	interval x_d = 2*(sqrt(49*n+25)-5)/49;

	cout << "x_d: " << x_d << endl;

	x = Intersect(x, interval(Inf(x_d), Sup(x)));

	interval n_d = x*(49*x+20)/4;

	cout << "n_d" << n_d << endl;

	n = Intersect(n, interval(Inf(n), Sup(n_d)));

	d = Intersect(d, interval(0, 1.0e8));

	interval y_new1 = (7*x-sqrt(d))/4;

	interval y_new2 = (7*x+sqrt(d))/4;

	cout << "y_new1: " << y_new1 << endl;

	cout << "y_new2: " << y_new2 << endl;

	interval y_new = hull(y_new1, y_new2);

	y = Intersect(y, y_new);

}

void example_Hansen_DAG() {

	dag<interval> ia_dag("hansen_max.nl");

	const int n = ia_dag.number_of_vars();

	const int m = ia_dag.number_of_cons();

	const interval vars[] = {
			interval(1, 10),
			interval(1, 10),
			interval(5.189197, 328.625),
			interval(1, 100),
			interval(1, 100),
			interval(1, 100),
			interval(15.567, 1050),
			interval(3, 202.4)
	};

	bool to_delete = false;

	ia_dag.set_variables(vars);

	cout << "z is set to " << ia_dag.get_variables()[2] << endl;

	for (int k=0; k<100; ++k) {

		to_delete = ia_dag.propagate();

		interval* var = ia_dag.get_variables();

		for (int i=0; i<n; ++i) {
			cout << "var[" << i << "] = " << var[i] << endl;
		}

		Hansen_HC_numerator(var[0], var[1], var[6]);

		const interval* con = ia_dag.get_constraints(to_delete);

		for (int i=0; i<m; ++i) {
			cout << "con[" << i << "] = " << con[i] << endl;
		}

	}


	const interval* r = ia_dag.get_constraints(to_delete);

	for (int i=0; i<n; ++i) {

		const interval val = r[i];

		if (AbsMin(val) > 1.0e-7) {
			cout << "Error at index: " << i << ", value " << val << endl;
		}
	}

	return;
}

void example_dummy() {

	dag<interval> ia_dag("dummy.nl");

	const int n = ia_dag.number_of_vars();

	const int m = ia_dag.number_of_cons();

//	const interval vars[] = {
//			interval(-1.0e8, 1.0e8),
//			interval(-1.0e8, 1.0e8),
//	};

	const interval vars[] = {
			interval(-1.0,-0.5),
			interval( 0.5, 1.0),
	};

	bool to_delete = false;

	ia_dag.set_variables(vars);

	for (int k=0; k<100; ++k) {

		to_delete = ia_dag.propagate();

		const interval* ia_var = ia_dag.get_variables();

		for (int i=0; i<n; ++i) {
			cout << "var[" << i << "] = " << ia_var[i] << endl;
		}

		const interval* con = ia_dag.get_constraints(to_delete);

		for (int i=0; i<m; ++i) {
			cout << "con[" << i << "] = " << con[i] << endl;
		}

	}


	const interval* r = ia_dag.get_constraints(to_delete);

	for (int i=0; i<n; ++i) {

		const interval val = r[i];

		if (AbsMin(val) > 1.0e-7) {
			cout << "Error at index: " << i << ", value " << val << endl;
		}
	}

	return;
}

int main() {

	example_Hansen_DAG();

	return 0;
}

