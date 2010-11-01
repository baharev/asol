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

int main() {

	dag<interval> ia_dag("hansen_max.nl");

	const int n = ia_dag.number_of_vars();

//			interval(5.189197, 328.695),

	const interval vars[] = {
			interval(1, 10),
			interval(1, 10),
			interval(15.567, 1050),
			interval(3, 202.4),
			interval(5.189197, 14),
			interval(1, 100),
			interval(1, 100),
			interval(1, 100)
	};

	bool to_delete = false;

	ia_dag.set_variables(vars);

	for (int i=0; i<100; ++i) {
		to_delete = ia_dag.propagate();
		cout << "x: " << ia_dag.get_variables()[0] << endl;
		cout << "y: " << ia_dag.get_variables()[1] << endl;
		cout << "z: " << ia_dag.get_variables()[3] << endl;
	}


	const interval* r = ia_dag.get_constraints(to_delete);

	for (int i=0; i<n; ++i) {

		const interval val = r[i];

		if (AbsMin(val) > 1.0e-7) {
			cout << "Error at index: " << i << ", value " << val << endl;
		}
	}

}

