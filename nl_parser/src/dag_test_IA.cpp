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

#include <iostream>
#include "interval_ext.hpp"
#include "dag.hpp"
using namespace std;
using namespace cxsc;

bool dag_test(const char* file_name) {

	cout << "Running test on file " << file_name << endl;

	bool passed = true;

	dag<interval> dd(file_name);

	const int n = dd.number_of_vars();

	bool to_delete = false;

	const interval* r = dd.get_constraints(to_delete);

	for (int i=0; i<n; ++i) {

		const interval val = r[i];

		if (AbsMin(val) > 1.0e-7) {
			cout << "Error at index: " << i << ", value " << val << endl;
			passed = false;
		}
	}

	if (passed)
		cout << "Test passed!" << endl;

	return passed;
}
