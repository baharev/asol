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

// A Double class should be implemented with ctor double(double, double)
#include <cmath>
using namespace std;

double sqr(double x) { return pow(x, 2); }
double  ln(double x) { return log(x); }
//double exp(double ) { return 0; }
//double sqrt(double x) { return ; }

double Inf(double ) { return 0; }
double Sup(double ) { return 0; }
bool in(const double, const double) { return false; }
bool ext_div(const double, const double, const double ) { return false; }
bool Disjoint(double, double) { return 0; }
double Intersect(const double, const double) { return 0; }

#include <iostream>
#include "dag.hpp"
using namespace std;

bool dag_test(const char* file_name) {

	cout << "Running test on file " << file_name << endl;

	bool passed = true;

	dag<double> dd(file_name);

	const int n = dd.number_of_vars();

	double* const r = new double[n];

	dd.get_constraints(r);

	for (int i=0; i<n; ++i) {

		const double val = r[i];

		if (val > 1.0e-7) {
			cout << "Error at index: " << i << ", value " << val << endl;
			passed = false;
		}
	}

	if (passed)
		cout << "Test passed!" << endl;

	return passed;
}
