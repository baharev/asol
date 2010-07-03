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
double sqr(double ) { return 0; }
double  ln(double ) { return 0; }
double exp(double ) { return 0; }
double sqrt(double ) { return 0; }
double Inf(double ) { return 0; }
double Sup(double ) { return 0; }
bool ext_div(const double, const double, const double ) { return false; }
bool Disjoint(double, double) { return 0; }
double Intersect(const double, const double) { return 0; }

#include "dag.hpp"

bool dag_test(const char* file_name) {

	dag<double> dd(file_name);

	return false;
}
