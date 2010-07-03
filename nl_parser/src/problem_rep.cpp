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

#include <sstream>
#include <cstdlib>
#include "problem_rep.hpp"
#include "utility.hpp"

using namespace std;

problem_rep::problem_rep(const string& name, const problem_statistics& statistics) {

	stats = statistics;
	problem_name = name;

	defined_var_index.reserve(stats.n_defined_variables());
	constraint_index.reserve(stats.n_equations());

	lower_bound.reserve(stats.n_variables());
	upper_bound.reserve(stats.n_variables());
	initial_point.reserve(stats.n_variables());
}

const string problem_rep::print(const char letter, const int index) {

	ostringstream s;
	s << letter << index << flush;
	return s.str();
}

const string problem_rep::add_primitive(const primitive_rep& prim) {

	const int index = n_primitives();

	primitives.push_back(prim);

	return print('t', index);
}

const int problem_rep::get_index(const string& arg) {

	const string val = arg.substr(1);

	if (arg.at(0) != 'n')
		return atoi( val.c_str() );

	const int index = n_numeric_constants();

	numeric_constants.push_back(val);

	return index;
}

const string problem_rep::add_unary_primitive(const string& op, const string& arg) {

	const int index = get_index(arg);
	const arg_t arg1(arg.at(0), index);

	return add_primitive( primitive_rep(op.c_str(), 1, arg1) );
}

const string problem_rep::add_binary_primitive(const string& op, const string& arg1, const string& arg2) {

	const int idx1 = get_index(arg1);
	const arg_t arg_1(arg1.at(0), idx1);

	const int idx2 = get_index(arg2);
	const arg_t arg_2(arg2.at(0), idx2);

	return add_primitive( primitive_rep(op.c_str(), 2, arg_1, arg_2) );
}

void problem_rep::add_defined_variable(const int defined_variable_index, const int primitive_index) {

	const int num_of_defined_vars = defined_variable_index - stats.n_variables();

	if ( num_of_defined_vars != n_defined_vars() )
		error("inconsistent number of defined variables");

	defined_var_index.push_back(primitive_index);
}

void problem_rep::add_constraint(const int index_of_constraint, const int primitive_index) {

	if ( index_of_constraint != n_constraints() )
		error("inconsistent number of constraints");

	constraint_index.push_back(primitive_index);
}

void problem_rep::push_back_variable_bounds(const string& lb, const string& ub) {

	const double l = atof(lb.c_str());
	const double u = atof(ub.c_str());

	if (l > u)
		error("inconsistent variable bounds");

	if (n_bounds() >= stats.n_variables())
		error("too many variable bounds");

	if (lower_bound.size() != upper_bound.size())
		error("inconsistent state of bound containers");

	lower_bound.push_back(lb);
	upper_bound.push_back(ub);
}

void problem_rep::push_back_initial_point(const int index, const string& x0) {

	const int size = n_initial_points();

	if (size >= stats.n_variables())
		error("too many initial points");

	if (size != index)
		error("unexpected variable index in x segment");

	initial_point.push_back(x0);
}

void problem_rep::check_initial_point_consistency() {

	const int n = n_bounds();

	if (n_bounds() != n_initial_points())
		error("conflicting number of bounds and initial points");

	for (int i=0; i<n; ++i) {
		const double lb = atof(lower_bound.at(i).c_str());
		const double ub = atof(upper_bound.at(i).c_str());
		const double x0 = atof(initial_point.at(i).c_str());

		if ((x0 < lb) || (x0 > ub)) {
			warning("initial point is out of bounds, resetting to midpoint");
			initial_point.at(i) = double_to_string( (lb+ub)/2.0 );
		}
	}
}

void problem_rep::check_consistency() {

	if (n_defined_vars() != stats.n_defined_variables())
		error("inconsistent number of defined variables");

	if (n_constraints() != stats.n_equations())
		error("inconsistent number of equations");

	if (n_initial_points() != 0) {
		check_initial_point_consistency();
	}
}
