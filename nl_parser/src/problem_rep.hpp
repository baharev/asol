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

#ifndef PROBLEM_REP_HPP_
#define PROBLEM_REP_HPP_

#include <vector>
#include <list>
#include <string>
#include "problem_statistics.hpp"
#include "primitive_rep.hpp"

using std::vector;
using std::list;
using std::string;

class problem_rep {

public:

	problem_rep(const string& name, const problem_statistics& statistics);

	problem_rep() { }

	const string add_unary_primitive(const string& op, const string& arg);

	const string add_binary_primitive(const string& op, const string& arg1, const string& arg2);

	void add_defined_variable(const int defined_variable_index, const int primitive_index);

	void add_constraint(const int index_of_constraint, const int primitive_index);

	void push_back_variable_bounds(const string& lb, const string& ub);

	void push_back_initial_point(const int index, const string& initial_point);

	void check_consistency();

	const string name() const { return problem_name; }

	int n_primitives()        const { return static_cast<int> (        primitives.size() ); }

	problem_statistics stats;

	friend class interpreter_C;
	friend class interpreter_DAG;

private:

	int n_numeric_constants() const { return static_cast<int> ( numeric_constants.size() ); }
	int n_defined_vars()      const { return static_cast<int> ( defined_var_index.size() ); }
	int n_constraints()       const { return static_cast<int> (  constraint_index.size() ); }
	int n_bounds()            const { return static_cast<int> (       lower_bound.size() ); }
	int n_initial_points()    const { return static_cast<int> (     initial_point.size() ); }

	const string print(const char letter, const int index);
	const string add_primitive(const primitive_rep& prim);
	const int get_index(const string& arg);

	void check_initial_point_consistency();

	string problem_name;

	list<primitive_rep> primitives;

	list<string> numeric_constants;
	vector<int>    defined_var_index;
	vector<int>    constraint_index;

	vector<string> lower_bound;
	vector<string> upper_bound;
	vector<string> initial_point;

};

#endif
