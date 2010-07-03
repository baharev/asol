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
#include "interpreter_DAG.hpp"
#include "problem_rep.hpp"
#include "utility.hpp"

using namespace std;

typedef list<string>::const_iterator li;
typedef list<primitive_rep>::const_iterator pli;

interpreter_DAG::interpreter_DAG(const problem_rep& problem, std::iostream& os) : prob(problem), out(os) {

	w_declarations();
	w_constraints();
	w_initial_point();
	w_logic();
}

void interpreter_DAG::w_declarations() {

	const int prm_sz = prob.n_primitives();
	const int var_sz = prob.stats.n_variables();
	const int con_sz = prob.n_constraints();
	const int nzeros = prob.stats.n_nonzeros();
	const int num_sz = prob.n_numeric_constants();


	const string problem_name = chop_extension_nl(prob.name());

	out << "problem name: " << problem_name << endl;
	out << "number of primitives -- t" << endl;
	out << prm_sz << endl;
	out << "number of variables -- v" << endl;
	out << var_sz << endl;
	out << "number of constraints -- c" << endl;
	out << con_sz << endl;
	out << "number of nonzeros" << endl;
	out << nzeros << endl;
	out << "number of numeric constants -- n" << endl;
	out << num_sz << endl;

	out << "numeric constant values" << endl;
	int i = 0;
	const li li_end( prob.numeric_constants.end() );
	for (li itr=prob.numeric_constants.begin(); i<num_sz; ++i, ++itr) {
		if (itr==li_end)
			error("unexpected end of container");
		out << (*itr) << endl;
	}

	out << "constraint indices in t" << endl;
	for (int i=0; i<con_sz; ++i) {
		out << prob.constraint_index.at(i) << endl;
	}
}

const string interpreter_DAG::arg_to_string(const arg_t& arg) {

	const char c = arg.operand();
	const int idx = arg.index();
	const string index = int_to_string(idx);

	if (idx < 0)
		error("incorrect index");

	string argument;

	if      (c == 'n') {
		if (idx >= prob.n_numeric_constants())
			error("incorrect n index");
		argument = "n\t" + index;
	}
	else if (c == 't') {
		if (idx >= prob.n_primitives())
			error("incorrect t index");
		argument = "t\t" + index;
	}
	else if (c == 'v') {
		const int n_vars = prob.stats.n_variables();

		if (idx < n_vars) {
			argument = "v\t" + index;
		}
		else {
			const int i = idx - n_vars;
			if (i >= prob.n_defined_vars())
				error("incorrect v index");
			const string def_idx = int_to_string(prob.defined_var_index.at(i));
			argument = "t\t" + def_idx;
		}
	}
	else {
		error("incorrect argument");
	}

	return argument;
}

void interpreter_DAG::w_unary_primitive(const int i, const primitive_rep& u_prim) {

	if (u_prim.arity()!=1)
		error("not an unary primitive");

	const arg_t arg1 = u_prim.arg_1();

	const string arg = arg_to_string(arg1);

	const string op = u_prim.operation();

	out << "1" << endl;
	out << op  << endl;
	out << i   << endl;
	out << arg << endl;
}

void interpreter_DAG::w_binary_primitive(const int i, const primitive_rep& b_prim) {

	if (b_prim.arity()!=2)
		error("not a binary primitive");

	const arg_t arg_1 = b_prim.arg_1();
	const arg_t arg_2 = b_prim.arg_2();

	const string arg1 = arg_to_string(arg_1);
	const string arg2 = arg_to_string(arg_2);

	const string op = b_prim.operation();

	out << "2" << endl;
	out << op  << endl;
	out << i   << endl;
	out << arg1 << endl;
	out << arg2 << endl;
}

void interpreter_DAG::w_constraints() {

	out << "constraints" << endl;

	const int n = prob.n_primitives();

	int i = 0;
	const pli li_end( prob.primitives.end() );
	for (pli itr = prob.primitives.begin(); i<n; ++i, ++itr) {

		if (itr==li_end)
			error("unexpected end of container");

		const primitive_rep& prim = *itr;

		if      (prim.arity()==1) {
			w_unary_primitive(i, prim);
		}
		else if (prim.arity()==2) {
			w_binary_primitive(i, prim);
		}
		else {
			error("incorrect arity");
		}
	}
}

void interpreter_DAG::w_initial_point() {

	const int n = prob.n_initial_points();

	if (n == 0) {
		error("testing requires an initial point");
	}

	out << "initial point x0" << endl;

	for (int i=0; i<n; ++i) {
		out << prob.initial_point.at(i) << endl;
	}
}

void interpreter_DAG::w_logic() {

	out << "EOF" << endl;
	out << endl;
}
