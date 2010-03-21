//=============================================================================
//
//  This code is part of ASOL (nonlinear system solver using affine arithmetic)
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

#include <fstream>
#include "dag.hpp"

namespace dag_builder {

void error(const char* msg) {
	std::cerr << std::endl << "Error: " << msg << "!" << std::endl;
	exit(-1);
}

//=============================================================================
//
// IMPORTANT!!!
//
// Everything from this point is subject to a COMPLETE REWRITE FROM SCRARTCH 
// by applying the interpreter pattern (see Gamma, Helm, Johnson, Vlissides:
// Design Patterns). More on this in the dag.hpp header file.
//
//=============================================================================

bool is_segment(const std::string& str) {

	const char c = str.at(0);
	return (c == 'F' ||
		    c == 'S' ||
			c == 'V' ||
			c == 'C' ||
			c == 'L' ||
			c == 'O' ||
			c == 'd' ||
			c == 'x' ||
			c == 'r' ||
			c == 'b' ||
			c == 'k' ||
			c == 'J' ||
			c == 'G' );
}

void read_file(const char* file_name, std::vector<std::string>& content, counters& s) {

	using std::ifstream;
	using std::string;

	int& n_vars   = s.n_vars;
	int& n_nzeros = s.n_nzeros;
	int& n_dfvs   = s.n_dfvs;
	int& n_node   = s.n_node;
	int& tmp_idx  = s.tmp_idx;
	int& num_idx  = s.num_idx;
	int& num_cnt  = s.num_cnt;

	n_vars = 0;
	n_dfvs = 0;
	n_node = 0;
	tmp_idx= 0;
	num_idx= 0;
	num_cnt= 0;

	ifstream in(file_name);

	if (!in.good())
		error("cannot read input from file");

	// Skip the first line
	string line;
	getline(in, line);
	line.clear();

	n_vars = 0;

	in >> n_vars;

	if (n_vars <=0)
		error("number of variables should be positive");

	int n_cons(-1);

	in >> n_cons;

	if (n_cons <=0)
		error("number of constraints should be positive");

	if (n_vars != n_cons)
		error("the number of variables should equal the number of constraints");

	int n_objs(-1);

	in >> n_objs;

	if (n_objs != 0)
		error("cannot handle objectives");

	// ???
	int n_rngs(-1);

	in >> n_rngs;

	if (n_rngs != 0)
		error("cannot handle ranges");

	int n_eqns(-1);

	in >> n_eqns;

	if (n_vars != n_eqns)
		error("the number of variables should equal the number of equations");

	// Skip the rest of the line and the next four lines
	for (int i=1; i<=5; ++i) {

		getline(in, line);
		line.clear();
	}

	// Check for discrete variables
	// exit if found any
	for (int i=1; i<=5; ++i) {

		int n_discr_var(-1);

		in >> n_discr_var;

		if (n_discr_var != 0)
			error("cannot handle discrete variables");
	}

	// Skip the rest of the line
	getline(in, line);
	line.clear();

	// Get the number of nonzeros in the Jacobian
	in >> n_nzeros;

	// Skip the rest of the line
	getline(in, line);
	line.clear();

	// Skip the next line
	getline(in, line);
	line.clear();

	n_dfvs = 0;

	for (int i=1; i<=5; ++i) {

		int buff(0);
		in >> buff;
		if (i==2 || i==4)
			n_dfvs += buff;
	}

	// Skip the rest of the line
	getline(in, line);

	do {
		line.clear();
		getline(in, line);
		if (line == "")
			break;
		content.push_back(line);

	} while ( !in.eof() );

	in.good();
	in.clear();
	in.close();
}

const std::string resolve_op(const std::string& op, int& n_args) {

	std::string res_op;
	n_args = -1;

	if      ( op == "o0") {
		res_op = "plus";
		n_args = 2;
	}
	else if ( op == "o1") {
		res_op = "minus";
		n_args = 2;
	}
	else if ( op == "o2") {
		res_op = "mult";
		n_args = 2;
	}
	else if ( op == "o3") {
		res_op = "div";
		n_args = 2;
	}
	else if ( op == "o5") {
		res_op = "pow";
		n_args = 2;
	}
	else if ( op == "o16") {
		res_op = "neg";
		n_args = 1;
	}
	else if ( op == "o43") {
		res_op = "ln";
		n_args = 1;
	}
	else if ( op == "o44") {
		res_op = "exp";
		n_args = 1;
	}
	else if ( op == "o54") {
		res_op = "sumlist";
		n_args = 0;
	}
	else {
		std::string msg(op);
		msg += " not implemented yet";
		error(msg.c_str());
	}

	assert ( n_args != -1);
	return res_op;
}

}
