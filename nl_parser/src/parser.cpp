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

#include <fstream>
#include <iostream>
#include <cctype>
#include <algorithm>
#include "parser.hpp"
#include "utility.hpp"

using namespace std;

namespace {

void skip_line(ifstream& in, const int lines_to_skip = 1) {

	string buff;

	for (int i=1; i<=lines_to_skip; ++i) {
		getline(in, buff);
		buff.clear();
	}
}

problem_statistics parse_header(ifstream& in) {

	string line;
	getline(in, line);

	if (line[0] != 'g')
		error("only ASCII format is handled, write with option g");

	line.clear();

	//  n_ stands for "number of"
	int n_vars = -1;

	in >> n_vars;

	if (n_vars <=0)
		error("number of variables should be positive");

	int n_cons = -1;

	in >> n_cons;

	if (n_cons <=0)
		error("number of constraints should be positive");

	if (n_vars != n_cons)
		warning("the number of variables should equal the number of constraints");

	int n_objs(-1);

	in >> n_objs;

	if (n_objs != 0)
		error("cannot handle objectives");

	// ??? What are these?
	int n_rngs(-1);

	in >> n_rngs;

	if (n_rngs != 0)
		error("cannot handle ranges");

	int n_eqns(-1);

	in >> n_eqns;

	if (n_cons != n_eqns)
		error("only equality constraints are handled");

	skip_line(in, 5);

	// Check for discrete variables
	for (int i=1; i<=5; ++i) {

		int n_discr_var(-1);

		in >> n_discr_var;

		if (n_discr_var != 0)
			error("cannot handle discrete variables");
	}

	skip_line(in);

	// Get the number of nonzeros in the Jacobian
	int n_nzeros = -1;

	in >> n_nzeros;

	if (n_nzeros <=0)
		error("the Jacobian is empty");

	skip_line(in, 2);

	int n_dfvs = 0;

	for (int i=1; i<=5; ++i) {
		// It is not clear how the number of defined variables is coded :(
		// Making assumptions...
		int buff(0);
		in >> buff;

		if (i==2 || i==4) {

			if (buff < 0)
				error("unexpected error");

			n_dfvs += buff;
		}
		else if (buff != 0) {

			error("unexpected error");
		}
	}

	skip_line(in);

	return problem_statistics(n_vars, n_eqns, n_nzeros, n_dfvs);
}

void chop_comment(string& line) {

	const string::size_type pos = line.find('#');

	if (pos == string::npos)
		return;

	// loop back from # to first non-white space
	string::size_type i = pos-1;

	while ( isspace(line.at(i)) )
		--i;

	line.erase(i+1);

	return;
}

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

bool is_J_segment(const std::string& str) {

	const char c = str.at(0);
	return (c == 'J');
}

bool is_r_segment(const std::string& str) {

	const char c = str.at(0);
	return (c == 'r');
}

bool is_operator(const string& op) {

	return op.at(0) == 'o';
}

const string resolve_operator(const string& op, int& arity) {

	string res_op;

	const int invalid = -1;
	const int unknown  = 0;
	arity = invalid;

	if      ( op == "o0") {
		res_op = "+";
		arity  = 2;
	}
	else if ( op == "o1") {
		res_op = "-";
		arity  = 2;
	}
	else if ( op == "o2") {
		res_op = "*";
		arity  = 2;
	}
	else if ( op == "o3") {
		res_op = "/";
		arity  = 2;
	}
	else if ( op == "o5") {
		res_op = "pow";
		arity  = 2;
	}
	else if ( op == "o16") {
		res_op = "-";
		arity  = 1;
	}
	else if ( op == "o43") {
		res_op = "log";
		arity  = 1;
	}
	else if ( op == "o44") {
		res_op = "exp";
		arity  = 1;
	}
	else if ( op == "o54") {
		res_op = "sumlist";
		arity  = unknown;
	}
	else if ( !is_operator(op) ) {

		error("input is not an operator");
	}
	else {
		string msg(op);
		msg += " not implemented yet";
		error(msg.c_str());
	}

	if( arity == invalid)
		error("unexpected error");

	return res_op;
}

bool is_nl_part_zero(const string& nl_part) {

	if (nl_part.at(0) != 'n')
		return false;

	const string val = nl_part.substr(1);

	const double value = atof(val.c_str());

	if (value == 0.0)
		return true;
	else
		return false;
}

typedef list<string>::iterator li;

}

bool parser::get_next_elem(li& i, string& arg) {

	const bool NOT_ARG = false;
	const bool  IS_ARG = true;

	++i;

	arg = *i;

	if (is_operator(arg))
		return NOT_ARG;

	return IS_ARG;
}

const li parser::try_to_add_unary_primitive(const string& op, const li& beg) {

	li i(beg);

	string arg;

	bool is_arg = get_next_elem(i, arg);

	if (!is_arg)
		return i;

	*beg = prob.add_unary_primitive(op, arg);

	return i;
}

const li parser::try_to_add_binary_primitive(const string& op, const li& beg) {

	li i(beg);

	string arg1;

	bool is_arg1 = get_next_elem(i, arg1);

	if (!is_arg1)
		return i;

	string arg2;

	bool is_arg2 = get_next_elem(i, arg2);

	if (!is_arg2)
		return i;

	*beg = prob.add_binary_primitive(op, arg1, arg2);

	return i;
}

const string parser::add_as_additions(vector<string>& args) {

	const int length = static_cast<int> ( args.size() );

	string t_prev = args.at(0);

	for (int k=1, new_index=prob.n_primitives(); k<length; ++k, ++new_index) {

		t_prev = prob.add_binary_primitive("+", t_prev, args.at(k));
	}

	return t_prev;
}

const li parser::add_sumlist(const li& beg) {

	li i(beg);

	if (*i != "o54")
		error("not a sumlist");

	++i;
	const int length = atoi( (*i).c_str() );

	if (length < 3)
		error("length of sumlist is incorrect");

	vector<string> args(length, "");

	for (int k=0; k<length; ++k) {
		++i;
		i = parse_nl_part(i);
		args.at(k) = *i;
	}

	*beg = add_as_additions(args);

	return i;
}

const string parser::parse_line_of_linear_part(const string& s) {

	const string::size_type space = s.find(" ");

	const string idx(s.substr(0, space));
	const string val(s.substr(space+1));

	const double value = atof(val.c_str());

	string arg = "v"+idx;

	if (value==-1) {
		// t = -v
		arg = prob.add_unary_primitive("-", arg);
	}
	else if (value==0) {
		// skip
		arg = "";
	}
	else if (value !=1 ) {
		// t = c*v
		arg = prob.add_binary_primitive("*", "n"+val, arg);
	}
	// else --> value == 1, just return the var

	return arg;
}

const string parser::parse_linear_part(li& i, const int length) {

	string arg = "";

	for (int k=0; k<length; ++k, ++i) {

		const string temp = parse_line_of_linear_part(*i); // 21 2 -> temp = n2*v21

		if (temp.length() == 0)
			continue; // 21 0 -> temp = 0*v21 -> skipping

		if (arg.length() == 0) {
			arg = temp;
		}
		else {
			arg = prob.add_binary_primitive("+", arg, temp); // t23 = t21 + t22
		}
	}

	return arg;
}

const li parser::try_to_add_primitive(const string& op, const int arity, const li& beg) {

	li i;

	if      (arity == 0) {

		i = add_sumlist(beg);
	}
	else if (arity == 1) {

		i = try_to_add_unary_primitive(op, beg);
	}
	else if (arity == 2) {

		i = try_to_add_binary_primitive(op, beg);
	}
	else {

		error("impossible arity");
	}

	return i;
}

const li parser::parse_nl_part(const li& nl_beg) {

	li i = nl_beg;

	string expr = *i;

	while (expr.at(0) == 'o') {

		li beg = i;

		int arity = -1;

		const string op = resolve_operator(expr, arity);

		const int size = prob.n_primitives();

		i = try_to_add_primitive(op, arity, beg);

		if (prob.n_primitives() != size) {
			++beg;
			++i;
			content.erase(beg, i);
			i = nl_beg;
		}

		expr = *i;
	}

	return i;
}

int parser::get_index(const string& part) {

	if (part.at(0) != 't')
		error("trivial cases are not supported");

	const string idx = part.substr(1);

	const int index = atoi( idx.c_str() );

	if (index != prob.n_primitives()-1)
		error("inconsistent indices");

	return index;
}

int parser::assemble_linear_nonlinear_parts(const string& nl_part, const string& linear_part) {

	bool nl_part_empty = is_nl_part_zero(nl_part);

	bool linear_part_empty = (linear_part.length() == 0);

	if (nl_part_empty && linear_part_empty)
		error("unexpected error");

	int index = -1;

	if (linear_part_empty) {

		index = get_index(nl_part);
	}
	else if (nl_part_empty) {

		index = get_index(linear_part);
	}
	else {
		// t32 = nl_part + linear_part
		index = prob.n_primitives();
		prob.add_binary_primitive("+", nl_part, linear_part);
	}

	return index;
}

void parser::parse_defined_variable(const li& segm_beg, const li& segm_end) {

	const string s(*segm_beg); // s is like "V21 3 2"

	const string::size_type first_space  = s.find(" ");
	const string::size_type second_space = s.find(" ", first_space+1);

	const string idx(s.substr(1, first_space-1)); // "21"
	const string lin(s.substr(first_space+1, second_space-first_space-1)); // "3"

	const int defined_var_index  = atoi( idx.c_str() );
	const int linear_part_length = atoi( lin.c_str() );

	li i = segm_beg;

	++i;

	const string linear_part = parse_linear_part(i, linear_part_length);

	li last = parse_nl_part(i);

	const string nl_part = *last;

	if (++last != segm_end)
		error("unexpected error");

	const int primitive_index = assemble_linear_nonlinear_parts(nl_part, linear_part);

	prob.add_defined_variable(defined_var_index, primitive_index);

	return;
}

const string parser::parse_J_segment(const li& pos, const int index) {

	li i(pos);

	for (int k=0; k<=index; ++k) {

		++i;
		i = find_if(i, content.end(), is_J_segment);

		if (i==content.end())
			error("failed to find J segment");
	}

	li J_beg(i);

	++i;
	const li J_end = find_if(i, content.end(), is_J_segment);

	const string s(*J_beg); // s is like "J23 4"

	const string::size_type space  = s.find(" ");

	const string idx(s.substr(1, space-1)); // "23"
	const string lin(s.substr(space+1));    // "4"

	const int J_index  = atoi( idx.c_str() );
	const int J_length = atoi( lin.c_str() );

	if (index != J_index)
		error("inconsistent index in J segment");

	++J_beg;

	const string linear_part = parse_linear_part(J_beg, J_length);

	return linear_part;
}

const string parser::parse_rhs(li& i, const int constraint_index) {

	i = find_if(i, content.end(), is_r_segment);

	if (i==content.end())
		error("failed to find r segment");

	advance(i, constraint_index+1);

	const string s(*i);

	const string::size_type space  = s.find(" ");

	if (space == string::npos)
		error("unexpected error while parsing r segment");

	const string type(s.substr(0, space));
	const string val( s.substr(space+1));

	if (type != "4")
		error("sorry, only equations are supported");

	const double value = -atof( val.c_str() );

	if (value==0.0)
		return "";

	const char first_char = val.at(0);

	const bool is_digit      = (isdigit(first_char) != 0);
	const bool is_minus_sign = (first_char == '-');

	string arg;

	if (is_digit) {
		arg = "n-" + val;
	}
	else if (is_minus_sign) {
		arg = 'n' + val.substr(1);
	}
	else {
		error("unexpected error while parsing r segment");
	}

	return arg;
}

const string parser::assemble_J_part_constant_term(const string& J_part, const string& constant_term) {

	if ((J_part.length() == 0) || (constant_term.length() == 0))
		return J_part + constant_term;

	return prob.add_binary_primitive("+", J_part, constant_term);
}

void parser::parse_constraint(const li& segm_beg, const li& segm_end) {

	const string s(*segm_beg);

	const string idx(s.substr(1));

	const int constraint_index = atoi( idx.c_str() );

	li i = segm_beg;

	++i;

	li last = parse_nl_part(i);

	const string nl_part = *last;

	if (++last != segm_end)
		error("unexpected error");

	i = segm_end;

	const string constant_term = parse_rhs(i, constraint_index);

	const string J_part = parse_J_segment(i, constraint_index);

	const string linear_part = assemble_J_part_constant_term(J_part, constant_term);

	const int primitive_index = assemble_linear_nonlinear_parts(nl_part, linear_part);

	prob.add_constraint(constraint_index, primitive_index);

	return;

}

void parser::parse_variable_bounds(const li& segm_beg, const li& segm_end) {

	int variable_counter = 0;

	li i(segm_beg);
	++i;

	const string double_bounded_variable("0 ");

	while (i != segm_end) {

		const string line(*i);

		if (line.substr(0, 2) != double_bounded_variable)
			error("only double bounded variables are supported");

		string::size_type pos = line.find(" ", 2);

		if (pos == string::npos || pos == 2)
			error("unexpected error while parsing b segment");

		const string lower_bound(line.substr(2, pos-2));
		const string upper_bound(line.substr(pos+1)    );

		prob.push_back_variable_bounds(lower_bound, upper_bound);

		++i;
		++variable_counter;
	}

	if (variable_counter != prob.stats.n_variables())
		error("inconsistent number of variables");
}

void parser::parse_initial_points(const li& segm_beg, const li& segm_end) {

	const string number_of_lines = (*segm_beg).substr(1);

	const int n_lines = atoi( number_of_lines.c_str() );

	if (prob.stats.n_variables() != n_lines)
		error("provide initial estimate for all variables or for none");

	int variable_counter = 0;

	li i(segm_beg);

	++i;

	while (i != segm_end) {

		const string line(*i);

		string::size_type pos = line.find(" ");

		if (pos == string::npos || pos == 0)
			error("unexpected error while parsing x segment");

		const string variable_index(line.substr(0, pos));
		const string initial_estimate(line.substr(pos+1));

		const int    idx = atoi(   variable_index.c_str() );

		prob.push_back_initial_point(idx, initial_estimate);

		++i;
		++variable_counter;
	}

	if (variable_counter != prob.stats.n_variables())
		error("inconsistent number of variables");
}

void parser::parse_segments() {

	li pos = content.begin();

	while (pos != content.end()) {

		li segm_beg = pos;

		++pos;

		li segm_end = find_if(pos, content.end(), is_segment);

		const char segment_type = segm_beg->at(0);

		if      (segment_type == 'V') {

			parse_defined_variable(segm_beg, segm_end);
		}
		else if (segment_type == 'C') {
			// parses the corresponding J segment and r value as well
			parse_constraint(segm_beg, segm_end);
		}
		else if (segment_type == 'b') {

			parse_variable_bounds(segm_beg, segm_end);
		}
		else if (segment_type == 'x') {

			parse_initial_points(segm_beg, segm_end);
		}
		else {

		}

		pos = segm_end;
	}

	prob.check_consistency();

}

void parser::grab_content(ifstream& in) {

	string line;

	do {
		line.clear();
		getline(in, line);

		if (line.length() == 0)
			break;

		chop_comment(line);
		content.push_back(line);

	} while ( !in.eof() );

}

parser::parser(const char* const file_name) {

	ifstream in(file_name);

	if (!in.good())
		error("cannot read input from file");

	cout << "Reading input file " << file_name << " ..." << endl;

	prob = problem_rep(file_name, parse_header(in));

	grab_content(in);

	in.clear();
	in.close();

	const int HEADER_LENGTH = 10;

	cout << "Successfully read " << (content.size()+HEADER_LENGTH) << " lines" << endl;

	cout << "Parsing segments ... " << endl;

	parse_segments();

	content.clear();

	cout << "Finished parsing segments" << endl;

}
