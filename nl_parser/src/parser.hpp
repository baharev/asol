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

#ifndef PARSER_HPP_
#define PARSER_HPP_

#include <list>
#include "problem_rep.hpp"

using std::list;
using std::ifstream;

class parser {

public:

	parser(const char* const file_name);

	const problem_rep& get_problem_representation() const { return prob; }

private:

	typedef list<string>::iterator li;

	bool get_next_elem(li& i, string& arg);

	const li try_to_add_unary_primitive(const string& op, const li& beg) ;

	const li try_to_add_binary_primitive(const string& op, const li& beg);

	const li parse_nl_part(const li& nl_beg);

	const string add_as_additions(vector<string>& args);

	const li add_sumlist(const li& beg);

	const string parse_line_of_linear_part(const string& s);

	const string parse_linear_part(li& i, const int length);

	const string parse_rhs(li& i, const int constraint_index);

	const li try_to_add_primitive(const string& op, const int arity, const li& beg);

	int assemble_linear_nonlinear_parts(const string& nl_part, const string& linear_part);

	void parse_defined_variable(const li& segm_beg, const li& segm_end);

	void parse_segments();

	int get_index(const string& part);

	const string parse_J_segment(const li& pos, const int index);

	const string assemble_J_part_constant_term(const string& J_part, const string& constant_term);

	void parse_constraint(const li& segm_beg, const li& segm_end);

	void parse_variable_bounds(const li& segm_beg, const li& segm_end);

	void parse_initial_points(const li& segm_beg, const li& segm_end);

	void grab_content(ifstream& in);

	problem_rep prob;

	list<string> content;

};

#endif
