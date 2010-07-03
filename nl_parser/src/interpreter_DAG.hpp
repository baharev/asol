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

#ifndef INTERPRETER_DAG_HPP_
#define INTERPRETER_DAG_HPP_

#include <iosfwd>
#include <string>

class problem_rep;
class primitive_rep;
class arg_t;

class interpreter_DAG {

public:

	interpreter_DAG(const problem_rep& problem, std::iostream& os);

private:

	void w_declarations();

	const std::string arg_to_string(const arg_t& arg);

	void w_unary_primitive(const int primitive_index, const primitive_rep& u_prim);

	void w_binary_primitive(const int primitive_index, const primitive_rep& b_prim);

	void w_constraints();

	void w_initial_point();

	void w_logic();

	const problem_rep& prob;

	std::iostream& out;
};

#endif
