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

#include <cstring>
#include "primitive_rep.hpp"
#include "utility.hpp"

using namespace std;

arg_t::arg_t() : c('?'), i(-1) {

}

arg_t::arg_t(const char operand, const int index) : c(operand), i(index) {

	if (!((operand == 'n') || (operand == 'v') || (operand == 't')))
		error("incorrect operand in arg_t ctor");

	if (index < 0)
		error("incorrect index in arg_t ctor");
}

bool arg_t::is_empty() const {

	return c == '?';
}

primitive_rep::primitive_rep(const char* operation, const int arity_of_op, const arg_t& arg_1, const arg_t& arg_2)
: n_args(arity_of_op), arg1(arg_1), arg2(arg_2)
{

	memset(op_symbol, '\0', sizeof(op_symbol));
	strncpy(op_symbol, operation, sizeof(op_symbol));

	if (op_symbol[sizeof(op_symbol)-1] != '\0')
		error("max length of the operator symbol is 3");

	if ((n_args < 1) || (n_args > 2))
		error("incorrect number of arguments");

	if ( (n_args == 1) && !arg2.is_empty() )
		error("inconsistent arguments in unary primitive");

	if ( (n_args == 2) && (arg1.is_empty() || arg2.is_empty()) )
		error("inconsistent arguments in binary primitive");
}
