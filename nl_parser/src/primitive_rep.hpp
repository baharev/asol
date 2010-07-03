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

#ifndef PRIMITIVE_REP_HPP_
#define PRIMITIVE_REP_HPP_

class arg_t {

public:

	arg_t();
	arg_t(const char operand, const int index);

	char operand() const { return c; }
	int index() const { return i; }

	bool is_empty() const;

private:

	char c;
	int  i;
};

class primitive_rep {

public:

	primitive_rep(const char* operation, int arity_of_op, const arg_t& arg_1, const arg_t& arg_2 = arg_t());
	primitive_rep( ) : n_args(-1) { }

	const char* operation() const { return op_symbol; }
	int             arity() const { return n_args; }
	const arg_t     arg_1() const { return arg1; }
	const arg_t     arg_2() const { return arg2; }

private:

	char op_symbol[4];
	int    n_args;
	arg_t arg1;
	arg_t arg2;
};

#endif
