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

#ifndef PROBLEM_STATISTICS_HPP_
#define PROBLEM_STATISTICS_HPP_

class problem_statistics {

public:

	problem_statistics(int variables, int equations, int nonzeros, int defined_variables) {
		vars = variables;
		eqns = equations;
		nzeros = nonzeros;
		dfvs = defined_variables;
	}

	problem_statistics() { vars = eqns = nzeros = dfvs = -1; }

	int n_variables() const { return vars; }
	int n_equations() const { return eqns; }
	int n_nonzeros()  const { return nzeros; }
	int n_defined_variables() const { return dfvs; }

private:

	int vars;
	int eqns;
	int nzeros;
	int dfvs;

};

#endif
