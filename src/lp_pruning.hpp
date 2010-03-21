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

#ifndef __LP_PRUNING_HPP
#define __LP_PRUNING_HPP

namespace lp_solver {

	class lp_impl;
}

namespace aa_lib {

class affine;

class LP_op {

	public:

		LP_op(const int n_vars, const int n_nzeros);

		bool prune(const affine r[], affine vars[], bool& toDelete, const bool only_feas_check = false);

		int get_simplex_iteration_count();

		~LP_op();

	private:

		LP_op(const LP_op& );
		LP_op& operator=(const LP_op& );

		bool choose_candidate(int& index, int& direction);
		bool check_if_narrow(const affine vars[]);
		bool lp_tobias();
		void tiny_coef(const affine r[]);
		void build_lp(const affine r[]);

		const int n_vars;
		double* const lo;
		double* const up;

		bool* const is_solved_min_x;
		bool* const is_solved_max_x;

		int* const ia;
		int* const ja;
		double* const ar;
		double* const aij_max;
		
		lp_solver::lp_impl* lp;

};

}

#endif
