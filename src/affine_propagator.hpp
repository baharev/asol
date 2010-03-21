//==============================================================================
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
//==============================================================================
//
// Version: 0.01
// Last updated on: 20 Mar 2010
//

#ifndef __AFFINE_PROPAGATOR_HPP
#define __AFFINE_PROPAGATOR_HPP

namespace aa_lib {

class affine;

class affine_propagator {

	public:

		affine_propagator(const int num_of_vars);

		void prune(const affine lin_cons[], affine vars[], const int num_of_vars, bool& toDelete, double& maxProgress);

		~affine_propagator();

	private:

		affine_propagator(const affine_propagator& );
		affine_propagator& operator=(const affine_propagator& );

		void reset_arrays();

		bool is_var_narrow(const affine vars[], const int i);

		bool check_for_convergence(const affine vars[]);

		int get_max_var_index(const affine& aa);

		double sum_aux_noise_vars(const affine& aa, const int max_var_idx);

		void sum_except_j(const affine& aa, const int j, const int n, const double sum_2, double& sum_inf, double& sum_sup);

		bool is_intersection_empty(const int idx, const double a_j, const double sum_inf, const double sum_sup);

		bool check_for_zero(const affine r[], const int begin, const int end);

		bool process_vars_in_eq(const affine& aa, const affine vars[]);

		bool process_all_eqs(const affine r[], const affine vars[], double& maxProgress);

		bool write_back(affine vars[], double& maxProgress);

		const int n_vars;

		double* const lo;
		double* const up;

		bool* const is_narrow;
};

}

#endif
