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

#include "affine.hpp"
#include "affine_propagator.hpp"
#include <iostream>
#include <cmath>
//#include <iomanip>
#include <assert.h>

using namespace std;

namespace aa_lib {

affine_propagator::affine_propagator(const int num_of_vars) :

	n_vars(num_of_vars),
	lo(new double[1+num_of_vars]),
	up(new double[1+num_of_vars]),
	is_narrow(new bool[1+num_of_vars])
{

}

affine_propagator::~affine_propagator() {

	delete[] lo;
	delete[] up;
	delete[] is_narrow;
}

void affine_propagator::reset_arrays() {

	for (int i=0; i<=n_vars; ++i) {
		lo[i] = -1.0;
		up[i] =  1.0;
	}

	for (int i=0; i<=n_vars; ++i)
		is_narrow[i] = false;
}

bool affine_propagator::is_var_narrow(const affine vars[], const int i) {

	const affine& x = vars[i-1];

	const bool prev_val = is_narrow[i];

	const bool narrow = aa_lib::is_narrow(x.inf(), x.sup());

	assert (!(prev_val == true && narrow == false));

	is_narrow[i] = narrow;

	return narrow;
}

bool affine_propagator::check_for_convergence(const affine vars[]) {

	bool has_converged = true;

	for (int i=1; i<=n_vars; ++i) {

		const affine& x = vars[i-1];

		assert(x.n == 2);
		assert(x.index[0] == 0);
		assert(x.index[1] == i);

		const bool narrow = is_var_narrow(vars, i);

		if (!narrow)
			has_converged = false;
	}

	return has_converged;
}

int affine_propagator::get_max_var_index(const affine& aa){

		assert(aa.n > 2);
		assert(aa.index[0] == 0);
		assert(aa.n <= aa.nmax);
		assert(aa.lb < aa.ub );

		const int m = aa.n-1; // number of noise symbols

		int max_var_idx = -1;
		for (int k=1; k<=m; ++k) {
			if (aa.index[k] <= n_vars)
				max_var_idx = k;
		}

		assert( max_var_idx >= 1 && max_var_idx <= n_vars);

		return max_var_idx;
}

double affine_propagator::sum_aux_noise_vars(const affine& aa, const int n) {

	const int m = aa.n-1; // number of all noise symbols
	//cout << "m = " << m << endl;
	//cout << "n = " << n << endl;
	double sum_2(0.0);

	for (int k=n+1; k<=m; ++k)
		sum_2 += fabs(aa.value[k]);

	//cout << "sum_2 = " << sum_2 << endl;
	return sum_2;
}

void affine_propagator::sum_except_j(const affine& aa, 
									   const int j, 
									   const int n, 
									   const double sum_2, 
									   double& sum_inf, 
									   double& sum_sup)
{
	//cout << endl << "j = " << j << endl << endl;

	sum_inf = aa.value[0];
	sum_sup = sum_inf;

	// sum of all non-aux noise symbols
	for (int i=1; i<=n; ++i) {

		//cout << "i = " << i << endl;
		if (i!=j) {

			const int index = aa.index[i];
			const double a_i = aa.value[i];

			double s_lo = 0.0;
			double s_up = 0.0;

			mult(a_i, a_i, lo[index], up[index], s_lo, s_up);

			sum_inf += s_lo;
			sum_sup += s_up;
		}
	}

	sum_inf -= sum_2;
	sum_sup += sum_2;
}

bool affine_propagator::is_intersection_empty(const int idx, const double a_j, const double sum_inf, const double sum_sup) {
		
	double eps_lo = -sum_inf/a_j;
	double eps_up = -sum_sup/a_j;

	if (eps_lo > eps_up)
		swap(eps_lo, eps_up);

	double chk_lo(0.0), chk_up(0.0);
	division(-sum_sup, -sum_inf, a_j, a_j, chk_lo, chk_up);
	
	assert ( chk_lo == eps_lo && chk_up == eps_up );

	//cout << "eps_j = " << "[ " << eps_lo << ", " << eps_up << "]" << endl;

	if (eps_lo > lo[idx]) {
		lo[idx] = eps_lo;
	}

	if (eps_up < up[idx]) {
		up[idx] = eps_up;
	}

	if (lo[idx] > up[idx]) {
		affine::isValid = false;
		return true;
	}

	return false;
}

// returns toDelete
bool affine_propagator::process_vars_in_eq(const affine& aa, const affine vars[]) {

	//const int m = aa.n-1; // number of noise symbols
	const int n = get_max_var_index(aa);

	const double sum_2 = sum_aux_noise_vars(aa, n);

	// for all noise vars in that eq
	for (int j=1; j<=n; ++j) {

		const int idx = aa.index[j];

		assert(idx>=1 && idx<=n_vars);

		if (is_var_narrow(vars, idx))
			continue;

		const double a_j = aa.value[j];

		if (a_j == 0.0) {
			//cout << endl << "Warning: a_j == 0.0" << endl;
			continue;
		}

		double sum_inf(0.0), sum_sup(0.0);

		sum_except_j(aa, j, n, sum_2, sum_inf, sum_sup);

		if (is_intersection_empty(idx, a_j, sum_inf, sum_sup)) {
			return true;
		}
	}

	return false;
}

bool affine_propagator::check_for_zero(const affine r[], const int begin, const int end) {

	assert(begin <= end);

	bool encloses_zero = true;

	for (int i=begin; i<=end; ++i) {

		const affine& aa = r[i];
		
		const double* const val = aa.value;
		const int   * const ind = aa.index;
		
		const int n = aa.n-1;
		
		double aa_lb = val[0];
		double aa_ub = aa_lb;
		
		for (int k=1; k<=n; ++k) {

			const double a_i = val[k];
			const int    idx = ind[k];

			double s_lo = 0.0;
			double s_up = 0.0;

			if (idx <= n_vars) {
				mult(a_i, a_i, lo[idx], up[idx], s_lo, s_up);
			}
			else {
				const double absval = fabs(a_i);
				s_lo = -absval;
				s_up =  absval;
			}

			aa_lb += s_lo;
			aa_ub += s_up;
		}
		// Should write back ???
		//if (aa_lb > aa.lb) aa.lb = aa_lb;
		//if (aa_ub < aa.ub) aa.ub = aa_ub;

		if (!((aa_lb <= 0.0) && (0.0 <= aa_ub))) {
			encloses_zero = false;
			break;
		}
	}			
			
	return encloses_zero;
}

bool affine_propagator::process_all_eqs(const affine r[], const affine vars[], double& maxProgress) {

	for (int p=0; p<n_vars; ++p) {

		bool has_converged = check_for_convergence(vars);

		if (has_converged) {
			maxProgress = 2.0;
			return false;
		}
		
		const affine& aa = r[p];

		//cout << scientific;
		//cout << endl << "Eq. #" << p << endl << aa << endl;

		const bool to_discard = process_vars_in_eq(aa, vars);

		if (to_discard) {
			affine::isValid = false;
			return true;
		}


		if (!check_for_zero(r, p, p)) {
			//cout << "Deleted" << endl;
			affine::isValid = false;
			return true;
		}
	}

	return false;
}

bool affine_propagator::write_back(affine vars[], double& maxProgress) {

	bool has_changed = false;

	for (int i=1; i<=n_vars; ++i) {

		const double progress = 1.0 - ((up[i]-lo[i])/2.0);

		if (progress > maxProgress)
			maxProgress = progress;

		// duplicate code, ouch!!!

		assert(lo[i] <= up[i]);

		if ((lo[i] > -0.99) || (up[i] < 0.99)) {

			has_changed = true;

			affine& x = vars[i-1];

			assert(x.n == 2);
			assert(x.index[0] == 0);
			assert(x.index[1] == i);

			double inf = 0.0;
			double sup = 0.0;

			const double a_i = x.value[1];

			mult(a_i, a_i, lo[i], up[i], inf, sup);

			const double a_0 = x.value[0];

			inf += a_0;
			sup += a_0;

			x.value[0] = (inf+sup)/2.0;
			x.value[1] = (sup-inf)/2.0;
			x.lb = inf;
			x.ub = sup;
		}
	}

	return has_changed;
}

// maxProgress >=0, <= 1; ==2; 
void affine_propagator::prune(const affine r[],
							  affine vars[],
							  const int num_of_vars,
							  bool& toDelete,
							  double& maxProgress)
{
	assert ( num_of_vars == n_vars ) ;

	maxProgress = 0.0;

	if (!affine::isValid) { 
		//cout << "Deleted" << endl;
		toDelete = true;
		return;
	}

	reset_arrays();

	bool to_discard = process_all_eqs(r, vars, maxProgress);

	if (to_discard) {
		affine::isValid = false;
		toDelete = true;
		return;
	}

	if (!check_for_zero(r, 0, n_vars-1)) {
		//cout << "Deleted" << endl;
		affine::isValid = false;
		toDelete = true;
		return;
	}

	bool has_changed = write_back(vars, maxProgress);

	if (!has_changed)
		maxProgress = 0.0;

}


}
