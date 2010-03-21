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
// Last updated on: 21 Mar 2010
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <assert.h>
#include "affine.hpp"
#include "lp_pruning.hpp"
#include "lp_impl.hpp"

using namespace std;
using namespace lp_solver;

//-----------------------------------------------------------------------------
//
// LP_op (the _op suffix stands for operator) implements the LP pruning
// procedure, based on Achterberg�s heuristic. For further details see 
//
// http://dx.doi.org/10.1002/aic.11777  or 
//
// http://reliablecomputing.eu/baharevCP09.pdf
//
//-----------------------------------------------------------------------------
//
// If the LP problem is badly scaled GLPK tends to get into an infinite loop.
// I followed the advice of the GLPK developer, Andrew Makhorin; the tiny_coef
// function aims to improve the scaling of the LP problem by replacing the tiny
// coefficients with exact zeros.
// Unfortunately it is still not enough to avoid getting into an infinite loop.
// There is no official solution to this issue but a dirty workaround, see:
//
// http://lists.gnu.org/archive/html/help-glpk/2009-10/msg00101.html
//
// /* glpspx01.c (primal simplex method) */
//
// ...
//
// static int restart_counter = 0;
// static int restart_limit = 0;
//
// void set_restart_limit(const int limit) {
//	restart_limit = limit;
// }
//
// static int restart_counter = 0;
// static int restart_limit = 0;
//
// void set_restart_limit(const int limit) {
//	restart_limit = limit;
// }
//
// int spx_primal(glp_prob *lp, const glp_smcp *parm)
// {
// 	  ...
//
//       if (parm->msg_lev >= GLP_MSG_ERR)
//          xprintf("Warning: numerical instability (primal simplex,"
//             " phase %s)\n", csa->phase == 1 ? "I" : "II");
//
//       ++restart_counter;
//
//       if (restart_counter > restart_limit) {
//     	  store_sol(csa, lp, GLP_UNDEF, GLP_UNDEF, 0);
//     	  ret = GLP_EFAIL;
//     	  goto done;
//       }
//
//      /* restart the search */
//       csa->phase = 0;
//
// If you patch your version of GLPK according to the code above, you can
// #define HACKED_GLPK. Otherwise be prepared to get into an infinite loop. :(
//
// The GLPK developer, Andrew Makhorin, informed me that he is developing a
// more stable version of the simplex method (bound flipping ratio test - bound
// swapping dual simplex algorithm). This will eliminate the need for these
// ugly fixes (set_restart_limit in lp_impl.cpp and tiny_coef here).
//
//-----------------------------------------------------------------------------
//
// A degenerate bridge pattern is applied here. It is degenerate because there
// is only one implementor (lp_impl; "Cheshire Cat" class). For more details
// see
//
// Bridge pattern in Gamma, Helm, Johnson, Vlissides: Design Patterns;
//
// Minimizing compile time dependencies in:
// Item 26-28 (pimpl idiom) in Sutter: Exceptional C++;
// Item 31 in Meyers: Effective C++, Third Edition.
//
//-----------------------------------------------------------------------------

namespace aa_lib {

LP_op::LP_op(const int num_of_vars, const int n_nzeros) :

	n_vars(num_of_vars),
	lo(new double[1+num_of_vars]),
	up(new double[1+num_of_vars]),

	is_solved_min_x(new bool[1+num_of_vars]),
	is_solved_max_x(new bool[1+num_of_vars]),

	ia(new int[1+n_nzeros]),
	ja(new int[1+n_nzeros]),
	ar(new double[1+n_nzeros]),
	aij_max(new double[num_of_vars]),
	lp(new lp_impl(num_of_vars))
{

}

LP_op::~LP_op() {

	delete[] lo;
	delete[] up;

	delete[] is_solved_min_x;
	delete[] is_solved_max_x;

	delete[] ia;
	delete[] ja;
	delete[] ar;
	delete[] aij_max;

	delete lp;
}


bool LP_op::choose_candidate(int& index, int& direction) {

	bool hasCandidate = false;

	double exctr_min =  10.0;
	double exctr_max =  10.0;
	int    min_idx   = -1;
	int    max_idx   = -1;

	for (int i=1; i<=n_vars; ++i) {

		double tmp_val = lp->get_col_prim(i);

		const double val = (tmp_val>1.0)?1.0:((tmp_val<-1.0)?-1.0:tmp_val);		

		if (fabs(val+1.0) < tol_solved_lp)
			is_solved_min_x[i] = true;

		if (fabs(1.0-val) < tol_solved_lp)
			is_solved_max_x[i] = true;

		if (!is_solved_min_x[i]) {

			const double temp = fabs(val+1.0);

			if (temp < exctr_min) {
				exctr_min = temp;
				min_idx   = i;
				hasCandidate = true;
			}
		}

		if (!is_solved_max_x[i]) {

			const double temp = fabs(1.0-val);

			if (temp < exctr_max) {
				exctr_max = temp;
				max_idx   = i;
				hasCandidate = true;
			}
		}

	}

	if (hasCandidate) {

		if (exctr_min < exctr_max) {
			index     = min_idx;
			// MIN
			direction = 1;
		}

		else {
			index     = max_idx;
			// MAX
			direction = -1;
		}

		assert((1<=index)&&(index<=n_vars));
	}

	return hasCandidate;
}

bool LP_op::check_if_narrow(const affine vars[]) {

	bool isFinished(true);

	for (int i=0; i<n_vars; ++i) {

		const affine& tmp = vars[i];

		if ( is_narrow(tmp.inf(), tmp.sup()) ) {
			is_solved_min_x[i+1] = true;
			is_solved_max_x[i+1] = true;
		}
		else {
			is_solved_min_x[i+1] = false;
			is_solved_max_x[i+1] = false;
			isFinished = false;
		}

	}

	return isFinished;
}

// Algorithm of Tobias Achterberg
bool LP_op::lp_tobias() {

	bool is_ok = true;

	int index     = -1;
	int direction = -1;

	while (choose_candidate(index, direction) && is_ok) {

		lp->set_obj_coef(index, direction);

		lp->simplex();

		if (OPT==lp->get_status()) {

		  const double Z = lp->get_obj_val();

		  double z1(0.0), z2(0.0);

		  if (direction == 1.0) {

			  is_solved_min_x[index] = true;
			  z1 = Z;

			  if (z1 > -1.0) {
				  lo[index] = z1;
			  }
			  else if ((z1+1.0) < -1.0e-4)			  
				 is_ok = false;
		  }
		  else {

			  is_solved_max_x[index] = true;
			  z2 = -Z;

			  if (z2 < 1.0) {
				  up[index] = z2;
			  }
			  else if ((z2-1.0) > 1.0e-4)
				 is_ok = false;				  
		  }

		}
		else {
		  is_ok = false;
		}

		lp->set_obj_coef(index, 0.0);
	}

	return is_ok;
}

const bool DEBUG_TINY = false;

void LP_op::tiny_coef(const affine r[]) {

	for (int i=0; i<n_vars; ++i)
		aij_max[i] = 0.0;

	double min_aij = numeric_limits<double>::max();
	double max_aij = 0.0;

	for (int i=0; i<n_vars; ++i) {

		const int n = r[i].n;

		const int*    const idx = r[i].index;
		const double* const val = r[i].value;

		double b_lb = -val[0];
		double b_ub = b_lb;

		for (int j=1; j<n; ++j) {

			const int col = idx[j];

			if (col <= n_vars) {

				const double aij = fabs(val[j]);

				if (aij > aij_max[i])
					aij_max[i] = aij;

				if ((aij != 0.0) && (aij < min_aij))
					min_aij = aij;
			}
			else {
				const double tmp = fabs(val[j]);
				b_lb -= tmp;
				b_ub += tmp;
			}
		}

		//assert(aij_max[i] > 0.0);

		if (max_aij < aij_max[i])
			max_aij = aij_max[i];

		b_lb = fabs(b_lb);
		b_ub = fabs(b_ub);

		const double b_max = (b_ub > b_lb)? b_ub : b_lb ;

		if (b_max > aij_max[i])
			aij_max[i] = b_max;
	}

	assert(max_aij > 0.0);
	assert(min_aij < numeric_limits<double>::max());

	const double ratio(max_aij/min_aij);

	if (DEBUG_TINY) {
		cout << setprecision(3) << scientific << flush;
		cout << endl << "min_aij:       " << min_aij << ", max_aij:   " << max_aij << flush;
		cout << "  ratio:  " << ratio << endl;
	}

	for (int i=0; i<n_vars; ++i)
		aij_max[i] *= 1.0e-7;
}

void LP_op::build_lp(const affine r[]) {

	tiny_coef(r);

	//int tiny_aij(0);
	//int tiny_rbnd(0);

	int k = 1;
	
	for (int i = 0; i < n_vars; ++i) {

		const int n = r[i].n;

		const int*    const idx = r[i].index;
		const double* const val = r[i].value;

		double row_lb = -val[0];
		double row_ub = row_lb;

		const double tiny = aij_max[i];

		for (int j = 1; j < n; ++j) {

			const int col = idx[j];

			if (col <= n_vars) {

				const double val_j = val[j];

				if (fabs(val_j) >= tiny) {

					ia[k] = i+1;
					ja[k] = col;
					ar[k] = val_j;
					++k;
				}
				else if (DEBUG_TINY) {
					//++tiny_aij;
					if (val_j == 0.0)
						cout << "removed tiny a[" << i+1 << "][" << col << "] = " << val_j << endl;
				}
			}
			else {
				const double tmp = fabs(val[j]);
				row_lb -= tmp;
				row_ub += tmp;
			}
		}

		if (fabs(row_lb) < tiny) {
			//++tiny_rbnd;
			if (DEBUG_TINY)
				cout << "removed tiny row_lb #" << i+1 << " = " << row_lb << endl;
			row_lb = 0.0;
		}

		if (fabs(row_ub) < tiny) {
			//++tiny_rbnd;
			if (DEBUG_TINY)
				cout << "removed tiny row_ub #" << i+1 << " = " << row_ub << endl;
			row_ub = 0.0;
		}

		int row_type = FX;

		if (row_lb != row_ub)
			row_type = DB;

		lp->set_row_bnds(i+1, row_type, row_lb, row_ub);

	}

	for (int j=1; j<=n_vars; ++j) {
		lp->set_col_bnds(j, DB, -1.0, 1.0);
		lp->set_obj_coef(j, 0.0);
	}

	lp->load_matrix(k-1, ia, ja, ar);

	//if (tiny_aij || tiny_rbnd) {
	//	cout << fixed;
	//	cout << "Removed a_ij: " << tiny_aij << ", row_bnd: " << tiny_rbnd << endl;
	//}
}

bool LP_op::prune(const affine r[], affine vars[], bool& toDelete, const bool only_feas_check) {

	lp->set_max_restart(20);

	lp->term_out(OFF);

	if (DEBUG_TINY) {
		lp->term_out(ON);
	}

	build_lp(r);

	lp->scale_prob(EQ);

	lp->term_out(ON);

	//----------------------------
	//glp_bfcp bfcp;
	//glp_get_bfcp(lp, &bfcp);
	//bfcp.type = GLP_BF_GR;
	//bfcp.piv_tol = 0.70;
	//bfcp.nrs_max = 10;
	//glp_set_bfcp(lp, &bfcp);
	//----------------------------
	//glp_term_out(GLP_OFF);
	//glp_cpx_basis(lp);
	//glp_term_out(GLP_ON);
	//----------------------------
	//glp_write_lp(lp, 0, "dump_lp");
	//glp_write_mps(lp, GLP_MPS_FILE, 0, "dump_mps");
	//----------------------------

	lp->std_basis();

	const int had_problems = lp->warm_up();

	if (had_problems) {
		cout << "Warning: numerical problems building the LP object!" << endl;
		return false;
	}

	lp->set_msg_lev(MSG_OFF);

	if (DEBUG_TINY) {
		lp->set_msg_lev(MSG_ERR);
	}

	lp->simplex();

	const int status = lp->get_status();

	if (status == NOFEAS) {

		toDelete = true;
		return true;
	}
	else if (status == OPT) {

		if (only_feas_check)
			return true;

		const bool is_too_narrow = check_if_narrow(vars);

		if (is_too_narrow)
			return true;

		for (int i=0; i<=n_vars; ++i) {
			lo[i] = -1.0;
			up[i] =  1.0;
		}

		const bool is_ok = lp_tobias();

		bool is_consistent(true);

		for (int i=0; i<=n_vars && is_consistent; ++i) {
			if (lo[i] >= up[i])
				is_consistent = false;
		}

		if (!is_ok || !is_consistent) {
			cout << "Warning: numerical problems during lp pruning" << endl;
			return false;
		}

	}
	else {

		cout << "Warning: numerical problems in lp feasibilty test" << endl;
		return false;
	}

	for (int i=1; i<=n_vars; ++i) {

		assert(lo[i] <= up[i]);

		if ((lo[i] > -0.99) || (up[i] < 0.99)) {

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
	return true;
}

int LP_op::get_simplex_iteration_count() {
	return lp->get_it_cnt();
}

}
