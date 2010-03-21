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

#ifndef __LP_IMPL_HPP
#define __LP_IMPL_HPP

#include "glpk.h"

// See the excellent documentation of GLPK for further details on the
// individual functions below

namespace lp_solver {

class lp_impl {

	public:

		explicit lp_impl(const int n_vars);

		double get_col_prim(int j) { return glp_get_col_prim(lp, j); }

		double get_obj_val()       { return glp_get_obj_val(lp); }

		int    get_it_cnt()        { return lpx_get_int_parm(lp, LPX_K_ITCNT); }

		int    get_status()        { return glp_get_status(lp); }

		void load_matrix(int nz, const int i[], const int j[], const double a[]) {
			
			glp_load_matrix(lp, nz, i, j, a);
		}

		void scale_prob(int flag) { glp_scale_prob(lp, flag); }

		void set_col_bnds(int j, int type, double lb, double ub) {
			
			glp_set_col_bnds(lp, j, type, lb, ub);
		}

		void set_row_bnds(int i, int type, double lb, double ub) {
			
			glp_set_row_bnds(lp, i, type, lb, ub);
		}

		void set_obj_coef(int j, double val) { glp_set_obj_coef(lp, j, val); }

		void set_msg_lev(const int level)    { parm->msg_lev = level; }

		void set_max_restart(const int);

		int  simplex()                       { return glp_simplex(lp, parm); }

		void std_basis()                     { glp_std_basis(lp); }

		int  term_out(int flag)              { return glp_term_out(flag); }

		int  warm_up()                       { return glp_warm_up(lp); }

		~lp_impl();

	private:

		lp_impl(const lp_impl& );

		lp_impl& operator=(const lp_impl& );

		glp_prob* lp;

		glp_smcp* parm;
};

const int DB = GLP_DB;
const int EQ = GLP_SF_EQ;
const int FX = GLP_FX;
const int MSG_ERR = GLP_MSG_ERR;
const int MSG_OFF = GLP_MSG_OFF;
const int MSG_ON  = GLP_MSG_ON;
const int NOFEAS  = GLP_NOFEAS;
const int OFF = GLP_OFF;
const int ON  = GLP_ON;
const int OPT = GLP_OPT;

}

#endif
