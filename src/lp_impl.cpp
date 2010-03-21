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

#include "lp_impl.hpp"

// See lp_pruning.cpp for further details on HACKED_GLPK

#define HACKED_GLPK
#ifdef HACKED_GLPK
extern "C" void set_restart_limit(const int limit);
#else
static void set_restart_limit(const int limit) { }
#endif

namespace lp_solver {

lp_impl::lp_impl(const int n_vars) {

	lp = glp_create_prob();

	glp_add_rows(lp, n_vars);

	glp_add_cols(lp, n_vars);

	parm = new glp_smcp;

	glp_init_smcp(parm);
}

lp_impl::~lp_impl() {

	delete parm;

	glp_delete_prob(lp);

	glp_free_env();
}



void lp_impl::set_max_restart(const int limit) {

	set_restart_limit(limit);
}

}

