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
//
//==============================================================================
//
// The code is a mess at the moment. It has to be refactored. I will follow the
// guidelines written in the following books:
//
// M Fowler; Refactoring
// WC Wake; Refactoring Workbook
// RC Martin; Agile Software Development
//
//
#include <iostream>
#include <iomanip>
#include <vector>
#include <deque>
#include <cmath>
#include <ctime>
#include <assert.h>
//#include <signal.h>

#define HAVE_CXSC
#ifdef  HAVE_CXSC
#include "interval_ext.hpp"
using namespace cxsc;
#endif

#include "affine.hpp"
#include "dag.hpp"

//#include "affine_propagator.hpp"

#define HAVE_GLPK
#ifdef  HAVE_GLPK
#include "lp_pruning.hpp"

#endif

using namespace std;
using namespace aa_lib;
using namespace dag_builder;

namespace {

	//--------------------------------------
	// The local solver needs these

	vector<double*> solutions;

	int Argc;
	char** Argv;
	double* pLb;
	double* pUb;

	//--------------------------------------

	int nbBox;
	int nbSplit;
	int depth;
	time_t tstart;
	deque<affine*> pending;

	affine* Res;

	interval* Ia_var;

	double* pWidth;
}

struct helper_classes {

	dag_builder::dag<affine>* aa_dag;
	dag_builder::dag<interval>* ia_dag;
	LP_op* LP_operator;
	//affine_propagator* aa_prop;
};

class width {

	public:

		explicit width(const double initial_width[],
					   const int n_vars,
					   const double tolerance);

		bool set_wprev(const double vars[]);

		bool rate_of_progress(const double current_width[],
							  const double ps,
							  const double pe,
							  bool& is_sufficient,
							  bool& is_excellent,
							  bool& has_converged,
							  int& index);

		const double* get_w0() const { return w; }

		double get_max_scaled_w() const { return max_prev; }

		double get_tol() const { return tol; }

		~width() { delete[] w; delete[] w_prev; }

	private:

		width();
		width(const width& );
		width& operator=(const width& );

		double* const w;

		double* const w_prev;

		double max_prev;

		const int n;

		const double tol;
};

width::width(const double initial_width[], const int n_vars, const double tolerance)
	: w(new double[n_vars]),
	  w_prev(new double[n_vars]),
	  max_prev(1.0),
	  n(n_vars),
	  tol(tolerance)
{
	for (int i=0; i<n_vars; ++i) {
		w[i] = initial_width[i];
		w_prev[i] = 1.0;
	}
}

bool width::set_wprev(const double current_width[]) {

	double max = 0.0;

	for (int i=0; i<n; ++i) {

		const double width = current_width[i]/w[i];

		assert((0.0 < width) && (width <= 1.0));

		w_prev[i] = width;

		if (width > max)
			max = width;
	}

	max_prev = max;

	return (max <= tol);
}

bool width::rate_of_progress(const double current_width[],
							 const double ps,
							 const double pe,
							 bool& is_sufficient,
							 bool& is_excellent,
							 bool& has_converged,
							 int& index )
{
	double max = 0.0;
	int idx(-2);
	is_sufficient = false;
	is_excellent = false;
	index = -1;

	for (int i=0; i<n; ++i){

		const double width = current_width[i]/w[i];

		//assert((0.0 < width) && (width <= 1.0));

		// favor higher indices ??? -> jaconsen.nl !!!
		if (  width > max) {
			max = width;
			idx = i;
		}

		const double w_diff = w_prev[i] - width;

		if ( w_diff > ps*max_prev) {
			is_sufficient = true;
		}

		w_prev[i] = width;
	}

	if ( max <= pe*max_prev )
		is_excellent = true;

	max_prev = max;

	has_converged = true;

	if (max > tol) {
		has_converged = false;
		index = idx+1;
	}

	assert( has_converged || (1 <=index && index<=n ) );

	return is_sufficient;
}

double elapsed_time() {

	return difftime(time(0), tstart);
}

void write_back(affine vars[], const dag_builder::dag<interval>& ia_dag) {

	affine::set_max_used_index(0);

	const int n_vars = ia_dag.number_of_vars();
	const interval* const ivars = ia_dag.get_variables();

	for (int i=0; i<n_vars; ++i) {
		vars[i] = affine(_double(Inf(ivars[i])), _double(Sup(ivars[i])));
	}
}

const double* comp_current_w(const affine vars[], const int n_vars) {

	double* const current_width = pWidth;

	for (int i=0; i<n_vars; ++i) {

		const double lb = vars[i].inf();
		const double ub = vars[i].sup();

		current_width[i] = ub-lb;
	}

	return current_width;
}

const double* comp_current_w(const interval ia_var[], const int n_vars) {

	double* const current_width = pWidth;

	for (int i=0; i<n_vars; ++i) {

		const double lb = _double(Inf(ia_var[i]));
		const double ub = _double(Sup(ia_var[i]));

		current_width[i] = ub-lb;
	}

	return current_width;
}

void ia_propagation(dag_builder::dag<interval>& ia_dag,
					affine vars[],
					bool& toDelete,
					int& index,
					width& w )
{
	const int n_vars = ia_dag.number_of_vars();

	// ??? makes an uncessary copy of the widths ...
	bool has_converged = w.set_wprev( comp_current_w(vars, n_vars) );

	if (has_converged) {
		index = -1;
		return;
	}

	interval* const ia_var = Ia_var;

	for (int i=0; i<n_vars; ++i)
		ia_var[i] = interval(vars[i].inf(), vars[i].sup());

	ia_dag.set_variables(ia_var);

	bool is_sufficient(true);
	bool is_excellent(false);

	for (int i=1; i<=n_vars && is_sufficient; ++i) {

		toDelete = ia_dag.propagate();

		if (toDelete)
			return;

		w.rate_of_progress( comp_current_w(ia_dag.get_variables(), n_vars), 0.025, 0.90, is_sufficient, is_excellent, has_converged, index);

	}

	write_back(vars, ia_dag);
}

//void aa_propagation(dag_builder::dag<affine>& aa_dag,
//					affine_propagator& aa_prop,
//					affine vars[],
//					bool& toDelete,
//					int& index,
//					width& w)
//{
//	double maxProgress(0.0);
//
//	const int n_vars = aa_dag.number_of_vars();
//
//	bool is_sufficient(true);
//	bool is_excellent(false);
//	bool has_converged(false);
//	affine* const r = ::Res;
//	//interval* const ia_var = ::Ia_var;
//
//	for (int i=1; i<=n_vars && is_sufficient && !has_converged; ++i) {
//
//		aa_dag.set_variables(vars);
//
//		if (toDelete = aa_dag.get_constraints(r))
//			return;
//
//		aa_prop.prune(r, vars, n_vars, toDelete, maxProgress);
//
//		if (toDelete)
//			return;
//
//		if (maxProgress <= 0.02) {
//			is_sufficient = false;
//		}
//		else if (maxProgress == 2.0) {
//			has_converged = true;
//		}
//		else { // maxProgress > 0.02
//
//			//ia_propagation(ia_dag, vars, toDelete, index, w);
//
//			//if (toDelete)
//			//	return;
//
//			//if (index == -1)
//			//	return;
//
//			w.rate_of_progress( comp_current_w(vars, n_vars), 0.05, 0.10, is_sufficient, is_excellent, has_converged, index);
//		}
//	}
//}

void g(affine vars[],
	   bool& toDelete,
	   int& index,
	   helper_classes& dags_lp,
	   width& w,
	   const bool only_feas_check = false)
{

	//bool is_sufficient(false);
	bool is_sufficient(true);
	bool is_excellent(false);
	//bool has_converged(true);
	bool has_converged(false);
	affine* const r = ::Res;
	//interval* const ia_var = ::Ia_var;

	dag_builder::dag<affine>&   aa_dag = *(dags_lp.aa_dag);
	dag_builder::dag<interval>& ia_dag = *(dags_lp.ia_dag);
	LP_op& LP_operator = *(dags_lp.LP_operator);

	const int n_vars = aa_dag.number_of_vars();

	//do {

	for (int k=1; ((k<=10 && is_sufficient)|| is_excellent)&& !has_converged; ++k) {

		is_sufficient = false;
		has_converged = true;

		ia_propagation(ia_dag, vars, toDelete, index, w);

		if (toDelete)
			return;

		if (index == -1)
			return;

		//===
		//aa_propagation(aa_dag, *(dags_lp.aa_prop), vars, toDelete, index, w);

		//if (toDelete)
		//	return;

		//if (index == -1)
		//	return;
		//===

		aa_dag.set_variables(vars);

		if (toDelete = aa_dag.get_constraints(r))
			return;

		const bool is_ok = LP_operator.prune(r, vars, toDelete, only_feas_check);
		if (!is_ok)
			break;

		if (toDelete)
			return;

		w.rate_of_progress( comp_current_w(vars, n_vars), 0.05, 0.90, is_sufficient, is_excellent, has_converged, index);
	}

	//} while (is_sufficient && (!has_converged));

	//-------------------------------------------------

	ia_propagation(ia_dag, vars, toDelete, index, w);
}

const bool APPLY_PROBING(true);

#define AMPL_IPOPT_LIB

#ifdef AMPL_IPOPT_LIB

extern int local_search(const double lb[], const double ub[], const double** x);
extern void init_ampl_ipopt_lib(char** argv);
extern void write_sol_to_file(const int n);

bool loc_solv(const affine vars[], const int n_vars)
{

	double* const lb = pLb;
	double* const ub = pUb;

	for (int i=0; i<n_vars; ++i) {
		lb[i] = vars[i].inf();
		ub[i] = vars[i].sup();
	}

	const double* x = 0;

	int ret_val = local_search(lb, ub, &x);

	bool new_solution_found(false);

	// Found a solution
	if (ret_val == 0) {

		bool push_back(true);

		const unsigned int size = solutions.size();

		for (unsigned int k=0; k<size && push_back; ++k) {

			double* y = solutions.at(k);

			bool is_different = false;

			for (int i=0; i<n_vars && !is_different; ++i) {
				// ??? change to rel and abs
				if (fabs(x[i]-y[i]) > 1.0e-4)
					is_different = true;
			}

			// found a solution already stored
			if (!is_different) {
				//cout << "This solution is already stored!" << endl;
				push_back = false;
			}
		}

		if (push_back) {

			new_solution_found = true;
			cout << setprecision(1) << fixed;
			cout << "A new solution found at box #" << nbBox << " ( depth: " << depth << ", " << flush;
			cout << " sol. #" << (solutions.size()+1) << "), " << elapsed_time() << " s" << endl;
			cout << setprecision(5) << fixed;
			cout.setf(ios_base::showpos);
			for (int i=0; i<n_vars; ++i)
				cout << x[i] << endl;
			cout.unsetf(ios_base::showpos);
			cout << setprecision(1) << fixed;
			cout << endl;

			double* sol = new double[n_vars];

			for (int i=0; i<n_vars; ++i)
				sol[i] = x[i];

			solutions.push_back(sol);

			write_sol_to_file(solutions.size());
		}
	}

	return new_solution_found;
}

#endif

bool is_too_narrow(const affine vars[], const int n_vars) {

	bool narrow(true);

	for (int i=0; i<n_vars && narrow; ++i) {

		const double lb = vars[i].inf();
		const double ub = vars[i].sup();

		if (!is_narrow(lb, ub))
			narrow = false;
	}

	return narrow;
}

void try_to_discard(affine vars[],
				   helper_classes& dags_lp,
				   width& w,
				   bool& toDelete,
				   char& status,
				   const bool narrow)
{
	dag_builder::dag<affine>&   aa_dag = *(dags_lp.aa_dag);
	dag_builder::dag<interval>& ia_dag = *(dags_lp.ia_dag);
	LP_op& LP_operator = *(dags_lp.LP_operator);

	// If the box is not too narrow, apply LP pruning
	const int n_vars = aa_dag.number_of_vars();

	if (!narrow) {

		affine* const r = Res;
		aa_dag.set_variables(vars);

		if (toDelete = aa_dag.get_constraints(r))
			return;

		LP_operator.prune(r, vars, toDelete);

	}

	// It is safe to apply IA propagation on narrow boxes
	if (!toDelete) {

		interval* const ia_var = Ia_var;
		for (int i=0; i<n_vars; ++i)
			ia_var[i] = interval(vars[i].inf(), vars[i].sup());

		ia_dag.set_variables(ia_var);

		bool is_sufficient(true);
		bool is_excellent(false);
		bool has_converged(true);
		int index = -1;

		for (int i=1; i<=n_vars && is_sufficient && !toDelete; ++i) {

			toDelete = ia_dag.propagate();

			w.rate_of_progress(comp_current_w(ia_dag.get_variables(), n_vars), 0.001, 0.90, is_sufficient, is_excellent, has_converged, index);

		}
	}

	if (!toDelete) {
		// Likely to be a parasite box
		status = '#';

		write_back(vars, ia_dag);
	}
	else {
		// Managed to get rid of the parasite box
		status = 'd';
	}
}

// return: contains
bool check_for_strict(const affine vars[], const double x[], const int n_vars, char& status) {

	bool contains = true;

	status = 'S';

	for (int i=0; i<n_vars; ++i) {

		const double lb = vars[i].inf();
		const double ub = vars[i].sup();

		const double fl = fabs(lb);
		const double fu = fabs(ub);

		const double xi = x[i];

		// Check easy containment first
		// Check abs tol
		const double threshold = 1.0e-6;
		if ((lb - threshold  <= xi ) && (xi <= ub + threshold)) {

			// Check strict containment
			if (xi < lb || xi > ub )
				status = 'E';
		}
		// Check rel tol
		else if ((lb - threshold*fl  <= xi ) && (xi <= ub + threshold*fu)) {

			status = 'E';
		}
		else {
			status = '#';
			contains = false;
			break;
		}
	}

	return contains;
}

// Check if the box contains a solution already stored
bool check_for_containment(const affine vars[], const int n_vars, char& status) {

	bool already_stored(false);

	status = '#';

	const int sol_sz = solutions.size();

	for (int j=0; j<sol_sz; ++j)	{

		const double* const x = solutions.at(j);

		bool contains = check_for_strict(vars, x, n_vars, status);

		if (contains) {
			assert ( (status == 'S') || (status == 'E') );
			already_stored = true;
			break;
		}
		else
			status = '#';
	}

	return already_stored;
}

void check_for_sol(affine vars[],
				   helper_classes& dags_lp,
				   width& w,
				   bool& toDelete,
				   char& status,
				   const bool narrow)
{

	const int n_vars = (dags_lp.aa_dag)->number_of_vars();

	const bool already_stored = check_for_containment(vars, n_vars, status);

	// This box seems not to contain an already stored solution
	// Perform local search first
	if (!already_stored) {

		bool new_solution = false;

#ifdef AMPL_IPOPT_LIB
		new_solution = loc_solv(vars, n_vars);
#endif

		if (!new_solution) {
			// Seems like there is no solution in the box
			// or we do not have local solver...
			try_to_discard(vars, dags_lp, w, toDelete, status, narrow);
		}
		else {
			cout << endl << "Found a new solution in check_for_sol()" << endl;
			check_for_strict(vars, *(solutions.rbegin()), n_vars, status);
			if (status == '#')
				cout << "Warning: IPOPT returned a strange solution!" << endl;
		}
	}

	assert (status != '?');
}

void init(const int n) {

	affine::set_valid();

	pLb = new double[n];
	pUb = new double[n];

	Res = new affine[n];

	Ia_var = new interval[n];

	pWidth = new double[n];
}

void destroy() {

	delete[] pLb; pLb = 0;
	delete[] pUb; pUb = 0;
	delete[] Res; Res = 0;
	delete[] Ia_var; Ia_var = 0;
	delete[] pWidth; pWidth = 0;
}

//double find_l_max(const affine vars[], const int index_set[], const int n_idx) {
//
//	assert (n_idx >= 2);
//
//	double l_max = vars[index_set[0]].inf();
//
//	for (int i=1; i<n_idx; ++i) {
//		const double l = vars[index_set[i]].inf();
//		if (l > l_max)
//			l_max = l;
//	}
//
//	return l_max;
//}

void check_index(const affine vars[], const int index, const double l_max, int& idx, double& w_max) {

	const double l = vars[index].inf();
	const double u = vars[index].sup();
	double       w = (u - l);
	//if (l >= l_max - 0.05 && w >= 0.05)
	//	w*=1.05;

	cout << "v" << index << ": " << flush;
	cout << " [" << l << ", " << u << "] dia: " << w << endl;

	if ( w > w_max ) {
		w_max = w;
		idx = index;
	}
}

void select_custom_index(const affine vars[], const int n, const width& w, int& index) {
	
	assert (n == 140);
	const int index_set[] = { 19, 39, 59 };
	const int n_index_set = sizeof(index_set)/sizeof(int);
	const int first_T_index = 120;

	int idx = -1;
	double w_max = -100.0;

	cout << setprecision(5) << fixed;

	//const double l_max = find_l_max(vars, index_set, n_index_set);

	for (int i=0; i<n_index_set; ++i) {

		check_index(vars, index_set[i], 0.0, idx, w_max);
	}

	if (w_max <= 0.05) {

		double maxw(-100.0);
		int maxi(-1);

		const double* const w0 = w.get_w0();

		for (int i=n-1; i>=0; --i) {

			const double w_i = (vars[i].sup() - vars[i].inf())/(i>=first_T_index?vars[i].sup():w0[i]);

			if (w_i > maxw) {

				maxw = w_i;
				maxi = i + 1;
			}
		}

		assert( (1<=maxi) && (maxi<=n));

		if (maxw > w.get_tol())
			index = maxi;
		else
			index = -1;
	}
	else
		index = idx + 1; // + 1 was missing !!!

	assert( (index == -1) || ((1<=index) && (index<=n)));

	if (index > 0) {
		const affine& t = vars[index-1];
		const double slb = t.inf();
		const double sub = t.sup();
		const double wi0 = w.get_w0()[index-1];
		const double swi = (sub-slb)/wi0;

		cout << endl << "Selected component: " << index << flush;
		cout << setprecision(5) << fixed;
		cout << " [" << slb << ", " << sub << "] scl: " << swi << endl;
		cout << setprecision(1) << fixed;
	}
}

void probing(affine box_orig[],
			 bool& toDelete,
			 int& index,
			 helper_classes& dags_lp,
			 width& w,
			 const int parts)
{
	assert(parts >= 2);

	//cout << endl << "=============================================================" << endl;
	//cout << "Probing starts with box: " << endl;
	//for (int i=0; i<N_VARS; ++i) {
	//	cout << setprecision(5) << fixed;
	//	cout << "[ " << (box_orig[i]).inf() << ", " << (box_orig[i]).sup() << "]" << endl;
	//}

	//dag_builder::dag<affine>&   aa_dag = *(dags_lp.aa_dag);
	//dag_builder::dag<interval>& ia_dag = *(dags_lp.ia_dag);
	//LP_op& LP_operator = *(dags_lp.LP_operator);

	const int n = (dags_lp.ia_dag)->number_of_vars();

	affine* const x = new affine[n];

	interval* const x_hull = new interval[n];

	bool has_changed(false);

	const int index_orig = index;

	for (int i=n-1; i>=0; --i) {
	//for (int i=0; i<n; ++i) {

		//cout << endl << "Probing component " << (i+1) << " of " << n << endl;

		const double lb = box_orig[i].inf();
		const double ub = box_orig[i].sup();

		bool is_initialized(false);

		for (int k=0; k<parts; ++k) {

			affine::set_valid();

			for (int j=0; j<n; ++j)
				x[j] = box_orig[j];

			const double p_lb = (((double) (parts-k  ))/((double) parts))*lb + (((double) (k  ))/((double) parts))*ub;
			const double p_ub = (((double) (parts-k-1))/((double) parts))*lb + (((double) (k+1))/((double) parts))*ub;

			x[i].set_bounds(p_lb, p_ub);

			bool to_delete = false;

			int idx = -1;

			//cout << endl << "Calling g(), component: " << (i+1) << ", part: " << (k+1)  << endl;

			//for (int i=0; i<N_VARS; ++i) {
			//	cout << setprecision(5) << fixed;
			//	cout << "[ " << (x[i]).inf() << ", " << (x[i]).sup() << "]" << endl;
			//}

			g(x, to_delete, idx, dags_lp, w, false);

			//cout << "After g()" << endl;

			//for (int i=0; i<N_VARS; ++i) {
			//	cout << setprecision(5) << fixed;
			//	cout << "[ " << (x[i]).inf() << ", " << (x[i]).sup() << "]" << endl;
			//}

			//cout << "to_delete: " << to_delete << endl << endl;

			if (to_delete) {
				//cout << "Dropped part " << (k+1) << " of " << parts << endl;
			}
			else {

				if (!is_initialized) {

					for (int j=0; j<n; ++j)
						x_hull[j] = interval(x[j].inf(), x[j].sup());

					is_initialized = true;
				}
				else {

					for (int j=0; j<n; ++j) {

						const double x_L = x[j].inf();
						const double x_U = x[j].sup();
						const real   h_L = Inf(x_hull[j]);
						const real   h_U = Sup(x_hull[j]);

						const real lb_hull = ( x_L < h_L )?( x_L ):( h_L );
						const real ub_hull = ( x_U > h_U )?( x_U ):( h_U );

						x_hull[j] = interval(lb_hull, ub_hull);
					}

				} // if (!is_initialized) {
			} // if (to_delete)
		} // for (int k=0; k<parts; ++k) {

		if (is_initialized) {

			//cout << setprecision(5) << fixed;
			//cout << endl << "After computing the hull: " << endl;
			//for (int j=0; j<n; ++j)
			//	cout << x_hull[j] << endl;
			//cout << endl;

			// write back x_hull to box_orig

			// check if the reduction is sufficient
			for (int j=0; j<n; ++j) {

				const double x_L = _double(Inf(x_hull[j]));
				const double x_U = _double(Sup(x_hull[j]));

				const double o_L = box_orig[j].inf();
				const double o_U = box_orig[j].sup();

				// ???
				assert(x_L >= o_L - 1.0e-6);
				assert(x_U <= o_U + 1.0e-6);

				const double suff_red = 0.01 * (o_U - o_L);

				if ((x_L >= o_L + suff_red)||(x_U + suff_red <= o_U)) {

					// ???
					has_changed = true;
					box_orig[j].set_bounds(x_L, x_U);

					const double excel_red = 0.499 * (o_U - o_L);

					if ((x_L >= o_L + excel_red)||(x_U + excel_red <= o_U)) {
						cout << setprecision(3) << fixed;
						cout << "Component " << (j+1) << " reduced by " << flush;
						cout << (100.0-(x_U-x_L)/(o_U-o_L)*100.0) << " % " << flush;
						cout << "on probing " << (i+1) << endl;
					}
				}
			}
		}
		else {
			cout << "Discarded during probing!" << endl;
			toDelete = true;
			break;
		}
	} // for (int i=0; i<n; ++i) {

	if (!toDelete && has_changed) {

		g(box_orig, toDelete, index, dags_lp, w);
		//ia_propagation(ia_dag, box_orig, toDelete, index, w);

	}
	else
		index = index_orig;

	//===================================================================================

	// Custom index selection rule for mss20
	//if (!toDelete && index != -1 ) {
	//	select_custom_index(box_orig, n, w, index);
	//}

	//==========================================================================

	//cout << endl << "Box on exiting probing:" << endl;
	//cout << setprecision(5) << fixed;
	//for (int j=0; j<n; ++j)
	//	cout << box_orig[j].inf() << '\t' << box_orig[j].sup() << endl;
	//cout << endl;

	delete[] x;
	delete[] x_hull;
}

void print_if_sol(affine vars[],
				  helper_classes& dags_lp,
				  width& w,
				  const bool narrow)
{

	bool toDelete(false);
	char status = '?';

	check_for_sol( vars, dags_lp, w, toDelete, status, narrow);

	if (!toDelete) {

		cout << setprecision(1) << fixed;
		cout << endl << "Not discarded (" << status << ") box #" << nbBox << " (" << depth << "), t: " << elapsed_time() << " s " << endl;

		const int n_vars = (dags_lp.aa_dag)->number_of_vars();

		for (int i=0; i<n_vars; ++i) {
			cout << setprecision(5) << scientific;
			cout.setf(ios_base::showpos);
			//cout << "[ " << (vars[i]).inf() << ", " << (vars[i]).sup() << "]" << endl;
			cout << (vars[i].inf()+vars[i].sup())/2.0 << '\t' << vars[i].inf() << '\t' << vars[i].sup() << endl;
			cout.unsetf(ios_base::showpos);
			cout << setprecision(1) << fixed;
		}
	}

	delete[] vars;
}

void index_sanity_check(const affine vars[], const int index, const int n_vars, const width& w) {

	assert((1<=index) && (index<= n_vars));

	double maxw = 0.0;
	int    maxi = -1;

	const double* const w_0 = w.get_w0();

	for (int i=0; i<n_vars; ++i) {
		const double wi = (vars[i].sup()-vars[i].inf())/w_0[i];
		if ( wi >= maxw) {
			maxw = wi;
			maxi =  i+1;
		}
	}

	assert( (0.0 < maxw) && (maxw <= 1.0) );
	assert((1<=maxi) && (maxi<= n_vars));

	const affine& t = vars[index-1];
	const double slb = t.inf();
	const double sub = t.sup();
	const double wi0 = w.get_w0()[index-1];
	const double swi = (sub-slb)/wi0;

	assert( fabs(maxw-swi) < 1.0e-14 );

	// Would fail: probing sets different values!!!
	// incorrect -> assert(w.get_max_scaled_w() == (sub-slb)/w_0);

	cout << endl << "Splitting component: " << index << flush;
	cout << setprecision(3) << fixed;
	cout << " [" << slb << ", " << sub << "] scl: " << swi << endl;
	cout << setprecision(1) << fixed;
}



//bool to_dump(false);
void nlss() {

	assert(Argc == 2 || Argc == 3);

	affine::set_max_used_index(0);

	cout << "Parsing input: " << Argv[1] << endl;
	cout << "Tolerance: " << tol_solved_lp << endl << endl;

	dag_builder::dag<affine> aa_dag(Argv[1]);
	dag_builder::dag<interval> ia_dag(Argv[1]);

	const int n_vars   = aa_dag.number_of_vars();

	const int n_nzeros = aa_dag.number_of_nzeros();

	LP_op LP_operator(n_vars, n_nzeros);
	//affine_propagator aa_prop(n_vars);

	helper_classes dags_lp;

	dags_lp.aa_dag = &aa_dag;
	dags_lp.ia_dag = &ia_dag;
	dags_lp.LP_operator = &LP_operator;
	//dags_lp.aa_prop = &aa_prop;

	assert(n_vars == ia_dag.number_of_vars());

	init(n_vars);

#ifdef AMPL_IPOPT_LIB
	init_ampl_ipopt_lib(Argv);
#endif

	affine* vars = new affine[n_vars];

	const affine* const aa_var(aa_dag.get_variables());

	for (int i=0; i<n_vars; ++i)
		vars[i] = affine(aa_var[i]);

	cout << setprecision(5) << fixed;

	cout.setf(ios_base::showpos);

	for (int i=0; i<n_vars; ++i)
		cout << "[ " << (vars[i]).inf() << ", " << (vars[i]).sup() << "]" << endl;

	cout.unsetf(ios_base::showpos);
	cout << setprecision(1) << fixed;

	time(&tstart);

	//=================================================

	bool toDelete = false;

	toDelete = ia_dag.propagate();
	if (toDelete)
		return;


	write_back(vars, ia_dag);

	affine* const r = Res;

	aa_dag.set_variables(vars);

	if (toDelete = aa_dag.get_constraints(r))
		return;

	LP_operator.prune(r, vars, toDelete);

	if (toDelete)
		return;

	interval* const ia_var = Ia_var;
	for (int i=0; i<n_vars; ++i)
		ia_var[i] = interval(vars[i].inf(), vars[i].sup());

	ia_dag.set_variables(ia_var);

	toDelete = ia_dag.propagate();
	if (toDelete)
		return;

	write_back(vars, ia_dag);

#ifdef AMPL_IPOPT_LIB

	loc_solv(vars, n_vars);

#endif

	width w( comp_current_w(vars, n_vars), n_vars, 1.0e-3);

	pending.push_back(vars);

	deque<int> bindices;

	bindices.push_back(-1);

	//===============================================

	while (!pending.empty()) {

		++depth;

		const unsigned int n_elems = pending.size();

		for (unsigned int k=1; k<=n_elems; ++k) {

			affine::set_max_used_index(n_vars);
			affine::set_valid();

			bool toDelete = false;
			int index = 0;

			++nbBox;

			affine* vars(pending.front());

			pending.pop_front();

			const int bisected_index_was(bindices.front());

			bindices.pop_front();

			cout << endl << "---------------------------------------------------------" << endl;
			cout << "Box number: " << nbBox << " (depth: "<< depth << "), previous index split was: " << flush;
			cout << bisected_index_was << endl;
			//if (nbBox == 3270) to_dump = true;
			//for (int i=0; i<n_vars; ++i) {
			//	cout << setprecision(5) << scientific;
			//	cout << "[ " << (vars[i]).inf() << ", " << (vars[i]).sup() << "]" << endl;
			//}

			g(vars, toDelete, index, dags_lp, w);

			//cout << endl;

			//for (int i=0; i<n_vars; ++i) {
			//	cout << setprecision(5) << scientific;
			//	cout << "[ " << (vars[i]).inf() << ", " << (vars[i]).sup() << "]" << endl;
			//}

			//cout << endl << "toDelete: " << toDelete << endl;
			//cout << "index: " << index << endl;

			if (toDelete) {
				delete[] vars;
				//return;
			}
			else {

				bool narrow = is_too_narrow(vars, n_vars);

				if ((index==-1) || narrow) {
					// the box is too narrow
					// we are not splitting it anymore
					print_if_sol(vars, dags_lp, w, narrow);
				}

				else {

					//======================================================

					if (APPLY_PROBING) {

						probing(vars, toDelete, index, dags_lp, w, 8);

						if (toDelete) {
							delete[] vars;
							continue;
						}

						narrow = is_too_narrow(vars, n_vars);

						if ((index==-1) || narrow) {
							// the box is too narrow
							// we are not splitting it anymore
							print_if_sol(vars, dags_lp, w, narrow);
							continue;
						}

					}

					//======================================================

					assert(1<=index && index<= n_vars);

					assert(!narrow);

#ifdef 	AMPL_IPOPT_LIB
					loc_solv(vars, n_vars);
#endif
					// disable it on custom splitting / scaling
					index_sanity_check(vars, index, n_vars, w);

					//=====================================================

					affine* const new_box = new affine[n_vars];

					//for (int i=0; i<n_vars; ++i)
					//	cout << fixed << "[ " << (vars[i]).inf() << ", " << (vars[i]).sup() << "]" << endl;

					bisect(vars, new_box, index, n_vars);

					pending.push_back(vars );
					pending.push_back(new_box);
					++nbSplit;

					bindices.push_back(index);
					bindices.push_back(index);
				} // if (index == -1)
			} // if (toDelete) {
		} //for (int k=1; k<=n_elems; ++k) {
	} // while (!pending.empty()) {

	cout << setprecision(1) << fixed;
	cout << endl << "Time: " << elapsed_time() << " s  (boxes: " << flush;
	cout << nbBox << ", splits: " <<  nbSplit << ", depth: " << depth << ", " << flush;
	cout << "simplex iterations: " << LP_operator.get_simplex_iteration_count() << ")" << endl;

	if (!solutions.empty()) {

		cout << endl << "Dumping the local solver\'s solutions..." << endl << endl;
		cout << setprecision(5) << fixed;
		for (int i=0; i<solutions.size(); ++i) {
				double* sol = solutions.at(i);
				cout << "Sol #" << i+1 << endl;

				cout.setf(ios_base::showpos);

				for (int j=0; j<n_vars; ++j)
					cout << sol[j] << endl;

				cout.unsetf(ios_base::showpos);

				delete[] sol;
				cout << endl;
		}
	}

	destroy();
}

// If you uncomment the lines below, you will get a broken program.
/*
void dump_info(int msg) {

	cout << endl << "Stopping on user request ... " << endl;
	cout << "No. of splits: " << nbSplit << endl;
	cout << "No. of examined boxes: " << nbBox << endl;
	cout << "No. of boxes pending: " << pending.size() << endl;
	cout << "at depth: " << depth << endl;
	cout << "Elapsed time : " << elapsed_time() << " s" << endl;

	time_t rawtime;
	time ( &rawtime );
	struct tm* timeinfo = localtime( &rawtime );

	cout << "Exited at: " << asctime(timeinfo) << endl;
	signal(msg, SIG_DFL);
	raise(msg);
}
*/

int main(int argc, char** argv) {

	//signal(SIGINT, dump_info);

	cout << "ASOL -- Affine solver. Version 0.01. Copyright (C) 2010 Ali Baharev" << endl;
	cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
	cout << "This is free software, you are welcome to redistribute it under the terms of" << endl;
	cout << "the GNU General Public License as published by the Free Software Foundation," << endl;
	cout << "either version 3 of the License, or (at your option) any later version." << endl;
	cout << "This program depends on:" << endl;
	cout << "GNU GLPK 4.39  (license: GNU  GPLv3)" << endl;
	cout << "C-XSC    2.4.0 (license: GNU LGPLv2)" << endl;
	cout << "IPOPT    3.6.1 (license: CPL)" << endl;
	cout << endl << "Built on " << __DATE__ << " " << __TIME__ << endl << endl;

	Argc = argc;
	Argv = argv;

	time_t rawtime;
	time( &rawtime );
	struct tm* timeinfo = localtime( &rawtime );
	cout << endl << "Started at: " << asctime(timeinfo) << endl;

	nlss();

	time ( &rawtime );
	timeinfo = localtime( &rawtime );
	cout << endl << "Finished at: " << asctime(timeinfo) << endl;

	cerr << endl << "SUCCESSFULLY FINISHED!" << endl << endl;

	//cin.get();

	return 0;
}

