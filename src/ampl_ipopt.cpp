//
//  This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
//  Version: 0.01
//
// 	Last updated on: 20 Mar 2010
//
//  The original code came with the following license.
//
//  The original source has been altered.
//
//==============================================================================
//
// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: ampl_ipopt.cpp 1293 2008-08-25 16:06:13Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "AmplTNLP.hpp"
#include "IpIpoptApplication.hpp"

#include "IpoptConfig.h"
#ifdef HAVE_CSTRING
# include <cstring>
#else
# ifdef HAVE_STRING_H
#  include <string.h>
# else
#  error "don't have header file for string"
# endif
#endif

// for printf
#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

/* AMPL includes */
#include "asl.h"
#include "asl_pfgh.h"

namespace {
	Ipopt::SmartPtr<Ipopt::IpoptApplication> app = new Ipopt::IpoptApplication();
	Ipopt::SmartPtr<Ipopt::AmplSuffixHandler> suffix_handler = new Ipopt::AmplSuffixHandler();
	Ipopt::SmartPtr<Ipopt::TNLP> ampl_tnlp;
}

void init_ampl_ipopt_lib(char** argv) {

	using namespace Ipopt;

	ApplicationReturnStatus retval;

	// Call Initialize the first time to create a journalist, but ignore
	// any options file

	retval = app->Initialize("");

	if (retval != Solve_Succeeded) {
		printf("ampl_ipopt.cpp: Error in first Initialize!!!!\n");
		exit(-100);
	}

	// Add the suffix handler for scaling
	suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
	suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
	suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Objective_Source, AmplSuffixHandler::Number_Type);
	// Modified for warm-start from AMPL
	suffix_handler->AddAvailableSuffix("ipopt_zL_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
	suffix_handler->AddAvailableSuffix("ipopt_zU_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
	suffix_handler->AddAvailableSuffix("ipopt_zL_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
	suffix_handler->AddAvailableSuffix("ipopt_zU_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);

	ampl_tnlp = new AmplTNLP(ConstPtr(app->Jnlst()), app->Options(), argv, suffix_handler);

	// Call Initialize again to process output related options
	retval = app->Initialize();

	if (retval != Solve_Succeeded) {
		printf("ampl_ipopt.cpp: Error in second Initialize!!!!\n");
		exit(-101);
	}

	SmartPtr<OptionsList> options = app->Options();

	options->SetIntegerValue("print_level", 0);

	//options->SetStringValue("start_with_resto", "yes");
	//options->SetStringValue("nlp_scaling_method", "equilibration-based");
	//options->SetStringValue("linear_scaling_on_demand", "no");

	app->Initialize("");
}

int local_search(const double lb[], const double ub[], const double** x)
{

	static bool called_1st_time(true);
	
	using namespace Ipopt;

	ApplicationReturnStatus retval;

	AmplTNLP* raw_ptr = static_cast<AmplTNLP*>(GetRawPtr(ampl_tnlp));

	ASL_pfgh* asl = raw_ptr->AmplSolverObject();

	for (int i=0; i<n_var; ++i)
		havex0[i] = 1;

	for (int i=0; i<n_var; ++i)
		X0[i] = (lb[i]+ub[i])/2.0;

	for (int i=0; i<n_var; ++i) {
		LUv[2*i  ] = lb[i];
		LUv[2*i+1] = ub[i];
	}

	//===========================================================

	if (called_1st_time) {
		called_1st_time = false;
		retval = app->OptimizeTNLP(ampl_tnlp);
	}
	else
		retval = app->ReOptimizeTNLP(ampl_tnlp);

	if (retval == Solve_Succeeded || 
		retval == Infeasible_Problem_Detected || 
		retval == Solved_To_Acceptable_Level ||
		retval == Feasible_Point_Found )
	{

		(*x) = raw_ptr->get_sol();
	}

    if (retval == Solve_Succeeded) {
      //message = "Optimal Solution Found";

	  return 0;
    }
    else if (retval == Solved_To_Acceptable_Level) {
      //message = "Solved To Acceptable Level."; ??? Shouldn't we accept such solutions???

	  return 1;
    }
    else if (retval == Infeasible_Problem_Detected) {
      //message = "Converged to a locally infeasible point. Problem may be infeasible.";

	  return 2;
    }
	else if (retval == Feasible_Point_Found) {
		//std::cout << "Feasible point found" << std::endl;
	  return 0; // Accept feasible points as well -> stewgou40 !!!
	}
	else {
      //throw "Unknown Error";
	  std::cout << "Unknown error, code: " << retval << std::endl;
	  return -1;
    }

}

#include <sstream>

void write_sol_to_file(const int n) {

	using namespace Ipopt;
	using namespace std;

	ostringstream s;

	s << n << flush;

	const string sol_number = s.str();

	const string msg = "Solution #" + sol_number;

	AmplTNLP* raw_ptr = static_cast<AmplTNLP*>(GetRawPtr(ampl_tnlp));

	raw_ptr->write_solution_file(msg);

	ASL_pfgh* asl = raw_ptr->AmplSolverObject();

	string stub_name( filename );
	
	stub_name = stub_name.substr(0, stub_name.size()-4);

	const string old_name = stub_name + ".sol";

	const string new_name(stub_name + "_" + sol_number + ".sol");

	int result = rename(old_name.c_str(), new_name.c_str());
	if ( result != 0 )
		cerr << "Error: missing the -AMPL flag or the .nl extension?" << endl;

}

