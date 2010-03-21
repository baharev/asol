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

#ifndef __DAG_HPP
#define __DAG_HPP

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <assert.h>

namespace dag_builder {

template <typename T>
class node;

template <typename T>
class dag {

	public:

		explicit dag(const char* file_name);

		void set_variables (const T var[]);

		const T* get_variables() const { return var; }

		bool get_constraints(T* r);

		bool propagate();

		int number_of_vars() const { return num_of_vars; }

		int number_of_nzeros() const { return num_of_nzeros; }

		~dag();

	private:

		dag<T> (const dag<T>&);
		dag<T>& operator=(const dag<T>& );

		void evaluate_all();
		bool propagate_all();

		node<T>** Node;
		T* var;
		T* tmp;
		T* num;

		int* dfv;
		int* con;

		int num_of_vars;
		int num_of_nzeros;
		int num_of_nodes;
};


//=============================================================================
//
//               DO NOT USE ANYTHING DECLARED BELOW THIS POINT
//
//=============================================================================
//
// You might be tempted to implement the directed acyclic graph (DAG) as a 
// node-based container (nodes storing pointers to the parent and child nodes
// besides the actual data). In my opinion, there is a much simpler (and
// probably more efficient) way to do this. Think about it: in this particular
// application, after building the DAG in the memory, we never insert or delete
// a node; we just perform evaluation or propagation (a particular traversal).
// We only need to have enough information to perform the traversal. We can
// perform both evaluation and propagation in a simple for loop if the nodes
// are placed in an appropriately chosen order in a memory block. In other
// words: we do not need any bookkeeping how the nodes are linked together
// provided the nodes are appropriately placed in a memory block.
//
// The nodes responsibility is to know how their evaluation / propagation has
// to be performed (which function to call) and what argument(s) should be 
// used. The DAG does not know the actual type of the contained nodes. The
// DAG does not store any explicit information on how the nodes are linked
// together: it is implicitly ensured by the way the DAG is built.
//
// The evaluation of the expression x^2+y^2-1 with bound constraints 
// -1.0 <= x <= -0.5 and 0.5 <= y <= 1.0 is used to illustrate the above
// described idea. Once the .nl file obtained from the AMPL modeling
// environment is processed, a similar procedure is performed in the
// constructor of the dag.
//
// private members of dag
// var: variables
// num: numerical constants
// tmp: temporary values
//
// Building the DAG:
//
// T var[] = { T(-1.0, -0.5), T( 0.5,  1.0) };
//
// T num[] = { T(2.0), T(-1.0) };
//
// T tmp[4];
//
// node<T>** Node = new node<T>* [4];
//
// We are evaluating: var[0]^2 + var[1]^2 - 1
//
// Step 0: tmp[0] = var[0]^2
// Node[0] = new PowInt<T>(tmp  , var  , num  );
//
// Step 1: tmp[1] = var[1]^2
// Node[1] = new PowInt<T>(tmp+1, var+1, num  );
//
// Step 2: tmp[2] = tmp[0] + tmp[1]
// Node[2] = new Plus<T>  (tmp+2, tmp  , tmp+1);
//
// Step 3: tmp[3] = tmp[2] - 1
// Node[3] = new Plus<T>  (tmp+3, tmp+2, num+1);
//
// Now the DAG of the equation is built, we are ready to evaluate the
// expression with respect to the variable bounds. Evaluation is just a 
// for loop:
//
// for (int i=0; i<4; ++i)
//    Node[i]->evaluate();
//
// The value of the expression is:
//    Node[3]->value();
//
// The trick is to get the offsets right, for example 3, 2, 1 in the call
// in Step 3:
// new Plus<T>  (tmp+3, tmp+2, num+1);
//
// The major part of the .nl file processing is about calculating these
// offsets.
//
// Although in a very different context but with a similar idea, see 
// Item 23 in Meyers: Effective STL.
//
//-----------------------------------------------------------------------------

template <typename T>
class node {

	public:

		virtual void evaluate() = 0;

		virtual bool propagate() = 0;

		explicit node(T* const Val) : val(Val) { }

		const T value() const { return *val; }

		virtual ~node() { };

	protected:

		T* const val;

	private:

		node<T> (const node&);
		node<T>& operator=(const node& );

};

// Unary operators
//------------------

template <typename T>
class Ln : public node<T> {

public:

	Ln<T> (T* const Val, T* const Arg)
		  : node<T>(Val),     arg(Arg)  { };

	virtual void evaluate();

	virtual bool propagate();

	~Ln<T>() { } ;

private:

	T* const arg;
};

template <typename T>
class Exp : public node<T> {

public:

	Exp<T> (T* const Val, T* const Arg)
		   : node<T>(Val),     arg(Arg)  { };

	virtual void evaluate();

	virtual bool propagate();

	~Exp<T>() { } ;

private:

	T* const arg;
};

template <typename T>
class Neg : public node<T> {

public:

	Neg<T> (T* const Val, T* const Arg)
		   : node<T>(Val),     arg(Arg)  { };

	virtual void evaluate();

	virtual bool propagate();

	~Neg<T>() { } ;

private:

	T* const arg;
};

// Implementation of unary operators
//------------------------------------

template <typename T> 
void Ln<T>::evaluate() {
	
	*(this->val) = ln(*arg);
}

template <typename T> 
void Exp<T>::evaluate()  {

	*(this->val) = exp(*arg); 
}

template <typename T> 
void Neg<T>::evaluate()  {
	
	*(this->val) = -(*arg);
}

// Inverse operations of the unary operators
//-----------------------------------------------------------------------------
//
// You only need the propagate() functions if backward propagation is 
// applicable to your type and you wish to use it.
// See backward propagation in: H. Schichl, A. Neumaier; 
// Interval Analysis on Directed Acyclic Graphs for Global Optimization;
// J. Global Optimization, 2005, 33, 541�562; Section 5.
//
// If you are using an interval data type and its interface does not 
// match the interface of cxsc::interval then apply the adapter pattern
// (see: Gamma, Helm, Johnson, Vlissides: Design Patterns) to get the necessary
// functions.
// The operator& is used to compute the intersection of intervals as defined by
// the C-XSC developers; ext_div stands for EXTended interval DIVision followed
// by the intersection of the input interval and the resulting (extended)
// interval. The rest of the functions are hopefully self-explanatory, for the
// documentation of the cxsc::interval see:
//
// http://www.math.uni-wuppertal.de/~xsc/xsc/cxsc/apidoc/html/classcxsc_1_1interval.html
//
// Documentation of the type requirements to use the DAG class is undoubtedly
// needed. I will provide an example for the HansenType (C-XSC implementation
// of Hansen's generalized interval arithmetic) using the DAG class to perform
// function evaluations. But till then, the only support comes from the compile
// error messages, sorry... :(
//
// The propagate() should return true if the box cannot contain a solution,
// or in other words: if the input and the image interval are disjoint.
//
// It is far from trivial how backward propagation should be applied to the
// affine data type (i.e. what is the intersection of two linear polynomials
// with bound constraints on the variables?), backward propagation is not
// implemented for the affine data type. If backward propagation is not
// applicable to your type either (for example you are using Hansen's
// generalized intervals) or you simply do not wish to apply backward
// propagation then provide a dummy implementation for the necessary functions
// (by either throwing an exception or calling exit; you may get warnings due
// to unreachable code). For an example see affine.cpp at the end of that file.
//
//-----------------------------------------------------------------------------

template <typename T>
bool Ln<T>::propagate() {

	      T& x = *arg;
	const T& y = *(this->val);

	// Inverse operator
	const T x_new( exp(y) );
	if ( Disjoint(x, x_new) )
		return true;
	else {
		x = x & x_new;
		return false;
	}
}

template <typename T>
bool Exp<T>::propagate() {

	      T& x = *arg;
	const T& y = *(this->val);

	// Inverse operator
	const T x_new( ln(y) );
	if ( Disjoint(x, x_new) )
		return true;
	else {
		x = x & x_new;
		return false;
	}
}

template <typename T>
bool Neg<T>::propagate() {

	      T& x = *arg;
	const T& y = *(this->val);

	// Inverse operator
	const T x_new( -y );
	if ( Disjoint(x, x_new) )
		return true;
	else {
		x = x & x_new;
		return false;
	}
}

// Binary operators
//-------------------

template <typename T>
class Plus : public node<T> {

public:

	Plus<T> (T* const Val, T* const Arg1, T* const Arg2)
		    : node<T>(Val),    arg1(Arg1),    arg2(Arg2)  { };

	virtual void evaluate();

	virtual bool propagate();

	~Plus<T>() { } ;

private:

	T* const arg1;
	T* const arg2;
};

template <typename T>
class Minus : public node<T> {

public:

	Minus<T> (T* const Val, T* const Arg1, T* const Arg2)
		     : node<T>(Val),    arg1(Arg1),    arg2(Arg2)  { };

	virtual void evaluate();
	
	virtual bool propagate();

	~Minus<T>() { } ;

private:

	T* const arg1;
	T* const arg2;
};

template <typename T>
class Mult : public node<T> {

public:

	Mult<T> (T* const Val, T* const Arg1, T* const Arg2)
		    : node<T>(Val),    arg1(Arg1),    arg2(Arg2)  { };

	virtual void evaluate();

	virtual bool propagate();

	~Mult<T>() { } ;

private:

	T* const arg1;
	T* const arg2;
};

template <typename T>
class Div : public node<T> {

public:

	Div<T> (T* const Val, T* const Arg1, T* const Arg2)
		   : node<T>(Val),    arg1(Arg1),    arg2(Arg2)  { };

	virtual void evaluate();

	virtual bool propagate();

	~Div<T>() { } ;

private:

	T* const arg1;
	T* const arg2;
};

template <typename T>
class PowInt : public node<T> {

public:

	PowInt<T> (T* const Val, T* const Arg1, T* const Arg2)
		      : node<T>(Val),    arg1(Arg1),    arg2(Arg2)  { };

	virtual void evaluate();

	virtual bool propagate();

	~PowInt<T>() { } ;

private:

	T* const arg1;
	T* const arg2;
};

// Implementation of binary operators
//------------------------------------

template <typename T> 
void Plus<T>::evaluate()  {
	
	*(this->val) = ((*arg1)+(*arg2));
}

template <typename T> 
void Minus<T>::evaluate() {
	
	*(this->val) = ((*arg1)-(*arg2));
}

template <typename T> 
void Mult<T>::evaluate()  {
	
	*(this->val) = ((*arg1)*(*arg2));
}

template <typename T> 
void Div<T>::evaluate()   {
	
	*(this->val) = ((*arg1)/(*arg2));
}

//-----------------------------------------------------------------------------
// Only x^2 is currently implemented. The application exits with the 
// corresponding error message if PowInt is called with a different 
// argument than ^2. 
// What makes for example x^3 function tricky? If zero is contained in the 
// argument interval then the second derivative changes sign inside that
// interval. As a consequence, Theorem 2 of Stolfi and Figueiredo; 
// Self-Validated Numerical Methods and Applications (1997); p. 57. is not
// applicable. There are ways to solve this issue but it requires further
// studies and I am busy... :(

template <typename T> 
void PowInt<T>::evaluate()  {
	
	*(this->val) = sqr(*arg1);
}

//-----------------------------------------------------------------------------
// As with unary operators, you only need the propagate() functions if 
// backward propagation is applicable to your type and you wish to use it.

template <typename T>
bool Plus<T>::propagate() {

	      T& x = *arg1;
	      T& y = *arg2;
	const T& z = *(this->val);

	// Inverse operator
	const T x_new( z-y );
	if ( Disjoint(x, x_new) )
		return true;
	else {
		x = x & x_new;
	}

	const T y_new( z-x );
	if ( Disjoint(y, y_new) )
		return true;
	else {
		y = y & y_new;
	}
	
	return false;
}

template <typename T>
bool Minus<T>::propagate() {

	      T& x = *arg1;
	      T& y = *arg2;
	const T& z = *(this->val);

	// Inverse operator
	const T x_new( z+y );
	if ( Disjoint(x, x_new) )
		return true;
	else {
		x = x & x_new;
	}

	const T y_new( x-z );
	if ( Disjoint(y, y_new) )
		return true;
	else {
		y = y & y_new;
	}

	return false;
}

template <typename T>
bool Mult<T>::propagate() {

	      T& x = *arg1;
	      T& y = *arg2;
	const T& z = *(this->val);

	// Inverse operator
	if (ext_div(x, z, y))
		return true;
	
	return ext_div(y, z, x);
}

template <typename T>
bool Div<T>::propagate() {

	      T& x = *arg1;
	      T& y = *arg2;
	const T& z = *(this->val);

	// Inverse operator
	const T x_new( z*y );
	if ( Disjoint(x, x_new) )
		return true;
	else {
		x = x & x_new;
	}
	
	return ext_div(y, x, z);
}

template <typename T>
bool PowInt<T>::propagate() {

	      T& x = *arg1;
	const T& y = *(this->val);

	// Inverse operator
	T x_new( sqrt(y) );

	if      (Inf(x) >= 0.0)
		;
	else if (Sup(x) <= 0.0)
		x_new = T(-Sup(x_new), -Inf(x_new));
	else
		x_new = T(-Sup(x_new), Sup(x_new));

	if ( Disjoint(x, x_new) )
		return true;
	else {
		x = x & x_new;
	}
	
	return false;
}

// n-ary operators
//-------------------

template <typename T>
class Sum_cx : public node<T> {

public:

	Sum_cx<T> (T* const Val, T* Arg[], const int n_args);

	virtual void evaluate();

	virtual bool propagate();

	~Sum_cx<T>() ;

private:

	Sum_cx<T> (const Sum_cx&);
	Sum_cx<T>& operator=(const Sum_cx& );

	T** const arg;
	const int n;
};

template <typename T>
Sum_cx<T>::Sum_cx(T* const Val, T* Arg[], const int n_args )
	: node<T>(Val), arg(new T*[n_args]), n(n_args)
{
	assert ( (n_args%2)==0 );

	for (int i=0; i<n_args; ++i)
		arg[i] = Arg[i];
}

//-----------------------------------------------------------------------------
// Sum_cx has to take ownership of some of its pointer arguments. It is both
// hideous and error-prone. The code should be valgrind clean, though. This
// issue will (also) be resolved by refactoring the code and moving the messy
// .nl file processing to another class (or even to a standalone application).
// Handling the the .nl file should not be the responsability of the dag class.
// More on this below, in the next long comment block.

template <typename T>
Sum_cx<T>::~Sum_cx() {
	const int n2 = n/2;
	for (int i=0; i<n2; ++i)
		delete arg[i];
	delete[] arg;
}

template <typename T> 
void Sum_cx<T>::evaluate()
{
	const int n2 = n/2;
	T res( (*(arg[0])) * (*(arg[n2])) );
	for (int i=1; i<n2; ++i)
		res += (*(arg[i])) * (*(arg[i+n2]));
	*(this->val) = res;
}

template <typename T> 
bool Sum_cx<T>::propagate()
{
	const T& sum_cx = *(this->val);

	const int n2 = n/2;

	for (int j=0; j<n2; ++j) {

		T scx_j(0.0);

		for (int i=0; i<n2; ++i) {
			if (i==j)
				continue;
			scx_j += (*(arg[i])) * (*(arg[i+n2]));
		}

		const T x_j_new( (sum_cx - scx_j)/(*(arg[j])) );
		T& x_j = *(arg[j+n2]);

		if ( Disjoint(x_j, x_j_new) )
			return true;
		else {
			x_j = x_j & x_j_new;
		}
	}

	return false;
}

template <typename T>
class Sum_x : public node<T> {

public:

	Sum_x<T> (T* const Val, T* Arg[], const int n_args);

	virtual void evaluate();

	virtual bool propagate();

	~Sum_x<T>() ;

private:

	Sum_x<T> (const Sum_x&);
	Sum_x<T>& operator=(const Sum_x& );

	T** const arg;
	const int n;
};

template <typename T>
Sum_x<T>::Sum_x(T* const Val, T* Arg[], const int n_args )
	: node<T>(Val), arg(new T*[n_args]), n(n_args)
{
	for (int i=0; i<n_args; ++i)
		arg[i] = Arg[i];
}

template <typename T>
Sum_x<T>::~Sum_x() {

	delete[] arg;
}

template <typename T> 
void Sum_x<T>::evaluate()
{
	T res(*(arg[0]));
	for (int i=1; i<n; ++i)
		res += *(arg[i]);
	*(this->val) = res;
}

template <typename T> 
bool Sum_x<T>::propagate()
{
	const T& sum = *(this->val);

	for (int j=0; j<n; ++j) {

		T s_j(0.0);

		for (int i=0; i<n; ++i) {
			if (i==j)
				continue;
			s_j += *(arg[i]);
		}

		const T x_j_new(sum - s_j);
		T& x_j = *(arg[j]);

		if ( Disjoint(x_j, x_j_new) )
			return true;
		else {
			x_j = x_j & x_j_new;
		}
	}

	return false;
}

//=============================================================================
//
// IMPORTANT!!!
//
// Everything from this point to the "end of parsing" comment is subject to a
// COMPLETE REWRITE FROM SCRARTCH by applying the interpreter pattern (see: 
// Gamma, Helm, Johnson, Vlissides: Design Patterns). It will make the code 
// easier to understand, maintain and extend.
//
// As discussed in the first comment block at the top of this file, the major
// part of the .nl file processing is about calculating the offsets (see also 
// the example in the first comment block).
//
//-----------------------------------------------------------------------------
//
// The code is very messy in its current form, partly because the .nl file was
// never meant to be processed directly. The AMPL developers made publicly
// available the AMPL/Solver Interface Library which provides the necessary
// routines to handle this output. This makes it relatively easy to attach a
// solver to AMPL if the solver uses real arithmetic (arithmetic on type
// double). This is almost always the case. As a consequence, you are not
// supposed to read the .nl file directly and there is no need to know how the
// .nl file contains the necessary information. Quite unfortunately, ASL can
// only work with the data type double and ASL is incompatible with abstract
// data types.
//
// As it is now, the dag class is responsible for processing the .nl file which 
// is a clear violation of the one class - one responsability rule.
//
// Many headers are now included (sstream, vector, algorithm) that has nothing 
// to do with the client's code, unfortunately the dag.hpp dumps those headers
// on the user. This creates a bunch of problems that could be avoided 
// otherwise (see for example: minimizing compile time dependencies at
// Item 26-28 in Sutter: Exceptional C++; or Item 31 in Meyers: Effective C++,
// Third Edition).
//
// Moving the .nl file parsing and transformation into another class is needed
// but making a separate, standalone application seems a much better solution. 
// Other researchers may wish to use that application even if they are not 
// interested in the dag class or in interval methods.
//
//-----------------------------------------------------------------------------

// To produce debug info set it to true
const bool VERBOSE_OUTPUT(false);

struct counters {
	int n_vars;
	int n_nzeros;
	int n_dfvs;
	int n_node;
	int tmp_idx;
	int num_idx;
	int num_cnt;
};

// These functions are independent of type T and they are implemented in dag.cpp
void error(const char* msg);

bool is_segment(const std::string& str);

void read_file(const char* file_name, std::vector<std::string>& content, counters& s);

const std::string resolve_op(const std::string& op, int& n_args);

template <typename T>
T* set_arg(const std::string& str,
		   counters& s,
			T* tmp,
			T* var,
			T* num,
			int* dfv)
{
	int& n_vars = s.n_vars;
	int& n_dfvs = s.n_dfvs;
	int& n_node = s.n_node;
	//int& tmp_idx = s.tmp_idx;
	int& num_idx = s.num_idx;
	int& num_cnt = s.num_cnt;

	const char first_char = str.at(0);

	T* arg = 0;

	if (first_char == 'v') {
		int var_idx = atoi(str.substr(1).c_str());
		
		if (var_idx < n_vars) {
			arg = var+var_idx;
		}
		else {
			// defined variable
			assert(var_idx-n_vars < n_dfvs);
			arg = tmp + dfv[var_idx-n_vars];
			
		}
	}
	else if (first_char == 't') {
		int temp_idx = atoi(str.substr(1).c_str());
		assert (temp_idx < n_node);
		arg = tmp+temp_idx;
	}
	else if (first_char == 'n') {
		double num_val = atof(str.substr(1).c_str());
		assert(num_idx < num_cnt);
		num[num_idx]= T(num_val);
		arg = num+num_idx;
		++num_idx;
	}
	else
		error("parsing argument");

	return arg;
}

// An std::pair could do the same job...
struct ind_val {
	ind_val(const int i, const double v) : ind(i), val(v) {};
	int    ind;
	double val;
};

template <typename T>
std::string parse_nl(std::vector<std::string>& expr, counters& cnt, const int beg, node<T>** Node, T* var, T* tmp, T* num, int* dfv) {

	using std::vector;
	using std::string;
	using std::cout;
	using std::endl;
	using std::flush;

	//int& n_vars = cnt.n_vars;
	int& tmp_idx = cnt.tmp_idx;
	//int& num_idx = cnt.num_idx;

	int i = beg;

	vector<string> args;

	int n_args;

	while (expr.size()>1) {

		if( expr.at(i).at(0)!='o' )
			return expr.at(i);

		const string op = expr.at(i);

		const string res_op = resolve_op(op, n_args);

		int j, k;

		if (n_args) {
			j = i+1;
			k = i+n_args;
		}
		else {

			if (op != "o54")
				error("unimplemented operator found");
			
			n_args = atoi((expr.at(i+1)).c_str());

			if (n_args < 3)
				error("unexpected error parsing n-ary operator");

			for (int m=1; m<=n_args; ++m) {

				string arg;

				if (expr.at(i+1+m).at(0) != 'o')

					arg = expr.at(i+1+m);
				else

					arg = parse_nl<T>(expr, cnt, i+1+m, Node, var, tmp, num, dfv);

				args.push_back(arg);
			}
			k= i+1+n_args;
			j=k+1;
		}

		bool operator_resolved(true);

		for (int n=j; n<=k; ++n) {

			if (expr.at(n).at(0)=='o') {
				operator_resolved = false;
				// set i to point to the next operator
				// and start over
				// args must be cleared
				i  =  n;
				break;
			}
			
			args.push_back(expr.at(n));
		}

		if (operator_resolved) {
			
			if (VERBOSE_OUTPUT) {
				cout << "t" << tmp_idx << " = " << res_op << flush;
				for (int pos=0; pos<n_args; ++pos) {
					cout << " " << args.at(pos) << flush;
				}
				cout << endl;
			}

			std::ostringstream s;

			s << tmp_idx;

			if (n_args==1) {

				T* arg = set_arg(args.at(0), cnt, tmp, var, num, dfv);

				if ( op == "o16") {
					Node[tmp_idx] = new Neg<T>( tmp+tmp_idx, arg);
				}
				else if ( op == "o43") {
					Node[tmp_idx] = new Ln<T>( tmp+tmp_idx, arg);
				}
				else if ( op == "o44") {
					Node[tmp_idx] = new Exp<T>( tmp+tmp_idx, arg);
				}
				else
					error("implementation not updated properly");
			}
			else if (n_args==2) {

				T* arg1 = set_arg(args.at(0), cnt, tmp, var, num, dfv);
				T* arg2 = set_arg(args.at(1), cnt, tmp, var, num, dfv);

				if      ( op == "o0") {
					Node[tmp_idx] = new Plus<T>( tmp+tmp_idx, arg1, arg2);
				}
				else if ( op == "o1") {
					Node[tmp_idx] = new Minus<T>( tmp+tmp_idx, arg1, arg2);
				}
				else if ( op == "o2") {
					Node[tmp_idx] = new Mult<T>( tmp+tmp_idx, arg1, arg2);
				}
				else if ( op == "o3") {
					Node[tmp_idx] = new Div<T>( tmp+tmp_idx, arg1, arg2);
				}
				else if ( op == "o5") {
					const string ns = (args.at(1).substr(1));
					const double nf = atof(ns.c_str());
					const int ni    = static_cast<int>(nf);
					if (!((nf==ni)&&(ni==2)))
						error("sorry, only the square function is implemented");

					Node[tmp_idx] = new PowInt<T>( tmp+tmp_idx, arg1, arg2);
				}
				else
					error("implementation not updated properly");
			}
			else if (n_args > 2) {

				assert( op == "o54");
				
				T* * arg = new T* [n_args];

				for (int p=0; p<n_args; ++p)
					arg[p] = set_arg(args.at(p), cnt, tmp, var, num, dfv);

				Node[tmp_idx] = new Sum_x<T>( tmp+tmp_idx, arg, n_args);

				delete[] arg;
			}
			else
				error("unexpected error");

			++tmp_idx;

			expr.at(i) = string("t"+ s.str());

			// erase() is inefficient for large vectors and should be
			// avoided; this issue will be resolved by the rewrite;
			// it is NOT a perfomance bottleneck, vectors are tipically
			// tiny in this application
			expr.erase(expr.begin()+(i+1), expr.begin()+(k+1));

			i = beg;
			
		}

		args.clear();

	}

	return (expr.at(0));
}

template <typename T>
void parse_lin(std::vector<ind_val>& linp, counters& cnt, const bool is_linear, const double r, node<T>** Node, T* var, T* tmp, T* num, int* dfv) {

	using std::cout;
	using std::endl;
	using std::flush;

	int& n_vars = cnt.n_vars;
	//int& n_dfvs = cnt.n_dfvs;
	//int& n_node = cnt.n_node;
	int& tmp_idx = cnt.tmp_idx;
	//int& num_idx = cnt.num_idx;
	//int& num_cnt = cnt.num_cnt;

	int lin_args(linp.size());

	if (!is_linear)
		++lin_args;

	if (r!= T(0.0))
		++lin_args;

	T* * arg = new T* [2*lin_args];

	int n_args(0);

	//if ( lin_args == 2 )
	//	cout << "t" << tmp_idx << " = plus" << flush;
	//else
	if (VERBOSE_OUTPUT) {
		cout << "t" << tmp_idx << " = sum" << flush;
	}

	if (!is_linear) {

		if (VERBOSE_OUTPUT) {
			cout << " 1*t" << tmp_idx-1 << flush;
		}
		arg[n_args] = new T(1.0);

		std::ostringstream s;
		s << 't' << (tmp_idx-1);
		arg[n_args+lin_args] = set_arg<T>(s.str(), cnt, tmp, var, num, dfv);

		++n_args;
	}

	for (int i=0; i<linp.size(); ++i) {

		ind_val p = linp.at(i);

		assert( p.val != 0.0 );
		assert( p.ind < n_vars);
		
		if (VERBOSE_OUTPUT) {
			cout << " c" << p.ind << "*v" << p.ind << flush;
		}

		arg[n_args] = new T(p.val);

		std::ostringstream s;
		s << 'v' << p.ind;
		arg[n_args+lin_args] = set_arg(s.str(), cnt, tmp, var, num, dfv);

		++n_args;
	}

	if ( r!= T(0.0) ) {

		if (VERBOSE_OUTPUT) {
			cout << " (-1)*r" << flush;
		}

		arg[n_args] = new T(-1.0);
		std::ostringstream s;
		s << 'n' << r ;
		arg[n_args+lin_args] = set_arg(s.str(), cnt, tmp, var, num, dfv);
		++n_args;
	}

	if (VERBOSE_OUTPUT) {
		cout << endl;
	}

	assert(n_args == lin_args);

	Node[tmp_idx] = new Sum_cx<T>( tmp+tmp_idx, arg, 2*lin_args);

	delete[] arg;

	++tmp_idx;
}

template <typename T>
void parse_V_C(std::vector<std::string>& expr, counters& s, std::vector<ind_val>& linp, const double r, node<T>** Node, T* var, T* tmp, T* num, int* dfv) {

	bool is_linear(false);

	if (expr.size() == 1)
		is_linear = true;
	else
		parse_nl<T>(expr, s, 0, Node, var, tmp, num, dfv);

	// ???
	if (linp.size() != 0 || r!=0.0) {
		parse_lin<T>(linp, s, is_linear, r, Node, var, tmp, num, dfv);

	}
}

template <typename T>
void parse_exprs(std::vector<std::string>& content, counters& cnt, node<T>** Node, T* var, T* tmp, T* num, int* dfv, int* con) {

	using std::vector;
	using std::string;
	using std::cout;
	using std::endl;

	int& n_vars = cnt.n_vars;
	int& n_dfvs = cnt.n_dfvs;
	int& tmp_idx = cnt.tmp_idx;

	while (!content.empty()) {

		vector<string>::iterator segm_beg = content.begin();
		
		assert (segm_beg != content.end());

		vector<string>::iterator segm_end = std::find_if(segm_beg+1, content.end(), is_segment);

		string s(*segm_beg);

		const char c = s.at(0);

		if (c == 'V') {

			const string::size_type a = s.find(" ");
			const string::size_type b = s.find(" ", a+1);

			const string idx(s.substr(1, a-1));
			const string lin(s.substr(a+1, b-a-1));

			const int v_idx = atoi( idx.c_str() );
			const int l_prt = atoi( lin.c_str() );

			vector<ind_val> linp;

			int i=1;

			while (i<=l_prt) {

				const string e(*(segm_beg+i));
				const string::size_type spc = e.find(" ");

				const string idx(e.substr(0, spc));
				const string val(e.substr(spc+1));

				const int    index = atoi(idx.c_str());
				const double value = atof(val.c_str());

				ind_val p(index, value);

				linp.push_back(p);

				++i;
			}

			vector<string> expr(segm_beg+i, segm_end);

			if (VERBOSE_OUTPUT) {
				cout << "Parsing defined variable #" << v_idx << endl;
			}

			double r = 0.0;
			if (expr.size() == 1) {
				r = -atof(((expr.at(0)).substr(1)).c_str());
			}

			parse_V_C<T>(expr, cnt, linp, r, Node, var, tmp, num, dfv);

			if (VERBOSE_OUTPUT) {
				cout << "Indices stored for this variable are: " << v_idx-n_vars << "   " << tmp_idx-1 << endl;
			}

			assert( v_idx-n_vars < n_dfvs );

			dfv[v_idx-n_vars] = tmp_idx-1;
			
			if (VERBOSE_OUTPUT) {
				cout << endl;			
			}
		}
		else if (c == 'C') {

			// nonlinear part
			//---------------------------------------
			const string idx(s.substr(1));

			const int c_idx = atoi( idx.c_str() );

			vector<string> expr(segm_beg+1, segm_end);


			// linear part
			//---------------------------------------

			string lin_name("J"+idx);

			const int off = segm_end-segm_beg-1;
			int pos = 1;
			while (true) {
				const string line = content.at(off+pos++);
				if (line.at(0) == 'J') {
					string curr( line.substr(0, line.find(' ')) );
					if (curr == lin_name) {
						break;
					}
				}
			}

			vector<string>::iterator j_beg = segm_beg+(off+pos-1);

			assert( j_beg != content.end() );

			vector<string>::iterator j_end = find_if(j_beg+1, content.end(), is_segment);

			const string j(*j_beg);

			const string::size_type a = j.find(" ");

			const string lin(j.substr(a+1));

			const int l_prt = atoi( lin.c_str() );

			vector<ind_val> linp;

			int i=1;

			while (i<=l_prt) {

				const string e(*(j_beg+i));
				const string::size_type spc = e.find(" ");

				const string idx(e.substr(0, spc));
				const string val(e.substr(spc+1));

				const int    index = atoi(idx.c_str());
				const double value = atof(val.c_str());

				if (value != 0.0) {

					ind_val p(index, value);
					linp.push_back(p);
				}

				++i;
			}

			// rhs of constraint
			//-------------------------------------------------

			vector<string>::iterator r_beg = find(content.begin(), content.end(), "r");

			assert (r_beg != content.end());

			string rhs(*(r_beg+1+c_idx));

			assert (rhs.at(0) == '4');

			const double r_val = atof((rhs.substr(2)).c_str());

			//const T r(r_val);

			//------------------------------------------------

			if (VERBOSE_OUTPUT) {
				cout << "Parsing constraint #" << c_idx << endl;
			}

			parse_V_C<T>(expr, cnt, linp, r_val, Node, var, tmp, num, dfv);

			if (VERBOSE_OUTPUT) {
				cout << "Indices stored for this constraint are: " << c_idx << "   " << tmp_idx-1 << endl;
			}

			assert( c_idx < n_vars );

			con[c_idx] = tmp_idx-1;

			if (VERBOSE_OUTPUT) {
				cout << endl;	
			}
			//content.erase(j_beg, j_end);


		}
		else if (c == 'b') {

			assert(n_vars);
			
			for (int i=1; i<=n_vars; ++i) {
				const string line(content.at(i));
				if (line.at(0) != '0')
					error("variables should have proper lower and upper bounds");

				string::size_type pos = line.find(" ", 3);

				const string val_lo ( line.substr(2, pos-2) );
				const string val_up ( line.substr(pos+1)    );

				var[i-1] = T(atof(val_lo.c_str()), atof(val_up.c_str()));
				//var[i-1] = T(atof(val_up.c_str()));
			}
		}
		else {
			if (VERBOSE_OUTPUT) {
				cout << "Found segment " << c << ", skipping..." << endl;
			}
		}

		// erase() is inefficient for large vectors and should be
		// avoided; this issue will be resolved by the rewrite;
		// it is NOT a perfomance bottleneck, vectors are tipically
		// tiny in this application
		content.erase(segm_beg, segm_end);

	}
}

//
// end of parsing and processing the .nl file
//=============================================================================
//
// After the rewriting the above bloody .nl file parsing and transformation
// (the new version will be based on the interpreter pattern), the dag 
// constructor will become cleaner too.
//

template <typename T>
dag<T>::dag(const char* file_name) {

	using std::vector;
	using std::string;

	vector<string> content;
	counters s;

	read_file(file_name, content, s);

	int& n_vars = s.n_vars;
	int& n_nzeros = s.n_nzeros;
	int& n_dfvs = s.n_dfvs;
	int& n_node = s.n_node;
	int& tmp_idx = s.tmp_idx;
	//int& num_idx = s.num_idx;
	int& num_cnt = s.num_cnt;

	int op_cnt(0);
	int nu_cnt(0);

	for (vector<string>::iterator i(content.begin()); i!=content.end(); ++i) {

		const char c = (*i).at(0);

		if (c == 'o')
			++op_cnt;
		else if (c == 'n')
			++nu_cnt;
	}

	nu_cnt += n_vars;
	n_node = op_cnt+n_vars+n_dfvs;
	num_cnt = nu_cnt;

	//------------------------------

	Node = new node<T>*[n_node];

	for (int i=0; i<n_node; ++i)
		Node[i] = 0;
	
	var = new T[n_vars];
	tmp = new T[n_node];
	num = new T[nu_cnt];
	
	dfv = new int[n_dfvs];
	con = new int[n_vars];

	num_of_vars = n_vars;
	num_of_nzeros = n_nzeros;
	num_of_nodes = 0;

	//---------------------------

	parse_exprs<T>(content, s, Node, var, tmp, num, dfv, con);

	assert (tmp_idx <= n_node);
	
	num_of_nodes = tmp_idx;

}

template <typename T>
dag<T>::~dag() {

	for (int i=0; i<num_of_nodes; ++i) {
		delete Node[i];
		Node[i] = 0;
	}

	delete[] Node;
	Node = 0;
	delete[] var;
	var = 0;
	delete[] tmp;
	tmp = 0;
	delete[] num;
	num = 0;

	delete[] dfv;
	dfv = 0;
	delete[] con;
	con = 0;
}

template <typename T>
void dag<T>::evaluate_all() {

	if (!num_of_nodes)
		error("number of actual nodes is not set");	
	
	if (VERBOSE_OUTPUT)
		std::cout << std::endl;

	for (int i=0; i<num_of_nodes; ++i) {
		Node[i]->evaluate();

		if (VERBOSE_OUTPUT)
			std::cout << "tmp[" << i << "] = " << Node[i]->value() << std::endl;
	}
}

template <typename T>
bool dag<T>::propagate_all() {

	if (VERBOSE_OUTPUT)
		std::cout << std::endl;

	for (int i=num_of_nodes-1; i>=0; --i) {
		const bool toDelete = Node[i]->propagate();
		if (toDelete)
			return true;
		if (VERBOSE_OUTPUT)
			std::cout << "tmp[" << i << "] = " << Node[i]->value() << std::endl;
	}

	return false;
}

template <typename T>
void dag<T>::set_variables (const T vars[]) {

	const int n = num_of_vars;
	
	for (int i=0; i<n; ++i)
		var[i] = vars[i];
}

template <typename T>
bool dag<T>::get_constraints(T r[]) {

	evaluate_all();
	
	const int n = num_of_vars;

	bool to_delete = false;

	for (int i=0; i<n; ++i) {

		r[i] = tmp[con[i]];

		if (!in(0.0, r[i]))
			to_delete = true;

		if (VERBOSE_OUTPUT) {
			std::cout << "constraint #" << i << " = " << std::endl << tmp[con[i]] << std::endl;
		}
	}

	return to_delete;
}

template <typename T>
bool dag<T>::propagate() {

	evaluate_all();
	
	const int n = num_of_vars;

	for (int i=0; i<n; ++i) {
		if (!(in(0.0, tmp[con[i]])))
			return true;
		else
			tmp[con[i]] = T(0.0);
	}

	const bool toDelete = propagate_all();

	return toDelete;
}

}

#endif
