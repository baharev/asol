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

#ifndef DAG_HPP_
#define DAG_HPP_

#include <iosfwd>
#include <sstream>
#include <string>
#include "dag_reader.hpp"
#include "utility.hpp"

template <typename T>
class node;

template <typename T>
class dag {

	public:

		explicit dag(const char* file_name);

		void set_variables (const T var[]);

		const T* get_variables() const { return var; }

		const T* get_constraints(bool& to_delete);

		bool propagate();

		int number_of_vars() const { return num_of_vars; }

		int number_of_nzeros() const { return num_of_nzeros; }

		~dag();

	private:

		dag<T> (const dag<T>&);
		dag<T>& operator=(const dag<T>& );

		void evaluate_all();
		bool propagate_all();

		T* get_arg(std::stringstream& in, const int index);
		void add_unary_primitive(const std::string& op, int i, T* arg);
		void add_binary_primitive(const std::string& op, int i, T* arg1, T* arg2);

		node<T>** Node;
		T* var;
		T* tmp;
		T* num;

		// FIXME dfv is unused
		//int* dfv;
		T** con;
		T* constraint_buffer;

		int num_of_vars;
		int num_of_nzeros;
		int num_of_prims;
		int num_of_cons;
		int num_of_nums;
};

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

// FIXME Make them inline!

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

template <typename T>
bool Ln<T>::propagate() {

	      T& x = *arg;
	const T& y = *(this->val);

	// Inverse operator
	const T x_new( exp(y) );
	if ( Disjoint(x, x_new) )
		return true;
	else {
		x = Intersect(x, x_new);
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
		x = Intersect(x, x_new);
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
		x = Intersect(x, x_new);
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
		x = Intersect(x, x_new);
	}

	const T y_new( z-x );
	if ( Disjoint(y, y_new) )
		return true;
	else {
		y = Intersect(y, y_new);
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
		x = Intersect(x, x_new);
	}

	const T y_new( x-z );
	if ( Disjoint(y, y_new) )
		return true;
	else {
		y = Intersect(y, y_new);
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
		x = Intersect(x, x_new);
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
		x = Intersect(x, x_new);
	}

	return false;
}

void skip_text(std::stringstream& in, const char* content, const bool check_content = true) {

	std::string s("");

	// TODO Understand why...
	do {
		getline(in, s);

		if (!in.good())
			error("problems with reading from the stream");

	} while (s=="");

	if ( check_content && (s != content) ) {
		error("unexpected content");
	}
}

void check_index(int i, int ub) {

	if ((i<0) || (i>=ub)) {
		// FIXME Define dbg() with Bug:
		error("offset is out of range");
	}
}

template <typename T>
T* dag<T>::get_arg(std::stringstream& in, const int i) {

	char a;
	in >> a;

	int offset;
	in >> offset;

	T* arg = 0;
	int ub = 0;

	if (a == 'n') {
		ub  = num_of_nums;
		arg = num;
	}
	else if (a == 'v') {
		ub  = num_of_vars;
		arg = var;
	}
	else if (a == 't') {
		ub  = i;
		arg = tmp;
	}
	else {
		error("incorrect argument");
	}

	check_index(offset, ub);

	arg += offset;

	return arg;
}

template <typename T>
void dag<T>::add_unary_primitive(const std::string& op, int i, T* arg) {

	// *t = f(*arg);
	T* t = tmp + i;

	if      (op == "exp") {
		Node[i] = new Exp<T> (t, arg);
	}
	else if (op == "log") {
		Node[i] = new Ln<T>  (t, arg);
	}
	else if (op == "-"  ) {
		Node[i] = new Neg<T> (t, arg);
	}
	else {
		error("implementation not updated properly");
	}
}

template <typename T>
void dag<T>::add_binary_primitive(const std::string& op, int i, T* arg1, T* arg2) {

	// *t = f(*arg1, *arg2);
	T* t = tmp + i;

	if     (op == "+") {
		Node[i] = new Plus<T>  (t, arg1, arg2);
	}
	else if (op == "-") {
		Node[i] = new Minus<T> (t, arg1, arg2);
	}
	else if (op == "*") {
		Node[i] = new Mult<T>  (t, arg1, arg2);
	}
	else if (op == "/") {
		Node[i] = new Div<T>   (t, arg1, arg2);
	}
	else if (op == "pow") {
		// TODO Only the square function is implemented
		if (*arg2 != T(2.0))
			error("sorry, only the square function is implemented");
		Node[i] = new PowInt<T>(t, arg1, arg2);
	}
	else {
		error("implementation not updated properly");
	}
}

template <typename T>
dag<T>::dag(const char* file_name) {

	std::stringstream in;

	grab_content(file_name, in);

	skip_text(in, "problem name", false);

	skip_text(in, "number of primitives");
	in >> num_of_prims;

	skip_text(in, "number of variables");
	in >> num_of_vars;

	skip_text(in, "number of constraints");
	in >> num_of_cons;

	skip_text(in, "number of nonzeros");
	in >> num_of_nzeros;

	skip_text(in, "number of numeric constants");
	in >> num_of_nums;

	// FIXME Error checking >=0

	//--------------------------------------------------------------------------

	Node = new node<T>*[num_of_prims];

	for (int i=0; i<num_of_prims; ++i)
		Node[i] = 0;

	var = new T[num_of_vars];
	tmp = new T[num_of_prims];
	num = new T[num_of_nums];

	// FIXME Defined variables are unused
	//dfv = new int[n_dfvs];

	con = new T*[num_of_cons];
	constraint_buffer = new T[num_of_cons];

	//--------------------------------------------------------------------------

	skip_text(in, "numeric constant values");

	for (int i=0; i<num_of_nums; ++i) {
		double value;
		in >> value;
		num[i] = T(value);
	}

	//--------------------------------------------------------------------------

	skip_text(in, "constraint indices in primitives");

	for (int i=0; i<num_of_cons; ++i) {
		int index;
		in >> index;
		check_index(index, num_of_prims);
		con[i] = tmp + index;
	}

	//--------------------------------------------------------------------------

	skip_text(in, "primitives");

	for (int i=0; i<num_of_prims; ++i) {

		int arity;
		in >> arity;

		std::string op;
		in >> op;

		int index;
		in >> index;

		if (index != i)
			error("unexpected error");

		if (arity == 1) {
			T* arg = get_arg(in, i);
			add_unary_primitive(op, i, arg);
		}
		else if (arity == 2) {
			T* arg1 = get_arg(in, i);
			T* arg2 = get_arg(in, i);
			add_binary_primitive(op, i, arg1, arg2);
		}
		else {
			error("incorrect arity");
		}
	}

	//--------------------------------------------------------------------------

	// FIXME Make it optional!

	skip_text(in, "initial point");

	for (int i=0; i<num_of_vars; ++i) {
		double value;
		in >> value;
		var[i] = T(value);
	}

	//--------------------------------------------------------------------------

	skip_text(in, "EOF");
}

template <typename T>
dag<T>::~dag() {

	for (int i=0; i<num_of_prims; ++i) {
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

	// FIXME Update if necessary
	//delete[] dfv;
	//dfv = 0;
	delete[] con;
	con = 0;

	delete[] constraint_buffer;
	constraint_buffer = 0;
}

template <typename T>
void dag<T>::evaluate_all() {

	for (int i=0; i<num_of_prims; ++i)
		Node[i]->evaluate();
}

template <typename T>
bool dag<T>::propagate_all() {

	for (int i=num_of_prims-1; i>=0; --i) {
		const bool toDelete = Node[i]->propagate();
		if (toDelete)
			return true;
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
const T* dag<T>::get_constraints(bool& to_delete) {

	evaluate_all();

	const int n = num_of_cons;

	to_delete = false;

	for (int i=0; i<n; ++i) {

		constraint_buffer[i] = *(con[i]);

		if (!in(0.0, constraint_buffer[i]))
			to_delete = true;
	}

	return constraint_buffer;
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

#endif
