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

#include <iostream>
#include <sstream>
#include "utility.hpp"
using namespace std;

namespace {

void print(string& msg, const char* message) {
	msg += message;
	msg += "!";
	cerr << endl << msg << endl;
}

}

void error(const char* message) {
	string msg("Error: ");
	print(msg, message);
	throw msg.c_str();
}

void warning(const char* message) {
	string msg("Warning: ");
	print(msg, message);
}

const string int_to_string(const int i)  {
	ostringstream s;
	s << i << flush;
	return s.str();
}

const string double_to_string(const double d)  {
	ostringstream s;
	s << d << flush;
	return s.str();
}

const string chop_extension_nl(const string& file_name) {

	const string::size_type length = file_name.size();

	if (length < 4)
		error("file name seems to be invalid");

	const string extension = file_name.substr(length-3);

	if (extension != ".nl")
		error("extension should be .nl");

	const string problem_name = file_name.substr(0, length-3);

	return problem_name;
}

