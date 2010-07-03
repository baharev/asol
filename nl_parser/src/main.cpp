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
#include <fstream>
#include "parser.hpp"
#include "utility.hpp"
#include "interpreter_C.hpp"
// FIXME Keep only necessary includes
#include "interpreter_DAG.hpp"

using namespace std;

void read_nl_file_list(const char* infile, vector<string>& tests) {

	ifstream in(infile);

	if (!in.good())
		error("cannot read input from file");

	string line;

	do {
		line.clear();
		getline(in, line);

		if (line.length() == 0)
			break;

		const string problem_name = chop_extension_nl(line);

		tests.push_back(problem_name);

	} while ( !in.eof() );

	in.clear();
	in.close();
}

void generate_source(const vector<string>& tests) {

	const int n = static_cast<int> (tests.size());

	for (int i=0; i<n; ++i) {

		const string problem_name = tests.at(i);

		const string nl_file_name = problem_name + ".nl";

		parser p(nl_file_name.c_str());

		const string outfile_name = problem_name + ".cpp";

		fstream out(outfile_name.c_str());

		if (!out.good())
			error("failed to create output file");

		interpreter_C ipt(p.get_problem_representation(), out);
	}
}

void generate_runner(const vector<string>& tests) {

	ofstream out("main.cpp");

	if (!out.good())
		error("failed to create output file");

	const int n = static_cast<int> (tests.size());

	for (int i=0; i<n; ++i)
		out << "bool " << tests.at(i) << "();" << endl;

	out << endl;
	out << "int main() {" << endl;
	for (int i=0; i<n; ++i)
		out << "  " << tests.at(i) << "();" << endl;
	out << "  return 0;" << endl;
	out << "}" << endl << endl;
}

void write_DAG(const vector<string>& tests) {

	const int n = static_cast<int> (tests.size());

	for (int i=0; i<n; ++i) {

		const string problem_name = tests.at(i);

		const string nl_file_name = problem_name + ".nl";

		parser p(nl_file_name.c_str());

		const string outfile_name = problem_name + ".dag";

		fstream out(outfile_name.c_str(), ios_base::out);

		if (!out.good())
			error("failed to create output file");

		interpreter_DAG ipt(p.get_problem_representation(), out);
	}
}

int main(int argc, char** argv) {

	if (argc != 2)
		error("incorrect number of arguments");

	vector<string> tests;

	read_nl_file_list(argv[1], tests);

	//generate_source(tests);

	//generate_runner(tests);

	write_DAG(tests);

	return 0;
}
