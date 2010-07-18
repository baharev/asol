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

#include <iostream>
#include "affine.hpp"

using namespace std;
using namespace aa_lib;

namespace {

const affine R = affine(1.9858775);

const affine V[] = { affine(0.0), affine(58.39), affine(89.57), affine(18.05)};

const affine k[][4] = {   { affine(0.0), affine(0.0), affine(0.0), affine(0.0)},
		            { affine(0.0), affine(0.0), affine(694.0825), affine(393.1971)},
		            { affine(0.0), affine(-149.7978), affine(0.0), affine(6811.3433)},
		            { affine(0.0), affine(926.2630), affine(1888.8509), affine(0.0)},  };


affine Lambda[4][4];

void compute_Lambda(const affine& T) {

	for (int i=1; i<=3; ++i) {
		for (int j=1; j<=3; ++j) {
			Lambda[i][j] = V[j]/V[i]*exp(-k[i][j]/(R*T));
			cout << "Lambda[" << i << "][" << j << "]: " << Lambda[i][j] << endl;
		}
	}
}

}

int main() {

	affine T(330.0, 380.0);
	affine x[] = { affine(0.0), affine(0.0,1.0), affine(0.0,1.0), affine(0.0,1.0) };

	compute_Lambda(T);
	return 0;
}
