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
#include <cassert>
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

const affine one(1.0);

const affine zero(0.0);

void compute_Lambda(const affine& T) {

	affine rT = reciprocal(R*T);

	for (int i=1; i<=3; ++i) {
		for (int j=1; j<=3; ++j) {
			Lambda[i][j] = V[j]/V[i]*exp(-k[i][j]*rT);
			cout << "Lambda[" << i << "][" << j << "]: " << endl;
			cout <<  Lambda[i][j] << endl;
		}
	}
}

void compute_s(const affine x[], affine s[]) {

	for (int i=1; i<=3; ++i)
		cout << endl << "x[" << i << "]" << endl << x[i] << endl;

	for (int i=1; i<=3; ++i) {

		s[i] = affine(0.0);

		for (int j=1; j<=3; ++j) {

			s[i] = s[i] + x[j]*Lambda[i][j];
		}
	}
}

void run() {

	affine::set_max_used_index(0);

	affine x[] = { affine(0.0), affine(0.0,1.0), affine(0.0,1.0), affine(0.0,1.0) };

	affine T(330.0, 380.0);

	affine s[] = { affine(0.0), affine(0.0,1.0), affine(0.0,1.0), affine(0.0,1.0) };

	compute_Lambda(T);

	affine sg[4];

	compute_s(x, sg);

	for (int i=1; i<=3; ++i) {
		const bool ok = s[i].intersect_domain(sg[i].inf(), sg[i].sup());
		assert(ok);
		//cout << endl << "s[" << i << "]" << endl << s[i] << endl;
	}

	affine xg[4];

	for (int j=1; j<=3; ++j) {

		for (int i=1; i<=3; ++i)
			xg[i] = x[i];

		xg[j] = one;

		for (int i=1; i<=3; ++i) {
			if (i==j)
				continue;
			xg[j] = xg[j] - x[i];
		}

		xg[j].set_range(x[j].inf(), x[j].sup());

		compute_s(xg, sg);

		for (int i=1; i<=3; ++i) {
			const bool ok = s[i].intersect_domain(sg[i].inf(), sg[i].sup());
			assert(ok);
			cout << endl << "s[" << i << "]" << endl << s[i] << endl;
		}
	}
}

}

void example___() {

	affine::set_max_used_index(0);

	affine x[] = { affine(0.0), affine(0.0,1.0), affine(0.0,1.0), affine(0.0,1.0) };

	affine T(330.0, 380.0);

	compute_Lambda(T);

	x[1].inf();

}

void example_1() {

	affine::set_max_used_index(0);

	affine x(2.0, 4.0);

	affine y = (x-affine(1.0))/(sqr(x)+affine(2.0));

	cout << endl << "y = (x-1)/(x^2+2)" << endl << y << endl;

}

void example_2() {

	affine::set_max_used_index(0);

	affine x(1.0, 2.0);

	affine y = (sqr(x)+x)/(affine(16.0)*x-affine(9.0));

	cout << endl << "y = (x^2+x)/(16*x-9)" << endl << y << endl;
}

void example_3() {

	affine::set_max_used_index(0);

	const affine x(1.0, 10.0);
	const affine y(1.0, 10.0);

	const affine xy = x*y;

	const affine numerator   = affine(5.0)*x-affine(4.0)*sqr(y)+affine(14.0)*xy;
	const affine denominator = sqr(x)+y+xy;

	cout << endl << "xy" << endl << xy << endl;
	cout << endl << "5x-4y^2+14xy" << endl << numerator << endl;
	cout << endl << "x^2+y+xy" << endl << denominator << endl;

	affine z = numerator/denominator;

	cout << endl << "z = (5x-4y^2+14xy)/(x^2+y+xy)" << endl << z << endl;
}

void example_4() {

	affine::set_max_used_index(0);

	affine one(1.0);

	affine x[] = {
			affine(0.0),
			affine(0.6, 0.7),
			affine(0.2, 0.3),
			affine(0.1, 0.2)
	};

	affine y[] = {
			affine(0.0),
			affine(0.1, 0.4),
			affine(0.2, 0.4),
			affine(0.3, 0.5)
	};

	affine z[7];

	z[0] = x[1]*y[1] + x[1]*y[2] + x[2]*y[2] + x[3]*y[1];

	z[1] = (x[1] + x[3])*y[1] + (x[1] + x[2])*y[2];

	z[2] = (one - x[2])*y[1] + (x[1] + x[2])*y[2];

	z[3] = (x[1] + x[3])*y[1] + (one - x[3])*y[2];

	z[4] = (one - x[2])*y[1] + (one - x[3])*y[2];

	z[5] = x[1]*(y[1] + y[2]) + x[2]*y[2] + x[3]*y[1];

	z[6] = x[1]*(one - y[3]) + x[2]*y[2] + x[3]*y[1];

	cout << endl << "z = x1*y1 + x1*y2 + x2*y2 + x3*y1" << endl << endl;

	for (int i=0; i<7; ++i) {
		cout << "n: " << z[i].n_elem() << endl;
		cout << z[i] << endl;
	}

}

void example_5() {

	affine::set_max_used_index(0);

	affine one(1.0);

	affine x[] = {
			affine(0.0),
			affine(0.6, 0.7),
			affine(0.2, 0.3),
			affine(0.1, 0.2)
	};

	affine y[] = {
			affine(0.0),
			affine(0.1, 0.4),
			affine(0.2, 0.4),
			affine(0.3, 0.5)
	};

	affine x1, x2, x3, y1, y2, y3;

	y1 = y[1];
	y2 = y[2];
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = x[3];

	affine z = x1*y1 + x1*y2 + x2*y2 + x3*y1;

	cout << "=========================================================" << endl;
	cout << "Procedure 2" << endl;
	cout << "Initial z" << endl << z << endl;

	cout << "---" << endl;
	cout << "y1" << endl;

	y1 = one - y[2] - y[3];
	y1.set_range(y[1].inf(), y[1].sup());
	y2 = y[2];
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = x1*y1 + x1*y2 + x2*y2 + x3*y1;

	cout << "z (x3 y1): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = x1*y1 + x1*y2 + x2*y2 + x3*y1;

	cout << "z (x2 y1): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = x1*y1 + x1*y2 + x2*y2 + x3*y1;

	cout << "z (x1 y1): " << endl << z << endl;

	cout << "---" << endl;
	cout << "y2" << endl;

	y1 = y[1];
	y2 = one - y[1] - y[3];
	y2.set_range(y[2].inf(), y[2].sup());
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = x1*y1 + x1*y2 + x2*y2 + x3*y1;

	cout << "z (x3 y2): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = x1*y1 + x1*y2 + x2*y2 + x3*y1;

	cout << "z (x2 y2): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = x1*y1 + x1*y2 + x2*y2 + x3*y1;

	cout << "z (x1 y2): " << endl << z << endl;

	cout << "---" << endl;
	cout << "y3" << endl;

	y1 = y[1];
	y2 = y[2];
	y3 = one - y[1] - y[2];
	y3.set_range(y[3].inf(), y[3].sup());

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = x1*y1 + x1*y2 + x2*y2 + x3*y1;

	cout << "z (x3 y3): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = x1*y1 + x1*y2 + x2*y2 + x3*y1;

	cout << "z (x2 y3): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = x1*y1 + x1*y2 + x2*y2 + x3*y1;

	cout << "z (x1 y3): " << endl << z << endl;

}

void example_6() {

	affine::set_max_used_index(0);

	affine one(1.0);

	affine x[] = {
			affine(0.0),
			affine(0.6, 0.7),
			affine(0.2, 0.3),
			affine(0.1, 0.2)
	};

	affine y[] = {
			affine(0.0),
			affine(0.1, 0.4),
			affine(0.2, 0.4),
			affine(0.3, 0.5)
	};

	affine x1, x2, x3, y1, y2, y3;

	y1 = y[1];
	y2 = y[2];
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = x[3];

	affine z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "=========================================================" << endl;
	cout << "Procedure 2" << endl;
	cout << "Initial z" << endl << z << endl;

	cout << "---" << endl;
	cout << "y1" << endl;

	y1 = one - y[2] - y[3];
	y1.set_range(y[1].inf(), y[1].sup());
	y2 = y[2];
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x3 y1): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x2 y1): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x1 y1): " << endl << z << endl;

	cout << "---" << endl;
	cout << "y2" << endl;

	y1 = y[1];
	y2 = one - y[1] - y[3];
	y2.set_range(y[2].inf(), y[2].sup());
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x3 y2): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x2 y2): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x1 y2): " << endl << z << endl;

	cout << "---" << endl;
	cout << "y3" << endl;

	y1 = y[1];
	y2 = y[2];
	y3 = one - y[1] - y[2];
	y3.set_range(y[3].inf(), y[3].sup());

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x3 y3): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x2 y3): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x1 y3): " << endl << z << endl;

}

void example_7() {

	affine::set_max_used_index(0);

	affine one(1.0);

	affine x[] = {
			affine(0.0),
			affine(0.6, 0.7),
			affine(0.2, 0.3),
			affine(0.1, 0.2)
	};

	affine y[] = {
			affine(0.0),
			affine(0.1, 0.4),
			affine(0.2, 0.4),
			affine(0.3, 0.5)
	};

	affine x1, x2, x3, y1, y2, y3;

	y1 = y[1];
	y2 = y[2];
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = x[3];

	affine z = (one-x2)*y1+(one-x3)*y2;

	cout << "=========================================================" << endl;
	cout << "Procedure 2" << endl;
	cout << "Initial z" << endl << z << endl;

	cout << "---" << endl;
	cout << "y1" << endl;

	y1 = one - y[2] - y[3];
	y1.set_range(y[1].inf(), y[1].sup());
	y2 = y[2];
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = (one-x2)*y1+(one-x3)*y2;

	cout << "z (x3 y1): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = (one-x2)*y1+(one-x3)*y2;

	cout << "z (x2 y1): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = (one-x2)*y1+(one-x3)*y2;

	cout << "z (x1 y1): " << endl << z << endl;

	cout << "---" << endl;
	cout << "y2" << endl;

	y1 = y[1];
	y2 = one - y[1] - y[3];
	y2.set_range(y[2].inf(), y[2].sup());
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = (one-x2)*y1+(one-x3)*y2;

	cout << "z (x3 y2): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = (one-x2)*y1+(one-x3)*y2;

	cout << "z (x2 y2): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = (one-x2)*y1+(one-x3)*y2;

	cout << "z (x1 y2): " << endl << z << endl;

	cout << "---" << endl;
	cout << "y3" << endl;

	y1 = y[1];
	y2 = y[2];
	y3 = one - y[1] - y[2];
	y3.set_range(y[3].inf(), y[3].sup());

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = (one-x2)*y1+(one-x3)*y2;

	cout << "z (x3 y3): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = (one-x2)*y1+(one-x3)*y2;

	cout << "z (x2 y3): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = (one-x2)*y1+(one-x3)*y2;

	cout << "z (x1 y3): " << endl << z << endl;

}

void example_8() {

	affine::set_max_used_index(0);

	affine one(1.0);

	affine x[] = {
			affine(0.0),
			affine(0.6, 0.7),
			affine(0.2, 0.3),
			affine(0.1, 0.2)
	};

	affine y[] = {
			affine(0.0),
			affine(0.1, 0.3),
			affine(0.2, 0.4),
			affine(0.3, 0.5)
	};

	affine x1, x2, x3, y1, y2, y3;

	y1 = y[1];
	y2 = y[2];
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = x[3];

	affine z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "=========================================================" << endl;
	cout << "y1 = [0.1, 0.3]" << endl;
	cout << "=========================================================" << endl;

	cout << "Initial z" << endl << z << endl;

	cout << "---" << endl;
	cout << "y1" << endl;

	y1 = one - y[2] - y[3];
	y1.set_range(y[1].inf(), y[1].sup());
	y2 = y[2];
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x3 y1): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x2 y1): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x1 y1): " << endl << z << endl;

	cout << "---" << endl;
	cout << "y2" << endl;

	y1 = y[1];
	y2 = one - y[1] - y[3];
	y2.set_range(y[2].inf(), y[2].sup());
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x3 y2): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x2 y2): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x1 y2): " << endl << z << endl;

	cout << "---" << endl;
	cout << "y3" << endl;

	y1 = y[1];
	y2 = y[2];
	y3 = one - y[1] - y[2];
	y3.set_range(y[3].inf(), y[3].sup());

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x3 y3): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x2 y3): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x1 y3): " << endl << z << endl;

}

void example_9() {

	affine::set_max_used_index(0);

	affine one(1.0);

	affine x[] = {
			affine(0.0),
			affine(0.6, 0.7),
			affine(0.2, 0.3),
			affine(0.1, 0.2)
	};

	affine y[] = {
			affine(0.0),
			affine(0.3, 0.4),
			affine(0.2, 0.4),
			affine(0.3, 0.5)
	};

	affine x1, x2, x3, y1, y2, y3;

	y1 = y[1];
	y2 = y[2];
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = x[3];

	affine z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "=========================================================" << endl;
	cout << "y1 = [0.3, 0.4]" << endl;
	cout << "=========================================================" << endl;

	cout << "Initial z" << endl << z << endl;

	cout << "---" << endl;
	cout << "y1" << endl;

	y1 = one - y[2] - y[3];
	y1.set_range(y[1].inf(), y[1].sup());
	y2 = y[2];
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x3 y1): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x2 y1): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x1 y1): " << endl << z << endl;

	cout << "---" << endl;
	cout << "y2" << endl;

	y1 = y[1];
	y2 = one - y[1] - y[3];
	y2.set_range(y[2].inf(), y[2].sup());
	y3 = y[3];

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x3 y2): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x2 y2): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x1 y2): " << endl << z << endl;

	cout << "---" << endl;
	cout << "y3" << endl;

	y1 = y[1];
	y2 = y[2];
	y3 = one - y[1] - y[2];
	y3.set_range(y[3].inf(), y[3].sup());

	x1 = x[1];
	x2 = x[2];
	x3 = one - x[1] - x[2];
	x3.set_range(x[3].inf(), x[3].sup());

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x3 y3): " << endl << z << endl;

	x1 = x[1];
	x2 = one - x[1] - x[3];
	x2.set_range(x[2].inf(), x[2].sup());
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x2 y3): " << endl << z << endl;

	x1 = one - x[2] -x[3];
	x1.set_range(x[1].inf(),x[1].sup());
	x2 = x[2];
	x3 = x[3];

	z = (x1+x3)*y1+(x1+x2)*y2;

	cout << "z (x1 y3): " << endl << z << endl;

}

int main() {

	run();

	/*
	example_1();

	example_2();

	example_3();

	example_4();

	example_5();

	example_6();

	example_7();

	example_8();

	example_9();
	*/

	return 0;
}
