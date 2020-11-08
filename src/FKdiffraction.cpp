#include "FKdiffraction.h"

complex<double> FKdiffraction::f(vector& x1, vector& x2)
{
	vector s = NULL_VEC;
	la.vec.initialize(s, 3);
	s.value[2] = -s0;
	vector r = NULL_VEC;
	la.vec.initialize(r, 3);
	r.value[2] = r0;
	vector n = NULL_VEC;
	la.vec.initialize(n, 3);
	n.value[2] = 1.0;
	vector temp = NULL_VEC;
	la.vec.minus(x1, s, temp);
	double S = la.vec.norm(temp);
	la.vec.minus(r, x1, temp);
	la.vec.plus(temp, x2, temp);
	double R = la.vec.norm(temp);
	la.vec.dispose(temp);
	const double k = 2.0 * pi / lambda;
	complex<double> f0val{ -sin(k * (R + S)),cos(k * (R + S)) };
	f0val /= R * S;
	f0val *= la.vec.cosAngle(n, r) - la.vec.cosAngle(n, s);
	la.vec.dispose(s);
	la.vec.dispose(r);
	la.vec.dispose(n);
	double Ac = -A0 / (2.0 * lambda);
	return Ac * f0val;
}