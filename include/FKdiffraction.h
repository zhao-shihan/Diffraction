#pragma once
#include <complex>
#include "LinearAlgebra.h"
using namespace std;

class FKdiffraction
{
public:
	double A0 = 1;
	double lambda = 632.8e-9;
	double s0 = 1;
	double r0 = 1;
	matrix intensity = NULL_MAT;
	complex<double> f(vector& x1, vector& x2);
private:
	LinearAlgebra la;
	const double pi = 3.1415926535897932;
};