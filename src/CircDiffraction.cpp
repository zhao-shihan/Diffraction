#include "CircDiffraction.h"

void CircDiffraction::initialize()
{
	x1 = new vector * [dimr];
	dS = new double[dimr];
	for (int i = 0; i < dimr; i++)
	{
		x1[i] = new vector[dimt];
		dS[i] = 0;
		for (int j = 0; j < dimt; j++)
		{
			x1[i][j] = NULL_VEC;
			la.vec.initialize(x1[i][j], 3);
		}
	}
	const double A = a / 2.0;
	const double B = b / 2.0;
	const double dr = 1.0 / (double)dimr;
	const double dt = 2.0 * pi / (double)dimt;
	const double drdt = dr * dt;
	double r = dr / 2.0;
	double t = dt / 2.0;
	for (int i = 0; i < dimr; i++)
	{
		dS[i] = A * B * r * drdt;
		for (int j = 0; j < dimt; j++)
		{
			x1[i][j].value[0] = A * r * cos(t);
			x1[i][j].value[1] = B * r * sin(t);
			t += dt;
		}
		t = dt / 2.0;
		r += dr;
	}
}

complex<double> CircDiffraction::integral(vector& x2)
{
	complex<double> I{ 0,0 };
	for (int i = 0; i < dimr; i++)
	{
		for (int j = 0; j < dimt; j++)
		{
			I += fk.f(x1[i][j], x2) * dS[i];
		}
	}
	return I;
}

void CircDiffraction::dispose()
{
	for (int i = 0; i < dimr; i++)
	{
		for (int j = 0; j < dimt; j++)
		{
			la.vec.dispose(x1[i][j]);
		}
		delete[]x1[i];
	}
	delete[]dS;
	delete[]x1;
}