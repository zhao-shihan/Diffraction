#include "RectDiffraction.h"

void RectDiffraction::initialize()
{
	x1s = new vector * [dimb];
	for (int i = 0; i < dimb; i++)
	{
		x1s[i] = new vector[dima];
		for (int j = 0; j < dima; j++)
		{
			x1s[i][j] = NULL_VEC;
			la.vec.initialize(x1s[i][j], 3);
		}
	}
	const double dx = a / (double)dima;
	const double dy = b / (double)dimb;
	double x = -a / 2.0 + dx / 2.0;
	double y = -b / 2.0 + dy / 2.0;
	for (int i = 0; i < dimb; i++)
	{
		for (int j = 0; j < dima; j++)
		{
			x1s[i][j].value[0] = x;
			x1s[i][j].value[1] = y;
			x += dx;
		}
		x = -a / 2.0 + dx / 2.0;
		y += dy;
	}
}

complex<double> RectDiffraction::integral(vector& x2)
{
	const double dx = a / (double)dima;
	const double dy = b / (double)dimb;
	const double dS = dx * dy;
	complex<double> I{ 0,0 };
	for (int i = 0; i < dimb; i++)
	{
		for (int j = 0; j < dima; j++)
		{
			I += fk.f(x1s[i][j], x2) * dS;
		}
	}
	return I;
}

void RectDiffraction::dispose()
{
	for (int i = 0; i < dimb; i++)
	{
		for (int j = 0; j < dima; j++)
		{
			la.vec.dispose(x1s[i][j]);
		}
		delete[]x1s[i];
	}
	delete[]x1s;
}