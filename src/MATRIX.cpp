#include "MATRIX.h"

bool MATRIX::initialize(matrix& A, int m, int n)
{
	dispose(A);
	A = NULL_MAT;
	if (m <= 0 || n <= 0)
	{
		return false;
	}
	try
	{
		A.value = new double* [m];
	}
	catch (std::bad_alloc)
	{
		A = NULL_MAT;
		return false;
	}
	for (int i = 0; i < m; i++)
	{
		try
		{
			A.value[i] = new double[n];
		}
		catch (std::bad_alloc)
		{
			for (int j = 0; j < i - 1; j++)
			{
				delete[] A.value[j];
			}
			delete[] A.value;
			A = NULL_MAT;
			return false;
		}
	}
	A.m = m;
	A.n = n;
	setZero(A);
	return true;
}

void MATRIX::dispose(matrix& A)
{
	if (A.value != NULL)
	{
		for (int i = 0; i < A.m; i++)
		{
			delete[] A.value[i];
			A.value[i] = NULL;
		}
		delete[] A.value;
		A.value = NULL;
	}
}

void MATRIX::print(matrix& A)
{
	for (int i = 0; i < A.m; i++)
	{
		for (int j = 0; j < A.n; j++)
		{
			std::cout << A.value[i][j] << ' ';
		}
		std::cout << std::endl;
	}
}

void MATRIX::fprint(matrix& A)
{
	std::ofstream fs;
	fs.open("data.csv");
	if (fs)
	{
		for (int i = 0; i < A.m; i++)
		{
			for (int j = 0; j < A.n; j++)
			{
				fs << A.value[i][j] << ',';
			}
			fs << std::endl;
		}
		fs.close();
	}
}

void MATRIX::setZero(matrix& A)
{
	if (A.value != NULL)
	{
		for (int i = 0; i < A.m; i++)
		{
			for (int j = 0; j < A.n; j++)
			{
				A.value[i][j] = 0;
			}
		}
	}
}

void MATRIX::copy(matrix& A, matrix& B)
{
	if (A.value != NULL)
	{
		bool flag = true;
		if (B.value == NULL || A.m != B.m || A.n != B.n)
		{
			flag = initialize(B, A.m, A.n);
		}
		if (flag)
		{
			for (int i = 0; i < A.m; i++)
			{
				for (int j = 0; j < A.n; j++)
				{
					B.value[i][j] = A.value[i][j];
				}
			}
		}
	}
}

void MATRIX::plus(matrix& A, matrix& B, matrix& C)
{
	if (A.value != NULL && B.value != NULL && A.m == B.m && A.n == B.n)
	{
		bool flag = true;
		if (C.value == NULL || A.m != C.m || A.n != C.n)
		{
			flag = initialize(C, A.m, A.n);
		}
		if (flag)
		{
			for (int i = 0; i < A.m; i++)
			{
				for (int j = 0; j < A.n; j++)
				{
					C.value[i][j] = A.value[i][j] + B.value[i][j];
				}
			}
		}
	}
}

void MATRIX::scalarProduct(double a, matrix& A, matrix& B)
{
	if (A.value != NULL)
	{
		bool flag = true;
		if (B.value == NULL || A.n != B.n || A.m != B.m)
		{
			flag = initialize(B, A.m, A.n);
		}
		if (flag)
		{
			for (int i = 0; i < A.m; i++)
			{
				for (int j = 0; j < A.n; j++)
				{
					B.value[i][j] = a * A.value[i][j];
				}
			}
		}
	}
}

void MATRIX::matrixProduct(matrix& A, matrix& B, matrix& C)
{
	if (A.value != NULL && B.value != NULL && A.n == B.m)
	{
		bool flag = true;
		if (C.value == NULL || A.m != C.m || B.n != C.n)
		{
			flag = initialize(C, A.m, B.n);
		}
		if (flag)
		{
			for (int i = 0; i < A.m; i++)
			{
				for (int j = 0; j < B.n; j++)
				{
					for (int k = 0; k < A.n; k++)
					{
						C.value[i][j] = A.value[i][k] * B.value[k][j];
					}
				}
			}
		}
	}
}

double MATRIX::det(matrix& A)
{
	return 0;
}