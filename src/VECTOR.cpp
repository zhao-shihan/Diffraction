#include "VECTOR.h"

bool VECTOR::initialize(vector& v, int n)
{
	dispose(v);
	v = NULL_VEC;
	if (n <= 0)
	{
		return false;
	}
	try
	{
		v.value = new double[n];
	}
	catch (std::bad_alloc)
	{
		v = NULL_VEC;
		return false;
	}
	v.n = n;
	setZero(v);
	return true;
}

void VECTOR::dispose(vector& v)
{
	if (v.value != NULL)
	{
		delete[] v.value;
		v.value = NULL;
	}
}

void VECTOR::print(vector& v)
{
	for (int i = 0; i < v.n; i++)
	{
		std::cout << v.value[i] << std::endl;
	}
}

void VECTOR::setZero(vector& v)
{
	if (v.value != NULL)
	{
		for (int i = 0; i < v.n; i++)
		{
			v.value[i] = 0;
		}
	}
}

void VECTOR::copy(vector& v1, vector& v2)
{
	if (v1.value != NULL)
	{
		bool flag = true;
		if (v2.value == NULL || v1.n != v2.n)
		{
			flag = initialize(v2, v1.n);
		}
		if (flag == true)
		{
			for (int i = 0; i < v1.n; i++)
			{
				v2.value[i] = v1.value[i];
			}
		}
	}
}

void VECTOR::plus(vector& v1, vector& v2, vector& v3)
{
	if (v1.value != NULL && v2.value != NULL && v1.n == v2.n)
	{
		bool flag = true;
		if (v3.value == NULL || v1.n != v3.n)
		{
			flag = initialize(v3, v1.n);
		}
		if (flag == true)
		{
			for (int i = 0; i < v1.n; i++)
			{
				v3.value[i] = v1.value[i] + v2.value[i];
			}
		}
	}
}

void VECTOR::minus(vector& v1, vector& v2, vector& v3)
{
	if (v1.value != NULL && v2.value != NULL && v1.n == v2.n)
	{
		bool flag = true;
		if (v3.value == NULL || v1.n != v3.n)
		{
			flag = initialize(v3, v1.n);
		}
		if (flag == true)
		{
			for (int i = 0; i < v1.n; i++)
			{
				v3.value[i] = v1.value[i] - v2.value[i];
			}
		}
	}
}

void VECTOR::scalarProduct(double a, vector& v1, vector& v2)
{
	if (v1.value != NULL)
	{
		bool flag = true;
		if (v2.value == NULL || v1.n != v2.n)
		{
			flag = initialize(v2, v1.n);
		}
		if (flag == true)
		{
			for (int i = 0; i < v1.n; i++)
			{
				v2.value[i] = a * v1.value[i];
			}
		}
	}
}

double VECTOR::dotProduct(vector& v1, vector& v2)
{
	if (v1.value != NULL && v2.value != NULL && v1.n == v2.n)
	{
		double I = 0;
		for (int i = 0; i < v1.n; i++)
		{
			I += v1.value[i] * v2.value[i];
		}
		return I;
	}
	else
	{
		return 0;
	}
}

double VECTOR::norm(vector& v)
{
	double norm2 = dotProduct(v, v);
	return sqrt(norm2);
}

void VECTOR::crossProduct(vector& v1, vector& v2, vector& v3)
{
	if (v1.value != NULL && v2.value != NULL && v1.n != v2.n)
	{
		bool flag = true;
		if (v3.value == NULL || v3.n != v1.n)
		{
			flag = initialize(v3, v1.n);
		}
		if (flag == true)
		{
			//!!!
		}
	}
}

double VECTOR::cosAngle(vector& v1, vector& v2)
{
	double dot = dotProduct(v1, v2);
	double norm1 = norm(v1);
	double norm2 = norm(v2);
	return dot / (norm1 * norm2);
}

double VECTOR::angle(vector& v1, vector& v2)
{
	double cosAngl = cosAngle(v1, v2);
	return acos(cosAngl);
}