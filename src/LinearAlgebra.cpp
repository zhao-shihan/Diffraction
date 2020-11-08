#include "LinearAlgebra.h"

void LinearAlgebra::product(matrix& A, vector& v, vector& v1)
{
	if (A.value != NULL && v.value != NULL && v.n == A.n)
	{
		bool flag = true;
		if (v1.value == NULL && v1.n != A.m)
		{
			flag = vec.initialize(v1, A.m);
		}
		if (flag == true)
		{
			for (int i = 0; i < A.m; i++)
			{
				for (int j = 0; j < A.n; j++)
				{
					v1.value[i] += A.value[i][j] * v.value[j];
				}
			}
		}
	}
}