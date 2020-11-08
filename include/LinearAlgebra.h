#pragma once
#include"VECTOR.h"
#include"MATRIX.h"

class LinearAlgebra
{
public:
	VECTOR vec;
	MATRIX mat;
	void product(matrix& A, vector& v, vector& v1);
};

