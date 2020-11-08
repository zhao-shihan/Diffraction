#pragma once
#include<iostream>
#include<cmath>
#define NULL_VEC { 0,NULL }

typedef struct
{
	int n;
	double* value;
} vector;

class VECTOR
{
public:
	bool initialize(vector& v, int n);
	void print(vector& v);
	void dispose(vector& v);
	void setZero(vector& v);
	void copy(vector& v1, vector& v2);
	void plus(vector& v1, vector& v2, vector& v3);
	void minus(vector& v1, vector& v2, vector& v3);
	void scalarProduct(double a, vector& v1, vector& v2);
	double dotProduct(vector& v1, vector& v2);
	double norm(vector& v);
	//Under construction...
	void crossProduct(vector& v1, vector& v2, vector& v3);
	double cosAngle(vector& v1, vector& v2);
	double angle(vector& v1, vector& v2);
};

