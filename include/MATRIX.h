#pragma once
#include<iostream>
#include<fstream>
#define NULL_MAT { 0,0,NULL }

typedef struct
{
	int m;
	int n;
	double** value;
} matrix;

class MATRIX
{
public:
	/// <summary>
	/// Initialize a matrix, set it's elements to 0.
	/// </summary>
	/// <param name="A">Matrix to be initialized.</param>
	/// <param name="m">Rows</param>
	/// <param name="n">Columns</param>
	/// <returns></returns>
	bool initialize(matrix& A, int m, int n);
	/// <summary>
	/// Print a matrix.
	/// </summary>
	/// <param name="A">Matrix to be print.</param>
	void print(matrix& A);
	/// <summary>
	/// Print a matrix into data.txt.
	/// </summary>
	/// <param name="A">Matrix to be print.</param>
	void fprint(matrix& A);
	/// <summary>
	/// Delete a matrix.
	/// </summary>
	/// <param name="A">Matrix to be disposed.</param>
	void dispose(matrix& A);
	/// <summary>
	/// Set a matrix to 0.
	/// </summary>
	/// <param name="A">Matrix to be set to 0.</param>
	void setZero(matrix& A);
	/// <summary>
	/// Copy a matrix.
	/// </summary>
	/// <param name="A">Source matrix.</param>
	/// <param name="B">Distination matrix.</param>
	void copy(matrix& A, matrix& B);
	/// <summary>
	/// Matrix plus A+B=C.
	/// </summary>
	/// <param name="A">Matrix A</param>
	/// <param name="B">Matrix B</param>
	/// <param name="C">Matrix C</param>
	void plus(matrix& A, matrix& B, matrix& C);
	/// <summary>
	/// Matrix scalar multiplication aA=B.
	/// </summary>
	/// <param name="a">Scalar a</param>
	/// <param name="A">Matrix A</param>
	/// <param name="B">Matrix B</param>
	void scalarProduct(double a, matrix& A, matrix& B);
	/// <summary>
	/// Matrix multiplication AB=C.
	/// </summary>
	/// <param name="A">Matrix A</param>
	/// <param name="B">Matrix B</param>
	/// <param name="C">Matrix C</param>
	void matrixProduct(matrix& A, matrix& B, matrix& C);
	double det(matrix& A);
};