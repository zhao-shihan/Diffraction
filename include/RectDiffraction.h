#pragma once
#include "FKdiffraction.h"
using namespace std;

// 单方孔衍射。
class RectDiffraction
{
public:
	FKdiffraction fk;
	//孔的宽度（沿x方向）
	double a = 1;
	//孔的高度（沿y方向）
	double b = 1;
	//孔宽的维数
	int dima = 1;
	//孔高的维数
	int dimb = 1;
	//屏的宽度（沿x方向）
	double a0 = 1;
	//屏的高度（沿y方向）
	double b0 = 1;
	//屏宽的维数
	int dima0 = 1;
	//屏高的维数
	int dimb0 = 1;
	/// <summary>
	/// 根据参数初始化x1s。
	/// </summary>
	void initialize();
	/// <summary>
	/// 计算方孔的基尔霍夫衍射积分。
	/// </summary>
	/// <param name="x2">屏上需要计算的点的坐标向量。</param>
	/// <returns></returns>
	complex<double> integral(vector& x2);
	/// <summary>
	/// 释放x1s。
	/// </summary>
	void dispose();
private:
	LinearAlgebra la;
	//根据输入的参数创建的在孔上的每一个格点的坐标向量。
	vector** x1s = NULL;
};