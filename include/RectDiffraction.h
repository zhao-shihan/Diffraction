#pragma once
#include "FKdiffraction.h"
using namespace std;

// ���������䡣
class RectDiffraction
{
public:
	FKdiffraction fk;
	//�׵Ŀ�ȣ���x����
	double a = 1;
	//�׵ĸ߶ȣ���y����
	double b = 1;
	//�׿��ά��
	int dima = 1;
	//�׸ߵ�ά��
	int dimb = 1;
	//���Ŀ�ȣ���x����
	double a0 = 1;
	//���ĸ߶ȣ���y����
	double b0 = 1;
	//�����ά��
	int dima0 = 1;
	//���ߵ�ά��
	int dimb0 = 1;
	/// <summary>
	/// ���ݲ�����ʼ��x1s��
	/// </summary>
	void initialize();
	/// <summary>
	/// ���㷽�׵Ļ�������������֡�
	/// </summary>
	/// <param name="x2">������Ҫ����ĵ������������</param>
	/// <returns></returns>
	complex<double> integral(vector& x2);
	/// <summary>
	/// �ͷ�x1s��
	/// </summary>
	void dispose();
private:
	LinearAlgebra la;
	//��������Ĳ����������ڿ��ϵ�ÿһ����������������
	vector** x1s = NULL;
};