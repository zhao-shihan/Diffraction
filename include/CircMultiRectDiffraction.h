#pragma once
#include "FKdiffraction.h"
using namespace std;

// ���η������е����䡣
class CircMultiRectDiffraction
{
public:
	FKdiffraction fk;
	//�׵�����
	int apCount = 1;
	//�������ߵľ���
	double d = 1;
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
	/// ���㻷�η������еĻ�������������֡�
	/// </summary>
	/// <param name="x2">������Ҫ����ĵ������������</param>
	/// <returns></returns>
	complex<double> integral(vector& x2);
	/// <summary>
	/// �ͷ�x1s��
	/// </summary>
	void dispose();
private:
	const double pi = 3.1415926535897932;
	LinearAlgebra la;
	//��������Ĳ��������Ĳ�ͬ�׵��ڿ��ϵ�ÿһ����������������
	vector*** x1s = NULL;
};