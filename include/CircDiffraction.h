#pragma once
#include "FKdiffraction.h"
using namespace std;

//����Բ�����䡣
class CircDiffraction
{
public:
	FKdiffraction fk;
	//�׵ĳ��ᣨ��x����
	double a = 1;
	//�׵Ķ��ᣨ��y����
	double b = 1;
	//����ά��
	int dimr = 1;
	//����ά��
	int dimt = 1;
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
	const double pi = 3.1415926535897932;
	//��������Ĳ����������ڿ��ϵ�ÿһ����������������
	vector** x1 = NULL;
	//��������Ĳ����������ڿ��ϵ�ÿһ���������Ԫ��
	double* dS = NULL;
};

