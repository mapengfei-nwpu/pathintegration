#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"structs.h"


#define PI 3.1415926

/*int main() {


	//���ȳ�ʼ���õ�����ĸ����ܶȺ;�
	//���ݸ����ܶȽ������ݲ���
	//����approximation���������µĵ�ľغ͸����ܶ�
	//���þغ͸����ܶȼ����µĸ����ܶ�
	//�����µĸ����ܶ��ٴ����ݲ���
	system("pause");
}*/

void fun() {
	float MEAN1 = -1.0;
	float MEAN2 = -1.0;
	float VAR1 = sqrtf(0.1);
	float VAR2 = sqrtf(0.1);
	float ROU12 = 0.01;
	float2D points[100];
	for (size_t i = 0; i < 100; i++)
	{
		/*PM(1, i) = points[i].x;
		PM(2, i) = points[i].y;
		PM(3, i) = VAR2**2.0*(PM(1, I) - MEAN1)**2.0 - 2.0*VAR1*VAR2*ROU12*(PM(1, I) - MEAN1)*(PM(2, I) - MEAN2) + VAR1**2.0*(PM(2, I) - MEAN2)**2.0;
		PM(3, i) = DEXP(-PM(3, I) / (2.0*VAR1**2.0*VAR2**2.0*(1.0 - ROU12**2.0)));
		PM(3, i) = PM(3, I) / (2.0*PI*VAR1*VAR2*DSQRT(1.0 - ROU12**2.0));*/
	}
}

//�ط��̵��Ҷ��
void RHS() {
	int a = 1;
	if (a == 1)printf("1");
	if (a == 2)printf("2");
	if (a == 3)printf("3");
	if (a == 4)printf("4");
	if (a == 5)printf("5");
}

//moment�Ǿء�ÿ��ֻҪ���������С�ÿ������5*n��Ԫ�ء�
//��ǰ��ָ�ʽ����ء�
void forwardDifference(float moment) {

}
//��䶯���Ժ󣬹����µĵ��ϵ�ֵ��
//�������֣��أ������ܶȺ���������Ҫ�õ����������
void approximation(float3Ds origin, float3Ds result) {

}
//����֡��þغ;ɵĸ����ܶȺ�������µĸ����ܶȺ�����
void integral() {
	//����㼯�����������ʷ֡�
	//Ȼ������֣��õ��µĸ����ܶȺ���
}
