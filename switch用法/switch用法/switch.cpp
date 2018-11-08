#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"structs.h"


#define PI 3.1415926

/*int main() {


	//首先初始化得到最初的概率密度和矩
	//根据概率密度进行泡泡布点
	//利用approximation函数计算新的点的矩和概率密度
	//利用矩和概率密度计算新的概率密度
	//根据新的概率密度再次泡泡布点
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

//矩方程的右端项。
void RHS() {
	int a = 1;
	if (a == 1)printf("1");
	if (a == 2)printf("2");
	if (a == 3)printf("3");
	if (a == 4)printf("4");
	if (a == 5)printf("5");
}

//moment是矩。每次只要更新它就行。每个矩有5*n个元素。
//向前差分格式计算矩。
void forwardDifference(float moment) {

}
//点变动了以后，估计新的点上的值。
//两个部分，矩，概率密度函数，都需要用到这个函数。
void approximation(float3Ds origin, float3Ds result) {

}
//求积分。用矩和旧的概率密度函数求出新的概率密度函数。
void integral() {
	//输入点集，进行三角剖分。
	//然后求积分，得到新的概率密度函数
}
