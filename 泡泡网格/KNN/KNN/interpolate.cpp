#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include"_head_.h"


#define WIDTH 50		//X轴方向区域数量
#define HEIGHT 50		//Y轴方向区域数量
#define SIDE  0.2		//区域边长
#define STARTX -5		//区域起始点
#define STARTY -5		//区域起始点

double weight[HEIGHT][WIDTH];

void pointInWhichArea(double2D point,int &i,int &j) {
	i = ((point.y - STARTY) / SIDE);
	j = ((point.x - STARTX) / SIDE);
	if (i < 0) i = 0;
	if (i >(WIDTH - 1)) i = WIDTH - 1;
	if (j < 0) j = 0;
	if (j > HEIGHT - 1) j = HEIGHT - 1;
}
void pointInWhichArea_2(double2D point, int &i, int &j) {
	i = ((point.y - STARTY + SIDE / 2.0) / SIDE);
	j = ((point.x - STARTX + SIDE / 2.0) / SIDE);
	if (i < 0) i = 0;
	if (i >WIDTH) i = WIDTH;
	if (j < 0) j = 0;
	if (j > HEIGHT) j = HEIGHT;
}
void func(double2D *points, double *value, int n) {
	int num[HEIGHT][WIDTH] = { 0 };
	int p, q;
	for (int i = 0; i < n; i++) {
		pointInWhichArea(points[i], p, q);
		++num[p][q];
		weight[p + 1][q + 1] += value[i];
	}
	for (int i = 0; i < HEIGHT; i++)
		for (int j = 0; j < WIDTH; j++)
			if (num[i][j] != 0)	weight[i + 1][j + 1] /= num[i][j];
}
double bilinearInterpolationInRectangle(double2D pp, double2D p1, double2D p2, double *f) {
	return ((p2.x - pp.x)*(p2.y - pp.y)*f[0] + (pp.x - p1.x)*(p2.y - pp.y)*f[1]
		+ (p2.x - pp.x)*(pp.y - p1.y)*f[2] + (pp.x - p1.x)*(pp.y - p1.y)*f[3]) / (p2.x - p1.x) / (p2.y - p1.y);
}

float ideal2D(double2D point) {
	int p, q;
	pointInWhichArea_2(point, p, q);
	point.x += SIDE / 2.0;
	point.y += SIDE / 2.0;

	double2D dot1, dot2;
	dot1.x = STARTX + q*SIDE;
	dot1.y = STARTY + p*SIDE;
	dot2.x = STARTX + q*SIDE + SIDE;
	dot2.y = STARTY + p*SIDE + SIDE;


	double f[4];
	f[0] = weight[p][q];
	f[1] = weight[p + 1][q];
	f[2] = weight[p][q + 1];
	f[3] = weight[p + 1][q + 1];

	double a;
	a = bilinearInterpolationInRectangle(point, dot1, dot2, f);
	
	if (a < 0) a = 0;
	if (a > 0.1) {
		a = 0.05;
	}
	else a = 0.1-a*0.5;
	return a;
}


float ideal3D(float x,float y) {
	double2D temp;
	temp.x = x;
	temp.y = y;
	return ideal2D(temp);
}
void idealConstruction(double2D *points, double *value, int n){
	func(points, value, n);
}
