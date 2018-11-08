#include <iostream>
#include <stdlib.h>
#include <cmath>
#include<string>
#include "structs.h"
using namespace std;
//初始化点和概率密度
void initialPoints(double2Ds &points, int nv) {
	points.elem = new double2D[nv];
	points.n = nv;
};
void initialProbability(double2Ds points, double *value) {
	double u1, u2, s1, s2, r12;
	//u1 = -1.0; u2 = -1.0;
	//s1 = s2 = sqrt(0.1);
	//r12 = 0.01;
	u1 = -2; u2 = -1.8; s1 = s2 = 1;
	r12 = 0.0;

	for (int i = 0; i < points.n; i++)
	{
		double x = points.elem[i].x;
		double y = points.elem[i].y;
		double temp;
		temp = s2*s2*(x - u1)*(x - u1) + s1*s1*(y - u2)*(y - u2) - 2 * r12*s1*s2*(y - u2)*(x - u1);
		temp = exp(-0.5*temp / ((1.0 - r12*r12)*s1*s1*s2*s2));
		value[i] = temp / (2.0*PI*s1*s2*sqrt(1.0 - r12*r12));
	}
}
void initialMoment(double2Ds points, Moment *moments) {
	for (int i = 0; i < points.n; i++)
	{
		moments[i]._1 = points.elem[i].x;
		moments[i]._2 = points.elem[i].y;
		moments[i]._3 = points.elem[i].x*points.elem[i].y;
		moments[i]._4 = points.elem[i].x*points.elem[i].x;
		moments[i]._5 = points.elem[i].y*points.elem[i].y;
	}
}
//通过矩计算概率密度
double p(double2D point, Moment moment) {
	double x, y, u1, u2, s1, s2, r12;
	x = point.x; y = point.y;
	u1 = moment._1; u2 = moment._2;
	s1 = sqrt(moment._4 - moment._1*moment._1);
	s2 = sqrt(moment._5 - moment._2*moment._2);
	r12 = (moment._3 - moment._1*moment._2) / (s1*s2);
	double temp, ret;
	temp = s2*s2*(x - u1)*(x - u1) + s1*s1*(y - u2)*(y - u2) - 2 * r12*s1*s2*(y - u2)*(x - u1);
	temp = exp(-0.5*temp / ((1.0 - r12*r12)*s1*s1*s2*s2));
	ret = temp / (2.0*PI*s1*s2*sqrt(1.0 - r12*r12));
	if (isnan(ret)) ret = 0.0;
	return ret;
}
//路径积分
void pathIntegration(double2Ds points, Moment *moments, double *area, double *old_probability, double *new_probability) {

	for (int i = 0; i < points.n; i++)
	{
		for (int j = 0; j < points.n; j++)
		{
			double pp = p(points.elem[i], moments[j]);

			new_probability[i] += area[j] * pp*old_probability[j];
		}
	}
}
//点的排序与去重
int double3DCompare(const void *v1, const void *v2) {
	double2D *x, *y;
	x = (double2D*)v1;
	y = (double2D*)v2;
	if (x->x < y->x)
		return(-1);
	else if (x->x > y->x)
		return(1);
	else
		return(0);
}
void de_duplicatingPoints(double2Ds &points) {
	int n = points.n;
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			if (points.elem[i].x == points.elem[j].x  &&  points.elem[i].y == points.elem[j].y) {
				points.elem[i] = points.elem[n - 1];
				--n; --i;
				break;
			}
		}
	}
	points.n = n;
}
void sortAndDeduplicatePointsWithValue(double2Ds &points) {
	de_duplicatingPoints(points);
	qsort(points.elem, points.n, sizeof(double2D), double3DCompare);
}
float ideal3D(float x, float y);
double* bubblePoint(double l, float(*ideal)(float,float), int &n);
void idealConstruction(double2D *points, double *value, int n);
void BubblePoint(double2Ds b_old, double2Ds &b_new, double *value) {
	// 构建一个分片线性插值函数
	idealConstruction(b_old.elem, value, b_old.n);
	int n = 0;
	// 进行泡泡布点
	double *data = bubblePoint(5.0, ideal3D,n);
	b_new.n = n;
	b_new.elem = (double2D*)malloc(n * sizeof(double2D));
	// 返回数据
	for (size_t i = 0; i < n; i++)
	{
		b_new.elem[i].x = data[i];
		b_new.elem[i].y = data[i + n];
	}
}
//主函数
int main() {
	double dt = 2.0*PI / 1.6;
	double t = dt;
	double2Ds points1;
	double2Ds points2;
	int3Ds triangles;
	double *area = NULL;
	double *value1 = NULL;
	double *value2 = NULL;
	Moment *moments1 = NULL;
	Moment *moments2 = NULL;
	initialPoints(points1, 4000);
	value1 = new double[points1.n]();
	value2 = new double[points1.n]();
	moments1 = new Moment[points1.n]();
	area = new double[points1.n]();
	manyRandPositions(5.0, points1.elem, points1.n);
	sortAndDeduplicatePointsWithValue(points1);
	initialMoment(points1, moments1);
	oderk4(moments1, points1.n, dt, dt);
	initialProbability(points1, value1);
	Triangulate(points1, triangles);
	calculateAreaOfEveryPoint(area, points1, triangles);
	pathIntegration(points1, moments1, area, value1, value2);
	delete[] value1;
	value1 = value2;
	value2 = NULL;
	char title[40];
	sprintf(title, "%.0f.txt", t);
	outputData3D(title, points1, value1);
	sprintf(title, "%.0f.jpg", t);
	paintTriangles(title, points1, triangles);
	t = t + dt;
	while (t<400) {
		BubblePoint(points1, points2, value1);
		char title[40];
		sprintf(title, "points%.0f.txt", t);
		outputData(title, points2);
		sortAndDeduplicatePointsWithValue(points2);
		if (value2 != NULL) delete[] value2;
		value2 = new double[points2.n];
		for (int i = 0; i < points2.n; i++) {
			//value2[i] = value1[i];
		}
		interpolation(points2, value2, points1, value1, triangles);
		delete[] value1; value1 = value2; value2 = NULL;
		if (moments2 != NULL) delete[] moments2;
		moments2 = new Moment[points2.n]();
		//for (int i = 0; i < points2.n; i++) { moments2[i] = moments1[i]; }
		interpolationMoment(points2, moments2, points1, moments1, triangles);
		delete[] moments1; moments1 = moments2; moments2 = NULL;
		delete[] points1.elem;
		points1 = points2;
		delete[] triangles.elem;
		Triangulate(points2, triangles);
		delete[] area;
		area = new double[points1.n]();
		calculateAreaOfEveryPoint(area, points1, triangles);

		double sum_area = 0.0;
		for (int i = 0; i < points1.n; i++) { sum_area += area[i]; }
		printf("面积：%lf\n", sum_area);


		//oderk4(moments1, points1.n, dt, dt);
		value2 = new double[points1.n]();
		pathIntegration(points1, moments1, area, value1, value2);
		delete[] value1; value1 = value2; value2 = NULL;
		sprintf(title, "%.0f.txt", t);
		outputData3D(title, points1, value1);
		sprintf(title, "%.0f.jpg", t);
		paintTriangles(title, points1, triangles);
		t = t + dt;
	}
	return 0;
}