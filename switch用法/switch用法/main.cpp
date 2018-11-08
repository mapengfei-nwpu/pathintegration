#include <iostream>
#include <stdlib.h>
#include <cmath>
#include<string>
#include "structs.h"
using namespace std;
//初始化点和概率密度
void initialPoints(double2Ds &points,int nv) {
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
void pathIntegration(double2Ds points,Moment *moments,double *area, double *old_probability, double *new_probability) {
	//               double* points_x, double* points_y, double *area, double *old
	for (int i = 0; i < points.n; i++)
	{
		for (int j = 0; j < points.n; j++)
		{
			double pp = p(points.elem[i], moments[j]);
			
			new_probability[i] += area[j]*pp*old_probability[j];
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
//主函数
int main() {
	double dt = 2.0*PI/1.6;
	double t = dt;
	double2Ds points1;
	double2Ds points2;
	double2Ds points3;
	double2Ds points4;
	int3Ds triangles;
	double *area = NULL;
	double *value1 = NULL; 
	double *value2 = NULL;
	Moment *moments1 = NULL;
	Moment *moments2 = NULL;
	initialPoints(points1, 51*51);
	value1		= new double[points1.n]();
	value2		= new double[points1.n]();
	moments1	= new Moment[points1.n]();
	area		= new double[points1.n]();
	//manyRandPositions(5.0, points1.elem, points1.n);
	//sortAndDeduplicatePointsWithValue(points1);
	for (size_t i = 0; i < 51; i++)
	{
		for (size_t j = 0; j < 51; j++)
		{
			points1.elem[i *51 + j].x = -5.0 + j*0.2;
			points1.elem[i * 51 + j].y = -5.0 + i*0.2;
		}
	}
	initialMoment(points1, moments1);
	
	

	oderk4(moments1, points1.n, dt, dt);
	
	initialProbability(points1, value1);
	for (size_t i = 0; i < 441; i++)
	{
		//printf("%lf  %lf  %lf  %lf  %lf  \n", moments1[i]._1, moments1[i]._2, moments1[i]._3, moments1[i]._4, moments1[i]._5);
		//printf("%lf\n", value1[i]);
	}
	Triangulate(points1, triangles);
	calculateAreaOfEveryPoint(area, points1, triangles);
	pathIntegration(points1, moments1, area, value1, value2);

	delete[] value1;
	value1 = value2;
	value2 = NULL;
	char title[40];
	sprintf(title, "%.0f.txt", t);
	outputData3D(title, points1, value2);
	sprintf(title, "%.0f.jpg", t);
	paintTriangles(title, points1, triangles);
	initialPoints(points4, 1500);
	manyRandPositions(5.0, points4.elem, points4.n);
	t = t + dt;
	while(t<400){
		initialPoints(points2,3000);
		initialPoints(points3, 1500);
		
		acceptanceRejectionSamplingMethod(points3, points1, value1, triangles);

		for (int i = 0; i < 1500; i++) { 
			points2.elem[i+1500] = points4.elem[i];
		}
		for (int i = 0; i < 1500; i++) {
			points2.elem[i] = points3.elem[i];
		}

		delete[]points3.elem;
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


double RHS(double M10, double M01, double M11, double M20, double M02, double TM, int MARK)
{
	//double ALPHA = 1.0, BETA = 0.3, YETA = 0.1, FF = 0.2, OMIGA = 1.2, SEGMA2 = 0.01;
	double GAMA = 0.2, BETA = 0.0, EPSON = 0.3, OMIGA1 = 1.0, OMIGA = 1.6, SIGMA1 = 1.0, SIGMA2 = 0.4;
	double ret;
	if (MARK == 1)
		ret = M01;
	if (MARK == 2)
		ret = -(GAMA*BETA / OMIGA1 / OMIGA1)*(3.0*M01*M02 - 2.0*M01*M01*M01) - OMIGA1*OMIGA1*EPSON*(3.0*M10*M20 - 2.0*M10*M10*M10) - GAMA*M01 - OMIGA1*OMIGA1*M10 + OMIGA1*SIGMA1*sin(OMIGA*TM);
	//ret = -YETA*M01 - ALPHA*M10 - 3.0*BETA*M10*M20 + 2.0*BETA*M10*M10*M10 + FF*cos(OMIGA*TM);

	if (MARK == 3)
		ret = M02 + (GAMA*BETA / OMIGA1 / OMIGA1)*(2.0*M01*M01*M01*M10 - 3.0*M02*M11) + OMIGA1*OMIGA1*EPSON*(2.0*M10*M10*M10*M10 - 3.0*M20*M20) - GAMA*M11 - OMIGA1*OMIGA1*M20 + OMIGA1*SIGMA1*M10*sin(OMIGA*TM);
	//ret = M02 - YETA*M11 - ALPHA*M20 - 3.0*BETA*M20*M20 + 2.0*BETA*M10*M10*M10*M10 + M10*FF*cos(OMIGA*TM);
	if (MARK == 4)
		ret = 2.0*M11;
	if (MARK == 5)
		ret = (GAMA*BETA / OMIGA1 / OMIGA1)*(4.0*M01*M01*M01*M01 - 6.0*M02*M02) + OMIGA1*OMIGA1*EPSON*(4.0*M10*M10*M10*M01 - 6.0*M20*M11) - 2.0*GAMA*M02 - 2.0*OMIGA1*OMIGA1*M11 + 2.0*OMIGA1*SIGMA1*M01*sin(OMIGA*TM) + OMIGA1*OMIGA1*SIGMA2*SIGMA2;
	//ret = -2.0*YETA*M02 - 2.0*ALPHA*M11 - 6.0*BETA*M20*M11 + 4.0*BETA*M10*M10*M10*M01 + SEGMA2 + 2.0*M01*FF*cos(OMIGA*TM);
	return ret;
}

void oderk4(Moment *moments, int n, double t, double dt) {

	double M10_1, M01_1, M11_1, M20_1, M02_1;
	double M10_2, M01_2, M11_2, M20_2, M02_2;
	double M10_3, M01_3, M11_3, M20_3, M02_3;
	for (int i = 0; i < n; i++) {
		double M10 = moments[i]._1;
		double M01 = moments[i]._2;
		double M11 = moments[i]._3;
		double M20 = moments[i]._4;
		double M02 = moments[i]._5;
		double DT0 = dt / 100;
		double TM0 = t - dt;
		while (abs(TM0 - t) > 1.0e-3) {
			//STEP 1
			M10_1 = M10 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 1);
			M01_1 = M01 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 2);
			M11_1 = M11 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 3);
			M20_1 = M20 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 4);
			M02_1 = M02 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 5);
			//STEP 2
			M10_2 = 0.75*M10 + 0.25*M10_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 1);
			M01_2 = 0.75*M01 + 0.25*M01_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 2);
			M11_2 = 0.75*M11 + 0.25*M11_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 3);
			M20_2 = 0.75*M20 + 0.25*M20_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 4);
			M02_2 = 0.75*M02 + 0.25*M02_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 5);
			//STEP 3
			M10_3 = M10 / 3.0 + 2.0*M10_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 1) / 3.0;
			M01_3 = M01 / 3.0 + 2.0*M01_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 2) / 3.0;
			M11_3 = M11 / 3.0 + 2.0*M11_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 3) / 3.0;
			M20_3 = M20 / 3.0 + 2.0*M20_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 4) / 3.0;
			M02_3 = M02 / 3.0 + 2.0*M02_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 5) / 3.0;
			//UPDATE
			M10 = M10_3;
			M01 = M01_3;
			M11 = M11_3;
			M20 = M20_3;
			M02 = M02_3;
			//LOOP CONDITIONS
			TM0 = TM0 + DT0;
			//cout << i << " ";
		}
		moments[i]._1 = M10;
		moments[i]._2 = M01;
		moments[i]._3 = M11;
		moments[i]._4 = M20;
		moments[i]._5 = M02;
	}
}