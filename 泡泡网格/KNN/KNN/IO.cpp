#include<stdio.h>
#include"structs.h"
void inputData(char *title, double2Ds points) {
	FILE *fp;
	fp = fopen(title, "r");
	for (int i = 0; i < points.n; ++i)
	{
		fscanf(fp, "%lf %lf", &(points.elem[i].x), &(points.elem[i].y));
	}
	fclose(fp);
}
void inputData3D(char *title, double2Ds points, double *value) {
	FILE *fp;
	fp = fopen(title, "r");
	for (int i = 0; i < points.n; ++i)
	{
		fscanf(fp, "%lf %lf %lf", &(points.elem[i].x), &(points.elem[i].y), &(value[i]));
	}
	fclose(fp);
}
void outputData(char *title, double2Ds points) {
	FILE *fp;
	fp = fopen(title, "w");
	for (int i = 0; i < points.n; ++i)
	{
		fprintf(fp, "%lf %lf\n", points.elem[i].x, points.elem[i].y);
	}
	fclose(fp);
}
void outputData3D(char *title, double2Ds points, double *value) {
	FILE *fp;
	fp = fopen(title, "w");
	for (int i = 0; i < points.n; ++i)
	{
		fprintf(fp, "%lf %lf %lf\n", points.elem[i].x, points.elem[i].y, value[i]);
	}
	fclose(fp);
}
