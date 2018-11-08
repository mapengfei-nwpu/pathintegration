#pragma once
#define PI 3.1415926535898
#define EPSILON  0.0001
#define NUM 2000
typedef struct tag_double2D {
	double x, y;
}double2D;
typedef struct tag_double3D {
	double x, y, z;
}double3D;
typedef struct tag_double2Ds {
	int n;
	double2D *elem;
}double2Ds;
typedef struct tag_double3Ds {
	int n;
	double3D *elem;
}double3Ds;
typedef struct tag_int3D {
	int x, y, z;
}int3D;
typedef struct tag_int3Ds {
	int n;
	int3D *elem;
}int3Ds;
typedef struct tag_int2D {
	int x, y;
}int2D;
typedef struct tag_int2Ds {
	int n;
	int2D *elem;
}int2Ds;
typedef struct tag_float2D {
	float x, y;
}float2D;
typedef struct tag_float3D {
	float x, y, z;
}float3D;
typedef struct tag_float2Ds {
	int n;
	float2D *elem;
}float2Ds;
typedef struct tag_float3Ds {
	int n;
	float3D *elem;
}float3Ds;
typedef struct tag_Moment {
	double _1, _2, _3, _4, _5;
}Moment;
typedef struct tag_Moments {
	int n;
	Moment *elem;
}Moments;
//IO.cpp
void outputData3D(char *title, double2Ds points, double *value);
void outputData(char *title, double2Ds points);
void inputData3D(char *title, double2Ds points, double *value);
void inputData(char *title, double2Ds points);
//random.cpp
double2D randPosition();
void manyRandPositions(double limits, double2D *a, int n);
void acceptanceRejectionSamplingMethod(double2Ds points_new, double2Ds points_old, double *value_old, int3Ds triangles);
//triangle.cpp
int Triangulate(double2Ds &points, int3Ds &triangles);
void interpolationMoment(double2Ds points_new, Moment *moments_new, double2Ds points_old, Moment *moments_old, int3Ds triangles);
void interpolation(double2Ds points_new, double *value_new, double2Ds points_old, double *value_old, int3Ds triangles);
void calculateAreaOfEveryPoint(double *area, double2Ds points, int3Ds triangles);
//oderk.cpp
void oderk4(Moment *moments, int n, double t, double dt);
double RHS(double M10, double M01, double M11, double M20, double M02, double TM, int MARK);
//bubblepoint.cpp
void bubblePoint(double2Ds points1, double2Ds points2, double *value);
//paint.cpp
void paintTriangles(char *title, double2Ds points, int3Ds triangles);
