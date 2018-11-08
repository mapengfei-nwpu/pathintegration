#include "structs.h"
#include <stdlib.h> 
#include<time.h>
#include<stdio.h>
double randomScalar() {
	return (rand() % RAND_MAX) / (float)(RAND_MAX / 2) - 1.0;
}

double2D randPosition() {

	double2D a;
	a.x = (rand() % RAND_MAX) / (float)(RAND_MAX / 2) - 1.0;
	a.y = (rand() % RAND_MAX) / (float)(RAND_MAX / 2) - 1.0;
	return a;
}
void manyRandPositions(double limits, double2D *a, int n) {
	int i = 0;
	srand((int)time(0));
	for (i = 0; i < n; i++) {
		double2D temp = randPosition();
		temp.x *= limits;
		temp.y *= limits;
		a[i] = temp;
	}
}
void acceptanceRejectionSamplingMethod(double2Ds points_new,double2Ds points_old,double *value_old,int3Ds triangles) {
	int i = 0;
	double c = 0;
	for (int j = 0; j < points_old.n; j++) { if (c < value_old[j]) c = value_old[j]; }
	
	while (i < points_new.n)
	{
		double2Ds p;
		p.n = 1; p.elem = new double2D[1];
		p.elem[0] = randPosition();
		p.elem[0].x *= 5.0;
		p.elem[0].y *= 5.0;
		double r = c*(randomScalar()+1.0)/2.0;
		double b;
		interpolation(p, &b, points_old, value_old, triangles);
		if (r <= b) {
			points_new.elem[i++] = p.elem[0];
		}
		delete[]p.elem;
	}
	printf("%lf\n", c);
}