
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include "_head_.h"

float3D randPosition() {

	float3D a;
	a.x = (rand() % RAND_MAX) / (float)(RAND_MAX / 2) - 1.0;
	a.y = (rand() % RAND_MAX) / (float)(RAND_MAX / 2) - 1.0;
	a.z = 0.0;//(rand() % RAND_MAX) / (float)(RAND_MAX / 2) - 1.0;
	return a;
}

void manyRandPositions(float limits, float3D *a, int n) {
	int i = 0;
	srand((int)time(0));
	for (i = 0; i < n; i++) {
		float3D temp = randPosition();
		temp.x *= limits;
		temp.y *= limits;
		temp.z *= limits;
		a[i] = temp;
	}
}