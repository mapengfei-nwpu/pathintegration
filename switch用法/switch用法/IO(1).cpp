#include <stdlib.h>
#include <stdio.h>
#include "_head_.h"


void inputData3D(char *title, float3D *positions, int n);
void outputData3D(char *title, float3D *positions, int n);

void inputData3D(char *title, float3D *positions, int n) {
	int i;
	int junk;
	FILE *fp;
	fp = fopen(title, "r");
	for (i = 0; i < n; ++i)
	{
		fscanf(fp, "%d %f %f %f",
			&junk,
			&(positions[i].x),
			&(positions[i].y),
			&(positions[i].z)
		);
	}
	fclose(fp);
}
void outputData3D(char *title, float3D *positions, int n) {
	int i;
	FILE *fp;
	fp = fopen(title, "w");
	for (i = 0; i < n; ++i)
	{
		fprintf(fp, "%d %f %f\n",
			i + 1,
			positions[i].x,
			positions[i].y
		);
	}
	fclose(fp);
}