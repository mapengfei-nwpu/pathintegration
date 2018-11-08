#pragma once
#pragma once
typedef struct tag_float3D {
	float x, y;
}float3D;
typedef struct tag_float2D {
	float x, y;
}float2D;
typedef struct tag_double2D {
	double x, y;
}double2D;

typedef struct tag_Particle {
	float3D position;
	float3D velocity;
	float3D force;
	float   overlap;
	float   ideal; //ÀíÏë¾àÀë
	int     tag;
}Particle;

typedef struct tag_Particles {
	Particle *particle;
	float    *dist;
	int      *index;
	int		 n;
}Particles;