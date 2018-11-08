#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>

typedef struct tag_float2D {
	float x, y;
}float2D;

typedef struct tag_Particle {
	float2D position;
	float2D velocity;
	float2D force;
	float   overlap;
	float   ideal; //理想距离
	int     tag;
}Particle;

typedef struct tag_Particles {
	Particle *particle;
	float    *dist;
	int      *index;
	int		 n;
}Particles;


//#include"cuda_knn.cuh"
bool knn_cuda_texture(const float * ref,
	int           ref_nb,
	const float * query,
	int           query_nb,
	int           dim,
	int           k,
	float *       knn_dist,
	int *         knn_index);

double TIME_STEP = 0.2;
double limits = 5;

// points placed on edges.
float2D *fixed_points;
int      fixed_points_nb;


void updateParticles(Particles particles, float(*ideal)(float,float));
Particle updateForce(int one_i, Particles b);
Particle updateVelocity(Particle particle);
Particle updatePosition(Particle particle);
float2D oneForce(Particle a, Particle b);
float molecularForce(float distance);
int adjustParticleNumber(Particles p_old, Particles p_new, float(*ideal)(float, float));

// euclidian distance.
float dist(float2D a, float2D b) {
	return sqrtf((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

// 
void updateIdeal(Particles particles, float(*ideal)(float, float)) {
	for (int i = 0; i < particles.n; i++)
		particles.particle[i].ideal = ideal(particles.particle[i].position.x,particles.particle[i].position.y);
}
void cudaNeighborSearch(Particles particles) {
	float *ref = (float*)malloc(2 * particles.n * sizeof(float));
	for (size_t i = 0; i < particles.n; i++)
	{
		ref[i] = particles.particle[i].position.x;
		ref[i + particles.n] = particles.particle[i].position.y;
	}
	knn_cuda_texture(ref, particles.n, ref, particles.n, 2, 20, particles.dist, particles.index);
	free(ref);
}
void updateParticles(Particles particles, float(*ideal)(float, float)) {
	updateIdeal(particles, ideal);
	cudaNeighborSearch(particles);
	for (int i = fixed_points_nb; i < particles.n; i++)
	{
		float2D  tempPosition = particles.particle[i].position;
		particles.particle[i] = updateForce(i, particles);
		particles.particle[i] = updateVelocity(particles.particle[i]);
		particles.particle[i] = updatePosition(particles.particle[i]);
		//边界控制
		if (particles.particle[i].position.x >  limits) particles.particle[i].position.x = limits;
		if (particles.particle[i].position.x < -limits) particles.particle[i].position.x = -limits;
		if (particles.particle[i].position.y >  limits) particles.particle[i].position.y = limits;
		if (particles.particle[i].position.y < -limits) particles.particle[i].position.y = -limits;
	}
}

Particle updateForce(int one_i, Particles particles) {
	int i = 0;
	Particle one = particles.particle[one_i];
	one.force = { 0,0};
	for (i = 0; i < 20; i++)
	{
		if (one_i != particles.index[i*particles.n + one_i]) {
			float2D tempForce = oneForce(one, particles.particle[particles.index[i*particles.n + one_i]]);
			one.force.x += tempForce.x;
			one.force.y += tempForce.y;
		}
	}
	return one;
}
Particle updateVelocity(Particle particle) {
	float m = 1;
	float c = 3.8429;
	float2D v = particle.velocity;
	float2D f = particle.force;
	particle.velocity.x = (f.x - c * v.x) / m*TIME_STEP + v.x;
	particle.velocity.y = (f.y - c * v.y) / m*TIME_STEP + v.y;
	return particle;
}
Particle updatePosition(Particle particle) {
	particle.position.x = particle.ideal*particle.velocity.x * TIME_STEP + particle.position.x;
	particle.position.y = particle.ideal*particle.velocity.y * TIME_STEP + particle.position.y;
	return particle;
}
float calculateEnergy(Particles particles) {
	float energy = (float)0;
	for (size_t i = 0; i < particles.n; i++)
	{
		float2D temp = particles.particle[i].velocity;
		energy += sqrtf(temp.x*temp.x + temp.y*temp.y );
	}
	return energy / particles.n;
}
void updateOverlap(Particles particles) {
	for (size_t i = 0; i < particles.n; i++)
	{
		particles.particle[i].overlap = (float)0.0;
		float ri = particles.particle[i].ideal;
		for (size_t j = 0; j < particles.n; j++)
		{
			float rj = particles.particle[j].ideal;
			float ij = dist(particles.particle[i].position, particles.particle[j].position);
			float temp = (float)2 * ri + rj - ij;
			if (temp>0) particles.particle[i].overlap += temp;
		}
		particles.particle[i].overlap /= ri;
		particles.particle[i].overlap -= 3;
	}
}
void initialParticles(Particles &p, int n) {
	p.particle = (Particle*)malloc(n * sizeof(Particle));
	p.dist = (float*)malloc(20 * n * sizeof(float));
	p.index = (int*)malloc(20 * n * sizeof(int));
	p.n = n;
}
void finishParticles(Particles &p) {
	free(p.particle);
	free(p.index);
	free(p.dist);
	p.particle = NULL;
	p.index = NULL;
	p.dist = NULL;
	p.n = 0;
}
float2D oneForce(Particle a, Particle b) {
	float2D force;
	float x = a.position.x - b.position.x;
	float y = a.position.y - b.position.y;
	float ls = dist(a.position, b.position);
	float li = a.ideal + b.ideal;
	float r = ls / li;
	float absForce = molecularForce(r);
	force.x = absForce * x / ls;
	force.y = absForce * y / ls;
	//avoid the repeat of particles.
	if (ls < 0.00001) force = { 10,10};
	return force;
}
float molecularForce(float distance) {
	//calculate the van der waals force with the distance of two particles
	float force;
	if (distance < 1.5)
		force = (float) 1.25 * distance * distance * distance
		- (float) 2.375 * distance * distance + (float) 1.125;
	else force = (float) 0.0;
	return force;
}
int adjustParticleNumber(Particles p_old, Particles p_new, float(*ideal)(float, float)) {
	// delete: -1, unknow: 0, reserve: 1, double: 2.
	updateIdeal(p_old, ideal);
	int *tag = (int*)calloc(p_old.n, sizeof(int));
	updateOverlap(p_old);
	// delete particles who are near the edges.
	for (int i = fixed_points_nb; i < p_old.n; i++) {
		double r = p_old.particle[i].ideal*0.5;
		double a = fabs(fabs(p_old.particle[i].position.x) - limits);
		double b = fabs(fabs(p_old.particle[i].position.y) - limits);
		if (b < a)a = b;//get the smaller one.
		if (a < r) {
			if (tag[i] != 1) tag[i] = -1;
		}
	}
	for (int i = fixed_points_nb; i < p_old.n; i++) {
		if (tag[i] == -1) continue;
		if (p_old.particle[i].overlap > 8) {

			if (tag[i] != 1) {
				tag[i] = -1; //delete
				double ri = p_old.particle[i].ideal;
				for (int j = 0; j < p_old.n; j++){

					double rj = p_old.particle[j].ideal;
					double r = 1.2*(ri + rj);
					double l = dist(p_old.particle[i].position, p_old.particle[j].position);
					if (l < r && tag[j] != -1) tag[j] = 1;
				}
			}
		}
	}
	// double the point . 
	for (int i = fixed_points_nb; i < p_old.n; i++)
		if ((fabs(p_old.particle[i].position.x) != limits &&
			fabs(p_old.particle[i].position.y) != limits) &&
			p_old.particle[i].overlap < 5.5) tag[i] = 2;
	// transfer data .
	int count1 = 0, count2 = 0;
	while (count1<p_old.n && count2 < p_new.n) {
		if (tag[count1] == 0 || tag[count1] == 1)
			p_new.particle[count2++] = p_old.particle[count1];
		if (tag[count1] == 2) {
			p_new.particle[count2++] = p_old.particle[count1];
			p_new.particle[count2++] = p_old.particle[count1];
			p_new.particle[count2 - 1].velocity = { 0,0 };
		}
		++count1;
	}
	free(tag);
	return count2;
}
int countPointOnOneEdge(float2D a, float2D b, float(*ideal)(float, float)) {
	float2D c;
	c.x = 0.5*(a.x + b.x);
	c.y = 0.5*(a.y + b.y);
	if (ideal(a.x, a.y) + ideal(b.x, b.y) < 0.75* dist(a, b)) {
		return countPointOnOneEdge(a, c, ideal) + countPointOnOneEdge(b, c, ideal) - 1;
	}
	else return 2;
}
int countPointOnEdges(float2D *apex, int apex_nb, float(*ideal)(float, float)) {
	int n = 0;
	for (int i = 0, j = apex_nb - 1; i < apex_nb; j = i++)
		n += countPointOnOneEdge(apex[i], apex[j], ideal);
	return n - apex_nb;
}
//when called outside this function, i should zero .
int ith = 0;
void placePointOnOneEdge(float2D a, float2D b, float(*ideal)(float, float)) {
	float2D c;
	c.x = 0.5*(a.x + b.x);
	c.y = 0.5*(a.y + b.y);
	if (ideal(a.x,a.y) + ideal(b.x,b.y) < 0.75*dist(a, b)) {
		fixed_points[ith++] = c;
		placePointOnOneEdge(a, c, ideal);
		placePointOnOneEdge(c, b, ideal);
	}
}
void placePointOnEdges(float2D *apex, int apex_nb, float(*ideal)(float, float)) {
	if (ith != 0) ith = 0;
	for (int i = 0, j = apex_nb - 1; i < apex_nb; j = i++) {
		fixed_points[ith++] = apex[i];
		placePointOnOneEdge(apex[i], apex[j], ideal);
	}
	
	//renew the global variable.
	ith = 0;
}


double* bubblePoint(double l,float(*ideal)(float ,float),int &n) {

	// limit the area of particles .
	limits = l;

	// alloc the memory .
	Particles particles;
	Particles particles_new;
	initialParticles(particles, 900);
	initialParticles(particles_new, 20000);


	// edges
	float2D apex[4];
	apex[0] = { (float)limits, (float)limits};
	apex[1] = { (float)limits,-(float)limits};
	apex[2] = { -(float)limits,-(float)limits};
	apex[3] = { -(float)limits, (float)limits};

	// on edges ...
	fixed_points_nb = countPointOnEdges(apex, 4, ideal);
	fixed_points = (float2D*)malloc(fixed_points_nb * sizeof(float2D));
	placePointOnEdges(apex, 4, ideal);
	
	// initial data 
	float2D a[900];
	for (size_t i = 0; i < 900; i++) {
		a[i].x = ((rand() % RAND_MAX) / (float)(RAND_MAX / 2) - 1.0)*limits;
		a[i].y = ((rand() % RAND_MAX) / (float)(RAND_MAX / 2) - 1.0)*limits;
	}
	for (size_t i = 0; i < fixed_points_nb; i++) a[i] = fixed_points[i];
	for (size_t i = 0; i < 900; i++){
		particles.particle[i].ideal = 0;
		particles.particle[i].force = { 0,0};
		particles.particle[i].velocity = { 0,0 };
		particles.particle[i].position = a[i];
	}

	// the main logical part !

	// number of particles after simulations .
	int new_nb = 0;
	int new_nbs[5] = { 0 };

	//outer and inner loop times .
	int outer = 100;
	int inner = 40;

	//pre-simulation .
	for (size_t j = 0; j < inner; j++)updateParticles(particles, ideal);
	printf("%d\n",fixed_points_nb);
	//real start .
	for (size_t i = 0; i < outer; i++) {

		//adjust particle number.
		updateIdeal(particles, ideal);
		new_nb = adjustParticleNumber(particles, particles_new, ideal);
		finishParticles(particles);
		initialParticles(particles, new_nb);
		for (size_t j = 0; j < new_nb; j++)	particles.particle[j] = particles_new.particle[j];

		//simulation.
		for (size_t j = 0; j < inner; j++)	updateParticles(particles, ideal);
		
		//judge the statement.
		int    good_particles_nb   = 0;
		double good_particles_rate = 0;
		updateOverlap(particles);
		for (size_t j = fixed_points_nb; j < particles.n; j++) {
			if (particles.particle[j].overlap > 5.5 && particles.particle[j].overlap < 8.0)
				++good_particles_nb;
		}
		good_particles_rate = (double)good_particles_nb / (double)(new_nb-fixed_points_nb);
		
		// judge how the new_nb changed. 
		double new_nb_average = 0;
		new_nbs[i % 5] = new_nb;
		for (size_t j = 0; j < 5; j++) new_nb_average += ((double)new_nbs[j] / 5.0);

		// stop conditions.
		printf("%lf  %lf\n", good_particles_rate, fabs(new_nb_average - new_nb));
		if (good_particles_rate > 0.99 && fabs(new_nb_average-new_nb) < 3)break;
		
	}

	// the last simulation
	for (size_t j = 0; j < 100; j++) updateParticles(particles, ideal);

	// transform data structure .
	double * data = (double*)malloc(2 * particles.n * sizeof(double));
	for (size_t j = 0; j < particles.n; j++) {
		data[j] = particles.particle[j].position.x;
		data[j + particles.n] = particles.particle[j].position.y;
	}
	n = particles.n;
	return data;
}
