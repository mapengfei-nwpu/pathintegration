
#include"structs.h"
#include<stdio.h>
#include<cmath>
int CircumCircle(double xp, double yp, double x1, double y1, double x2,
	double y2, double x3, double y3, double &xc, double &yc, double &r) {
	double m1, m2, mx1, mx2, my1, my2;
	double dx, dy, rsqr, drsqr;
	if (abs(y1 - y2) < EPSILON && abs(y2 - y3) < EPSILON)
		return(false);
	if (abs(y2 - y1) < EPSILON) {
		m2 = -(x3 - x2) / (y3 - y2);
		mx2 = (x2 + x3) / 2.0;
		my2 = (y2 + y3) / 2.0;
		xc = (x2 + x1) / 2.0;
		yc = m2 * (xc - mx2) + my2;
	}
	else if (abs(y3 - y2) < EPSILON) {
		m1 = -(x2 - x1) / (y2 - y1);
		mx1 = (x1 + x2) / 2.0;
		my1 = (y1 + y2) / 2.0;
		xc = (x3 + x2) / 2.0;
		yc = m1 * (xc - mx1) + my1;
	}
	else {
		m1 = -(x2 - x1) / (y2 - y1);
		m2 = -(x3 - x2) / (y3 - y2);
		mx1 = (x1 + x2) / 2.0;
		mx2 = (x2 + x3) / 2.0;
		my1 = (y1 + y2) / 2.0;
		my2 = (y2 + y3) / 2.0;
		xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
		yc = m1 * (xc - mx1) + my1;
	}
	dx = x2 - xc;
	dy = y2 - yc;
	rsqr = dx * dx + dy * dy;
	r = sqrt(rsqr);
	dx = xp - xc;
	dy = yp - yc;
	drsqr = dx * dx + dy * dy;
	return((drsqr <= rsqr) ? true : false);
}
int Triangulate(double2Ds &points, int3Ds &triangles) {

	int ntri = 0;
	int nv = points.n;
	int3D *v = new int3D[3 * nv];
	double2D *p = new double2D[nv + 3];							//the origin points and three points of the big triangle
	for (int i = 0; i < nv; i++) p[i] = points.elem[i];			//do not need origin data anymore！

	int *complete = NULL;
	int2D *edges = NULL;
	int2D *p_EdgeTemp;
	int nedge = 0;
	int trimax, emax = 200;
	int status = 0;
	int inside;
	int i, j, k;
	double xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r;
	double xmin, xmax, ymin, ymax, xmid, ymid;
	double dx, dy, dmax;

	/* Allocate memory for the completeness list, fag for each triangle */
	trimax = 4 * nv;
	complete = new int[trimax];
	/* Allocate memory for the edge list */
	edges = new int2D[emax];
	/*
	Find the maximum and minimum vertex bounds.
	This is to allow calculation of the bounding triangle
	*/
	xmin = p[0].x;
	ymin = p[0].y;
	xmax = xmin;
	ymax = ymin;
	for (i = 1; i < nv; i++) {
		if (p[i].x < xmin) xmin = p[i].x;
		if (p[i].x > xmax) xmax = p[i].x;
		if (p[i].y < ymin) ymin = p[i].y;
		if (p[i].y > ymax) ymax = p[i].y;
	}
	dx = xmax - xmin;
	dy = ymax - ymin;
	dmax = (dx > dy) ? dx : dy;
	xmid = (xmax + xmin) / 2.0;
	ymid = (ymax + ymin) / 2.0;
	/*
	Set up the supertriangle
	his is a triangle which encompasses all the sample points.
	The supertriangle coordinates are added to the end of the
	vertex list. The supertriangle is the first triangle in
	the triangle list.
	*/
	p[nv + 0].x = xmid - 20 * dmax;
	p[nv + 0].y = ymid - dmax;
	p[nv + 1].x = xmid;
	p[nv + 1].y = ymid + 20 * dmax;
	p[nv + 2].x = xmid + 20 * dmax;
	p[nv + 2].y = ymid - dmax;
	v[0].x = nv;
	v[0].y = nv + 1;
	v[0].z = nv + 2;
	complete[0] = false;
	ntri = 1;
	for (i = 0; i < nv; i++) {
		xp = p[i].x;
		yp = p[i].y;
		nedge = 0;
		for (j = 0; j < ntri; j++) {
			if (complete[j])
				continue;
			x1 = p[v[j].x].x;
			y1 = p[v[j].x].y;
			x2 = p[v[j].y].x;
			y2 = p[v[j].y].y;
			x3 = p[v[j].z].x;
			y3 = p[v[j].z].y;
			inside = CircumCircle(xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r);
			if (xc + r + EPSILON < xp)
				// Suggested
				// if (xc + r + EPSILON < xp)
				complete[j] = true;
			if (inside) {
				/* Check that we haven't exceeded the edge list size */
				if (nedge + 3 >= emax) {
					emax += 100;
					p_EdgeTemp = new int2D[emax];
					for (int i = 0; i < nedge; i++) { // Fix by John Bowman
						p_EdgeTemp[i] = edges[i];
					}
					delete[]edges;
					edges = p_EdgeTemp;
				}
				edges[nedge + 0].x = v[j].x;
				edges[nedge + 0].y = v[j].y;
				edges[nedge + 1].x = v[j].y;
				edges[nedge + 1].y = v[j].z;
				edges[nedge + 2].x = v[j].z;
				edges[nedge + 2].y = v[j].x;
				nedge += 3;
				v[j] = v[ntri - 1];
				complete[j] = complete[ntri - 1];
				ntri--;
				j--;
			}
		}
		for (j = 0; j < nedge - 1; j++) {
			for (k = j + 1; k < nedge; k++) {
				if ((edges[j].x == edges[k].y) && (edges[j].y == edges[k].x)) {
					edges[j].x = -1;
					edges[j].y = -1;
					edges[k].x = -1;
					edges[k].y = -1;
				}
				/* Shouldn't need the following, see note above */
				if ((edges[j].x == edges[k].x) && (edges[j].y == edges[k].y)) {
					edges[j].x = -1;
					edges[j].y = -1;
					edges[k].x = -1;
					edges[k].y = -1;
				}
			}
		}
		for (j = 0; j < nedge; j++) {
			if (edges[j].x < 0 || edges[j].y < 0)
				continue;
			v[ntri].x = edges[j].x;
			v[ntri].y = edges[j].y;
			v[ntri].z = i;
			complete[ntri] = false;
			ntri++;
		}
	}
	for (i = 0; i < ntri; i++) {
		if (v[i].x >= nv || v[i].y >= nv || v[i].z >= nv) {
			v[i] = v[ntri - 1];
			ntri--;
			i--;
		}
	}
	triangles.n = ntri;
	triangles.elem = new int3D[ntri];
	for (int i = 0; i < ntri; i++) triangles.elem[i] = v[i];
	delete[] v;
	delete[] p;
	delete[] edges;
	delete[] complete;
	return 0;
}
double triangleArea(double2D *points) {
	double x1 = points[0].x;
	double y1 = points[0].y;
	double x2 = points[1].x;
	double y2 = points[1].y;
	double x3 = points[2].x;
	double y3 = points[2].y;
	return abs(0.5*(x1*y2 + x2*y3 + x3*y1 - x1*y3 - x2*y1 - x3*y2));
}
int isPointInTriangle(double2D point, double2D *triangle) {
	double2D triangle1[3];
	double2D triangle2[3];
	double2D triangle3[3];
	triangle1[0] = triangle2[0] = triangle3[0] = point;
	triangle1[2] = triangle3[2] = triangle[0];
	triangle1[1] = triangle2[2] = triangle[2];
	triangle2[1] = triangle3[1] = triangle[1];
	if (abs(triangleArea(triangle) - triangleArea(triangle1) - triangleArea(triangle2) - triangleArea(triangle3))<EPSILON)
		return 1;
	else return 0;
}
int3D findTriangleContainsPoint(double2Ds points, int3Ds triangles, double2D point) {
	double2D temp[3];
	for (int i = 0; i < triangles.n; i++) {
		temp[0] = points.elem[triangles.elem[i].x];
		temp[1] = points.elem[triangles.elem[i].y];
		temp[2] = points.elem[triangles.elem[i].z];
		if (isPointInTriangle(point, temp))
			return triangles.elem[i];
	}
	int3D wrong = { -1,-1,-1 };
	return wrong;
}
double triangleAreaCoordinateInterpolation(double2D p, double2D *points, double *value) {
	//三角形面积坐标插值
	//具体推导：https://en.wikipedia.org/wiki/Barycentric_coordinate_system
	double2D p1 = points[0]; double a1 = value[0];
	double2D p2 = points[1]; double a2 = value[1];
	double2D p3 = points[2]; double a3 = value[2];
	//p=(1-c1-c2)p1+c1*p2+c2*p3;
	//p=p1+c1(p2-p1)+c2(p3-p1);
	double x1 = p2.x - p1.x; double y1 = p2.y - p1.y;
	double x2 = p3.x - p1.x; double y2 = p3.y - p1.y;
	double x = p.x - p1.x; double y = p.y - p1.y;
	double c1, c2 = (x*y1 - y*x1) / (x2*y1 - x1*y2);
	if (abs(x1) > EPSILON) c1 = (x - c2*x2) / x1;
	else c1 = (y - c2*y2) / y1;
	//返回函数值
	double ret = (1 - c1 - c2) * a1 + c1 * a2 + c2 * a3;
	return ret;
	//return (a1+a2+a3)/3.0;
}
void interpolation(double2Ds points_new, double *value_new, double2Ds points_old, double *value_old, int3Ds triangles)
{
	for (size_t i = 0; i < points_new.n; i++)
	{
		int3D temp_triangle;
		//find which triangle contains the points_old[i] at first. 
		temp_triangle = findTriangleContainsPoint(points_old, triangles, points_new.elem[i]);
		//if the point is not in any triangle
		if (temp_triangle.x == -1) {
			value_new[i] = 0.0;
			continue;
		}
		//the temp_triangle have the index of the three points which make up the triangle.
		//the index are of points_old.
		double2D temp_points[3];
		temp_points[0] = points_old.elem[temp_triangle.x];
		temp_points[1] = points_old.elem[temp_triangle.y];
		temp_points[2] = points_old.elem[temp_triangle.z];
		double2D temp_point = points_new.elem[i];
		double temp_values[3];
		temp_values[0] = value_old[temp_triangle.x];
		temp_values[1] = value_old[temp_triangle.y];
		temp_values[2] = value_old[temp_triangle.z];
		value_new[i] = triangleAreaCoordinateInterpolation(points_new.elem[i], temp_points, temp_values);
	}
}
void interpolationMoment(double2Ds points_new, Moment *moments_new, double2Ds points_old, Moment *moments_old, int3Ds triangles) {
	for (size_t i = 0; i < points_new.n; i++)
	{
		int3D temp_triangle;
		//find which triangle contains the points_old[i] at first. 
		temp_triangle = findTriangleContainsPoint(points_old, triangles, points_new.elem[i]);
		//if the point is not in any triangle
		if (temp_triangle.x == -1) {
			continue;
		}
		//the temp_triangle have the index of the three points which make up the triangle.
		//the index are of points_old.
		double2D temp_points[3];
		temp_points[0] = points_old.elem[temp_triangle.x];
		temp_points[1] = points_old.elem[temp_triangle.y];
		temp_points[2] = points_old.elem[temp_triangle.z];
		double2D temp_point = points_new.elem[i];
		double temp_values[3];
		temp_values[0] = moments_old[temp_triangle.x]._1;
		temp_values[1] = moments_old[temp_triangle.y]._1;
		temp_values[2] = moments_old[temp_triangle.z]._1;
		moments_new[i]._1 = triangleAreaCoordinateInterpolation(points_new.elem[i], temp_points, temp_values);
		temp_values[0] = moments_old[temp_triangle.x]._2;
		temp_values[1] = moments_old[temp_triangle.y]._2;
		temp_values[2] = moments_old[temp_triangle.z]._2;
		moments_new[i]._2 = triangleAreaCoordinateInterpolation(points_new.elem[i], temp_points, temp_values);
		temp_values[0] = moments_old[temp_triangle.x]._3;
		temp_values[1] = moments_old[temp_triangle.y]._3;
		temp_values[2] = moments_old[temp_triangle.z]._3;
		moments_new[i]._3 = triangleAreaCoordinateInterpolation(points_new.elem[i], temp_points, temp_values);
		temp_values[0] = moments_old[temp_triangle.x]._4;
		temp_values[1] = moments_old[temp_triangle.y]._4;
		temp_values[2] = moments_old[temp_triangle.z]._4;
		moments_new[i]._4 = triangleAreaCoordinateInterpolation(points_new.elem[i], temp_points, temp_values);
		temp_values[0] = moments_old[temp_triangle.x]._5;
		temp_values[1] = moments_old[temp_triangle.y]._5;
		temp_values[2] = moments_old[temp_triangle.z]._5;
		moments_new[i]._5 = triangleAreaCoordinateInterpolation(points_new.elem[i], temp_points, temp_values);
	}
}
double2D pedalOfPointToLine(double x1, double y1, double x2, double y2, double xc, double yc) {
	double2D p;
	double dxc1 = xc - x1;
	double dyc1 = yc - y1;
	double dx12 = x1 - x2;
	double dy12 = y1 - y2;
	double sqxy = dx12*dx12 + dy12*dy12;
	p.x = x1 + dx12*dy12 / sqxy*dyc1 + dx12*dx12 / sqxy*dxc1;
	p.y = y1 + dy12*dy12 / sqxy*dyc1 + dy12*dx12 / sqxy*dxc1;
	//p.x = (x1 + x2) / 2.0;
	//p.y = (y1 + y2) / 2.0;
	return p;
}
double3D calculateThreePartsAreaOfATriangle(double2D *triangle) {
	//将一个三角形分成三个多边形，每个四边形所占的面积
	double junk1, junk2, junk3, junk4;
	junk1 = junk2 = junk3 = junk4 = 0;
	double2D p1, p2, p3, pc;
	double xc, yc;
	double3D s;
	double2D points[3];

	CircumCircle(junk1, junk2,
		triangle[0].x, triangle[0].y,
		triangle[1].x, triangle[1].y,
		triangle[2].x, triangle[2].y,
		xc, yc, junk3);
	pc.x = xc; pc.y = yc;
	if (isPointInTriangle(pc, triangle)) {
		p1 = pedalOfPointToLine(triangle[0].x, triangle[0].y, triangle[1].x, triangle[1].y, xc, yc);
		p2 = pedalOfPointToLine(triangle[1].x, triangle[1].y, triangle[2].x, triangle[2].y, xc, yc);
		p3 = pedalOfPointToLine(triangle[2].x, triangle[2].y, triangle[0].x, triangle[0].y, xc, yc);
		points[1].x = xc; points[1].y = yc;
		points[0] = triangle[0];
		points[2] = p1;
		double a = triangleArea(points);
		points[2] = p3;
		a += triangleArea(points);
		points[0] = triangle[1];
		points[2] = p1;
		double b = triangleArea(points);
		points[2] = p2;
		b += triangleArea(points);
		points[0] = triangle[2];
		points[2] = p2;
		double c = triangleArea(points);
		points[2] = p3;
		c += triangleArea(points);
		s.x = a;
		s.y = b;
		s.z = c;
		s.x = s.y = s.z = triangleArea(triangle) / 3;
		return s;
	}
	else {
		s.x = s.y = s.z = triangleArea(triangle) / 3;
		return s;
	}
}
void calculateAreaOfEveryPoint(double *area, double2Ds points, int3Ds triangles) {
	//计算所有三角形的面积
	double sum = 0.0;
	for (int i = 0; i < triangles.n; i++)
	{
		double2D temp[3];
		double3D S;
		temp[0] = points.elem[triangles.elem[i].x];
		temp[1] = points.elem[triangles.elem[i].y];
		temp[2] = points.elem[triangles.elem[i].z];
		S = calculateThreePartsAreaOfATriangle(temp);
		area[triangles.elem[i].x] += S.x;
		area[triangles.elem[i].y] += S.y;
		area[triangles.elem[i].z] += S.z;
		double d = triangleArea(temp) - S.x - S.y - S.z;
	}

}