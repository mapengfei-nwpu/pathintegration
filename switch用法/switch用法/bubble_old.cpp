#include<stdlib.h>
#include <math.h>
#include "structs.h"
#define h (double)0.01
#define TERMS 6			//多项式项数
#define AREAS 3600		//区域个数
#define WIDTH 60		//X轴方向区域数量
#define HEIGHT 60		//Y轴方向区域数量
#define SIDE  0.1		//区域边长
#define STARTX -3		//区域起始点
#define STARTY -3		//区域起始点
typedef struct tagNode {
	int a;
	struct tagNode *next;
} Node, *List;
typedef struct tag_Point	//单个点的数据结构
{
	double x;
	double y;
}Point;
double distance(Point a, Point b);
int filt_one(int n, Point *a, Point b);
Point force(List *list, Point *point, Point *velocity, int n, double *interpolation);
void listCreate(List *list, Point *point, int n, double *interpolation);
void listAdd(List *list, int data);
void listDelete(List *list);
void listReCreate(List *list, Point *point, int n, double *interpolation);
int pointInWhichArea(Point point) {
	int i = ((point.y - STARTY) / SIDE);
	int j = ((point.x - STARTX) / SIDE);
	if (i < 0) i = 0;
	if (i >(WIDTH - 1)) i = WIDTH - 1;
	if (j < 0) j = 0;
	if (j > HEIGHT - 1) j = HEIGHT - 1;
	int area = i * WIDTH + j;
	return area;
}
double* func(Point *points, double *value, int n) {
	int *num = new int[AREAS]();
	double *res = new double[AREAS]();
	for (int i = 0; i < n; i++) {
		int area = pointInWhichArea(points[i]);
		num[area]++;
		res[area] += value[i];
	}
	for (int i = 0; i < AREAS; i++) {
		if (num[i] != 0)res[i] /= num[i];
	}
	return res;
}
double ideal(Point point, double *value) {
	double a = value[pointInWhichArea(point)];
	if (a > 0.5) a = 0.05;
	else a = 0.15 - 0.1*a;
	return a;
}
void range_kutta(Point *point, int n, Point *edge, int m, double *interpolation)
{
	int i, j;
	Point *velocity;
	List *list;

	velocity = (Point*)malloc(n * sizeof(Point));
	list = (List*)malloc(n * sizeof(List));

	listCreate(list, point, n, interpolation);
	for (i = 0; i < n; ++i) {
		velocity[i].x = 0.0;
		velocity[i].y = 0.0;
	}
	Point *k1, *k2, *k3, *k4;
	Point *l1, *l2, *l3, *l4;
	Point *point_1, *point_2, *point_3;
	Point *velocity_1, *velocity_2, *velocity_3;
	Point *point_swap;
	Point *velocity_swap;

	k1 = (Point*)malloc(n * sizeof(Point));
	k2 = (Point*)malloc(n * sizeof(Point));
	k3 = (Point*)malloc(n * sizeof(Point));
	k4 = (Point*)malloc(n * sizeof(Point));

	l1 = (Point*)malloc(n * sizeof(Point));
	l2 = (Point*)malloc(n * sizeof(Point));
	l3 = (Point*)malloc(n * sizeof(Point));
	l4 = (Point*)malloc(n * sizeof(Point));

	point_1 = (Point*)malloc(n * sizeof(Point));
	point_2 = (Point*)malloc(n * sizeof(Point));
	point_3 = (Point*)malloc(n * sizeof(Point));

	velocity_1 = (Point*)malloc(n * sizeof(Point));
	velocity_2 = (Point*)malloc(n * sizeof(Point));
	velocity_3 = (Point*)malloc(n * sizeof(Point));

	point_swap = (Point*)malloc(n * sizeof(Point));
	velocity_swap = (Point*)malloc(n * sizeof(Point));
	for (j = 0; j < 500; ++j) {
		for (i = 0; i < n; ++i) {
			k1[i] = velocity[i];
			l1[i] = force(list, point, velocity, i, interpolation);
			point_1[i].x = point[i].x + h * k1[i].x / (double) 2.0;
			point_1[i].y = point[i].y + h * k1[i].y / (double) 2.0;
			velocity_1[i].x = velocity[i].x + h * l1[i].x / (double) 2.0;
			velocity_1[i].y = velocity[i].y + h * l1[i].y / (double) 2.0;
		}
		for (i = 0; i < n; ++i) {
			k2[i].x = velocity[i].x + h * l1[i].x / (double) 2.0;
			k2[i].y = velocity[i].y + h * l1[i].y / (double) 2.0;
			l2[i] = force(list, point_1, velocity_1, i, interpolation);
			point_2[i].x = point[i].x + h * k2[i].x / (double) 2.0;
			point_2[i].y = point[i].y + h * k2[i].y / (double) 2.0;
			velocity_2[i].x = velocity[i].x + h * l2[i].x / (double) 2.0;
			velocity_2[i].y = velocity[i].y + h * l2[i].y / (double) 2.0;
		}
		for (i = 0; i < n; ++i) {
			k3[i].x = velocity[i].x + h * l2[i].x / (double) 2.0;
			k3[i].y = velocity[i].y + h * l2[i].y / (double) 2.0;
			l3[i] = force(list, point_2, velocity_2, i, interpolation);
			point_3[i].x = point[i].x + h * k3[i].x;
			point_3[i].y = point[i].y + h * k3[i].y;
			velocity_3[i].x = velocity[i].x + h * l3[i].x;
			velocity_3[i].y = velocity[i].y + h * l3[i].y;
		}
		for (i = 0; i < n; ++i) {
			k4[i].x = velocity[i].x + h * l2[i].x / (double) 2.0;
			k4[i].y = velocity[i].y + h * l2[i].y / (double) 2.0;
			l4[i] = force(list, point_3, velocity_3, i, interpolation);
		}
		for (i = 0; i < n; ++i) {
			point_swap[i].x =
				point[i].x +
				h * (k1[i].x + (double) 2.0 * k2[i].x + (double) 2.0 * k3[i].x + k4[i].x) / (double) 6.0;
			point_swap[i].y =
				point[i].y +
				h * (k1[i].y + (double) 2.0 * k2[i].y + (double) 2.0 * k3[i].y + k4[i].y) / (double) 6.0;
			velocity_swap[i].x =
				velocity[i].x +
				h * (l1[i].x + (double) 2.0 * l2[i].x + (double) 2.0 * l3[i].x + l4[i].x) / (double) 6.0;
			velocity_swap[i].y =
				velocity[i].y +
				h * (l1[i].y + (double) 2.0 * l2[i].y + (double) 2.0 * l3[i].y + l4[i].y) / (double) 6.0;
			if (!filt_one(m, edge, point_swap[i])) {
				if (j < 200) {
					double2D temp_point = randPosition();
					point_swap[i].x = temp_point.x*1.0;
					point_swap[i].y = temp_point.y*1.0;
				}
				else point_swap[i] = point[i];
			}
		}
		for (i = 0; i < n; ++i) {
			point[i] = point_swap[i];
			velocity[i] = velocity_swap[i];
		}
		listReCreate(list, point, n, interpolation);
	}
}
double distance(Point a, Point b) {
	double x = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
	double y = sqrt(x);
	return y;
}
Point force(List *list, Point *point, Point *velocity, int n, double *interpolation) {
	Point f;
	f.x = 0.0;
	f.y = 0.0;
	List p = list[n];
	while (p != NULL) {
		double x = point[p->a].x - point[n].x;
		double y = point[p->a].y - point[n].y;
		double ls = sqrt(x * x + y * y);
		double li = (ideal(point[p->a], interpolation) + ideal(point[n], interpolation));
		double r = ls / li;
		double fr;
		if (r < 1.5)
			fr = (double) 1.25 * r * r * r - (double) 2.375 * r * r + (double) 1.125;
		else fr = (double) 0.0;
		double fx = fr * x / ls;
		double fy = fr * y / ls;
		f.x += fx;
		f.y += fy;
		p = p->next;
	}
	f.x = -f.x - (double) 3.8429 * velocity[n].x;
	f.y = -f.y - (double) 3.8429 * velocity[n].y;
	return f;
}
void listCreate(List *list, Point *point, int n, double* interpolation) {
	int i, j;
	for (i = 0; i < n; ++i) {
		list[i] = NULL;
		for (j = 0; j < n; ++j) {
			double x = point[i].x - point[j].x;
			double y = point[i].y - point[j].y;
			double ls = sqrt(x * x + y * y);
			double li = (ideal(point[i], interpolation) + ideal(point[j], interpolation));
			double r = ls / li;
			if (r < 1.5 && ls >0.00001) {
				listAdd(&list[i], j);
			}
		}
	}
}
void listAdd(List *list, int data) {
	List node = (List)malloc(sizeof(Node));
	node->a = data;
	List temp = *list;
	*list = node;
	node->next = temp;
}
void listDelete(List *list) {
	List p = *list;
	List temp;
	while (p != NULL) {
		temp = p->next;
		free(p);
		p = temp;
	}
	*list = NULL;
}
void listReCreate(List *list, Point *point, int n, double *interpolation) {
	int i;
	for (i = 0; i < n; ++i) {
		listDelete(&list[i]);
	}
	listCreate(list, point, n, interpolation);
}
int filt_one(int n, Point *a, Point b) {
	int i, j, c = 0;
	for (i = 0, j = n - 1; i < n; j = i++) {
		if (((a[i].y > b.y) != (a[j].y > b.y)) &&
			(b.x < (a[j].x - a[i].x) * (b.y - a[i].y) / (a[j].y - a[i].y) + a[i].x))
		{
			c = !c;
		}
	}
	return c;
}
void bubblePoint(double2Ds points1, double2Ds points2, double *value) {
	Point *points = new Point[NUM];
	for (int i = 0; i < NUM; i++) {
		points[i].x = points1.elem[i].x;
		points[i].y = points1.elem[i].y;
	}
	double *value2 = func(points, value, NUM);
	int num2 = 4;
	Point edge[4];
	edge[0].x = edge[0].y = edge[1].x = edge[3].y = 3.0;
	edge[1].y = edge[2].x = edge[2].y = edge[3].x = -3.0;
	range_kutta(points, NUM, edge, num2, value2);
	for (int i = 0; i < NUM; i++) {
		points1.elem[i].x = points[i].x;
		points1.elem[i].y = points[i].y;
	}
	delete[]points;  delete value2;
}
