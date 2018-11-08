#include "structs.h"
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace cv;

//»­Í¼º¯Êý
void paintPolygon(Mat img, double2Ds points) {

	Point2d *a = new Point2d[points.n + 1];
	for (size_t i = 0; i < points.n; i++)
	{
		a[i].x = 1500 + 200 * points.elem[i].x; a[i].y = 1500 + 200 * points.elem[i].y;
	}
	a[points.n].x = 1500 + 200 * points.elem[0].x; a[points.n].y = 1500 + 200 * points.elem[0].y;
	for (size_t i = 0; i < points.n; i++)
	{
		line(img, a[i], a[i + 1], Scalar(255, 255, 255), 1, CV_AA, 0);
	}
}
void paintTriangles(char *title, double2Ds points, int3Ds triangles) {

	double2Ds temp_points; temp_points.n = 3;
	temp_points.elem = new double2D[3];
	Mat img(3800, 3800, CV_8UC3, Scalar::all(0));
	for (int i = 0; i < triangles.n; i++)
	{

		temp_points.elem[0] = points.elem[triangles.elem[i].x];
		temp_points.elem[1] = points.elem[triangles.elem[i].y];
		temp_points.elem[2] = points.elem[triangles.elem[i].z];
		paintPolygon(img, temp_points);
	}
	//imshow(title, img);
	imwrite(title, img);
}