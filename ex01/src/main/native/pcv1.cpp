//============================================================================
// Name        : pcv1.cpp
// Author      : Ronny Haensch
// Version     : 1.0
// Copyright   : -
// Description : Load and display an image; do some basic linear algebra.
// Group       : Marcus Grum (), Robin Vobruba (343773), Jens Jawer () und Marcus Pannwitz (343479)
//============================================================================

#include <iostream>
#ifdef WIN32
	#include "stdafx.h"
#endif
#include "opencv2/opencv.hpp"

using namespace std;
using namespace cv;

// some function headers (see below)
bool isPointOnLine(Mat& point, Mat& line, double eps = pow((double)10, (double)-5));
void applyH(Mat& geomObj, Mat& H, Mat& result, string type);
void getH(Mat& T, Mat& R, Mat& S, Mat& H);
void getScaleMatrix(double lambda, Mat& S);
void getRotMatrix(double phi, Mat& R);
void getTranslMatrix(double dx, double dy, Mat& T);
void getConnectingLine(Mat& p1, Mat& p2, Mat& line);

/**
 * usage: path to image in argv[1]
 */
int main(int argc, char** argv) {

	// check if image path was defined
	if (argc != 2) {
		cerr << "Usage: pcv1 <path_to_image>" << endl;
		return -1;
	}

	// window names
	string win1 = string ("Image");

	// load image as gray-scale, path in argv[1]
	cout << "Load image: start" << endl;
	Mat inputImage;

	inputImage = imread( argv[1]);
	if (!inputImage.data) {
		cerr << "ERROR: image could not be loaded from " << argv[1] << endl;
	} else {
		cout << "Load image: done ( " << inputImage.rows << " x " << inputImage.cols << " )" << endl;
	}

	// show input image
	namedWindow( win1.c_str(), CV_WINDOW_AUTOSIZE );
	imshow( win1.c_str(), inputImage );
	waitKey(0);

	// ----------------------------------------

	// the two given points as OpenCV matrices
	Mat x(2, 1, CV_32FC1);
	x.at<float>(0, 0) = 2;
	x.at<float>(1, 0) = 3;

	Mat y(2, 1, CV_32FC1);
	y.at<float>(0, 0) = -4;
	y.at<float>(1, 0) = 5;

	// same points in homogeneous coordinates
	Mat v1(3, 1, CV_32FC1);
	v1.at<float>(0, 0) = x.at<float>(0, 0);
	v1.at<float>(1, 0) = x.at<float>(1, 0);
	v1.at<float>(2, 0) = 1;
	Mat v2(3, 1, CV_32FC1);
	// define v2 as homogeneous version of y
	v2.at<float>(0, 0) = y.at<float>(0, 0);
	v2.at<float>(1, 0) = y.at<float>(1, 0);
	v2.at<float>(2, 0) = 1;

	// print points
	cout << "point 1: " << v1.t() << "^T" << endl;
	cout << "point 2: " << v2.t() << "^T" << endl;
	cout << endl;

	// the connecting line between those points in homogeneous coordinates
	Mat line;
	getConnectingLine(v1, v2, line);

	// print line
	cout << "joining line: " << line << "^T" << endl;
	cout << endl;

	// the parameters of the transformation
	int dx = 6; // translation in x
	int dy = -7; // translation in y
	double phi = 15; // rotation angle in degree
	double lambda = 8; // scaling factor

	// matrices for transformation
	Mat T, R, S, H; // translation, rotation, scaling, homography
	// calculate translation matrix
	getTranslMatrix(dx, dy, T);
	// calculate rotation matrix
	getRotMatrix(phi, R);
	// calculate scale matrix
	getScaleMatrix(lambda, S);
	// combine individual transformations to a homography
	getH(T, R, S, H);

	// print calculated matrices
	cout << "Translation matrix: " << endl;
	cout << T << endl;
	cout << endl;
	cout << "Rotation matrix: " << endl;
	cout << R << endl;
	cout << endl;
	cout << "Scaling matrix: " << endl;
	cout << S << endl;
	cout << endl;
	cout << "Homography: " << endl;
	cout << H << endl;
	cout << endl;

	// transform first point x (and print it)
	Mat v1_new;
	applyH(v1, H, v1_new, "point");
	cout << "new point 1: " << v1_new << "^T" << endl;

	// transform second point y (and print it)
	Mat v2_new;
	applyH(v2, H, v2_new, "point");
	cout << "new point 2: " << v2_new << "^T" << endl;
	cout << endl;

	// transform joining line (and print it)
	Mat line_new;
	applyH(line, H, line_new, "line");
	cout << "new line: " << line_new << "^T" << endl;
	cout << endl;

	// check if transformed points are still on transformed line
	bool xOnLine = isPointOnLine(v1_new, line_new);
	bool yOnLine = isPointOnLine(v2_new, line_new);
	if (xOnLine) {
		cout << "first point lies still on the line *yay*" << endl;
	} else {
		cout << "first point does not lie on the line *oh oh*" << endl;
	}

	if (yOnLine) {
		cout << "second point lies still on the line *yay*" << endl;
	} else {
		cout << "second point does not lie on the line *oh oh*" << endl;
	}

	waitKey(0);

	return 0;
}

/**
 * calculates joining line between two points
 * @param p1  the first point in homogeneous coordinates
 * @param p2  the second point in homogeneous coordinates
 * @param line  the joining line in homogeneous coordinates
 */
void getConnectingLine(Mat& p1, Mat& p2, Mat& line) {

	// init line
	line = Mat::zeros(3, 1, CV_32FC1);

	// Calculate the cross product to get the line, which intersects both points
	line.at<float>(0, 0) = p1.at<float>(1, 0) * p2.at<float>(2, 0) - p1.at<float>(2, 0) * p2.at<float>(1, 0);
	line.at<float>(1, 0) = p1.at<float>(2, 0) * p2.at<float>(0, 0) - p1.at<float>(0, 0) * p2.at<float>(2, 0);
	line.at<float>(2, 0) = p1.at<float>(0, 0) * p2.at<float>(1, 0) - p1.at<float>(1, 0) * p2.at<float>(0, 0);
}

/**
 * generates translation matrix T defined by translation (dx, dy)^T
 * @param dx  the translation in x-direction
 * @param dy  the translation in y-direction
 * @param T  the resulting translation matrix
 */
void getTranslMatrix(double dx, double dy, Mat& T) {

	// init T as identity matrix
	T = Mat::eye(3, 3, CV_32FC1);

	// The translation is represented in the third column
	T.at<float>(0, 2) = dx;
	T.at<float>(1, 2) = dy;
}

/**
 * generates rotation matrix R defined by angle phi
 * @param phi  the rotation angle in degree
 * @param R  the resulting rotation matrix
 */
void getRotMatrix(double phi, Mat& R) {

	// init T as identity matrix
	R = Mat::eye(3, 3, CV_32FC1);

	// transform degree to radian
	phi = phi / 180 * CV_PI;

	// Calculate the rotation values
	R.at<float>(0, 0) = R.at<float>(1, 1) = (float)cos(phi);
	R.at<float>(1, 0) = (float)sin(phi);
	R.at<float>(0, 1) = -R.at<float>(1, 0);
}

/**
 * generates scaling matrix S defined by scaling factor lambda
 * @param lambda the scaling parameter
 * @param S the resulting scaling matrix
 */
void getScaleMatrix(double lambda, Mat& S) {

	// init S
	S = Mat::eye(3, 3, CV_32FC1);

	// The scaling is represented on main diagonal of S, except the last one (2, 2)
	S.at<float>(0, 0) = S.at<float>(1, 1) = (float)lambda;
}

/**
 * combines translation-, rotation-, and scaling-matrices to a single transformation matrix H
 * @param T	: translation matrix
 * @param R	: rotation matrix
 * @param S	: scaling matrix
 * @param H	: resulting homography
 */
void getH(Mat& T, Mat& R, Mat& S, Mat& H) {

	// init H
	H = Mat::eye(3, 3, CV_32FC1);

	// Combine the matrices
	// 1st - Rotation and Scaling
	for (int i = 0; i < 4; i++) {
		H.at<float>(i / 2, i % 2) = R.at<float>(i / 2, i % 2) * S.at<float>(i / 2, i % 2);
	}

	// 2nd - Translation
	H.at<float>(0, 2) = T.at<float>(0, 2);
	H.at<float>(1, 2) = T.at<float>(1, 2);
}

/**
 * transforms a geometric object by a given homography
 * @param geomObj the geometric object to be transformed
 * @param H the homography defining the transformation
 * @param result the transformed object
 * @param type the type of the geometric object (for now: only point and line)
 */
void applyH(Mat& geomObj, Mat& H, Mat& result, string type) {

	// Init the result matrix
	result = Mat::zeros(3, 1, CV_32FC1);

	// Store temporary the transformation matrix
	Mat Temp_H;

	if (type.compare("point") == 0) {
		// if object is a point
		Temp_H = H;
	} else if (type.compare("line") == 0) {
		// if object is a line
		cv::transpose(H.inv(), Temp_H);
	} else {
		cerr << "ERROR: Do not know how to transform " << type << endl;
		return;
	}

	result.at<float>(0, 0) = geomObj.at<float>(0, 0) * Temp_H.at<float>(0, 0) + geomObj.at<float>(1, 0) * Temp_H.at<float>(0, 1) + geomObj.at<float>(2, 0) * Temp_H.at<float>(0, 2);
	result.at<float>(1, 0) = geomObj.at<float>(0, 0) * Temp_H.at<float>(1, 0) + geomObj.at<float>(1, 0) * Temp_H.at<float>(1, 1) + geomObj.at<float>(2, 0) * Temp_H.at<float>(1, 2);
	result.at<float>(2, 0) = geomObj.at<float>(0, 0) * Temp_H.at<float>(2, 0) + geomObj.at<float>(1, 0) * Temp_H.at<float>(2, 1) + geomObj.at<float>(2, 0) * Temp_H.at<float>(2, 2);
}

/**
 * checks if a point is on a line
 * @param point the given point
 * @param line the given line
 * @param eps the used accuracy (set to 10^-5 by default (see header))
 */
bool isPointOnLine(Mat& point, Mat& line, double eps) {

	// Calculate the distance between line and point
	double Value = (double)(line.at<float>(0, 0) * point.at<float>(0, 0) + line.at<float>(1, 0) * point.at<float>(1, 0) + line.at<float>(2, 0) * point.at<float>(2, 0));

	// Output of calculation result for debug
	cout << "Check for \"Point in Line\":  Value (" << Value << ") | EPS (" << eps << ")" << endl;

	// if distance is to large return false, else true
	return (Value > eps) ? false : true;
}
