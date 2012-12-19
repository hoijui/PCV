//============================================================================
// Name        : pcv3.cpp
// Author      : Ronny Haensch
// Version     : 1.0
// Copyright   : -
// Description : Loads calibration image and calibrates camera
// Group       : Marcus Grum (), Robin Vobruba (), Jens Jawer () und Marcus Pannwitz (343479)
//============================================================================

#ifdef _WIN32
	#include "stdafx.h"
#endif
#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

// global variable to save point pairs
vector<Point2f> pointList;

// some function headers (see below)
// functions to be implemented
Mat calibrate(Mat& points2D, Mat& points3D);
Mat solve_dlt(Mat& A);
void decondition(Mat& T_2D, Mat& T_3D, Mat& P);
Mat getDesignMatrix_camera(Mat& points2D, Mat& points3D);
Mat transform(Mat& H, Mat& p);
Mat getCondition2D(Mat& p);
Mat getCondition3D(Mat& p);
void interprete(Mat& P);
// given functions
int getPoints(struct winInfo& calib, string fname, Mat& points2D, Mat& points3D);
void getPoints(int event, int x, int y, int flags, void* param);

// structure to fuse image data and window title
struct winInfo {Mat img; string name;};

// usage: path to calibration image in argv[1], path to object-point file in argv[2]
// main function
int main(int argc, char** argv) {

	// check if file paths were defined
	if (argc != 3) {
		cerr << "Usage: pcv3 <path to calibration image> <path to object-point file>" << endl;
		return -1;
	}

	// titles of point-selection windows
	string winName1 = string ("Calibration image");

	// load image calibration image, paths in argv[1]
	Mat calibImage = imread(argv[1]);
	if (!calibImage.data) {
		cerr << "ERROR: Could not load image " << argv[1] << endl;
		return -1;
	}
	// fuse image data and window title
	struct winInfo calib;
	calib.img = calibImage.clone();
	calib.name = winName1;

	// get corresponding points within the image
	Mat points2D, points3D;
	string fname = argv[2];
	int numberOfPointPairs = getPoints(calib, fname, points2D, points3D);

	// just some putput
	cout << "Number of defined point pairs: " << numberOfPointPairs << endl;
	cout << endl << "2D Points in image:" << endl;
	cout << points2D << endl;
	cout << endl << "3D Points at object:" << endl;
	cout << points3D << endl;

	// calculate projection matrix
	Mat P = calibrate(points2D, points3D);

	// decompose P to get camera parameter
	interprete(P);

	return 0;
}

// Extract and print information about interior and exterior orientation from camera (all 11 parameters)
/*
P	The 3x4 projection matrix
*/
void interprete(Mat& P) {

	// see pcv8_WS1213_cameramodel.pdf page 11+
	// see pcv9_WS1213_projectionmatrix.pdf page 9

	Mat M(P, Rect(0, 0, 3, 3));
	Mat Q_x = Mat::eye(3, 3, CV_32FC1);
	Mat Q_y = Mat::eye(3, 3, CV_32FC1);
	Mat Q_z = Mat::eye(3, 3, CV_32FC1);
	
	float t = sqrt( pow(M.at<float>(2, 1), 2) + pow(M.at<float>(2, 2), 2) );
	float s = -M.at<float>(2, 1) / t;
	float c =  M.at<float>(2, 2) / t;
	Q_x.at<float>(1, 1) = Q_x.at<float>(2, 2) = c;
	Q_x.at<float>(1, 2) = -s;
	Q_x.at<float>(2, 1) =  s;
	
	Mat M_1 = M * Q_x;
	M_1.at<float>(2, 1) = 0;
	
	t = sqrt( pow(M_1.at<float>(2, 0), 2) + pow(M_1.at<float>(2, 2), 2) );
	s = M.at<float>(2, 0) / t;
	c = M.at<float>(2, 2) / t;
	Q_y.at<float>(0, 0) = Q_x.at<float>(2, 2) = c;
	Q_y.at<float>(0, 2) =  s;
	Q_y.at<float>(2, 0) = -s;

	Mat M_2 = M_1 * Q_y;
	M_2.at<float>(2, 0) = M_2.at<float>(2, 1) = 0;
	
	t = sqrt( pow(M_1.at<float>(1, 0), 2) + pow(M_1.at<float>(1, 1), 2) );
	s = -M.at<float>(1, 0) / t;
	c =  M.at<float>(1, 1) / t;
	Q_z.at<float>(0, 0) = Q_x.at<float>(1, 1) = c;
	Q_z.at<float>(0, 1) = -s;
	Q_z.at<float>(1, 0) =  s;

	Mat M_3 = M_2 * Q_z;
	M_2.at<float>(1, 0) = M_2.at<float>(2, 0) = M_2.at<float>(2, 1) = 0;

	Mat R = M * Q_x * Q_y * Q_z;
	Mat Q = Q_z.t() * Q_y.t() * Q_x.t();

	cout << "Interpretation of the 3x4 projection matrix" << endl;
	cout << "Matrix R:" << endl;
	cout << R << endl;
	cout << endl << "Matrix Q:" << endl;
	cout << Q << endl;
}

static Mat createConditioningMatrix(const Mat& points) {

	const int d = points.rows;
	const int n = points.cols;

	// center
	Mat center = Mat::zeros(d, 1, points.type());
	for (int c = 0; c < n; ++c) {
		center += points.col(c);
	}
	center /= n;
	Mat centeringTranslation = Mat::zeros(d, d, points.type());
	centeringTranslation.col(d - 1) = -center;

	// set std. deviation to 1
	Scalar averageSqDist = 0.0;
	for (int c = 0; c < n; ++c) {
		const Scalar distSq = sum(points.col(c) - center);
		averageSqDist += distSq;
	}
	averageSqDist /= (Scalar)n;
	Scalar invAverageSqDist = (Scalar)1 / averageSqDist;
	Mat stdDeviationNormalization = Mat::eye(d, d, points.type());
	//stdDeviationNormalization = stdDeviationNormalization / averageSqDist;

	//return stdDeviationNormalization * centeringTranslation;
}

// estimate projection matrix
/*
points2D	set of 2D points within the image
points3D	set of 3D points at the object
return		the projection matrix to be computed
*/
Mat calibrate(Mat& points2D, Mat& points3D) {

	// see pcv9_WS1213_projectionmatrix.pdf page 7

	Mat ConMat2D = getCondition2D(points2D);
	Mat ConMat3D = getCondition3D(points3D);
	Mat TransfPoints2D = transform(ConMat2D, points2D);
	Mat TransfPoints3D = transform(ConMat3D, points3D);
	Mat DesMat = getDesignMatrix_camera(TransfPoints2D, TransfPoints3D);
	Mat EstProjMat = solve_dlt(DesMat);
	decondition(ConMat2D, ConMat3D, EstProjMat);

	return EstProjMat;
}

// solve homogeneous equation system by usage of SVD
/*
A		the design matrix
return		the estimated projection matrix
*/
Mat solve_dlt(Mat& A) {

	Mat X;
	SVD::solveZ(A, X);
	return X.reshape(1, 3);
}

// decondition a projection matrix that was estimated from conditioned point clouds
/*
T_2D	conditioning matrix of set of 2D image points
T_3D	conditioning matrix of set of 3D object points
P	conditioned projection matrix that has to be un-conditioned (in-place)
*/
void decondition(Mat& T_2D, Mat& T_3D, Mat& P) {

	// Decondition of homography matrix P
	P = T_2D.inv() * P * T_3D;
}

// define the design matrix as needed to compute projection matrix
/*
points2D	set of 2D points within the image
points3D	set of 3D points at the object
return		the design matrix to be computed
*/
Mat getDesignMatrix_camera(Mat& points2D, Mat& points3D) {

	Mat base = points2D;
	Mat attach = points3D;

	// The design matrix has at least 9 rows in case of 4 points. If there are more points, the number of rows equals two times the number of points
	Mat T = (base.cols < 5) ? Mat::zeros(9, 12, CV_32FC1) : Mat::zeros(base.cols * 2, 12, CV_32FC1);

	// Loop through all corresponding points
	for(int i = 0; i < base.cols; ++i)
	{
		int y = i * 2;

		for(int j = 0; j < 4; ++j)
		{
			// Two times: -w' * x
			T.at<float>(y + 1, j + 4) = T.at<float>(y, j) = -base.at<float>(2, i) * attach.at<float>(j, i);

			// u' * x
			T.at<float>(y, 8 + j) = base.at<float>(0, i) * attach.at<float>(j, i);

			// v' * x
			T.at<float>(y + 1, 8 + j) = base.at<float>(1, i) * attach.at<float>(j, i);
		}
	}

	return T;
}

// apply transformation to set of points
/*
H		matrix representing the transformation
p		input points
return		transformed points
*/
Mat transform(Mat& H, Mat& p) {

	Mat p_trans = Mat(p.size(), p.type());

	for (int i = 0; i < p.cols; i++)
	{
		Mat tPoint = H * p.col(i);

		for (int j = 0; j < tPoint.rows; j++)
		{
			p_trans.at<float>(j, i) = tPoint.at<float>(j, 0);
		}
	}

	return p_trans;
}

// get the conditioning matrix of given points
/*
p		the points as matrix
return		the condition matrix
*/
Mat getCondition2D(Mat& p) {

	// calculate center
	float transX = 0, transY = 0;

	// for each point
	for(int i = 0; i < p.cols; ++i)
	{
		transX += p.at<float>(0, i) / p.at<float>(2, i);
		transY += p.at<float>(1, i) / p.at<float>(2, i);
	}

	transX /= p.cols;
	transY /= p.cols;

	// calculate scale
	float scaleX = 0, scaleY = 0;

	// for each point
	for(int i = 0; i < p.cols; ++i)
	{
		scaleX += abs((p.at<float>(0, i) / p.at<float>(2, i)) - transX);
		scaleY += abs((p.at<float>(1, i) / p.at<float>(2, i)) - transY);
	}

	scaleX /= p.cols;
	scaleY /= p.cols;

	// build condition matrix
	Mat cond = Mat::eye(3, 3, CV_32FC1);
	cond.at<float>(0, 0) = 1 / scaleX;
	cond.at<float>(0, 2) = -transX / scaleX;
	cond.at<float>(1, 1) = 1 / scaleY;
	cond.at<float>(1, 2) = -transY / scaleY;

	return cond;
}

// get the conditioning matrix of given points
/*
	p		the points as matrix
	return	the condition matrix 
*/
Mat getCondition3D(Mat& p)
{
	// calculate center
	float transX = 0, transY = 0, transZ = 0;

	// for each point
	for(int i = 0; i < p.cols; ++i)
	{
		transX += p.at<float>(0, i) / p.at<float>(3, i);
		transY += p.at<float>(1, i) / p.at<float>(3, i);
		transZ += p.at<float>(2, i) / p.at<float>(3, i);
	}

	transX /= p.cols;
	transY /= p.cols;
	transZ /= p.cols;

	// calculate scale
	float scaleX = 0, scaleY = 0, scaleZ = 0;

	// for each point
	for(int i = 0; i < p.cols; ++i)
	{
		scaleX += abs((p.at<float>(0, i) / p.at<float>(3, i)) - transX);
		scaleY += abs((p.at<float>(1, i) / p.at<float>(3, i)) - transY);
		scaleZ += abs((p.at<float>(2, i) / p.at<float>(3, i)) - transZ);
	}

	scaleX /= p.cols;
	scaleY /= p.cols;
	scaleZ /= p.cols;

	// build condition matrix
	Mat cond = Mat::eye(4, 4, CV_32FC1);
	cond.at<float>(0, 0) = 1 / scaleX;
	cond.at<float>(0, 3) = -transX / scaleX;
	cond.at<float>(1, 1) = 1 / scaleY;
	cond.at<float>(1, 3) = -transY / scaleY;
	cond.at<float>(2, 2) = 1 / scaleY;
	cond.at<float>(2, 3) = -transZ / scaleZ;

	return cond;
}

/* *****************************
GIVEN FUNCTIONS
***************************** */

// display image and catch the point marked by left mouse clicks
// points have to be clicked in same order as in the file
// points will be in homogeneous coordinates
/*
calib		structure containing calibration image
fname		path to file that contains 3D real world (image points have to be defined in same order)
points2D	points within the image (to be defined by this method)
points3D	points at the object (to be read from file)
*/
int getPoints(struct winInfo& calib, string fname, Mat& points2D, Mat& points3D) {

	// show input images and install mouse callback
	namedWindow(calib.name.c_str(), 0);
	imshow(calib.name.c_str(), calib.img);
	setMouseCallback(calib.name.c_str(), getPoints, (void*)(&calib));
	// wait until any key was pressed
	waitKey(0);

	/*
	// for debugging purposes (points from slides)
	pointList.clear();
	pointList.push_back(Point2f(18.5,   46.8));
	pointList.push_back(Point2f(99.1,  146.5));
	pointList.push_back(Point2f(13.8,  221.8));
	pointList.push_back(Point2f(242.1,  52.5));
	pointList.push_back(Point2f(151.1, 147.1));
	pointList.push_back(Point2f(243.1, 224.5));
	*/

	// allocate memory for point-lists (represented as matrix)
	points2D = Mat(3, pointList.size(), CV_32FC1);
	points3D = Mat(4, pointList.size(), CV_32FC1);
	int n = 0; // number of point pairs

	fstream calibFile(fname.c_str(), ios::in);

	string buffer;
	int e;
	// read points from global variable, transform them into homogeneous coordinates
	for (vector<Point2f>::iterator p = pointList.begin(); p != pointList.end(); p++, n++) {
		// define image point
		points2D.at<float>(0, n) = p->x;
		points2D.at<float>(1, n) = p->y;
		points2D.at<float>(2, n) = 1;

		// read object points from file 
		getline(calibFile, buffer);
		e = buffer.find(" ");
		points3D.at<float>(0, n) = atof(buffer.substr(0,e).c_str());
		buffer = buffer.substr(e+1, string::npos);
		e = buffer.find(" ");
		points3D.at<float>(1, n) = atof(buffer.substr(0,e).c_str());
		buffer = buffer.substr(e+1, string::npos);
		points3D.at<float>(2, n) = atof(buffer.c_str());
		points3D.at<float>(3, n) = 1;
	}

	// be tidy
	destroyWindow(calib.name.c_str());
	calibFile.close();

	return n;
}

// mouse call back to get points and draw circles
/*
event	specifies encountered mouse event
x,y	position of mouse pointer
flags	not used here
param	a struct containing used IplImage and window title
*/
void getPoints(int event, int x, int y, int flags, void* param) {

	// cast to structure
	struct winInfo* win = (struct winInfo*)param;

	switch(event) {
		// if left mouse button was pressed
		case CV_EVENT_LBUTTONDOWN:{
			// create point representing mouse position
			Point2f p = Point2f(x,y);
			// draw green point
			circle(win->img, p, 2, Scalar(0, 255, 0), 2);
			// draw green circle
			circle(win->img, p, 15, Scalar(0, 255, 0), 2);
			// update image
			imshow(win->name.c_str(), win->img);
			// add point to point list
			pointList.push_back(p);
		} break;
	}
}
