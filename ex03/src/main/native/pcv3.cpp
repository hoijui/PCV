//============================================================================
// Name        : pcv3.cpp
// Author      : Ronny Haensch
// Version     : 1.0
// Copyright   : -
// Description : loads calibration image and calibrates camera
//============================================================================

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

// extract and print information about interior and exterior orientation from camera (all 11 parameters)
/*
P	The 3x4 projection matrix
*/
void interprete(Mat& P) {

	// see pcv8_WS1213_cameramodel.pdf page 11+
	// see pcv9_WS1213_projectionmatrix.pdf page 9

	Mat M = P.colRange(0, 2);
	Mat upper(3, 3, P.type());
	Mat rest(3, 3, P.type());
	RQDecomp3x3(M, upper, rest);
	Mat rotationMatrix = rest;
	Mat calibrationMatrix = upper;

	Mat P234 = P.colRange(1, 3);
	Mat P134(3, 3, P.type());
	P134.row(0) = P.row(0);
	P134.row(1) = P.row(2);
	P134.row(2) = P.row(3);
	Mat P124(3, 3, P.type());
	P134.row(0) = P.row(0);
	P134.row(1) = P.row(1);
	P134.row(2) = P.row(3);
	Mat P123 = P.colRange(0, 2);

	Mat projectionCenter(3, 1, CV_32FC1);
	projectionCenter.at<float>(0, 0) = determinant(P234);
	projectionCenter.at<float>(1, 0) = determinant(P134);
	projectionCenter.at<float>(2, 0) = determinant(P124);
	const float W = determinant(P123);
	projectionCenter = projectionCenter / W;

	Mat principlePoint(2, 1, CV_32FC1);
	principlePoint.at<float>(0, 0) = calibrationMatrix.at<float>(0, 3);
	principlePoint.at<float>(1, 0) = calibrationMatrix.at<float>(1, 3);

	const float principleDistance = calibrationMatrix.at<float>(0, 0);

	const float skew = calibrationMatrix.at<float>(0, 1);

	float omega = 0.0f; // TODO extract from rotation matrix
	float phi = 0.0f; // TODO extract from rotation matrix
	float kappa = 0.0f; // TODO extract from rotation matrix

	// Interior orientation
	cout << "C (projection center): " << projectionCenter << endl;
	cout << "R (rotation matrix): " << rotationMatrix << endl;
	cout << "Rotation angles (omega, phi, kappa): " << omega << ", " << phi << ", " << kappa << endl;

	// Exterior orientation
	cout << "alpha_x (principle distance) [px]: " << principleDistance << endl;
	cout << "s (skew): " << skew << endl;
	cout << "principle point (x_0, y_0) [px]: " << principlePoint << endl;
	const float aspectRatio = 0.0f;
	cout << "aspect ratio (gamme = alpha_y / alpha_x) [<no unit>]: " << aspectRatio << endl;

	// TODO
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
	stdDeviationNormalization = stdDeviationNormalization / averageSqDist;

	return stdDeviationNormalization * centeringTranslation;
}

// estimate projection matrix
/*
points2D	set of 2D points within the image
points3D	set of 3D points at the object
return		the projection matrix to be computed
*/
Mat calibrate(Mat& points2D, Mat& points3D) {

	// see pcv9_WS1213_projectionmatrix.pdf page 7

	const int n = points2D.cols;

	const Mat conditioner2D = createConditioningMatrix(points2D); // T'
	const Mat conditioner3D = createConditioningMatrix(points3D); // T

	// condition
	const Mat points2DConditioned = conditioner2D * points2D;
	const Mat points3DConditioned = conditioner3D * points3D;

	Mat designMatrix = Mat::zeros(2*n, 12, points2D.type());
	//vector<Mat> designMatrices;
	for (int c = 0; c < n; ++c) {
		const float u = points2DConditioned.at<float>(0, c);
		const float v = points2DConditioned.at<float>(1, c);
		const float w = points2DConditioned.at<float>(2, c);
		const Mat& X = points3DConditioned.col(c);
		//Mat designMatrix = Mat::zeros(2, 12, points2D.type());
		// -w * X^T
		designMatrix.at<float>(2*c + 0, 0)  = -w * X.at<float>(0);
		designMatrix.at<float>(2*c + 0, 1)  = -w * X.at<float>(1);
		designMatrix.at<float>(2*c + 0, 2)  = -w * X.at<float>(2);
		designMatrix.at<float>(2*c + 0, 3)  = -w * X.at<float>(3);
		// -w * X^T
		designMatrix.at<float>(2*c + 1, 4)  = -w * X.at<float>(0);
		designMatrix.at<float>(2*c + 1, 5)  = -w * X.at<float>(1);
		designMatrix.at<float>(2*c + 1, 6)  = -w * X.at<float>(2);
		designMatrix.at<float>(2*c + 1, 7)  = -w * X.at<float>(3);
		// u * X^T
		designMatrix.at<float>(2*c + 0, 8)  = u * X.at<float>(0);
		designMatrix.at<float>(2*c + 0, 9)  = u * X.at<float>(1);
		designMatrix.at<float>(2*c + 0, 10) = u * X.at<float>(2);
		designMatrix.at<float>(2*c + 0, 11) = u * X.at<float>(3);
		// v * X^T
		designMatrix.at<float>(2*c + 1, 8)  = v * X.at<float>(0);
		designMatrix.at<float>(2*c + 1, 9)  = v * X.at<float>(1);
		designMatrix.at<float>(2*c + 1, 10) = v * X.at<float>(2);
		designMatrix.at<float>(2*c + 1, 11) = v * X.at<float>(3);

		//designMatrices.push_back(designMatrix);
	}

	// DEPRECATED TODO make one big design matrix with columns merged?

	Mat p;
	SVD::solveZ(designMatrix, p); // designMatrix * p = 0

	Mat Ptilde = p.reshape(3, 4);

	Mat P = conditioner2D * Ptilde * conditioner3D;

	return P;
}

// solve homogeneous equation system by usage of SVD
/*
A		the design matrix
return		the estimated projection matrix
*/
Mat solve_dlt(Mat& A) {

	// TODO
}

// decondition a projection matrix that was estimated from conditioned point clouds
/*
T_2D	conditioning matrix of set of 2D image points
T_3D	conditioning matrix of set of 3D object points
P	conditioned projection matrix that has to be un-conditioned (in-place)
*/
void decondition(Mat& T_2D, Mat& T_3D, Mat& P) {

	// TODO
}

// define the design matrix as needed to compute projection matrix
/*
points2D	set of 2D points within the image
points3D	set of 3D points at the object
return		the design matrix to be computed
*/
Mat getDesignMatrix_camera(Mat& points2D, Mat& points3D) {

	// TODO
}

// apply transformation to set of points
/*
H		matrix representing the transformation
p		input points
return		transformed points
*/
Mat transform(Mat& H, Mat& p) {

	// TODO
}

// get the conditioning matrix of given points
/*
p		the points as matrix
return		the condition matrix
*/
Mat getCondition2D(Mat& p) {

	// TODO
}

// get the conditioning matrix of given points
/*
p		the points as matrix
return		the condition matrix 
*/
Mat getCondition3D(Mat& p) {

	// TODO
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
