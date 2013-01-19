//============================================================================
// Name        : pcv6.cpp
// Author      : Ronny Haensch
// Version     : 1.0
// Copyright   : -
// Description : projective and euclidian reconstruction
//============================================================================

#include <iostream>
#include <fstream>

#include <list>

#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

// global variable to save point pairs
vector<Point> pointList;

// some function headers (see below)
Mat getCondition2D(Mat& p);
Mat getCondition3D(Mat& p);
Mat transform(Mat& H, Mat& p);
Mat solve_dlt(Mat& A);
Mat getDesignMatrix_fundamental(Mat& n1, Mat& n2);
Mat getDesignMatrix_homography3D(Mat& fst, Mat& snd);
Mat getFundamentalMatrix(Mat& p1, Mat& p2);
void forceSingularity(Mat& F);
void decondition_fundamental(Mat& T1, Mat& T2, Mat& F);
void decondition_homography3D(Mat& T1, Mat& T2, Mat& H);
Mat homography3D(Mat& Xp, Mat& Xm);
void defineCameras(Mat& F, Mat& P1, Mat& P2);
Mat linearTriangulation(Mat& P1, Mat& P2, Mat& x1, Mat& x2);
void getEpipols(Mat& F, Mat& e1, Mat& e2);
Mat makeSkewMatrix(Mat& e);
// given functions
int readMatchingPoints(string fname, Mat& x1, Mat& x2);
int readControlPoints(string fname, Mat& x1, Mat& x2, Mat& Xm);
void savePointList(string fname, Mat& points);

// structure to fuse image data and window title
struct winInfo {Mat img; string name;};

// usage: path to image points in argv[1], path to control points in argv[2]
// main function
int main(int argc, char** argv) {

	// check if image paths were defined
	if (argc != 3) {
		cerr << "Usage: pcv6 <path to image points> <path to control points>" << endl;
		return -1;
	}

	// get corresponding points within the two images
	Mat x1, x2;
	int numberOfPointPairs = readMatchingPoints(argv[1], x1, x2);

	// just some putput
	cout << "Number of defined point pairs: " << numberOfPointPairs << endl;

	// calculate fundamental matrix
	Mat F = getFundamentalMatrix(x1, x2);

	cout << endl << "Fundamental matrix: " << F << endl;

	// define projection matrices
	Mat P1, P2;
	defineCameras(F, P1, P2);

	cout << endl << "Camera 1: " << P1 << endl;
	cout << endl << "Camera 2: " << P2 << endl;

	// linear triangulation of image points
	// resulting projective reconstruction
	Mat X = linearTriangulation(P1, P2, x1, x2);

	// save reconstructed points to file
	savePointList("projectiveReconstruction.asc", X);

	// read control points
	// clean re-use of x1- and x2-matrices
	Mat Xm;
	numberOfPointPairs = readControlPoints(argv[2], x1, x2, Xm);

	cout << endl << "Control points:" << endl;
	cout << "x1: " << x1 << endl;
	cout << "x2: " << x2 << endl;
	cout << "Xm: " << Xm << endl;

	// linear triangulation of image points
	Mat Xp = linearTriangulation(P1, P2, x1, x2);
	cout << "Reconstructed control points: " << Xp << endl;

	// Transform projective reconstructed points to euclidian reconstruction
	Mat H = homography3D(Xp, Xm);
	Mat X_final = transform(H, X);

	// save reconstructed points to file
	savePointList("euclidianReconstruction.asc", X_final);

	return 0;
}

// calculates epipols from fundamental matrix
/*
F	the fundamental matrix
e1	first epipol
e2	second epipol
*/
void getEpipols(Mat& F, Mat& e1, Mat& e2) {

	// TODO
}

// generates skew matrix from vector
/*
e		given vector
return		skew matrix
*/
Mat makeSkewMatrix(Mat& e) {

	// TODO
}

// generates 2 projection matrices by using fundamental matrix
/*
F	the fundamental matrix
P1	first projection matrix (standard P matrix)
P2	second projection matrix based on F
*/
void defineCameras(Mat& F, Mat& P1, Mat& P2) {

	// TODO
}

// triangulates given set of image points based on projection matrices
/*
P1	projection matrix of first image
P2	projection matrix of second image
x1	image point set of first image
x2	image point set of second image
return	triangulated object points
*/
Mat linearTriangulation(Mat& P1, Mat& P2, Mat& x1, Mat& x2) {

	// TODO
}

// computes 3D homography
/*
X1	first set of points
X2	second set of points
H	computed homography
*/
Mat homography3D(Mat& X1, Mat& X2) {

	// TODO
}

// decondition a homography that was estimated from conditioned point clouds
/*
T_to	conditioning matrix of first set of points
T_from	conditioning matrix of second set of points
H	conditioned homography that has to be un-conditioned (in-place)
*/
void decondition_homography3D(Mat& T_to, Mat& T_from, Mat& H) {

	// TODO
}

// compute fundamental matrix
/*
fst	first set of points
snd	second set of points
return	the estimated fundamental matrix
*/
Mat getFundamentalMatrix(Mat& fst, Mat& snd) {

	// TODO
}

// solve homogeneous equation system by usage of SVD
/*
A		the design matrix
return		the estimated fundamental matrix
*/
Mat solve_dlt(Mat& A) {

	// TODO
}

// decondition a fundamental matrix that was estimated from conditioned point clouds
/*
T_fst	conditioning matrix of first set of points
T_snd	conditioning matrix of second set of points
F	conditioned fundamental matrix that has to be un-conditioned (in-place)
*/
void decondition_fundamental(Mat& T_fst, Mat& T_snd, Mat& F) {

	// TODO
}

// define the design matrix as needed to compute fundamental matrix
/*
fst		first set of points
snd		second set of points
return		the design matrix to be computed
*/
Mat getDesignMatrix_fundamental(Mat& fst, Mat& snd) {

	// TODO
}

// define the design matrix as needed to compute fundamental matrix
/*
fst	first set of points
snd	second set of points
A	the design matrix to be computed
*/
Mat getDesignMatrix_homography3D(Mat& fst, Mat& snd) {

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
return		the condition matrix (already allocated)
*/
Mat getCondition2D(Mat& p) {

	// TODO
}

// get the conditioning matrix of given 3D points
/*
p	the points as matrix
T	the condition matrix (already allocated)
*/
Mat getCondition3D(Mat& p) {

	// TODO
}


// enforce rank of 2 on fundamental matrix
/*
F	the matrix to be changed
*/
void forceSingularity(Mat& F) {

	// TODO
}

/* ***********************
   *** Given Functions ***
   *********************** */

// saves point list to file
/*
fname		file name
points		matrix of points
*/
void savePointList(string fname, Mat& points) {

	// open file for write
	fstream file(fname.c_str(), ios::out);
	if (!file) {
		cerr << "ERROR: cannot open file " << fname << endl;
		return;
	}

	// if homogeneous points: norm and write points
	if (points.rows == 4) {
		for (int i=0; i < points.cols; i++) {
			file
					<< points.at<float>(0, i) / points.at<float>(3, i) << ","
					<< points.at<float>(1, i) / points.at<float>(3, i) << ","
					<< points.at<float>(2, i) / points.at<float>(3, i) << endl;
		}
	}
	// if euclidian points: write points
	if (points.rows == 3) {
		for (int i=0; i < points.cols; i++) {
			file
					<< points.at<float>(0, i) << ","
					<< points.at<float>(1, i) << ","
					<< points.at<float>(2, i) << endl;
		}
	}

	// close file
	file.close();
}

// read homologeous points from file
/*
fname	path of file containing point list
x1	pointer to matrix containing points of first image
x2	pointer to matrix containing points of second image
*/
int readMatchingPoints(string fname, Mat& x1, Mat& x2) {

	// open file
	fstream file(fname.c_str(), ios::in);
	if (!file) {
		cerr << "ERROR: Cannot open file " << fname << endl;
		exit(-1);
	}

	// read points into list
	list<Scalar> points;
	int end, numberOfPointPairs = 0;
	string buffer;
	// read line by line
	while (getline(file, buffer)) {
		// counting
		numberOfPointPairs++;

		// get semicolon
		end = buffer.find(';');
		// first point before semicolon
		string fst = buffer.substr(0, end);
		// second point after semicolon
		string snd = buffer.substr(end+1);

		// get x and y
		Scalar cur;
		end = fst.find(',');
		cur.val[0] = atof(fst.substr(0, end).c_str());
		cur.val[1] = atof(fst.substr(end+1).c_str());

		// get x and y
		end = snd.find(',');
		cur.val[2] = atof(snd.substr(0, end).c_str());
		cur.val[3] = atof(snd.substr(end+1).c_str());

		// push point pair to list
		points.push_back(cur);
	}

	// allocate memory for point matrices
	x1 = Mat(3, numberOfPointPairs, CV_32FC1);
	x2 = Mat(3, numberOfPointPairs, CV_32FC1);

	// fill point matrices
	int i=0;
	for (list<Scalar>::iterator p = points.begin(); p != points.end(); p++,i++) {
		// x1
		x1.at<float>(0, i) = (*p).val[0];
		x1.at<float>(1, i) = (*p).val[1];
		x1.at<float>(2, i) = 1;
		// x2
		x2.at<float>(0, i) = (*p).val[2];
		x2.at<float>(1, i) = (*p).val[3];
		x2.at<float>(2, i) = 1;
	}

	return numberOfPointPairs;
}

// read control points from file
/*
fname	path of file containing point list
x1	pointer to matrix containing points of first image
x2	pointer to matrix containing points of second image
Xm	pointer to matrix containing object points
*/
int readControlPoints(string fname, Mat& x1, Mat& x2, Mat& Xm) {

	// open file
	fstream file(fname.c_str(), ios::in);
	if (!file) {
		cerr << "ERROR: Cannot open file " << fname << endl;
		exit(-1);
	}

	// read points into list
	list< vector<double> > points;
	int pos1, pos2, numberOfPointPairs = 0;
	string buffer;
	// read line by line
	while (getline(file, buffer)) {
		// counting
		numberOfPointPairs++;

		// get semicolons
		pos1 = buffer.find(';');
		pos2 = buffer.rfind(';');
		// first point before first semicolon
		string fst = buffer.substr(0, pos1);
		// second point between semicolons
		string snd = buffer.substr(pos1+1, pos2);
		// third point after last semicolon
		string thd = buffer.substr(pos2+1);

		// current point triplet
		vector<double> cur;

		// get x and y of first point
		pos1 = fst.find(',');
		cur.push_back( atof(fst.substr(0, pos1).c_str()) );
		cur.push_back( atof(fst.substr(pos1+1).c_str()) );
		// get x and y of second point
		pos1 = snd.find(',');
		cur.push_back( atof(snd.substr(0, pos1).c_str()) );
		cur.push_back( atof(snd.substr(pos1+1).c_str()) );
		// get x, y, and z of object point
		pos1 = thd.find(',');
		pos2 = thd.rfind(',');
		cur.push_back( atof(thd.substr(0, pos1).c_str()) );
		cur.push_back( atof(thd.substr(pos1+1, pos2).c_str()) );
		cur.push_back( atof(thd.substr(pos2+1).c_str()) );

		// push point triplet to list
		points.push_back(cur);
	}

	// allocate memory for point matrices
	x1 = Mat(3, numberOfPointPairs, CV_32FC1);
	x2 = Mat(3, numberOfPointPairs, CV_32FC1);
	Xm = Mat(4, numberOfPointPairs, CV_32FC1);

	// fill point matrices
	int i=0;
	for (list< vector<double> >::iterator p = points.begin(); p != points.end(); p++,i++) {
		// x1
		x1.at<float>(0, i) = (*p).at(0);
		x1.at<float>(1, i) = (*p).at(1);
		x1.at<float>(2, i) = 1;
		// x2
		x2.at<float>(0, i) = (*p).at(2);
		x2.at<float>(1, i) = (*p).at(3);
		x2.at<float>(2, i) = 1;
		// Xm
		Xm.at<float>(0, i) = (*p).at(4);
		Xm.at<float>(1, i) = (*p).at(5);
		Xm.at<float>(2, i) = (*p).at(6);
		Xm.at<float>(3, i) = 1;
	}

	return numberOfPointPairs;
}
