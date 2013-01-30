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
void setSmallValuesToZero(Mat& mat);
void savePointListOFF(string fname, Mat& points);
float sumMembers(const Mat& mat);
float absDifferenceSum(const Mat& mat1, const Mat& mat2);

// structure to fuse image data and window title
struct winInfo {Mat img; string name;};

/**
 * usage: path to image points in argv[1], path to control points in argv[2]
 */
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
	for (int di = 0; di < X.cols; ++di) {
		cout << endl << "X(" << di << ")=" << X.col(di) << endl;
	}

	// save reconstructed points to file
	savePointList("projectiveReconstruction.asc", X);
	savePointListOFF("projectiveReconstruction.off", X);

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
	cout << endl << endl;
	float reconstructionErrorSum = absDifferenceSum(Xp, Xm);
	cout << "Reconstruction error: " << reconstructionErrorSum << endl;

	// Transform projective reconstructed points to euclidian reconstruction
	Mat H = homography3D(Xp, Xm);
	Mat X_final = transform(H, X);

	// save reconstructed points to file
	savePointList("euclidianReconstruction.asc", X_final);
	savePointListOFF("euclidianReconstruction.off", X);

	return 0;
}

/**
 * calculates epipols from fundamental matrix
 * @param F  the fundamental matrix
 * @param e1  first epipol
 * @param e2  second epipol
 */
void getEpipols(Mat& F, Mat& e1, Mat& e2) {

	// # dimensions
	const int D = F.rows; // 3

	Mat U(D, D, CV_32FC1);
	Mat diag(D, 1, CV_32FC1);
	Mat Vt(D, D, CV_32FC1);

	// decompose
	SVD::compute(F, diag, U, Vt);

	// extract epipoles
	// see pcv12_WS1213_relor.pdf page 8
	e1 = Vt.col(D-1); // last column
	e2 = U.col(D-1); // last column
}

/**
 * generates skew matrix from vector
 * @param e  given vector
 * @return  skew matrix
 */
Mat makeSkewMatrix(Mat& e) {

	// # dimensions
	const int D = e.rows; // 3

	Mat skew = Mat::zeros(D, D, e.type());
	skew.at<float>(0, 1) = -e.at<float>(2, 0);
	skew.at<float>(1, 0) =  e.at<float>(2, 0);
	skew.at<float>(0, 2) =  e.at<float>(1, 0);
	skew.at<float>(2, 0) = -e.at<float>(1, 0);
	skew.at<float>(1, 2) = -e.at<float>(0, 0);
	skew.at<float>(2, 1) =  e.at<float>(0, 0);

	return skew;
}

/**
 * generates 2 projection matrices by using fundamental matrix
 * @param F  the fundamental matrix
 * @param P1  first projection matrix (standard P matrix)
 * @param P2  second projection matrix based on F
 */
void defineCameras(Mat& F, Mat& P1, Mat& P2) {

	// # dimensions
	const int D = F.rows; // 3

	P1 = Mat::eye(D, D + 1, F.type());

	P2 = Mat::zeros(D, D + 1, F.type()); // 3 x 4 matrix
	Mat e1;
	Mat e2;
	getEpipols(F, e1, e2);
	P2.colRange(0, 3) = makeSkewMatrix(e2) * F;
	// P2.col(3) = e2;
	for (int i = 0; i < e2.rows; ++i) {
		P2.at<float>(i, 3) = e2.at<float>(i, 0);
	}
}

/**
 * triangulates given set of image points based on projection matrices
 * @param P1  projection matrix of first image
 * @param P2  projection matrix of second image
 * @param x1  image point set of first image
 * @param x2  image point set of second image
 * @return  triangulated object points
 */
Mat linearTriangulation(Mat& P1, Mat& P2, Mat& x1, Mat& x2) {

	// number of point sets
	const int n = x1.cols;

	const int D = P1.cols;

	Mat X_Os(D, n, P1.type());

	// see pcv13_WS1213_triangulation.pdf page 10
	Mat A(D, D, P1.type());
	for (int pi = 0; pi < n; ++pi) {
		const float x  = x1.at<float>(0, pi);
		const float y  = x1.at<float>(1, pi);
		const float x_ = x2.at<float>(0, pi);
		const float y_ = x2.at<float>(1, pi);

		for (int di = 0; di < D; ++di) {
			// x *  (p3 - p1)
			A.at<float>(0, di) = (x  * P1.at<float>(2, di)) - P1.at<float>(0, di);

			// y *  (p3 - p2)
			A.at<float>(1, di) = (y  * P1.at<float>(2, di)) - P1.at<float>(1, di);

			// x_ * (p_3 - p_1)
			A.at<float>(2, di) = (x_ * P2.at<float>(2, di)) - P2.at<float>(0, di);

			// y_ * (p_3 - p_2)
			A.at<float>(3, di) = (y_ * P2.at<float>(2, di)) - P2.at<float>(1, di);
		}

		Mat X_O = Mat::zeros(D, 1, P1.type());
		SVD::solveZ(A, X_O);
		// X_Os.col(pi) = X_O;
		for (int i = 0; i < X_O.rows; ++i) {
			X_Os.at<float>(i, pi) = X_O.at<float>(i, 0);
		}
	}

	return X_Os;
}

/**
 * computes 3D homography
 * @param X1  first set of points
 * @param X2  second set of points
 * @param H  computed homography
 */
Mat homography3D(Mat& X1, Mat& X2) {

	Mat& base = X1;
	Mat& attach = X2;

	// coordinate transformation matrices for conditioning
	Mat tBase = getCondition3D(base);
	Mat tAttach = getCondition3D(attach);

	// create condition matrices
	Mat cBase = transform(tBase, base);
	Mat cAttach = transform(tAttach, attach);

	// create design matrix, for projection attach -> base
	Mat design = getDesignMatrix_homography3D(cBase, cAttach);

	// SVD and reshaping
	Mat H = solve_dlt(design);

	// creating final homography
	decondition_homography3D(tBase, tAttach, H);

	return H;
}

/**
 * decondition a homography that was estimated from conditioned point clouds
 * @param T_to	conditioning matrix of first set of points
 * @param T_from	conditioning matrix of second set of points
 * @param H	conditioned homography that has to be un-conditioned (in-place)
 */
void decondition_homography3D(Mat& T_to, Mat& T_from, Mat& H) {

	H = T_from.t() * H * T_to;
}

/**
 * compute fundamental matrix
 * @param fst  first set of points
 * @param snd  second set of points
 * @return  the estimated fundamental matrix
 */
Mat getFundamentalMatrix(Mat& fst, Mat& snd) {

	// coordinate transformation matrices for conditioning
	Mat tFirst = getCondition2D(fst);
	Mat tSecond = getCondition2D(snd);

	// create conditioned matrices
	Mat cFirst = transform(tFirst, fst);
	Mat cSecond = transform(tSecond, snd);

	// create design matrix, for fundamental matrix between fst and snd
	Mat design = getDesignMatrix_fundamental(cFirst, cSecond);

	// svd and reshaping
	Mat F = solve_dlt(design);

	// enforce singularity of F: SVD, kill one rank, recompose F
	forceSingularity(F);

	// creating final fundamental matrix
	decondition_fundamental(tFirst, tSecond, F);

	return F;
}

/**
 * solve homogeneous equation system by usage of SVD
 * @param A  the design matrix
 * @return  the estimated fundamental matrix
 */
Mat solve_dlt(Mat& A) {

	const int n = sqrt((float)A.cols);
	Mat f = Mat::zeros(1, n*n, CV_32FC1);
	SVD::solveZ(A, f);
	return f.reshape(0, n);
}

/**
 * decondition a fundamental matrix that was estimated from conditioned point clouds
 * @param T_fst  conditioning matrix of first set of points
 * @param T_snd  conditioning matrix of second set of points
 * @param F  conditioned fundamental matrix that has to be un-conditioned (in-place)
 */
void decondition_fundamental(Mat& T_fst, Mat& T_snd, Mat& F) {

	F = T_snd.t() * F * T_fst;
}

/**
 * define the design matrix as needed to compute fundamental matrix
 * @param fst  first set of points
 * @param snd  second set of points
 * @return  the design matrix to be computed
 */
Mat getDesignMatrix_fundamental(Mat& fst, Mat& snd) {

	// The design matrix has at least 9 rows in case of 8 points. If there are more points, the number of rows equals the number of points
	Mat design = (fst.cols < 8) ? Mat::zeros(9, 9, CV_32FC1) : Mat::zeros(fst.cols, 9, CV_32FC1);

	// Loop through all corresponding points
	for (int i = 0; i < fst.cols; ++i) {
		float x_fst = fst.at<float>(0, i),
		x_snd = snd.at<float>(0, i),
		y_fst = fst.at<float>(1, i),
		y_snd = snd.at<float>(1, i);

		design.at<float>(i, 0) = x_fst * x_snd;
		design.at<float>(i, 1) = y_fst * x_snd;
		design.at<float>(i, 2) = x_snd;

		design.at<float>(i, 3) = x_fst * y_snd;
		design.at<float>(i, 4) = y_fst * y_snd;
		design.at<float>(i, 5) = y_snd;

		design.at<float>(i, 6) = x_fst;
		design.at<float>(i, 7) = y_fst;
		design.at<float>(i, 8) = 1;
	}

	// "Speed up for large rectangular matrices (i.e. m > 2n) by using A^T*A instead of A"; see pcv5_WS1213_DLT.pdf, page 25
	if (design.rows > 2 * design.cols) {
		design = design.t() * design;
	}

	return design;
}

/**
 * define the design matrix as needed to compute fundamental matrix
 * @param fst  first set of points
 * @param snd  second set of points
 * @param A  the design matrix to be computed
 */
Mat getDesignMatrix_homography3D(Mat& fst, Mat& snd) {

	const Mat& base = fst;
	const Mat& attach = snd;

	// design matrix: at least 5 points required ->
	// size at least (16 x 16), if more points selected, then (3*#points x 16)
	// see pcv7_WS1213_3Dhomo.pdf, page 14
	Mat designMat = Mat::zeros(base.cols <= 5 ? 16 : base.cols * 3, 16, CV_32FC1);

	for (int i = 0; i < base.cols; ++i) {
		const int r1 = i * 2;
		const int r2 = r1 + 1;
		const int r3 = r1 + 2;

		// two times: -t' * x
		designMat.at<float>(r3, 0 + 8) = designMat.at<float>(r2, 0 + 4) = designMat.at<float>(r1, 0) = -base.at<float>(3, i) * attach.at<float>(0, i);  // -t' * x
		designMat.at<float>(r3, 1 + 8) = designMat.at<float>(r2, 1 + 4) = designMat.at<float>(r1, 1) = -base.at<float>(3, i) * attach.at<float>(1, i);  // -t' * y
		designMat.at<float>(r3, 2 + 8) = designMat.at<float>(r2, 2 + 4) = designMat.at<float>(r1, 2) = -base.at<float>(3, i) * attach.at<float>(2, i);  // -t' * z
		designMat.at<float>(r3, 3 + 8) = designMat.at<float>(r2, 3 + 4) = designMat.at<float>(r1, 3) = -base.at<float>(3, i) * attach.at<float>(3, i);  // -t' * w

		// u' * x
		designMat.at<float>(r1, 12) = base.at<float>(0, i) * attach.at<float>(0, i);
		designMat.at<float>(r1, 13) = base.at<float>(0, i) * attach.at<float>(1, i);
		designMat.at<float>(r1, 14) = base.at<float>(0, i) * attach.at<float>(2, i);
		designMat.at<float>(r1, 15) = base.at<float>(0, i) * attach.at<float>(3, i);

		// v' * x
		designMat.at<float>(r2, 12) = base.at<float>(1, i) * attach.at<float>(0, i);
		designMat.at<float>(r2, 13) = base.at<float>(1, i) * attach.at<float>(1, i);
		designMat.at<float>(r2, 14) = base.at<float>(1, i) * attach.at<float>(2, i);
		designMat.at<float>(r2, 15) = base.at<float>(1, i) * attach.at<float>(3, i);

		// w' * x
		designMat.at<float>(r3, 12) = base.at<float>(2, i) * attach.at<float>(0, i);
		designMat.at<float>(r3, 13) = base.at<float>(2, i) * attach.at<float>(1, i);
		designMat.at<float>(r3, 14) = base.at<float>(2, i) * attach.at<float>(2, i);
		designMat.at<float>(r3, 15) = base.at<float>(2, i) * attach.at<float>(3, i);
	}

	return designMat;
}

/**
 * apply transformation to set of points
 * @param H  matrix representing the transformation
 * @param p  input points
 * @return  transformed points
 */
Mat transform(Mat& H, Mat& p) {

	return H * p;
}

static Mat getConditionXD(Mat& p) {

	// # dimensions
	const int D = p.rows;

	// calculate center
	Mat center = Mat::zeros(D, 1, p.type());
	for (int pi = 0; pi < p.cols; ++pi) {
		const float w = p.at<float>(D - 1, pi);
		for (int di = 0; di < (D - 1); ++di) {
			center.at<float>(di, 0) += p.at<float>(di, pi) / w;
		}
	}
	center /= p.cols;

	// calculate scale
	Mat scale = Mat::zeros(D, 1, p.type());
	for (int pi = 0; pi < p.cols; ++pi) {
		const float w = p.at<float>(D - 1, pi);
		for (int di = 0; di < (D - 1); ++di) {
			scale.at<float>(di, 0) += abs((p.at<float>(di, pi) / w) - center.at<float>(di, 0));
		}
	}
	scale /= p.cols;

	// build condition matrix
	Mat cond = Mat::eye(D, D, p.type());
	for (int di = 0; di < (D - 1); ++di) {
		cond.at<float>(di, di) = 1.0f / scale.at<float>(di, 0);
		cond.at<float>(di, D - 1) = - center.at<float>(di, 0) / scale.at<float>(di, 0);
	}

	return cond;
}

/**
 * get the conditioning matrix of given points
 * @param p  the points as matrix
 * @return  the condition matrix (already allocated)
 */
Mat getCondition2D(Mat& p) {

	return getConditionXD(p);
}

/**
 * get the conditioning matrix of given 3D points
 * @param p  the points as matrix
 * @param T  the condition matrix (already allocated)
 */
Mat getCondition3D(Mat& p) {

	return getConditionXD(p);
}

/**
 * enforce rank of 2 on fundamental matrix
/*
F  the matrix to be changed
*/
void forceSingularity(Mat& F) {

	Mat U(3, 3, CV_32FC1);
	Mat d(3, 1, CV_32FC1);
	Mat Vt(3, 3, CV_32FC1);

	// decompose
	SVD::compute(F, d, U, Vt);

	// set d3 to 0  (enforce rank 2)
	d.at<float>(2, 0) = 0;

	// recombine
	F = U * Mat::diag(d) * Vt;
}

void setSmallValuesToZero(Mat& mat) {

	for (int ri = 0; ri < mat.rows; ++ri) {
		for (int ci = 0; ci < mat.cols; ++ci) {
			if (abs(mat.at<float>(ri, ci)) < 1.0e-5) {
				mat.at<float>(ri, ci) = 0.0f;
			}
		}
	}
}

void savePointListOFF(string fname, Mat& points) {

	// open file for write
	fstream file(fname.c_str(), ios::out);
	if (!file) {
		cerr << "ERROR: cannot open file " << fname << endl;
		return;
	}

	// if homogeneous points: normalize
	if (points.rows == 4) {
		for (int i=0; i < points.cols; i++) {
			points.at<float>(0, i) /= points.at<float>(3, i);
			points.at<float>(1, i) /= points.at<float>(3, i);
			points.at<float>(2, i) /= points.at<float>(3, i);
			points.at<float>(3, i) /= points.at<float>(3, i);
		}
	}
	setSmallValuesToZero(points);
	file << "OFF" << endl;
	file << points.cols << " 0" << " 0" << endl;
	// write euclidian vertices
	for (int i=0; i < points.cols; i++) {
		file
				<< points.at<float>(0, i) << ","
				<< points.at<float>(1, i) << ","
				<< points.at<float>(2, i) << endl;
	}
	// write faces (vertex indices)
	/*for (int i=0; i < points.cols; i++) {
		file
				<< i << " "
				<< (i + 1) % points.cols << " "
				<< (i + 2) % points.cols << endl;
	}*/

	// close file
	file.close();
}

float sumMembers(const Mat& mat) {

	return Mat(
			Mat::ones(1, mat.rows, mat.type())
			* mat
			* Mat::ones(mat.cols, 1, mat.type())
			).at<float>(0, 0);
}

float absDifferenceSum(const Mat& mat1, const Mat& mat2) {
	return sumMembers(abs(mat1 - mat2));
}

/* ***********************
   *** Given Functions ***
   *********************** */

/**
 * saves point list to file
 * @param fname  file name
 * @param points  matrix of points
 */
void savePointList(string fname, Mat& points) {

	// open file for write
	fstream file(fname.c_str(), ios::out);
	if (!file) {
		cerr << "ERROR: cannot open file " << fname << endl;
		return;
	}

	// if homogeneous points: normalize
	if (points.rows == 4) {
		for (int i=0; i < points.cols; i++) {
			points.at<float>(0, i) /= points.at<float>(3, i);
			points.at<float>(1, i) /= points.at<float>(3, i);
			points.at<float>(2, i) /= points.at<float>(3, i);
			points.at<float>(3, i) /= points.at<float>(3, i);
		}
	}
	setSmallValuesToZero(points);
	// write euclidian points
	for (int i=0; i < points.cols; i++) {
		file
				<< points.at<float>(0, i) << ","
				<< points.at<float>(1, i) << ","
				<< points.at<float>(2, i) << endl;
	}

	// close file
	file.close();
}

/**
 * read homologeous points from file
 * @param fname  path of file containing point list
 * @param x1  pointer to matrix containing points of first image
 * @param x2  pointer to matrix containing points of second image
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

/**
 * read control points from file
 * @param fname  path of file containing point list
 * @param x1  pointer to matrix containing points of first image
 * @param x2  pointer to matrix containing points of second image
 * @param Xm  pointer to matrix containing object points
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
