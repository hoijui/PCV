//============================================================================
// Name        : pcv5.cpp
// Author      : Ronny Haensch
// Version     : 1.0
// Copyright   : -
// Description : loads two images, computes relative orientation
// Group       : Marcus Grum (), Robin Vobruba (), Jens Jawer () und Marcus Pannwitz (343479)
//============================================================================

#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

// global variable to save point pairs
vector<Point2f> pointList;

// some function headers (see below)
Mat getFundamentalMatrix(Mat& p1, Mat& p2);
Mat getCondition2D(Mat& p);
Mat transform(Mat& H, Mat& p);
Mat getDesignMatrix_fundamental(Mat& n1, Mat& n2);
void forceSingularity(Mat& F);
Mat solve_dlt(Mat& A);
void decondition(Mat& T1, Mat& T2, Mat& F);
void visualize(struct winInfo& img1, struct winInfo& img2, Mat& p_fst, Mat& p_snd, Mat& F);
double getError(Mat& p_fst, Mat& p_snd, Mat& F);
// given functions
int getPoints(struct winInfo& img1, struct winInfo& img2, Mat& p1, Mat& p2);
void getPoints(int event, int x, int y, int flags, void* param);
void drawEpiLine(Mat& img, double a, double b, double c);

// structure to fuse image data and window title
struct winInfo {Mat img; string name;};

// usage: path to first image in argv[1], path to second image in argv[2]
// main function
int main(int argc, char** argv) {

    // check if image paths were defined
    if (argc != 3){
	cerr << "Usage: pcv5 <path to 1st image> <path to 2nd image>" << endl;
	return -1;
    }

    // titles of some windows
    string winName1 = string ("First image");
    string winName2 = string ("Second image");

    // load images, paths in argv[1] and argv[2]
    Mat fstImage = imread(argv[1]);
    Mat sndImage = imread(argv[2]);
    
    if ( (!fstImage.data) || (!sndImage.data)){
	cerr << "ERROR: Could not load images" << endl;
	return -2;
    }

    // fuse image data and window title
    struct winInfo fst;
    fst.img = fstImage.clone();
    fst.name = winName1;    
    struct winInfo snd;
    snd.img = sndImage.clone();
    snd.name = winName2;
   
    // get corresponding points within the two images
    // start with one point within the first image, then click on corresponding point in second image
    Mat p_fst, p_snd;

// comment out if you want to click
#define PREDEFINED_POINTS 1

#ifdef PREDEFINED_POINTS
	p_fst = (Mat_<float>(3, 8, CV_32FC1) <<
        84, 61, 64, 140, 190, 139, 63, 43,
        19, 31, 149, 136, 160, 207, 245, 201,
        1, 1, 1, 1, 1, 1, 1, 1);

	p_snd = (Mat_<float>(3, 8, CV_32FC1) <<
        98, 84, 80, 156, 197, 152, 77, 45,
        5, 19, 133, 129, 150, 193, 235, 187,
        1, 1, 1, 1, 1, 1, 1, 1);

    namedWindow( fst.name.c_str(), 0 );
    imshow( fst.name.c_str(), fst.img );
    namedWindow( fst.name.c_str(), 0 );
    imshow( snd.name.c_str(), snd.img );

    int numberOfPointPairs = 8;
#else
	int numberOfPointPairs = getPoints(fst, snd, p_fst, p_snd);
#endif
    
    // just some putput
    cout << "Number of defined point pairs: " << numberOfPointPairs << endl;
    cout << endl << "Points in first image:" << endl;
    cout << p_fst << endl;
    cout << endl << "Points in second image:" << endl;
    cout << p_snd << endl;
    
    // calculate fundamental matrix
    Mat F = getFundamentalMatrix(p_fst, p_snd);

    // visualize epipolar lines
    visualize(fst, snd, p_fst, p_snd, F);

    // calculate geometric error
    double err = getError(p_fst, p_snd, F);
    cout << "Geometric error: " << err << endl;
    
    return 0;

}

// compute fundamental matrix
/*
fst	first set of points
snd	second set of points
return	the estimated fundamental matrix
*/
Mat getFundamentalMatrix(Mat& fst, Mat& snd){

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
	decondition(tFirst, tSecond, F);

	return F;
}

// solve homogeneous equation system by usage of SVD
/*
A		the design matrix
return		the estimated fundamental matrix
*/
Mat solve_dlt(Mat& A) {
    
	Mat f = Mat::zeros(1, 9, CV_32FC1);
	SVD::solveZ(A, f);
	f = f.reshape(1, 3);
//f = f.reshape(0, 3);
	return f;
}

// decondition a fundamental matrix that was estimated from conditioned point clouds
/*
T_fst	conditioning matrix of first set of points
T_snd	conditioning matrix of second set of points
F	conditioned fundamental matrix that has to be un-conditioned (in-place)
*/
void decondition(Mat& T_fst, Mat& T_snd, Mat& F){
  
    //F = T_fst.inv() * F * T_snd;
    F = T_snd.t() * F * T_fst;
}

// define the design matrix as needed to compute fundamental matrix
/*
fst		first set of points
snd		second set of points
return		the design matrix to be computed
*/
Mat getDesignMatrix_fundamental(Mat& fst, Mat& snd){

	// The design matrix has at least 9 rows in case of 8 points. If there are more points, the number of rows equals the number of points
	Mat design = (fst.cols < 8) ? Mat::zeros(9, 9, CV_32FC1) : Mat::zeros(fst.cols, 9, CV_32FC1);

	// Loop through all corresponding points
	for (int i = 0; i < fst.cols; ++i)
	{
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

	return design;
}

// apply transformation to set of points
/*
H		matrix representing the transformation
p		input points
return		transformed points
*/
Mat transform(Mat& H, Mat& p) {
	return H * p;
}

// get the conditioning matrix of given points
/*
p		the points as matrix
return		the condition matrix (already allocated)
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

// enforce rank of 2 on fundamental matrix
/*
F	the matrix to be changed
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

// draws epipolar lines in image
/*
img     target window
points  points in other image
F       fundamental matrix between both images
*/
void visualizeHelper(struct winInfo& win, Mat& points, Mat F) {

	// for each point
	for (int i = 0; i < points.cols; ++i)
	{
		// l' = F * x or l = F^T * x'
		const Mat& x = points.col(i);
		const Mat line = F * x;
		double a = line.at<float>(0, 0);
		double b = line.at<float>(1, 0);
		double c = line.at<float>(2, 0);

		drawEpiLine(win.img, a, b, c);
	}

	// update image in windows
	imshow(win.name.c_str(), win.img);
}

// draws epipolar lines into images
/*
img1	structure containing fst image
img2	structure containing snd image
p_fst	first point set
p_snd	second point set
F	fundamental matrix
*/
void visualize(struct winInfo& img1, struct winInfo& img2, Mat& p_fst, Mat& p_snd, Mat& F){

    // draw epipolar lines for p_fst in img2
    visualizeHelper(img2, p_fst, F);

    // draw epipolar lines for p_snd in img1
    visualizeHelper(img1, p_snd, F.t());

    // save Windows for Doc
    imwrite("EpipolarLinesL.png", img1.img);
    imwrite("EpipolarLinesR.png", img2.img);
    // wait until any key was pressed
    waitKey(0);
}

// calculate geometric error of estimated fundamental matrix
/*
p_fst		first set of points
p_snd		second set of points
F		fundamental matrix
return		geometric error
*/
double getError(Mat& p_fst, Mat& p_snd, Mat& F) {

	double sum = 0;
	for(int i = 0; i < p_fst.cols; ++i)
	{
		const Mat& x = p_fst.col(i);
		const Mat& x_ = p_snd.col(i);
		const Mat F_x = F*x;
		const Mat F_T_x_ = F.t()*x_;

		sum = sum + pow(Mat(x_.t() * F_x).at<float>(0, 0), 2) /
				( pow(F_x.at<float>(0, 0), 2)
					+ pow(F_x.at<float>(1, 0), 2)
					+ pow(F_T_x_.at<float>(0, 0), 2)
					+ pow(F_T_x_.at<float>(1, 0), 2)
				);
	}
	int N=p_fst.cols;

	return 1.0/N*sum;
}

/* ***********************
   *** Given Functions ***
   *********************** */

// draw line given in homogeneous representation into image
/*
img	the image
a,b,c	the line parameters
*/
void drawEpiLine(Mat& img, double a, double b, double c){

  // calculate intersection with image borders
  double x_0, x_1, y_0, y_1;
  Point p1 = Point(-c/a, 0);						// schnittpunkt mit unterer bildkante (x-achse)
  Point p2 = Point(0, -c/b);						// schnittpunkt mit linker bildkante (y-achse)
  Point p3 = Point((-b*(img.rows-1)-c)/a, img.rows-1);		// schnittpunkt mit oberer bildkante
  Point p4 = Point(img.cols-1, (-a*(img.cols-1)-c)/b);		// schnittpunkt mit rechter bildkante
  
  // check start and end points
  Point startPoint, endPoint, cur_p;
  startPoint.x = startPoint.y = endPoint.x = endPoint.y = 0;
  bool set_start = false, set_end = false;
  for(int p=0; p<4; p++){
    switch(p){
      case 0: cur_p = p1; break;
      case 1: cur_p = p2; break;
      case 2: cur_p = p3; break;
      case 3: cur_p = p4; break;
    }
    if ( (cur_p.x >= 0) && (cur_p.x < img.cols) && (cur_p.y >= 0) && (cur_p.y < img.rows) ){
      if (!set_start){
	startPoint = cur_p;
	set_start = true;
      }else{
	endPoint = cur_p;
	set_end = true;
      }
    }
  }
  
  // draw line
  line(img, startPoint, endPoint, Scalar(0,0,255), 1);

}

// display two images and catch the point pairs marked by left mouse clicks
// points will be in homogeneous coordinates
/*
fst		structure containing fst image
snd		structure containing image that has to be snded
p_fst		points within the fst image (to be defined by this method)
p_snd		points within the second image (to be defined by this method)
*/
int getPoints(struct winInfo& fst, struct winInfo& snd, Mat& p_fst, Mat& p_snd){

    // show input images and install mouse callback
    namedWindow( fst.name.c_str(), 0 );
    imshow( fst.name.c_str(), fst.img );
    setMouseCallback(fst.name.c_str(), getPoints, (void*)(&fst));
    namedWindow( snd.name.c_str(), 0 );
    imshow( snd.name.c_str(), snd.img );
    setMouseCallback(snd.name.c_str(), getPoints, (void*)(&snd));
    // wait until any key was pressed
    waitKey(0);

    // allocate memory for point-lists (represented as matrix)
    p_fst = Mat(3, pointList.size()/2, CV_32FC1);
    p_snd = Mat(3, pointList.size()/2, CV_32FC1);
    int n=0;	// number of point pairs
    // read points from global variable, transform them into homogeneous coordinates
    for(vector<Point2f>::iterator p = pointList.begin(); p != pointList.end(); p++){
	p_fst.at<float>(0, n) = p->x;
	p_fst.at<float>(1, n) = p->y;
	p_fst.at<float>(2, n) = 1;
	p++;
	p_snd.at<float>(0, n) = p->x;
	p_snd.at<float>(1, n) = p->y;
	p_snd.at<float>(2, n) = 1;
	n++;
    }

    // save Windows for Doc
    imwrite("PointCorrespondencesL.png", fst.img);
    imwrite("PointCorrespondencesR.png", snd.img);

    // close windows
    destroyWindow( fst.name.c_str() );
    destroyWindow( snd.name.c_str() );

    return n;

}

// mouse call back to get points and draw circles
/*
event	specifies encountered mouse event
x,y	position of mouse pointer
flags	not used here
param	a struct containing used IplImage and window title
*/
void getPoints(int event, int x, int y, int flags, void* param){

  // cast to structure
  struct winInfo* win = (struct winInfo*)param;
  
  switch(event){
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
    }break;
  }
}
