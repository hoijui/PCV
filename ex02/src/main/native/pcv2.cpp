//============================================================================
// Name        : pcv2.cpp
// Author      : Ronny Haensch
// Version     : 1.0
// Copyright   : -
// Description : loads three images, computes homographies,
//               and stiches the images together
//============================================================================

#include <iostream>
#include "opencv2/opencv.hpp"

using namespace std;
using namespace cv;

// global variable to save point pairs
vector<Point2f> pointList;

// some function headers (see below)
// given functions
int getPoints(struct winInfo& img1, struct winInfo& img2, Mat& p1, Mat& p2);
void getPoints(int event, int x, int y, int flags, void* param);
Mat stich(Mat& img1, Mat& img2, Mat& H);
// functions to be implemented
Mat homography2D(Mat& p1, Mat& p2);
void decondition(Mat& T1, Mat& T2, Mat& H);
Mat getDesignMatrix_homography2D(Mat& n1, Mat& n2);
Mat transform(Mat& H, Mat& p);
Mat getCondition2D(Mat& p);
Mat solve_dlt(Mat& A);

/**
 * structure to fuse image data and window title
 */
struct winInfo {Mat img; string name;};

/**
 * usage: path to base image in argv[1], path to other images in argv[2] and argv[3]
 */
int main(int argc, char** argv) {

    // check if image paths were defined
	if (argc != 4) {
	cerr << "Usage: pcv2 <path to base image> <path to 2nd image> <path to 3rd image>" << endl;
	return -1;
    }

    // titles of some windows
    string winName1 = string("Base image");
    string winName2 = string("Image to attach");
    string winName3 = string("Panorama");

    // load image first two images, paths in argv[1] and argv[2]
    Mat baseImage = imread(argv[1]);
    Mat attachImage = imread(argv[2]);

    // fuse image data and window title
    struct winInfo base;
    base.img = baseImage.clone();
	base.name = winName1;
    struct winInfo attach;
    attach.img = attachImage.clone();
    attach.name = winName2;

    // get corresponding points within the two image
    // start with one point within the attached image, then click on corresponding point in base image
    Mat p_basis, p_attach;
    int numberOfPointPairs = getPoints(base, attach, p_basis, p_attach);

    // just some putput
    cout << "Number of defined point pairs: " << numberOfPointPairs << endl;
    cout << endl << "Points in base image:" << endl;
    cout << p_basis << endl;
    cout << endl << "Points in second image:" << endl;
    cout << p_attach << endl;

    // calculate homography
    Mat H = homography2D(p_basis, p_attach);

    // create panorama
    Mat panorama = stich(baseImage, attachImage, H);

    // display panorama (resizeable)
	namedWindow(winName3.c_str(), 0);
	imshow(winName3.c_str(), panorama);
    waitKey(0);
    destroyWindow( winName3.c_str());

    // clear everything for the next image
    pointList.clear();

    // panorama is new base image, third image is the image to attach
    baseImage = panorama;
    // load third image
    attachImage = imread(argv[3]);

    // fuse image data and window title
    base.img = baseImage.clone();
    attach.img = attachImage.clone();

    // get corresponding points within the two image
    // start with one point within the attached image, then click on corresponding point in base image
    numberOfPointPairs = getPoints(base, attach, p_basis, p_attach);

    // just some putput
    cout << "Number of defined point pairs: " << numberOfPointPairs << endl;
    cout << endl << "Points in base image:" << endl;
    cout << p_basis << endl;;
    cout << endl << "Points in second image:" << endl;
    cout << p_attach << endl;;

    // calculate homography
    H = homography2D(p_basis, p_attach);

    // create panorama
    panorama = stich(baseImage, attachImage, H);

    // display panorama (resizeable)
    namedWindow( winName3.c_str(), 0 );
    imshow( winName3.c_str(), panorama );
    waitKey(0);
    destroyWindow( winName3.c_str());

    imwrite("panorama.png", panorama);

	return 0;
}

/**
 * compute homography
 * @param to  first set of points
 * @param from  second set of points
 * @param H  the homography to be computed
*/
Mat homography2D(Mat& base, Mat& attach) {

    // coordinate transformation matrices for conditioning
    Mat tBase = getCondition2D(base);
    Mat tAttach = getCondition2D(attach);

    // create condition matrices
    Mat cBase = transform(tBase, base);
    Mat cAttach = transform(tAttach, attach);

    // create design matrix, for projection attach -> base
    Mat design = getDesignMatrix_homography2D(cBase, cAttach);

    // svd and reshaping
    Mat H = solve_dlt(design);
    
    // creating final homography
    decondition(tBase, tAttach, H);

    return H;
}

/**
 * solve homogeneous equation system by usage of SVD
 * @param A  the design matrix
 * @param H  the estimated homography
 */
Mat solve_dlt(Mat& A) {
    
    Mat d, u, vt;
    SVD::compute(A, d, u, vt, 0); // gibts noch einen kürzeren befehl?!?

    Mat h = vt.row(8);

    // reshape h -> H
    Mat H = Mat(3, 3, CV_32FC1);
    for(int x = 0; x < 3; ++x)
    {
        for(int y = 0; y < 3; ++y)
        {
            H.at<float>(y, x) = h.at<float>(0, y * 3 + x);
        }
    }

    return H;
}

/**
 * decondition a homography that was estimated from conditioned point clouds
 * @param T_to  conditioning matrix of first set of points
 * @param T_from  conditioning matrix of second set of points
 * @param H  conditioned homography that has to be un-conditioned (in-place)
 */
void decondition(Mat& T_base, Mat& T_attach, Mat& H) {

    H = T_base.inv() * H * T_attach;
}

/**
 * define the design matrix as needed to compute 2D-homography
 * @param base  first set of points:	attach = H * base
 * @param base  second set of points:	attach = H * base
 * @param A  the design matrix to be computed
 */
Mat getDesignMatrix_homography2D(Mat& base, Mat& attach) {

    // design matrix: at least 9x9, if more points selected, then 2*pointsx9
    Mat designMat = Mat::zeros(base.cols == 4 ? 9 : base.cols * 2, 9, CV_32FC1);

    for(int i = 0; i < base.cols; ++i)
    {
        int y = i * 2;

        // zweimal: -w' * x
        designMat.at<float>(y + 1, 0 + 3) = designMat.at<float>(y, 0) = -base.at<float>(2, i) * attach.at<float>(0, i);  // -w' * x
        designMat.at<float>(y + 1, 1 + 3) = designMat.at<float>(y, 1) = -base.at<float>(2, i) * attach.at<float>(1, i);  // -w' * y
        designMat.at<float>(y + 1, 2 + 3) = designMat.at<float>(y, 2) = -base.at<float>(2, i) * attach.at<float>(2, i);  // -w' * w

        // u' * x
        designMat.at<float>(y, 6) = base.at<float>(0, i) * attach.at<float>(0, i);
        designMat.at<float>(y, 7) = base.at<float>(0, i) * attach.at<float>(1, i);
        designMat.at<float>(y, 8) = base.at<float>(0, i) * attach.at<float>(2, i);

        // v' * x
        designMat.at<float>(y + 1, 6) = base.at<float>(1, i) * attach.at<float>(0, i);
        designMat.at<float>(y + 1, 7) = base.at<float>(1, i) * attach.at<float>(1, i);
        designMat.at<float>(y + 1, 8) = base.at<float>(1, i) * attach.at<float>(2, i);
    }

    return designMat;
}

/**
 * apply transformation to set of points
 * @param H  matrix representing the transformation
 * @param p  input points
 * @param n  transformed points
 */
Mat transform(Mat& H, Mat& p) {

    Mat transformed(3, p.cols, CV_32FC1);

    for(int i = 0; i < p.cols; ++i)
    {
        Mat tPoint = H * p.col(i);

        transformed.at<float>(0, i) = tPoint.at<float>(0, 0);
        transformed.at<float>(1, i) = tPoint.at<float>(1, 0);
        transformed.at<float>(2, i) = tPoint.at<float>(2, 0);
    }

   return transformed;
}

/**
 * get the conditioning matrix of given points
 * @param p  the points as matrix
 * @param T  the condition matrix (already allocated)
 */
Mat getCondition2D(Mat& p) {

    // calculate center
    float transX = 0,
    transY = 0;

    for(int i = 0; i < p.cols; ++i)
    {
        transX += p.at<float>(0, i) / p.at<float>(2, i);
        transY += p.at<float>(1, i) / p.at<float>(2, i);
    }

    transX /= p.cols;
    transY /= p.cols;

    // calculate scale
    float scaleX = 0,
        scaleY = 0;

    for(int i = 0; i < p.cols; ++i)
    {
        scaleX += abs((p.at<float>(0, i) / p.at<float>(2, i)) - transX);
        scaleY += abs((p.at<float>(0, i) / p.at<float>(2, i)) - transY);
    }

    scaleX /= 4;
    scaleY /= 4;

    // build condition matrix
    Mat cond = Mat::zeros(3, 3, CV_32FC1);
    cond.at<float>(0, 0) = 1 / scaleX;
    cond.at<float>(0, 2) = -transX / scaleX;
    cond.at<float>(1, 1) = 1 / scaleY;
    cond.at<float>(1, 2) = -transY / scaleY;
    cond.at<float>(2, 2) = 1;

    return cond;
}

/* *************************
   **** Given Functions ****
   ************************* */

/**
 * stitch two images together by transforming one of them by a given homography
 * @param base		the base image
 * @param attach		the image to be attached
 * @param H		the homography to warp the second image
 * @param panorama	the resulting image
 */
Mat stich(Mat& base, Mat& attach, Mat& H) {

    // compute corners of warped image
    Mat corners(1, 4, CV_32FC2);
    corners.at<Vec2f>(0, 0) = Vec2f(0,0);
    corners.at<Vec2f>(0, 1) = Vec2f(0,attach.rows);
    corners.at<Vec2f>(0, 2) = Vec2f(attach.cols,0);
    corners.at<Vec2f>(0, 3) = Vec2f(attach.cols,attach.rows);
	perspectiveTransform(corners, corners, H);

	// compute size of resulting image and allocate memory
    float x_start = min( min( corners.at<Vec2f>(0, 0)[0], corners.at<Vec2f>(0, 1)[0]), (float)0);
    float x_end   = max( max( corners.at<Vec2f>(0, 2)[0], corners.at<Vec2f>(0, 3)[0]), (float)base.cols);
    float y_start = min( min( corners.at<Vec2f>(0, 0)[1], corners.at<Vec2f>(0, 2)[1]), (float)0);
    float y_end   = max( max( corners.at<Vec2f>(0, 1)[1], corners.at<Vec2f>(0, 3)[1]), (float)base.rows);

    // create translation matrix in order to copy both images to correct places
    Mat T = Mat::zeros(3,3,CV_32FC1);
    T.at<float>(0, 0) = 1;
    T.at<float>(1, 1) = 1;
    T.at<float>(2, 2) = 1;
    T.at<float>(0, 2) = -x_start;
    T.at<float>(1, 2) = -y_start;

    // change homography to take necessary translation into account
    T = T * H;
    // warp second image and copy it to output image
    Mat panorama;
    warpPerspective(attach, panorama, T, Size(x_end - x_start + 1, y_end - y_start + 1), CV_INTER_LINEAR);

    // copy base image to correct position within output image
    Mat roi(panorama, Rect(-x_start,-y_start,base.cols, base.rows));
    base.copyTo(roi, base);

    return panorama;
}

/**
 * display two images and catch the point pairs marked by left mouse clicks
 * points will be in homogeneous coordinates
 * @param base		structure containing base image
 * @param attach		structure containing image that has to be attached
 * @param p_base		points within the base image (to be defined by this method)
 * @param p_attach	points within the second image (to be defined by this method)
 */
int getPoints(struct winInfo& base, struct winInfo& attach, Mat& p_base, Mat& p_attach) {

    cout << endl;
    cout << "Please select at least four points by clicking at the corresponding image positions:" << endl;
    cout << "Firstly click at the point that shall be transformed (within the image to be attached), followed by a click on the corresponding point within the base image" << endl;
    cout << "Continue until you have collected as many point pairs as you wish" << endl;
    cout << "Stop the point selection by pressing any key" << endl << endl;

    // show input images and install mouse callback
    namedWindow( base.name.c_str(), 0 );
    imshow( base.name.c_str(), base.img );
    setMouseCallback(base.name.c_str(), getPoints, (void*) &base);
    namedWindow( attach.name.c_str(), 0 );
    imshow( attach.name.c_str(), attach.img );
    setMouseCallback(attach.name.c_str(), getPoints, (void*) &attach);
    // wait until any key was pressed
    waitKey(0);

    destroyWindow( base.name.c_str() );
    destroyWindow( attach.name.c_str() );

    // allocate memory for point-lists (represented as matrix)
    p_base = Mat(3, pointList.size()/2, CV_32FC1);
    p_attach = Mat(3, pointList.size()/2, CV_32FC1);
    int n=0; // number of point pairs
    // read points from global variable, transform them into homogeneous coordinates
	for (vector<Point2f>::iterator p = pointList.begin(); p != pointList.end(); p++) {
	p_attach.at<float>(0, n) = (*p).x;
	p_attach.at<float>(1, n) = (*p).y;
	p_attach.at<float>(2, n) = 1;
	p++;
	p_base.at<float>(0, n) = (*p).x;
	p_base.at<float>(1, n) = (*p).y;
	p_base.at<float>(2, n) = 1;
	n++;
    }

    return n;
}

/**
 * mouse call back to get points and draw circles
 * @param event	specifies encountered mouse event
 * @param x,y	position of mouse pointer
 * @param flags	not used here
 * @param param	a struct containing used IplImage and window title
 */
void getPoints(int event, int x, int y, int flags, void* param) {

  // cast to structure
  struct winInfo* win = (struct winInfo*)param;

	switch (event) {
    // if left mouse button was pressed
		case CV_EVENT_LBUTTONDOWN: {
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
