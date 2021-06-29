#pragma once

#ifndef SNAKE_H
#define SNAKE_H

#include <QImage>
//#include "/media/Khue/Dev. Software/C++libraries/boost/boost/utility.hpp"
#include <list>
#include "point.h"
#include <math.h>
#include <complex>
#include <iostream>
#include "stuffs.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2_algorithms.h>

//Include libraries for image processing
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <itkGradientMagnitudeImageFilter.h>
#include "itkRescaleIntensityImageFilter.h"
#include "itkApproximateSignedDistanceMapImageFilter.h"
#include "itkGradientVectorFlowImageFilter.h"
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkFastChamferDistanceImageFilter.h>
#include "itkThresholdImageFilter.h"

typedef unsigned char						UnsignedCharPixelType;
typedef float								FloatPixelType;
typedef itk::Image<FloatPixelType, 2>		FloatImageType;
typedef itk::Image<UnsignedCharPixelType,2>	UnsignedCharImageType;
typedef FloatImageType						ImageType;

typedef std::list<Point> PointList;
typedef std::complex<double> complex;

enum DistanceType{EUCLIDEAN, DTW, PRATT,GDM,LARSSON};

//complex* getFourierCoeffiecients2(PointList snake);
//PointList getShapeFromFourierCoeffiecients(complex* input,int N);
complex* DFT(complex* input, int N);		//Discrete Fourier Transform
complex* nDFT(PointList snake,bool SPI, bool RotInvariant);	//Normalized Discrete Fourier Transform
complex* iDFT(complex* input, int N);		//Inverse Discrete Fourier Transform	
PointList iDFT2Points(complex* input,int N);	//Convert iDFT to a list of points
complex* PointsToDFT(PointList points,bool SPInvariant, bool RotInvariant);
complex* points2complex(PointList points);
PointList complex2points(complex *z,int N);
double dtwDistance(complex *s, complex *t, int N);
double dtwDistance(complex *s, complex *t, int N, int wLenght);
double dtwDistance2(complex *s, complex *t, int N);
double LB_Keogh(complex *q, complex *c, int n, int r);
double Pratt(complex* S, int sizeS, complex* T, int sizeT);
double PartitionDistance(complex* S, int sizeS, complex* T, int sizeT);

complex* fillPoints(complex* c, int n, int N);
double maxR12(complex* S, int sizeS, complex* T, int sizeT);
complex max(complex *q, int first, int last);
double interArea(complex* S, int sizeS, complex* T, int sizeT);

class Snake 
{
/* Private Attributes */
private:
	PointList snake;
	PointList shape;
	PointList snakeU;
	PointList snakeL;
	int width;
	int height;
	double snakelength;
    double** gradient;
    double** flow;
	double** e_uniformity ;
	double** e_curvature  ;
	double** e_image    ;
	double** e_balloon;
	double** e_flow;
	double** e_prior;
	
	/*complex* shapeCoeff1;
	complex* snakeCoeff1;
	complex* shapeCoeff;
	complex* snakeCoeff;*/


    //the DFT and iDFT of the shape
//    complex* shWithSPI;
//    complex* ShWithSPI;
//    complex* shWithoutSPI;
//    complex* ShWithoutSPI;

public:
	bool useSPI;
    bool useRotInvar;
	bool useDTW;
	int windowLength;
    DistanceType distanceType;
    DistanceType distanceTypeSP;
/*Public Methods*/
public:
	//Constructors and destructor
	Snake();
    Snake(int width, int height, double** gradient, double** flow, PointList points, PointList shape);
    Snake(const char* fileName);
    Snake(const char* fileName, PointList shape);
    Snake(const char* fileName, double sigma, double threshold);
    Snake(const char* fileName, PointList shape, double sigma, double threshold);
	~Snake();

	//Operators
	Snake& operator=(const Snake &s); 
	Snake replace(Point oldpoint,Point newpoint);

    void initGradientFlow(const char* fileName);
    void initGradientFlow(const char* fileName, double sigma, double threshold);
    void setGradientFlow(QImage imageGrad, QImage imageFlow);
	
	bool step1(double continuityCoeff,double curvatureCoeff,double imageCoeff,double flowCoeff,double balloonCoeff); 
	bool step2(double continuityCoeff,double curvatureCoeff,double imageCoeff,double flowCoeff,double balloonCoeff,double priorCoeff); 
	
	double getsnakelength();
	
	complex* DFTSnake();
	complex* DFTShape();
	inline PointList getSnake(){return snake;}
	inline PointList getShape(){return shape;}
	PointList getSnakeU();
	PointList getSnakeL();
    inline double** getGradient(){return gradient;}
	inline void setSnake(PointList list){this->snake.clear() ; this->snake=list; }    
    inline void setShape(PointList s){ this->shape.clear(); this->shape=s; /*sh=this->DFTShape(); Sh=iDFT(sh,(int)this->shape.size());*/}
	//inline void enableSPI(bool b){useSPI = b;}
	inline void enableDTW(bool b){useDTW = b;}
	Point get(int i);

	//Auto-adaptation
	void rebuild(int space);
    void rebuild2(int n);
	void removeOverlappingPoints(int minlen);
    int removeOverlappingPoints2(int minlen);
	void addMissingPoints(int maxlen);
	void moveAwayPoints(int maxlen); 
    int getOptimalStartingPoint();
    double getOptimalRotation(double angleStep);
	PointList shiftSnake(int startingIndex);
    PointList rotateSnake(double angle);
    void shift(int index);
    void rotate(double angle);    //rotate all points of the snake by ¨angle¨ degrees
    void reSample(int n);

    double EuclideanModified();

/*Private Methods*/
private:	
	double f_uniformity(Point prev, Point next, Point p);
	double f_curvature(Point prev, Point p, Point next);
	//double f_gflow(Point cur, Point p);
	double f_image(Point cur, Point p);
	double f_balloon(Point prev, Point p, Point cur, Point next);
	double f_gflow(Point cur, Point p);
	double f_prior(); 
	void dynamicAllocation();
	void freeAllocation();
	
};

#endif
