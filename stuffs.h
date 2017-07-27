#pragma once

#ifndef STUFFS_H
#define STUFFS_H

#include<iostream>
#include <iomanip>
//#include <QImage>
#include <math.h>
#include "point.h"
#include <list>
#include <QTextStream>
#include <QFile>

typedef std::list<Point> PointList;

////Include libraries for image processing
//#include "itkImage.h"
//#include "itkImageFileReader.h"
//#include "itkImageFileWriter.h"
//#include <itkGradientMagnitudeImageFilter.h>
//#include "itkRescaleIntensityImageFilter.h"
//#include "itkApproximateSignedDistanceMapImageFilter.h"
//#include "itkGradientVectorFlowImageFilter.h"
//#include "itkGradientImageFilter.h"
//#include "itkVector.h"

//typedef unsigned char						UnsignedCharPixelType;
//typedef float								FloatPixelType;
//typedef itk::Image<FloatPixelType, 2>		FloatImageType;
//typedef itk::Image<UnsignedCharPixelType,2>	UnsignedCharImageType;
//typedef FloatImageType						ImageType;

void print(const char * title, double** f, int width, int height);
void print3D(const char * title, double*** f, int width, int height);
double** newArray(int width, int height);
double*** newArray3D(int x, int y, int z);
double*** GVFC(double** f, int w, int h, int ITER, double mu);
double*** GVF(double** f, int w, int h, int ITER, double mu);
double** GVFMagnitude(double*** gvf, int w, int h);
void save(double** im, int width, int height, const char* file);
double** normalizedImage(const char* file, int &width, int &height);
double** normalizedEdgeMap(const char* file, int &width, int &height);
void saveVectorImage(double*** gvf, int w, int h);
double distance2D(Point A, Point B);
PointList getSegmentation(PointList l);
Point* listToPointArray(PointList l, int &N);
double distance(int i, Point* S, int sizeS, Point* T, int sizeT);
double Pratt2(Point* S, int sizeS, Point* T, int sizeT);
void saveToFile(PointList l, QString file);
void saveToFile(QString l, QString fileName);

#endif
