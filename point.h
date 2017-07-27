#pragma once
#ifndef POINT_H
#define POINT_H

#include <math.h>


class Point
{
private:
	int x;
	int y;
public:
	Point();
	~Point();
	Point(int x,int y);
	int getX() {return x;}
	int getY() {return y;}
	void setX(int X){x=X;}
	void setY(int Y){y=Y;}
	void setLocation(int x, int y);
    double norm(){return sqrt((double)x*x+y*y);}

	Point& operator=(const Point &p);
	friend Point operator+(const Point &p1,const Point &p2);
	friend Point operator-(const Point &p1,const Point &p2);
	friend bool operator== (const Point &P1, const Point &P2);
    friend bool operator!= (const Point &P1, const Point &P2);
};


#endif
