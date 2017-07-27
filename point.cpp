#include "point.h"


Point::Point()
{
   
}


Point::~Point()
{
}

Point::Point(int x,int y)
{
	this->x=x;
	this->y=y;
}


void Point::setLocation(int x, int y)
{
	this->x=x;
	this->y=y;
}

 Point& Point::operator=(const Point &p) 
 {
	 this->x=p.x;
	 this->y=p.y;

    return *this;  // Return a reference to myself.
 }

 Point operator+(const Point &p1,const Point &p2)
 {
	 return Point(p1.x+p2.x,p1.y+p2.y);
 }
 Point operator-(const Point &p1,const Point &p2)
 {
	 return Point(p1.x-p2.x,p1.y-p2.y);
 }
 
bool operator== (const Point &P1, const Point &P2)
{
	return (P1.x == P2.x && P1.y == P2.y);
}
 
bool operator!= (const Point &P1, const Point &P2)
{
    return !(P1 == P2);
}
