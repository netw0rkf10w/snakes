#ifndef SCENE_H
#define SCENE_H

#include "point.h"
#include <list>
#include <QtGui>
#include <iostream>
#include "/media/Khue/Dev. Software/C++libraries/boost/boost/utility.hpp"

typedef std::list<Point> PointList;

class Scene : public QGraphicsScene
{
	Q_OBJECT
public:	
	Scene(QObject *parent = 0);
	~Scene();
private:
	//QGraphicsScene* scene;
	int numOfPoints;
	Point first;
    Point last;
	PointList points;
	bool mouseInteraction;
private slots:
	void mouseMoveEvent(QGraphicsSceneMouseEvent * e);
	void mousePressEvent(QGraphicsSceneMouseEvent * e);
	void mouseReleaseEvent(QGraphicsSceneMouseEvent* e); 
public:
	void setNumOfPoints(int n){ numOfPoints = n; }
    void enableMouseIteration(bool b){ mouseInteraction = b ;}
    void draw(PointList points, QPen pen);
	void drawLine(PointList points);
    void drawLineWithSP(PointList points);
    //void drawLine(PointList points, bool displayIndex);
    void drawLine(PointList points,QColor color);
    inline PointList getPoints(){ return points; }
    inline void setPoints(PointList l){ points = l; }
};

#endif // SCENE_H
