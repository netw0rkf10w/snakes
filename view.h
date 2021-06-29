#ifndef VIEW_H
#define VIEW_H

#include "point.h"
#include <list>
#include <QtGui>
#include <QGraphicsView>
#include <iostream>


class View : public QGraphicsView
{
	//Q_OBJECT

public:	
	View(QWidget *parent = 0);
	~View();
private:
	//QGraphicsScene* scene;
    QString imageFile;
    QImage image;

protected:
	//Holds the current centerpoint for the view, used for panning and zooming
    QPointF CurrentCenterPoint;
 
    //From panning the view
    QPoint LastPanPoint;
 
    //Set the current centerpoint in the
    void SetCenter(const QPointF& centerPoint);
    QPointF GetCenter() { return CurrentCenterPoint; }
 
    //Take over the interaction
    /*virtual void mousePressEvent(QMouseEvent* event);
    virtual void mouseReleaseEvent(QMouseEvent* event);
    virtual void mouseMoveEvent(QMouseEvent* event);*/
    virtual void wheelEvent(QWheelEvent* event);
    virtual void resizeEvent(QResizeEvent* event);
    virtual void drawBackground ( QPainter * painter, const QRectF & rect );
public:
	//void setNumOfPoints2(int n){numOfPoints=n;}
    void setImageFile(QString f){ imageFile=f; }
    void setImage(QImage im){ image=im; }
};

#endif // VIEW_H
