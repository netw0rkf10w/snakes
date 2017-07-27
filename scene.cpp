#include "scene.h"

const double PI = 3.1415926536; 
PointList createPoints(Point begin, Point end, int num);
PointList shiftPoints(PointList points, Point begin, Point end);

double distance(Point A, Point B) 
{
	int ux = A.getX()-B.getX();
	int uy = A.getY()-B.getY();
	double un = ux*ux+uy*uy;
    return sqrt(un);
}

Scene::Scene(QObject *parent):QGraphicsScene(parent)
{
    numOfPoints=40;
	mouseInteraction=true;
}

Scene::~Scene()
{
	points.clear();
}

void Scene::draw(PointList points, QPen pen)
{
    PointList::iterator i;
	for(i=points.begin();i!=points.end();i++)
	{
        //addLine(i->getX(),i->getY(),i->getX(),i->getY(),pen);
        addRect(i->getX(),i->getY(),1,1,pen);
	}		
}

void Scene::drawLine(PointList points)
{
    QPen pen(QBrush(Qt::red),3);
    PointList::iterator i;
	for(i=points.begin();i!=points.end();i++)
	{
        addRect(i->getX()-2,i->getY()-2,4,4,pen);//addEllipse(i->getX()-2,i->getY()-2,4,4,pen);
        if(i==boost::prior(points.end()))
        {
            this->addLine(i->getX(),i->getY(),points.begin()->getX(),points.begin()->getY(),pen);
        }
        else
        {
            this->addLine(boost::next(i)->getX(),boost::next(i)->getY(),i->getX(),i->getY(),pen);
        }
	}
}

void Scene::drawLineWithSP(PointList points)
{
    QPen pen(QBrush(Qt::red),3);
    PointList::iterator i;
    for(i=points.begin();i!=points.end();i++)
    {
        if(i==points.begin())
        {
            addRect(i->getX()-4,i->getY()-4,8,8,pen);  //addEllipse(i->getX()-4,i->getY()-4,8,8,pen);
        }
        else
        {
            addRect(i->getX()-2,i->getY()-2,4,4,pen);//addEllipse(i->getX()-2,i->getY()-2,4,4,pen);
            this->addLine(boost::prior(i)->getX(),boost::prior(i)->getY(),i->getX(),i->getY(),pen);
        }

    }
}

//void Scene::drawLine(PointList points, bool displayIndex)
//{
//    QPen pen1(QBrush(Qt::red),3);
//    QPen pen2(QBrush(Qt::blue),3);
//    QPen pen;
//    PointList::iterator i;
//    for(i=points.begin();i!=points.end();i++)
//    {
//        if(i->isOriginal()) pen=pen1;
//        else pen=pen2;
//        if(i==points.begin())
//        {
//            addRect(i->getX()-4,i->getY()-4,8,8,pen);  //addEllipse(i->getX()-4,i->getY()-4,8,8,pen);
//        }
//        else
//        {
//            addRect(i->getX()-2,i->getY()-2,4,4,pen);//addEllipse(i->getX()-2,i->getY()-2,4,4,pen);
//            this->addLine(boost::prior(i)->getX(),boost::prior(i)->getY(),i->getX(),i->getY(),pen);
//        }
//        if(displayIndex)
//        {
//            if(i->getIndex() > 0)
//            {
//                QGraphicsTextItem *text=addText(QString::number(i->getIndex()));
//                text->setPos(i->getX(),i->getY());
//            }
//        }
//
//    }
//}


void Scene::drawLine(PointList points,QColor color)
{
    QPen pen(QBrush(color),3);
    PointList::iterator i;
    for(i=points.begin();i!=points.end();i++)
    {
        if(i==points.begin()) addEllipse(i->getX()-4,i->getY()-4,8,8,pen);
        else
        {
            addEllipse(i->getX()-2,i->getY()-2,4,4,pen);
            this->addLine(boost::prior(i)->getX(),boost::prior(i)->getY(),i->getX(),i->getY(),pen);
        }

    }
}

void Scene::mouseMoveEvent(QGraphicsSceneMouseEvent* e)
{
	
	if(mouseInteraction)
	{
		last = Point(e->scenePos().x(),e->scenePos().y());
		if(e->buttons() & Qt::LeftButton)
		{
			//this->addLine(this->first.getX(),this->first.getY(),e->scenePos().x(),e->scenePos().y(),pen);
			points.clear();
			points=createPoints(this->first,last,this->numOfPoints);	
		}	
		if(e->buttons() & Qt::RightButton)
		{
			
			PointList s=shiftPoints(points,first,last);
			if((int)s.size()>0)
			{			
				points.clear();
				points=s;	
			}
			
			first=last;
		}
		this->clear();
        QPen pen(QBrush(Qt::red),3);
        draw(points,pen);
	}
	
}

void Scene::mousePressEvent(QGraphicsSceneMouseEvent* e)
{
	first.setX(e->scenePos().x());
	first.setY(e->scenePos().y());
	//if(e->button()==Qt::LeftButton) std::cout<<first.getX()<<"\t"<<first.getX()<<std::endl;
	//if(e->button()==Qt::LeftButton) addEllipse(e->scenePos().x()-5,e->scenePos().y()-5,10,10,pen);
}

void Scene::mouseReleaseEvent(QGraphicsSceneMouseEvent* e)
{
	last.setX(e->pos().x());
	last.setY(e->pos().y());
	//std::cout<<last.x()<<"\t"<<last.y()<<std::endl;
}

PointList createPoints(Point begin, Point end, int num)
{
	double angle=2*PI/num;
	int radius=distance(begin,end);	
	Point p; 
	Point center((end.getX()+begin.getX())/2,(end.getY()+begin.getY())/2);
	PointList s;
	for(int i=0;i<num;i++)
	{
		p.setLocation(radius*cos(angle*i)+center.getX(),radius*sin(angle*i)+center.getY());
		s.push_back(p);
	}
	return s;
}

PointList shiftPoints(PointList points, Point begin, Point end)
{
	PointList s;
	PointList::iterator i;
	int dx=end.getX()-begin.getX();
	int dy=end.getY()-begin.getY();
	Point p; 
	for(i=points.begin();i!=points.end();i++)
	{
		p.setLocation(i->getX()+dx,i->getY()+dy);
		s.push_back(p);
	}
	return s;
}
