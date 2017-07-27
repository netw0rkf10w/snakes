#include "snake.h"
#include <algorithm>
#include <fstream>
#include <float.h>
#include "stuffs.h"

#define	NBHD_RADIUS 1

const double PI = 3.1415926536; 



void normalize(double** array3x3); 
void normalize2(double** array3x3);
void normalize3(double** array3x3);

Snake Snake::replace(Point oldpoint,Point newpoint)
//Replace the i-th point of the snake by the point p
{
	Snake news=*this;
	PointList list=news.getSnake();
	std::replace(list.begin(),list.end(),oldpoint,newpoint);
	news.setSnake(list);
	list.clear();
	return news;
}


//Default constructor
Snake::Snake()
{
	//SnakeSize=0;
	//SnakePixels=std::list<Pixel>::;
	//SnakeImage=NULL;
	//dynamicAllocation();
}

//Destructor
Snake::~Snake()
{
	this->snake.clear();
	this->shape.clear();
    snakeU.clear();
    snakeL.clear();
	//freeAllocation();
	/*int i;
	for (i = 0; i < width; ++i)
	{
		delete[] gradient[i];
		delete[] flow[i];
	}
	delete[] gradient;
	delete[] flow;*/
}



Snake& Snake::operator=(const Snake &s)
{
	this->e_balloon=s.e_balloon;
	this->e_curvature=s.e_curvature;
	this->e_flow=s.e_flow;
	this->e_image=s.e_image;
	this->e_prior=s.e_prior;
	this->e_uniformity=s.e_uniformity;
	this->flow=s.flow;
	this->gradient=s.gradient;
	this->height=s.height;
	this->shape=s.shape;
	this->snake=s.snake;
	this->snakelength=s.snakelength;
	this->useDTW=s.useDTW;
	this->useSPI=s.useSPI;
	this->width=s.width;
	return *this;
}


void Snake::dynamicAllocation()
{
	e_uniformity = new double*[3];
	e_curvature = new double*[3];
	e_image = new double*[3];
	e_flow = new double*[3];
	e_balloon = new double*[3];
	e_prior = new double*[3];
	for(int i=0;i<3;i++) 
	{
		e_uniformity[i]=new double[3];
		e_curvature[i]=new double[3];
		e_image[i]=new double[3];
		e_flow[i]=new double[3];
		e_balloon[i]=new double[3];
		e_prior[i]=new double[3];
	}	
}

void Snake::freeAllocation()
{
	for(int i=0;i<3;i++) 
	{
		delete[] e_uniformity[i];
		delete[] e_curvature[i];
		delete[] e_image[i];
		delete[] e_flow[i];
		delete[] e_balloon[i];
		delete[] e_prior[i];
	}
	delete[] e_uniformity;	
	delete[] e_curvature;
	delete[] e_image;
	delete[] e_flow;
	delete[] e_balloon;
	delete[] e_prior;
}


Snake::Snake(int width, int height, double** gradient, double** flow, PointList points, PointList shape)
{
	this->snake = points;
	this->shape=shape;
	this->gradient = gradient;
	this->flow=flow;
	this->width = width;
	this->height = height;

	// VERY IMPORTANT !!! Dynamic allocation
	dynamicAllocation();	
}



void Snake::initGradientFlow(const char* fileName)
{
    /* Read the image */
    typedef itk::ImageFileReader<ImageType>		ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fileName);
    reader->Update();

    /* Create the gradient image */
//    typedef itk::GradientMagnitudeImageFilter< ImageType, ImageType > GradientFilter;
//    GradientFilter::Pointer fGrad = GradientFilter::New();
//    fGrad->SetInput(reader->GetOutput());
//    fGrad->Update();

    typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< ImageType, ImageType > GradientFilter;
    GradientFilter::Pointer fGrad = GradientFilter::New();
    fGrad->SetInput(reader->GetOutput());
    fGrad->SetSigma(1.0);
    fGrad->Update();
    ImageType::Pointer imageGrad = fGrad->GetOutput();

    /* Rescale to the range 0-255 */
    typedef itk::RescaleIntensityImageFilter<ImageType,ImageType>    RescaleFilterType2;
    RescaleFilterType2::Pointer gradient0255 = RescaleFilterType2::New();
    gradient0255->SetOutputMinimum(   0 );
    gradient0255->SetOutputMaximum( 255 );
    gradient0255->SetInput(fGrad->GetOutput());
    gradient0255->Update();

    /* Create the vector flow image*/
    typedef  itk::ApproximateSignedDistanceMapImageFilter< ImageType, ImageType  > ApproximateSignedDistanceMapImageFilterType;
    ApproximateSignedDistanceMapImageFilterType::Pointer fFlow = ApproximateSignedDistanceMapImageFilterType::New();
    fFlow->SetInput(gradient0255->GetOutput());
    fFlow->SetInsideValue(255);
    fFlow->SetOutsideValue(0);
    fFlow->Update();
    ImageType::Pointer imageFlow = fFlow->GetOutput();

    //Get the width and the height of the image
	width=imageGrad->GetLargestPossibleRegion().GetSize()[0]; 
	height=imageGrad->GetLargestPossibleRegion().GetSize()[1];

    gradient = newArray(width,height);
    flow = newArray(width,height);
	//
	ImageType::IndexType index;
	for(int i=0;i<width;i++)
	{
		for(int j=0;j<height;j++)
        {
			index[0]=i;
			index[1]=j;
			gradient[i][j]=imageGrad->GetPixel(index);
            flow[i][j]=255.0-imageFlow->GetPixel(index);
		}
	}


    //Save the flow image
//    typedef itk::ImageFileWriter<UnsignedCharImageType>  WriterType;
//    typedef itk::RescaleIntensityImageFilter<FloatImageType,UnsignedCharImageType>    FilterType;

//    FilterType::Pointer filter = FilterType::New();
//    filter->SetOutputMinimum(   0 );
//    filter->SetOutputMaximum( 255 );

//    WriterType::Pointer writer = WriterType::New();
//    filter->SetInput(fFlow->GetOutput());
//    filter->Update();
//    writer->SetInput(filter->GetOutput());
//     writer->SetFileName("1.png");
//    try
//    {
//        writer->Update();
//    }
//    catch( itk::ExceptionObject & err )
//    {
//        std::cerr << "ExceptionObject caught !" << std::endl;
//        std::cerr << err << std::endl;
//        return;
//    }


}



void Snake::initGradientFlow(const char* fileName, double sigma, double threshold)
{
    /* Read the image */
    typedef itk::ImageFileReader<ImageType>		ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fileName);
    reader->Update();

    typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< ImageType, ImageType > GradientFilter;
    GradientFilter::Pointer fGrad = GradientFilter::New();
    fGrad->SetInput(reader->GetOutput());
    fGrad->SetSigma(sigma);
    fGrad->Update();
    //save this as the gradient image for snake processing later
    ImageType::Pointer imageGrad = fGrad->GetOutput();


    /* Rescale to the range 0-255 */
    typedef itk::RescaleIntensityImageFilter<ImageType,ImageType>    RescaleFilterType2;
    RescaleFilterType2::Pointer gradient0255 = RescaleFilterType2::New();
    gradient0255->SetOutputMinimum(   0 );
    gradient0255->SetOutputMaximum( 255 );
    gradient0255->SetInput(fGrad->GetOutput());
    gradient0255->Update();


    /* Add a threshold filter to adjust the edges */
    typedef itk::ThresholdImageFilter <ImageType>        ThresholdImageFilterType;
    ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();
    thresholdFilter->SetInput(gradient0255->GetOutput());
    thresholdFilter->ThresholdOutside(0, threshold);
    thresholdFilter->SetOutsideValue(255);

    /* Create the vector flow image from the thresholded image */
    typedef  itk::ApproximateSignedDistanceMapImageFilter< ImageType, ImageType  > ApproximateSignedDistanceMapImageFilterType;
    ApproximateSignedDistanceMapImageFilterType::Pointer fFlow = ApproximateSignedDistanceMapImageFilterType::New();
    fFlow->SetInput(thresholdFilter->GetOutput());
    fFlow->SetInsideValue(0);
    fFlow->SetOutsideValue(255);
    fFlow->Update();
    ImageType::Pointer imageFlow = fFlow->GetOutput();

    //Get the width and the height of the image
    width=imageGrad->GetLargestPossibleRegion().GetSize()[0];
    height=imageGrad->GetLargestPossibleRegion().GetSize()[1];

    gradient = newArray(width,height);
    flow = newArray(width,height);
    //
    ImageType::IndexType index;
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {
            index[0]=i;
            index[1]=j;
            gradient[i][j]=imageGrad->GetPixel(index);
            flow[i][j]=imageFlow->GetPixel(index);
        }
    }
}



void Snake::setGradientFlow(QImage imageGrad, QImage imageFlow)
{
    //Get the width and the height of the image
    int width=imageGrad.width();
    int height=imageGrad.height();
    gradient = newArray(width,height);
    flow = newArray(width,height);

    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {
            gradient[i][j]=imageGrad.pixelIndex(i,j);
            flow[i][j]=imageFlow.pixelIndex(i,j);
        }
    }
}


Snake::Snake(const char* fileName)
{
    initGradientFlow(fileName);
	int dimensionMin=std::min(width,height);
	int radius=dimensionMin/2-5;	
	Point p; 
	Point center(width/2,height/2);
	PointList s;

	double angle=0.314159/2; //18.0 degrees;
	for(int i=0;i<40;i++)
	{
		p.setLocation(radius*cos(angle*i)+center.getX(),radius*sin(angle*i)+center.getY());
		s.push_back(p);
	}
	this->snake=s;
    this->shape=s;
	
	// VERY IMPORTANT !!! Dynamic allocation
	dynamicAllocation();	
}

Snake::Snake(const char* fileName, double sigma, double threshold)
{
    initGradientFlow(fileName,sigma,threshold);
    int dimensionMin=std::min(width,height);
    int radius=dimensionMin/2-5;
    Point p;
    Point center(width/2,height/2);
    PointList s;

    double angle=0.314159/2; //18.0 degrees;
    for(int i=0;i<40;i++)
    {
        p.setLocation(radius*cos(angle*i)+center.getX(),radius*sin(angle*i)+center.getY());
        s.push_back(p);
    }
    this->snake=s;
    this->shape=s;

    // VERY IMPORTANT !!! Dynamic allocation
    dynamicAllocation();
}

Snake::Snake(const char* fileName, PointList shape)
{
    initGradientFlow(fileName);
	int dimensionMin=std::min(width,height);
	int radius=dimensionMin/2-5;	
	Point p; 
	Point center(width/2,height/2);
	PointList s;
	
	double angle=2*PI/(int)shape.size();
	for(int i=0;i<(int)shape.size();i++)
	{
		p.setLocation(radius*cos(angle*i)+center.getX(),radius*sin(angle*i)+center.getY());
		s.push_back(p);
	}
	this->snake=s;
	this->shape=shape;
	// VERY IMPORTANT !!! Dynamic allocation
	dynamicAllocation();	
}

Snake::Snake(const char* fileName, PointList shape, double sigma, double threshold)
{
    initGradientFlow(fileName,sigma,threshold);
    int dimensionMin=std::min(width,height);
    int radius=dimensionMin/2-5;
    Point p;
    Point center(width/2,height/2);
    PointList s;

    double angle=2*PI/(int)shape.size();
    for(int i=0;i<(int)shape.size();i++)
    {
        p.setLocation(radius*cos(angle*i)+center.getX(),radius*sin(angle*i)+center.getY());
        s.push_back(p);
    }
    this->snake=s;
    this->shape=shape;
    // VERY IMPORTANT !!! Dynamic allocation
    dynamicAllocation();
}

Point Snake::get(int j)
{
	PointList::iterator i=snake.begin();
	std::advance(i,j);
    Point p(i->getX(),i->getY());
	return p;
}

 

bool Snake::step1(double continuityCoeff,double curvatureCoeff,double imageCoeff,double flowCoeff,double balloonCoeff) 
{
	bool changed=false;
	Point p(0,0);
	int x,y;
 
	// compute length of original snake (used by method: f_uniformity)
	this->snakelength = getsnakelength();
	PointList newsnake;
	// for each point of the previous snake
	for(int i=0;i<(int)snake.size();i++) 
	{
		Point prev = get((i+(int)snake.size()-1)%(int)snake.size());
		Point cur  = get(i);
		Point next = get((i+1)%(int)snake.size());

        // compute all energies
		for(int dy=-1;dy<=1;dy++) 
		{
			for(int dx=-1;dx<=1;dx++) 
			{
				p.setLocation(cur.getX()+dx, cur.getY()+dy);
				e_uniformity[1+dx][1+dy] = f_uniformity(prev,next,p);
				e_curvature[1+dx][1+dy]  = f_curvature(prev,p,next);
				e_flow[1+dx][1+dy]       = f_gflow(cur,p);
				e_balloon[1+dx][1+dy]    = f_balloon(prev,p,cur,next);
				e_image[1+dx][1+dy]    = f_image(cur,p);
			}
		}
 
		// normalize energies
		/*normalize(e_uniformity);
		normalize(e_curvature);
		normalize(e_flow);
		normalize(e_balloon);*/
        normalize3(e_uniformity);
        normalize3(e_curvature);
        normalize3(e_flow);
        normalize3(e_image);
        normalize3(e_balloon);
 
		// find the point with the minimum sum of energies
		double emin, e;
        emin=0;
        x=cur.getX();
        y=cur.getY();

		if(continuityCoeff	>0) emin+= continuityCoeff * e_uniformity[1][1]; // internal energy
		if(curvatureCoeff	>0) emin+= curvatureCoeff  * e_curvature[1][1];  // internal energy
		if(imageCoeff		>0) emin+= imageCoeff * e_image[1][1];    // external energy
		if(flowCoeff		>0) emin+= flowCoeff * e_flow[1][1];       // external energy				
		if(balloonCoeff		>0) emin+= balloonCoeff * e_balloon[1][1]; 

		for(int dy=-1;dy<=1;dy++) 
		{
			for(int dx=-1;dx<=1;dx++) 
			{
				e = 0;
				if(continuityCoeff	>0) e+= continuityCoeff * e_uniformity[1+dx][1+dy]; // internal energy
				if(curvatureCoeff	>0) e+= curvatureCoeff  * e_curvature[1+dx][1+dy];  // internal energy
				if(imageCoeff		>0) e+= imageCoeff * e_image[1+dx][1+dy];    // external energy
				if(flowCoeff		>0) e+= flowCoeff * e_flow[1+dx][1+dy];       // external energy				
				if(balloonCoeff		>0) e+= balloonCoeff * e_balloon[1+dx][1+dy]; 
 
				if (e<emin) 
				{ 
					emin=e; 
					x=cur.getX()+dx; 
					y=cur.getY()+dy; 
				}
			}
		}
 
		// boundary check
		if (x<1) x=1;
		if (x>=(this->width-1)) x=this->width-2;
		if (y<1) y=1;
		if (y>=(this->height-1)) y=this->height-2;
 
		// compute the returned value
		
		if (x!=cur.getX() || y!=cur.getY()) 
		{
			changed=true;			
		}
			
 
		// create the point in the new snake
		Point p1(x,y);
		newsnake.push_back(p1);			
	}
 
	// new snake becomes current
	this->snake.clear();
	this->snake=newsnake;
 
	return changed;
}



 
bool Snake::step2(double continuityCoeff,double curvatureCoeff,double imageCoeff,double flowCoeff,double balloonCoeff,double priorCoeff) 
{
	bool changed=false;
	Point p(0,0);
	double emin, e;
	int x,y;
 
	// compute length of original snake (used by method: f_uniformity)
	this->snakelength = getsnakelength();
	
	
	// for each point of the previous snake
	for(int i=0;i<(int)snake.size();i++) 
	{
		Point prev = get((i+(int)snake.size()-1)%(int)snake.size());
		Point cur  = get(i);
		Point next = get((i+1)%(int)snake.size());

		// compute all energies
		for(int dy=-1;dy<=1;dy++) 
		{
			for(int dx=-1;dx<=1;dx++) 
			{
				p.setLocation(cur.getX()+dx, cur.getY()+dy);
				e_uniformity[1+dx][1+dy] = f_uniformity(prev,next,p);
				e_curvature[1+dx][1+dy]  = f_curvature(prev,p,next);
				e_flow[1+dx][1+dy]       = f_gflow(cur,p);
				e_balloon[1+dx][1+dy]    = f_balloon(prev,p,cur,next);
				e_image[1+dx][1+dy]    = f_image(cur,p);
				if(priorCoeff>0) 
				{
					Snake s2=this->replace(cur,p);
					e_prior[1+dx][1+dy]  = s2.f_prior();
				}
			}
		}
 
		// normalize energies
		normalize3(e_uniformity);
        normalize3(e_curvature);
        normalize3(e_flow);
		normalize3(e_image);
        normalize3(e_balloon);
		if(priorCoeff>0) normalize3(e_prior);
 
		// find the point with the minimum sum of energies
		
		emin=0;
		x=cur.getX(); 
		y=cur.getY(); 
		if(continuityCoeff	!=0) emin+= continuityCoeff * e_uniformity[1][1]; // internal energy
		if(curvatureCoeff	!=0) emin+= curvatureCoeff  * e_curvature[1][1];  // internal energy
		if(imageCoeff		!=0) emin+= imageCoeff * e_image[1][1];    // external energy
		if(flowCoeff		!=0) emin+= flowCoeff * e_flow[1][1];       // external energy				
		if(balloonCoeff		!=0) emin+= balloonCoeff * e_balloon[1][1]; 
		if(priorCoeff		!=0)	emin+= priorCoeff * e_prior[1][1]; 
		for(int dy=-1;dy<=1;dy++) 
		{
			for(int dx=-1;dx<=1;dx++) 
			{
				e = 0;
				if(continuityCoeff	!=0) e+= continuityCoeff * e_uniformity[1+dx][1+dy]; // internal energy
				if(curvatureCoeff	!=0) e+= curvatureCoeff  * e_curvature[1+dx][1+dy];  // internal energy
				if(imageCoeff		!=0) e+= imageCoeff * e_image[1+dx][1+dy];    // external energy
				if(flowCoeff		!=0) e+= flowCoeff * e_flow[1+dx][1+dy];       // external energy				
				if(balloonCoeff		!=0) e+= balloonCoeff * e_balloon[1+dx][1+dy]; 
				if(priorCoeff		!=0)	e+= priorCoeff * e_prior[1+dx][1+dy]; 
				if (e<emin) 
				{ 
					emin=e; 
					x=cur.getX()+dx; 
					y=cur.getY()+dy; 
				}
			}
		}
	
		// boundary check
		if (x<1) x=1;
		if (x>=(this->width-1)) x=this->width-2;
		if (y<1) y=1;
		if (y>=(this->height-1)) y=this->height-2;

		if (x!=cur.getX() || y!=cur.getY()) 
		{
			changed=true;
			// create the point in the new snake
			PointList newsnake;
            Point p1(x,y);
            newsnake=this->replace(cur,p1).getSnake();
			this->snake.clear();
			this->snake=newsnake;
			newsnake.clear();
		}
		// compute the returned value
		/*Point p1(x,y);
		newsnake.push_back(p1);	*/
	}
	//this->snake=newsnake;
	
	return changed;
}



// normalize energy matrix
void normalize(double** array3x3) 
{
	double sum=0;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			sum+=std::abs(array3x3[i][j]);
 
	if (sum==0) return;
 
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			array3x3[i][j]/=sum;
}

void normalize2(double** array3x3) 
{
	double max=0;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			if(std::abs(array3x3[i][j])>max)
				max=std::abs(array3x3[i][j]);
		}
	}		
 
	if (max==0) return;
 
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			array3x3[i][j]/=max;
}

void normalize3(double** array3x3) 
{
	double max=0;
	double min=DBL_MAX;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			if(std::abs(array3x3[i][j])>max)
				max=std::abs(array3x3[i][j]);
			if(std::abs(array3x3[i][j])<min)
				min=std::abs(array3x3[i][j]);
		}
	}		
 
	if (max==min) return;
 
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			array3x3[i][j]=(array3x3[i][j]-min)/(max-min);
}

 
double Snake::getsnakelength() 
{
	// total length of snake
	double length=0;
	for(int i=0;i<(int)snake.size();i++) 
	{
		Point cur   = get(i);
		Point next  = get((i+1)%(int)snake.size());
		length+=distance2D(cur, next);
	}
	return length;
}
 


double Snake::f_uniformity(Point prev, Point next, Point p) 
{
 
	// length of previous segment
	double un = distance2D(prev, p);
 
	// mesure of uniformity
	double avg = snakelength/(int)snake.size();
	double dun = std::abs(un-avg);
 
	// elasticity energy
	return dun*dun;
}
 
double Snake::f_curvature(Point prev, Point p, Point next) 
{
	int ux = p.getX()-prev.getX();
	int uy = p.getY()-prev.getY();
	double un = std::sqrt((double)ux*ux+uy*uy);
 
	int vx = p.getX()-next.getX();
	int vy = p.getY()-next.getY();
	double vn = std::sqrt((double)vx*vx+vy*vy);
 
	if (un==0 || vn==0) return 0;
 
	double cx = (vx+ux)/(un*vn);
	double cy = (vy+uy)/(un*vn);
 
	// curvature energy
	double cn = cx*cx+cy*cy;
	return cn;
}
 
double Snake::f_gflow(Point cur, Point p) 
{
	// gradient flow
//	int dcur = this->flow[cur.getX()][cur.getY()];
//	int dp   = this->flow[p.getX()][p.getY()];
//	double d = dp-dcur;
//	return d;
    double dp   = -this->flow[p.getX()][p.getY()];
    return dp;
}

double Snake::f_image(Point cur, Point p) 
{
//	double d = distance2D(cur, p);
//	double g = this->gradient[cur.getX()][cur.getY()];
//	double e = g*d;
//	return e;
    double g = this->gradient[p.getX()][p.getY()];
    return g;
}

double Snake::f_balloon(Point prev, Point p, Point cur, Point next)
{ 
	// 
	Point t=next-prev;
    double tx=(double)t.getX()/t.norm();
    double ty=(double)t.getY()/t.norm();
	Point s=p-cur;
    double u=(double)s.getX()+(2*NBHD_RADIUS+1)*ty;
    double v=(double)s.getY()-(2*NBHD_RADIUS+1)*tx;
	return u*u+v*v;
}

double Snake::f_prior() 
{
    complex* sh=this->DFTShape();
    complex* sn=this->DFTSnake();
    complex* Sh=iDFT(sh,(int)this->shape.size());
    complex* Sn=iDFT(sn,(int)this->snake.size());
    double d;
    switch(distanceType)
    {
    case EUCLIDEAN:
        d=0;
        for(int i=0;i<(int)snake.size();i++)
        {
            d=d+std::norm(sh[i]-sn[i]);
        }
        //d=this->EuclideanModified();
        break;
    case DTW:
        //d=dtwDistance(Sh,Sn,(int)this->shape.size(),this->windowLength);
        d=dtwDistance(Sh,Sn,(int)this->shape.size());
        break;
    case PRATT:
        d=1-Pratt(Sn,(int)this->snake.size(),Sh,(int)this->shape.size());
        break;
    case GDM:
        d = interArea(Sn,(int)this->snake.size(),Sh,(int)this->shape.size());
        break;
    case LARSSON:
        d=1-maxR12(sn,(int)this->snake.size(),sh,(int)this->shape.size());
        break;
    }
    delete[] sn;
    delete[] Sn;
    return d;
}

double Snake::EuclideanModified()
{
    complex* sh=this->DFTShape();
    complex* Sh=iDFT(sh,(int)this->shape.size());
    //complex* sh=PointsToDFT(this->shape,this->useSPI,this->useRotInvar);
    complex* sn;
    complex* Sn;
    int N=snake.size();
    double min=DBL_MAX;
    double d;
    for(int i=0;i<N;i++)
    {
        //sn=this->DFTSnake();
        PointList l=this->shiftSnake(i);
        sn=PointsToDFT(l,this->useSPI,this->useRotInvar);
        Sn=iDFT(sn,(int)this->snake.size());
        d=0;
        for(int i=0;i<(int)snake.size();i++)
        {
            d=d+std::norm(sh[i]-sn[i]);
        }
        //Get the index where the distance is minimal
        if(d < min)
        {
            min=d;
        }
        l.clear();
        delete[] sn;
        delete[] Sn;
    }
    delete[] sh;
    delete[] Sh;
    return min;
}

//double interArea(complex* S, int sizeS, complex* T, int sizeT)
//{
//    typedef CGAL::Simple_cartesian<double> K;
//    typedef K::Point_2 PointCGAL;
//    typedef CGAL::Polygon_2<K> Polygon_2;
//    typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;

//    PointCGAL *points=new PointCGAL[sizeS];
//    PointCGAL *points2=new PointCGAL[sizeT];
//    for(int i=0;i<sizeS;i++)
//    {
//        points[i]=PointCGAL(100*S[i].real(),100*S[i].imag());
//        //points[i]=PointCGAL(S[i].real()/500,S[i].imag()/500);
//    }
//    for(int i=0;i<sizeT;i++)
//    {
//        points2[i]=PointCGAL(100*T[i].real(),100*T[i].imag());
//        //points2[i]=PointCGAL(T[i].real()/500,T[i].imag()/500);
//    }
//    Polygon_2 poly1(points, points+sizeS);
//    Polygon_2 poly2(points2, points2+sizeT);
//    //CGAL::General_polygon_with_holes_2<K> poly3;
//    std::list<Polygon_with_holes_2> polyI;

//    CGAL::intersection(poly1, poly2, std::back_inserter(polyI));
//    //CGAL::difference(poly1, poly2, std::back_inserter(polyI));
//    //CGAL::symmetric_difference(poly1, poly2, std::back_inserter(polyI));
//    double totalArea = 0;
//    typedef std::list<Polygon_with_holes_2>::iterator LIT;
//    for(LIT lit = polyI.begin(); lit!=polyI.end(); lit++)
//    {
//        totalArea+=lit->outer_boundary().area();
//    }
//    double a1=poly1.area();
//    double a2=poly2.area();
//    totalArea = a1+a2-2*totalArea;
//    delete[] points;
//    delete[] points2;
//    polyI.clear();
//    return totalArea;
//}

double interArea(complex* S, int sizeS, complex* T, int sizeT)
{
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_2 PointCGAL;
    typedef CGAL::Polygon_2<K> Polygon_2;
    typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;

    PointCGAL *points=new PointCGAL[sizeS];
    PointCGAL *points2=new PointCGAL[sizeT];
    for(int i=0;i<sizeS;i++)
    {
        points[i]=PointCGAL(100*S[i].real(),100*S[i].imag());
        //points[i]=PointCGAL(S[i].real()/500,S[i].imag()/500);
    }
    for(int i=0;i<sizeT;i++)
    {
        points2[i]=PointCGAL(100*T[i].real(),100*T[i].imag());
        //points2[i]=PointCGAL(T[i].real()/500,T[i].imag()/500);
    }
    Polygon_2 poly1(points, points+sizeS);
    Polygon_2 poly2(points2, points2+sizeT);
    //CGAL::General_polygon_with_holes_2<K> poly3;
    std::list<Polygon_with_holes_2> polyI;

    CGAL::intersection(poly1, poly2, std::back_inserter(polyI));
    //CGAL::difference(poly1, poly2, std::back_inserter(polyI));
    //CGAL::symmetric_difference(poly1, poly2, std::back_inserter(polyI));
    double totalArea = 0;
    typedef std::list<Polygon_with_holes_2>::iterator LIT;
    for(LIT lit = polyI.begin(); lit!=polyI.end(); lit++)
    {
        totalArea+=lit->outer_boundary().area();
    }
    double a1=poly1.area();
    double a2=poly2.area();
    totalArea = a1+a2-2*totalArea;
    delete[] points;
    delete[] points2;
    polyI.clear();
    return totalArea;
}


double maxR12(complex* S, int sizeS, complex* T, int sizeT)
//compute iDFT(/S*T) and return the max coefficient
{
    double r;
    complex* P = new complex[sizeS];
    for(int i=0;i<sizeS;i++)
    {
        P[i]=conj(S[i])*T[i];
    }
    r=std::norm(max(iDFT(P,sizeS),0,sizeS-1));
    delete[] P;
    return r;
}

double distance(int i, complex* S, complex* T, int sizeT)
//Distance between the i-th pixel of S and its closest pixel in T
{
    double min=std::norm(S[i]-T[0]);
    for(int j=0;j<sizeT;j++)
    {
        if(min>std::norm(S[i]-T[j])) min=std::norm(S[i]-T[j]);
    }
    return min;
}

complex* fillPoints(complex* c, int n, int N)
//c is the coordinates of the contour, n is the size and N-1 is the number of points added between two points of the contour (N segments)
{
    double x1,y1,x2,y2,stepx,stepy;
    complex* v=new complex[n*N];
    int j=0;
    for(int i=0;i<n;i++)
    {
        x1=c[i].real();
        y1=c[i].imag();
        if(i<n-1)
        {
            x2=c[i+1].real();
            y2=c[i+1].imag();
        }
        else
        {
            x2=c[0].real();
            y2=c[0].imag();
        }
        v[j]=c[i];
        j++;
        stepx=(x2-x1)/N;
        stepy=(y2-y1)/N;
        for(int i=1;i<N;i++)
        {
            v[j]=complex(x1+i*stepx,y1+i*stepy);
            j++;
        }
    }
    return v;
}

double Pratt(complex* S, int sizeS, complex* T, int sizeT)
{
    double sum=0;
//    int N=3;
//    complex* nS=fillPoints(S,sizeS,N);
//    complex* nT=fillPoints(T,sizeT,N);
//    for(int i=0;i<sizeS*N;i++)
//    {
//        sum=sum+(double)1/(1+distance(i,nS,nT,sizeT*N)*distance(i,nS,nT,sizeT*N));
//    }
    for(int i=0;i<sizeS;i++)
    {
        sum=sum+(double)1/(1+distance(i,S,T,sizeT)*distance(i,S,T,sizeT));
    }
    sum=sum/std::max(sizeS,sizeT);
//    delete[] nS;
//    delete[] nT;
    return sum;
}

double PolygonArea(complex* p,int N)
{
    int i,j;
    double area = 0;
    for (i=0;i<N;i++)
    {
      j = (i + 1) % N;
      area += p[i].real() * p[j].imag();
      area -= p[i].imag() * p[j].real();
    }

    area /= 2;
    return(area < 0 ? -area : area);
}

double PartitionDistance(complex* S, int sizeS, complex* T, int sizeT)
{
    double d=PolygonArea(S,sizeS)-PolygonArea(T,sizeT);
    return(d < 0 ? -d : d);
}


complex* getFourierCoeffiecients(PointList snake)
{
	int N=(int)snake.size();
	complex* v=new complex[N];
	
	PointList::iterator j;
	int i=0;
	for(j=snake.begin();j!=snake.end();j++)
	{
		v[i]=complex (j->getX(),j->getY());
		i++;
	}

	// Get the DFT
	complex* C=DFT(v,N);

	// Normalize the Fourier coefficients
	double norm;
	double phase;
	double teta1=std::arg(C[N/2-1]);//std::cout<<"teta1="<<teta1<<std::endl;
	double teta2=std::arg(C[1+N/2]);//std::cout<<"teta2="<<teta2<<std::endl;std::cout<<std::endl;
	
	complex* C2=new complex[N];
	//complex x;
	C2[N/2]=0;
	for(int i=0;i<N;i++)
	{
		if(i!= N/2)
		{
			norm=std::abs(C[i])/std::abs(C[1+N/2]);
			phase=std::arg(C[i])-(teta2+teta1)/2+((double)i-(double)N/2)*(teta1-teta2);
			//phase=std::arg(C[i])-teta2;
			C2[i]=std::polar(norm,phase);
		}
	}
	delete[] C;
	delete[] v;
	return C2;	
}

int Snake::getOptimalStartingPoint()
{
    complex* sh=this->DFTShape();
    complex* Sh=iDFT(sh,(int)this->shape.size());
    //complex* sh=PointsToDFT(this->shape,this->useSPI,this->useRotInvar);
	complex* sn;
    complex* Sn;
    int N=snake.size();
    double min=DBL_MAX;
    double d;
	int startingIndex=0;
	for(int i=0;i<N;i++)
	{
		//sn=this->DFTSnake();
		PointList l=this->shiftSnake(i);
        sn=PointsToDFT(l,this->useSPI,this->useRotInvar);
        Sn=iDFT(sn,(int)this->snake.size());
        switch(distanceTypeSP)
        {
        case EUCLIDEAN:
            d=0;
            for(int i=0;i<(int)snake.size();i++)
            {
                d=d+std::norm(sh[i]-sn[i]);
            }
            break;
        case DTW:
            d=dtwDistance(Sh,Sn,(int)this->shape.size(),this->windowLength);
            break;
        case PRATT:
            d=1-Pratt(Sn,(int)this->snake.size(),Sh,(int)this->shape.size());
            break;
        case GDM:
            d = interArea(Sn,(int)this->snake.size(),Sh,(int)this->shape.size());
            break;
        case LARSSON:
            d=1-maxR12(sn,(int)this->snake.size(),sh,(int)this->shape.size());
            break;
        }

        //Get the index where the distance is minimal
        if(d < min)
		{
            min=d;
			startingIndex=i;
		} 
		l.clear();
		delete[] sn;
        delete[] Sn;
	} 
	delete[] sh;
    delete[] Sh;
	return startingIndex;
}

void Snake::reSample(int n)
//N is the new size of the snake. n < snake.size()
{
    int N=snake.size();
    if(N>n)
    {
        int q=N/n;
        int r=N%n;
        PointList list;
        for(int i=1;i<=n;i++)
        {
            if(i<n-r+1) list.push_back(this->get(i*q-q));
            else list.push_back(this->get(i*q+i+r-n-q));
        }
        snake.clear();
        snake=list;
        int a=snake.size();

        //list.clear();
    }

}


double Snake::getOptimalRotation(double angleStep)
//Return the angle by which the initial snake has to be rotated by
{
    complex* sh=this->DFTShape();
    complex* Sh=iDFT(sh,(int)this->shape.size());
    complex* sn;
    complex* Sn;
    int N=snake.size();
    int n=shape.size();
    double min=DBL_MAX;
    double d;
    double rad=angleStep*PI/180; //angle step in radian
    double bestAngle=0;
    int it=360/angleStep;
    for(int i=0;i<it;i++)
    {
        PointList l=this->rotateSnake(i*angleStep);
        sn=PointsToDFT(l,this->useSPI,this->useRotInvar);
        Sn=iDFT(sn,(int)this->snake.size());
        switch(distanceTypeSP)
        {
        case EUCLIDEAN:
            d=0;
            for(int i=0;i<(int)snake.size();i++)
            {
                d=d+std::norm(sh[i]-sn[i]);
            }
            break;
        case DTW:
            d=dtwDistance(Sh,Sn,(int)this->shape.size(),this->windowLength);
            break;
        case PRATT:
            d=1-Pratt(Sn,(int)this->snake.size(),Sh,(int)this->shape.size());
            break;
        case GDM:
            d = interArea(Sn,(int)this->snake.size(),Sh,(int)this->shape.size());
            break;
        case LARSSON:
            d=1-maxR12(sn,(int)this->snake.size(),sh,(int)this->shape.size());
            break;
        }

        //Get the index where the distance is minimal
        if(d < min)
        {
            min=d;
            bestAngle=i*angleStep;
        }
        l.clear();
        delete[] sn;
        delete[] Sn;
    }
    delete[] sh;
    delete[] Sh;
    return bestAngle;
}

double dtwDistance2(complex *s, complex *t, int N)
{
	int i,j;
	double cost;
	double **DTW=new double*[N];
	for(i=0;i<N;i++)
	{
		DTW[i]=new double[N];
	}
	
	for(i=1;i<N;i++)
	{
		DTW[0][i]= DBL_MAX;
		DTW[i][0]= DBL_MAX;
	}
	DTW[0][0]=0;
	for(i=1;i<N;i++)
	{
		for(j=1;j<N;j++)
		{
			cost=std::norm(s[i]-t[j]);
			DTW[i][j]=cost + std::min(DTW[i-1][j-1],std::min(DTW[i][j-1],DTW[i-1][j]));
		}
	}
	//return std::sqrt(DTW[N-1][N-1]);
	cost = DTW[N-1][N-1];
	for(i=0;i<N;i++)
	{
		delete []DTW[i];
	}
	delete []DTW;
	return cost;
}

double dtwDistance(complex *s, complex *t, int N, int wLenght)
{
	int i,j;
	double cost;
	double **DTW=new double*[N];
	for(i=0;i<N;i++)
	{
		DTW[i]=new double[N];
	}
	
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			DTW[i][j]= DBL_MAX;
		}
	}
	DTW[0][0]=0;
	for(i=1;i<N;i++)
	{
		for(j=std::max(1,i-wLenght);j<std::min(N,i+wLenght);j++)
		{
			cost=std::norm(s[i]-t[j]);
			DTW[i][j]=cost + std::min(DTW[i-1][j-1],std::min(DTW[i][j-1],DTW[i-1][j]));
		}
	}
	//return std::sqrt(DTW[N-1][N-1]);
	cost = DTW[N-1][N-1];
	for(i=0;i<N;i++)
	{
		delete []DTW[i];
	}
	delete []DTW;
	return cost;
}

double dtwDistance(complex *s, complex *t, int N)
{
	int i,j;
	double cost;
	double **DTW=new double*[N];
	double **C=new double*[N];
	for(i=0;i<N;i++)
	{
		DTW[i]=new double[N];
		C[i]=new double[N];
	}
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			C[i][j]=std::norm(s[i]-t[j]);
		}
	}
	DTW[0][0]=0;
	for(i=1;i<N;i++)
	{
		DTW[i][1]= DTW[i-1][1]+C[i][1];
		DTW[1][i]= DTW[1][i-1]+C[1][i];
	}
	for(i=1;i<N;i++)
	{
		for(j=1;j<N;j++)
		{
			DTW[i][j]=C[i][j] + std::min(DTW[i-1][j-1],std::min(DTW[i][j-1],DTW[i-1][j]));
		}
	}
	cost = DTW[N-1][N-1];
	for(i=0;i<N;i++)
	{
		delete[] DTW[i];
		delete[] C[i];
	}
	delete[] DTW;
	delete[] C;
	return cost;
}

complex max(complex *q, int first, int last)
{
	int iMax=0;
	double max=0;
	for(int i=first;i<=last;i++)
	{
		if(max < std::norm(q[i])) 
		{
				max=std::norm(q[i]);
				iMax=i;
		}
	}
	return q[iMax];
}

complex min(complex *q, int first, int last)
{
	int iMin=0;
	double min=DBL_MAX;
	for(int i=first;i<=last;i++)
	{
		if(std::norm(q[i]) < min) 
		{
				min=std::norm(q[i]);
				iMin=i;
		}
	}
	return q[iMin];
}

double LB_Keogh(complex *q, complex *c, int n, int r)
//Keogh, Exact indexing of dynamic time warping
{
	complex* U=new complex[n]; 
	complex* L=new complex[n]; 
	int first,last,i;
	for(i=0;i<n;i++)
	{
		first=i-r;
		if(first < 0) first=0;
		last=i+r;
		if(last > n-1) last=n-1;
		U[i]=max(q,first,last);
		L[i]=min(q,first,last);
	}
	double sum=0;
	for(i=0;i<n;i++)
	{
		if (std::norm(c[i]) > std::norm(U[i])) 
		{
			sum = sum+std::norm(c[i]-U[i]);
		}
		else
		{
			if (std::norm(c[i]) < std::norm(L[i])) 
			{
				sum = sum+std::norm(c[i]-L[i]);
			}
		}
	}
	return std::sqrt(sum);
}

PointList complex2points(complex *z,int N)
{
	PointList list;
    Point p(0,0);
	for(int i=0;i<N;i++)
	{
        p.setLocation((int)(z[i].real()),(int)(z[i].imag()));
		list.push_back(p);
	}
	delete[] z;
	return list;
}

complex* points2complex(PointList points)
{
	int N=(int)points.size();
	complex* v=new complex[N];
	
	PointList::iterator j;
	int i=0;
	for(j=points.begin();j!=points.end();j++)
	{
		v[i]=complex (j->getX(),j->getY());
		i++;
	}
	return v;
}

PointList Snake::getSnakeU()
// upper snake
{
	int n=snake.size();
	complex* U=new complex[n];
	int i;
	int first,last;
	complex* c=points2complex(snake);
	for(i=0;i<n;i++)
	{
		first=i - this->windowLength;
		if(first < 0) first=0;
		last=i + this->windowLength;
		if(last > n-1) last=n-1;
		U[i]=max(c,first,last);
	}
	PointList p=complex2points(U,n);
	delete[] U;
	delete[] c;
	return p;
}

complex* DFT(complex* input, int N)
{
	complex* C=new complex[N];
	complex x;
	complex sum;
	for(int m=0;m<N;m++)	//corresponding to m = -N/2,...,N/2
	{
		sum=0;
		for(int l=0;l<N;l++)//corresponding to l = 0,...,N-1
		{
			x=complex(0.0,-2.0*PI*l*(m-N/2)/N);
			//x=complex(0,-2*PI*l*m/N);
			sum=sum+(input[l])*(std::exp(x));
		}
		//C[m]=((complex)1/(complex)N)*sum;
        C[m]=sum/(complex)N;
	}
	return C;
}

complex* nDFT(complex* input, int N, bool SPInvariant, bool RotInvariant)
{
	complex* C=DFT(input,N);

	// Normalize the Fourier coefficients
	double norm;
	double phase;
	double teta1=std::arg(C[N/2-1]);
	double teta2=std::arg(C[1+N/2]);
	
	complex* C2=new complex[N];
	//complex x;
	C2[N/2]=0;
	for(int i=0;i<N;i++)
	{
		if(i!= N/2)
		{
			norm=std::abs(C[i])/std::abs(C[1+N/2]);
            if(RotInvariant)
            {
                if(SPInvariant)
                {
                    phase=std::arg(C[i])-(teta2+teta1)/2+((double)i-(double)N/2)*(teta1-teta2);
                }
                else
                {
                    phase=std::arg(C[i])-teta2;
                }
            }
            else
            {
                phase=std::arg(C[i]);
            }
            C2[i]=std::polar(norm,phase);
		}
	}	
	delete[] C;
	return C2;	
}


complex* iDFT(complex* input, int N)
{
	complex* z=new complex[N];
	complex sum;
	complex x;
	for(int l=0;l<N;l++) //corresponding to l = 0,...,N-1
	{
		sum=0;
		for(int m=0;m<N;m++) //corresponding to m = -N/2,...,N/2
		{
			x=complex(0,2*PI*l*(m-N/2)/N);
            //x=complex(0,-2*PI*l*m/N);
			sum=sum+(input[m])*(std::exp(x));
		}
        //z[l]=sum/(complex)N;
        z[l]=sum;
	}
	return z;
}


//complex* getFourierCoeffiecients(PointList snake)
//{
//	int N=(int)snake.size();
//	complex* C=new complex[N];
//	PointList::iterator j;
//	int i=0;
//	for(j=snake.begin();j!=snake.end();j++)
//	{
//		C[i]=complex (j->getX(),j->getY());
//		i++;
//	}
//	CFFT::Forward(C,N);
//	// Normalize the coefficients. See Farred's paper
//	double norm;
//	double phase;
//	complex* C2=new complex[N];
//	C2[N/2]=0;
//	for(int i=0;i<N;i++)
//	{
//		if(i!= N/2)
//		{
//			norm=std::abs(C[i])/std::abs(C[1+N/2]);
//			phase=std::arg(C[i])-(1/2)*(std::arg(C[1+N/2])+std::arg(C[-1+N/2])) + ((i-N/2)/2)*(std::arg(C[-1+N/2])-std::arg(C[1+N/2]));
//			C2[i]=std::polar(norm,phase);
//		}
//	}
//	return C2;
//}


PointList iDFT2Points(complex* input,int N)
// Reconstruct the form from the iDFT
{
	PointList list;
	complex* z=iDFT(input,N);
	Point p(0,0);
    double max = 0;
    double norm;
    for(int i=0;i<N;i++)
    {
        norm = std::abs(z[i]);
        if(max < norm) max = norm;
    }

	for(int i=0;i<N;i++)
	{
        p.setLocation(100*z[i].real()/max,100*z[i].imag()/max);
		list.push_back(p);
	}
	delete[] z;
	return list;
}

complex* PointsToDFT(PointList points, bool SPInvariant, bool RotInvariant)
{
	int N=(int)points.size();
	complex* v=new complex[N];
	
	PointList::iterator j;
	int i=0;
	for(j=points.begin();j!=points.end();j++)
	{
		v[i]=complex (j->getX(),j->getY());
		i++;
	}
    complex* C= nDFT(v,N,SPInvariant,RotInvariant);
	delete[] v;
	return C;
}

complex* Snake::DFTSnake()
{
	int N=(int)snake.size();
	complex* v=new complex[N];
	
	PointList::iterator j;
	int i=0;
	for(j=snake.begin();j!=snake.end();j++)
	{
		v[i]=complex (j->getX(),j->getY());
		i++;
	}

    complex* C= nDFT(v,N,this->useSPI,this->useRotInvar);
	delete[] v;
	return C;
}

complex* Snake::DFTShape()
{
	int N=(int)shape.size();
	complex* v=new complex[N];
	
	PointList::iterator j;
	int i=0;
	for(j=shape.begin();j!=shape.end();j++)
	{
		v[i]=complex (j->getX(),j->getY());
		i++;
	}

    complex* C= nDFT(v,N,this->useSPI,this->useRotInvar);
	delete[] v;
	return C;
}


// rebuild the snake using cubic spline interpolation
void Snake::rebuild2(int nmb)
{
    // precompute length(i) = length of the snake from start to point #i
    double* clength = new double[(int)snake.size()+1];
    clength[0]=0;
    for(int i=0;i<(int)snake.size();i++)
    {
        Point cur   = get(i);
        Point next  = get((i+1)%(int)snake.size());
        clength[i+1]=clength[i]+distance2D(cur, next);
    }

    // compute number of points in the new snake
    double total = clength[(int)snake.size()];
    //int nmb = (int)(0.5+total/space);

    // build a new snake
    PointList newsnake;
    for(int i=0,j=0;j<nmb;j++) {
        // current length in the new snake
        double dist = (j*total)/nmb;

        // find corresponding interval of points in the original snake
        while(! (clength[i]<=dist && dist<clength[i+1])) i++;

        // get points (P-1,P,P+1,P+2) in the original snake
        Point prev  = get((i+(int)snake.size()-1)%(int)snake.size());
        Point cur   = get(i);
        Point next  = get((i+1)%(int)snake.size());
        Point next2  = get((i+2)%(int)snake.size());

        // do cubic spline interpolation
        double t =  (dist-clength[i])/(clength[i+1]-clength[i]);
        double t2 = t*t, t3=t2*t;
        double c0 =  1*t3;
        double c1 = -3*t3 +3*t2 +3*t + 1;
        double c2 =  3*t3 -6*t2 + 4;
        double c3 = -1*t3 +3*t2 -3*t + 1;
        double x = prev.getX()*c3 + cur.getX()*c2 + next.getX()* c1 + next2.getX()*c0;
        double y = prev.getY()*c3 + cur.getY()*c2 + next.getY()* c1 + next2.getY()*c0;
        Point newpoint( (int)(0.5+x/6), (int)(0.5+y/6) );

        // add computed point to the new snake
        newsnake.push_back(newpoint);
    }
    this->snake = newsnake;
}


// rebuild the snake using cubic spline interpolation
void Snake::rebuild(int space) 
{
 	// precompute length(i) = length of the snake from start to point #i
	double* clength = new double[(int)snake.size()+1];
	clength[0]=0;
	for(int i=0;i<(int)snake.size();i++) 
	{
		Point cur   = get(i);
		Point next  = get((i+1)%(int)snake.size());
		clength[i+1]=clength[i]+distance2D(cur, next);
	}
 
	// compute number of points in the new snake
    double total = clength[(int)snake.size()];
	int nmb = (int)(0.5+total/space);
 
	// build a new snake
	PointList newsnake;
	for(int i=0,j=0;j<nmb;j++) {
		// current length in the new snake
		double dist = (j*total)/nmb;
 
		// find corresponding interval of points in the original snake
		while(! (clength[i]<=dist && dist<clength[i+1])) i++;
 
		// get points (P-1,P,P+1,P+2) in the original snake
		Point prev  = get((i+(int)snake.size()-1)%(int)snake.size());
		Point cur   = get(i);
		Point next  = get((i+1)%(int)snake.size());
		Point next2  = get((i+2)%(int)snake.size());
 
		// do cubic spline interpolation
		double t =  (dist-clength[i])/(clength[i+1]-clength[i]);
		double t2 = t*t, t3=t2*t;
		double c0 =  1*t3;
		double c1 = -3*t3 +3*t2 +3*t + 1;
		double c2 =  3*t3 -6*t2 + 4;
		double c3 = -1*t3 +3*t2 -3*t + 1;
		double x = prev.getX()*c3 + cur.getX()*c2 + next.getX()* c1 + next2.getX()*c0;
		double y = prev.getY()*c3 + cur.getY()*c2 + next.getY()* c1 + next2.getY()*c0;
		Point newpoint( (int)(0.5+x/6), (int)(0.5+y/6) );
 
		// add computed point to the new snake
		newsnake.push_back(newpoint);
	}
	this->snake = newsnake;
}
 
 
void Snake::removeOverlappingPoints(int minlen) 
{
	// for each point of the snake
	for(int i=0;i<(int)snake.size();i++) 
	{
		Point cur = get(i);
 
		// check the other points (right half)
		for(int di=1+(int)snake.size()/2;di>0;di--) 
		{
			Point end  = get((i+di)%(int)snake.size());
			double dist = distance2D(cur,end);
 
			// if the two points are to close...
			if ( dist>minlen ) continue;
 
			// ... cut the "loop" part og the snake
            for(int k=0;k<di;k++) snake.remove( get((i+1) %(int)snake.size()) );
			break;
		}
	}
}

int Snake::removeOverlappingPoints2(int minlen)
{
    // for each point of the snake
    int n=0;
    for(int i=0;i<(int)snake.size();i++)
    {
        Point cur = get(i);

        // check the other points (right half)
        for(int di=1+(int)snake.size()/2;di>0;di--)
        {
            Point end  = get((i+di)%(int)snake.size());
            double dist = distance2D(cur,end);

            // if the two points are to close...
            if ( dist>minlen ) continue;

            // ... cut the "loop" part og the snake
            for(int k=0;k<di;k++)
            {
                snake.remove( get((i+1) %(int)snake.size()) );
                n++;
            }
            break;
        }
    }
    return n;
}

void Snake::moveAwayPoints(int maxlen) 
// return the number of points removed
{
	// for each point of the snake
	for(int i=0;i<(int)snake.size();i++) 
	{
		Point prev = get((i+(int)snake.size()-1)%(int)snake.size());
		Point cur  = get(i);
		Point next = get((i+1)%(int)snake.size());
		//remove the point that is far away from its neighbors
		if(distance2D(cur,prev)>maxlen && distance2D(cur,next)>maxlen)
		{
			snake.remove(cur);
			PointList::iterator j=snake.begin();
			std::advance(j,i);
			//Point newpoint(
			//snake.insert(j, newpoint); 
		}
			
		
		//add a new point 
		
	}
}
 
void Snake::addMissingPoints(int maxlen) 
{
	PointList::iterator j=snake.begin();
	// for each point of the snake
	for(int i=0;i<(int)snake.size();i++)
	{
		Point prev  = get((i+(int)snake.size()-1)%(int)snake.size());
		Point cur   = get(i);
		Point next  = get((i+1)%(int)snake.size());
		Point next2  = get((i+2)%(int)snake.size());
 
		// if the next point is too far then add a new point
		if ( distance2D(cur,next)>maxlen ) 
		{ 
			// precomputed Uniform cubic B-spline for t=0.5
			double c0=0.125/6.0, c1=2.875/6.0, c2=2.875/6.0, c3=0.125/6.0;
			double x = prev.getX()*c3 + cur.getX()*c2 + next.getX()* c1 + next2.getX()*c0;
			double y = prev.getY()*c3 + cur.getY()*c2 + next.getY()* c1 + next2.getY()*c0;
			Point newpoint( (int)(0.5+x), (int)(0.5+y) );
			snake.insert(boost::next(j), newpoint); 
			i--;
		}
		else
		{
			std::advance(j,1);
		}		
	}
}


void Snake::shift(int index)
{
    PointList newpoints;
    Point p;
    int N=snake.size();
    for(int i=0;i<N;i++)
    {
        p=get((index+i)%N);
        newpoints.push_back(p);
    }
    snake.clear();
    snake=newpoints;
}


void Snake::rotate(double angle)    //rotate all points of the snake by ¨angle¨ degrees
{
    double x,y;
    double rad=angle*PI/180;
    PointList newpoints;
    PointList::iterator i;
    Point p;

    //Get the barycenter of the contour
    Point Center;
    x=0;
    y=0;
    for(i=snake.begin();i!=snake.end();i++)
    {
        x=x+i->getX();
        y=y+i->getY();
    }
    Center.setLocation(x/(int)snake.size(),y/(int)snake.size());

    //Get the snake after the rotation by rad
    for(i=snake.begin();i!=snake.end();i++)
    {
        x=Center.getX() + (i->getX()-Center.getX())*cos(rad)-(i->getY()-Center.getY())*sin(rad);
        y=Center.getX() + (i->getX()-Center.getX())*sin(rad)+(i->getY()-Center.getY())*cos(rad);
        p.setLocation(x,y);
        newpoints.push_back(p);
    }
    snake.clear();
    snake=newpoints;
}


PointList Snake::shiftSnake(int startingIndex)
{
    PointList l;
    Point p;
    int N=snake.size();
    for(int i=0;i<N;i++)
    {
        p=get((startingIndex+i)%N);
        l.push_back(p);
    }
    return l;
}

PointList Snake::rotateSnake(double angle)
{
    double x,y;
    double rad=angle*PI/180;
    PointList newpoints;
    PointList::iterator i;
    Point p;

    //Get the barycenter of the contour
    Point Center;
    x=0;
    y=0;
    for(i=snake.begin();i!=snake.end();i++)
    {
        x=x+i->getX();
        y=y+i->getY();
    }
    Center.setLocation(x/(int)snake.size(),y/(int)snake.size());

    //Get the snake after the rotation by rad
    for(i=snake.begin();i!=snake.end();i++)
    {
        x=Center.getX() + (i->getX()-Center.getX())*cos(rad)-(i->getY()-Center.getY())*sin(rad);
        y=Center.getX() + (i->getX()-Center.getX())*sin(rad)+(i->getY()-Center.getY())*cos(rad);
        p.setLocation(x,y);
        newpoints.push_back(p);
    }
    return newpoints;
}





















