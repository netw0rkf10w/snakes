#include "stuffs.h"


using namespace std;


void print(const char * title, double** f, int width, int height)
{
    cout<<"******* "<<title<<"********"<<endl;
    for(int i=0;i<width;i++)
    {
        cout<<endl;
        for(int j=0;j<height;j++)
        {
            cout<< fixed << setprecision (2) <<f[i][j]<<"\t";
        }
    }
    cout<<endl;
}


void print3D(const char * title, double*** f, int width, int height)
{
    cout<<"******* "<<title<<"********"<<endl;
    for(int i=0;i<width;i++)
    {
        cout<<endl;
        for(int j=0;j<height;j++)
        {
            cout<< fixed << setprecision (2) <<"("<<f[i][j][0]<<","<<f[i][j][1]<<")"<<"\t";
        }
    }
    cout<<endl;
}


double** newArray(int width, int height)
{
    double** gradient;
    // VERY IMPORTANT !!!Dynamic allocations for the gradient and the flow matrix
    int dim_allouee1 = 0; // nombre d'éléments alloués avec succès sur la dimension 2
    try
    {
        // dimension 1 : tableau de 10 pointeurs vers des tableaux d'entiers
        gradient = new double * [ width ];
        // initialiser les 10 pointeurs à 0 (NULL)
        std::fill_n( gradient, width, static_cast<double*>( 0 ) );

        // dimension 2 : les tableaux d'entiers
        for ( dim_allouee1 = 0; dim_allouee1 < width; ++dim_allouee1)
        {
            gradient[ dim_allouee1 ] = new double[ height ];
        }

    }
    catch ( const std::bad_alloc & ) // erreur d'allocation
    {
        // désallouer tout ce qui a été alloué avec succès
        for ( int i = 0; i < dim_allouee1; ++i )
        {
            delete [] gradient[ i ];
        }
        delete [] gradient;
    }
    return gradient;
}


double*** newArray3D(int x, int y, int z)
{
    double*** tab;
    // VERY IMPORTANT !!!Dynamic allocations for the gradient and the flow matrix
    int assignedy = 0; //number of allocated elements of 2nd dimension
    int assignedz = 0; //number of allocated elements of 3rd dimension
    try
    {
        // dimension 1
        tab = new double**[ x ];
        std::fill_n( tab, x, static_cast<double**>( 0 ) );

        // dimension 2
        for( assignedy = 0; assignedy < x; ++assignedy)
        {
            tab[assignedy] = new double*[ y ];
            std::fill_n( tab[ assignedy ], y, static_cast<double*>( 0 ) );

            //dimension 3
            for( assignedz = 0; assignedz < y; ++assignedz)
            {
                tab[assignedy][assignedz] = new double[ z ];
            }

        }

    }
    catch ( const std::bad_alloc & ) // In case of allocation error
    {
        // Free those successfully allocated
        for(int i = 0; i < assignedy; ++i)
        {
            for(int j = 0; j < assignedz; ++j)
            {
                delete[] tab[i][j];
            }
            delete[] tab[i];
        }
        delete [] tab;
    }
    return tab;
}



double*** GVFC(double** f, int w, int h, int ITER, double mu)
{
    double mag2, temp, tempx, tempy;
    int count, x, y, XN_1, XN_2, YN_1, YN_2, XN, YN;

    /* define constants and create row-major double arrays */
    XN = w;
    YN = h;
    XN_1 = XN - 1;
    XN_2 = XN - 2;
    YN_1 = YN - 1;
    YN_2 = YN - 2;

    double** u = newArray(XN,YN);
    double** v = newArray(XN,YN);
    double** fx = newArray(XN,YN);
    double** fy = newArray(XN,YN);
    double** Lu = newArray(XN,YN);
    double** Lv = newArray(XN,YN);
    double** g  = newArray(XN,YN);
    double** c1 = newArray(XN,YN);
    double** c2 = newArray(XN,YN);
    double** b  = newArray(XN,YN);

    double*** gvf = newArray3D(w,h,2);

    /************** I: Normalize the edge map to [0,1] **************/
//	fmax = -1e10;
//	fmin = 1e10;
//	for (x=0; x<=YN*XN-1; x++)
//    {
//		fmax = std::max(fmax,f[x]);
//		fmin = std::min(fmin,f[x]);
//	}
//	assert(fmax!=fmin);


//	for (x=0; x<=YN*XN-1; x++)
//		f[x] = (f[x]-fmin)/(fmax-fmin);



    /**************** II: Compute edge map gradient *****************/
    /* I.1: Neumann boundary condition:
    *      zero normal derivative at boundary
    */
    /* Deal with corners */
    fx[0][0] = fy[0][0] = fx[XN_1][0] = fy[XN_1][0] = 0;
    fx[XN_1][YN_1] = fy[XN_1][YN_1] = fx[0][YN_1] = fy[0][YN_1] = 0;

    /* Deal with left and right column */
    for (y=1; y<=YN_2; y++)
    {
        fx[0][y] = fx[XN_1][y] = 0;
        fy[0][y] = 0.5 * (f[0][y+1] - f[0][y-1]);
        fy[XN_1][y] = 0.5 * (f[XN_1][y+1] - f[XN_1][y-1]);
    }

    /* Deal with top and bottom row */
    for (x=1; x<=XN_2; x++)
    {
        fy[x][0]= fy[x][YN_1] = 0;
        fx[x][0] = 0.5 * (f[x+1][0] - f[x-1][0]);
        fx[x][YN_1] = 0.5 * (f[x+1][YN_1] - f[x-1][YN_1]);
    }

    /* I.2: Compute interior derivative using central difference */
    for (y=1; y<=YN_2; y++)
    {
        for (x=1; x<=XN_2; x++)
        {
            /* NOTE: f is stored in column major */
            fx[x][y] = 0.5*(f[x+1][y]-f[x-1][y]);
            fy[x][y] = 0.5*(f[x][y+1]-f[x][y-1]);
        }
    }

    //For testing, compute and show the gradient image
//    double** gradient = newArray(w,h);
//    for (int y=0;y<h;y++)
//    {
//        for (int x=0;x<h;x++)
//        {
//            gradient[x][y] = sqrt(fx[x][y]*fx[x][y]+fy[x][y]*fy[x][y]);
//        }
//    }
//    print("Gradient",gradient,w,h);

    /******* III: Compute parameters and initializing arrays **********/
    temp = -1.0/(mu*mu);
    for (y=0; y<=YN_1; y++)
    {
        for (x=0; x<=XN_1; x++)
        {
            tempx = fx[x][y];
            tempy = fy[x][y];
            /* initial GVF vector */
            u[x][y] = tempx;
            v[x][y] = tempy;
            /* gradient magnitude square */
            mag2 = tempx*tempx + tempy*tempy;

            g[x][y] = mu;
            b[x][y] = mag2;

            c1[x][y] = b[x][y] * tempx;
            c2[x][y] = b[x][y] * tempy;
        }
    }

    /* free memory of fx and fy */
    for(int i=0;i<XN;i++)
    {
        delete[] fx[i];
        delete[] fy[i];
    }
    delete[] fx;
    delete[] fy;


    /************* Solve GVF = (u,v) iteratively ***************/
    for (count=1; count<=ITER; count++)
    {
        /* IV: Compute Laplace operator using Neuman condition */
        /* IV.1: Deal with corners */
        Lu[0][0] = (u[0][1] + u[1][0])*0.5 - u[0][0];
        Lv[0][0] = (v[0][1] + v[1][0])*0.5 - v[0][0];
        Lu[XN_1][0] = (u[XN_2][0] + u[XN_1][1])*0.5 - u[XN_1][0];
        Lv[XN_1][0] = (v[XN_2][0] + v[XN_1][1])*0.5 - v[XN_1][0];
        Lu[0][YN_1] = (u[1][YN_1] + u[0][YN_2])*0.5 - u[0][YN_1];
        Lv[0][YN_1] = (v[1][YN_1] + v[0][YN_2])*0.5 - v[0][YN_1];
        Lu[XN_1][YN_1] = (u[XN_2][YN_1] + u[XN_1][YN_2])*0.5 - u[XN_1][YN_1];
        Lv[XN_1][YN_1] = (v[XN_2][YN_1] + v[XN_1][YN_2])*0.5 - v[XN_1][YN_1];

        /* IV.2: Deal with left and right columns */
        for (y=1; y<=YN_2; y++)
        {
            Lu[0][y] = (2*u[1][y] + u[0][y-1] + u[0][y+1])*0.25 - u[0][y];
            Lv[0][y] = (2*v[1][y] + v[0][y-1] + v[0][y+1])*0.25 - v[0][y];
            Lu[XN_1][y] = (2*u[XN_2][y] + u[XN_1][y-1]
            + u[XN_1][y+1])*0.25 - u[XN_1][y];
            Lv[XN_1][y] = (2*v[XN_2][y] + v[XN_1][y-1]
            + v[XN_1][y+1])*0.25 - v[XN_1][y];
        }

        /* IV.3: Deal with top and bottom rows */
        for (x=1; x<=XN_2; x++)
        {
            Lu[x][0] = (2*u[x][1] + u[x-1][0] + u[x+1][0])*0.25 - u[x][0];
            Lv[x][0] = (2*v[x][1] + v[x-1][0] + v[x+1][0])*0.25 - v[x][0];
            Lu[x][YN_1] = (2*u[x][YN_2] + u[x-1][YN_1]
            + u[x+1][YN_1])*0.25 - u[x][YN_1];
            Lv[x][YN_1] = (2*v[x][YN_2] + v[x-1][YN_1]
            + v[x+1][YN_1])*0.25 - v[x][YN_1];
        }

        /* IV.4: Compute interior */
        for (y=1; y<=YN_2; y++)
        {
            for (x=1; x<=XN_2; x++)
            {
                Lu[x][y] = (u[x-1][y] + u[x+1][y]
                + u[x][y-1] + u[x][y+1])*0.25 - u[x][y];
                Lv[x][y] = (v[x-1][y] + v[x+1][y]
                + v[x][y-1] + v[x][y+1])*0.25 - v[x][y];
            }
        }

        /******** V: Update GVF ************/
        for (y=0; y<=YN_1; y++)
        {
            for (x=0; x<=XN_1; x++)
            {
                u[x][y] = (1- b[x][y])*u[x][y] + g[x][y]*Lu[x][y]*4 + c1[x][y];
                v[x][y] = (1- b[x][y])*v[x][y] + g[x][y]*Lv[x][y]*4 + c2[x][y];
            }
        }

//                /* print iteration number */
//                printf("%5d",count);
//                if (count%15 == 0)
//                    printf("\n");
    }
    //printf("\n");

    /* copy u,v to the output in column major order */
//    for (y=0; y<=YN_1; y++)
//    {
//        for (x=0; x<=XN_1; x++)
//        {
//            ou[x*YN+y] = u[x][y];
//            ov[x*YN+y] = v[x][y];
//        }
//    }
    for (y=0; y<=YN_1; y++)
    {
        for (x=0; x<=XN_1; x++)
        {
            gvf[x][y][0] = u[x][y];
            gvf[x][y][1] = v[x][y];
        }
    }

//    for(int i=0;i<w;i=i+30)
//    {
//        cout<<endl;
//        for(int j=0;j<h;j=j+30)
//        {
//            cout<< fixed << setprecision (5) <<"("<<gvf[i][j][0]<<","<<gvf[i][j][1]<<")"<<"\t";
//            //gradient[i][j]=(double)imageGrad->GetPixel(index);
//        }
//        cout<<endl;
//    }



    /* free all the array memory */
    for(int i=0;i<XN;i++)
    {
        delete[] u[i];
        delete[] v[i];
        delete[] Lu[i];
        delete[] Lv[i];
        delete[] g[i];
        delete[] c1[i];
        delete[] c2[i];
        delete[] b[i];
    }
    delete[] u;
    delete[] v;
    delete[] Lu;
    delete[] Lv;
    delete[] g;
    delete[] c1;
    delete[] c2;
    delete[] b;

    // return U and V arrays
    return gvf;
}



/*
 * @param f : image normalized in [0,1]
 * @param w : width of image
 * @param h : height of image
 * @param ITER : number of iterations
 * @param mu : iteration step
 * @return u[x,y] and v[x,y] arrays
 */
double*** GVF(double** f, int w, int h, int ITER, double mu)
{
    // create empty arrays
    double** u = newArray(w,h);
    double** v = newArray(w,h);
    double** fx = newArray(w,h);
    double** fy = newArray(w,h);
    double** Lu = newArray(w,h);
    double** Lv = newArray(w,h);
    double*** gvf = newArray3D(w,h,2);

    //double** gradient = newArray(w,h);

    for (int y=0;y<h;y++)
    {
        for (int x=0;x<w;x++)
        {
            fx[x][y] = 0;
            fy[x][y] = 0;
            u[x][y] = 0;
            v[x][y] = 0;
            Lu[x][y] = 0;
            Lv[x][y] = 0;
        }
    }

    // precompute edge-map (gradient)
    for (int y=1;y<(h-1);y++) {
        for (int x=1;x<(w-1);x++) {
            fx[x][y] = f[x+1][y]-f[x-1][y];
            fy[x][y] = f[x][y+1]-f[x][y-1];
            //gradient[x][y] = 255*sqrt(fx[x][y]*fx[x][y]+fy[x][y]*fy[x][y]);
        }
    }
    //print("Gradient",gradient,w,h);
    //save(gradient,w,h,"gradient.png");
    // iterative diffusion
    for(int loop=0;loop<ITER;loop++)
    {
        // compute laplacian of U and V
        for (int y=1;y<(h-1);y++)
        {
            for (int x=1;x<(w-1);x++)
            {
                Lu[x][y] = -u[x][y] + 0.25*(u[x-1][y]+u[x+1][y]+u[x][y-1]+u[x][y+1]);
                Lv[x][y] = -v[x][y] + 0.25*(v[x-1][y]+v[x+1][y]+v[x][y-1]+v[x][y+1]);
            }
        }

        // update U and V
        for (int y=0;y<h;y++)
        {
            for (int x=0;x<w;x++)
            {
                double gnorm2 = fx[x][y]*fx[x][y] + fy[x][y]*fy[x][y];

                u[x][y] += mu*4*Lu[x][y] - (u[x][y]-fx[x][y])*gnorm2;
                v[x][y] += mu*4*Lv[x][y] - (v[x][y]-fy[x][y])*gnorm2;
                gvf[x][y][0]=u[x][y];
                gvf[x][y][1]=v[x][y];
            }
        }
    }

    for(int i=0;i<w;i++)
    {
        delete[] u[i];
        delete[] v[i];
        delete[] fx[i];
        delete[] fy[i];
        delete[] Lu[i];
        delete[] Lv[i];
    }
    delete[] u;
    delete[] v;
    delete[] fx;
    delete[] fy;
    delete[] Lu;
    delete[] Lv;

    // return U and V arrays
    return gvf;
}



double** GVFMagnitude(double*** gvf, int w, int h)
{
    double** gvfm = newArray(w,h);
    for (int y=0;y<h;y++)
    {
        for (int x=0;x<w;x++)
        {
            gvfm[x][y] = (double)255*sqrt(gvf[x][y][0]*gvf[x][y][0] + gvf[x][y][1]*gvf[x][y][1]);
        }
    }
    return gvfm;
}


/**
 * Draws a line between two points p1(p1x,p1y) and p2(p2x,p2y).
 * This function is based on the Bresenham's line algorithm and is highly
 * optimized to be able to draw lines very quickly. There is no floating point
 * arithmetic nor multiplications and divisions involved. Only addition,
 * subtraction and bit shifting are used.
 *
 * Note that you have to define your own customized setPixel(x,y) function,
 * which essentially lights a pixel on the screen.
 */
PointList lineBresenham(int p1x, int p1y, int p2x, int p2y)
{
    int F, x, y;
    PointList l;

    if (p1x > p2x)  // Swap points if p1 is on the right of p2
    {
        swap(p1x, p2x);
        swap(p1y, p2y);
    }

    // Handle trivial cases separately for algorithm speed up.
    // Trivial case 1: m = +/-INF (Vertical line)
    if (p1x == p2x)
    {
        if (p1y > p2y)  // Swap y-coordinates if p1 is above p2
        {
            swap(p1y, p2y);
        }

        x = p1x;
        y = p1y;
        while (y <= p2y)
        {
            l.push_back(Point(x,y));
            y++;
        }
        return l;
    }
    // Trivial case 2: m = 0 (Horizontal line)
    else if (p1y == p2y)
    {
        x = p1x;
        y = p1y;

        while (x <= p2x)
        {
            l.push_back(Point(x,y));
            x++;
        }
        return l;
    }


    int dy            = p2y - p1y;  // y-increment from p1 to p2
    int dx            = p2x - p1x;  // x-increment from p1 to p2
    int dy2           = (dy << 1);  // dy << 1 == 2*dy
    int dx2           = (dx << 1);
    int dy2_minus_dx2 = dy2 - dx2;  // precompute constant for speed up
    int dy2_plus_dx2  = dy2 + dx2;


    if (dy >= 0)    // m >= 0
    {
        // Case 1: 0 <= m <= 1 (Original case)
        if (dy <= dx)
        {
            F = dy2 - dx;    // initial F

            x = p1x;
            y = p1y;
            while (x <= p2x)
            {
                l.push_back(Point(x,y));
                if (F <= 0)
                {
                    F += dy2;
                }
                else
                {
                    y++;
                    F += dy2_minus_dx2;
                }
                x++;
            }
        }
        // Case 2: 1 < m < INF (Mirror about y=x line
        // replace all dy by dx and dx by dy)
        else
        {
            F = dx2 - dy;    // initial F

            y = p1y;
            x = p1x;
            while (y <= p2y)
            {
                l.push_back(Point(x,y));
                if (F <= 0)
                {
                    F += dx2;
                }
                else
                {
                    x++;
                    F -= dy2_minus_dx2;
                }
                y++;
            }
        }
    }
    else    // m < 0
    {
        // Case 3: -1 <= m < 0 (Mirror about x-axis, replace all dy by -dy)
        if (dx >= -dy)
        {
            F = -dy2 - dx;    // initial F

            x = p1x;
            y = p1y;
            while (x <= p2x)
            {
                l.push_back(Point(x,y));
                if (F <= 0)
                {
                    F -= dy2;
                }
                else
                {
                    y--;
                    F -= dy2_plus_dx2;
                }
                x++;
            }
        }
        // Case 4: -INF < m < -1 (Mirror about x-axis and mirror
        // about y=x line, replace all dx by -dy and dy by dx)
        else
        {
            F = dx2 + dy;    // initial F

            y = p1y;
            x = p1x;
            while (y >= p2y)
            {
                l.push_back(Point(x,y));
                if (F <= 0)
                {
                    F += dx2;
                }
                else
                {
                    x++;
                    F += dy2_plus_dx2;
                }
                y--;
            }
        }
    }
    return l;
}


PointList getSegmentation(PointList snake)
{
    PointList seg;
    PointList::iterator it;
    PointList l;
    int N = snake.size();
    int x1,y1,x2,y2;
    it=snake.begin();
    for(int i=0; i< N-1; i++)
    {
        x1=it->getX();
        y1=it->getY();
        it++;
        x2=it->getX();
        y2=it->getY();
        l=lineBresenham(x1,y1,x2,y2);
        PointList::iterator itseg = seg.end();
        seg.splice(itseg,l);
    }
    seg.unique();
    return seg;
}

Point* listToPointArray(PointList l, int &N)
{
    N=(int)l.size();
    Point* p=new Point[N];
    PointList::iterator it = l.begin();
    for(int i=0; i< N; i++)
    {
        int x=it->getX();
        int y=it->getY();
        p[i]=Point(x,y);
        it++;
    }
    return p;
}

double distance2D(Point A, Point B)
{
    int ux = A.getX()-B.getX();
    int uy = A.getY()-B.getY();
    double un = ux*ux+uy*uy;
    return sqrt(un);
}

double distance(int i, Point* S, int sizeS, Point* T, int sizeT)
//S, T: arrays of pixels. This function compute the distance between the i-th pixel of S and its closest pixel in T
{
    if(i > sizeS) return 0;
    double min=Point(S[i]-T[0]).norm();
    for(int j=0;j<sizeT;j++)
    {
        if(min>Point(S[i]-T[j]).norm()) min=Point(S[i]-T[j]).norm();
    }
    return min;
}


double Pratt2(Point* S, int sizeS, Point* T, int sizeT)
{
    double sum=0;
    for(int i=0;i<sizeS;i++)
    {
        sum=sum+(double)1/(1+distance(i,S,sizeS,T,sizeT)*distance(i,S,sizeS,T,sizeT));
    }
    sum=sum/std::max(sizeS,sizeT);
    return sum;
}

void saveToFile(PointList l, QString fileName)
{
    QFile file(fileName);
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    //out << "Copyright 2012 Huu Dien Khue Le \n huudienkhue.le@gmail.com \n";
    if((int)l.size() > 0)
    {
        //out<< "N \t"<<(int)l.size()<<"\n";
        //out<< "x \t y \n";
        PointList::iterator j;
        for(j=l.begin();j!=l.end();j++)
        {
            out<<j->getX()<< "\t"<<j->getY()<<"\n";
        }
    }
    else out<< "The snake is null."<<"\n";
    // optional, as QFile destructor will already do it:
    file.close();
}

void saveToFile(QString l, QString fileName)
{
    QFile file(fileName);
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    //out << "Copyright 2012 Huu Dien Khue Le \n huudienkhue.le@gmail.com \n";
    if(!l.isNull())
    {
        out<<l<<"\n";
    }
    file.close();
}


//void save(double** im, int width, int height, const char* file)
//{
//    QImage image(width,height,QImage::Format_RGB32);
//    QRgb pixel;
//    for(int i=0;i<width;i++)
//    {
//        for(int j=0;j<height;j++)
//        {
//            pixel = qRgb(im[i][j],im[i][j],im[i][j]);
//            image.setPixel(i,j,pixel);
//        }
//    }
//    image.save(QString(file),"PNG");
//}




//double** normalizedImage(const char* file, int &width, int &height)
//{
////    QImage image = new QImage(file);
////    int width = image.width();
////    int height = image.height();

//    /* Read the image file*/
//    typedef itk::ImageFileReader<ImageType>		ReaderType;
//    ReaderType::Pointer reader = ReaderType::New();
//    reader->SetFileName(file);
//    reader->Update();

//    ImageType::Pointer imageGrad = reader->GetOutput();

//    //Get the width and the height of the image
//    width=imageGrad->GetLargestPossibleRegion().GetSize()[0];
//    height=imageGrad->GetLargestPossibleRegion().GetSize()[1];

//    /* Get the maximum and the minimum values for normalization */
//    typedef itk::MinimumMaximumImageCalculator <ImageType> ImageCalculatorFilterType;
//    ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
//    imageCalculatorFilter->SetImage(imageGrad);
//    imageCalculatorFilter->Compute();
//    double max = imageCalculatorFilter->GetMaximum();
//    double min = imageCalculatorFilter->GetMinimum();
//    //cout<<max<<"  "<<min<<endl;

//    /* Create the normalized gradien image, i.e. edge map */
//    double** gradient = newArray(width, height);

//    ImageType::IndexType index;
//    for(int i=0;i<width;i++)
//    {
//        for(int j=0;j<height;j++)
//        {
//            index[0]=i;
//            index[1]=j;
//            gradient[i][j]=(double)(imageGrad->GetPixel(index)-min)/(max-min);
//            //gradient[i][j]=(double)imageGrad->GetPixel(index);
//        }
//    }

//    return gradient;
//}




//double** normalizedEdgeMap(const char* file, int &width, int &height)
//{
////    QImage image = new QImage(file);
////    int width = image.width();
////    int height = image.height();

//    /* Read the image file*/
//    typedef itk::ImageFileReader<UnsignedCharImageType>		ReaderType;
//    ReaderType::Pointer reader = ReaderType::New();
//    reader->SetFileName(file);
//    reader->Update();

//    /* Create the gradient image */
//    typedef itk::GradientMagnitudeImageFilter< UnsignedCharImageType, ImageType > GradientFilter;
//    GradientFilter::Pointer fGrad = GradientFilter::New();
//    fGrad->SetInput(reader->GetOutput());
//    fGrad->Update();
//    ImageType::Pointer imageGrad = fGrad->GetOutput();


//        typedef itk::RescaleIntensityImageFilter<ImageType,UnsignedCharImageType>    RescaleFilterType;

//        RescaleFilterType::Pointer gradBeforeThresholdingRescale = RescaleFilterType::New();
//        gradBeforeThresholdingRescale->SetOutputMinimum(   0 );
//        gradBeforeThresholdingRescale->SetOutputMaximum( 255 );
//        gradBeforeThresholdingRescale->SetInput(fGrad->GetOutput());
//        gradBeforeThresholdingRescale->Update();
//        typedef itk::ImageFileWriter< UnsignedCharImageType >  GVFWriterType;
//        GVFWriterType::Pointer writer = GVFWriterType::New();
//        writer->SetFileName("outitk.png");
//        writer->SetInput( gradBeforeThresholdingRescale->GetOutput() );
//        writer->Update();

//    //Get the width and the height of the image
//    width=imageGrad->GetLargestPossibleRegion().GetSize()[0];
//    height=imageGrad->GetLargestPossibleRegion().GetSize()[1];

//    /* Get the maximum and the minimum values for normalization */
//    typedef itk::MinimumMaximumImageCalculator <ImageType> ImageCalculatorFilterType;
//    ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
//    imageCalculatorFilter->SetImage(imageGrad);
//    imageCalculatorFilter->Compute();
//    double max = imageCalculatorFilter->GetMaximum();
//    double min = imageCalculatorFilter->GetMinimum();

//    /* Create the normalized gradien image, i.e. edge map */
//    double** gradient;
//    // VERY IMPORTANT !!!Dynamic allocations for the gradient and the flow matrix
//    int dim_allouee1 = 0; // nombre d'éléments alloués avec succès sur la dimension 2
//    try
//    {
//        // dimension 1 : tableau de 10 pointeurs vers des tableaux d'entiers
//        gradient = new double * [ width ];
//        // initialiser les 10 pointeurs à 0 (NULL)
//        std::fill_n( gradient, width, static_cast<double*>( 0 ) );

//        // dimension 2 : les tableaux d'entiers
//        for ( dim_allouee1 = 0; dim_allouee1 < width; ++dim_allouee1)
//        {
//            gradient[ dim_allouee1 ] = new double[ height ];
//        }

//    }
//    catch ( const std::bad_alloc & ) // erreur d'allocation
//    {
//        // désallouer tout ce qui a été alloué avec succès
//        for ( int i = 0; i < dim_allouee1; ++i )
//        {
//            delete [] gradient[ i ];
//        }
//        delete [] gradient;
//    }

//    ImageType::IndexType index;
//    for(int i=0;i<width;i++)
//    {
//        for(int j=0;j<height;j++)
//        {
//            index[0]=i;
//            index[1]=j;
//            gradient[i][j]=(double)(imageGrad->GetPixel(index)-min)/(max-min);
//            //gradient[i][j]=(double)imageGrad->GetPixel(index);
//        }
//    }
//    return gradient;
//}




//void saveVectorImage(double*** gvf, int w, int h)
//{
//    typedef itk::Vector< float, 2 > PixelType;
//    typedef itk::Image< PixelType, 2 > VectorImageType;

//    // Then the image object can be created
//    VectorImageType::Pointer image = VectorImageType::New();
//    // The image region should be initialized
//    VectorImageType::IndexType start;
//    VectorImageType::SizeType size;
//    size[0] = w; // size along X
//    size[1] = h; // size along Y
//    start[0] = 0; // first index on X
//    start[1] = 0; // first index on Y
//    VectorImageType::RegionType region;
//    region.SetSize( size );
//    region.SetIndex( start );
//    // Pixel data is allocated
//    image->SetRegions( region );
//    image->Allocate();
//    // The image buffer is initialized to a particular value
//    VectorImageType::PixelType initialValue;
//    // A vector can initialize all its components to the
//    // same value by using the Fill() method.
//    initialValue.Fill( 0.0 );
//    // Now the image buffer can be initialized with this
//    // vector value.
//    image->FillBuffer( initialValue );

//    VectorImageType::IndexType pixelIndex;
//    VectorImageType::PixelType pixelValue;
//    for (int y=0;y<h;y++)
//    {
//        for (int x=0;x<w;x++)
//        {
//            pixelIndex[0] = x; // x position
//            pixelIndex[1] = y; // y position
//            //pixelValue[0] = gvf[x][y][0]; // x component
//            //pixelValue[1] = gvf[x][y][1]; // y component
//            pixelValue[0] = 1; // x component
//            pixelValue[1] = 1; // y component
//            image->SetPixel( pixelIndex, pixelValue );
//        }
//    }

//    typedef itk::ImageFileWriter< VectorImageType >  WriterType;
//    WriterType::Pointer writer = WriterType::New();
//    writer->SetFileName("vector.mhd");

//    // Software Guide : BeginCodeSnippet
//    writer->SetInput(image);
//    writer->Update();
//    // Software Guide : BeginLatex
//    //
//    // The Vector class inherits the operator \code{[]} from the
//    // \doxygen{FixedArray} class. This makes it possible to access the
//    // Vector's components using index notation.
//    //
//    // Software Guide : EndLatex
//    // Software Guide : BeginCodeSnippet

//    // Software Guide : EndCodeSnippet
//    // Software Guide : BeginLatex
//    //
//    // We can now store this vector in one of the image pixels by defining an
//    // index and invoking the \code{SetPixel()} method.
//    //
//    // Software Guide : EndLatex
//    // Software Guide : BeginCodeSnippet
//}




