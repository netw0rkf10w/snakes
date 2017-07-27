#include "application.h"
#include "image.h"
#include "itkThresholdImageFilter.h"
#include "stuffs.h"



Application::Application(QWidget *parent, Qt::WFlags flags):QMainWindow(parent, flags)
{
	ui.setupUi(this);

    ui.continuityBegin->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.curvatureBegin->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.imageBegin->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.flowBegin->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.balloonBegin->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.priorBegin->setValidator(new QDoubleValidator(-100,100,10,this));
    //Set default values for the coefficents
    ui.continuityBegin->setText("1.0");
    ui.curvatureBegin->setText("0.0");
    ui.imageBegin->setText("0.0");
    ui.flowBegin->setText("0.0");
    ui.balloonBegin->setText("0.0");
    ui.priorBegin->setText("1.0");

    ui.continuityStep->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.curvatureStep->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.imageStep->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.flowStep->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.balloonStep->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.priorStep->setValidator(new QDoubleValidator(-100,100,10,this));
    //Set default values for the coefficents
    ui.continuityStep->setText("0.02");
    ui.curvatureStep->setText("0.02");
    ui.imageStep->setText("0.02");
    ui.flowStep->setText("0.02");
    ui.balloonStep->setText("0.02");
    ui.priorStep->setText("0.02");

    ui.continuityEnd->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.curvatureEnd->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.imageEnd->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.flowEnd->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.balloonEnd->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.priorEnd->setValidator(new QDoubleValidator(-100,100,10,this));
    //Set default values for the coefficents
    ui.continuityEnd->setText("1.0");
    ui.curvatureEnd->setText("0.0");
    ui.imageEnd->setText("0.0");
    ui.flowEnd->setText("1.0");
    ui.balloonEnd->setText("0.0");
    ui.priorEnd->setText("1.0");

	ui.maxIterations->setValidator(new QIntValidator(0,100000,this));
	ui.autoAdaptLoop->setValidator(new QIntValidator(0,100,this));
	ui.minLength->setValidator(new QIntValidator(1,100,this));
	ui.maxLength->setValidator(new QIntValidator(1,100,this));

    ui.maxIterations->setText("250");
	ui.autoAdaptation->setChecked(false);	
	ui.autoAdaptLoop->setText("10");
	ui.minLength->setText("5");
	ui.maxLength->setText("10");
    ui.autoAdaptLoop->setEnabled(ui.autoAdaptation->isChecked());
    ui.minLength->setEnabled(ui.autoAdaptation->isChecked());
    ui.maxLength->setEnabled(ui.autoAdaptation->isChecked());

    ui.enableSPI->setChecked(false);
    ui.enableRotInvariant->setChecked(true);

	ui.wLength->setValidator(new QDoubleValidator(0,100,10,this));
	ui.wLength->setText("1.0");
    ui.wLength->setVisible(false);
    ui.wLengthLabel->setVisible(false);

    ui.angle->setValidator(new QDoubleValidator(-500,500,10,this));
    ui.angle->setText("0.0");
    ui.indexShift->setValidator(new QIntValidator(0,500,this));
    ui.indexShift->setText("0");

    ui.comboDistance->addItem("Euclidean");
    ui.comboDistance->addItem("Dynamic Time Warping");
    ui.comboDistance->addItem("Pratt");
    ui.comboDistance->addItem("Generic Discrepancy");
    ui.comboDistance->addItem("Larsson");

    ui.comboDistShifting->addItem("Euclidean");
    ui.comboDistShifting->addItem("DTW");
    ui.comboDistShifting->addItem("Pratt");
    ui.comboDistShifting->addItem("GDM");
    ui.comboDistShifting->addItem("Larsson");


    //Set default coefficient values for the 1st step of the auto-run scheme
    ui.continuityCoeff_2->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.curvatureCoeff_2->setValidator(new QDoubleValidator(-100,100,10,this));
    ui.flowCoeff_2->setValidator(new QDoubleValidator(-100,100,10,this));

    ui.maxIterations_2->setValidator(new QIntValidator(0,100000,this));
    ui.autoAdaptLoop_2->setValidator(new QIntValidator(0,100,this));
    ui.minLength_2->setValidator(new QIntValidator(1,100,this));
    ui.maxLength_2->setValidator(new QIntValidator(1,100,this));

    ui.continuityCoeff_2->setText("2.0");
    ui.curvatureCoeff_2->setText("0.0");
    ui.flowCoeff_2->setText("1.0");
    ui.maxIterations_2->setText("250");
    ui.autoAdaptation_2->setChecked(false);
    ui.autoAdaptLoop_2->setText("10");
    ui.minLength_2->setText("5");
    ui.maxLength_2->setText("10");
    ui.autoAdaptLoop_2->setEnabled(ui.autoAdaptation_2->isChecked());
    ui.minLength_2->setEnabled(ui.autoAdaptation_2->isChecked());
    ui.maxLength_2->setEnabled(ui.autoAdaptation_2->isChecked());

    ui.spaceForRebuilding->setValidator(new QIntValidator(0,100000,this));
    ui.spaceForRebuilding->setText("2");


	this->theFileName=QString::null;
	stop=false;

	pen= QPen(QBrush(Qt::red),3);
    //penblue= QPen(QBrush(Qt::blue),1);

    /* Initialize the main view */
	sceneMain = new Scene();
	sceneMain->update();
	sceneMain->setSceneRect(ui.viewMain->rect());
	//sceneMain->addLine(0,0,0,100,pen);
	//sceneMain->addLine(0,0,100,0,pen);
	ui.viewMain->setScene(sceneMain);
	ui.viewMain->show();	
	//ui.viewMain->setAlignment(Qt::AlignLeft);
	//ui.viewMain->setMouseTracking(true);
	//ui.viewMain->

//	sceneShape = new Scene();
//	sceneShape->enableMouseIteration(false);
//	sceneShape->update();
//	ui.viewShape->setScene(sceneShape);
//	ui.viewShape->show();

    isInPhase1 = false;
    isInPhase2 = false;
    isInPhaseChangeSP = false;

    /* Initialize the constructed view */
	sceneSnake = new Scene();	
	sceneSnake->enableMouseIteration(false);
	sceneSnake->update();
	ui.viewSnake->setScene(sceneSnake);
	ui.viewSnake->show();


    /* Initialize the gradient image view */
    scene1 = new Scene();
    scene1->enableMouseIteration(false);
    scene1->update();
    scene1->setSceneRect(ui.view1->rect());
    ui.view1->setScene(scene1);
    ui.view1->show();

    scene2 = new Scene();
    scene2->enableMouseIteration(false);
    scene2->update();
    scene2->setSceneRect(ui.view2->rect());
    ui.view2->setScene(scene2);
    ui.view2->show();

    scene3 = new Scene();
    scene3->enableMouseIteration(false);
    scene3->update();
    scene3->setSceneRect(ui.view3->rect());
    ui.view3->setScene(scene3);
    ui.view3->show();

    scene4 = new Scene();
    scene4->enableMouseIteration(false);
    scene4->update();
    scene4->setSceneRect(ui.view4->rect());
    ui.view4->setScene(scene4);
    ui.view4->show();


    scenePostProcessing = new Scene();
    scenePostProcessing->enableMouseIteration(false);
    scenePostProcessing->update();
    ui.viewPostProcessing->setScene(scenePostProcessing);
    ui.viewPostProcessing->show();

    ui.sigmaSlider->setMinimum(0);
    ui.sigmaSlider->setMaximum(100);
    ui.sigmaSlider->setValue(10);
    sigma=(double)ui.sigmaSlider->value()/10;
    ui.sigma->setText(QString::number(sigma));
    ui.sigmaMin->setText("0");
    ui.sigmaMax->setText("10");

    ui.thresholdSlider->setMinimum(0);
    ui.thresholdSlider->setMaximum(255);
    ui.thresholdSlider->setValue(230);
    threshold=(double)ui.thresholdSlider->value();
    ui.threshold->setText(QString::number(threshold));

	// Set up action signals and slots
	connect(ui.actionOpen, SIGNAL(triggered()), this, SLOT(slotOpen()));
	connect(ui.actionOpen_Shape, SIGNAL(triggered()), this, SLOT(slotOpenShape()));
	connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));
	connect(ui.actionSave, SIGNAL(triggered()), this, SLOT(slotSave()));
    connect(ui.actionRun, SIGNAL(triggered()), this, SLOT(slotGo()));
	connect(ui.actionStop, SIGNAL(triggered()), this, SLOT(slotStop()));
    connect(ui.actionContinue, SIGNAL(triggered()), this, SLOT(slotContinue()));
    connect(ui.autoAdaptation, SIGNAL(clicked()), this, SLOT(slotUpdateParameters()));
    connect(ui.autoAdaptation_2, SIGNAL(clicked()), this, SLOT(slotUpdateParameters()));
    connect(ui.enableSPI, SIGNAL(clicked()), this, SLOT(slotUpdateParameters()));
    connect(ui.enableRotInvariant, SIGNAL(clicked()), this, SLOT(slotUpdateParameters()));
    connect(ui.rotateButton, SIGNAL(clicked()), this, SLOT(slotRotateSnake()));
    connect(ui.shiftButton, SIGNAL(clicked()), this, SLOT(slotShiftSnake()));
    connect(ui.comboDistance, SIGNAL(currentIndexChanged(int)),this, SLOT(slotUpdateParameters()));
    connect(ui.comboDistShifting, SIGNAL(currentIndexChanged(int)),this, SLOT(slotUpdateParameters()));
    connect(ui.defineSnakeButton, SIGNAL(clicked()), this, SLOT(slotDefineInitSnake()));
    connect(ui.sigmaSlider, SIGNAL(valueChanged(int)), this, SLOT(slotUpdatePreprocessingParameters()));
    connect(ui.thresholdSlider, SIGNAL(valueChanged(int)), this, SLOT(slotUpdatePreprocessingParameters()));
    connect(ui.browseGroundTruthButton, SIGNAL(clicked()), this, SLOT(slotBrowseGroundTruth()));
    connect(ui.browseSnakeFileButton, SIGNAL(clicked()), this, SLOT(slotBrowseSnakeFile()));
    connect(ui.computeButton, SIGNAL(clicked()), this, SLOT(slotCompute()));
    connect(ui.rebuildButton, SIGNAL(clicked()), this, SLOT(slotRebuild()));

}

Application::~Application()
{
	//theImage->Delete();
}

void Application::slotUpdatePreprocessingParameters()
{

    sigma = (double)ui.sigmaSlider->value()/10;
    ui.sigma->setText(QString::number(sigma));


    threshold=(double)ui.thresholdSlider->value();
    ui.threshold->setText(QString::number(threshold));

    if(!theFileName.isNull())    display();
}

void Application::display()
{

    /* Create the gradient image */

    //    typedef itk::GradientMagnitudeImageFilter< ImageType, ImageType > GradientFilter;
    //    GradientFilter::Pointer fGrad = GradientFilter::New();
    //    fGrad->SetInput(reader->GetOutput());
    //    fGrad->Update();

    typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< ImageType, ImageType > GradientFilter;
    GradientFilter::Pointer fGrad = GradientFilter::New();
    fGrad->SetInput(reader->GetOutput());
    fGrad->SetSigma(sigma);
    fGrad->Update();
    //save this as the gradient image for snake processing later
    floatGradient = fGrad->GetOutput();


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
    floatFlow = fFlow->GetOutput();


    /* Rescale the created images for displaying */
    typedef itk::RescaleIntensityImageFilter<FloatImageType,UnsignedCharImageType>    RescaleFilterType;

    RescaleFilterType::Pointer gradRescale1 = RescaleFilterType::New();
    gradRescale1->SetOutputMinimum(   0 );
    gradRescale1->SetOutputMaximum( 255 );
    gradRescale1->SetInput(fGrad->GetOutput());
    gradRescale1->Update();

    RescaleFilterType::Pointer gradRescale2 = RescaleFilterType::New();
    gradRescale2->SetOutputMinimum(   0 );
    gradRescale2->SetOutputMaximum( 255 );
    gradRescale2->SetInput(thresholdFilter->GetOutput());
    gradRescale2->Update();


    RescaleFilterType::Pointer flowRescale = RescaleFilterType::New();
    flowRescale->SetOutputMinimum(   0 );
    flowRescale->SetOutputMaximum( 255 );
    flowRescale->SetInput(fFlow->GetOutput());
    flowRescale->Update();


    //Save the flow image
//        typedef itk::ImageFileWriter<UnsignedCharImageType>  WriterType;
//        WriterType::Pointer writer = WriterType::New();
//        //writer->SetInput(fFlow->GetOutput());
//        writer->SetFileName("1.png");
//        writer->SetInput(imageFlow);
//        try
//        {
//            writer->Update();
//        }
//        catch( itk::ExceptionObject & err )
//        {
//            std::cerr << "ExceptionObject caught !" << std::endl;
//            std::cerr << err << std::endl;
//            return;
//        }


    //UnsignedCharImageType::Pointer imageGrad=filter->GetOutput();

    UnsignedCharImageType::Pointer imageGrad1=gradRescale1->GetOutput();
    UnsignedCharImageType::Pointer imageGrad2=gradRescale2->GetOutput();
    UnsignedCharImageType::Pointer imageFlow=flowRescale->GetOutput();

    ImageType::IndexType index;
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {
            index[0]=i;
            index[1]=j;
            //gradImage.setPixel(i,j,imageGrad->GetPixel(index));
            pixel = qRgb(imageGrad1->GetPixel(index),imageGrad1->GetPixel(index),imageGrad1->GetPixel(index));
            gradImageBeforeThresholding.setPixel(i,j,pixel);
            pixel = qRgb(imageGrad2->GetPixel(index),imageGrad2->GetPixel(index),imageGrad2->GetPixel(index));
            gradImageAfterThresholding.setPixel(i,j,pixel);
            pixel = qRgb(imageFlow->GetPixel(index),imageFlow->GetPixel(index),imageFlow->GetPixel(index));
            flowImage.setPixel(i,j,pixel);
            //gradient[i][j]=imageGrad->GetPixel(index);
            //flow[i][j]=imageFlow->GetPixel(index);
        }
    }
    scene2->clear();
    scene3->clear();
    scene4->clear();
    ui.view2->setImage(gradImageBeforeThresholding);
    ui.view3->setImage(gradImageAfterThresholding);
    ui.view4->setImage(flowImage);
    scene2->update();
    scene3->update();
    scene4->update();
}


void Application::slotRotateSnake()
{
    double angle=ui.angle->text().toDouble();
    s.rotate(angle);
    this->updateAllViews();
}

void Application::slotShiftSnake()
{
    int i=ui.indexShift->text().toInt();
    if(i!=0)
    {
        s.shift(i);
        this->updateAllViews();
    }
}

void Application::slotDefineInitSnake()
{
    if((int)s.getSnake().size() > 0) this->sceneMain->setPoints(s.getSnake());
}

void Application::getParameters()
{
    continuityBegin=ui.continuityBegin->text().toDouble();
    curvatureBegin=ui.curvatureBegin->text().toDouble();
    imageBegin=ui.imageBegin->text().toDouble();
    flowBegin=ui.flowBegin->text().toDouble();
    balloonBegin=ui.balloonBegin->text().toDouble();
    priorBegin=ui.priorBegin->text().toDouble();

    continuityStep=ui.continuityStep->text().toDouble();
    curvatureStep=ui.curvatureStep->text().toDouble();
    imageStep=ui.imageStep->text().toDouble();
    flowStep=ui.flowStep->text().toDouble();
    balloonStep=ui.balloonStep->text().toDouble();
    priorStep=ui.priorStep->text().toDouble();

    continuityEnd=ui.continuityEnd->text().toDouble();
    curvatureEnd=ui.curvatureEnd->text().toDouble();
    imageEnd=ui.imageEnd->text().toDouble();
    flowEnd=ui.flowEnd->text().toDouble();
    balloonEnd=ui.balloonEnd->text().toDouble();
    priorEnd=ui.priorEnd->text().toDouble();

    maxIterations=ui.maxIterations->text().toInt();
	autoAdaptation=ui.autoAdaptation->isChecked();
	autoAdaptLoop=ui.autoAdaptLoop->text().toInt();
	minLength=ui.minLength->text().toInt();
	maxLength=ui.maxLength->text().toInt();
	enableSPI=ui.enableSPI->isChecked();
    autoShifting=ui.autoShifting->isChecked();
	s.useSPI=ui.enableSPI->isChecked();
    s.useRotInvar=ui.enableRotInvariant->isChecked();
	s.windowLength=ui.wLength->text().toInt();

    s.distanceType=(DistanceType)ui.comboDistance->currentIndex();
    s.distanceTypeSP=(DistanceType)ui.comboDistShifting->currentIndex();

    //Get parameters for the 1st step of the auto-run scheme
    continuityCoeff_2=ui.continuityCoeff_2->text().toDouble();
    curvatureCoeff_2=ui.curvatureCoeff_2->text().toDouble();
    flowCoeff_2=ui.flowCoeff_2->text().toDouble();
    maxIterations_2=ui.maxIterations_2->text().toInt();
    autoAdaptation_2=ui.autoAdaptation_2->isChecked();
    autoAdaptLoop_2=ui.autoAdaptLoop_2->text().toInt();
    minLength_2=ui.minLength_2->text().toInt();
    maxLength_2=ui.maxLength_2->text().toInt();


}


void Application::slotStop()
{
	stop=true;
}

void Application::updateMainView()
{
    if((int)s.getSnake().size() > 0)
    {
        sceneMain->clear();
        sceneMain->drawLine(s.getSnake());
        //sceneMain->drawLine(s.getSnake());
    }
}

void Application::updateSnakeView()
{
    if((int)s.getSnake().size() > 0)
    {
        sceneSnake->clear();
        sceneSnake->drawLine(iDFT2Points(s.DFTSnake(),s.getSnake().size()),Qt::blue);
    }
    if((int)theShape.size() > 0)
        sceneSnake->drawLine(iDFT2Points(s.DFTShape(),theShape.size()),Qt::red);
}

//void Application::updateShapeView()
//{
//    if((int)theShape.size() > 0)
//    {
//        sceneShape->clear();
//        sceneShape->drawLine(iDFT2Points(s.DFTShape(),s.getShape().size()));
//    }

//}

void Application::updateAllViews()
{
    updateMainView();
    updateSnakeView();
//    updateShapeView();
}

void Application::slotUpdateParameters()
{
    s.useRotInvar=ui.enableRotInvariant->isChecked();
    ui.enableSPI->setEnabled(s.useRotInvar);
    s.useSPI=ui.enableSPI->isChecked();

    ui.autoAdaptLoop->setEnabled(ui.autoAdaptation->isChecked());
    ui.minLength->setEnabled(ui.autoAdaptation->isChecked());
    ui.maxLength->setEnabled(ui.autoAdaptation->isChecked());

    ui.autoAdaptLoop_2->setEnabled(ui.autoAdaptation_2->isChecked());
    ui.minLength_2->setEnabled(ui.autoAdaptation_2->isChecked());
    ui.maxLength_2->setEnabled(ui.autoAdaptation_2->isChecked());


    ui.wLengthLabel->setVisible((DistanceType)ui.comboDistance->currentIndex()==DTW);
    ui.wLength->setVisible((DistanceType)ui.comboDistance->currentIndex()==DTW);

	//and display the the constructed snake from its Fourier coefficients
//	if(s.getShape().size()>0)
//	{
//        updateShapeView();
//	}
	if(s.getSnake().size()>0)
	{
        updateSnakeView();
	}	
}

void Application::disallowChangingParameters()
{
    ui.groupGeneralParameters->setEnabled(false);
    ui.groupFirstStepParameters->setEnabled(false);
}

void Application::allowChangingParameters()
{
    ui.groupGeneralParameters->setEnabled(true);
    ui.groupFirstStepParameters->setEnabled(true);
}

void Application::slotTest()
{
	/*PointList l=s.shiftSnake(3);
	QGraphicsScene* scene2 = new QGraphicsScene();		
	scene2->update();
	ui.graphicsView_2->setScene(scene2);
	ui.graphicsView_2->show();	
	scene2->clear();
	addPointListToScene(scene2,l,pen);*/
}

void Application::slotOpenShape()
{
    QString fileName = QFileDialog::getOpenFileName(this,tr("Open Shape Image"), QDir::currentPath(),tr("Text files (*.txt);;Image files (*.png *.bmp *.jpg);;All files (*.*)"));
    if (fileName.isEmpty())
        return;
    QFileInfo fi(fileName);
    QString ext = fi.completeSuffix();
    if(ext == "txt")
    {
        QFile file(fileName);
        if(!file.open(QIODevice::ReadOnly))
        {
            QMessageBox::information(0, "error", file.errorString());
        }

        QTextStream in(&file);

        theShape.clear();
        while(!in.atEnd())
        {
            QString line = in.readLine();
            //std::cout<<line.toUtf8().constData()<<std::endl;
            QStringList fields = line.split("\t");
            Point p;
            if(!fields.isEmpty())
            {
                std::cout<<fields.first().toUtf8().constData()<<" "<<fields.last().toUtf8().constData()<<std::endl;
                p.setX(atoi(fields.first().toUtf8().constData()));
                p.setY(atoi(fields.last().toUtf8().constData()));
                theShape.push_back(p);
                //std::cout<<p.getX()<<" "<<p.getY()<<std::endl;
            }
        }

        file.close();
    }
    else
    {
        QImage image(fileName);
        if (image.isNull())
        {
            QMessageBox::information(this, tr("Image Viewer"),tr("Cannot load %1.").arg(fileName));
            return;
        }
        QGraphicsBlurEffect *effect = new QGraphicsBlurEffect(this);
        setGraphicsEffect(true? effect : 0);
        ShapeDlg *dlg=new ShapeDlg;
        dlg->setTheFileName(fileName);
        dlg->exec();
        delete effect;

        if((int)dlg->getShape().size()>0)
        {
            //Retrieve the shape reference
            theShape.clear();
            theShape = dlg->getShape();
        }
    }

    if(theShape.size() > 0)
    {
        s.getShape().clear();
        s.setShape(theShape);
        this->sceneMain->setNumOfPoints(theShape.size());

        //and display the the constructed snake from its Fourier coefficients
        sceneSnake->clear();
        sceneSnake->drawLine(iDFT2Points(s.DFTShape(),theShape.size()));
    }
}


void Application::slotOpen()
{
	
    QString fileName = QFileDialog::getOpenFileName(this,tr("Open Image"), QDir::currentPath(),tr("Image files (*.png *.bmp *.jpg);;All files (*.*)"));
    if (!fileName.isEmpty()) 
	{
        QImage image(fileName);
        if (image.isNull()) 
		{
            QMessageBox::information(this, tr("Image Viewer"),tr("Cannot load %1.").arg(fileName));
            return;
        }	
        theFileName=fileName;
        /*QGraphicsScene* scene=new QGraphicsScene;
        scene->addPixmap(QPixmap::fromImage(QImage(theFileName)));
        ui.viewMain->setScene(scene);
        ui.viewMain->show();*/
        //image = QImage(theFileName);
        ui.viewMain->setImage(image);
        ui.view1->setImage(image);
        width = image.width();
        height = image.height();
        sceneMain->setSceneRect(0,0,width,height);
        scene1->setSceneRect(0,0,width,height);
        scene2->setSceneRect(0,0,width,height);
        scene3->setSceneRect(0,0,width,height);
        scene4->setSceneRect(0,0,width,height);

        gradImageBeforeThresholding= QImage(width,height,QImage::Format_RGB32);
        gradImageAfterThresholding= QImage(width,height,QImage::Format_RGB32);
        flowImage= QImage(width,height,QImage::Format_RGB32);

        /* Create the gradient image */
        reader = ReaderType::New();
        reader->SetFileName(theFileName.toUtf8().constData());
        reader->Update();
        display();
	}
}


void Application::run()
{
    statusBar()->showMessage("Processing...");
    if(isInPhase1)
    {
        //Run without prior energy to get the contour and then deduce the optimal starting point
        statusBar()->showMessage("Phase 1: Get the optimal starting point position...");
        while(stop==false && loop< maxIterations_2 && s.step2(continuityCoeff_2, curvatureCoeff_2, 0.0, flowCoeff_2, 0.0,0.0))
        {
            if (autoAdaptation_2 && loop % autoAdaptLoop_2==0)
            {
                s.removeOverlappingPoints(minLength_2);
                s.addMissingPoints(maxLength_2);
            }

            updateMainView();
            updateSnakeView();
            loop++;
            ui.iterations->setText(QString::number(loop));
            qApp->processEvents();
        }
        // if the snake is stopped, and this was not caused by clicking on the Stop button, then the process passes to Phase 1.5.
        // Otherwise, it is still in phase 1
        if(stop==false)
        {
            isInPhase1 = false;
            isInPhaseChangeSP = true;
        }
    }

    if(isInPhaseChangeSP)
    {
        statusBar()->showMessage("Phase 1: Get the optimal starting point position...");
        if (autoAdaptation_2) s.reSample((int)theShape.size());
        //Get the optimal starting point position
        int start=s.getOptimalStartingPoint();
        ui.OptimalStart->setText(QString::number(start));

        //Reset the snake
        if((int)this->sceneMain->getPoints().size()>0)
            s.setSnake(this->sceneMain->getPoints());

        //Change the starting position to the optimal one
        s.shift(start);
        isInPhase1 = false;
        isInPhaseChangeSP = false;
    }

    isInPhase2 = true;
    loop=0;
    if(isInPhase2)
    {
        if (autoAdaptation) // If auto adaptation,  no prior energy
        {
            while(stop==false && loop< maxIterations && s.step2(continuityCoeff, curvatureCoeff, imageCoeff, flowCoeff, balloonCoeff,0))
            {
                if (loop%autoAdaptLoop==0)
                {
                    s.removeOverlappingPoints(minLength);
                    s.addMissingPoints(maxLength);
                }
                updateMainView();
                updateSnakeView();

                qApp->processEvents();
                loop++;
                ui.iterations->setText(QString::number(loop));

                if((continuityStep > 0 && continuityCoeff < continuityEnd) || (continuityStep < 0 && continuityCoeff > continuityBegin)) continuityCoeff+= continuityStep;
                if((curvatureStep > 0 && curvatureCoeff < curvatureEnd) || (curvatureStep < 0 && curvatureCoeff > curvatureBegin)) curvatureCoeff+= curvatureStep;
                if((imageStep > 0 && imageCoeff < imageEnd) || (imageStep < 0 && imageCoeff > imageBegin)) imageCoeff+= imageStep;
                if((flowStep > 0 && flowCoeff < flowEnd) || (flowStep < 0 && flowCoeff > flowBegin)) flowCoeff+= flowStep;
                if((balloonStep > 0 && balloonCoeff < balloonEnd) || (balloonStep < 0 && balloonCoeff > balloonBegin)) balloonCoeff+= balloonStep;
                if((priorStep > 0 && priorCoeff < priorEnd) || (priorStep < 0 && priorCoeff > priorBegin)) priorCoeff+= priorStep;

                ui.continuityCoeff->setText(QString::number(continuityCoeff));
                ui.curvatureCoeff->setText(QString::number(curvatureCoeff));
                ui.imageCoeff->setText(QString::number(imageCoeff));
                ui.flowCoeff->setText(QString::number(flowCoeff));
                ui.balloonCoeff->setText(QString::number(balloonCoeff));
                ui.priorCoeff->setText(QString::number(priorCoeff));
            }
        }
        else
        {
            //If no adaptation and the number of points of the snake is greater than the one of the reference, then resample the snake
            /*int n=s.getShape().size();
            s.reSample(n);*/
            while(stop==false && loop< maxIterations && s.step2(continuityCoeff, curvatureCoeff, imageCoeff, flowCoeff, balloonCoeff,priorCoeff))
            {
                //Make the snake look better
                int n = s.getSnake().size();
                if (loop%autoAdaptLoop==0)
                {
                    if(s.removeOverlappingPoints2(minLength) > 0)
                        s.rebuild2(n);
                }

                updateMainView();
                updateSnakeView();
                loop++;
                ui.iterations->setText(QString::number(loop));
                qApp->processEvents();

                ////automatically change the starting point
                if(autoShifting)
                {
                    int start=s.getOptimalStartingPoint();
                    PointList l=s.shiftSnake(start);
                    s.getSnake().clear();
                    s.setSnake(l);
                }

                if((continuityStep > 0 && continuityCoeff < continuityEnd) || (continuityStep < 0 && continuityCoeff > continuityBegin)) continuityCoeff+= continuityStep;
                if((curvatureStep > 0 && curvatureCoeff < curvatureEnd) || (curvatureStep < 0 && curvatureCoeff > curvatureBegin)) curvatureCoeff+= curvatureStep;
                if((imageStep > 0 && imageCoeff < imageEnd) || (imageStep < 0 && imageCoeff > imageBegin)) imageCoeff+= imageStep;
                if((flowStep > 0 && flowCoeff < flowEnd) || (flowStep < 0 && flowCoeff > flowBegin)) flowCoeff+= flowStep;
                if((balloonStep > 0 && balloonCoeff < balloonEnd) || (balloonStep < 0 && balloonCoeff > balloonBegin)) balloonCoeff+= balloonStep;
                if((priorStep > 0 && priorCoeff < priorEnd) || (priorStep < 0 && priorCoeff > priorBegin)) priorCoeff+= priorStep;

                ui.continuityCoeff->setText(QString::number(continuityCoeff));
                ui.curvatureCoeff->setText(QString::number(curvatureCoeff));
                ui.imageCoeff->setText(QString::number(imageCoeff));
                ui.flowCoeff->setText(QString::number(flowCoeff));
                ui.balloonCoeff->setText(QString::number(balloonCoeff));
                ui.priorCoeff->setText(QString::number(priorCoeff));
            }
        }
    }
    statusBar()->showMessage("Done !");
}

void Application::slotGo()
{	
//    //PointList l=this->sceneMain->getPoints();
//    int n=this->sceneMain->getPoints().size();
//    complex* c = points2complex(this->sceneMain->getPoints());
//    complex* f = fillPoints(c, n, 4);
//    PointList l = complex2points(f,n*4);
//    sceneSnake->clear();
//    //sceneSnake->drawLine(iDFT2Points(s.DFTSnake(),s.getSnake().size()));
//    sceneSnake->drawLine(l);
//    //qApp->processEvents();

    if( theFileName.isNull() )
    {
        QMessageBox msgBox;
        msgBox.setText("Please open an image!");
        msgBox.exec();
    }
    else
    {
        if((int)theShape.size()>0)
        {
            s= Snake(theFileName.toUtf8().constData(),theShape,sigma,threshold);
        }
        else
        {
            //QMessageBox::information(this, tr("No shape reference"),tr("No shape reference. The prior coefficient will be set to 0."));
            ui.priorCoeff->setText("0.0");
            priorCoeff=ui.priorCoeff->text().toDouble();
            s= Snake(theFileName.toUtf8().constData(),sigma,threshold);
        }
        if((int)this->sceneMain->getPoints().size()>0)
        {
            //makeSureAllPointsAreInTheImage and all the poinst are set to original (i.e. points of the initial snake)
            PointList::iterator i;
            PointList l = this->sceneMain->getPoints();
            QImage im(theFileName);
            for(i=l.begin();i!=l.end();i++)
            {
                if (i->getX() < 1) i->setX(1);
                if (i->getX() >= (im.width()-1)) i->setX(im.width()-2);
                if (i->getY() < 1) i->setY(1);
                if (i->getY() >= (im.height()-1)) i->setY(im.height()-2);
            }
            s.setSnake(l);
        }
        getParameters();
        continuityCoeff=ui.continuityBegin->text().toDouble();
        curvatureCoeff=ui.curvatureBegin->text().toDouble();
        imageCoeff=ui.imageBegin->text().toDouble();
        flowCoeff=ui.flowBegin->text().toDouble();
        balloonCoeff=ui.balloonBegin->text().toDouble();
        priorCoeff=ui.priorBegin->text().toDouble();
        //disallowChangingParameters();

        isInPhase1 = ui.autoRun->isChecked();
        isInPhaseChangeSP = false;

        loop=0;
        stop=false;
        run();
//        statusBar()->showMessage("Processing...");
//        loop=0;
//        stop=false;
//        if(!ui.autoRun->isChecked())
//        {
//            if (autoAdaptation) // If auto adaptation,  no prior energy
//            {
//                while(stop==false && loop< maxIterations && s.step2(continuityCoeff, curvatureCoeff, imageCoeff, flowCoeff, balloonCoeff,0))
//                {
//                    if (loop%autoAdaptLoop==0)
//                    {
//                        s.removeOverlappingPoints(minLength);
//                        s.addMissingPoints(maxLength);
//                    }
//                    updateMainView();
//                    updateSnakeView();

//                    qApp->processEvents();
//                    loop++;
//                    ui.iterations->setText(QString::number(loop));
//                }
//            }
//            else
//            {
//                //If no adaptation and the number of points of the snake is greater than the one of the reference, then resample the snake
//                /*int n=s.getShape().size();
//                s.reSample(n);*/
//                while(stop==false && loop< maxIterations && s.step2(continuityCoeff, curvatureCoeff, imageCoeff, flowCoeff, balloonCoeff,priorCoeff))
//                {
//                    updateMainView();
//                    updateSnakeView();
//                    loop++;
//                    ui.iterations->setText(QString::number(loop));
//                    qApp->processEvents();

//                    ////automatically change the starting point
//                    if(autoShifting)
//                    {
//                        int start=s.getOptimalStartingPoint();
//                        PointList l=s.shiftSnake(start);
//                        s.getSnake().clear();
//                        s.setSnake(l);
//                    }
//                }
//            }

//        }
//        else //If the auto-run mode is set
//        {
//            //Run without prior energy to get the contour and then deduce the optimal starting point
//            statusBar()->showMessage("Get the optimal starting point position...");
//            while(stop==false && loop< maxIterations_2 && s.step2(continuityCoeff_2, curvatureCoeff_2, 0.0, flowCoeff_2, 0.0,0.0))
//            {
//                if (autoAdaptation_2 && loop % autoAdaptLoop_2==0)
//                {
//                    s.removeOverlappingPoints(minLength_2);
//                    s.addMissingPoints(maxLength_2);
//                }

//                updateMainView();
//                updateSnakeView();
//                loop++;
//                ui.iterations->setText(QString::number(loop));
//                qApp->processEvents();
//            }
//            if (autoAdaptation_2) s.reSample((int)theShape.size());
//            //Get the optimal starting point position
//            int start=s.getOptimalStartingPoint();
//            ui.OptimalStart->setText(QString::number(start));

//            //Reset the snake
//            if((int)this->sceneMain->getPoints().size()>0)
//                s.setSnake(this->sceneMain->getPoints());

//            //Change the starting position to the optimal one
//            s.shift(start);

//            //The parameters may have been modified by the user
//            getParameters();

//            //Allow to run
//            stop=false;
//            loop=0;
//            //Process the second step, add prior energy
//            statusBar()->showMessage("Processing... Please wait...");
//            while(stop==false && loop< maxIterations)
//            {
//                s.step2(continuityCoeff, curvatureCoeff, imageCoeff, flowCoeff, balloonCoeff,priorCoeff);
//                updateMainView();
//                updateSnakeView();
//                loop++;
//                ui.flow->setText(QString::number(flowCoeff));
//                if(flowCoeff<1) flowCoeff+=flowStep;
//                ui.iterations->setText(QString::number(loop));
//                qApp->processEvents();
//            }
//        }

//        statusBar()->showMessage("Done !");
//        allowChangingParameters();
    }
}


void Application::slotContinue()
{

	int n=s.getShape().size();
    s.reSample(n);

	stop=false;
    getParameters();
    if((continuityStep > 0 && continuityCoeff < continuityEnd) || (continuityStep < 0 && continuityCoeff > continuityBegin)) continuityCoeff+= continuityStep;
    if((curvatureStep > 0 && curvatureCoeff < curvatureEnd) || (curvatureStep < 0 && curvatureCoeff > curvatureBegin)) curvatureCoeff+= curvatureStep;
    if((imageStep > 0 && imageCoeff < imageEnd) || (imageStep < 0 && imageCoeff > imageBegin)) imageCoeff+= imageStep;
    if((flowStep > 0 && flowCoeff < flowEnd) || (flowStep < 0 && flowCoeff > flowBegin)) flowCoeff+= flowStep;
    if((balloonStep > 0 && balloonCoeff < balloonEnd) || (balloonStep < 0 && balloonCoeff > balloonBegin)) balloonCoeff+= balloonStep;
    if((priorStep > 0 && priorCoeff < priorEnd) || (priorStep < 0 && priorCoeff > priorBegin)) priorCoeff+= priorStep;

    ui.continuityCoeff->setText(QString::number(continuityCoeff));
    ui.curvatureCoeff->setText(QString::number(curvatureCoeff));
    ui.imageCoeff->setText(QString::number(imageCoeff));
    ui.flowCoeff->setText(QString::number(flowCoeff));
    ui.balloonCoeff->setText(QString::number(balloonCoeff));
    ui.priorCoeff->setText(QString::number(priorCoeff));
    run();
    
//	if(autoShifting && !ui.autoRun->isChecked())
//	{
//        int start2=s.getOptimalStartingPoint();
//        ui.OptimalStart->setText(QString::number(start2));
//		PointList l=s.shiftSnake(start2);
//		s.getSnake().clear();
//		s.setSnake(l);
//	}

//	int loop=0;
//    disallowChangingParameters();
//	statusBar()->showMessage("Processing...");
//	while(stop==false && loop< maxIterations && s.step2(continuityCoeff, curvatureCoeff, imageCoeff, flowCoeff, balloonCoeff,priorCoeff))
//	{
		
//        updateMainView();
//        updateSnakeView();
//		loop++;
//        ui.iterations->setText(QString::number(loop));
//		qApp->processEvents();
		
//		//automatically change the starting point
//        if(autoShifting && !ui.autoRun->isChecked())
//		{
//            int start=s.getOptimalStartingPoint();
//			PointList l=s.shiftSnake(start);
//			s.getSnake().clear();
//			s.setSnake(l);
//		}
//	}
//	statusBar()->showMessage("Done !");
//    allowChangingParameters();
}


void Application::slotExit()
{
	//ui.widget->GetRenderWindow()->GetInteractor()->TerminateApp(); //Close interactor before closing app. If not, app's windows process will not be killed.
	qApp->exit();
}

void Application::slotSave()
{
	//ui.viewMain->scene()->addLine(0, 0, -100, 200, QPen(QBrush(Qt::black),1));
	/*QString fileName = QFileDialog::getSaveFileName(this,tr("Save Image"),QDir::currentPath(),tr("PNG (*.png);; JPG (*.jpg) ;; BMP (*.bmp) ;; TIFF(*.TIFF) ;; All Files (*.*)") );
    if( !fileName.isNull() )
    {
        writeFloatImageToUnsignedCharFile(getTheItkImage(), fileName.toUtf8().constData());
    }*/
    //QString fileName = QFileDialog::getSaveFileName(this,tr("Export data"),QDir::currentPath(),tr("TEXT (*.txt)") );
    QString dir = QFileDialog::getExistingDirectory(this, tr("Choose directory"),QDir::currentPath(),QFileDialog::ShowDirsOnly|QFileDialog::DontResolveSymlinks);
    if( !dir.isNull() )
    {
        //Save the image with the final snake
        QString imagefile = dir + "/result.png";
        QPixmap pixMap = QPixmap::grabWidget(ui.viewMain);
        pixMap.save(imagefile);

        QString fileName = dir + "/result_snake.txt";
        saveToFile(s.getSnake(),fileName);
        //PointList seglist = getSegmentation(s.getSnake());
        PointList seglist = s.getSnake();
        fileName = dir + "/result_segmentation.txt";
        saveToFile(seglist,fileName);
        fileName = dir + "/result_groundtruth.txt";
        saveToFile(groundTruth,fileName);
        fileName = dir + "/result_Pratt.txt";
        saveToFile(ui.pratt->text(),fileName);
    }

}


void Application::slotBrowseGroundTruth()
{
    QString fileName = QFileDialog::getOpenFileName(this,tr("Open ground truth data file"), QDir::currentPath(),tr("Text files (*.txt)"));
    if (fileName.isEmpty())
        return;
    QFileInfo fi(fileName);
    QString ext = fi.completeSuffix();
    if(ext == "txt")
    {
        ui.groundTruthFile->setText(fileName);
        QFile file(fileName);
        if(!file.open(QIODevice::ReadOnly))
        {
            QMessageBox::information(0, "Error", file.errorString());
        }

        QTextStream in(&file);
        PointList l;
        while(!in.atEnd())
        {

            QString line = in.readLine();
            //std::cout<<line.toUtf8().constData()<<std::endl;
            QStringList fields = line.split("\t");
            int x,y;
            //Point p;
            if(!fields.isEmpty())
            {
                //std::cout<<fields.first().toUtf8().constData()<<" "<<fields.last().toUtf8().constData()<<std::endl;
                x=atoi(fields.first().toUtf8().constData());
                y=atoi(fields.last().toUtf8().constData());
                l.push_back(Point(x,y));
            }
        }
        file.close();
        groundTruth.clear();
        groundTruth = getSegmentation(l);
        scenePostProcessing->setSceneRect(0,0,width,height);
        scenePostProcessing->update();
        scenePostProcessing->clear();
        QPen penblue(QBrush(Qt::blue),1);
        scenePostProcessing->draw(groundTruth,penblue);

    }
}

void Application::slotBrowseSnakeFile()
{
    QString fileName = QFileDialog::getOpenFileName(this,tr("Open snake data file"), QDir::currentPath(),tr("Text files (*.txt)"));
    if (fileName.isEmpty())
        return;
    QFileInfo fi(fileName);
    QString ext = fi.completeSuffix();
    if(ext == "txt")
    {
        ui.snakeFile->setText(fileName);
        QFile file(fileName);
        if(!file.open(QIODevice::ReadOnly))
        {
            QMessageBox::information(0, "Error", file.errorString());
        }

        QTextStream in(&file);
        PointList l;
        while(!in.atEnd())
        {

            QString line = in.readLine();
            //std::cout<<line.toUtf8().constData()<<std::endl;
            QStringList fields = line.split("\t");
            int x,y;
            //Point p;
            if(!fields.isEmpty())
            {
                //std::cout<<fields.first().toUtf8().constData()<<" "<<fields.last().toUtf8().constData()<<std::endl;
                x=atoi(fields.first().toUtf8().constData());
                y=atoi(fields.last().toUtf8().constData());
                l.push_back(Point(x,y));
            }
        }
        file.close();

        aSnake.clear();
        aSnake = getSegmentation(l);
        //scenePostProcessing->setSceneRect(0,0,width,height);
        //scenePostProcessing->update();
        //scenePostProcessing->clear();
        //QPen penblue(QBrush(Qt::blue),1);
        //scenePostProcessing->draw(groundTruth,penblue);

    }
}


void Application::slotCompute()
{
    PointList seglist;
    if(ui.snakeFile->text().isEmpty())
        seglist = getSegmentation(s.getSnake());
    else seglist = aSnake;

    if((int)seglist.size() <= 0 || (int)groundTruth.size() <= 0) return;

    QImage image = QImage(theFileName);
    ui.viewPostProcessing->setImage(image);
    QPen pen(QBrush(Qt::red),1);
    QPen pen2(QBrush(Qt::blue),1);
    scenePostProcessing->clear();
    scenePostProcessing->draw(seglist,pen);
    scenePostProcessing->draw(groundTruth,pen2);
    scenePostProcessing->update();

    int nS,nT;
    Point* pS = listToPointArray(seglist,nS);
    Point* pT = listToPointArray(groundTruth,nT);
    double pratt = Pratt2(pS,nS,pT,nT);
    ui.pratt->setText(QString::number(pratt));
}

void Application::slotRebuild()
{
    int space=ui.spaceForRebuilding->text().toInt();
    s.rebuild(space);
    updateMainView();
    updateSnakeView();
    //qApp->processEvents();
}
