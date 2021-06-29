#ifndef APPLICATION_H
#define APPLICATION_H

#include "point.h" 
#include "snake.h"
#include "ui_application.h"
#include "dialog.h"
#include "scene.h"
//#include "ui_dialog.h"
//#include "itkGradientImageFilter.h"

#include <QtGui>
#include <QFlags>


//typedef unsigned char						UnsignedCharPixelType;
//typedef float								FloatPixelType;
////typedef itk::RGBPixel<unsigned char>		RGBPixelType;
//typedef itk::Image<FloatPixelType, 2>		FloatImageType;
//typedef itk::Image<UnsignedCharPixelType,2>	UnsignedCharImageType;
//typedef FloatImageType						ImageType;
////typedef itk::ImageFileReader<ImageType>		ReaderType;

typedef itk::ImageFileReader<ImageType>		ReaderType;

class Application : public QMainWindow
{
	Q_OBJECT

public:	
    Application(QWidget *parent = 0, Qt::WindowFlags flags = 0);
	~Application();
	//void paintEvent(QPaintEvent*);
private slots:
	//void mouseMoveEvent(QMouseEvent *e);
public slots:
    void slotExit();
    void slotOpen();
    void slotOpenShape();
    void slotSave();
    void slotGo();
    void slotTest();
    void slotContinue();
    void slotStop();
    void slotUpdateParameters();
    void slotRotateSnake();
    void slotShiftSnake();
    void slotDefineInitSnake();
    void slotUpdatePreprocessingParameters();
    void slotCompute();
    void slotBrowseGroundTruth();
    void slotBrowseSnakeFile();
    void slotRebuild();

private:
    void display();
private:
    bool isInPhase1, isInPhase2, isInPhaseChangeSP;
	/*Private attributes*/
    double sigma,threshold;
    int width, height;
    QImage image, gradImageBeforeThresholding, gradImageAfterThresholding, flowImage;
    QRgb pixel;
    ReaderType::Pointer reader;

    QString groundTruthFile;

	//User interface parameters
	Ui::ApplicationClass ui;
	Scene* sceneMain;    
	Scene* sceneSnake;
    Scene* scene1;
    Scene* scene2;
    Scene* scene3;
    Scene* scene4;
    Scene* scenePostProcessing;

	QPen pen;
    //QPen penblue;

    ImageType::Pointer floatGradient;
    ImageType::Pointer floatFlow;

	bool stop;
	int loop;
	//
	QString theFileName;
	Snake s;
	PointList theShape;
    PointList groundTruth;
    PointList aSnake;

	double continuityCoeff;
	double curvatureCoeff;
	double imageCoeff;
	double flowCoeff;
	double balloonCoeff;
	double priorCoeff;

    double continuityBegin;
    double curvatureBegin;
    double imageBegin;
    double flowBegin;
    double balloonBegin;
    double priorBegin;

    double continuityEnd;
    double curvatureEnd;
    double imageEnd;
    double flowEnd;
    double balloonEnd;
    double priorEnd;

    double continuityStep;
    double curvatureStep;
    double imageStep;
    double flowStep;
    double balloonStep;
    double priorStep;


	int maxIterations;
	bool autoAdaptation;
	int autoAdaptLoop;
	int minLength;
	int maxLength;
	bool enableSPI;
    bool autoShifting;


    //Parameters for the 1st step of the auto-run scheme
    double continuityCoeff_2;
    double curvatureCoeff_2;
    double flowCoeff_2;
    int maxIterations_2;
    bool autoAdaptation_2;
    int autoAdaptLoop_2;
    int minLength_2;
    int maxLength_2;

	/*Private methods*/
	void getParameters();
    void updateMainView();
    void updateSnakeView();
    //void updateShapeView();
    void updateAllViews();
    void disallowChangingParameters();
    void allowChangingParameters();
	void run();

};

#endif // APPLICATION_H
