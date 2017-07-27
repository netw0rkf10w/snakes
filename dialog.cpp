#include "dialog.h"


ShapeDlg::ShapeDlg(QWidget* parent) : QDialog(parent)
{
	dialog.setupUi(this);

	dialog.continuityCoeff->setValidator(new QDoubleValidator(0,100,10,this));
	dialog.curvatureCoeff->setValidator(new QDoubleValidator(0,100,10,this));
	dialog.imageCoeff->setValidator(new QDoubleValidator(0,100,10,this));
	dialog.flowCoeff->setValidator(new QDoubleValidator(0,100,10,this));
    dialog.balloonCoeff->setValidator(new QDoubleValidator(0,100,10,this));

	dialog.maxIterations->setValidator(new QIntValidator(0,100000,this));
	dialog.autoAdaptLoop->setValidator(new QIntValidator(0,100,this));
	dialog.minLength->setValidator(new QIntValidator(1,100,this));
	dialog.maxLength->setValidator(new QIntValidator(1,100,this));

	//Set default values for the coefficents
	dialog.continuityCoeff->setText("1.0");
    dialog.curvatureCoeff->setText("0.0");
    dialog.imageCoeff->setText("0.0");
	dialog.flowCoeff->setText("1.0");
    dialog.balloonCoeff->setText("0.0");

	dialog.autoAdaptation->setChecked(true);	
	dialog.maxIterations->setText("300");
	dialog.autoAdaptLoop->setText("10");
	dialog.minLength->setText("15");
	dialog.maxLength->setText("30");

	pen=QPen(QBrush(Qt::red),3);

    scene = new Scene();
    scene->enableMouseIteration(true);
    scene->update();
    scene->setSceneRect(dialog.graphicsView->rect());

    dialog.graphicsView->setScene(scene);
    dialog.graphicsView->show();

	connect(dialog.rebuildButton,SIGNAL(clicked()), this, SLOT(slotBuild()));
	connect(dialog.cancelButton,SIGNAL(clicked()), this, SLOT(slotCancel()));
	connect(dialog.useButton,SIGNAL(clicked()), this, SLOT(slotUse()));

}

ShapeDlg::~ShapeDlg()
{
	shape.clear();
}

void ShapeDlg::slotUse()
{
	this->accept();	
}

void ShapeDlg::slotBuild()
{
	Build(theFileName);
}

void ShapeDlg::slotCancel()
{
	shape.clear();
	this->reject();
}


void ShapeDlg::Build(QString fileName)
{
    //theFileName=fileName;
	continuityCoeff=dialog.continuityCoeff->text().toDouble();
	curvatureCoeff=dialog.curvatureCoeff->text().toDouble();
	imageCoeff=dialog.imageCoeff->text().toDouble();
	flowCoeff=dialog.flowCoeff->text().toDouble();
    balloonCoeff=dialog.balloonCoeff->text().toDouble();

	maxIterations=dialog.maxIterations->text().toInt();
	autoAdaptation=dialog.autoAdaptation->isChecked();
	autoAdaptLoop=dialog.autoAdaptLoop->text().toInt();
	minLength=dialog.minLength->text().toInt();
	maxLength=dialog.maxLength->text().toInt();

	//ui.graphicsView_2->scale(0.5,0.5);
	//Snake s= getInitialSnake(fileName.toLatin1().constData());
    //const char* f = fileName.toUtf8().constData();
    Snake s(fileName.toUtf8().constData());
    if((int)scene->getPoints().size()>0)
    {
        //makeSureAllPointsAreInTheImage
        PointList::iterator i;
        PointList l = scene->getPoints();
        QImage im(fileName);
        for(i=l.begin();i!=l.end();i++)
        {
            if (i->getX() < 1) i->setX(1);
            if (i->getX() >= (im.width()-1)) i->setX(im.width()-2);
            if (i->getY() < 1) i->setY(1);
            if (i->getY() >= (im.height()-1)) i->setY(im.height()-2);
        }
        s.setSnake(l);
    }
    int loop=0;
    while(loop< maxIterations && s.step2(continuityCoeff, curvatureCoeff, imageCoeff, flowCoeff, balloonCoeff,0))
	{
		// auto adapt the number of points in the snake
        if (autoAdaptation && (loop % autoAdaptLoop)==0)
		{
			s.removeOverlappingPoints(minLength);
			s.addMissingPoints(maxLength);
		}
		scene->clear();
		scene->drawLine(s.getSnake());	
        dialog.size->setText(QString::number((int)s.getSnake().size()));
		qApp->processEvents();
		loop++;
	}
 
	// rebuild using spline interpolation
	//if (autoAdaptation) s.rebuild(maxLength);	
	shape = s.getSnake();
	//std::cout<<shape.size()<<std::endl;
}

int ShapeDlg::exec()
{
	if(!theFileName.isEmpty())
	{
        QImage im(theFileName);
        scene->setSceneRect(0,0,im.width(),im.height());
        scene->clear();
        dialog.graphicsView->setImageFile(theFileName);
        dialog.graphicsView->setImage(im);
        scene->update();
        Build(theFileName);
        return QDialog::exec();
	}
	return 0;
//    Q_D(QDialog);
//
//    if (d->eventLoop) {
//        qWarning("QDialog::exec: Recursive call detected");
//        return -1;
//    }
//
//    bool deleteOnClose = testAttribute(Qt::WA_DeleteOnClose);
//    setAttribute(Qt::WA_DeleteOnClose, false);
//
//    d->resetModalitySetByOpen();
//
//    bool wasShowModal = testAttribute(Qt::WA_ShowModal);
//    setAttribute(Qt::WA_ShowModal, true);
//    setResult(0);
//
////On Windows Mobile we create an empty menu to hide the current menu
//#ifdef Q_WS_WINCE_WM
//#ifndef QT_NO_MENUBAR
//    QMenuBar *menuBar = 0;
//    if (!findChild<QMenuBar *>())
//        menuBar = new QMenuBar(this);
//    if (qt_wince_is_smartphone()) {
//        QAction *doneAction = new QAction(tr("Done"), this);
//        menuBar->setDefaultAction(doneAction);
//        connect(doneAction, SIGNAL(triggered()), this, SLOT(_q_doneAction()));
//    }
//#endif //QT_NO_MENUBAR
//#endif //Q_WS_WINCE_WM
//
//    bool showSystemDialogFullScreen = false;
//#ifdef Q_OS_SYMBIAN
//    if (qobject_cast<QFileDialog *>(this) || qobject_cast<QFontDialog *>(this) ||
//        qobject_cast<QWizard *>(this)) {
//        showSystemDialogFullScreen = true;
//    }
//#endif // Q_OS_SYMBIAN
//
//    if (showSystemDialogFullScreen) {
//        setWindowFlags(windowFlags() | Qt::WindowSoftkeysVisibleHint);
//        setWindowState(Qt::WindowFullScreen);
//    }
//    show();
//
//#ifdef Q_WS_MAC
//    d->mac_nativeDialogModalHelp();
//#endif
//
//    QEventLoop eventLoop;
//    d->eventLoop = &eventLoop;
//    QPointer<QDialog> guard = this;
//    (void) eventLoop.exec(QEventLoop::DialogExec);
//    if (guard.isNull())
//        return QDialog::Rejected;
//    d->eventLoop = 0;
//
//    setAttribute(Qt::WA_ShowModal, wasShowModal);
//
//    int res = result();
//    if (deleteOnClose)
//        delete this;
//#ifdef Q_WS_WINCE_WM
//#ifndef QT_NO_MENUBAR
//    else if (menuBar)
//        delete menuBar;
//#endif //QT_NO_MENUBAR
//#endif //Q_WS_WINCE_WM
//    return res;
}

void ShapeDlg::slotChangeText()
{
	dialog.label->setText("hi hi");
}
