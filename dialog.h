#ifndef DIALOG_H
#define DIALOG_H

#include "ui_dialog.h"
#include "snake.h"
#include "scene.h"

#include <QtGui>


class ShapeDlg : public QDialog
{
	Q_OBJECT

public:	
	ShapeDlg(QWidget *parent = 0);
	~ShapeDlg();
	void setTheFileName(QString fileName){theFileName=fileName;}
	void Build(QString fileName);
	PointList getShape(){return shape;}
	
public slots:
	virtual void slotChangeText();
	virtual void slotBuild();
	virtual void slotCancel();
	virtual void slotUse();
    int exec();
private:
	Ui::Dialog dialog;
    Scene* scene;
	QString theFileName;
	PointList shape;
	double continuityCoeff;
	double curvatureCoeff;
	double imageCoeff;
	double flowCoeff;
	double balloonCoeff;
	double priorCoeff;
	bool autoAdaptation;
	int maxIterations;
	int autoAdaptLoop;
	int minLength;
	int maxLength;
	QPen pen;
};

#endif // DIALOG_H
