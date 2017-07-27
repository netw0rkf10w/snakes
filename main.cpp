
#include "application.h"
#include <QtGui>
//#include <qwt_plot.h>
//#include <qwt_plot_curve.h>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	Application w;
	w.setWindowFlags(Qt::Window);
	w.show();

	/*QwtPlot myPlot;
    myPlot.show();*/
	return a.exec();
}
