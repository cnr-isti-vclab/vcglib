#include <qapplication.h>
#include <qimage.h>
#include <segmentform.h>
#include <segmentator.h>
#include <qdir.h>
#include <qcolor.h>
#include <SimpleGLWidget.h>
#include <qtimer.h>

Segmentator *s;
QTimer *timer;

int main( int argc, char ** argv )
{
	s=new Segmentator();

	//s->LoadFromDir("./venacava/","prova.txt");//to chANGE

	//s->InitSegmentation(0.5,0.2,20,10.f);
	//
    QApplication a( argc, argv );

    SegmentForm w;
    w.show();

	//assign pointer to pricipal form
	w.simpleGLWidget1->w=&w;
	
	/*s=new Segmentator();*/

	//s->LoadFromDir("./venacava/","prova.txt");//to chANGE

	//s->InitSegmentation(0.5,0.2,20,10.f);
	

	timer = new QTimer(w.simpleGLWidget1 );
	QTimer::connect( timer, SIGNAL(timeout()), w.simpleGLWidget1, SLOT(Update()) );
    timer->start(0); // 2 seconds single-shot timer
	

    a.connect( &a, SIGNAL( lastWindowClosed() ), &a, SLOT( quit() ) );
    return a.exec();
}
