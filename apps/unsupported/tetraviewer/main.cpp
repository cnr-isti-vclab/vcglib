#include <qapplication.h>
#include "mainframe.h"

#include <wrap\io_tetramesh\import_ts.h>
#include <wrap\io_tetramesh\import_ply.h>
#include <vcg\complex\tetramesh\update\topology.h>
#include <vcg\complex\tetramesh\update\bounding.h>
#include <vcg\complex\tetramesh\update\normal.h>


//#include <tetrastats.h>
#include "myglwidget.h"


MyTetraMesh TM;
MyTetraMesh *tm;
TetraStats<MyTetraMesh> Stats;
typedef vcg::tetra::io::ImporterTS<MyTetraMesh> ImpTS;
typedef vcg::tetra::UpdateTetraTopology<MyTetraMesh::VertexContainer ,MyTetraMesh::TetraContainer > UT;
typedef vcg::tetra::UpdateNormals<MyTetraMesh> UN;
typedef vcg::tetra::UpdateBounding<MyTetraMesh> UB;


//MainFrame *wp;

void openTetraMesh(const char* filename)
{
//opening the tetrahedral mesh
	QString path=QString(filename);
	QString ext =path.right(3);
	TM=MyTetraMesh();

	if (ext==".ts")
		ImpTS::Open(TM,filename);
	else
		vcg::tetra::io::ImporterPLY<MyTetraMesh> ::Open(TM,filename);
	
	UT::TTTopology(TM.vert,TM.tetra);
	UT::ClearVTTopology(TM.vert,TM.tetra);
	UT::VTTopology(TM.vert,TM.tetra);
	UN::PerVertex(TM);
	UB::Box(TM);
	tm=&TM;
	Stats.SetTetraMesh(tm);
	Stats.Update();
}



int main( int argc, char ** argv )
{
	
	tm=0;

    QApplication a( argc, argv );
    MainFrame w;
//	wp=&w;
    w.show();	

    a.connect( &a, SIGNAL( lastWindowClosed() ), &a, SLOT( quit() ) );

    return a.exec();
}

