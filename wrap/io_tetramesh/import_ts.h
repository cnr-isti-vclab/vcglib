#ifndef __VCGLIB_IMPORTERTS
#define __VCGLIB_IMPORTERTS
#define NULL 0
#include <vcg/space/point3.h>

namespace vcg {
namespace tetra {
namespace io {

template <typename  MESHTYPE>
class ImporterTS{
  typedef typename MESHTYPE Tetramesh;
	typedef typename Tetramesh::VertexPointer VertexPointer;
	typedef typename Tetramesh::VertexType VertexType;
	typedef typename Tetramesh::TetraType FaceType;
	typedef typename Tetramesh::VertexIterator VertexIterator;
	typedef typename Tetramesh::TetraIterator FaceIterator;
	typedef typename Tetramesh::ScalarType ScalarType;
	typedef Point3<ScalarType> Point3x;

public:
static int Open( Tetramesh & m, const char * filename ){	
	int nvertex;
	int ntetra;
	float x;
	float y;
	float z;
	int tp0;
	int tp1;
	int tp2;
	int tp3;
	float mass;
	FILE *f;
	Tetramesh::VertexType p1;
	f = fopen(filename,"r");
	if(f == NULL ) 
		{
			printf( "The file was not opened\n" );
			return -1;
		}
   else
   {
		fscanf(f, "%i", &nvertex );
		fscanf(f, "%i", &ntetra );
		int j;
		for (j=0;j<nvertex;j++)
		{
			fscanf(f, "%f", &x );
			fscanf(f, "%f", &y );
			fscanf(f, "%f", &z );
		  //fscanf(f, "%f", &mass );
      p1.ClearFlags();
      p1.P()=Point3x(x, y,z );
			m.vert.push_back(p1);
		}
		m.tetra.reserve(ntetra);
    m.vert.reserve(nvertex);
		for (j=0;j<ntetra;j++)
		{
			fscanf(f, "%i", &tp0 );
			fscanf(f, "%i", &tp1 );
			fscanf(f, "%i", &tp2 );
			fscanf(f, "%i", &tp3 );
			
			Tetramesh::TetraType  newTetra;
			m.tetra.push_back(newTetra);
			m.tetra.back().Init(&m.vert[tp0],&m.vert[tp1],&m.vert[tp2],&m.vert[tp3]); 
		}
	 }
	 m.vn = nvertex;
	 m.tn = ntetra;

		return 0;
	 }
	};// end class
};// end of io
};// end of tri
};// end of vcg
#endif