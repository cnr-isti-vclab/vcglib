#include <vcg\simplex\vertex\with\atvn.h>
#include <vcg\simplex\tetrahedron\with\atav.h>
#include <vcg\complex\tetramesh\base.h>

class MyTetrahedron;
class DUMMYEDGETYPE;
class DUMMYFACETYPE;

class MyVertex: public vcg::VertexATVNd<DUMMYEDGETYPE,DUMMYFACETYPE,MyTetrahedron>{};

class	MyTetrahedron:	public vcg::TetraATAV<MyVertex,MyTetrahedron>{};

typedef vcg::tetra::Tetramesh< std::vector<MyVertex> ,std::vector<MyTetrahedron> > MyTetraMesh;
