
#ifndef _TRIMESHTYPE_H
#define _TRIMESHTYPE_H

#include<vcg/simplex/vertex/vertex.h>
#include<vcg/simplex/face/with/afav.h>

using namespace std;
using namespace vcg;

class	CFace;
class	CEdge;
class	CVertex : public	Vertex<float,CEdge,CFace>{};
class	CFace :public FaceAFAV<CVertex,CEdge,CFace>{};

class	CMesh:	public tri::TriMesh< vector<CVertex>,	vector<CFace > >{};

#endif