#include <vcg/traced_vector.h>
#include <vcg/simple_temporary_data.h>


#include<vcg/space/point3.h>
#include<vcg/space/color4.h>
using namespace vcg;
#include <vcg/simplex/vertexplus/base.h>
#include<vcg/math/base.h>


class EdgeProto{};
class Vertex0   : public VertexSimp1< Vertex0,EdgeProto,vcg::vert::Coord3fOpt>{};

typedef TVector<Vertex0> TVect;
int main(){
int i;
	
// partial test for CORE DATA ***********
	// put some vertex in a vector
	TVector<Vertex0> c1;

	for( i = 0; i < 10; i++)
		c1.push_back(Vertex0());

	c1.EnableAttribute<Vertex0::CoordType>();	

	// put some more vertex inthe vector
	for( i= 0; i < 10; i++)
		c1.push_back(Vertex0());

	c1[2].P()=Point3f(1.0,2,3);
	Point3f p = c1[2].P();

	// drop  the attributes
	c1.DisableAttribute<Vertex0::CoordType>();	
// *****************************************


// USER DATA
// si puo' fare in 2 modi: Se il container e'di tipo TVector (traced vector)
// si puo'usare 
//(1)	c.NewTempData<TIPO_ATTRIBUTO>(); che  si occupa di resizare
// i container di dati temporanei (cioe' riflette tutte le push_back,reserve,resize 
// eventualemnte fatte sul vettore). NewTempoData resituisce una handle per accedere
// al dato.

// (2) si usa SimpleTempData che e' come uno tranne che non supporta automaticamente
//	le variazioni di dimensione
//	del vettore (si puo' fare a mano pero')...s

// partial test for USER DATA ***********
// modo (1)

	for( i = 0; i < 10; i++)
		c1.push_back(Vertex0());

	CATEntry<TVect,EntryCATMulti<TVect> >::Insert(c1); // questa riga sparira'

	TempData<TVect,int> handle =  c1.NewTempData<int>();
	handle[&c1[3]] = 71;
	// put some more vertex inthe vector
	for( i = 0; i < 10; i++)
		c1.push_back(Vertex0());
	int h = handle[&c1[3]];
	c1.DeleteTempData(handle);
// *****************************************


// partial test for USER DATA ***********
// modo (2)

	std::vector<Vertex0> c;
	for( i = 0; i < 10; i++)
		c.push_back(Vertex0());

	SimpleTempData<std::vector<Vertex0>,int> tmp(c);
	tmp.Start();
	tmp[&c[1]] = 22;
	int hh = tmp[&c[1]];
	tmp.Stop();
// **************************************


	}



