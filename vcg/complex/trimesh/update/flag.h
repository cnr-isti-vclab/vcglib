/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.2  2004/03/10 00:46:10  cignoni
changed to the face::IsBorder() style

Revision 1.1  2004/03/05 10:59:24  cignoni
Changed name from plural to singular (normals->normal)

Revision 1.1  2004/03/04 00:37:56  cignoni
First working version!


****************************************************************************/
#ifndef __VCG_TRI_UPDATE_FLAGS
#define __VCG_TRI_UPDATE_FLAGS

namespace vcg {
namespace tri {
/** \addtogroup trimesh */
/*@{*/
/// Management, updating and computation of per-vertex and per-face flags (like border flags).
/// This class is used to compute or update some of the flags that can be stored in the mesh components. For now just Border flags (e.g. the flag that tells if a given edge of a face belong to a border of the mesh or not).

template <class UpdateMeshType>
class UpdateFlags
{

public:
typedef UpdateMeshType MeshType; 
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;


static void FaceBorderFromFF(MeshType &m)
{
	const int BORDERFLAG[3]={FaceType::BORDER0,FaceType::BORDER1,FaceType::BORDER2};
	FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)if(!(*fi).IsD())
		for(int j=0;j<3;++j)
		{
			//if(!(*fi).IsManifold(j)) (*fi).SetCF(j);
			//else 
      if(face::IsBorder(*fi,j)) (*fi).SetB(j);
					 else (*fi).ClearB(j);
		}
}

// versione minimale che non calcola i complex flag.
void FaceBorderFromNone()
{
	typedef PEdge<MCTYPE, MVTYPE, MFTYPE > MPEDGE;
	vector<MPEDGE> e;
	face_iterator pf;
	vector<MPEDGE>::iterator p;

	if( fn == 0 ) return;

	e.resize(fn*3);								// Alloco il vettore ausiliario
	p = e.begin();
	for(pf=face.begin();pf!=face.end();++pf)			// Lo riempio con i dati delle facce
		if( ! (*pf).IsDeleted() )
			for(int j=0;j<3;++j)
			{
				(*p).Set(&(*pf),j);
				(*pf).ClearB(j);
				++p;
			}
	assert(p==e.end());
	sort(e.begin(), e.end());							// Lo ordino per vertici
	
	vector<MPEDGE>::iterator pe,ps;
	for(ps = e.begin(), pe=e.begin(); pe<=e.end(); ++pe)	// Scansione vettore ausiliario
	{
		if( pe==e.end() || *pe != *ps )					// Trovo blocco di edge uguali
		{
			if(pe-ps==1) 	{	
					//++nborder;
					ps->f->SetB(ps->z);
			} else
			if(pe-ps!=2)  {  // Caso complex!!
				for(;ps!=pe;++ps)
						ps->f->SetB(ps->z); // Si settano border anche i complex.
			} 
			ps = pe;
//			++ne;										// Aggiorno il numero di edge
		}
	}
//	TRACE("found %i border (%i complex) on %i edges\n",nborder,ncomplex,ne);

}

	/// Bisogna carlcolare il border flag delle facce
void VertexBorderFromFace()
{
	vertex_iterator v;
	face_iterator f;

	for(v=vert.begin();v!=vert.end();++v)
		(*v).ClearB();

	for(f=face.begin();f!=face.end();++f)
	{
		for(int z=0;z<3;++z)
			if( (*f).IsB(z) )
			{
				(*f).V0(z)->SetB();
				(*f).V1(z)->SetB();
			}
	}
}


}; // end class

/*@}*/
}	// End namespace
}	// End namespace


#endif
