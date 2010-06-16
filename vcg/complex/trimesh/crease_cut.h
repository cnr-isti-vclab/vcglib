/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2008                                                \/)\/    *
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

#ifndef __VCG_CREASE_CUT
#define __VCG_CREASE_CUT
#include<vcg/simplex/face/jumping_pos.h>
#include<vcg/complex/trimesh/append.h>
#include<vcg/complex/trimesh/update/normal.h>
namespace vcg {
namespace tri {
	
/*
Crease Angle
Assume che:
la mesh abbia la topologia ff
la mesh non abbia complex (o se li aveva fossero stati detached)
Abbia le normali per faccia normalizzate!!


Prende una mesh e duplica tutti gli edge le cui normali nelle facce incidenti formano un angolo maggiore 
di <angle> (espresso in rad).
foreach face
 foreach unvisited vert vi
   scan the star of triangles around vi duplicating vi each time we encounter a crease angle.

the new (and old) vertexes are put in a std::vector that is swapped with the original one at the end.
 
Si tiene un vettore di interi 3 *fn che dice l'indice del vertice puntato da ogni faccia. 
quando si scandisce la stella intorno ad un vertici, per ogni wedge si scrive l'indice del vertice corrsipondente.
 
 
*/

template<class MESH_TYPE>
void CreaseCut(MESH_TYPE &m, float angleRad)
{
	typedef typename MESH_TYPE::CoordType				CoordType;
	typedef typename MESH_TYPE::ScalarType			ScalarType;
	typedef typename MESH_TYPE::VertexType			VertexType;
	typedef typename MESH_TYPE::VertexPointer		VertexPointer;
	typedef typename MESH_TYPE::VertexIterator	VertexIterator;
	typedef typename MESH_TYPE::FaceIterator		FaceIterator;
	typedef typename MESH_TYPE::FaceType				FaceType;
	typedef typename MESH_TYPE::FacePointer		FacePointer;
	
	tri::Allocator<MESH_TYPE>::CompactVertexVector(m);
	tri::Allocator<MESH_TYPE>::CompactFaceVector(m);
	
	tri::UpdateNormals<MESH_TYPE>::NormalizeFace(m);
	
	assert(m.HasFFTopology());
  typename MESH_TYPE::ScalarType cosangle=math::Cos(angleRad);

	tri::UpdateFlags<MESH_TYPE>::VertexClearV(m);
	std::vector<int> indVec(m.fn*3,-1);
	int newVertexCounter=m.vn;
	int creaseCounter=0;
	int startVn=m.vn;
	FaceIterator fi;
	//const FaceType * nextf;
	for(fi=m.face.begin();fi!=m.face.end();++fi) 
		 for(int j=0;j<3;++j) 
				if(!(*fi).V(j)->IsV() )  // foreach unvisited vertex we loop around it searching for creases.
					{
						(*fi).V(j)->SetV();
						
				    face::JumpingPos<FaceType> iPos(&*fi,j,(*fi).V(j));		
						size_t vertInd = Index(m,iPos.v);	 // 
						bool isBorderVertex = iPos.FindBorder();   // for border vertex we start from the border.
						face::JumpingPos<FaceType> startPos=iPos;
						if(!isBorderVertex)                        // for internal vertex we search the first crease and start from it  
								{
									do {
												ScalarType dotProd = iPos.FFlip()->cN().dot(iPos.f->N());
												iPos.NextFE();						
												if(dotProd<cosangle) break;							
									} while (startPos!=iPos);
									startPos=iPos;                       // the found crease become the new starting pos.
								}
						
						int locCreaseCounter=0;
						int curVertexCounter =vertInd;
						 
						do {																													// The real Loop
								ScalarType dotProd=iPos.FFlip()->cN().dot(iPos.f->N());			// test normal with the next face (fflip)
								size_t faceInd = Index(m,iPos.f);	
								indVec[faceInd*3+ iPos.VInd()] = curVertexCounter;				
								
								if(dotProd<cosangle)
											{ //qDebug("  Crease FOUND");
									    	++locCreaseCounter;									 
												curVertexCounter=newVertexCounter;
												newVertexCounter++;
											}														
								iPos.NextFE();						
						} while (startPos!=iPos);
						if(locCreaseCounter>0 && (!isBorderVertex) ) newVertexCounter--;
					}
	
	// A questo punto ho un vettore che mi direbbe per ogni faccia quale vertice devo mettere. Dopo che ho aggiunto i vertici necessari,
	// rifaccio il giro delle facce
	qDebug("adding %i vert for %i crease edges ",newVertexCounter-m.vn, creaseCounter);				
	tri::Allocator<MESH_TYPE>::AddVertices(m,newVertexCounter-m.vn);
 
	tri::UpdateFlags<MESH_TYPE>::VertexClearV(m);
	for(fi=m.face.begin();fi!=m.face.end();++fi)
		for(int j=0;j<3;++j) // foreach unvisited vertex
		{
			size_t faceInd = Index(m, *fi);	
			size_t vertInd = Index(m, (*fi).V(j));	
			int curVertexInd = indVec[faceInd*3+ j];
			assert(curVertexInd != -1);
			assert(curVertexInd < m.vn);
			if(curVertexInd < startVn) assert(size_t(curVertexInd) == vertInd); 
			if(curVertexInd >= startVn)
				{
					m.vert[curVertexInd].ImportData(*((*fi).V(j)));
					(*fi).V(j) = & m.vert[curVertexInd];
				}
		}		
		tri::UpdateNormals<MESH_TYPE>::PerVertexFromCurrentFaceNormal(m);
}

} // end namespace tri
} // end namespace vcg
#endif

/*
 face::Pos<FaceType> startPos=iPos;
 const FaceType * nextf;
 
 do{
	 size_t faceInd = Index(m,iPos.f);	
	 assert(iPos.f->V(iPos.VInd()) ==iPos.v);
	 qDebug("Face %2i:%i - v%2i -- %2i",faceInd,iPos.VInd(),Index(m,iPos.v),curVertexCounter);
	 indVec[faceInd*3+ iPos.VInd()] = curVertexCounter;
	 
	 nextf=iPos.FFlip();
	 ScalarType dotProd=nextf->cN() * iPos.f->N();
	 iPos.FlipF();
	 iPos.FlipE();
	 if(dotProd<cosangle)
	 {
		 ++creaseCounter;									 
		 if(iPos!=startPos)  
		 {
			 qDebug("Crease");
			 curVertexCounter=newVertexCounter;
			 newVertexCounter++;
		 }
	 }
 }
 while((iPos!=startPos) && (!iPos.IsBorder()) );
 qDebug("-- End Loop on v %i",Index(m,iPos.v));
				 }
	
	// A questo punto ho un vettore che mi direbbe per ogni faccia quale vertice devo mettere. Dopo che ho aggiunto i vertici necessari,
	// rifaccio il giro delle facce
	qDebug("adding %i vert for %i crease edges ",newVertexCounter-m.vn, creaseCounter);				
	tri::Allocator<MESH_TYPE>::AddVertices(m,newVertexCounter-m.vn);
 
	tri::UpdateFlags<MESH_TYPE>::VertexClearV(m);
	for(fi=m.face.begin();fi!=m.face.end();++fi)
 for(int j=0;j<3;++j) // foreach unvisited vertex
 {
			size_t faceInd = Index(m, *fi);	
			size_t vertInd = Index(m, (*fi).V(j));	
			int curVertexInd = indVec[faceInd*3+ j];
			qDebug("Setting Face %2i vert %2i : %2i -> %2i",faceInd,j,vertInd,curVertexInd);
			assert(curVertexInd != -1);
			assert(curVertexInd < m.vn);
			if(curVertexInd < startVn) assert(curVertexInd == vertInd); 
			if(curVertexInd >= startVn)
			{
				//m.vert[curVertexInd].ImportData(*((*fi).V(j)));
				tri::Append<CMeshO,CMeshO>::ImportVertex(m.vert[curVertexInd],*((*fi).V(j)));
				(*fi).V(j) = & m.vert[curVertexInd];
			}
 }		
 
*/
