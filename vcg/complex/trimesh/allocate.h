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
Revision 1.1  2004/02/19 13:11:06  cignoni
Initial commit


****************************************************************************/

namespace vcg {
namespace tri {

template <class AllocateMeshType>
class Allocator
{
 
public:
typedef AllocateMeshType MeshType; 
typedef typename MeshType::VertexPointer VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::FacePointer FacePointer;
typedef typename MeshType::FaceIterator FaceIterator;

/* This class is used to */ 
template<class SimplexPointerType>
class PointerUpdater
{
public:
  void Clear();
  void Update(SimplexPointerType &vp);
  bool NeedUpdate();
  
  SimplexPointerType oldBase;
  SimplexPointerType newBase;
};


/** Function to add n vertices to the mesh. The second parameter hold a vector of 
	pointers to pointer to elements of the mesh that should be updated after a 
	possible vector realloc. 
	@param n Il numero di vertici che si vuole aggiungere alla mesh.
	@param local_var Vettore di variabili locali che rappresentano puntatori a vertici. 
	restituisce l'iteratore al primo elemento aggiunto.
*/
static VertexIterator AddVertices(MeshType &m,int n, PointerUpdater<VertexPointer> &pu)
{
  VertexIterator last=m.vert.end();
  pu.Clear();
	if(m.vert.empty()) pu.oldBase=0;  // if the vector is empty we cannot find the last valid element
	else               pu.oldBase=&*m.vert.begin(); 
    
	for(int i=0; i<n; ++i)
	{
    m.vert.push_back(MeshType::VertexType());
		m.vert.back().ClearFlags();
	}

	m.vn+=n;

	pu.newBase = &*m.vert.begin();
	if(pu.NeedUpdate())
		{
			FaceIterator fi;
			for (fi=m.face.begin(); fi!=m.face.end(); ++fi)
        if(!(*fi).IsD())
        {
          pu.Update((*fi).V(0));
          pu.Update((*fi).V(1));
          pu.Update((*fi).V(2));
        }
		
		// e poiche' lo spazio e' cambiato si ricalcola anche last da zero  
		if(last!=0)  
			{ 
				last = m.vert.begin(); 
				advance(last,siz+1);
			}
			else last=m.vert.begin(); 
		}
 
	return last;// deve restituire l'iteratore alla prima faccia aggiunta;
}

static VertexIterator AddVertices(MeshType &m, int n)
{
    PointerUpdater<VertexPointer> pu;
    return AddVertices(m, n,pu);
}

/** Function to add n faces to the mesh.
	@param n Il numero di facce che si vuole aggiungere alla mesh
*/
static FaceIterator AddFaces(MeshType &m, int n)
{
	PointerUpdater<FacePointer> pu;
  return AddFaces(m,n,pu);
}
/** Function to add n faces to the mesh. 
  NOTA: Aggiorna fn;
	The second parameter hold a vector of 
	pointers to pointer to elements of the mesh that should be updated after a 
	possible vector realloc.
	@param n Facce da aggiungere
	@param local_var Vettore di variabili locali che rappresentano puntatori a facce, occorre, 
	perche' questi valori siano consistenti, aggiornarli ogni qual volta venga eseguito un resize
	del contenitore delle facce.
*/
static FaceIterator AddFaces(MeshType &m, int n, PointerUpdater<FacePointer> &pu)
{
  FaceIterator last=m.vert.end();
  pu.Clear();
	if(m.face.empty()) pu.oldBase=0;  // if the vector is empty we cannot find the last valid element
	else               pu.oldBase=&*m.face.begin(); 
    
	unsigned int siz=0;	for(int i=0; i<n; ++i)
	{
    m.face.push_back(MeshType::FaceType());
		m.face.back().ClearFlags();
	}

	m.fn+=n;
	
  pu.newBase = &*m.face.begin();

	FaceIterator oldbegin, newbegin;
	oldbegin = face.begin();
	FaceIterator last=face.end();
	if(face.empty()) last=0;
	            else last--;

	unsigned int siz=0;
	MFTYPE dum;
	dum.Supervisor_Flags()=0;
	for(int i=0; i<n; ++i)
		face.push_back(dum);
	
	fn+=n;
	newbegin = face.begin();
	if(newbegin != oldbegin)// se e' cambiato lo spazio (vector abbastanza grande o lista)
	{
		if(MFTYPE::OBJ_TYPE & MFTYPE::OBJ_TYPE_A) 
		{
			FaceIterator f;
			for (f=face.begin(); f!=face.end(); ++f)
				for(int k=0; k<(*f).size(); ++k)if(!(*f).IsD())
					(*f).F(k) = (*f).F(k)-&*oldbegin+&*newbegin;
		}
		vector<face_base **>::iterator jit;
		for(jit=local_var.begin(); jit!=local_var.end(); ++jit)
			if((**jit) !=0 ) **jit = **jit-&*oldbegin+&*newbegin;
		
		// deve restituire l'iteratore alla prima faccia aggiunta;
		if(last!=0) 
		{ 
			last = face.begin(); 
			advance(last,siz+1);
		}
		else last=face.begin(); 
	}
	else // 
	{ assert(newbegin == oldbegin);
		// se non e'cambiato lo spazio (vector abbastanza grande o lista)
		if(last==0) last = face.begin(); // se il vettore era vuoto si restituisce begin
		           else advance(last,1); // altrimenti il primo dopo quello che era in precedenza l'ultimo valido.
	}

	return last;

}

}; // end class
} // End Namespace TriMesh
} // End Namespace vcg
