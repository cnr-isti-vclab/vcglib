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
Revision 1.8  2007/05/21 11:12:10  cignoni
Corrected gcc compiling issues

Revision 1.7  2006/01/19 14:18:08  spinelli
fix bug end iterator++

Revision 1.6  2005/05/30 09:43:41  spinelli
vertexIterator sostituito con VertexIterator

Revision 1.5  2005/05/17 21:14:56  ganovelli
some typecast (crs4)

Revision 1.4  2004/10/28 00:47:42  cignoni
Better Doxygen documentation

Revision 1.3  2004/09/20 08:37:47  cignoni
Better Doxygen docs

Revision 1.2  2004/05/10 14:41:25  ganovelli
name of adhacency function updated



****************************************************************************/

#ifndef __VCGLIB_EDGEALLOCATOR
#define __VCGLIB_EDGEALLOCATOR

namespace vcg {
	namespace edg {
		/** \addtogroup edgemesh */
		/*@{*/
		/// Class to safely add vertexes and faces to a mesh updating all the involved pointers.
		/// It provides static memeber to add either vertex or faces to a edgemesh.
		template <class AllocateMeshType>
		class Allocator
		{

		public:
			typedef AllocateMeshType MeshType; 
			typedef typename MeshType::VertexType     VertexType;
			typedef typename MeshType::VertexPointer  VertexPointer;
			typedef typename MeshType::VertexIterator VertexIterator;
			typedef typename MeshType::EdgeType       EdgeType;
			typedef typename MeshType::EdgePointer    EdgePointer;
			typedef typename MeshType::EdgeIterator   EdgeIterator;

			/** This class is used when allocating new vertexes and faces to update 
			the pointers that can be changed when resizing the involved vectors of vertex or faces.
			It can also be used to prevent any update of the various mesh fields
			(e.g. in case you are building all the connections by hand as in a importer);
			*/ 
			template<class SimplexPointerType>
			class PointerUpdater
			{
			public:
				void Clear(){newBase=oldBase=newEnd=oldEnd=0;preventUpdateFlag=false;};
				void Update(SimplexPointerType &vp)
				{
					vp=newBase+(vp-oldBase);
				}
				bool NeedUpdate() {if(newBase!=oldBase && !preventUpdateFlag) return true; else return false;}

				SimplexPointerType oldBase;
				SimplexPointerType newBase;
				SimplexPointerType newEnd;
				SimplexPointerType oldEnd;
				bool preventUpdateFlag; /// when true no update is considered necessary.
			};


			/** Function to safely add n vertices to a mesh. 

			@param m The mesh to be expanded
			@param n the number of vertexes to be added
			@param pu A PointerUpdater that stores the relocation that can be happened.
			*/
			static VertexIterator AddVertices(MeshType &m,int n, PointerUpdater<VertexPointer> &pu)
			{
				VertexIterator last=m.vert.end();
				pu.Clear();
				if(m.vert.empty()) pu.oldBase=0;  // if the vector is empty we cannot find the last valid element
				else               pu.oldBase=&*m.vert.begin(); 


				for(int i=0; i<n; ++i)
				{
					m.vert.push_back(VertexType());
					m.vert.back().ClearFlags();
				}

				m.vn+=n;

				pu.newBase = &*m.vert.begin();
				if(pu.NeedUpdate())
				{
					EdgeIterator ei;
					for (ei=m.edges.begin(); ei!=m.edges.end(); ++ei)
						if(!(*ei).IsD())
						{
							pu.Update((*ei).V(0));
							pu.Update((*ei).V(1));
						}
				}
				// e poiche' lo spazio e' cambiato si ricalcola anche last da zero  
				unsigned int siz=(unsigned int)m.vert.size()-n;	
				//if(last!=(VertexIterator)0)  
				//{ 
				last = m.vert.begin(); 
				advance(last,siz);
				//}
				//else last=m.vert.begin(); 


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
			static EdgeIterator AddEdges(MeshType &m, int n)
			{
				PointerUpdater<EdgePointer> pu;
				return AddEdges(m,n,pu);
			}
			/** Function to add n faces to the mesh. 
			NOTA: Aggiorna fn;
			*/
			static EdgeIterator AddEdges(MeshType &m, int n, PointerUpdater<EdgePointer> &pu)
			{
				EdgeIterator last=m.edges.end();
				pu.Clear();
				if(m.edges.empty()) {
					pu.oldBase=0;  // if the vector is empty we cannot find the last valid element
					//last=0;
				}	else  {
					pu.oldBase=&*m.edges.begin(); 
					last=m.edges.end();
				} 

				m.edges.resize(m.edges.size()+n);
				/*for(int i=0; i<n; ++i)
				{
				m.edges.push_back(MeshType::EdgeType());
				m.edges.back().ClearFlags();
				}*/

				m.en+=n;

				pu.newBase = &*m.edges.begin();

				if(pu.NeedUpdate())
				{
					EdgeIterator ei;
					for (ei=m.edges.begin(); ei!=m.edges.end(); ++ei)
						if(!(*ei).IsD())
						{
							if(EdgeType::HasEEAdjacency())
							{
								pu.Update((*ei).EEp(0));
								pu.Update((*ei).EEp(1));
							}
						}
						VertexIterator vi;
						for (vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
							if(!(*vi).IsD())
							{
								if(VertexType::HasVEAdjacency())
									pu.Update((*vi).VEp());
							}
				}						
				// e poiche' lo spazio e' cambiato si ricalcola anche last da zero  
				unsigned int siz=(unsigned int)m.edges.size()-(unsigned int)n;	
				//if(last!=(EdgeIterator)0)  
				//	{ 
				last = m.edges.begin(); 
				advance(last,siz);
				//	}
				//else last=m.edges.begin(); 


				return last;
			}

		}; // end class
		/*@}*/
	} // End Namespace TriMesh
} // End Namespace vcg

#endif
