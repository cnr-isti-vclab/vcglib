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
Revision 1.43  2008/04/18 17:45:23  cignoni
fast return for compacting functions if no compacting is needed

Revision 1.42  2008/04/10 09:18:57  cignoni
moved Index function from append to the allocate

Revision 1.41  2008/04/03 22:47:10  cignoni
template the reorder functions on the vector types (for ocf)

Revision 1.40  2008/03/11 09:22:07  cignoni
Completed the garbage collecting functions CompactVertexVector and CompactFaceVector.

Revision 1.39  2007/12/11 20:18:55  cignoni
forgotten required std::

Revision 1.38  2007/12/11 11:35:50  cignoni
Added the CompactVertexVector garbage collecting function.

Revision 1.37  2007/10/16 16:46:53  cignoni
Added Allocator::DeleteFace and Allocator::DeleteVertex; Now the use of SetD() should be deprecated.

Revision 1.36  2007/01/11 10:24:25  cignoni
Added test in AddVertices  to do not update un-initalized vert references (for newly allocated faces)

Revision 1.35  2006/11/29 15:58:50  cignoni
Added check with the new end and avoided dangerous updating of already updated pointers

Revision 1.34  2006/11/28 22:34:28  cignoni
Added default constructor with null initialization to adjacency members.
AddFaces and AddVertices NEED to know if the topology is correctly computed to update it.

Revision 1.33  2006/11/13 13:12:27  ponchio
Removed a couple of useless assert.

Revision 1.32  2006/10/27 11:06:29  ganovelli
the calls to  HasFFAdjacency e HasVFAdjacency have been changed to override them for the optional attributes (see vcg/complex/trimesh/base.h)

Revision 1.31  2006/10/17 06:54:14  fiorin
Added #include <assert.h>

Revision 1.30  2006/10/02 09:31:47  ponchio
usual typename missing

Revision 1.29  2006/09/29 15:11:41  giec
Fixed a few bug.

Revision 1.28  2006/09/29 14:40:22  cignoni
Removed a useless, wrong version of AddFaces

Revision 1.27  2006/02/28 12:22:48  spinelli
fix bug end iterator++

Revision 1.26  2006/02/28 12:13:49  spinelli
fix bug end iterator++

Revision 1.25  2005/11/10 15:37:58  cignoni
Removed flags clearing (now it should be in the constructor of face and vertex)

Revision 1.24  2005/10/13 09:32:11  cignoni
Re-inserted the cFFp and cVFp access. If only the const version of the member function exists, the compiler will call it
when a non-const object invokes that function

Revision 1.23  2005/10/12 17:26:19  ponchio
cFFp doesn not exist -> FFp (there is the const version...)
same for cVFp.

Revision 1.22  2005/10/12 10:47:21  cignoni
Removed clearing of flags of added faces. Now the flag component has a constructor that clear it.
FF and VF adjacency are updated only if they are present and consistent (e.g. only if VFp(k) != 0 or FFp(k)!=0)

Revision 1.21  2005/07/01 11:22:00  cignoni
Corrected for the fourth time line a cast to Facetype at line 341.
Read the notes there before changing it again

Revision 1.20  2005/06/09 14:14:29  ganovelli
two warnings on type cast

Revision 1.19  2005/04/27 16:08:39  callieri
in addfaces, added casting for face* returned from vertex.VFp() [borland]

Revision 1.18  2005/03/23 13:22:57  turini
Wrong left parenthesis removed.

Revision 1.17  2005/03/23 11:29:49  ganovelli
cast int->iterator corrected

Revision 1.16  2005/02/19 10:43:11  ponchio
reverted tarini mod

Revision 1.15  2005/02/08 17:14:28  tarini
aggiunto un typecast a (FaceType*) per farlo compilare under Mingw comp

Revision 1.14  2004/10/14 15:08:04  pietroni
added #include <vector>

Revision 1.13  2004/09/07 07:36:32  fasano
Replaced some typename definitions

Revision 1.12  2004/08/25 15:15:26  ganovelli
minor changes to comply gcc compiler (typename's and stuff)

Revision 1.11  2004/08/07 17:38:00  pietroni
solved errors on AddFaces relative to VFp pointers of faces

Revision 1.10  2004/08/07 16:16:32  pietroni
corrected errors in AddFaces ( must be updated pointers to chain of faces of VFTopology)

Revision 1.9  2004/08/05 16:44:06  pietroni
added addafaces funtion with local values

Revision 1.8  2004/07/15 11:40:34  ganovelli
VFb to VFp

Revision 1.7  2004/05/11 14:12:13  ganovelli
general comment: minor modifications to compile with g++. Almost all
insertions of "typename" keyword and new line at the end of file

Revision 1.6  2004/05/10 13:24:21  cignoni
Updated names of adj functions and added ending newline

Revision 1.5  2004/04/21 14:06:10  ganovelli
#ifndef added

Revision 1.4  2004/03/31 14:43:56  cignoni
bug in update of VF adj

Revision 1.3  2004/03/12 15:25:29  cignoni
Corrected bug on the return of a wrong iterator

Revision 1.2  2004/03/03 15:35:52  cignoni
Yet another cr lf mismatch

Revision 1.1  2004/02/24 21:36:42  cignoni
grouped documentation, changed typenames and reflection mechanism

Revision 1.1  2004/02/19 13:11:06  cignoni
Initial commit


****************************************************************************/

#ifndef __VCGLIB_TRIALLOCATOR
#define __VCGLIB_TRIALLOCATOR

#include <vector>
#include <string>
#include <set>
#include <assert.h>
#include <vcg/complex/trimesh/base.h>
#include <vcg/container/simple_temporary_data.h>

namespace vcg {
	namespace tri {
		/** \addtogroup trimesh */
		
		template<class MeshType>
		size_t Index(MeshType &m, typename MeshType::VertexType &v) {return &v-&*m.vert.begin();}
		template<class MeshType>
		size_t Index(MeshType &m, typename MeshType::FaceType &f) {return &f-&*m.face.begin();}

		template<class MeshType>
		size_t Index(MeshType &m, const typename MeshType::VertexType *vp) {return vp-&*m.vert.begin();}
		template<class MeshType>
		size_t Index(MeshType &m, const typename MeshType::FaceType * fp) {return fp-&*m.face.begin();}

		// Placeholder. 
		// this one is called by the Compact and overridden by more specialized functions for OCF classes.
		// that manage also the additional types
		template <class face_type>
			void ReorderFace( std::vector<size_t> &newVertIndex, std::vector<face_type>  &vert)
		{}
		template <class vertex_type>
		void ReorderVert( std::vector<size_t> &newVertIndex, std::vector<vertex_type> &vert)
		{}
		
		template <class MeshType, class ATTR_CONT>
		void ReorderAttribute(ATTR_CONT &c,std::vector<size_t> & newVertIndex, MeshType &m){
			typename std::set<typename MeshType::PointerToAttribute>::iterator ai;	
				for(ai = c.begin(); ai != c.end(); ++ai)
					((typename MeshType::PointerToAttribute)(*ai)).Reorder(newVertIndex);
		}

		template <class MeshType, class ATTR_CONT>
		void ResizeAttribute(ATTR_CONT &c,const int & sz, MeshType &m){
			typename std::set<typename MeshType::PointerToAttribute>::iterator ai;	
				for(ai =c.begin(); ai != c.end(); ++ai)
					((typename MeshType::PointerToAttribute)(*ai)).Resize(m.vn);
		}

		/*@{*/
		/// Class to safely add vertexes and faces to a mesh updating all the involved pointers.
		/// It provides static memeber to add either vertex or faces to a trimesh.
		template <class AllocateMeshType>
		class Allocator
		{

		public:
			typedef AllocateMeshType MeshType;
			typedef typename MeshType::VertexType     VertexType;
			typedef typename MeshType::VertexPointer  VertexPointer;
			typedef typename MeshType::VertexIterator VertexIterator;
			typedef typename MeshType::VertContainer VertContainer;

			typedef typename MeshType::EdgeType     EdgeType;
			typedef typename MeshType::EdgePointer  EdgePointer;
			typedef typename MeshType::EdgeIterator EdgeIterator;
			typedef typename MeshType::EdgeContainer EdgeContainer;

			typedef typename MeshType::FaceType       FaceType;
			typedef typename MeshType::FacePointer    FacePointer;
			typedef typename MeshType::FaceIterator   FaceIterator;
			typedef typename MeshType::FaceContainer FaceContainer;
			typedef typename MeshType::PointerToAttribute PtrToAttr;
			typedef typename std::set<PtrToAttr>::iterator AttrIterator;
			typedef typename std::set<PtrToAttr>::const_iterator AttrConstIterator;
			typedef typename std::set<PtrToAttr >::iterator PAIte;

			/** This class is used when allocating new vertexes and faces to update 
			the pointers that can be changed when resizing the involved vectors of vertex or faces.
			It can also be used to prevent any update of the various mesh fields
			(e.g. in case you are building all the connections by hand as in a importer);
			*/ 
			template<class SimplexPointerType>
			class PointerUpdater
			{
			public:
				PointerUpdater(void) : newBase(0), oldBase(0), newEnd(0), oldEnd(0), preventUpdateFlag(false) { ; }
				void Clear(){newBase=oldBase=newEnd=oldEnd=0;};
				void Update(SimplexPointerType &vp)
				{
					if(vp>=newBase && vp<newEnd) return;
					assert(vp>=oldBase);
					assert(vp<oldEnd);
					vp=newBase+(vp-oldBase);
				}
				bool NeedUpdate() {if(oldBase && newBase!=oldBase && !preventUpdateFlag) return true; else return false;}

				SimplexPointerType newBase;
				SimplexPointerType oldBase;
				SimplexPointerType newEnd;
				SimplexPointerType oldEnd;
				bool preventUpdateFlag; /// when true no update is considered necessary.
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
				VertexIterator last;
				if(n == 0) return m.vert.end();
				pu.Clear();
				if(m.vert.empty()) pu.oldBase=0;  // if the vector is empty we cannot find the last valid element
        else {
          pu.oldBase=&*m.vert.begin(); 
  				pu.oldEnd=&m.vert.back()+1; 
        }

				m.vert.resize(m.vert.size()+n);
				m.vn+=n;

				typename std::set<PtrToAttr>::iterator ai;
				for(ai = m.vert_attr.begin(); ai != m.vert_attr.end(); ++ai)
					((PtrToAttr)(*ai)).Resize(m.vert.size());

				pu.newBase = &*m.vert.begin();
				pu.newEnd =  &m.vert.back()+1; 
		      if(pu.NeedUpdate())
					{
					FaceIterator fi;
					for (fi=m.face.begin(); fi!=m.face.end(); ++fi)
						if(!(*fi).IsD())
							for(int i=0; i < (*fi).VN(); ++i)
								if ((*fi).cV(i)!=0) pu.Update((*fi).V(i));
					EdgeIterator ei;
					for (ei=m.edge.begin(); ei!=m.edge.end(); ++ei)
						if(!(*ei).IsD())
						{
							if(HasEVAdjacency (m)) { pu.Update((*ei).V(0));pu.Update((*ei).V(1));}
							if(HasHEVAdjacency(m))   pu.Update((*ei).HEVp());
						}

						// e poiche' lo spazio e' cambiato si ricalcola anche last da zero  
				}
				unsigned int siz=(unsigned int)m.vert.size()-n;	
			
				last = m.vert.begin(); 
				advance(last,siz);
			
				return last;// deve restituire l'iteratore alla prima faccia aggiunta;
			}

			/** Function to add n vertices to the mesh.
			First wrapper, with no parameters
			*/
			static VertexIterator AddVertices(MeshType &m, int n)
			{
				PointerUpdater<VertexPointer> pu;
				return AddVertices(m, n,pu);
			}

			/** Function to add n vertices to the mesh.
			Second Wrapper, with a vector of vertex pointers to be updated.
			*/
			static VertexIterator AddVertices(MeshType &m, int n, std::vector<VertexPointer *> &local_vec)
			{
				PointerUpdater<VertexPointer> pu;
				VertexIterator v_ret =  AddVertices(m, n,pu);

				typename std::vector<VertexPointer *>::iterator vi;
            for(vi=local_vec.begin();vi!=local_vec.end();++vi)
               pu.Update(**vi);
				return v_ret;
			}

			/** Function to add n edges to the mesh. The second parameter hold a vector of 
			pointers to pointer to elements of the mesh that should be updated after a 
			possible vector realloc. 
			@param n number of edges to be added
			@param local_var vector of pointers to pointers to edges to be updated.
			return an iterator to the first element added
			*/
			static EdgeIterator AddEdges(MeshType &m,int n, PointerUpdater<EdgePointer> &pu)
			{
				EdgeIterator last;
				if(n == 0) return m.edge.end();
				pu.Clear();
				if(m.edge.empty()) pu.oldBase=0;  // if the vector is empty we cannot find the last valid element
				else {
					pu.oldBase=&*m.edge.begin(); 
					pu.oldEnd=&m.edge.back()+1; 
				}

				m.edge.resize(m.edge.size()+n);
				m.en+=n;

				typename std::set<typename MeshType::PtrToAttr>::iterator ai;
				for(ai = m.edge_attr.begin(); ai != m.edge_attr.end(); ++ai)
					((typename MeshType::PtrToAttr)(*ai)).Resize(m.edge.size());

				pu.newBase = &*m.edge.begin();
				pu.newEnd =  &m.edge.back()+1; 
		      if(pu.NeedUpdate())
					{
					int ii = 0;
					FaceIterator fi;
					for (fi=m.face.begin(); fi!=m.face.end(); ++fi){
						if(HasFHEAdjacency(m))
							pu.Update((*fi).FHEp());
						if(!(*fi).IsD())
							for(int i=0; i < (*fi).VN(); ++i)
								if ((*fi).cFEp(i)!=0) pu.Update((*fi).FEp(i));
					}

					VertexIterator vi;
					for (vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
						if(!(*vi).IsD())
								if ((*vi).cVEp()!=0) pu.Update((*vi).VEp());

					EdgeIterator ei = m.edge.begin();
					while(ii < m.en - n){// cycle on all the faces except the new ones
						if(!(*ei).IsD())
						{
							if(HasHENextAdjacency(m)) pu.Update((*ei).HENp());
							if(HasHEPrevAdjacency(m)) pu.Update((*ei).HEPp());
							if(HasHEOppAdjacency(m)) pu.Update((*ei).HEOp());
							++ii;
						}
						++ei;
					}

						// e poiche' lo spazio e' cambiato si ricalcola anche last da zero  
				}
				unsigned int siz=(unsigned int)m.edge.size()-n;	
			
				last = m.edge.begin(); 
				advance(last,siz);
			
				return last;// deve restituire l'iteratore alla prima faccia aggiunta;
			}

			/** Function to add n vertices to the mesh.
			First wrapper, with no parameters
			*/
			static EdgeIterator AddEdges(MeshType &m, int n)
			{
				PointerUpdater<EdgePointer> pu;
				return AddEdges(m, n,pu);
			}

      /** Function to add n vertices to the mesh.
			Second Wrapper, with a vector of vertex pointers to be updated.
			*/
			static EdgeIterator AddEdges(MeshType &m, int n, std::vector<EdgePointer*> &local_vec)
			{
				PointerUpdater<EdgePointer> pu;
				EdgeIterator v_ret =  AddEdges(m, n,pu);

				typename std::vector<EdgePointer *>::iterator ei;
            for(ei=local_vec.begin();ei!=local_vec.end();++ei)
               pu.Update(**ei);
				return v_ret;
			}


			/** Function to add n faces to the mesh.
			First wrapper, with no parameters
			*/
			static FaceIterator AddFaces(MeshType &m, int n)
			{
				PointerUpdater<FacePointer> pu;
				return AddFaces(m,n,pu);
			}

      /** Function to add n faces to the mesh.
			Second Wrapper, with a vector of face pointer to be updated.
			*/
			static FaceIterator AddFaces(MeshType &m, int n,std::vector<FacePointer *> &local_vec)
			{
				PointerUpdater<FacePointer> pu;
				FaceIterator f_ret= AddFaces(m,n,pu);

				typename std::vector<FacePointer *>::iterator fi;
            for(fi=local_vec.begin();fi!=local_vec.end();++fi)
               pu.Update(**fi);
				return f_ret;
			}

			/** Function to add n faces to the mesh. 
      This is the only full featured function that is able to manage correctly all the internal pointers of the mesh (ff and vf relations).
			NOTE: THIS FUNCTION ALSO UPDATE FN
			*/
			static FaceIterator AddFaces(MeshType &m, int n, PointerUpdater<FacePointer> &pu)
			{
				FaceIterator  last, fi;
				if(n == 0) return m.face.end();
				pu.Clear();
				if(m.face.empty()) {
					pu.oldBase=0;  // if the vector is empty we cannot find the last valid element
				}	else  {
					pu.oldBase=&*m.face.begin(); 
					pu.oldEnd=&m.face.back()+1; 
					last=m.face.end();
				} 

				m.face.resize(m.face.size()+n);
				m.fn+=n;


				typename std::set<PtrToAttr>::iterator ai;
				for(ai = m.face_attr.begin(); ai != m.face_attr.end(); ++ai)
					((PtrToAttr)(*ai)).Resize(m.face.size());

				pu.newBase = &*m.face.begin();
				pu.newEnd  = &m.face.back()+1;
				
				if(pu.NeedUpdate())
				{
					int ii = 0;
					FaceIterator fi = m.face.begin();
					while(ii<m.fn-n) // cycle on all the faces except the new ones
					{
						if(!(*fi).IsD())
						{
							if(HasFFAdjacency(m))
								for(int i  = 0; i < (*fi).VN(); ++i)
									if ((*fi).cFFp(i)!=0) pu.Update((*fi).FFp(i));
						
							if(HasVFAdjacency(m))
								for(int i = 0; i < (*fi).VN(); ++i)
									if ((*fi).cVFp(i)!=0) pu.Update((*fi).VFp(i));
						  ++ii;
						}
						++fi;
						}
						VertexIterator vi;
						for (vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
							if(!(*vi).IsD())
							{
								if(HasVFAdjacency(m))
									if ((*vi).cVFp()!=0)
										pu.Update((FaceType * &)(*vi).VFp());
								// Note the above cast is probably not useful if you have correctly defined 
								// your vertex type with the correct name of the facetype as a template argument; 
								//  pu.Update((FaceType*)(*vi).VFp()); compiles on old gcc and borland
								//  pu.Update((*vi).VFp());            compiles on .net and newer gcc
							}
						EdgeIterator ei;
						for (ei=m.edge.begin(); ei!=m.edge.end(); ++ei)
							if(!(*ei).IsD())
							{
								if(HasEFAdjacency(m))
									if ((*ei).cEFp()!=0)
										pu.Update((FaceType * &)(*ei).EFp());
								// Note the above cast is probably not useful if you have correctly defined 
								// your vertex type with the correct name of the facetype as a template argument; 
								//  pu.Update((FaceType*)(*vi).VFp()); compiles on old gcc and borland
								//  pu.Update((*vi).VFp());            compiles on .net and newer gcc
							}
	
							// e poiche' lo spazio e' cambiato si ricalcola anche last da zero  

				}
				unsigned int siz=(unsigned int)m.face.size()-n;	
				last = m.face.begin(); 
				advance(last,siz);
				return last;
			}

		/** Function to delete a face from the mesh. 
			NOTE: THIS FUNCTION ALSO UPDATE FN
		*/
		static void DeleteFace(MeshType &m, FaceType &f)
		{
			assert(!f.IsD());
			f.SetD();
			--m.fn;
		}

		/** Function to delete a vertex from the mesh. 
			NOTE: THIS FUNCTION ALSO UPDATE vn
		*/
		static void DeleteVertex(MeshType &m, VertexType &v)
		{
			assert(!v.IsD());
			v.SetD();
			--m.vn;
		}

		/** Function to delete an edge from the mesh. 
			NOTE: THIS FUNCTION ALSO UPDATE en
		*/
		static void DeleteEdge(MeshType &m, EdgeType &e)
		{
			assert(!e.IsD());
			e.SetD();
			--m.en;
		}
			
		/* 
		Function to compact all the vertices that have been deleted and put them to the end of the vector. 
		after this pass the isD test in the scanning of vertex vector, is no more strongly necessary.
		It should not be called when TemporaryData is active;
		*/
		
		static void CompactVertexVector( MeshType &m ) 
		{
			// If already compacted fast return please!
			if(m.vn==(int)m.vert.size()) return; 
			
			// newVertIndex [ <old_vert_position> ] gives you the new position of the vertex in the vector;
			std::vector<size_t> newVertIndex(m.vert.size(),std::numeric_limits<size_t>::max() );
			
			size_t pos=0;
			size_t i=0;
			
			for(i=0;i<m.vert.size();++i)
			{
				if(!m.vert[i].IsD())
				{
					if(pos!=i)
						m.vert[pos]=m.vert[i];
					newVertIndex[i]=pos;
					++pos;
				}
			}
			assert((int)pos==m.vn);
			
			// call a templated reordering function that manage any additional data internally stored by the vector 
			// for the default std::vector no work is needed (some work is typically needed for the OCF stuff) 
			ReorderVert<typename MeshType::VertexType>(newVertIndex,m.vert);
			
			// reorder the optional atttributes in m.vert_attr to reflect the changes 
			ReorderAttribute(m.vert_attr,newVertIndex,m);

			m.vert.resize(m.vn);

			// resize the optional atttributes in m.vert_attr to reflect the changes 
			ResizeAttribute(m.vert_attr,m.vn,m);

			FaceIterator fi;
			VertexPointer vbase=&m.vert[0];
			for(fi=m.face.begin();fi!=m.face.end();++fi)
				if(!(*fi).IsD())
					for(i=0;i<3;++i)
					{
						size_t oldIndex = (*fi).V(i) - vbase;
						assert(vbase <= (*fi).V(i) && oldIndex < newVertIndex.size());
						(*fi).V(i) = vbase+newVertIndex[oldIndex];
					}
				
		}

		/* 
		Function to compact all the vertices that have been deleted and put them to the end of the vector. 
		after this pass the isD test in the scanning of vertex vector, is no more strongly necessary.
		It should not be called when TemporaryData is active;
		*/
		
		static void CompactFaceVector( MeshType &m ) 
		{
		  // If already compacted fast return please!
			if(m.fn==(int)m.face.size()) return; 
			 
			// newFaceIndex [ <old_face_position> ] gives you the new position of the face in the vector;
			std::vector<size_t> newFaceIndex(m.face.size(),std::numeric_limits<size_t>::max() );
			
			size_t pos=0;
			size_t i=0;
			
			for(i=0;i<m.face.size();++i)
			{
				if(!m.face[i].IsD())
				{
					if(pos!=i)
						m.face[pos]=m.face[i];
					newFaceIndex[i]=pos;
					++pos;
				}
			}
			assert((int)pos==m.fn);
			
			// call a templated reordering function that manage any additional data internally stored by the vector 
			// for the default std::vector no work is needed (some work is typically needed for the OCF stuff) 
		  ReorderFace<typename MeshType::FaceType>(newFaceIndex,m.face);
					
			// reorder the optional atttributes in m.face_attr to reflect the changes 
			ReorderAttribute(m.face_attr,newFaceIndex,m);

			// Loop on the vertices to correct VF relation
			VertexIterator vi;
			FacePointer fbase=&m.face[0];
			for (vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
					if(!(*vi).IsD())
					{
						if(HasVFAdjacency(m))
									if ((*vi).cVFp()!=0)
									{
										size_t oldIndex = (*vi).cVFp() - fbase;
										assert(fbase <= (*vi).cVFp() && oldIndex < newFaceIndex.size());
										(*vi).VFp() = fbase+newFaceIndex[oldIndex];
									}
					}
			
			// Loop on the faces to correct VF and FF relations
			m.face.resize(m.fn);
			// resize the optional atttributes in m.face_attr to reflect the changes 
			ResizeAttribute(m.face_attr,m.vn,m);

			FaceIterator fi;
			for(fi=m.face.begin();fi!=m.face.end();++fi)
				if(!(*fi).IsD())
				{
					if(HasVFAdjacency(m))
								for(i=0;i<3;++i)
								if ((*fi).cVFp(i)!=0)
									{
										size_t oldIndex = (*fi).VFp(i) - fbase;
										assert(fbase <= (*fi).VFp(i) && oldIndex < newFaceIndex.size());
										(*fi).VFp(i) = fbase+newFaceIndex[oldIndex];
									}
					if(HasFFAdjacency(m))
										for(i=0;i<3;++i)
											if ((*fi).cFFp(i)!=0)
											{
												size_t oldIndex = (*fi).FFp(i) - fbase;
												assert(fbase <= (*fi).FFp(i) && oldIndex < newFaceIndex.size());
												(*fi).FFp(i) = fbase+newFaceIndex[oldIndex];
											}
				}
		}

public:

	/// Per Vertex Attributes
	template <class ATTR_TYPE>
	static 
	bool IsValidHandle( MeshType & m,  const typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE> & a){
		if(a._handle == NULL) return false;
		for(AttrIterator i = m.vert_attr.begin(); i!=m.vert_attr.end();++i)
			if ( (*i).n_attr == a.n_attr ) return true;
		return false;
	}

	template <class ATTR_TYPE> 
	static
	typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>
	 AddPerVertexAttribute( MeshType & m, std::string name){
		PAIte i;
		PtrToAttr h; 
		h._name = name;
		if(!name.empty()){
			i = m.vert_attr.find(h);
			assert(i ==m.vert_attr.end() );// an attribute with this name exists
		}
		h._sizeof = sizeof(ATTR_TYPE);
		h._padding = 0;
		h._handle = (void*) new SimpleTempData<VertContainer,ATTR_TYPE>(m.vert);
		m.attrn++;
		h.n_attr = m.attrn;
		std::pair < AttrIterator , bool> res =  m.vert_attr.insert(h);
		return typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>(res.first->_handle,res.first->n_attr );
	 }

	template <class ATTR_TYPE> 
	static
	typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>
	 AddPerVertexAttribute( MeshType & m){
		 return AddPerVertexAttribute<ATTR_TYPE>(m,std::string(""));
	 }

	template <class ATTR_TYPE> 
	static
		typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>
	 GetPerVertexAttribute( MeshType & m, const std::string & name){
		assert(!name.empty());
		PtrToAttr h1; h1._name = name;
		typename std::set<PtrToAttr > :: iterator i;

		i =m.vert_attr.find(h1);
		if(i!=m.vert_attr.end()){
			if(	(*i)._padding != 0 ){
					PtrToAttr attr = (*i);						// copy the PointerToAttribute
					m.vert_attr.erase(i);						// remove it from the set
					FixPaddedPerVertexAttribute<ATTR_TYPE>(m,attr);				
					std::pair<AttrIterator,bool> new_i = m.vert_attr.insert(attr);	// insert the modified PointerToAttribute
					assert(new_i.second);
					i = new_i.first;				
				}

			return typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>((*i)._handle,(*i).n_attr);

		}
			else
				return typename MeshType:: template PerVertexAttributeHandle<ATTR_TYPE>(NULL,0);
		
	}

	template <class ATTR_TYPE> 
	static
		void
	DeletePerVertexAttribute( MeshType & m,typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE> & h){
		typename std::set<PtrToAttr > ::iterator i;
		for( i = m.vert_attr.begin(); i !=  m.vert_attr.end(); ++i)
			if( (*i)._handle == h._handle ){
				delete ((SimpleTempData<VertContainer,ATTR_TYPE>*)(*i)._handle);
				m.vert_attr.erase(i); 
				return;}
			assert(0);
	}

	static
		void	DeletePerVertexAttribute( MeshType & m,  std::string name){
		AttrIterator i;
		PtrToAttr h1; h1._name = name;
		i = m.vert_attr.find(h1);
		assert(i!=m.vert_attr.end());
		delete ((SimpleTempDataBase<VertContainer>*)(*i)._handle);
		m.vert_attr.erase(i);
	}



	/// Per Edge Attributes
	template <class ATTR_TYPE>
	static 
	bool IsValidHandle( MeshType & m,  const typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE> & a){
		if(a._handle == NULL) return false;
		for(AttrIterator i = m.edge_attr.begin(); i!=m.edge_attr.end();++i)
			if ( (*i).n_attr == a.n_attr ) return true;
		return false;
	}

	template <class ATTR_TYPE> 
	static
	typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE>
	 AddPerEdgeAttribute( MeshType & m, std::string name){
		PAIte i;
		PtrToAttr h; 
		h._name = name;
		if(!name.empty()){
			i = m.edge_attr.find(h);
			assert(i ==m.edge_attr.end() );// an attribute with this name exists
		}
		h._handle = (void*) new SimpleTempData<EdgeContainer,ATTR_TYPE>(m.edge);
 		m.attrn++;
		h.n_attr = m.attrn;
		std::pair < AttrIterator , bool> res =  m.edge_attr.insert(h);
		return typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE>(res.first->_handle,res.first->n_attr);
	 }

	template <class ATTR_TYPE> 
	static
	typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>
	 AddPerEdgeAttribute( MeshType & m){
		 return AddPerEdgeAttribute<ATTR_TYPE>(m,std::string(""));
	 }
	
	template <class ATTR_TYPE> 
	static
		typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE>
	 GetPerEdgeAttribute( const MeshType & m, const std::string & name){
		assert(!name.empty());
		PtrToAttr h1; h1._name = name;
		typename std::set<PtrToAttr > ::const_iterator i;

		i =m.edge_attr.find(h1);
		if(i!=m.edge_attr.end())
				return typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>((*i)._handle,(*i).n_attr);
			else
				return typename MeshType:: template PerFaceAttributeHandle<ATTR_TYPE>(NULL,0);
	}

	template <class ATTR_TYPE> 
	static
		void
	DeletePerEdgeAttribute( MeshType & m,typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE> & h){
		typename std::set<PtrToAttr > ::iterator i;
		for( i = m.edge_attr.begin(); i !=  m.edge_attr.end(); ++i)
			if( (*i)._handle == h._handle ){
				delete ((SimpleTempData<FaceContainer,ATTR_TYPE>*)(*i)._handle);
				m.edge_attr.erase(i); 
				return;}
			assert(0);
	}

	static
		void	DeletePerEdgeAttribute( MeshType & m,  std::string name){
		AttrIterator i;
		PtrToAttr h1; h1._name = name;
		i = m.edge_attr.find(h1);
		assert(i!=m.edge_attr.end());
		delete ((SimpleTempDataBase<EdgeContainer>*)(*i)._handle);
		m.edge_attr.erase(i);
	}

	/// Per Face Attributes
	template <class ATTR_TYPE>
	static 
	bool IsValidHandle( MeshType & m,  const typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE> & a){
		if(a._handle == NULL) return false;
		for(AttrIterator i = m.face_attr.begin(); i!=m.face_attr.end();++i)
			if ( (*i).n_attr == a.n_attr ) return true;
		return false;
	}

	template <class ATTR_TYPE> 
	static
	typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>
	 AddPerFaceAttribute( MeshType & m, std::string name){
		PAIte i;
		PtrToAttr h; 
		h._name = name;
		if(!name.empty()){
			i = m.face_attr.find(h);
			assert(i ==m.face_attr.end() );// an attribute with this name exists
		}
		h._handle = (void*) new SimpleTempData<FaceContainer,ATTR_TYPE>(m.face);
		m.attrn++;
		h.n_attr = m.attrn;
		std::pair < AttrIterator , bool> res =  m.face_attr.insert(h);
		return typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>(res.first->_handle,res.first->n_attr);
	 }

	template <class ATTR_TYPE> 
	static
	typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>
	 AddPerFaceAttribute( MeshType & m){
		 return AddPerFaceAttribute<ATTR_TYPE>(m,std::string(""));
	 }
		
	template <class ATTR_TYPE> 
	static
		typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>
	 GetPerFaceAttribute( const MeshType & m, const std::string & name){
		assert(!name.empty());
		PtrToAttr h1; h1._name = name;
		typename std::set<PtrToAttr > ::const_iterator i;

		i =m.face_attr.find(h1);
		if(i!=m.face_attr.end())
				return typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>((*i)._handle,(*i).n_attr);
			else
				return typename MeshType:: template PerFaceAttributeHandle<ATTR_TYPE>(NULL,0);
	}

	template <class ATTR_TYPE> 
	static
		void
	DeletePerFaceAttribute( MeshType & m,typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE> & h){
		typename std::set<PtrToAttr > ::iterator i;
		for( i = m.face_attr.begin(); i !=  m.face_attr.end(); ++i)
			if( (*i)._handle == h._handle ){
				delete ((SimpleTempData<FaceContainer,ATTR_TYPE>*)(*i)._handle);
				m.face_attr.erase(i); 
				return;}
			assert(0);
	}

	static
		void	DeletePerFaceAttribute( MeshType & m,  std::string name){
		AttrIterator i;
		PtrToAttr h1; h1._name = name;
		i = m.face_attr.find(h1);
		assert(i!=m.face_attr.end());
		delete ((SimpleTempDataBase<FaceContainer>*)(*i)._handle);
		m.face_attr.erase(i);
	}

	/// Per Mesh Attributes
	template <class ATTR_TYPE>
	static 
	bool IsValidHandle( MeshType & m,  const typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE> & a){
		if(a._handle == NULL) return false;
		for(AttrIterator i = m.mesh_attr.begin(); i!=m.mesh_attr.end();++i)
			if ( (*i).n_attr == a.n_attr ) return true;
		return false;
	}

	template <class ATTR_TYPE> 
	static
	typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE>
	 AddPerMeshAttribute( MeshType & m, std::string name){
		PAIte i;
		PtrToAttr h; 
		h._name = name;
		if(!name.empty()){
			i = m.mesh_attr.find(h);
			assert(i ==m.mesh_attr.end() );// an attribute with this name exists
		}
		h._handle = (void*) new Attribute<ATTR_TYPE>();
		m.attrn++;
		h.n_attr = m.attrn;
		std::pair < AttrIterator , bool> res =  m.mesh_attr.insert(h);
		return typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE>(res.first->_handle,res.first->n_attr);
	 }
	
	template <class ATTR_TYPE> 
	static
		typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE>
	 GetPerMeshAttribute( const MeshType & m, const std::string & name){
		assert(!name.empty());
		PtrToAttr h1; h1._name = name;
		typename std::set<PtrToAttr > ::const_iterator i;

		i =m.mesh_attr.find(h1);
		if(i!=m.mesh_attr.end())
				return typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE>((*i)._handle,(*i).n_attr);
			else
				return typename MeshType:: template PerMeshAttributeHandle<ATTR_TYPE>(NULL,0);
	}

	template <class ATTR_TYPE> 
	static
		void
	DeletePerMeshAttribute( MeshType & m,typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE> & h){
		typename std::set<PtrToAttr > ::iterator i;
		for( i = m.mesh_attr.begin(); i !=  m.mesh_attr.end(); ++i)
			if( (*i)._handle == h._handle ){
				delete (( Attribute<ATTR_TYPE> *)(*i)._handle);
				m.mesh_attr.erase(i); 
				return;}
			assert(0);
	}

	static
		void	DeletePerMeshAttribute( MeshType & m,  std::string name){
		AttrIterator i;
		PtrToAttr h1; h1._name = name;
		i = m.mesh_attr.find(h1);
		assert(i!=m.mesh_attr.end());
		delete ((AttributeBase  *)(*i)._handle);
		m.mesh_attr.erase(i);
	}

	template <class ATTR_TYPE>
	static 
		void FixPaddedPerVertexAttribute ( MeshType & m,PtrToAttr & pa){

			// create the container of the right type
			SimpleTempData<VertContainer,ATTR_TYPE>* _handle =  new SimpleTempData<VertContainer,ATTR_TYPE>(m.vert);

			// copy the padded container in the new one
			_handle->Resize(m.vert.size());
			for(int i  = 0; i < m.vert.size(); ++i){
				ATTR_TYPE * dest = &(*_handle)[i];
				char * ptr = (char*)( ((SimpleTempDataBase<VertContainer> *)pa._handle)->DataBegin());
				memcpy((void*)dest ,  
				(void*) &(ptr[i * pa._sizeof + pa._padding]) ,sizeof(ATTR_TYPE));
			}

			// remove the padded container 
			delete ((SimpleTempDataBase<VertContainer>*) pa._handle);

			// update the pointer to data
			pa._handle = _handle;

			// zero the padding 
			pa._padding = 0;
	}
}; // end class
		

		/*@}*/
	} // End Namespace TriMesh
} // End Namespace vcg

#endif
