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

#include <vcg/complex/local_optimization.h>
#include <vcg/simplex/face/topology.h>

namespace vcg
{
	namespace tri
	{
		/** \addtogroup trimesh */
		/* @{ */

		/*!
		*	This Class is specialization of LocalModification for the edge flip
		*	It wraps the atomic operation EdgeFlip to be used in a optimizatin routine.
		* Note that it has knowledge of the heap of the class LocalOptimization because
		* it is responsible of updating it after a flip has been performed
		*/
		template <class TRIMESH_TYPE, class MYTYPE>
		class TriEdgeFlip : public LocalOptimization< TRIMESH_TYPE >::LocModType
		{
		protected:
			typedef	typename TRIMESH_TYPE::FaceType				FaceType;
			typedef typename TRIMESH_TYPE::FacePointer		FacePointer;
			typedef typename TRIMESH_TYPE::FaceIterator		FaceIterator;
			typedef	typename TRIMESH_TYPE::VertexType			VertexType;
      typedef	typename TRIMESH_TYPE::VertexPointer	VertexPointer;
			typedef	typename TRIMESH_TYPE::ScalarType			ScalarType;
			typedef	typename TRIMESH_TYPE::CoordType			CoordType;
			typedef vcg::face::Pos<FaceType>							PosType;
			typedef typename LocalOptimization<TRIMESH_TYPE>::HeapElem HeapElem;
			typedef typename LocalOptimization<TRIMESH_TYPE>::HeapType HeapType;

			/*! 
			*	the pos of the flipping
			*/
			PosType _pos;

			/*!
			* priority in the heap
			*/
			ScalarType _priority;

			/*!
			* Mark for updating
			*/
			int _localMark;

			/*!
			*	mark for up_dating
			*/
			static int& GlobalMark()
			{ 
				static int im = 0; 
				return im;
			};


		public:
			/*!
			*	static data to gather statistical information 
			*	about the reasons of collapse failures
			*/
			struct FailStat 
			{
				static int &Volume()           {static int vol=0; return vol;}
				static int &LinkConditionFace(){static int lkf=0; return lkf;}
				static int &LinkConditionEdge(){static int lke=0; return lke;}
				static int &LinkConditionVert(){static int lkv=0; return lkv;}
				static int &OutOfDate()        {static int ofd=0; return ofd;}
				static int &Border()           {static int bor=0; return bor;}
				static void Init() 
				{
					Volume()           =0;
					LinkConditionFace()=0;
					LinkConditionEdge()=0;
					LinkConditionVert()=0;
					OutOfDate()        =0;
					Border()           =0;
				}
			};

			/*!
			*	Default constructor
			*/
			inline TriEdgeFlip()
			{};

			/*!
			*	Constructor with <I>pos</I> type
			*/
      inline TriEdgeFlip(PosType pos, int mark)
			{
				_pos = pos;
				_localMark = mark;
				_priority = ComputePriority();
			};

			/*!
			*/
			~TriEdgeFlip()
			{
			};

			/*!
			*	Return the LocalOptimization type
			*/
			ModifierType IsOfType()
			{ 
				return TriEdgeFlipOp;
			};

			/*!
			* Check if the pos is updated
			*/
			bool IsUpToDate()
			{
				VertexPointer v0 = _pos.V(0);
				VertexPointer v1 = _pos.V(1);

				if ((v0->IsD()) || (v1->IsD()) || _localMark < vcg::math::Min< int >( v0->IMark(), v1->IMark()))
				{
					++FailStat::OutOfDate();
					return false;
				}
				return true;
			};

			/*!
			* Check if this flipping operation can be performed
			*/
			bool IsFeasible()
			{
				return vcg::face::CheckFlipEdge(*_pos.f, _pos.z);
			};

			/*!
			*	Compute the priority of this optimization
			*/
			ScalarType ComputePriority()
			{
				FacePointer f = _pos.f;
				int z					= _pos.z;

				const ScalarType RatioThr = 20;
				const ScalarType AngleThr = (ScalarType)(M_PI/3600.0);
				
				int z1 = (*f).FFi(z);
				FacePointer f1 = (*f).FFp(z);
				VertexType *vx = (*f).FFp(z)->V2(z1); // ... ->V2((*f).FFi(z));
				CoordType &n0 = (*f).N();
				CoordType &n1 = (*f1).N();
				CoordType &n0d = (*f).FFp1(z)->N();
				CoordType &n0u = (*f).FFp2(z)->N();
				CoordType &n1d = (*f1).FFp1(z1)->N();
				CoordType &n1u = (*f1).FFp2(z1)->N();

				ScalarType a01 = AngleN((*f).N(),(*f1).N());
				ScalarType e01		= vcg::Distance((*f).V(z)->cP(), (*f).V1(z)->cP() );
				ScalarType e01f	= vcg::Distance((*f).V2(z)->cP(), vx->cP() );

				//Compute Edge Lenght Note that border edges are lenght 0
				ScalarType e0d = (vcg::face::IsBorder(*f, 1)) ? 0 : Distance((*f).V1(z)->cP(), (*f).V2(z)->cP() ) ;
				ScalarType e0u = (vcg::face::IsBorder(*f, 2)) ? 0 : Distance((*f).V(z)->cP() , (*f).V2(z)->cP() );
				ScalarType e1d = (vcg::face::IsBorder(*f1, (z1+1)%3)) ? 0 : Distance((*f).V1(z)->cP(), vx->cP() );
				ScalarType e1u = (vcg::face::IsBorder(*f1, (z1+2)%3)) ? 0 : Distance((*f).V(z)->cP() , vx->cP() ); 

				CoordType  n01u = ((vx->cP() - (*f).V(z)->cP())  ^ ((*f).V2(z)->cP() - (*f).V(z)->cP())).Normalize();
				CoordType  n01d = (((*f).V1(z)->cP() - vx->cP()) ^ ((*f).V2(z)->cP() - vx->cP())).Normalize();
				ScalarType a01f = vcg::AngleN(n01u,n01d); 
				ScalarType af0u = vcg::AngleN(n01u,n0u);
				ScalarType af0d = vcg::AngleN(n01d,n0d);
				ScalarType af1u = vcg::AngleN(n01u,n1u);
				ScalarType af1d = vcg::AngleN(n01d,n1d);

				e01 = e01f = e0d = e1d = e0u = e1u = 1; //pezza per pesare solo gli angoli!!!

				ScalarType OldCurvature = math::Max<ScalarType>(e01*a01, math::Max<ScalarType>(e0u*vcg::AngleN(n0,n0u), math::Max<ScalarType>( e0d*vcg::AngleN(n0,n0d), math::Max<ScalarType>(e1u*vcg::AngleN(n1,n1u) , e1d*vcg::AngleN(n1,n1d)))));
				ScalarType NewCurvature = math::Max<ScalarType>(e01f*a01f, math::Max<ScalarType>(e0u*AngleN(n01u,n0u), math::Max<ScalarType>( e0d*vcg::AngleN(n01d,n0d), math::Max<ScalarType>(e1u*vcg::AngleN(n01u,n1u), e1d*vcg::AngleN(n01d,n1d)))));

				_priority = (NewCurvature+AngleThr) - OldCurvature;
				return _priority;
			};

			/*!
			* Return the priority of this optimization
			*/
			ScalarType Priority() const
			{
				return _priority;
			};

			/*!
			* Execute the flipping of the edge
			*/
			void Execute(TRIMESH_TYPE &m)
			{
				vcg::face::FlipEdge(*_pos.f, _pos.z);
			};

			/*!
			*/
			const char* Info(TRIMESH_TYPE &m)
			{
				static char dump[60];
				sprintf(dump,"%i -> %i %g\n", _pos.V(0)-&m.vert[0], _pos.V(1)-&m.vert[0],-_priority);
				return dump;
			};

			/*!
			*/
			static void Init(TRIMESH_TYPE &mesh, HeapType &heap)
			{
				heap.clear();
				FaceIterator f_iter;
				for (f_iter = mesh.face.begin(); f_iter!=mesh.face.end(); ++f_iter)
				{
					if (! (*f_iter).IsD() )
					{
						for (unsigned int i=0; i<3; i++)
						{
							VertexPointer v0 = (*f_iter).V(i);
							VertexPointer v1 = (*f_iter).V((i+1)%3);
							if (v1-v0 > 0)
								heap.push_back( HeapElem( new MYTYPE(PosType(&*f_iter, i), mesh.IMark()) ) );
						} // endfor
					} // endif
				} //endfor
			};

			/*!
			*/
			void UpdateHeap(HeapType &heap)
			{
				GlobalMark()++;
				PosType pos(_pos.f, _pos.z);
				pos.FlipF();
				
				_pos.f->V1(_pos.z)->IMark() = GlobalMark();
				_pos.f->V2(_pos.z)->IMark() = GlobalMark();
				pos.f->V1(pos.z)->IMark() = GlobalMark();
				pos.f->V2(pos.z)->IMark() = GlobalMark();

				
				if (_pos.f->V2(_pos.z) - _pos.f->V1(_pos.z) > 0)
					heap.push_back( HeapElem( new MYTYPE( PosType(_pos.f, (_pos.z+1)%3, _pos.f->V1(_pos.z)), GlobalMark() ) ) );
				
				if (_pos.f->V(_pos.z) - _pos.f->V2(_pos.z) >0 )
					heap.push_back( HeapElem( new MYTYPE( PosType(_pos.f, (_pos.z+2)%3, _pos.f->V2(_pos.z)), GlobalMark() ) ) );

				if (pos.f->V2(pos.z) - pos.f->V1(pos.z) > 0)
					heap.push_back( HeapElem( new MYTYPE( PosType(pos.f, (pos.z+1)%3, pos.f->V1(pos.z)), GlobalMark() ) ) );

				if (pos.f->V(pos.z) - pos.f->V2(pos.z) > 0)
					heap.push_back( HeapElem( new MYTYPE( PosType(pos.f, (pos.z+2)%3, pos.f->V2(pos.z)), GlobalMark() ) ) );

				std::push_heap(heap.begin(),heap.end());
			};
		}; // end of TriEdgeFlip class

		/*! @} */
	}; // end of namespace tri
}; // end of namespace vcg