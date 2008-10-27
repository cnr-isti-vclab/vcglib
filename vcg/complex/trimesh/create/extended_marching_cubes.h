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
/***************************************************************************/


#ifndef __VCG_EXTENDED_MARCHING_CUBES
#define __VCG_EXTENDED_MARCHING_CUBES

#include <float.h>
#include <assert.h>
#include <vector>
#include <vcg/math/base.h>
#include <vcg/math/matrix.h>
#include <vcg/math/lin_algebra.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/complex/trimesh/update/edges.h>
#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/update/topology.h>
#include <vcg/complex/trimesh/allocate.h>
#include <vcg/space/point3.h>
#include "emc_lookup_table.h"

namespace vcg
{
	namespace tri
	{
		// Doxygen documentation
		/** \addtogroup trimesh */
		/*@{*/

			/*
		* Cube description:
		*         3 ________ 2           _____2__     
		*         /|       /|         / |       /|    
		*       /  |     /  |      11/  3   10/  |    
		*   7 /_______ /    |      /__6_|__ /    |1   
		*    |     |  |6    |     |     |     |    
		*    |    0|__|_____|1    |     |__|_0|__|    
		*    |    /   |    /      7   8/   5    /     
		*    |  /     |  /        |  /     |  /9      
		*    |/_______|/          |/___4___|/         
		*   4          5                  
		*/

		//! This class implements the Extended Marching Cubes algorithm.
		/*!
		*	The implementation is enough generic: this class works only on one volume cell for each 
		*	call to <CODE>ProcessCell</CODE>. Using the field value at the cell corners, it adds to the
		*	mesh the triangles set approximating the surface that cross that cell. 
		*	@param TRIMESH_TYPE (Template parameter) the mesh type that will be constructed
		*	@param WALKER_TYPE	(Template parameter) the class that implements the traversal ordering of the volume.
		**/
		template<class TRIMESH_TYPE, class WALKER_TYPE>
		class ExtendedMarchingCubes
		{
		public:
#if defined(__GNUC__)
			typedef unsigned int				size_t;
#else
#ifdef			_WIN64
			typedef unsigned __int64    size_t;
#else
			typedef _W64 unsigned int   size_t;
#endif
#endif 
			typedef typename vcg::tri::Allocator< TRIMESH_TYPE > AllocatorType;
			typedef typename TRIMESH_TYPE::ScalarType			ScalarType;
			typedef typename TRIMESH_TYPE::VertexType			VertexType;
			typedef typename TRIMESH_TYPE::VertexPointer	VertexPointer;
			typedef typename TRIMESH_TYPE::VertexIterator	VertexIterator;
			typedef typename TRIMESH_TYPE::FaceType				FaceType;
			typedef typename TRIMESH_TYPE::FacePointer		FacePointer;
			typedef typename TRIMESH_TYPE::FaceIterator		FaceIterator;
			typedef typename TRIMESH_TYPE::CoordType			CoordType;
			typedef typename TRIMESH_TYPE::CoordType*			CoordPointer;

			typedef struct
			{
				size_t face, edge;
			} LightEdge;

			/*!
			*	Constructor
			*	\param mesh		The mesh that will be constructed
			*	\param volume	The volume describing the field
			*	\param walker	The class implementing the traversal policy
			*	\param angle	The feature detection threshold misuring the sharpness of a feature(default is 30 degree)
			*/
			ExtendedMarchingCubes(TRIMESH_TYPE &mesh, WALKER_TYPE &walker, ScalarType angle=30)
			{
				_mesh		= &mesh;
				_walker = &walker;
				_featureAngle = vcg::math::ToRad(angle);
				_initialized = _finalized = false;
			};

			/*!
			*	Execute the initialiazation. 
			*	This method must be executed before the first call to <CODE>ApplyEMC</CODE>
			*/
			void Initialize()
			{
				assert(!_initialized && !_finalized);
				_featureFlag	= VertexType::NewBitFlag();
				_initialized	= true;
			};

			/*!
			*	
			*	This method must be executed after the last call to <CODE>ApplyEMC</CODE>
			*/
			void Finalize()
			{
				assert(_initialized && !_finalized);
				FlipEdges();

				VertexIterator v_iter = _mesh->vert.begin();
				VertexIterator v_end	= _mesh->vert.end();
				for ( ; v_iter!=v_end; v_iter++)
					v_iter->ClearUserBit( _featureFlag );
				VertexType::DeleteBitFlag( _featureFlag	);
				_featureFlag = 0;
				_mesh = NULL;
				_walker = NULL;
				_finalized = true;
			};

			/*!
			* Apply the <I>extended marching cubes</I> algorithm to the volume cell identified by the two points <CODE>min</CODE> and <CODE>max</CODE>.
			*	All the three coordinates of the first point must be smaller than the respectives three coordinatas of the second point.
			*	\param min	the first point
			*	\param max	the second point
			*/
		void ProcessCell(const vcg::Point3i &min, const vcg::Point3i &max)
			{
				assert(_initialized && !_finalized);
				assert(min[0]<max[0] && min[1]<max[1] && min[2]<max[2]);
				_corners[0].X()=min.X();		_corners[0].Y()=min.Y();		_corners[0].Z()=min.Z();
				_corners[1].X()=max.X();		_corners[1].Y()=min.Y();		_corners[1].Z()=min.Z();
				_corners[2].X()=max.X();		_corners[2].Y()=max.Y();		_corners[2].Z()=min.Z();
				_corners[3].X()=min.X();		_corners[3].Y()=max.Y();		_corners[3].Z()=min.Z();
				_corners[4].X()=min.X();		_corners[4].Y()=min.Y();		_corners[4].Z()=max.Z();
				_corners[5].X()=max.X();		_corners[5].Y()=min.Y();		_corners[5].Z()=max.Z();
				_corners[6].X()=max.X();		_corners[6].Y()=max.Y();		_corners[6].Z()=max.Z();
				_corners[7].X()=min.X();		_corners[7].Y()=max.Y();		_corners[7].Z()=max.Z();

				unsigned char cubetype = 0;
				if ((_field[0] = _walker->V(_corners[0].X(), _corners[0].Y(), _corners[0].Z())) >= 0) cubetype+=  1;
				if ((_field[1] = _walker->V(_corners[1].X(), _corners[1].Y(), _corners[1].Z())) >= 0) cubetype+=  2;
				if ((_field[2] = _walker->V(_corners[2].X(), _corners[2].Y(), _corners[2].Z())) >= 0) cubetype+=  4;
				if ((_field[3] = _walker->V(_corners[3].X(), _corners[3].Y(), _corners[3].Z())) >= 0) cubetype+=  8;
				if ((_field[4] = _walker->V(_corners[4].X(), _corners[4].Y(), _corners[4].Z())) >= 0) cubetype+= 16;
				if ((_field[5] = _walker->V(_corners[5].X(), _corners[5].Y(), _corners[5].Z())) >= 0) cubetype+= 32;
				if ((_field[6] = _walker->V(_corners[6].X(), _corners[6].Y(), _corners[6].Z())) >= 0) cubetype+= 64;
				if ((_field[7] = _walker->V(_corners[7].X(), _corners[7].Y(), _corners[7].Z())) >= 0) cubetype+=128;

				if (cubetype==0 || cubetype==255)
					return;

				size_t vertices_idx[12];
				memset(vertices_idx, -1, 12*sizeof(size_t));
				int code = EMCLookUpTable::EdgeTable(cubetype);
				VertexPointer vp = NULL;
				if (   1&code ) { _walker->GetXIntercept(_corners[0], _corners[1], vp); vertices_idx[ 0] = vp - &_mesh->vert[0]; }
				if (   2&code ) { _walker->GetYIntercept(_corners[1], _corners[2], vp); vertices_idx[ 1] = vp - &_mesh->vert[0]; }
				if (   4&code ) { _walker->GetXIntercept(_corners[3], _corners[2], vp); vertices_idx[ 2] = vp - &_mesh->vert[0]; }
				if (   8&code ) { _walker->GetYIntercept(_corners[0], _corners[3], vp); vertices_idx[ 3] = vp - &_mesh->vert[0]; }
				if (  16&code ) { _walker->GetXIntercept(_corners[4], _corners[5], vp); vertices_idx[ 4] = vp - &_mesh->vert[0]; }
				if (  32&code ) { _walker->GetYIntercept(_corners[5], _corners[6], vp); vertices_idx[ 5] = vp - &_mesh->vert[0]; }
				if (  64&code ) { _walker->GetXIntercept(_corners[7], _corners[6], vp); vertices_idx[ 6] = vp - &_mesh->vert[0]; }
				if ( 128&code ) { _walker->GetYIntercept(_corners[4], _corners[7], vp); vertices_idx[ 7] = vp - &_mesh->vert[0]; }
				if ( 256&code ) { _walker->GetZIntercept(_corners[0], _corners[4], vp); vertices_idx[ 8] = vp - &_mesh->vert[0]; }
				if ( 512&code ) { _walker->GetZIntercept(_corners[1], _corners[5], vp); vertices_idx[ 9] = vp - &_mesh->vert[0]; }
				if (1024&code ) { _walker->GetZIntercept(_corners[2], _corners[6], vp); vertices_idx[10] = vp - &_mesh->vert[0]; }
				if (2048&code ) { _walker->GetZIntercept(_corners[3], _corners[7], vp); vertices_idx[11] = vp - &_mesh->vert[0]; }

				int m, n, vertices_num;
				int components = EMCLookUpTable::TriTable(cubetype, 1)[0]; //unsigned int components =  triTable[cubetype][1][0];
				int	*indices   = &EMCLookUpTable::TriTable(cubetype, 1)[components+1]; //int					 *indices   = &EMCLookUpTable::TriTable(cubetype, 1, components+1);

				std::vector< size_t > vertices_list;
				for (m=1; m<=components; m++)
				{
					// current sheet contains vertices_num vertices
					vertices_num = EMCLookUpTable::TriTable(cubetype, 1)[m]; //vertices_num = triTable[cubetype][1][m];

					// collect vertices
					vertices_list.clear();
					for (n=0; n<vertices_num; ++n)
						vertices_list.push_back( vertices_idx[ indices[n] ] );

					VertexPointer feature = FindFeature( vertices_list );
					if (feature != NULL) // i.e. is a valid vertex
					{
						// feature -> create triangle fan around feature vertex
						size_t feature_idx = feature - &_mesh->vert[0];
						size_t face_idx			= _mesh->face.size();
						vertices_list.push_back( vertices_list[0] );
						AllocatorType::AddFaces(*_mesh, (int) vertices_num);
						for (int j=0; j<vertices_num; ++j, face_idx++)
						{
							_mesh->face[face_idx].V(0) = &_mesh->vert[ vertices_list[j  ] ];
							_mesh->face[face_idx].V(1) = &_mesh->vert[ vertices_list[j+1] ];
							_mesh->face[face_idx].V(2) = &_mesh->vert[ feature_idx			 ];
						}
					}
					else
					{
						// no feature -> old marching cubes triangle table
						for (int j=0; EMCLookUpTable::PolyTable(vertices_num, j) != -1; j+=3) //for (int j=0; polyTable[vertices_num][j] != -1; j+=3)
						{
							size_t face_idx			= _mesh->face.size();
							AllocatorType::AddFaces(*_mesh, 1);
							//_mesh->face[ face_idx].V(0) = &_mesh->vert[ vertices_idx[ indices[ polyTable[vertices_num][j  ] ] ] ];
							//_mesh->face[ face_idx].V(1) = &_mesh->vert[ vertices_idx[ indices[ polyTable[vertices_num][j+1] ] ] ];
							//_mesh->face[ face_idx].V(2) = &_mesh->vert[ vertices_idx[ indices[ polyTable[vertices_num][j+2] ] ] ];
							_mesh->face[ face_idx].V(0) = &_mesh->vert[ vertices_idx[ indices[ EMCLookUpTable::PolyTable(vertices_num, j  ) ] ] ];
							_mesh->face[ face_idx].V(1) = &_mesh->vert[ vertices_idx[ indices[ EMCLookUpTable::PolyTable(vertices_num, j+1) ] ] ];
							_mesh->face[ face_idx].V(2) = &_mesh->vert[ vertices_idx[ indices[ EMCLookUpTable::PolyTable(vertices_num, j+2) ] ] ];
						}
					}
					indices += vertices_num;

				}
			}; // end of ApplyEMC

		private:
			/*!
			*/
			WALKER_TYPE		*_walker;
			/*!
			*/
			TRIMESH_TYPE	*_mesh;
			/*!
			*/
			bool	_initialized;;
			/*!
			*/
			bool	_finalized;
			/*!
			*	The feature detection threshold misuring the sharpness of a feature
			*/
			ScalarType _featureAngle;
			/*!
			*	The flag used for marking the feature vertices.
			*/
			int				_featureFlag;
			/*!
			*	Array of the 8 corners of the volume cell being processed
			*/
			vcg::Point3i _corners[8];
			/*!
			*	The field value at the cell corners
			*/
			ScalarType	 _field[8];


			/*!
			*	Tests if the surface patch crossing the current cell contains a sharp feature
			*	\param vertices_idx	The list of vertex indices intersecting the edges of the current cell
			*	\return							The pointer to the new Vertex if a feature is detected; NULL otherwise.
			*/
			VertexPointer FindFeature(const std::vector<size_t> &vertices_idx)
			{
				unsigned int i, j, rank;
				size_t vertices_num = (size_t) vertices_idx.size();

				CoordType *points		= new CoordType[ vertices_num ];
				CoordType *normals	= new CoordType[ vertices_num ];

				for (i=0; i<vertices_num; i++)
				{
					points[i]		= _mesh->vert[ vertices_idx[i] ].P();
					normals[i]	= _mesh->vert[ vertices_idx[i] ].N();
				}

				// move barycenter of points into (0, 0, 0)
				CoordType center((ScalarType) 0.0, (ScalarType) 0.0, (ScalarType) 0.0);
				for (i=0; i<vertices_num; ++i)  
					center += points[i];
				center /= (ScalarType) vertices_num;
				for (i=0; i<vertices_num; ++i)  
					points[i] -= center;

				// normal angle criterion
				double c, minC, maxC; 
				CoordType axis;
				for (minC=1.0, i=0; i<vertices_num-1; ++i)
				{
					for (j=i+1; j<vertices_num; ++j)
					{
						c = normals[i]*normals[j];
						if (c < minC)
						{
							minC = c;
							axis = normals[i] ^ normals[j];
						}
					}
				} //end for (minC=1.0, i=0; i<vertNumber; ++i)

				if (minC > cos(_featureAngle)) 
					return NULL;  // invalid vertex

				// ok, we have a feature: is it edge or corner, i.e. rank 2 or 3 ?
				axis.Normalize();
				for (minC=1.0, maxC=-1.0, i=0; i<vertices_num; ++i)
				{
					c = axis * normals[i];
					if (c < minC)  minC = c;
					if (c > maxC)  maxC = c;
				}
				c = vcg::math::Max< double >(fabs(minC), fabs(maxC));
				c = sqrt(1.0-c*c);
				rank = (c > cos(_featureAngle) ? 2 : 3);

				// setup linear system (find intersection of tangent planes)
				vcg::ndim::Matrix<double>  A((unsigned int) vertices_num, 3);
				double *b = new double[ vertices_num ];
				for (i=0; i<vertices_num; ++i)
				{
					A[i][0] =  normals[i][0];
					A[i][1] =  normals[i][1];
					A[i][2] =  normals[i][2];
					b[i]		=  (points[i] * normals[i]);
				}

				// SVD of matrix A
				vcg::ndim::Matrix<double>  V(3, 3);
				double *w = new double[vertices_num];
				vcg::SingularValueDecomposition< typename vcg::ndim::Matrix<double> > (A, w, V, LeaveUnsorted, 100);

				// rank == 2 -> suppress smallest singular value
				if (rank == 2)
				{
					double        smin   = DBL_MAX; // the max value, as defined in <float.h>
					unsigned int  sminid = 0;
					unsigned int  srank  = vcg::math::Min< unsigned int >(vertices_num, 3u);

					for (i=0; i<srank; ++i)
					{
						if (w[i] < smin)
						{
							smin   = w[i];
							sminid = i;
						}
					}
					w[sminid] = 0.0;
				}

				// SVD backsubstitution -> least squares, least norm solution x
				double *x = new double[3];
				vcg::SingularValueBacksubstitution< vcg::ndim::Matrix<double> >(A, w, V, x, b);

				// transform x to world coords
				CoordType point((ScalarType) x[0], (ScalarType) x[1], (ScalarType) x[2]);
				point += center;

				// insert the feature-point 
				VertexPointer mean_point = &*AllocatorType::AddVertices( *_mesh, 1);
				mean_point->SetUserBit(_featureFlag);
				mean_point->P() = point;
				mean_point->N().SetZero();
				delete []x;
				delete []points;
				delete []normals;
				return mean_point;
			} // end of FindFeature

			/*!
			*	Postprocessing step performed during the finalization tha flip some of the mesh edges.
			*	The flipping criterion is quite simple: each edge is flipped if it will connect two 
			*	feature samples after the flip.
			*/
			void FlipEdges()
			{
				size_t i;
				std::vector< LightEdge > edges;
				FaceIterator f_iter = _mesh->face.begin();
				FaceIterator f_end  = _mesh->face.end();
				for (i=0; f_iter!=f_end; f_iter++, i++)
				{
					if (f_iter->V(1) > f_iter->V(0))
					{
						LightEdge le;
						le.face = i;
						le.edge = 0;
						edges.push_back( le );
					}
					if (f_iter->V(2) > f_iter->V(1)) 
					{
						LightEdge le;
						le.face = i;
						le.edge = 1;
						edges.push_back( LightEdge(le));
					}
					if (f_iter->V(0) > f_iter->V(2)) 
					{
						LightEdge le;
						le.face = i;
						le.edge = 2;
						edges.push_back( le );
					}
				}
				vcg::tri::UpdateTopology< TRIMESH_TYPE >::VertexFace( *_mesh );
				vcg::tri::UpdateTopology< TRIMESH_TYPE >::FaceFace( *_mesh );

				typename std::vector< LightEdge >::iterator e_it		= edges.begin();
				typename std::vector< LightEdge >::iterator e_end	= edges.end();

				FacePointer g, f;
				int			w, z;
				for( ; e_it!=e_end; e_it++)
				{
					f = &_mesh->face[e_it->face];
					z = (int) e_it->edge;

					if (vcg::face::CheckFlipEdge< FaceType >(*f, z))
					{
						VertexPointer v0, v1, v2, v3;
						v0 = f->V(z);
						v1 = f->V1(z);
						v2 = f->V2(z);
						g = f->FFp(z);
						w = f->FFi(z);
						v3 = g->V2(w);
						bool b0, b1, b2, b3;
						b0 = !v0->IsUserBit(_featureFlag);
						b1 = !v1->IsUserBit(_featureFlag);
						b2 =	v2->IsUserBit(_featureFlag);
						b3 =	v3->IsUserBit(_featureFlag);
						if( b0 && b1 && b2 && b3)
							vcg::face::FlipEdge< FaceType >(*f, z);

					} // end if (vcg::face::CheckFlipEdge< _Face >(*f, z))
				} // end for( ; e_it!=e_end; e_it++)
			}; //end of FlipEdges
		}; // end of class ExtendedMarchingCubes
		//	/*! @} */
		// end of Doxygen documentation

	} // end of namespace tri
}; // end of namespace vcg

#endif // __VCG_EXTENDED_MARCHING_CUBES
