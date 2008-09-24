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
#ifndef __VCG_MESH_RESAMPLER
#define __VCG_MESH_RESAMPLER

#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/complex/trimesh/update/edges.h>
#include <vcg/complex/trimesh/create/marching_cubes.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/trimesh/closest.h>
#include <vcg/space/box3.h>

namespace vcg {
namespace tri {


/** \addtogroup trimesh */
/*@{*/
/*@{*/
/** Class Resampler.
    This is class reasmpling a mesh using marching cubes methods
		@param OLD_MESH_TYPE (Template Parameter) Specifies the type of mesh to be resampled
		@param NEW_MESH_TYPE (Template Parameter) Specifies the type of output mesh.
 */

template <class OLD_MESH_TYPE,class NEW_MESH_TYPE, class FLT>
	class Resampler : public BasicGrid<FLT>
{
	typedef OLD_MESH_TYPE Old_Mesh;
	typedef NEW_MESH_TYPE New_Mesh;

	//template <class OLD_MESH_TYPE,class NEW_MESH_TYPE>
	class Walker : BasicGrid<float>
	{
	private:
		typedef int VertexIndex;
		typedef  OLD_MESH_TYPE Old_Mesh;
		typedef  NEW_MESH_TYPE New_Mesh;
		typedef typename New_Mesh::CoordType NewCoordType;
		typedef typename New_Mesh::VertexType* VertexPointer;
		typedef typename Old_Mesh::FaceContainer FaceCont;
		typedef typename vcg::GridStaticPtr<typename Old_Mesh::FaceType> GridType;

	protected:

		int SliceSize;
		int	CurrentSlice;
		typedef trimesh::FaceTmark<Old_Mesh> MarkerFace;
		MarkerFace markerFunctor;
		
	
		VertexIndex *_x_cs; // indici dell'intersezioni della superficie lungo gli Xedge della fetta corrente
		VertexIndex	*_y_cs; // indici dell'intersezioni della superficie lungo gli Yedge della fetta corrente
		VertexIndex *_z_cs; // indici dell'intersezioni della superficie lungo gli Zedge della fetta corrente
		VertexIndex *_x_ns; // indici dell'intersezioni della superficie lungo gli Xedge della prossima fetta 
		VertexIndex *_z_ns; // indici dell'intersezioni della superficie lungo gli Zedge della prossima fetta 

		//float *_v_cs;///values of distance fields for each direction in current slice
		//float *_v_ns;///values of distance fields for each direction in next slice

		typedef typename  std::pair<bool,float> field_value;
		field_value* _v_cs;
		field_value* _v_ns;

		New_Mesh	*_newM;
		Old_Mesh	*_oldM;
		GridType _g;
		
	public:
		float max_dim;
		float offset; // an offset value that is always added to the returned value. Useful for extrarcting isosurface  at a different threshold
		/*Walker(Volume_Dataset <short> *Vo,float in,const Box3i &bbox,vcg::Point3i &resolution)
		{*/
		/*	init=in;
			Vol=Vo;*/

	
	
		Walker(const Box3f &_bbox, Point3i _siz )
		{
			this->bbox= _bbox;
			this->siz=_siz;
			ComputeDimAndVoxel();
			
			SliceSize = (this->siz.X()+1)*(this->siz.Z()+1);
			CurrentSlice = 0;
			offset=0;

			_x_cs = new VertexIndex[ SliceSize ];
			_y_cs = new VertexIndex[ SliceSize ];
			_z_cs = new VertexIndex[ SliceSize ];
			_x_ns = new VertexIndex[ SliceSize ];
			_z_ns = new VertexIndex[ SliceSize ];

			_v_cs= new field_value[(this->siz.X()+1)*(this->siz.Z()+1)];
			_v_ns= new field_value[(this->siz.X()+1)*(this->siz.Z()+1)];
			
		};

		~Walker()
		{}

		
		float V(const Point3i &p)
		{
			return (V(p.V(0),p.V(1),p.V(2)));
		}

		float V(int x,int y,int z)
		{
			assert ((y==CurrentSlice)||(y==(CurrentSlice+1)));
			
			//test if it is outside the bb of the mesh
			//vcg::Point3f test=vcg::Point3f((float)x,(float)y,(float)z);
			/*if (!_oldM->bbox.IsIn(test))
				return (1.f);*/
			int index=GetSliceIndex(x,z);
			
			if (y==CurrentSlice)
			{
				return _v_cs[index].second+offset;
			}
			else
			{
				return _v_ns[index].second+offset;
			}
		}
		///return true if the distance form the mesh is less than maxdim and return distance
		bool DistanceFromMesh(int x,int y,int z,Old_Mesh *mesh, float &dist)
		{

			typename Old_Mesh::FaceType *f=NULL;
			const float max_dist = max_dim;
			vcg::Point3f testPt;
			this->IPToP(Point3i(x,y,z),testPt);

			vcg::Point3f closestNorm;
			vcg::Point3f closestPt;
			vcg::Point3f pip;

			// Note that PointDistanceBaseFunctor does not require the edge and plane precomptued.
			// while the PointDistanceFunctor requires them. 
			
			vcg::face::PointDistanceBaseFunctor<typename Old_Mesh::ScalarType > PDistFunct;
			f = _g.GetClosest(PDistFunct,markerFunctor,testPt,max_dist,dist,closestPt);
			
			if (f==NULL) return false;

			InterpolationParameters(*f,closestPt, pip[0], pip[1], pip[2]);
			closestNorm =  (f->V(0)->cN())*pip[0]+ (f->V(1)->cN())*pip[1] + (f->V(2)->cN())*pip[2] ;

			assert(!f->IsD());
				Point3f dir=(testPt-closestPt);
			/*	dist=dir.Norm();*/

				dir.Normalize();
				//direction of normal inside the mesh
				if ((dir*closestNorm)<0)
					dist=-dist;
				//the intersection exist
				return true;
		}

		///compute the values if an entire slice (per y) distances>dig of a cell are signed with double of
		/// the distance of the bb
		void ComputeSliceValues(int slice,field_value *slice_values)
		{
			float dist;
			for (int i=0; i<=this->siz.X(); i++)
			{
				for (int k=0; k<=this->siz.Z(); k++)
					{
						int index=GetSliceIndex(i,k);
						if (DistanceFromMesh(i,slice,k,_oldM,dist))///compute the distance,inside volume of the mesh is negative
						{
							//put computed values in the slice values matrix
							slice_values[index]=field_value(true,dist);
							//end putting values
						}
						else
							slice_values[index]=field_value(false,dist);
					}
			}
		}

		template<class EXTRACTOR_TYPE>
		void ProcessSlice(EXTRACTOR_TYPE &extractor)
		{
			for (int i=0; i<this->siz.X(); i++)
			{
				for (int k=0; k<this->siz.Z(); k++)
				{
						Point3i p1(i,CurrentSlice,k);
						Point3i p2=p1+Point3i(1,1,1);
						extractor.ProcessCell(p1, p2);
				}
			}
		}


		template<class EXTRACTOR_TYPE>
		void BuildMesh(Old_Mesh &old_mesh,New_Mesh &new_mesh,EXTRACTOR_TYPE &extractor,vcg::CallBackPos *cb)
		{
			_newM=&new_mesh;
			_oldM=&old_mesh;
			
			// the following two steps are required to be sure that the point-face distance without precomputed data works well.
			tri::UpdateNormals<Old_Mesh>::PerFaceNormalized(old_mesh);
			tri::UpdateFlags<Old_Mesh>::FaceProjection(old_mesh);
			
			_g.Set(_oldM->face.begin(),_oldM->face.end());
			markerFunctor.SetMesh(&old_mesh);

			_newM->Clear();

			Begin();
			extractor.Initialize();
			
			for (int j=0; j<=this->siz.Y(); j++)
			{
				cb((100*j)/this->siz.Y(),"Marching ");
				ProcessSlice<EXTRACTOR_TYPE>(extractor);//find cells where there is the isosurface and examine it			
				NextSlice();
			}
			extractor.Finalize();
			
			typename New_Mesh::VertexIterator vi;
			for(vi=new_mesh.vert.begin();vi!=new_mesh.vert.end();++vi) 
				if(!(*vi).IsD())
					{
						IPToP((*vi).cP(),(*vi).P());
					}
		}
		
		//return the index of a vertex in slide as it was stored
		int GetSliceIndex(int x,int z)
		{
			VertexIndex index = x+z*(this->siz.X()+1);
			return (index);
		}
	
		//swap slices , the initial value of distance fields ids set as double of bbox of space
		void NextSlice() 
		{
			
			memset(_x_cs, -1, SliceSize*sizeof(VertexIndex));
			memset(_y_cs, -1, SliceSize*sizeof(VertexIndex));
			memset(_z_cs, -1, SliceSize*sizeof(VertexIndex));
			

			std::swap(_x_cs, _x_ns);
			std::swap(_z_cs, _z_ns);

			std::swap(_v_cs, _v_ns);

			CurrentSlice ++;
			
			ComputeSliceValues(CurrentSlice + 1,_v_ns);
		}
		
		//initialize data strucures , the initial value of distance fields ids set as double of bbox of space
		void Begin()
		{

			CurrentSlice = 0;

			memset(_x_cs, -1, SliceSize*sizeof(VertexIndex));
			memset(_y_cs, -1, SliceSize*sizeof(VertexIndex));
			memset(_z_cs, -1, SliceSize*sizeof(VertexIndex));
			memset(_x_ns, -1, SliceSize*sizeof(VertexIndex));
			memset(_z_ns, -1, SliceSize*sizeof(VertexIndex));

			ComputeSliceValues(CurrentSlice,_v_cs);
			ComputeSliceValues(CurrentSlice+1,_v_ns);
		}

		
	
		
		bool Exist(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v)
		{ 
			int i = p1.X();// - _bbox.min.X())/_cell_size.X();
			int z = p1.Z();// - _bbox.min.Z())/_cell_size.Z();
			VertexIndex index = i+z*this->siz.X();

			//VertexIndex index =GetSliceIndex(//
			int v_ind = 0;
			if (p1.X()!=p2.X()) //intersezione della superficie con un Xedge
			{
				if (p1.Y()==CurrentSlice)
				{
					if (_x_cs[index]!=-1)
					{
						v_ind = _x_cs[index];
						v = &_newM->vert[v_ind];
						assert(!v->IsD());
						return true;
					}

				}
				else
				{
					if (_x_ns[index]!=-1)
					{
						v_ind = _x_ns[index];
						v = &_newM->vert[v_ind];
						assert(!v->IsD());
						return true;
					}
				}
				v = NULL;
				return false;
			}
			else if (p1.Y()!=p2.Y()) //intersezione della superficie con un Yedge
			{
				if (_y_cs[index]!=-1)
				{
					v_ind =_y_cs[index];
					v = &_newM->vert[v_ind];
					assert(!v->IsD());
					return true;
				}
				else
				{
					v = NULL;
					return false;
				}

			}
			else if (p1.Z()!=p2.Z())
			//intersezione della superficie con un Zedge
			{
				if (p1.Y()==CurrentSlice)
				{
					if ( _z_cs[index]!=-1)
					{
						v_ind = _z_cs[index];
						v = &_newM->vert[v_ind];
						assert(!v->IsD());
						return true;
					}

				}
				else
				{
					if (_z_ns[index]!=-1)
					{
						v_ind = _z_ns[index];
						v = &_newM->vert[v_ind];
						assert(!v->IsD());
						return true;
					}
				}
				v = NULL;
				return false;
			}
			assert (0);
			return false;
		}

		///interpolate 
		NewCoordType Interpolate(const vcg::Point3i &p1, const vcg::Point3i &p2,int dir)
		{	
			float f1 = (float)V(p1);
			float f2 = (float)V(p2);
			float u = (float) f1/(f1-f2);
			NewCoordType ret=vcg::Point3f((float)p1.V(0),(float)p1.V(1),(float)p1.V(2));
			ret.V(dir) = (float) p1.V(dir)*(1.f-u) + u*(float)p2.V(dir);
			return (ret);
		}

		///if there is a vertex in z axis of a cell return the vertex or create it
		void GetXIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
		{ 
      assert(p1.X()+1 == p2.X());
			assert(p1.Y()   == p2.Y());
			assert(p1.Z()   == p2.Z());
			
			int i = p1.X();// (p1.X() - _bbox.min.X())/_cell_size.X();
			int z = p1.Z();//(p1.Z() - _bbox.min.Z())/_cell_size.Z();
			VertexIndex index = i+z*this->siz.X();
			VertexIndex pos;
			if (p1.Y()==CurrentSlice)
			{
				if ((pos=_x_cs[index])==-1)
				{
					_x_cs[index] = (VertexIndex) _newM->vert.size();
					pos = _x_cs[index];
					Allocator<New_Mesh>::AddVertices( *_newM, 1 );
					v = &_newM->vert[pos];
					v->P()=Interpolate(p1,p2,0);
					return;
				}
			}
			if (p1.Y()==CurrentSlice+1)
			{
				if ((pos=_x_ns[index])==-1)
				{
					_x_ns[index] = (VertexIndex) _newM->vert.size();
					pos = _x_ns[index];
					Allocator<New_Mesh>::AddVertices( *_newM, 1 );
					v = &_newM->vert[pos];
					v->P()=Interpolate(p1,p2,0);
					return;
				}
			}
			v = &_newM->vert[pos];
		}
		
		///if there is a vertex in y axis of a cell return the vertex or create it
		void GetYIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
		{
      assert(p1.X()   == p2.X());
			assert(p1.Y()+1 == p2.Y());
			assert(p1.Z()   == p2.Z());
			
			int i = p1.X(); // (p1.X() - _bbox.min.X())/_cell_size.X();
			int z = p1.Z(); // (p1.Z() - _bbox.min.Z())/_cell_size.Z();
			VertexIndex index = i+z*this->siz.X();
			VertexIndex pos;
			if ((pos=_y_cs[index])==-1)
			{
				_y_cs[index] = (VertexIndex) _newM->vert.size();
				pos = _y_cs[index];
				Allocator<New_Mesh>::AddVertices( *_newM, 1);
				v = &_newM->vert[ pos ];
				v->P()=Interpolate(p1,p2,1);
			}
			v = &_newM->vert[pos];
		}
		
		///if there is a vertex in z axis of a cell return the vertex or create it
		void GetZIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
		{
      assert(p1.X()   == p2.X());
			assert(p1.Y()   == p2.Y());
			assert(p1.Z()+1 == p2.Z());
			
			int i = p1.X(); //(p1.X() - _bbox.min.X())/_cell_size.X();
			int z = p1.Z(); //(p1.Z() - _bbox.min.Z())/_cell_size.Z();
			VertexIndex index = i+z*this->siz.X();

			VertexIndex pos;
			if (p1.Y()==CurrentSlice)
			{
				if ((pos=_z_cs[index])==-1)
				{
					_z_cs[index] = (VertexIndex) _newM->vert.size();
					pos = _z_cs[index];
					Allocator<New_Mesh>::AddVertices( *_newM, 1 );
					v = &_newM->vert[pos];
					v->P()=Interpolate(p1,p2,2);
					return;
				}
			}
			if (p1.Y()==CurrentSlice+1)
			{
				if ((pos=_z_ns[index])==-1)
				{
					_z_ns[index] = (VertexIndex) _newM->vert.size();
					pos = _z_ns[index];
					Allocator<New_Mesh>::AddVertices( *_newM, 1 );
					v = &_newM->vert[pos];
					v->P()=Interpolate(p1,p2,2);
					return;
				}
			}
			v = &_newM->vert[pos];
		}

	};//end class walker

public:

typedef Walker  /*< Old_Mesh,New_Mesh>*/  MyWalker;

typedef vcg::tri::MarchingCubes<New_Mesh, MyWalker> MyMarchingCubes;

///resample the mesh using marching cube algorithm ,the accuracy is the dimension of one cell the parameter
static void Resample(Old_Mesh &old_mesh,New_Mesh &new_mesh,vcg::Point3<int> accuracy,float max_dist, float thr=0, vcg::CallBackPos *cb=0 )
{
	///be sure that the bounding box is updated
	vcg::tri::UpdateBounding<Old_Mesh>::Box(old_mesh);
	
  Box3f volumeBox = old_mesh.bbox;
	volumeBox.Offset(volumeBox.Diag()/10.0f);
	MyWalker	walker(volumeBox,accuracy);

	walker.max_dim=max_dist+fabs(thr);
	walker.offset = - thr; 
	MyMarchingCubes mc(new_mesh, walker);
	walker.BuildMesh(old_mesh,new_mesh,mc,cb);
}
 

};//end class resampler

};//end namespace trimesh
};//end namespace vcg
#endif