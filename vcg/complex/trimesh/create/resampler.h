#ifndef __VCG_MESH_RESAMPLER
#define __VCG_MESH_RESAMPLER

#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/complex/trimesh/update/edges.h>
#include <vcg/complex/trimesh/create/extended_marching_cubes.h>
#include <vcg/complex/trimesh/create/marching_cubes.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/trimesh/closest.h>
#include <vcg/space/box3.h>

namespace vcg {
namespace trimesh {

class RES
{
public:

	enum MarchMode  {MMarchingCubes,MExtendedMarchingCubes} ;
};

/** \addtogroup trimesh */
/*@{*/
/*@{*/
/** Class Resampler.
    This is class reasmpling a mesh using marching cubes methods
		@param OLD_MESH_TYPE (Template Parameter) Specifies the type of mesh to be resampled
		@param NEW_MESH_TYPE (Template Parameter) Specifies the type of output mesh.
		@param MARCHING_ALGORITHM (Template Parameter) Specifies the type of marching cube algorithm (extended or not).
 */

template <class OLD_MESH_TYPE,class NEW_MESH_TYPE>
class Resampler:RES
{
	typedef typename OLD_MESH_TYPE Old_Mesh;
	typedef typename NEW_MESH_TYPE New_Mesh;

	template <class OLD_MESH_TYPE,class NEW_MESH_TYPE>
	class Walker
	{
	private:
		typedef int VertexIndex;
		typedef typename OLD_MESH_TYPE Old_Mesh;
		typedef typename NEW_MESH_TYPE New_Mesh;
		typedef typename New_Mesh::CoordType NewCoordType;
		typedef typename New_Mesh::VertexType* VertexPointer;
		typedef typename Old_Mesh::FaceContainer FaceCont;
		typedef typename GridStaticPtr<FaceCont> Grid;
		typedef typename vcg::Box3<int> BoundingBox;
		typedef vcg::tri::Allocator< New_Mesh > Allocator;

	protected:
		BoundingBox		_bbox;
		vcg::Point3i	_resolution;
		vcg::Point3i	_cell_size;

		float max_dim;

		int _slice_dimension;
		int	_current_slice;
  
	
		VertexIndex *_x_cs; // indici dell'intersezioni della superficie lungo gli Xedge della fetta corrente
		VertexIndex	*_y_cs; // indici dell'intersezioni della superficie lungo gli Yedge della fetta corrente
		VertexIndex *_z_cs; // indici dell'intersezioni della superficie lungo gli Zedge della fetta corrente
		VertexIndex *_x_ns; // indici dell'intersezioni della superficie lungo gli Xedge della prossima fetta 
		VertexIndex *_z_ns; // indici dell'intersezioni della superficie lungo gli Zedge della prossima fetta 

		New_Mesh	*_newM;
		Old_Mesh	*_oldM;
		Grid _g;

	public:

		Walker(const BoundingBox &bbox,vcg::Point3i &resolution)
		{
			assert (resolution.V(0)<=bbox.DimX());
			assert (resolution.V(1)<=bbox.DimY());
			assert (resolution.V(2)<=bbox.DimZ());

			_bbox= bbox;
			_resolution = resolution;
			_cell_size.X() = _bbox.DimX()/_resolution.X();
			_cell_size.Y() = _bbox.DimY()/_resolution.Y();
			_cell_size.Z() = _bbox.DimZ()/_resolution.Z();
		
			///extend bb until the box - resolution and cell matches
			while ((_bbox.DimX()%_cell_size.X())!=0)
					_bbox.max.X()++;

			while ((_bbox.DimY()%_cell_size.Y())!=0)
					_bbox.max.Y()++;

			while ((_bbox.DimZ()%_cell_size.Z())!=0)
					_bbox.max.Z()++;
			
			//exetend bb to 1 cell for each side
			_bbox.max+=_cell_size;
			_bbox.min-=_cell_size;

			///resetting resolution values
			_resolution.X()=_bbox.DimX()/_cell_size.X();
			_resolution.Y()=_bbox.DimY()/_cell_size.Y();
			_resolution.Z()=_bbox.DimZ()/_cell_size.Z();

			///asserting values
			assert(_bbox.DimX()%_cell_size.X()==0);
			assert(_bbox.DimY()%_cell_size.Y()==0);
			assert(_bbox.DimZ()%_cell_size.Z()==0);

			assert(_cell_size.X()*_resolution.X()==_bbox.DimX());
			assert(_cell_size.Y()*_resolution.Y()==_bbox.DimY());
			assert(_cell_size.Z()*_resolution.Z()==_bbox.DimZ());

			_slice_dimension = _resolution.X()*_resolution.Z();
		
			Point3f diag=Point3f((float)_cell_size.V(0),(float)_cell_size.V(1),(float)_cell_size.V(2));
			max_dim=diag.Norm();///diagonal of a cell

			_x_cs = new VertexIndex[ _slice_dimension ];
			_y_cs = new VertexIndex[ _slice_dimension ];
			_z_cs = new VertexIndex[ _slice_dimension ];
			_x_ns = new VertexIndex[ _slice_dimension ];
			_z_ns = new VertexIndex[ _slice_dimension ];
			
		};

		~Walker()
		{}

		template<class EXTRACTOR_TYPE>
		void BuildMesh(Old_Mesh &old_mesh,New_Mesh &new_mesh,EXTRACTOR_TYPE &extractor)
		{
			_newM=&new_mesh;
			_oldM=&old_mesh;

			Point3f min=Point3f((float)_bbox.min.V(0),(float)_bbox.min.V(1),(float)_bbox.min.V(2));
			Point3f max=Point3f((float)_bbox.max.V(0),(float)_bbox.max.V(1),(float)_bbox.max.V(2));

			vcg::Box3<float> BBf=vcg::Box3<float>(min,max);

			_g.SetBBox(BBf);

			_g.Set(_oldM->face);

			_newM->Clear();
			vcg::Point3i p1, p2;

			Begin();
			extractor.Initialize();
			for (int j=_bbox.min.Y(); j<_bbox.max.Y()-_cell_size.Y(); j+=_cell_size.Y())
			{
				for (int i=_bbox.min.X(); i<_bbox.max.X()-_cell_size.X(); i+=_cell_size.X())
				{
					for (int k=_bbox.min.Z(); k<_bbox.max.Z()-_cell_size.Z(); k+=_cell_size.Z())
					{
						p1.X()=i;
						p1.Y()=j;
						p1.Z()=k;
						p2.X()=i+_cell_size.X();
						p2.Y()=j+_cell_size.Y();
						p2.Z()=k+_cell_size.Z();	
						if (ExistIntersection(p1,p2))
							extractor.ProcessCell(p1, p2);
					}
				}
				NextSlice();
			}
			extractor.Finalize();

			/*_newM= NULL;*/
		};
	
		float V(Point3i p)
		{
			return (V(p.V(0),p.V(1),p.V(2)));
		}
	
		///control if exist any intersection with the mesh using all points of the cell
		bool ExistIntersection(Point3i min,Point3i max)
		{
			vcg::Point3i _corners[8];

			_corners[0].X()=min.X();		_corners[0].Y()=min.Y();		_corners[0].Z()=min.Z();
			_corners[1].X()=max.X();		_corners[1].Y()=min.Y();		_corners[1].Z()=min.Z();
			_corners[2].X()=max.X();		_corners[2].Y()=max.Y();		_corners[2].Z()=min.Z();
			_corners[3].X()=min.X();		_corners[3].Y()=max.Y();		_corners[3].Z()=min.Z();
			_corners[4].X()=min.X();		_corners[4].Y()=min.Y();		_corners[4].Z()=max.Z();
			_corners[5].X()=max.X();		_corners[5].Y()=min.Y();		_corners[5].Z()=max.Z();
			_corners[6].X()=max.X();		_corners[6].Y()=max.Y();		_corners[6].Z()=max.Z();
			_corners[7].X()=min.X();		_corners[7].Y()=max.Y();		_corners[7].Z()=max.Z();

			////control if there is a face
			vcg::Point3f Norm;
			Old_Mesh::FaceType *f=NULL;
			float distm;
			Point3f pip;
			Point3f Target;
			vcg::Point3f test;
			for (int i=0;i<8;i++)
			{
				f=NULL;
				distm=max_dim;
				test=vcg::Point3f((float)_corners[i].X(),(float)_corners[i].Y(),(float)_corners[i].Z());
				vcg::trimesh::Closest<Old_Mesh,Grid,float>((*_oldM),test,_g,distm,Norm,Target,f,pip);
				if (f==NULL)
					return false;
			}
			return true;
		}
	
		///si potrebbe ottimizzare sfruttando valori che in realta' ho gia' calcolato
		float V(int pi, int pj, int pk)
		{
			vcg::Point3f Norm;
			Old_Mesh::FaceType *f=NULL;
			float dist=max_dim;;
			vcg::Point3f test=vcg::Point3f((float)pi,(float)pj,(float)pk);
	
			Point3f pip;
			Point3f Target;
			vcg::trimesh::Closest<Old_Mesh,Grid,float>((*_oldM),test,_g,dist,Norm,Target,f,pip);

			assert(f!=NULL);

			Point3f dir=(test-Target);

			//dist=dir.Norm();
			dir=dir.Normalize();

			//direction of normal inside or outside the mesh
 			if ((f->N()*dir)>0)
				return (dist);
			else
				return (-dist);
		}

		bool Exist(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v)
		{ 
			return false;
		//int i_idx = p1.X()-_bbox.min.X();
		//int k_idx = p2.Z()-_bbox.min.Z();
		//int index = i_idx+k_idx*_resolution.X();
		//if (p1.X()!=p2.X()) //intersezione della superficie con un Xedge
		//	return (p1.Y()==_current_slice)? _x_cs[index]!=-1 : _x_ns[index]!=-1;
		//else if (p1.Y()!=p2.Y()) //intersezione della superficie con un Yedge
		//	return _y_cs[index]!=-1;
		//else if (p1.Z()!=p2.Z()) //intersezione della superficie con un Zedge
		//	return (p1.Y()==_current_slice)? _z_cs[index]!=-1 : _z_ns[index]!=-1;
		}


		NewCoordType Interpolate(const vcg::Point3i &p1, const vcg::Point3i &p2,int dir)
		{	
			float f1 = V(p1);
			float f2 = V(p2);
			float u = (float) f1/(f1-f2);
			NewCoordType ret=vcg::Point3f((float)p1.V(0),(float)p1.V(1),(float)p1.V(2));
			ret.V(dir) = (float) p1.V(dir)*(1-u) + u*p2.V(dir);
			return (ret);
		}

		void GetXIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
		{ 
			int i = (p1.X() - _bbox.min.X())/_cell_size.X();
			int z = (p1.Z() - _bbox.min.Z())/_cell_size.Z();
			VertexIndex index = i+z*_resolution.X();
			VertexIndex pos;
			if (p1.Y()==_current_slice)
			{
				if ((pos=_x_cs[index])==-1)
				{
					_x_cs[index] = (VertexIndex) _newM->vert.size();
					pos = _x_cs[index];
					Allocator::AddVertices( *_newM, 1 );
					v = &_newM->vert[pos];
					v->P()=Interpolate(p1,p2,0);
					return;
				}
			}
			if (p1.Y()==_current_slice+_cell_size.Y())
			{
				if ((pos=_x_ns[index])==-1)
				{
					_x_ns[index] = (VertexIndex) _newM->vert.size();
					pos = _x_ns[index];
					Allocator::AddVertices( *_newM, 1 );
					v = &_newM->vert[pos];
					v->P()=Interpolate(p1,p2,0);
					return;
				}
			}
			v = &_newM->vert[pos];
		}

		void GetYIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
		{
			int i = (p1.X() - _bbox.min.X())/_cell_size.X();
			int z = (p1.Z() - _bbox.min.Z())/_cell_size.Z();
			VertexIndex index = i+z*_resolution.X();
			VertexIndex pos;
			if ((pos=_y_cs[index])==-1)
			{
				_y_cs[index] = (VertexIndex) _newM->vert.size();
				pos = _y_cs[index];
				Allocator::AddVertices( *_newM, 1);
				v = &_newM->vert[ pos ];
				v->P()=Interpolate(p1,p2,1);
			}
			v = &_newM->vert[pos];
		}

		void GetZIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
		{
			int i = (p1.X() - _bbox.min.X())/_cell_size.X();
			int z = (p1.Z() - _bbox.min.Z())/_cell_size.Z();
			VertexIndex index = i+z*_resolution.X();

			VertexIndex pos;
			if (p1.Y()==_current_slice)
			{
				if ((pos=_z_cs[index])==-1)
				{
					_z_cs[index] = (VertexIndex) _newM->vert.size();
					pos = _z_cs[index];
					Allocator::AddVertices( *_newM, 1 );
					v = &_newM->vert[pos];
					v->P()=Interpolate(p1,p2,2);
					return;
				}
			}
			if (p1.Y()==_current_slice+_cell_size.Y())
			{
				if ((pos=_z_ns[index])==-1)
				{
					_z_ns[index] = (VertexIndex) _newM->vert.size();
					pos = _z_ns[index];
					Allocator::AddVertices( *_newM, 1 );
					v = &_newM->vert[pos];
					v->P()=Interpolate(p1,p2,2);
					return;
				}
			}
			v = &_newM->vert[pos];
		}



		void NextSlice() 
		{
			memset(_x_cs, -1, _slice_dimension*sizeof(VertexIndex));
			memset(_y_cs,	-1, _slice_dimension*sizeof(VertexIndex));
			memset(_z_cs, -1, _slice_dimension*sizeof(VertexIndex));

			std::swap(_x_cs, _x_ns);
			std::swap(_z_cs, _z_ns);		
			
			_current_slice += _cell_size.Y();
		}

		void Begin()
		{
			_current_slice = _bbox.min.Y();

			memset(_x_cs, -1, _slice_dimension*sizeof(VertexIndex));
			memset(_y_cs, -1, _slice_dimension*sizeof(VertexIndex));
			memset(_z_cs, -1, _slice_dimension*sizeof(VertexIndex));
			memset(_x_ns, -1, _slice_dimension*sizeof(VertexIndex));
			memset(_z_ns, -1, _slice_dimension*sizeof(VertexIndex));
		
		}
	};//end class walker

public:

typedef typename  Walker< Old_Mesh,New_Mesh> MyWalker;

typedef typename vcg::tri::MarchingCubes<New_Mesh, MyWalker> MarchingCubes;
typedef typename vcg::tri::ExtendedMarchingCubes<New_Mesh, MyWalker> ExtendedMarchingCubes;

///resample the mesh using marching cube algorithm ,the accuracy is the dimension of one cell the parameter
template <RES::MarchMode mm>
static void Resample(Old_Mesh &old_mesh,New_Mesh &new_mesh,vcg::Point3<int> accuracy)
{
	new_mesh.Clear();
	if (Old_Mesh::HasPerFaceNormal())
		vcg::tri::UpdateNormals<Old_Mesh>::PerFaceNormalized(old_mesh);
	if (Old_Mesh::HasPerVertexNormal())
		vcg::tri::UpdateNormals<Old_Mesh>::PerVertexNormalized(old_mesh);
		///the mesh must have plane for ugrid
	if (!Old_Mesh::FaceType::HasEdgePlane())
		assert(0);
	else
		vcg::tri::UpdateEdges<Old_Mesh>::Set(old_mesh);

	///be sure that the bounding box is updated
	vcg::tri::UpdateBounding<Old_Mesh>::Box(old_mesh);
	

	// MARCHING CUBES CALLS
	Point3i min=Point3i((int)ceil(old_mesh.bbox.min.V(0)),(int)ceil(old_mesh.bbox.min.V(1)),(int)ceil(old_mesh.bbox.min.V(2)));
	Point3i max=Point3i((int)ceil(old_mesh.bbox.max.V(0)),(int)ceil(old_mesh.bbox.max.V(1)),(int)ceil(old_mesh.bbox.max.V(2)));

	///exetend out to BB for resample the limits of the mesh
	/*min-=Point3i(1,1,1);
	max+=Point3i(1,1,1);*/
	vcg::Box3<int> boxInt=Box3<int>(min,max);
	
	float rx=((float)boxInt.DimX())/(float)accuracy.X();
	float ry=((float)boxInt.DimY())/(float)accuracy.Y();
	float rz=((float)boxInt.DimZ())/(float)accuracy.Z();

	int rxi=(int)ceil(rx);
	int ryi=(int)ceil(ry);
	int rzi=(int)ceil(rz);
	
	Point3i res=Point3i(rxi,ryi,rzi);

	MyWalker	walker(boxInt,res);
	if (mm==MMarchingCubes)
	{
		MarchingCubes mc(new_mesh, walker);
		walker.BuildMesh<MarchingCubes>(old_mesh,new_mesh,mc);
	}
	else if (mm==MExtendedMarchingCubes)
	{
		ExtendedMarchingCubes mc(new_mesh, walker,30);
		walker.BuildMesh<ExtendedMarchingCubes>(old_mesh,new_mesh,mc);
	}
	
	
	
	if (New_Mesh::HasFFTopology())
		vcg::tri::UpdateTopology<New_Mesh>::FaceFace(new_mesh);
	if (New_Mesh::HasVFTopology())
		vcg::tri::UpdateTopology<New_Mesh>::VertexFace(new_mesh);	
	if (New_Mesh::HasPerFaceNormal())
		vcg::tri::UpdateNormals<New_Mesh>::PerFaceNormalized(new_mesh);
	if (New_Mesh::HasPerVertexNormal())
		vcg::tri::UpdateNormals<New_Mesh>::PerVertexNormalized(new_mesh);

}
 

};//end class resampler

};//end namespace trimesh
};//end namespace vcg
#endif