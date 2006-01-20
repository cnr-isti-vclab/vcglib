#ifndef __VCG_MESH_RESAMPLER
#define __VCG_MESH_RESAMPLER

#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/complex/trimesh/update/edges.h>
//#include <vcg/complex/trimesh/create/extended_marching_cubes.h>
#include <vcg/complex/trimesh/create/marching_cubes.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/trimesh/closest.h>
#include <vcg/space/box3.h>

//#include <volume_dataset.h>//debugghe

namespace vcg {
namespace trimesh {


/** \addtogroup trimesh */
/*@{*/
/*@{*/
/** Class Resampler.
    This is class reasmpling a mesh using marching cubes methods
		@param OLD_MESH_TYPE (Template Parameter) Specifies the type of mesh to be resampled
		@param NEW_MESH_TYPE (Template Parameter) Specifies the type of output mesh.
 */

template <class OLD_MESH_TYPE,class NEW_MESH_TYPE>
class Resampler
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
		typedef typename vcg::GridStaticPtr<typename Old_Mesh::FaceType> GridType;
		typedef typename vcg::Box3<int> BoundingBox;
		//typedef typename std::pair<vcg::Point3i,vcg::Point3i> PointPair;
		typedef vcg::tri::Allocator< New_Mesh > Allocator;

	protected:
		BoundingBox		_bbox;
		vcg::Point3i	_resolution;
		vcg::Point3i	_cell_size;

		
		float dim_diag;

		int _slice_dimension;
		int	_current_slice;
  
	
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
		/*Walker(Volume_Dataset <short> *Vo,float in,const BoundingBox &bbox,vcg::Point3i &resolution)
		{*/
		/*	init=in;
			Vol=Vo;*/

		void SetBBParameters()
		{
			_cell_size.X() =_bbox.DimX()/_resolution.X();
			_cell_size.Y() =_bbox.DimY()/_resolution.Y();
			_cell_size.Z() =_bbox.DimZ()/_resolution.Z();
		
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

			_slice_dimension = (_resolution.X()+1)*(_resolution.Z()+1);
		
			//Point3f diag=Point3f((float)_cell_size.V(0),(float)_cell_size.V(1),(float)_cell_size.V(2));
			//max_dim=diag.Norm();///diagonal of a cell
			//
			_current_slice = _bbox.min.Y();

			Point3f minD=Point3f((float)_bbox.min.V(0),(float)_bbox.min.V(1),(float)_bbox.min.V(2));
			Point3f maxD=Point3f((float)_bbox.max.V(0),(float)_bbox.max.V(1),(float)_bbox.max.V(2));
			/*Point3f d=(maxD-minD);
			dim_diag=d.Norm();*/
		}

		Walker(const BoundingBox &bbox,vcg::Point3i &resolution)
		{
			assert (resolution.V(0)<=bbox.DimX());
			assert (resolution.V(1)<=bbox.DimY());
			assert (resolution.V(2)<=bbox.DimZ());

			_bbox= bbox;

			_resolution = resolution;

			SetBBParameters();

			_x_cs = new VertexIndex[ _slice_dimension ];
			_y_cs = new VertexIndex[ _slice_dimension ];
			_z_cs = new VertexIndex[ _slice_dimension ];
			_x_ns = new VertexIndex[ _slice_dimension ];
			_z_ns = new VertexIndex[ _slice_dimension ];

			_v_cs= new field_value[(_resolution.X()+1)*(_resolution.Z()+1)];
			_v_ns= new field_value[(_resolution.X()+1)*(_resolution.Z()+1)];
			
		};

		~Walker()
		{}

		
		float V(Point3i p)
		{
			return (V(p.V(0),p.V(1),p.V(2)));
		}

		float V(int x,int y,int z)
		{
			assert ((y==_current_slice)||(y==(_current_slice+_cell_size.Y())));
			
			//test if it is outside the bb of the mesh
			//vcg::Point3f test=vcg::Point3f((float)x,(float)y,(float)z);
			/*if (!_oldM->bbox.IsIn(test))
				return (1.f);*/
			int index=GetSliceIndex(x,z);
			
			if (y==_current_slice)
			{
				//assert(_v_cs[index]<dim_diag);
				assert(_v_cs[index].first);
				return _v_cs[index].second;
			}
			else
			{
				//assert(_v_ns[index]<dim_diag);
				assert(_v_ns[index].first);
				return _v_ns[index].second;
			}
		}
		///return true if the distance form the mesh is less than maxdim and return distance
		bool DistanceFromMesh(int x,int y,int z,Old_Mesh *mesh,float &dist)
		{

			Old_Mesh::FaceType *f=NULL;
			//float distm=max_dim;
			dist=max_dim;
			vcg::Point3f test=vcg::Point3f((float)x,(float)y,(float)z);

			////test if it is outside the bb of the mesh
			/*if (!_oldM->bbox.IsIn(test))
			{
				dist=1.f;
				return true;
			}*/

			vcg::Point3f Norm;
			vcg::Point3f Target;
			vcg::Point3f pip;

			//vcg::tri::get<Old_Mesh,GridType,float>((*mesh),test,_g,dist,Norm,Target,f,pip);
			
			f= vcg::trimesh::GetClosestFace<Old_Mesh,GridType>( *mesh,_g,test,max_dim,dist,Target,Norm,pip);
		
			if (f==NULL)
					return false;
			else
			{
				assert(!f->IsD());
				Point3f dir=(test-Target);
			/*	dist=dir.Norm();*/

				dir.Normalize();
				//direction of normal inside the mesh
				if ((dir*Norm)<0)
					dist=-dist;
				//the intersection exist
				return true;
			}
		}

		///compute the values if an entire slice (per y) distances>dig of a cell are signed with double of
		/// the distance of the bb
		void CumputeSliceValues(int slice,field_value *slice_values)
		{
			float dist;
			for (int i=_bbox.min.X(); i<=_bbox.max.X(); i+=_cell_size.X())
			{
				for (int k=_bbox.min.Z(); k<=_bbox.max.Z(); k+=_cell_size.Z())
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
		void ProcessSlice(std::vector<vcg::Point3i> cells,EXTRACTOR_TYPE &extractor)
		{
			std::vector<vcg::Point3i>::iterator it;

			for (it=cells.begin();it<cells.end();it++)
			{
						assert((*it).Y()==_current_slice);
						assert(V(*it)<=max_dim);
						assert(_bbox.IsIn(*it));
						vcg::Point3i p1=(*it)+_cell_size;
						assert((*it)<_bbox.max);
						assert(p1<=_bbox.max);
						extractor.ProcessCell((*it), p1);
			}

		}

		void SetGrid()
		{	
			_g.Set(_oldM->face.begin(),_oldM->face.end());
		}

		template<class EXTRACTOR_TYPE>
		void BuildMesh(Old_Mesh &old_mesh,New_Mesh &new_mesh,EXTRACTOR_TYPE &extractor)
		{
			_newM=&new_mesh;
			_oldM=&old_mesh;
			
			SetGrid();

			_newM->Clear();

			vcg::Point3i p1, p2;
			Begin();
			extractor.Initialize();
			for (int j=_bbox.min.Y(); j<=_bbox.max.Y()-_cell_size.Y(); j+=_cell_size.Y())
			{
				ProcessSlice<EXTRACTOR_TYPE>(FindCells(),extractor);//find cells where there is the isosurface and examine it			
				NextSlice();
			}
			extractor.Finalize();
			/*_newM= NULL;*/
		}
		
		//return the index of a vertex in slide as it was stored
		int GetSliceIndex(int x,int z)
		{
			int ii = (x - _bbox.min.X())/_cell_size.X();
			int zz = (z - _bbox.min.Z())/_cell_size.Z();
			VertexIndex index = ii+zz*(_resolution.X()+1);
			return (index);
		}

		///return true if exist in the cell one value <0 and another one >0
		bool FindMinMax(vcg::Point3i min,vcg::Point3i max)
		{
			assert((min.X()<max.X())&&(min.Y()<max.Y())&&(min.Z()<max.Z()));

			vcg::Point3i _corners[8];
			///control for each corner of the 
			_corners[0].X()=min.X();		_corners[0].Y()=min.Y();		_corners[0].Z()=min.Z();
			_corners[1].X()=max.X();		_corners[1].Y()=min.Y();		_corners[1].Z()=min.Z();
			_corners[2].X()=max.X();		_corners[2].Y()=max.Y();		_corners[2].Z()=min.Z();
			_corners[3].X()=min.X();		_corners[3].Y()=max.Y();		_corners[3].Z()=min.Z();
			_corners[4].X()=min.X();		_corners[4].Y()=min.Y();		_corners[4].Z()=max.Z();
			_corners[5].X()=max.X();		_corners[5].Y()=min.Y();		_corners[5].Z()=max.Z();
			_corners[6].X()=max.X();		_corners[6].Y()=max.Y();		_corners[6].Z()=max.Z();
			_corners[7].X()=min.X();		_corners[7].Y()=max.Y();		_corners[7].Z()=max.Z();

			float min_value=max_dim;
			float max_value=-max_dim;
			field_value value;
			for (int i=0;i<8;i++)
			{
				//if one value is > that bbox.diag this value is not valid
				//that is the mark

				if (_corners[i].Y()==_current_slice)
					value=_v_cs[GetSliceIndex(_corners[i].X(),_corners[i].Z())];	
				else
					value=_v_ns[GetSliceIndex(_corners[i].X(),_corners[i].Z())];

				if (value.first==false)
					return false;

				//assign new values of min and max
				if (value.second<min_value)
					min_value=value.second;
				if (value.second>max_value)
					max_value=value.second;
			}	

			/////do not test with zero..
			if ((min_value<=0.f)&&(max_value>=0.f))
				return true;

			return false;
		}

		///filter the cells from to_hexamine vector to the ones that 
		/// min and max of the cell are <0 and >0
		std::vector<vcg::Point3i> FindCells()
		{
			std::vector<vcg::Point3i> res;
			for (int i=_bbox.min.X(); i<=_bbox.max.X()-_cell_size.X(); i+=_cell_size.X())
			{
				for (int k=_bbox.min.Z(); k<=_bbox.max.Z()-_cell_size.Z(); k+=_cell_size.Z())
				{
					int x0=i;
					int y0=_current_slice;
					int z0=k;
					int x1=x0+_cell_size.X();
					int y1=y0+_cell_size.Y();
					int z1=z0+_cell_size.Z();
					vcg::Point3i p0=Point3i(x0,y0,z0);
					vcg::Point3i p1=Point3i(x1,y1,z1);
					assert(p0<_bbox.max);
					if (FindMinMax(p0,p1))
						res.push_back(p0);
				}
			}
			return res;
		}

		//swap slices , the initial value of distance fields ids set as double of bbox of space
		void NextSlice() 
		{
			
			memset(_x_cs, -1, _slice_dimension*sizeof(VertexIndex));
			memset(_y_cs, -1, _slice_dimension*sizeof(VertexIndex));
			memset(_z_cs, -1, _slice_dimension*sizeof(VertexIndex));
			

			std::swap(_x_cs, _x_ns);
			std::swap(_z_cs, _z_ns);

			std::swap(_v_cs, _v_ns);

			_current_slice += _cell_size.Y();
			
			//memset(_v_ns, dim_diag*2.f, _slice_dimension*sizeof(float));
			//memset(_v_ns, field_value(false,0.f), _slice_dimension*sizeof(field_value));

			CumputeSliceValues(_current_slice+ _cell_size.Y(),_v_ns);
		}
		
		//initialize data strucures , the initial value of distance fields ids set as double of bbox of space
		void Begin()
		{

			_current_slice = _bbox.min.Y();

			memset(_x_cs, -1, _slice_dimension*sizeof(VertexIndex));
			memset(_y_cs, -1, _slice_dimension*sizeof(VertexIndex));
			memset(_z_cs, -1, _slice_dimension*sizeof(VertexIndex));
			memset(_x_ns, -1, _slice_dimension*sizeof(VertexIndex));
			memset(_z_ns, -1, _slice_dimension*sizeof(VertexIndex));

			/*memset(_v_cs, dim_diag*2.f, _slice_dimension*sizeof(float));
			memset(_v_ns, dim_diag*2.f, _slice_dimension*sizeof(float));*/
			
			/*memset(_v_cs, field_value(false,0.f), _slice_dimension*sizeof(field_value));
			memset(_v_ns, field_value(false,0.f), _slice_dimension*sizeof(field_value));*/

			CumputeSliceValues(_current_slice,_v_cs);
			CumputeSliceValues(_current_slice+_cell_size.Y(),_v_ns);
		
		}

		
	
		
		bool Exist(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v)
		{ 
			int i = (p1.X() - _bbox.min.X())/_cell_size.X();
			int z = (p1.Z() - _bbox.min.Z())/_cell_size.Z();
			VertexIndex index = i+z*_resolution.X();

			//VertexIndex index =GetSliceIndex(//
			int v_ind = 0;
			if (p1.X()!=p2.X()) //intersezione della superficie con un Xedge
			{
				if (p1.Y()==_current_slice)
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
				if (p1.Y()==_current_slice)
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
			assert ((p1.Y()==_current_slice)||(p1.Y()==(_current_slice+_cell_size.Y())));

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
		
		///if there is a vertex in y axis of a cell return the vertex or create it
		void GetYIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
		{
			assert ((p1.Y()==_current_slice)||(p1.Y()==(_current_slice+_cell_size.Y())));

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
		
		///if there is a vertex in z axis of a cell return the vertex or create it
		void GetZIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
		{
			assert ((p1.Y()==_current_slice)||(p1.Y()==(_current_slice+_cell_size.Y())));

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

	};//end class walker

public:

typedef typename  Walker< Old_Mesh,New_Mesh> MyWalker;

typedef typename vcg::tri::MarchingCubes<New_Mesh, MyWalker> MarchingCubes;

///resample the mesh using marching cube algorithm ,the accuracy is the dimension of one cell the parameter
static void Resample(Old_Mesh &old_mesh,New_Mesh &new_mesh,vcg::Point3<int> accuracy,float max_dist)
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

	
	vcg::Box3<int> boxInt=Box3<int>(min,max);


	float rx=((float)boxInt.DimX())/(float)accuracy.X();
	float ry=((float)boxInt.DimY())/(float)accuracy.Y();
	float rz=((float)boxInt.DimZ())/(float)accuracy.Z();

	int rxi=(int)ceil(rx);
	int ryi=(int)ceil(ry);
	int rzi=(int)ceil(rz);
	
	Point3i res=Point3i(rxi,ryi,rzi);

	MyWalker	walker(boxInt,res);

	walker.max_dim=max_dist;

	/*new_mesh.vert.reserve(old_mesh.vn*2);
	new_mesh.face.reserve(old_mesh.fn*2);*/

	/*if (mm==MMarchingCubes)
	{*/
		MarchingCubes mc(new_mesh, walker);
		walker.BuildMesh<MarchingCubes>(old_mesh,new_mesh,mc);
	/*}*/
	/*else if (mm==MExtendedMarchingCubes)
	{
		ExtendedMarchingCubes mc(new_mesh, walker,30);
		walker.BuildMesh<ExtendedMarchingCubes>(old_mesh,new_mesh,mc);
	}*/
	
	
	
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