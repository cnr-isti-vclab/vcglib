#include<vector>
#include<vcg/space/point3.h>
#include<vcg/space/plane3.h>
#include<vcg/space/segment3.h>
#include<vcg/space/intersection3.h>
#include<vcg/space/index/grid_static_ptr.h>



namespace vcg{

/** \addtogroup complex */
/*@{*/
/** 
    Function computing the intersection between  a grid and a plane. It returns all the cells intersected
*/
template < typename  GridType,typename ScalarType>
bool Intersect(   GridType & grid,Plane3<ScalarType> plane, vector<typename GridType::Cell *> &cells){					
	Point3d p,_d;
	Plane3d pl;
	_d.Import(plane.Direction());
	pl.SetDirection(_d);
	pl.SetOffset(plane.Offset());
	for( int ax = 0; ax <3; ++ax)
			{ int axis = ax;
				int axis0 = (axis+1)%3;
				int axis1 = (axis+2)%3;
				int i,j;
				Point3i pi;

				Segment3<double> seg;
				seg.P0().Import(grid.bbox.min);
				seg.P1().Import(grid.bbox.min);
				seg.P1()[axis] = grid.bbox.max[axis];

				for(i = 0 ; i <= grid.siz[axis0]; ++i){
					for(j = 0 ; j <= grid.siz[axis1]; ++j)
						{
							seg.P0()[axis0] = grid.bbox.min[axis0]+ (i+0.1) * grid.voxel[axis0] ;
							seg.P1()[axis0] = grid.bbox.min[axis0]+ (i+0.1) * grid.voxel[axis0];
							seg.P0()[axis1] = grid.bbox.min[axis1]+ (j+0.1) * grid.voxel[axis1];
							seg.P1()[axis1] = grid.bbox.min[axis1]+ (j+0.1) * grid.voxel[axis1];
							if ( Intersection(pl,seg,p))
								{
									pi[axis] =	min(max(0,floor((p[axis ]-grid.bbox.min[axis])/grid.voxel[axis])),grid.siz[axis]);
									pi[axis0] = i;
									pi[axis1] = j;
									grid.Grid(pi,axis,cells);
								}
						}
					}
			}
		sort(cells.begin(),cells.end());
		cells.erase(unique(cells.begin(),cells.end()),cells.end());
		
		return false;
	}

/*@}*/

/** \addtogroup complex */
/*@{*/
/** 
    Function computing the intersection between  a trimesh and a plane. It returns an EdgeMesh.
		Note: This version always returns a segment for each triangle of the mesh which intersects with the plane. In other
		words there are 2*n vertices where n is the number of segments fo the mesh. You can run vcg::edge::Unify to unify
		the vertices closer that a given value epsilon. Note that, due to subtraction error during triangle plane intersection,
		it is not safe to put epsilon to 0. 
*/
// TODO si dovrebbe considerare la topologia face-face della trimesh per derivare quella della edge mesh..
template < typename  TriMeshType, typename EdgeMeshType, class ScalarType>
bool Intersection( TriMeshType & m, Plane3<ScalarType>  pl,EdgeMeshType & em,double& ave_length,
								 typename GridStaticPtr<typename TriMeshType::FaceContainer> *grid,
									typename vector< typename GridStaticPtr<typename TriMeshType::FaceContainer>::Cell* >& cells){
		typedef typename TriMeshType::FaceContainer FaceContainer;
		typedef GridStaticPtr<FaceContainer> GridType;
		EdgeMeshType::VertexIterator vi;
		TriMeshType::FaceIterator fi;

		Intersect(*grid,pl,cells);
		Segment3<ScalarType> seg;
		ave_length = 0.0;
		vector<GridType::Cell*>::iterator ic;
		GridType::Cell fs,ls;
		for(ic = cells.begin(); ic != cells.end();++ic)
			{
				grid->Grid(*ic,fs,ls);
				GridType::Link * lk = fs;
				while(lk != ls){
					TriMeshType::FaceType & face = *(lk->Elem());
						if(!face.IsS())
						{
							face.SetS();
							if(vcg::Intersection(pl,face,seg))// intersezione piano triangolo
								{
									// add to em
									ave_length+=seg.Length();
									vcg::edge::Allocator<EdgeMeshType>::AddEdges(em,1);
									vi = vcg::edge::Allocator<EdgeMeshType>::AddVertices(em,2);
									(*vi).P() = seg.P0();
									em.edges.back().V(0) = &(*vi); 
									vi++;
									(*vi).P() = seg.P1();
									em.edges.back().V(1) = &(*vi); 
								}
					 	}//endif 
						lk++;
					}//end while
			}
		ave_length/=em.en;

		return true;
	}
/*@}*/
} // end namespace vcg