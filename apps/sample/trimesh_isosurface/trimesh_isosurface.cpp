#include <stdio.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/space/point3.h>
#include <vcg/space/box3.h>
#include <vcg/math/perlin_noise.h>
#include "trivial_walker.h"
#include <vcg/complex/trimesh/create/marching_cubes.h>
#include <vcg/complex/trimesh/create/extended_marching_cubes.h>

using namespace std;
using namespace vcg;

#include <vcg/simplex/vertex/vertex.h>
#include <vcg/simplex/face/face.h>
#include <vcg/complex/trimesh/base.h>
#include <vcg/complex/trimesh/allocate.h>

typedef float ScalarType;

class MyEdge;
class MyFace;
class MyVertex : public vcg::Vertex< ScalarType, MyEdge, MyFace > {};
class MyFace		: public vcg::Face< MyVertex, MyEdge, MyFace> {};
class MyMesh		: public vcg::tri::TriMesh< std::vector< MyVertex>, std::vector< MyFace > > {};

template <class VOX_TYPE>
class Volume
{
public:
  typedef VOX_TYPE VoxelType;
  vector<VoxelType> Vol;
  
  Point3i sz;   /// Dimensioni griglia come numero di celle per lato
  
  const Point3i &ISize() {return sz;};   /// Dimensioni griglia come numero di celle per lato
 	
  void Init(Point3i _sz)
  {
    sz=_sz;
    Vol.resize(sz[0]*sz[1]*sz[2]);
  }

  float Val(const int &x,const int &y,const int &z) const {
      return cV(x,y,z).V(); 
    //else return numeric_limits<float>::quiet_NaN( ); 
  }

  float &Val(const int &x,const int &y,const int &z) {
      return V(x,y,z).V(); 
    //else return numeric_limits<float>::quiet_NaN( ); 
  }

	VOX_TYPE &V(const int &x,const int &y,const int &z) {
		return Vol[x+y*sz[0]+z*sz[0]*sz[1]]; 
	}

	const VOX_TYPE &cV(const int &x,const int &y,const int &z) const {
		return Vol[x+y*sz[0]+z*sz[0]*sz[1]]; 
	}


enum { XAxis=0,YAxis=1,ZAxis=2} VolumeAxis;

template < class VertexPointerType, enum VolumeAxis AxisVal >
  void GetIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
{
			float f1 = Val(p1.X(), p1.Y(), p1.Z())-thr;
			float f2 = Val(p2.X(), p2.Y(), p2.Z())-thr;
			float u = (float) f1/(f1-f2);
			if(AxisVal==XAxis) v->P().X() = (float) p1.X()*(1-u) + u*p2.X();
                    else v->P().X() = (float) p1.X();
			if(AxisVal==YAxis) v->P().Y() = (float) p1.Y()*(1-u) + u*p2.Y();
			              else v->P().Y() = (float) p1.Y();
			if(AxisVal==ZAxis) v->P().Z() = (float) p1.Z()*(1-u) + u*p2.Z();
			              else v->P().Z() = (float) p1.Z();
}

template < class VertexPointerType >
  void GetXIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
{ GetIntercept<VertexPointerType,XAxis>(p1,p2,v,thr); }

template < class VertexPointerType >
  void GetYIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
{ GetIntercept<VertexPointerType,YAxis>(p1,p2,v,thr); }

template < class VertexPointerType >
  void GetZIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
{ GetIntercept<VertexPointerType,ZAxis>(p1,p2,v,thr); }
};
template <class VolumeType>
class RawVolumeImporter
{
public:
  enum DataType
{
		// Funzioni superiori
  UNDEF=0,
  BYTE=1;
  SHORT=2;
  FLOAT=3;
};

static bool Open(const char *filename, VolumeType &V, Point3i sz, DataType d)
{
return true;
}
};

class SimpleVoxel
{
private:
  float _v;
public:
  float &V() {return _v;};
  float V() const {return _v;};
};


typedef Volume<SimpleVoxel> MyVolume;

int main(int argc, char *argv[])
{
	MyVolume	volume;
  
  typedef vcg::tri::TrivialWalker<MyMesh,MyVolume>	MyWalker;
	typedef vcg::tri::MarchingCubes<MyMesh, MyWalker>	MyMarchingCubes;
	MyWalker walker;
	

  // Simple initialization of the volume with some cool perlin noise
	volume.Init(Point3i(64,64,64));
  for(int i=0;i<64;i++)
    for(int j=0;j<64;j++)
      for(int k=0;k<64;k++)
        volume.Val(i,j,k)=(j-32)*(j-32)+(k-32)*(k-32)  + i*10*math::Perlin::Noise(i*.2,j*.2,k*.2);


	// MARCHING CUBES
	MyMesh		mc_mesh;
	printf("[MARCHING CUBES] Building mesh...");
	MyMarchingCubes					mc(mc_mesh, walker);
	walker.BuildMesh<MyMarchingCubes>(mc_mesh, volume, mc, 20*20);
	vcg::tri::io::ExporterPLY<MyMesh>::Save( mc_mesh, "marching_cubes.ply");

	printf("OK!\n");
};