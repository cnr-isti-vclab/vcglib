
#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/vertex/component.h>
#include<vcg/simplex/face/base.h>
#include<vcg/simplex/face/component.h>
#include<vcg/complex/trimesh/base.h>
#include<vcg/complex/trimesh/allocate.h>


class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>		::AsVertexType,
											vcg::Use<MyFace>			::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f,vcg::vertex::Normal3f>{};
class MyFace    : public vcg::Face< MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f> {};

class MyMesh : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

 
float Irradiance(MyMesh::VertexType v){
	// .....
	return 1.0;
}

int main()
{


  MyMesh m;
  //...here m is filled
  
  // add a per-vertex attribute with type float named "Irradiance"
  MyMesh::PerVertexAttributeHandle<float> ih = vcg::tri::Allocator<MyMesh>::AddPerVertexAttribute<float> (m,std::string("Irradiance"));

  // add a per-vertex attribute with type float named "Radiosity"   
  vcg::tri::Allocator<MyMesh>::AddPerVertexAttribute<float> (m,std::string("Radiosity"));
 
  // add a per-vertex attribute with type bool and no name specified
  MyMesh::PerVertexAttributeHandle<bool> blocked_h = vcg::tri::Allocator<MyMesh>::AddPerVertexAttribute<bool> (m); 
  
  // add a per-vertex attribute with type bool and no name specified
  MyMesh::PerFaceAttributeHandle<bool> blocked_hf = vcg::tri::Allocator<MyMesh>::AddPerFaceAttribute<bool> (m); 

  MyMesh::VertexIterator vi; int i = 0;
  for(vi   = m.vert.begin(); vi != m.vert.end(); ++vi,++i){
   ih[vi]  = Irradiance(*vi);  // [] operator takes a iterator
   ih[*vi] = Irradiance(*vi);  // or a MyMesh::VertexType object
   ih[&*vi]= Irradiance(*vi);  // or a pointer to it
   ih[i]   = Irradiance(*vi);  // or an integer index
  }
    
  // Once created with AddPerVertexAttribute, an handle to the attribute can be obtained as follows
  MyMesh::PerVertexAttributeHandle<float> rh = vcg::tri::Allocator<MyMesh>::GetPerVertexAttribute<float>(m,"Radiosity");

  // you can query if an attribute is present or not
  bool hasRadiosity = vcg::tri::HasPerVertexAttribute(m,"Radiosity");

  // you can delete an attibute by name
  vcg::tri::Allocator<MyMesh>::DeletePerVertexAttribute(m,"Radiosity");

  // you can delete an attibute by handle
  vcg::tri::Allocator<MyMesh>::DeletePerVertexAttribute(m,blocked_h);

  bool res ;
 
  res = vcg::tri::Allocator<MyMesh>::IsValidHandle(m,ih);printf("%d\n",res);
  res = vcg::tri::Allocator<MyMesh>::IsValidHandle(m,blocked_hf);printf("%d\n",res);
  vcg::tri::Allocator<MyMesh>::DeletePerVertexAttribute(m,ih);
  vcg::tri::Allocator<MyMesh>::DeletePerFaceAttribute(m,blocked_hf);
  res = vcg::tri::Allocator<MyMesh>::IsValidHandle(m,ih);printf("%d\n",res);
  res = vcg::tri::Allocator<MyMesh>::IsValidHandle(m,blocked_hf);printf("%d\n",res);

 
}
