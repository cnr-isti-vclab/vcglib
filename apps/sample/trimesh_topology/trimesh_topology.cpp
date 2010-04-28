#include <vector>

#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/face/base.h>
#include<vcg/simplex/face/topology.h>
#include<vcg/complex/trimesh/base.h>
#include<vcg/complex/trimesh/create/platonic.h>

// topology computation
#include<vcg/complex/trimesh/update/topology.h>

// half edge iterators
#include<vcg/simplex/face/pos.h>



using namespace vcg;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>::AsVertexType,
																				Use<MyFace>::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::BitFlags  >{};
class MyFace    : public Face  < MyUsedTypes, face::VertexRef,face::FFAdj, face::Mark, face::BitFlags > {};
class MyMesh : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace > >{};


int main(int ,char ** ){

	MyMesh m;

	//generate a mesh
	vcg::tri::Icosahedron(m);

	//update the face-face topology 
	vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);

  // Now for each face the F() members are meaningful
  
  if(face::IsBorder(m.face[0],0)) printf("Edge 0 of face 0 is a border\n");
                       else printf("Edge 0 of face 0 is NOT a border\n"); // always this path!

  vcg::face::FFDetach<MyFace>(m.face[0],0);  // Detach the face [0] from the mesh
  vcg::face::FFDetach<MyFace>(m.face[0],1);
  vcg::face::FFDetach<MyFace>(m.face[0],2);

  if(face::IsBorder(m.face[0],0)) printf("Edge 0 of face 0 is a border\n"); // always this path!
                       else printf("Edge 0 of face 0 is NOT a border\n"); 

  m.face[0].SetD(); // deleting face [0] (i.e. marked as deleted) 


	// declare an iterator on the mesh
	vcg::face::Pos<MyMesh::FaceType> he, hei;

  // Now a simple search and trace of all the border of the mesh
  MyMesh::FaceIterator fi; 
  UnMarkAll(m);
  int BorderEdgeNum=0;
  int HoleNum=0;
  for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
			{	
				for(int j=0;j<3;j++)
				{
          if ( face::IsBorder(*fi,j) && tri::IsMarked(m,&*fi))
            {
              tri::Mark(m,&*fi);
              hei.Set(&*fi,j,fi->V(j));
							he=hei;
							do
							{
                BorderEdgeNum++;	
								he.NextB(); // next edge along a border 
                tri::Mark(m,he.f);
       				}
							while (he.f!=hei.f);
							HoleNum++;
            }			
				}
			}
		
  printf("Mesh has %i holes and %i border edges\n",HoleNum,BorderEdgeNum);
  return 0;
}

