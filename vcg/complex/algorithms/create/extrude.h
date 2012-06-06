#ifndef EXTRUDE_H
#define EXTRUDE_H

namespace vcg {
namespace tri {


template <class MeshType> class Extrude
{
  public:
  typedef typename MeshType::FacePointer FacePointer;
  typedef typename MeshType::EdgePointer EdgePointer;
  typedef typename MeshType::VertexPointer VertexPointer;
  typedef typename MeshType::FaceType FaceType;
  typedef typename MeshType::EdgeType EdgeType;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::EdgeIterator EdgeIterator;
  typedef typename MeshType::FaceIterator FaceIterator;

static void ProfileWithCap(MeshType &profile, MeshType &surface, const vcg::Similarityf &sim )
{
  surface.Clear();

  for(VertexIterator vi=profile.vert.begin();vi!=profile.vert.end();++vi)
  {
    VertexIterator vp=tri::Allocator<MeshType>::AddVertices(surface,2);
    vp->P()=vi->P();
    ++vp;
    vp->P()= sim*vi->P() ;
  }

  for(EdgeIterator ei=profile.edge.begin();ei!=profile.edge.end();++ei)
  {
    int i0=tri::Index(profile,ei->V(0));
    int i1=tri::Index(profile,ei->V(1));

    FaceIterator fp= tri::Allocator<MeshType>::AddFaces(surface,2);
    fp->V(0) = &surface.vert[i0*2];
    fp->V(1) = &surface.vert[i1*2];
    fp->V(2) = &surface.vert[i0*2+1];
    ++fp;
    fp->V(0) = &surface.vert[i1*2+1];
    fp->V(1) = &surface.vert[i0*2+1];
    fp->V(2) = &surface.vert[i1*2];
  }

  MeshType cap1;

  tri::CapEdgeMesh(profile,cap1);
  if(cap1.fn==0) CapEdgeMesh(profile,cap1,true);
  tri::Append<MeshType,MeshType>::Mesh(surface,cap1);

  for(VertexIterator vi=cap1.vert.begin();vi!=cap1.vert.end();++vi)
    vi->P() = sim*vi->P();
  tri::Append<MeshType,MeshType>::Mesh(surface,cap1);

  tri::Clean<MeshType>::RemoveDuplicateVertex(surface);
  bool oriented,orientable;
  tri::UpdateTopology<MeshType>::FaceFace(surface);
  tri::Clean<MeshType>::OrientCoherentlyMesh(surface,oriented,orientable);
}

static void ProfileWithCap(MeshType &profile, MeshType &surface, const Point3f offset)
{
  Similarityf tra;
  tra.SetTranslate(offset);
  ProfileWithCap(profile,surface,tra);
}

}; // end class

} // end namespace tri
} // end namespace vcg

#endif // EXTRUDE_H
