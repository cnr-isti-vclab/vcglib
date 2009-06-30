
namespace vcg{namespace tri{

/*
seeks and removes all doublets (a pair of quads sharing two consecutive edges)
by merging them into a single quad (thus removing one vertex and two tri faces)-
Returns number of removed Doublets
*/
template <class Mesh>
int BitQuadRemoveDoublets(Mesh &m)
{
  int res=0;
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::FaceType FaceType;
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) {
    fi->Q()=1;
    for (int k=0; k<3; k++) {
      if ( IsDoublet(*fi,k) ){
        res++;
        RemoveDoublet(*fi,k,m);
        if (fi->IsD()) break; // break wedge circle, if face disappeard
      }
    }
  }
  return res;
}

/*
marks (Quality=0) and approx. counts doublets (a pair of quads sharing two consecutive edges)
*/
template <class Mesh>
int BitQuadMarkDoublets(Mesh &m)
{
  int res=0;
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::FaceType FaceType;
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) {
    fi->Q()=1;
    for (int k=0; k<3; k++) {
      if ( IsDoublet(*fi,k) ){
        res++;
        if (fi->IsF((k+1)%3)) res++; // counts for a quad
        fi->Q()=0;
      }
    }
  }
  assert (res%2==0);
  return res/4; // return doublet pairs (approx, as a quad could be a part of many pairs)
}

/*
marks (Quality=0) and counts singlets (vertex B in an A-B-A-C quad)
*/
template <class Mesh>
int BitQuadMarkSinglets(Mesh &m)
{
  int res=0;
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::FaceType FaceType;
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) {
    fi->Q()=1;
    for (int k=0; k<3; k++) {
      if ( IsSinglet(*fi,k) ){
        res++;
        fi->Q()=0;
      }
    }
  }
  assert (res%2==0);
  return res/2; // return number of  singlet pairs
}

/*
deletes singlets, reutrns number of
*/
template <class Mesh>
int BitQuadRemoveSinglets(Mesh &m)
{
  int res=0;
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::FaceType FaceType;
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) {
    for (int k=0; k<3; k++) {
      if ( IsSinglet(*fi,k) ){
        res++;
        RemoveSinglet(*fi,k,m);
        return res;
        break;
      }
    }
  }
  return res; // return singlet pairs (approx, as a quad could be a part of many pairs)
}


/* returns average quad quality, and assigns it to triangle quality
*/
template <class Mesh>
typename Mesh::ScalarType MeasureBitQuadQuality(Mesh &m)
{
  assert(Mesh::HasPerFaceFlags());
  typename Mesh::ScalarType res = 0;
  int div = 0;
  typedef typename Mesh::FaceIterator FaceIterator;
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) {
    if (fi->IsAnyF()) {
      
      typename Mesh::ScalarType q = quadQuality( &*fi, FauxIndex(&*fi) );
      
      if (Mesh::HasPerFaceQuality()) fi->Q() = q;
      res += q;
      div++;
    }
  }
  if (!div) return 0; else return res / div;
}

}} // end namespace vcg::tri
