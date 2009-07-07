
namespace vcg{namespace tri{

// helper function: mark a quadface, setting Q at 0, and neight at .75, 0.5...
template <class Mesh>
void MarkFace(typename Mesh::FaceType* f, Mesh &m){
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::FaceType FaceType;

  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) {
     fi->Q() = 1; 
  }
  
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) f->FFp(i)->FFp(j)->Q() = 0.75;
  }
  for (int i=0; i<3; i++) {
    f->FFp(i)->Q() = 0.50;
  }
  f->Q() = 0;
  
}

// helper function: mark a quadface, setting Q at 0, and neight at .75, 0.5...
template <class Mesh>
void MarkVertex(typename Mesh::FaceType* f, int wedge, Mesh &m){
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::FaceType FaceType;
  typedef typename Mesh::VertexType VertexType;
  VertexType *v = f->V(wedge);
  
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) {
     if (fi->V0(0)==v || fi->V1(0)==v ||fi->V2(0)==v ) fi->Q() = 0; 
     // else fi->Q() = 1; 
  }
  
}

template <class Mesh>
bool MarkSmallestEdge(Mesh &m, bool perform)
{
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::FaceType FaceType;
  typedef typename Mesh::ScalarType ScalarType;
  ScalarType min = std::numeric_limits<ScalarType>::max();
  
  FaceType *fa=NULL; int w=0;
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD())
  for (int k=0; k<3; k++) {
    FaceType *f=&*fi;
    
    if (f->IsF(k)) continue;
    if (f->FFp(k) == f ) continue; // skip borders
    
    ScalarType score;

    score = (f->P0(k) - f->P1(k)).Norm();
    if (score<min) {
      min=score; 
      fa = f;
      w = k;
    }
        
  }
  if (fa) {
    if (perform) {
      return CollapseQuadEdge(*fa,w,m);
    } else {
      fa->Q()=0.0;
      fa->FFp(w)->Q()=0.0;
      return true;
    }
  }
  return false;
}

// returns: 0 if fail. 1 if edge. 2 if diag.
template <class Mesh>
int MarkSmallestEdgeOrDiag(Mesh &m, typename Mesh::ScalarType edgeMult, bool perform)
{
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::FaceType FaceType;
  typedef typename Mesh::ScalarType ScalarType;
  ScalarType min = std::numeric_limits<ScalarType>::max();
  
  FaceType *fa=NULL; int w=0; bool counterDiag = false;
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD())
  for (int k=0; k<3; k++) {
    FaceType *f=&*fi;
    
    if (f->FFp(k) >= f ) continue; // skip borders (==), and do it one per edge
    
    ScalarType score;

    score = (f->P0(k) - f->P1(k)).Norm();
    
    if (!f->IsF(k)) score*=edgeMult; // edges are supposed to be smaller!
    
    if (score<min) {
      min=score; 
      fa = f;
      w = k;
      counterDiag=false;
    }
    
    if (f->IsF(k)) { // for diag faces, test counterdiag too
      score = CounterDiag(f).Norm();
      if (score<min) {
        min=score; 
        fa = f;
        w = k;
        counterDiag=true;
      }
    }
    

        
  }
  if (fa) {
    if (perform) {
      if (fa->IsF(w)) {
        if (counterDiag) {
          CollapseQuadCounterDiag(*fa, PosOnDiag(*fa,true), m ); return 2;
        } else {
          CollapseQuadDiag(*fa, PosOnDiag(*fa,false), m ); return 2;
        }
      } else {
        if  (CollapseQuadEdge(*fa,w,m)) return 1;
      }
    } else {
      fa->Q()=0.0;
      fa->FFp(w)->Q()=0.0;
      if (fa->IsF(w)) return 2; else return 1;
    }
  }
  return 0;
}


template <class Mesh>
void MarkSmallestDiag(Mesh &m)
{
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::FaceType FaceType;
  typedef typename Mesh::ScalarType ScalarType;
  ScalarType min = std::numeric_limits<ScalarType>::max();
  
  FaceType *fa=NULL; 
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) {
    FaceType *f=&*fi;
    
    ScalarType score;

    score = Diag(f).Norm();
    if (score<min) {
      min=score; 
      fa = f;
    }
    
    score = CounterDiag(f).Norm();
    if (score<min) {
      min=score; 
      fa = f;
    }
    
  }
  if (fa) {
    fa->Q()=0.0;
    fa->FFp(FauxIndex(fa))->Q()=0.0;
  }

}

template <class Mesh>
bool IdentifyAndCollapseSmallestDiag(Mesh &m){
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::FaceType FaceType;
  typedef typename Mesh::ScalarType ScalarType;
  ScalarType min = std::numeric_limits<ScalarType>::max();
  
  FaceType *fa=NULL; bool flip;
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) {
    FaceType *f=&*fi;
    
    ScalarType score;
    
    score = Diag(f).Norm();
    if (score<min) {
      min=score; 
      fa = f;
      flip = false;
    }
    
    score = CounterDiag(f).Norm();
    if (score<min) {
      min=score; 
      fa = f;
      flip = true;
    }
    
  }
  if (!fa) return false;
  
  if (TestAndRemoveDoublet(*fa,0,m)) { return true; }
  if (TestAndRemoveDoublet(*fa,1,m)) { return true; }
  if (TestAndRemoveDoublet(*fa,2,m)) { return true; }
  int k = FauxIndex(fa);
  if (TestAndRemoveDoublet( *fa->FFp(k),(fa->FFi(k)+2)%3, m )) return true;

  if (flip) {
    if (!CheckFlipBitQuadDiag(*fa) ) {
      // I can't collapse (why?)
      MarkFace(fa,m);
      return false;
    } else 
    CollapseQuadCounterDiag(*fa, PosOnDiag(*fa,true), m );
  }
  else  {
    CollapseQuadDiag(*fa, PosOnDiag(*fa,false), m );
  }
  return true;
}



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
marks (Quality=0) and approx. counts profitable vertex rotations
(vertex rotations which make edge shorter
*/
template <class Mesh, bool perform>
int BitQuadMarkVertexRotations(Mesh &m)
{
  int res=0;
  typedef typename Mesh::VertexIterator VertexIterator;
  typedef typename Mesh::VertexType VertexType;
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::FaceType FaceType;
  for (VertexIterator vi = m.vert.begin();  vi!=m.vert.end(); vi++) if (!vi->IsD()) vi->ClearV();
  if (!perform)
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) fi->Q()=1.0;

  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) {
    
    for (int k=0; k<3; k++) {
      if (fi->V(k)->IsV()) continue;
      if (TestBitQuadVertexRotation(*fi,k)) {
        res++;
        fi->V(k)->SetV();
        if (!perform)
          MarkVertex(&*fi, k, m); //fi->Q()=0;
        else {
          RotateBitQuadVertex(*fi, k); //fi->Q()=0;
          return 1;
        }
      }
    }
  }
  return res;
}

// mark (and count) all edges that are worth rotating
// if perform == true, actually rotate them
template <class Mesh, bool perform>
int BitQuadMarkEdgeRotations(Mesh &m)
{
  int count = 0;
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::FaceType FaceType;
  
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) fi->Q()=1;
  
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) {
    //if (count>0) break;
    for (int k=0; k<3; k++) {
      if (fi->IsF(k)) continue;
      if (fi->FFp(k)<= &*fi) continue; // only once per real (non faux) edge, and only for non border ones
      int best = TestBitQuadEdgeRotation(*fi, k);
      if (perform) {
        if (best==+1) if (RotateBitQuadEdge<FaceType, true>(*fi, k)) count++;
        if (best==-1) if (RotateBitQuadEdge<FaceType,false>(*fi, k)) count++;
      }
      else {
        if (best!=0) { fi->Q()=0; fi->FFp(k)->Q()=0; count++; }
      }
    }
  }

  return count;
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
