#include <vector>
#include <vcg/complex/trimesh/subset.h>
#include <vcg/simplex/face/jumping_pos.h>
#include <vcg/space/planar_polygon_tessellation.h>

/** BIT-QUAD creation support:
    a few basic operations to work with bit-quads simplices
    (quads defined by faux edges over a tri mesh backbone)


   [ basic operations: ]

   bool IsDoublet(const FaceType& f, int wedge)
   void RemoveDoublet(FaceType &f, int wedge, MeshType& m)
    - identifies and removed "Doublets" (pair of quads sharing two consecutive edges)

   bool IsSinglet(const FaceType& f, int wedge)
   void RemoveSinglet(FaceType &f, int wedge, MeshType& m)

   void FlipDiag(FaceType &f)
    - rotates the faux edge of a quad (quad only change internally)

   bool RotateEdge(FaceType& f, int w0a);
    - rotate a quad edge (clockwise or counterclockwise, specified via template)
   
   bool RotateVertex(FaceType &f, int w0)
    - rotate around a quad vertex ("wind-mill" operation)

   void CollapseDiag(FaceType &f, ... p , MeshType& m)
    - collapses a quad on its diagonal. 
    - p identifies the pos of collapsed point 
      (as either the parametric pos on the diagonal, or a fresh coordtype)
     

   [ helper functions: ]

   ScalarType quadQuality( ... ); 
    - returns the quality for a given quad 
    - (should be made into a template parameter for methods using it)
    - currently measures how squared each angle is
   
   int FauxIndex(const FaceType* f);
    - returns index of the only faux edge of a quad (otherwise, assert)

   int CountBitPolygonInternalValency(const FaceType& f, int wedge)
    - returns valency of vertex in terms  of polygons (quads, tris...)
    
  
*/

// this must become a parameter in the corresponding class
#define DELETE_VERTICES 0
// let's not remove them after all...
// since TwoManyfold is weak, the vertex could still be used elsewhere...

namespace vcg{namespace tri{

/* simple geometric-interpolation mono-function class used 
as a default template parameter to BitQuad class */
template <class VertexType>
class GeometricInterpolator{
public:
  typedef typename VertexType::ScalarType ScalarType;
  static void Apply( const VertexType &a,  const VertexType &b, ScalarType t, VertexType &res){
    assert (&a != &b);
    res.P() = a.P()*(1-t) + b.P()*(t);
  }
};

template <
  // first template parameter: the tri mesh (with face-edges flagged)
  class _MeshType, 
  // second template parameter: used to define interpolations between points
  class Interpolator = GeometricInterpolator<typename _MeshType::VertexType> 
>
class BitQuad{
public:

typedef _MeshType MeshType;
typedef typename MeshType::ScalarType ScalarType;
typedef typename MeshType::CoordType CoordType;
typedef typename MeshType::FaceType FaceType;
typedef typename MeshType::FaceType* FaceTypeP;
typedef typename MeshType::VertexType VertexType;
typedef typename MeshType::FaceIterator FaceIterator;
typedef typename MeshType::VertexIterator VertexIterator;

static void MarkFaceF(FaceType *f){
  f->V(0)->SetS();
  f->V(1)->SetS();
  f->V(2)->SetS();
  int i=FauxIndex(f);
  f->FFp( i )->V2( f->FFi(i) )->SetS();

}


template <bool verse>
static bool RotateEdge(FaceType& f, int w0a){
  FaceType *fa = &f;
  assert(! fa->IsF(w0a) );

  VertexType *v0, *v1;
  v0= fa->V0(w0a);
  v1= fa->V1(w0a);
  
  int w1a = (w0a+1)%3;
  int w2a = (w0a+2)%3;

  FaceType *fb = fa->FFp(w0a);
  
  MarkFaceF(fa);
  MarkFaceF(fb);
  
  int w0b = fa->FFi(w0a);
  int w1b = (w0b+1)%3;
  int w2b = (w0b+2)%3;
    
  if (fa->IsF(w2a) == verse) {
    if (!CheckFlipDiag(*fa)) return false;
    FlipDiag(*fa);
    // recover edge index, so that (f, w0a) identifies the same edge as before
    FaceType *fc = fa->FFp(FauxIndex(fa));
    for (int i=0; i<3; i++){
      if ( v0==fa->V0(i) && v1==fa->V1(i) ) w0a = i;
      if ( v0==fc->V0(i) && v1==fc->V1(i) ) { fa = fc; w0a = i; }
    }
  }
  
  if (fb->IsF(w2b) == verse) {
    if (!CheckFlipDiag(*fb)) return false;
    FlipDiag(*fb);
  }
  
  if (!CheckFlipEdge(*fa,w0a)) return false;
  FlipEdge(*fa,w0a);
  return true;
}

/* small helper function which returns the index of the only
   faux index, assuming there is exactly one (asserts out otherwise)
*/
static int FauxIndex(const FaceType* f){
  if (f->IsF(0)) return 0;
  if (f->IsF(1)) return 1;
  assert(f->IsF(2));
  return 2;
}

// rotates the diagonal of a quad
static void FlipDiag(FaceType &f){
  int faux = FauxIndex(&f);
  FaceType* fa = &f;
  FaceType* fb = f.FFp(faux);
  vcg::face::FlipEdge(f, faux);
  // ripristinate faux flags
  fb->ClearAllF();
  fa->ClearAllF();
  for (int k=0; k<3; k++) {
    if (fa->FFp(k) == fb) fa->SetF(k);
    if (fb->FFp(k) == fa) fb->SetF(k);
  }
}


// given a vertex (i.e. a face and a wedge), 
// this function tells us how the average edge lenght around a vertex would change
// if that vertex is rotated
static ScalarType AvgEdgeLenghtVariationIfVertexRotated(const FaceType &f, int w0)
{
  assert(!f.IsD());
  
  ScalarType 
    before=0, // sum of quad edges (originating from v)
    after=0;  // sum of quad diag (orginating from v)
  int guard = 0;

  // rotate arond vertex
  const FaceType* pf = &f;
  int pi = w0;
  int n = 0; // vertex valency
  int na = 0; 
  do {
    ScalarType triEdge = (pf->P0(pi) - pf->P1(pi) ).Norm();
    if (pf->IsF(pi)) { after += triEdge; na++;}
    else { before+= triEdge; n++; }
    if ( pf->IsF((pi+1)%3)) { after += CounterDiag( pf ).Norm(); na++; }
    
    const FaceType *t = pf;
    t = pf->FFp( pi );
    if (pf == t ) return std::numeric_limits<ScalarType>::max(); // it's a mesh border! flee!
    pi = pf->cFFi( pi );
    pi = (pi+1)%3; // FaceType::Next( pf->FFi( pi ) );
    pf = t;
    assert(guard++<100);
  } while (pf != &f);
  assert (na == n);
  return (after-before)/n;
}

/*
  const FaceType* pf = &f;
  int pi = wedge;
  int res = 0, guard=0;
  do {
    if (!pf->IsAnyF()) return false; // there's a triangle!
    if (!pf->IsF(pi)) res++;
    const FaceType *t = pf;
    t = pf->FFp( pi );
    if (pf == t ) return false;
    pi = pf->cFFi( pi );
    pi = (pi+1)%3; // FaceType::Next( pf->FFi( pi ) );
    pf = t;
    assert(guard++<100);
  } while (pf != &f);
*/

// given a vertex (i.e. a face and a wedge), 
// this function tells us if it should be rotated or not
// (currently, we should iff it is shortened)
static bool TestVertexRotation(const FaceType &f, int w0)
{
  assert(!f.IsD());
  // rotate quad IFF this way edges become shorter:
  return AvgEdgeLenghtVariationIfVertexRotated(f,w0)<0;
}


static bool RotateVertex(FaceType &f, int w0)
{
  
  int guard = 0;

  FaceType* pf = &f;
  int pi = w0;
  int n = 0; // vertex valency
  
  if (pf->IsF((pi+2) % 3)) { 
      pi = (pi+2)%3;
      // do one step back
      int tmp = pf->FFi(pi); pf = pf->FFp(pi); pi = tmp;  // flipF 
  }
  
  const FaceType* stopA = pf;
  const FaceType* stopB = pf->FFp(FauxIndex(pf));

  // rotate around vertex, flipping diagonals if necessary,
  do {
    bool mustFlip;
    if (pf->IsF(pi)) { 
      // if next edge is faux, move on other side of quad
      int tmp = (pf->FFi(pi)+1)%3; pf = pf->FFp(pi); pi = tmp;  // flipF 
      mustFlip = false;
    }
    else {
      mustFlip = true;
    }

    FaceType *lastF = pf;
    
    int tmp = (pf->FFi(pi)+1)%3; pf = pf->FFp(pi); pi = tmp;  // flipF
        
    if (mustFlip) {
      if (!CheckFlipDiag(*lastF)) return false; // cannot flip??
      FlipDiag(*lastF);
    }
    MarkFaceF(pf);

    
  } while (pf != stopA && pf!= stopB);
  
  // last pass: rotate arund vertex again, changing faux status
  stopA=pf;
  do {
    int j = pi;
    if (pf->IsF(j)) 
      { pf->ClearF(j); IncreaseValency(pf->V(j)); } 
    else 
      { pf->SetF(j); DecreaseValency(pf->V(j)); }

    j = (j+2)%3;
    if (pf->IsF(j)) pf->ClearF(j); else pf->SetF(j);
    int tmp = (pf->FFi(pi)+1)%3; pf = pf->FFp(pi); pi = tmp;  // flipF flipV
  } while (pf != stopA );

  return true;
}





// flips the faux edge of a quad
static void FlipEdge(FaceType &f, int k){
  assert(!f.IsF(k));
  FaceType* fa = &f;
  FaceType* fb = f.FFp(k);
  assert(fa!=fb); // else, rotating a border edge

  // backup prev other-quads-halves
  FaceType* fa2 = fa->FFp( FauxIndex(fa) );
  FaceType* fb2 = fb->FFp( FauxIndex(fb) );

  IncreaseValency( fa->V2(k) );
  IncreaseValency( fb->V2(f.FFi(k)) );
  DecreaseValency( fa->V0(k) );
  DecreaseValency( fa->V1(k) );

  vcg::face::FlipEdge(*fa, k);
  
  // ripristinate faux flags
  fb->ClearAllF();
  fa->ClearAllF();
  for (int k=0; k<3; k++) {
    //if (fa->FFp(k) == fa2) fa->SetF(k);
    //if (fb->FFp(k) == fb2) fb->SetF(k);
    if (fa->FFp(k)->IsF( fa->FFi(k) )) fa->SetF(k);
    if (fb->FFp(k)->IsF( fb->FFi(k) )) fb->SetF(k);
  }
  
}

// check if a quad diagonal can be topologically flipped
static bool CheckFlipDiag(FaceType &f){
  return (vcg::face::CheckFlipEdge(f, FauxIndex(&f) ) );
}

// given a face (part of a quad), returns its diagonal
static CoordType Diag(const FaceType* f){
  int i = FauxIndex(f);
  return f->P1( i ) - f->P0( i );
}


// given a face (part of a quad), returns other diagonal
static CoordType CounterDiag(const FaceType* f){
  int i = FauxIndex(f);
  return f->cP2( i ) - f->cFFp( i )->cP2(f->cFFi(i) ) ;
}

/* helper function:
   collapses a single face along its faux edge. 
   Updates FF adj of other edges. */
static void _CollapseDiagHalf(FaceType &f, int faux, MeshType& m)
{
  int faux1 = (faux+1)%3;
  int faux2 = (faux+2)%3;
  
  DecreaseValency(f.V2(faux)); // update valency 
  
  FaceType* fA = f.FFp( faux1 );
  FaceType* fB = f.FFp( faux2 );
  
  MarkFaceF(fA);
  MarkFaceF(fB);

  int iA = f.FFi( faux1 );
  int iB = f.FFi( faux2 );
  
  if (fA==&f && fB==&f) { 
    // both non-faux edges are borders: tri-face disappears, just remove the vertex
    if (DELETE_VERTICES)
    if (GetValency(f.V(faux2))==0)
      Allocator<MeshType>::DeleteVertex(m,*(f.V(faux2))); 
  } else {
    if (fA==&f) {
      fB->FFp(iB) = fB;  fB->FFi(iB) = iB;
    } else {
      fB->FFp(iB) = fA;  fB->FFi(iB) = iA;
    }

    if (fB==&f) {
      fA->FFp(iA) = fA;  fA->FFi(iA) = iA;
    } else {
      fA->FFp(iA) = fB;  fA->FFi(iA) = iB;
    }
  }
  
  Allocator<MeshType>::DeleteFace(m,f);

}

static void RemoveDoublet(FaceType &f, int wedge, MeshType& m){
  if (f.IsF((wedge+1)%3) ) {
    VertexType *v = f.V(wedge);
    FlipDiag(f);
    // quick hack: recover wedge index after flip
    if (f.V(0)==v) wedge = 0;
    else if (f.V(1)==v) wedge = 1;
    else {
      assert(f.V(2)==v);
      wedge = 2;
    }
  }
  ScalarType k=(f.IsF(wedge))?1:0;
  CollapseDiag(f, k, m);
  VertexType *v = f.V(wedge);
}

static void RemoveSinglet(FaceType &f, int wedge, MeshType& m){
  FaceType *fa, *fb; // these will die
  FaceType *fc, *fd; // their former neight
  fa = & f;
  fb = fa->FFp(wedge);
  int wa0 = wedge;
  int wa1 = (wa0+1)%3 ;
  int wa2 = (wa0+2)%3 ;
  int wb0 = (fa->FFi(wa0)+1)%3;
  int wb1 = (wb0+1)%3 ; 
  int wb2 = (wb0+2)%3 ;
  assert (fb == fa->FFp( wa2 ) ); // otherwise, not a singlet
  
  // valency decrease
  DecreaseValency(fa->V(wa1));
  DecreaseValency(fa->V(wa2));
  if (fa->IsF(wa0)) {
    DecreaseValency(fa->V(wa2)); // double decrease of valency
  } else {
    DecreaseValency(fa->V(wa1)); // double decrease of valency
  }
  
  // no need to MarkFaceF !
  
  fc = fa->FFp(wa1);
  fd = fb->FFp(wb1);
  int wc = fa->FFi(wa1);
  int wd = fb->FFi(wb1);
  fc->FFp(wc) = fd;
  fc->FFi(wc) = wd;
  fd->FFp(wd) = fc;
  fd->FFi(wd) = wc;
  // faux status of survivors: unchanged
  assert( ! ( fc->IsF( wc) ) );
  assert( ! ( fd->IsF( wd) ) );
  
  Allocator<MeshType>::DeleteFace( m,*fa );
  Allocator<MeshType>::DeleteFace( m,*fb );
  
  DecreaseValency(fa->V(wedge) );
  if (DELETE_VERTICES) 
  if (GetValency(fa->V(wedge))==0) Allocator<MeshType>::DeleteVertex( m,*fa->V(wedge) );
}


static bool TestAndRemoveDoublet(FaceType &f, int wedge, MeshType& m){
  if (IsDoublet(f,wedge)) {
     RemoveDoublet(f,wedge,m);
     return true;
  }
  return false;
}

static bool TestAndRemoveSinglet(FaceType &f, int wedge, MeshType& m){
  if (IsSinglet(f,wedge)) {
     RemoveSinglet(f,wedge,m);
     return true;
  }
  return false;
}

// given a face and a wedge, counts its valency in terms of quads (and triangles)
// uses only FF, assumes twomanyfold
// returns -1 if border
static int CountBitPolygonInternalValency(const FaceType& f, int wedge){
  const FaceType* pf = &f;
  int pi = wedge;
  int res = 0;
  do {
    if (!pf->IsF(pi)) res++;
    const FaceType *t = pf;
    t = pf->FFp( pi );
    if (pf == t ) return -1;
    pi = (pi+1)%3; // FaceType::Next( pf->FFi( pi ) );
    pf = t;
  } while (pf != &f);
  return res;
}

// given a face and a wedge, returns if it host a doubet
// assumes tri and quad only. uses FF topology only.
static bool IsDoubletFF(const FaceType& f, int wedge){
  const FaceType* pf = &f;
  int pi = wedge;
  int res = 0, guard=0;
  do {
    if (!pf->IsAnyF()) return false; // there's a triangle!
    if (!pf->IsF(pi)) res++;
    const FaceType *t = pf;
    t = pf->FFp( pi );
    if (pf == t ) return false;
    pi = pf->cFFi( pi );
    pi = (pi+1)%3; // FaceType::Next( pf->FFi( pi ) );
    pf = t;
    assert(guard++<100);
  } while (pf != &f);
  return (res == 2);
}

// version that uses vertex valency
static bool IsDoublet(const FaceType& f, int wedge){
  return (GetValency( f.V(wedge)) == 2) && (!f.V(wedge)->IsB() ) ;
}

// given a face and a wedge, returns if it host a singlets
// assumes tri and quad only. uses FF topology only.
static bool IsSingletFF(const FaceType& f, int wedge){
  const FaceType* pf = &f;
  int pi = wedge;
  int res = 0, guard=0;
  do {
    if (!pf->IsAnyF()) return false; // there's a triangle!
    if (!pf->IsF(pi)) res++;
    const FaceType *t = pf;
    t = pf->FFp( pi );
    if (pf == t ) return false;
    pi = pf->cFFi( pi );
    pi = (pi+1)%3; // FaceType::Next( pf->FFi( pi ) );
    pf = t;
    assert(guard++<100);
  } while (pf != &f);
  return (res == 1);
}

// version that uses vertex valency
static bool IsSinglet(const FaceType& f, int wedge){
  return (GetValency( f.V(wedge) ) == 1) && (!f.V(wedge)->IsB() ) ;
}

static bool CollapseEdgeDirect(FaceType &f, int w0, MeshType& m){
  FaceType * f0 = &f;

  assert( !f0->IsF(w0) );
  
  VertexType *v0, *v1;
  v0 = f0->V0(w0);
  v1 = f0->V1(w0);
  
  if (!RotateVertex(*f0,w0)) return false;

  // quick hack: recover original wedge
  if      (f0->V(0) == v0) w0 = 0;
  else if (f0->V(1) == v0) w0 = 1;
  else if (f0->V(2) == v0) w0 = 2; 
  else assert(0);
  
  assert( f0->V1(w0) == v1 );
  assert( f0->IsF(w0) );
  
  CollapseDiag(*f0,PosOnDiag(*f0,false), m);
  return true;
}

static bool CollapseEdge(FaceType &f, int w0, MeshType& m){
  FaceTypeP f0 = &f;
  
  assert(!f0->IsF(w0)); // don't use this method to collapse diag.
  
  FaceTypeP f1 = f0->FFp(w0); 
  int w1 = f0->FFi(w0); 
  
  assert(f0!=f1); // can't collapse border edges!
  
  // choose: rotate around V0 or around V1?
  if (
    AvgEdgeLenghtVariationIfVertexRotated(*f0,w0)
    <
    AvgEdgeLenghtVariationIfVertexRotated(*f1,w1)
  )    return CollapseEdgeDirect(*f0,w0,m); 
  else return CollapseEdgeDirect(*f1,w1,m);
}



/** collapses a quad diagonal a-b 
  forming the new vertex in between the two old vertices.
   if k == 0, new vertex is in a
   if k == 1, new vertex is in b
   if k == 0.5, new vertex in the middle, etc
*/

static void CollapseCounterDiag(FaceType &f, ScalarType interpol, MeshType& m){
  //CoordType p;
  //int fauxa = FauxIndex(&f);
  //p = f.V(fauxa)->P()*(1-k) + f.V( (fauxa+1)%3 )->P()*(k);
  
  FlipDiag(f);
  CollapseDiag(f,interpol,m);
}

/*
static void CollapseCounterDiag(FaceType &f, ScalarType k, MeshType& m){
  CoordType p;
  int fauxa = FauxIndex(&f);
  p = f.P2(fauxa)*(1-k) + f.FFp( fauxa )->P2( f.FFi( fauxa ) )*(k);
  CollapseCounterDiag(f,p,m);
}
*/

//static void CollapseCounterDiag(FaceType &f, const CoordType &p, MeshType& m){
//  FlipDiag(f);
//  CollapseDiag(f,p,m);
//}


static void CollapseDiag(FaceType &f, ScalarType interpol, MeshType& m){
    
  FaceType* fa = &f;
  int fauxa = FauxIndex(fa);
  FaceType* fb = fa->FFp(fauxa);
  assert (fb!=fa);
  int fauxb = FauxIndex(fb);
  
  VertexType* va = fa->V(fauxa); // va lives
  VertexType* vb = fb->V(fauxb); // vb dies
  
  Interpolator::Apply( *(f.V0(fauxa)), *(f.V1(fauxa)), interpol, *va);

  // update FV...
  bool border = false;
  int pi = fauxb;
  FaceType* pf = fb; /* pf, pi could be a Pos<FaceType> p(pf, pi) */
  // rotate around vb, (same-sense-as-face)-wise
  do {
    pf->V(pi) = va;
    
    pi=(pi+2)%3;
    FaceType *t = pf->FFp(pi);
    if (t==pf) { border= true; break; }
    pi = pf->FFi(pi);
    pf = t;
  } while (pf!=fb);
  
  // rotate around va, (counter-sense-as-face)-wise
  if (border) {
    int pi = fauxa;
    FaceType* pf = fa; /* pf, pi could be a Pos<FaceType> p(pf, pi) */
    do {
      pi=(pi+1)%3;
      pf->V(pi) = va;
      FaceType *t = pf->FFp(pi);
      if (t==pf) break;
      pi = pf->FFi(pi);
      pf = t;
    } while (pf!=fb);
  }
    
  // update FF, delete faces
  _CollapseDiagHalf(*fb, fauxb, m);
  _CollapseDiagHalf(*fa, fauxa, m);
  
  SetValency(va, GetValency(vb)+GetValency(va)-2);
  SetValency(vb, GetValency(vb)+GetValency(va)-2);
  if (DELETE_VERTICES) Allocator<MeshType>::DeleteVertex(m,*vb);
  
}      




// helper function: find a good position on a diag to collapse a point
// currently, it is point in the middle, 
//    unless a mixed border-non border edge is collapsed, then it is an exreme
static ScalarType PosOnDiag(const FaceType& f, bool counterDiag){
  bool b0, b1, b2, b3; // which side of the quads are border
  
  const FaceType* fa=&f;
  int ia = FauxIndex(fa);
  const FaceType* fb=fa->cFFp(ia); 
  int ib = fa->cFFi(ia);
  
  b0 = fa->FFp((ia+1)%3) == fa;
  b1 = fa->FFp((ia+2)%3) == fa;
  b2 = fb->FFp((ib+1)%3) == fb;
  b3 = fb->FFp((ib+2)%3) == fb;
  
  if (counterDiag) {
    if (  (b0||b1) && !(b2||b3) ) return 1;
    if ( !(b0||b1) &&  (b2||b3) ) return 0;
  } else {
    if (  (b1||b2) && !(b3||b0) ) return 0;
    if ( !(b1||b2) &&  (b3||b0) ) return 1;
  }
  //if (f->FF( FauxIndex(f) )->IsB(
  return 0.5f;
}

// trick! hide valency in flags
typedef enum { VALENCY_FLAGS = 24 } ___; // this bit and the 4 successive one are devoted to store valency

static void SetValency(VertexType *v, int n){
  //v->Q() = n;
  assert(n>=0 && n<=31);
  v->Flags()&= ~(31<<VALENCY_FLAGS);
  v->Flags()|= n<<VALENCY_FLAGS;
}

static int GetValency(const VertexType *v){
  //return (int)(v->cQ());
  return ( v->Flags() >> (VALENCY_FLAGS) ) & 31;
}

static void IncreaseValency(VertexType *v, int dv=1){
#ifdef NDEBUG
  v->Flags() += dv<<VALENCY_FLAGS;
#else
  SetValency( v, GetValency(v)+dv );
#endif
}

static int DecreaseValency(VertexType *v, int dv=1){
#ifdef NDEBUG
  v->Flags() -= dv<<VALENCY_FLAGS;
#else
  SetValency( v, GetValency(v)-dv );
#endif
}

static void UpdateValencyInFlags(MeshType& m){
  for (VertexIterator vi = m.vert.begin();  vi!=m.vert.end(); vi++) if (!vi->IsD()) {
    SetValency(&*vi,0);
  }
  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) {
     for (int w=0; w<3; w++) 
     if (!fi->IsF(w))
       IncreaseValency( fi->V(w)); 
  }
}

static void UpdateValencyInQuality(MeshType& m){
  for (VertexIterator vi = m.vert.begin();  vi!=m.vert.end(); vi++) if (!vi->IsD()) {
     vi->Q() = 0; 
  }

  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++) if (!fi->IsD()) {
     for (int w=0; w<3; w++) 
         fi->V(w)->Q() += (fi->IsF(w)||fi->IsF((w+2)%3) )? 0.5f:1; 
  }
}

static bool HasConsistentValencyFlag(MeshType &m) {
  UpdateValencyInFlags(m);
  UpdateValencyInQuality(m);
  for (VertexIterator vi = m.vert.begin();  vi!=m.vert.end(); vi++) if (!vi->IsD()) {
    if (GetValency(&*vi)!=vi->Q()) return false;
  }
  return true;
}

private:
  
// helper function:
// returns quality of a quad formed by points a,b,c,d
// quality is computed as "how squared angles are"
static ScalarType quadQuality(const CoordType &a, const CoordType &b, const CoordType &c, const CoordType &d){
  ScalarType score = 0;
  score += 1 - math::Abs( Cos( a,b,c) );
  score += 1 - math::Abs( Cos( b,c,d) );
  score += 1 - math::Abs( Cos( c,d,a) );
  score += 1 - math::Abs( Cos( d,a,b) );
  return score / 4;
}

// helper function:
// returns quality of a given (potential) quad
static ScalarType quadQuality(FaceType *f, int edge){
  
  CoordType 
    a = f->V0(edge)->P(),
    b = f->FFp(edge)->V2( f->FFi(edge) )->P(),
    c = f->V1(edge)->P(),
    d = f->V2(edge)->P();
  
  return quadQuality(a,b,c,d);

}

/**
helper function:
given a quad edge, retruns:
   0 if that edge should not be rotated
  +1 if it should be rotated clockwise (+1)
  -1 if it should be rotated counterclockwise (-1)
Currently an edge is rotated iff it is shortened by that rotations
(shortcut criterion)
*/
static int TestEdgeRotation(const FaceType &f, int w0)
{
  const FaceType *fa = &f;
  assert(! fa->IsF(w0) );
  ScalarType q0,q1,q2;
  CoordType v0,v1,v2,v3,v4,v5;
  int w1 = (w0+1)%3;
  int w2 = (w0+2)%3;
  
  v0 = fa->P(w0);
  v3 = fa->P(w1);
  
  if (fa->IsF(w2) ) {
    v1 = fa->cFFp(w2)->V2( fa->cFFi(w2) )->P();
    v2 = fa->P(w2);
  } else {
    v1 = fa->P(w2);
    v2 = fa->cFFp(w1)->V2( fa->cFFi(w1) )->P();
  }
  
  const FaceType *fb = fa->cFFp(w0);
  w0 = fa->cFFi(w0);
  
  w1 = (w0+1)%3;
  w2 = (w0+2)%3;
  if (fb->IsF(w2) ) {
    v4 = fb->cFFp(w2)->V2( fb->cFFi(w2) )->P();
    v5 = fb->P(w2);
  } else {
    v4 = fb->P(w2);
    v5 = fb->cFFp(w1)->V2( fb->cFFi(w1) )->P();
  }
  
  /*
  //  max overall quality criterion:
  q0 = quadQuality(v0,v1,v2,v3) +  quadQuality(v3,v4,v5,v0); // keep as is?
  q1 = quadQuality(v1,v2,v3,v4) +  quadQuality(v4,v5,v0,v1); // rotate CW?
  q2 = quadQuality(v5,v0,v1,v2) +  quadQuality(v2,v3,v4,v5); // rotate CCW?
  
  if (q0>=q1 && q0>=q2) return 0;
  if (q1>=q2) return 1;*/
  
  // min distance (shortcut criterion)
  q0 = (v0 - v3).SquaredNorm();
  q1 = (v1 - v4).SquaredNorm();
  q2 = (v5 - v2).SquaredNorm();
  if (q0<=q1 && q0<=q2) return 0;
  if (q1<=q2) return 1;
  return -1;
}
private:
  
// helper function:
// cos of angle abc. This should probably go elsewhere
static ScalarType Cos(const CoordType &a, const CoordType &b, const CoordType &c )
{
  CoordType 
    e0 = b - a,
    e1 = b - c;
  ScalarType d =  (e0.Norm()*e1.Norm());
  if (d==0) return 0.0;
  return (e0*e1)/d;
}








};
}} // end namespace vcg::tri
