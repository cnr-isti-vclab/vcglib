/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
#ifndef _AUTOALIGN_4PCS_H_
#define _AUTOALIGN_4PCS_H_

/**
implementation of the 4PCS method from the paper:
"4-Points Congruent Sets for Robust Pairwise Surface Registration"
D.Aiger, N.Mitra D.Cohen-Or, SIGGRAPH 2008
ps: the name of the variables are out of vcg standard but like the one
used in the paper pseudocode.
*/

#include <vcg/complex/complex.h>
#include <vcg/space/point_matching.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/math/random_generator.h>
namespace vcg{
namespace tri{

template <class MeshType>
class FourPCS {
public:
    /* mesh only for using spatial indexing functions (to remove) */
  class PVertex;    // dummy prototype never used
  class PFace;

  class PUsedTypes: public vcg::UsedTypes < vcg::Use<PVertex>::template AsVertexType,
                                            vcg::Use<PFace  >::template AsFaceType >{};

  class PVertex : public vcg::Vertex< PUsedTypes,vcg::vertex::BitFlags,vcg::vertex::Coord3f ,vcg::vertex::Mark>{};
  class PFace   : public vcg::Face<   PUsedTypes> {};
  class PMesh   : public vcg::tri::TriMesh< std::vector<PVertex>, std::vector<PFace> > {};

  typedef typename MeshType::ScalarType ScalarType;
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::VertexPointer VertexPointer;
  typedef typename MeshType::VertexType VertexType;
  typedef vcg::Point4< vcg::Point3<ScalarType> > FourPoints;
  typedef vcg::GridStaticPtr<typename PMesh::VertexType, ScalarType > GridType;

  /* class for Parameters */
  struct Param
  {
    ScalarType delta;   // Approximation Level
    int feetSize;       // how many points in the neighborhood of each of the 4 points
    ScalarType overlap; // overlap estimation as a percentage
    int scoreFeet;      // how many of the feetsize points must match (max feetsize*4) to try an early interrupt
    int n_samples_on_Q; // number of samples on P
    int seed;
    ScalarType cosAngle; // max admittable angle that can be admitted between matching points in alignments (expressed as cos(ang) )

    void Default(){
      delta = 0.5;
      feetSize = 25;
      overlap = 0.5;
      scoreFeet = 50;
      n_samples_on_Q=500;
      seed =0;
      cosAngle = 0; // normals must differ more than 90 degree to be considered bad.
    }
  };

  class Couple
  {
  public:
    VertexPointer p0,p1;
    Couple(VertexPointer i, VertexPointer j, float d) : p0(i),p1(j),dist(d){}
    float dist;
    const bool operator < (const   Couple & o) const {return dist < o.dist;}
    VertexPointer operator[](const int &i){return (i==0)? this->p0 : this->p1;}
  };

  struct Candidate
  {
    Candidate(){}
    Candidate(FourPoints _p,vcg::Matrix44<ScalarType>_T):p(_p),T(_T){}
    FourPoints  p;
    vcg::Matrix44<ScalarType> T;
    int score;
    int base; // debug: for which base
    inline bool operator <(const Candidate & o) const {return score > o.score;}
  };


  // class for the point  'ei'
  struct EPoint{
    EPoint(vcg::Point3<ScalarType> _p, int _i):pos(_p),pi(_i){}
    vcg::Point3<ScalarType> pos;
    int pi;        //index to R[1|2]
    void GetBBox(vcg::Box3<ScalarType> & b){b.Add(pos);}
  };


  Param par;    /// parameters

  MeshType  *P;    // mesh from which the coplanar base is selected
  MeshType  *Q;    // mesh where to find the correspondences

  std::vector<VertexPointer> subsetQ;  // subset of the vertices in Q
  std::vector<VertexPointer> subsetP; // random selection on P

  ScalarType side;               // side

  PMesh     Invr;                // invariants
  math::MarsenneTwisterRNG rnd;
  std::vector< Candidate > U;    // the
  int iwinner;                   // winner == U[iwinner]
  FourPoints B;                  // coplanar base
  std::vector<FourPoints> bases; // used bases
  std::vector<VertexType*> ExtB[4]; // selection of vertices "close" to the four point
  ScalarType radius;

  std::vector<Couple > R1;
  ScalarType r1,r2;

    GridType *ugrid; // griglia
    vcg::GridStaticPtr<typename MeshType::VertexType, ScalarType > ugridQ;
    vcg::GridStaticPtr<typename MeshType::VertexType, ScalarType > ugridP;

 // the two main functions to be used
  void Init(MeshType &Mov, MeshType &Fix);
  bool Align( int   L, vcg::Matrix44f & result, vcg::CallBackPos * cb = NULL );

// ---- auxiliary functions
    bool SelectCoplanarBase(); // on P
    bool FindCongruent();      // of base B, on Q, with approximation delta
    void ComputeR1();
    bool IsTransfCongruent(FourPoints fp,vcg::Matrix44<ScalarType> & mat, float &  trerr);
    int EvaluateSample(Candidate & fp, CoordType & tp, CoordType & np);
    void EvaluateAlignment(Candidate & fp);
    void TestAlignment(Candidate & fp);

    /* returns the closest point between to segments x1-x2 and x3-x4.  */
    void IntersectionLineLine(const CoordType & x1,const CoordType & x2,const CoordType & x3,const CoordType & x4, CoordType&x)
    {
      CoordType a = x2-x1, b = x4-x3, c = x3-x1;
      x = x1 + a * ((c^b).dot(a^b)) / (a^b).SquaredNorm();
    }


    /* debug tools */
public:
    std::vector<vcg::Matrix44f> allTr;// tutte le trasformazioni provate
    FILE * db;
    char namemesh1[255],namemesh2[255];
    int n_base;
    void InitDebug(const char * name1, const char * name2){
        db = fopen("debugPCS.txt","w");
        sprintf(&namemesh1[0],"%s",name1);
        sprintf(&namemesh2[0],"%s",name2);
        n_base = 0;
    }

    void FinishDebug(){
        fclose(db);
    }
};

template <class MeshType>
void FourPCS<MeshType>:: Init(MeshType &_movP,MeshType &_fixQ)
{
        P = &_movP;
        Q = &_fixQ;
        tri::UpdateBounding<MeshType>::Box(*P);
        if(par.seed==0) rnd.initialize(time(0));
        else rnd.initialize(par.seed);
        ugridQ.Set(Q->vert.begin(),Q->vert.end());
        ugridP.Set(P->vert.begin(),P->vert.end());

        float radius=0;
        tri::PoissonPruning(*Q,subsetQ,radius,par.n_samples_on_Q,par.seed);
        tri::PoissonPruning(*P,subsetP,radius,par.n_samples_on_Q,par.seed);
        float ratio = std::min<int>(Q->vert.size(),par.n_samples_on_Q) / (float) Q->vert.size();

        // estimate neigh distance
        float avD = 0.0;
        for(int i = 0 ; i < 100; ++i){
            int ri = rnd.generate(Q->vert.size());
            std::vector< CoordType > samples;
            std::vector<ScalarType > dists;
            std::vector<VertexType* > ress;
            vcg::tri::GetKClosestVertex<
                    MeshType,
                    vcg::GridStaticPtr<typename MeshType::VertexType, ScalarType>,
                    std::vector<VertexType*>,
                    std::vector<ScalarType>,
                    std::vector< CoordType > >(*Q,ugridQ,2,Q->vert[ri].cP(),Q->bbox.Diag(), ress,dists, samples);
            assert(ress.size() == 2);
            avD+=dists[1];
        }
        avD    /=100;                        // average vertex-vertex distance
        avD /= sqrt(ratio);

        par.delta = avD * par.delta;
        side = P->bbox.Dim()[P->bbox.MaxDim()]*par.overlap; //rough implementation
    }

template <class MeshType>
bool FourPCS<MeshType>::SelectCoplanarBase()
{
  // choose the inter point distance
  ScalarType dtol = side*0.1; //rough implementation

  //choose the first two points

  // first point random
  B[0] = P->vert[ rnd.generate(P->vert.size())].P();

  // second a point at distance d+-dtol
  int i;
  for(i = 0; i < P->vert.size(); ++i){
    int id = rnd.generate(P->vert.size()-1);
    ScalarType dd = (P->vert[id].P() - B[0]).Norm();
    if(  ( dd < side + dtol) && (dd > side - dtol)){
      B[1] = P->vert[id].P();
      break;
    }
  }
  if(i ==  P->vert.size()) return false;

  // third point at distance side*0.8 from middle way between B[0] and B[1]
  const vcg::Point3f middle = (B[0]+B[1])/2.0;
  for(i = 0; i < P->vert.size(); ++i){
    int id = rnd.generate(P->vert.size()-1);
    ScalarType dd = (P->vert[id].P() - middle).Norm();
    if(  ( dd < side*0.8) ){
      B[2] = P->vert[id].P();
      break;
    }
  }
  if(i ==  P->vert.size()) return false;

        //fourth point
        float cpr = rnd.generate01();
        vcg::Point3f crossP = B[0] *(1-cpr)+B[1]*cpr;
        CoordType B4 = B[2]+(crossP-B[2]).Normalize()*side;
        CoordType n = ((B[0]-B[1]).normalized() ^ (B[2]-B[1]).normalized()).normalized();
        ScalarType radius = dtol;

        std::vector<typename MeshType::VertexType*> closests;
        std::vector<ScalarType> distances;
        std::vector<CoordType> points;

         vcg::tri::GetInSphereVertex<
                    MeshType,
                    vcg::GridStaticPtr<typename MeshType::VertexType, ScalarType >,
                    std::vector<typename MeshType::VertexType*>,
                    std::vector<ScalarType>,
                    std::vector<CoordType>
                >(*P,ugridP,B4,radius,closests,distances,points);

        if(closests.empty())
            return false;
        int bestInd = -1;   ScalarType bestv=std::numeric_limits<float>::max();
        for(i = 0; i <closests.size(); ++i){
         ScalarType dist_from_plane = fabs((closests[i]->P() - B[1]).normalized().dot(n));
            if( dist_from_plane < bestv){
                bestv = dist_from_plane;
                bestInd = i;
            }
        }
        if(bestv >dtol)
            return false;
        B[3] =  closests[bestInd]->P();

//printf("B[3] %d\n", (typename MeshType::VertexType*)closests[best] - &(*P->vert.begin()));

        // compute r1 and r2
        CoordType x;
//        std::swap(B[1],B[2]);
        IntersectionLineLine(B[0],B[1],B[2],B[3],x);

        r1 = (x - B[0]).dot(B[1]-B[0]) / (B[1]-B[0]).SquaredNorm();
        r2 = (x - B[2]).dot(B[3]-B[2]) / (B[3]-B[2]).SquaredNorm();

        if( ((B[0]+(B[1]-B[0])*r1)-(B[2]+(B[3]-B[2])*r2)).Norm() > par.delta )
            return false;

        radius  =side*0.5;
        std::vector< CoordType > samples;
        std::vector<ScalarType > dists;

        for(int i  = 0 ; i< 4; ++i){
            vcg::tri::GetKClosestVertex<
                MeshType,
                vcg::GridStaticPtr<typename MeshType::VertexType, ScalarType >,
                std::vector<VertexType*>,
                std::vector<ScalarType>,
                std::vector< CoordType > >(*P,ugridP, par.feetSize ,B[i],radius, ExtB[i],dists, samples);
        }

    //for(int i  = 0 ; i< 4; ++i)
 //        printf("%d ",ExtB[i].size());
    //    printf("\n");
return true;

}


template <class MeshType>
bool FourPCS<MeshType>::IsTransfCongruent(FourPoints fp, vcg::Matrix44<ScalarType> & mat, float &  trerr)
{
  std::vector<vcg::Point3<ScalarType> > fix(4);
  std::vector<vcg::Point3<ScalarType> > mov(4);
  for(int i = 0 ; i < 4; ++i) {
    mov[i]=B[i];
    fix[i]=fp[i];
  }

  if(fabs( Distance(fix[0],fix[1]) - Distance(mov[0],mov[1]) ) > par.delta) return false;
  if(fabs( Distance(fix[0],fix[2]) - Distance(mov[0],mov[2]) ) > par.delta) return false;
  if(fabs( Distance(fix[0],fix[3]) - Distance(mov[0],mov[3]) ) > par.delta) return false;
  if(fabs( Distance(fix[1],fix[2]) - Distance(mov[1],mov[2]) ) > par.delta) return false;
  if(fabs( Distance(fix[1],fix[3]) - Distance(mov[1],mov[3]) ) > par.delta) return false;
  if(fabs( Distance(fix[2],fix[3]) - Distance(mov[2],mov[3]) ) > par.delta) return false;

  /*
  vcg::Point3<ScalarType> n,p;
  n = (( B[1]-B[0]).normalized() ^ ( B[2]- B[0]).normalized())*( B[1]- B[0]).Norm();
  p =  B[0] + n;
  mov.push_back(p);
  n = (( fp[1]-fp[0]).normalized() ^ (fp[2]- fp[0]).normalized())*( fp[1]- fp[0]).Norm();
  p =  fp[0] + n;
  fix.push_back(p);
*/

  vcg::ComputeRigidMatchMatrix(fix,mov,mat);

  ScalarType err = 0.0;
  for(int i = 0; i < 4; ++i) err+= (mat * mov[i] - fix[i]).SquaredNorm();

  trerr = vcg::math::Sqrt(err);
  return  trerr  < par.delta;
}

template <class MeshType>
void
FourPCS<MeshType>::ComputeR1()
{
    R1.clear();
    for(int vi = 0; vi  < subsetQ.size(); ++vi)
      for(int vj = vi; vj < subsetQ.size(); ++vj){
//          ScalarType d = ((Q->vert[subsetQ[vi]]).P()-(Q->vert[subsetQ[vj]]).P()).Norm();
            ScalarType d = (subsetQ[vi]->P()- subsetQ[vj]->P()).Norm();
              if( (d < side+par.delta))
            {
                R1.push_back(Couple(subsetQ[vi],subsetQ[vj],d ));
                R1.push_back(Couple(subsetQ[vj],subsetQ[vi],d));
            }
    }

    std::sort(R1.begin(),R1.end());
}

template <class MeshType>
bool FourPCS<MeshType>::FindCongruent() { // of base B, on Q, with approximation delta
    bool done = false;
    std::vector<EPoint> R2inv;
    int n_closests = 0, n_congr = 0;
    int ac =0 ,acf = 0,tr = 0,trf =0;
    ScalarType d1,d2;
    d1 = (B[1]-B[0]).Norm();
    d2 = (B[3]-B[2]).Norm();

    typename PMesh::VertexIterator vii;
    typename std::vector<Couple>::iterator bR1,eR1,bR2,eR2,ite;
    bR1 = std::lower_bound<typename std::vector<Couple>::iterator,Couple>(R1.begin(),R1.end(),Couple(0,0,d1-par.delta));
    eR1 = std::lower_bound<typename std::vector<Couple>::iterator,Couple>(R1.begin(),R1.end(),Couple(0,0,d1+par.delta));
    bR2 = std::lower_bound<typename std::vector<Couple>::iterator,Couple>(R1.begin(),R1.end(),Couple(0,0,d2-par.delta));
    eR2 = std::lower_bound<typename std::vector<Couple>::iterator,Couple>(R1.begin(),R1.end(),Couple(0,0,d2+par.delta));

    // in  [bR1,eR1) there are all the pairs ad a distance d1 +- par.delta
    // in  [bR1,eR1) there are all the pairs ad a distance d2 +- par.delta

    if(bR1 == R1.end()) return false;// if there are no such pairs return
    if(bR2 == R1.end()) return false; // if there are no such pairs return

    // put [bR1,eR1) in a mesh to have the search operator for free (lazy me)
    Invr.Clear();
    int i = &(*bR1)-&(*R1.begin());
    for(ite = bR1; ite != eR1;++ite){
        vii = vcg::tri::Allocator<PMesh>::AddVertices(Invr,1);
//      (*vii).P() = Q->vert[R1[i][0]].P() + (Q->vert[R1[i][1]].P()-Q->vert[R1[i][0]].P()) * r1;
        (*vii).P() =         R1[i].p0->P() + (        R1[i].p1->P() -       R1[i].p0->P()) * r1;
        ++i;
    }
    if(Invr.vert.empty() ) return false;

    // index remaps a vertex of Invr to its corresponding point in R1
    typename PMesh::template PerVertexAttributeHandle<int> id = vcg::tri::Allocator<PMesh>::template AddPerVertexAttribute<int>(Invr,std::string("index"));
    i = &(*bR1)-&(*R1.begin());
    for(vii = Invr.vert.begin(); vii != Invr.vert.end();++vii,++i)  id[vii] = i;

    vcg::tri::UpdateBounding<PMesh>::Box(Invr);
    //    printf("Invr size %d\n",Invr.vn);

    ugrid = new GridType();
    ugrid->Set(Invr.vert.begin(),Invr.vert.end());

    i = &(*bR2)-&(*R1.begin());
    // R2inv contains all the points generated by the couples in R2 (with the reference to remap into R2)
    for(ite = bR2; ite != eR2;++ite){
//        R2inv.push_back( EPoint( Q->vert[R1[i][0]].P() + (Q->vert[R1[i][1]].P()-Q->vert[R1[i][0]].P()) * r2,i));
        R2inv.push_back( EPoint( R1[i].p0->P() + (R1[i].p1->P() - R1[i].p0->P()) * r2,i));
        ++i;
    }

    n_closests = 0; n_congr = 0; ac =0 ; acf = 0; tr = 0; trf = 0;
    printf("R2Inv.size  = %d \n",R2inv.size());
    for(uint i = 0 ; i < R2inv.size() ; ++i){

        std::vector<typename PMesh::VertexType*> closests;

        // for each point in R2inv get all the points in R1 closer than par.delta
        vcg::Matrix44<ScalarType> mat;
        vcg::Box3f bb;
        bb.Add(R2inv[i].pos+vcg::Point3f(par.delta,par.delta, par.delta));
        bb.Add(R2inv[i].pos-vcg::Point3f(par.delta,par.delta, par.delta));

        vcg::tri::GetInBoxVertex<PMesh,GridType,std::vector<typename PMesh::VertexType*> >
             (Invr,*ugrid,bb,closests);

        if(closests.size() > 5)
            closests.resize(5);

         n_closests+=closests.size();
         for(uint ip = 0; ip < closests.size(); ++ip){
                FourPoints p;
                p[0] = R1[id[closests[ip]]][0]->P();
                p[1] = R1[id[closests[ip]]][1]->P();
                p[2] = R1[ R2inv[i].pi][0]->P();
                p[3] = R1[ R2inv[i].pi][1]->P();

                float trerr;
              n_base++;
                    if(!IsTransfCongruent(p,mat,trerr)) {
                        trf++;
                        //char name[255];
                        //sprintf(name,"faileTR_%d_%f.aln",n_base,trerr);
                        //fprintf(db,"TransCongruent %s\n", name);
                        //SaveALN(name, mat);
                    }
                    else{
                        tr++;
                        n_congr++;
                        Candidate c(p,mat);
                        EvaluateAlignment(c);

                        if( c.score > par.scoreFeet)
                            U.push_back(c);

/*
                        EvaluateAlignment(U.back());
                        U.back().base = bases.size()-1;

                        if( U.back().score > par.scoreFeet){
                            TestAlignment(U.back());
                            if(U.back().score > par.scoreAln)
                                {
                                    done = true; break;
                                }
                            }
*/
                        //char name[255];
                        //sprintf(name,"passed_score_%5d_%d.aln",U.back().score,n_base);
                        //fprintf(db,"OK TransCongruent %s, score: %d \n", name,U.back().score);
                        //SaveALN(name, mat);
                    }
                }
     }

     delete ugrid;
     vcg::tri::Allocator<PMesh>::DeletePerVertexAttribute(Invr,id);
     printf("n_closests %5d = (An %5d ) + ( Tr %5d ) + (OK) %5d\n",n_closests,acf,trf,n_congr);

     return done;
//     printf("done n_closests %d congr %d in %f s\n ",n_closests,n_congr,(clock()-start)/(float)CLOCKS_PER_SEC);
//     printf("angle:%d %d, trasf %d %d\n",ac,acf,tr,trf);
}



template <class MeshType>
int FourPCS<MeshType>::EvaluateSample(Candidate & fp, CoordType & tp, CoordType & np)
{
  VertexType*   v=0;
  ScalarType   dist ;
  radius = par.delta;
  tp = fp.T * tp;

  vcg::Point4<ScalarType> np4;
  np4 = fp.T * vcg::Point4<ScalarType>(np[0],np[1],np[2],0.0);
  np[0] = np4[0]; np[1] = np4[1];     np[2] = np4[2];

  if(ugridQ.bbox.IsIn(tp))
   v = vcg::tri::GetClosestVertex(*Q, ugridQ, tp, radius,  dist  );

  if(v!=0)
  {
    if( v->N().dot(np) > par.cosAngle )  return 1;
    else return -1;
  }
  else return 0;
}


template <class MeshType>
void
FourPCS<MeshType>::EvaluateAlignment(Candidate  & fp){
        int n_delta_close = 0;
        for(int i  = 0 ; i< 4; ++i) {
            for(uint j = 0; j < ExtB[i].size();++j){
                CoordType np = ExtB[i][j]->cN();;
                CoordType tp = ExtB[i][j]->P();
                n_delta_close+=EvaluateSample(fp,tp,np);
            }
        }
        fp.score = n_delta_close;
}

template <class MeshType>
void FourPCS<MeshType>::TestAlignment(Candidate  & fp){
        radius = par.delta;
        int n_delta_close = 0;
         for(uint j = 0; j < subsetP.size();++j){
                CoordType np = subsetP[j]->N();
                CoordType tp = subsetP[j]->P();
                n_delta_close+=EvaluateSample(fp,tp,np);
             }
        fp.score =  n_delta_close;
}


template <class MeshType>
bool FourPCS<MeshType>::     Align(  int L, vcg::Matrix44f & result, vcg::CallBackPos * cb )
{
    int bestv = 0;
    bool found;
    int n_tries = 0;
    U.clear();

    if(L==0)
    {
        L = (log(1.0-0.9) / log(1.0-pow((float)par.overlap,3.f)))+1;
        printf("using %d bases\n",L);
    }

    ComputeR1();

    for(int t  = 0; t  < L; ++t )
    {
      do
      {
        n_tries = 0;
        do
        {
          n_tries++;
          found = SelectCoplanarBase();
        }
        while(!found && (n_tries <50));
        if(!found) {
          par.overlap*=0.9;
          side = P->bbox.Dim()[P->bbox.MaxDim()]*par.overlap; //rough implementation
          ComputeR1();
        }
      } while (!found && (par.overlap >0.1));

      if(par.overlap <0.1) {
        printf("FAILED");
        return false;
      }
      bases.push_back(B);
      if(cb) cb(t*100/L,"Trying bases");
      if(FindCongruent())
        break;
    }

    if(U.empty()) return false;

//    std::sort(U.begin(),U.end());
    if(cb) cb(90,"TestAlignment");
    bestv  = -std::numeric_limits<float>::max();
    iwinner = 0;

    for(int i = 0 ; i <  U.size() ;++i)
     {
        TestAlignment(U[i]);
        if(U[i].score > bestv){
            bestv = U[i].score;
            iwinner = i;
            }
    }

    result = U[iwinner].T;
    Invr.Clear();
    return true;
}

} // namespace tri
} // namespace vcg
#endif
