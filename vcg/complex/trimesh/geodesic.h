/*#***************************************************************************
 * Geodesic.h                                                         o o    *
 *                                                                  o     o  *
 * Visual Computing Group                                           _  O  _  *
 * IEI Institute, CNUCE Institute, CNR Pisa                          \/)\/   *
 *                                                                  /\/|     *
 * Copyright(C) 1999 by Paolo Cignoni,                                 |     *
 * All rights reserved.                                                \     *
 *                                                                           *
 * Permission  to use, copy, modify, distribute  and sell this  software and *
 * its documentation for any purpose is hereby granted without fee, provided *
 * that  the above copyright notice appear  in all copies and that both that *
 * copyright   notice  and  this  permission  notice  appear  in  supporting *
 * documentation. the author makes  no representations about the suitability *
 * of this software for any purpose. It is provided  "as is" without express *
 * or implied warranty.                                                      *
 *                                                                           *
 ***************************************************************************#*/
/*#**************************************************************************
  History
	$Log: not supported by cvs2svn $

 *#**************************************************************************/



#include <assert.h>
#include <vcg/container/simple_temporary_data.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/math/base.h>
#include <deque>
#include <vector>
#include <list>
#include <functional>

namespace vcg{
template <class MeshType>
class Geo{

	public:

	typedef typename MeshType::VertexPointer VertexPointer;

	template <class MeshType>
	struct TempData{
		TempData(){}
		TempData(const double & d_){d=d_;visited=false;}
		double d;
		bool visited;
		};

	typedef SimpleTempData<std::vector<typename MeshType::VertexType>, TempData<MeshType> >  TempDataType;
	static TempDataType  & TD(){ static TempDataType td; return td;}
	
	
struct pred: public std::binary_function<VertexPointer,VertexPointer,bool>{
			bool operator()(const VertexPointer& v0, const VertexPointer& v1) const
				{return (Geo<MeshType>::TD()[v0].d > Geo<MeshType>::TD()[v1].d);}
		};

		
static 	 typename MeshType::VertexPointer BuildSP( 
		 MeshType & m,
		 std::vector<typename MeshType::VertexPointer> & _frontier,
		 double & max_distance,
		 bool fartestOnBorder = false
		 )
		 {
		 TD().c = &m.vert;

		 bool isLeaf;
		 std::vector<typename MeshType::VertexPointer> frontier;
		 std::vector<typename MeshType::VertexPointer> :: iterator tmp;
		 frontier.clear();
		 //Requirements
		 assert(m.HasVFTopology());

		 if(m.vn==0) return NULL;

		 MeshType::VertexIterator ii;
		 std::list<typename MeshType::VertexPointer> children;
		 typename MeshType::VertexPointer curr,fartest,pw1;	
		 bool toQueue;

		 TD().Start(TempData<MeshType>(-1.0));

		 std::list<typename MeshType::VertexPointer>::iterator is;
		 std::deque<typename MeshType::VertexPointer> leaves;
		 std::vector <std::pair<typename MeshType::VertexPointer,typename MeshType::ScalarType> > expansion;

		 std::vector <typename MeshType::VertexPointer >::const_iterator ifr;
		 face::VFIterator<typename MeshType::FaceType> x;int k;
		 typename MeshType::VertexPointer pw;

		 for(ifr = _frontier.begin(); ifr != _frontier.end(); ++ifr){
			 TD()[*ifr].visited= true;
			 TD()[*ifr].d = 0.0;	
			 }

 for(ifr = _frontier.begin(); ifr != _frontier.end(); ++ifr)
	 {	
		 // determina la distanza dei vertici della fan
		 for( x.f = (*ifr)->VFp(), x.z = (*ifr)->VFi(); x.f!=0; ++x )
			 for(k=0;k<2;++k)
				 {
				 if(k==0) pw = x.f->V1(x.z);
				 else     pw = x.f->V2(x.z);

				if(TD()[pw].d ==-1){
					TD()[pw].d = Distance(pw->cP(),(*ifr)->cP());
					frontier.push_back(pw);				
					}
				}
	 }
    // initialize Heap
	 make_heap(frontier.begin(),frontier.end(),pred());	
	  double curr_d,d_curr = 0.0;
		max_distance=0.0;
		std::vector<typename MeshType::VertexPointer >:: iterator iv;

 		while(!frontier.empty())
			{ //printf("size: %d\n", frontier.size());
				expansion.clear();
				pop_heap(frontier.begin(),frontier.end(),pred());
				curr = frontier.back();
				frontier.pop_back();
				d_curr = TD()[curr].d;
				TD()[curr].visited = true;

				isLeaf = (!fartestOnBorder || curr->IsB());

			 face::VFIterator<typename MeshType::FaceType> x;int k;

			
				for( x.f = curr->VFp(), x.z = curr->VFi(); x.f!=0; ++x )
					for(k=0;k<2;++k)
						{
						if(k==0) {
							pw = x.f->V1(x.z);
							pw1=x.f->V2(x.z);
							}
						else {
							pw = x.f->V2(x.z);
							pw1=x.f->V1(x.z);
							}
				
				const double & d_pw1  = TD()[pw1].d;

				if((!TD()[pw1].visited ) || d_curr == 0.0)
					{
					if(TD()[pw].d == -1){
								curr_d = TD()[curr].d + (pw->P()-curr->P()).Norm();
								expansion.push_back(std::pair<typename MeshType::VertexPointer,typename MeshType::ScalarType>(pw,curr_d));
					 }
					 continue;
					}
	
		 assert( TD()[pw1].d != -1);
		 assert( (curr!=pw) && (pw!=pw1) && (pw1 != curr));				
		 assert(d_pw1!=-1.0);

		 //************** calcolo della distanza di pw in base alle distanze note di pw1 e curr
		 //************** sapendo che (curr,pw,pw1) e'una faccia della mesh
		 //************** (vedi figura in file distance.gif)
		 Point3<MeshType::ScalarType> w_c = pw->cP()- curr->cP();
		 Point3<MeshType::ScalarType> w_w1 = pw->cP()- pw1->cP();
		 Point3<MeshType::ScalarType> w1_c = pw1->cP()- curr->cP();

		 double ew_c  = (w_c).Norm();
		 double ew_w1 = (w_w1).Norm();
		 double ec_w1 = (w1_c).Norm();
		 double	alpha,alpha_, beta,beta_,theta_c,theta,h,delta,s,a,b;

		 alpha = acos((w_c*w1_c)/(ew_c*ec_w1));
		 s = (d_curr + d_pw1+ec_w1)/2;
		 a = s/ec_w1;
		 b = a*s;
		 alpha_ = 2*acos ( math::Min(1.0,sqrt(  (b- a* d_pw1)/d_curr)));

		 if ( alpha+alpha_ > M_PI){
			 curr_d = d_curr + ew_c;		
			 }else
			 {
				 beta_ = 2*acos ( math::Min(1.0,sqrt(  (b- a* d_curr)/d_pw1)));
				 beta  = acos((w_w1)*(-w1_c)/(ew_w1*ec_w1));

				 if ( beta+beta_ > M_PI)
					 curr_d = d_pw1  + ew_w1;
				 else 
					 {
					 theta	= M_PI-alpha-alpha_;
					 delta	= cos(theta)* ew_c;
					 h		= sin(theta)* ew_c;
					 curr_d = sqrt( pow(h,2)+ pow(d_curr + delta,2));
					 }
			 }
		 //**************************************************************************************
		 toQueue = (TD()[(pw)].d==-1);

		 if(toQueue){// se non e'gia' in coda ce lo mette
			 expansion.push_back(std::pair<typename MeshType::VertexPointer,typename MeshType::ScalarType>(pw,curr_d));
			 }else
			 {
				 if(  TD()[(pw)].d > curr_d )
						TD()[(pw)].d = curr_d;
			 }
				
		 if(isLeaf){
			 if(d_curr > max_distance){
				 max_distance = d_curr;
				 fartest = curr;
				 }
			 }


				}
		 std::vector <std::pair<typename MeshType::VertexPointer,typename MeshType::ScalarType> > ::iterator i;
		 for(i = expansion.begin(); i!= expansion.end(); ++i)
				{
					TD()[(*i).first].d = (*i).second;
					frontier.push_back((*i).first);
					push_heap(frontier.begin(),frontier.end(),pred());
						} // end for
}// end while

// scrivi le distanze sul campo qualita' (nn: farlo parametrico)
 MeshType::VertexIterator vi;
 for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
	 (*vi).Q() = TD()[&(*vi)].d; 
 
  
  TD().Stop();	

 return fartest;

 }
	

public:
static void FartestPoint( MeshType & m,
									std::vector<typename MeshType::VertexPointer> & fro,//insieme di vertici da cui trovare le distanze
									typename MeshType::VertexPointer & fartest,		//punto piu'lontano
									double & distance){					//distaza geodesica
									fartest = BuildSP(m,fro,distance,false); 
					  }
static void FartestBPoint(
									 MeshType & m,
									 std::vector<typename MeshType::VertexPointer> & fro, //insieme di vertici da cui trovare le distanze
									 typename MeshType::VertexPointer & fartest,	    //punto piu'lontano
									 double & distance){
									 fartest =  BuildSP(m,fro,distance,true); 
					  }

static void DistanceFromBorder(	MeshType & m,
									  typename MeshType::VertexPointer & v0,	//ritorna il vertice piu'lontano da ogni punto sul bordo
										typename MeshType::VertexPointer & v1, // ritorna il vertice di bordo piu'vicino a v0
										double & distance			// distanza geodesica tra v0 e v1
					)
					{
					std::vector<typename MeshType::VertexPointer> fro;
					MeshType::VertexIterator vi;
					MeshType::VertexPointer fartest;
					for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
						if( (*vi).IsB())
							fro.push_back(&(*vi));
					FartestPoint(m,fro,fartest,distance);
					}

static void MostInternal(	MeshType & m,
									  typename MeshType::VertexPointer & v0,	//ritorna il vertice piu'lontano da ogni punto sul bordo
										typename MeshType::VertexPointer & v1, // ritorna il vertice di bordo piu'vicino a v0
										double & distance			// distanza geodesica tra v0 e v1
					)
					{
					std::vector<typename MeshType::VertexPointer> fro;
					MeshType::VertexIterator vi;
					MeshType::VertexPointer fartest;
					for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
						if( (*vi).IsB())
							fro.push_back(&(*vi));
					FartestPoint(m,fro,fartest,distance);
					fro.clear();
					fro.push_back(fartest);
					FartestBPoint(m,fro,fartest,distance);
					}

};
};// end namespace