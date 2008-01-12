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
/*#**************************************************************************
  History
	$Log: not supported by cvs2svn $
	Revision 1.6  2008/01/12 19:07:05  ganovelli
	Recompiled from previous out of date version. Still to revise but working
	
	Revision 1.5  2005/12/13 17:17:19  ganovelli
	first importing from old version. NOT optimized! It works with VertexFace Adjacency even over non manifolds
	

 *#**************************************************************************/



#include <assert.h>
#include <vcg/container/simple_temporary_data.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/math/base.h>
#include <deque>
#include <vector>
#include <list>
#include <functional>

/*
class for computing approximated geodesic distances on a mesh.

basic example: farthest vertex from a specified one
	MyMesh m;
	MyMesh::VertexPointer seed,far;
	MyMesh::ScalarType dist;

  vcg::Geo<MyMesh> g;
	g.FarthestVertex(m,seed,far,d);

*/
namespace vcg{
template <class MeshType>
class Geo{

	public:

	typedef typename MeshType::VertexPointer VertexPointer;
	typedef typename MeshType::ScalarType  ScalarType;

	/* Auxiliary class for keeping the heap of vertices to visit and their estimated distance
	*/
	struct VertDist{
		VertDist(){}
		VertDist(VertexPointer _v, ScalarType _d):v(_v),d(_d){}
		VertexPointer v;
		ScalarType d;
	};


	/* Temporary data to associate to all the vertices: estimated distance and boolean flag
	*/
	template <class MeshType>
	struct TempData{
		TempData(){}
		TempData(const ScalarType & d_){d=d_;visited=false;}
		ScalarType d;
		bool visited;
		};

	typedef SimpleTempData<std::vector<typename MeshType::VertexType>, TempData<MeshType> >  TempDataType;
	TempDataType  * TD;
	
	
	struct pred: public std::binary_function<VertDist,VertDist,bool>{
		pred(){};
			bool operator()(const VertDist& v0, const VertDist& v1) const
				{return (v0.d > v1.d);}
		};

	/*
	  starting from the seeds, it assign a distance value to each vertex. The distance of a vertex is its 
	  approximated geodesic distance to the closest seeds.
	  This is function is not meant to be called (although is not prevented). Instead, it is invoked by 
	  wrapping function.
	*/
 	 typename MeshType::VertexPointer Visit( 
		 MeshType & m,
		 std::vector<VertDist> & _frontier,
		 ScalarType & max_distance,
		 bool fartestOnBorder = false
		 )
{
	bool isLeaf,toQueue;
	std::vector<VertDist> frontier;
	MeshType::VertexIterator ii;
	std::list<typename MeshType::VertexPointer> children;
	typename MeshType::VertexPointer curr,farthest,pw1;	
	std::list<typename MeshType::VertexPointer>::iterator is;
	std::deque<typename MeshType::VertexPointer> leaves;
	std::vector <std::pair<typename MeshType::VertexPointer,typename MeshType::ScalarType> > expansion;
	std::vector <VertDist >::iterator ifr;
	face::VFIterator<typename MeshType::FaceType> x;int k;
	typename MeshType::VertexPointer pw;

	//Requirements
	assert(m.HasVFTopology());
	assert(!_frontier.empty());

	TD = new TempDataType(m.vert);
	TD->Start(TempData<MeshType>(-1.0));

	for(ifr = _frontier.begin(); ifr != _frontier.end(); ++ifr){
		(*TD)[(*ifr).v].visited= true;
		(*TD)[(*ifr).v].d = 0.0;	
		(*ifr).d = 0.0;
		}

	for(ifr = _frontier.begin(); ifr != _frontier.end(); ++ifr)
		{	
			// determina la distanza dei vertici della fan
			for( x.f = (*ifr).v->VFp(), x.z = (*ifr).v->VFi(); x.f!=0; ++x )
				for(k=0;k<2;++k)
					{
					if(k==0) pw = x.f->V1(x.z);
					else     pw = x.f->V2(x.z);

					if((*TD)[pw].d ==-1){
						(*TD)[pw].d = Distance(pw->cP(),(*ifr).v->cP());
						frontier.push_back(VertDist(pw,(*TD)[pw].d));				
						}
					}
		}
		// initialize Heap
		make_heap(frontier.begin(),frontier.end(),pred());	
		ScalarType curr_d,d_curr = 0.0;
		max_distance=0.0;
		std::vector<VertDist >:: iterator iv;

 		while(!frontier.empty())
			{ //printf("size: %d\n", frontier.size());
				expansion.clear();
				pop_heap(frontier.begin(),frontier.end(),pred());
				curr = (frontier.back()).v;
				frontier.pop_back();
				d_curr =  (*TD)[curr].d;
				(*TD)[curr].visited = true;


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
				
				const ScalarType & d_pw1  =  (*TD)[pw1].d;

				if((! (*TD)[pw1].visited ) || d_curr == 0.0)
					{
					if( (*TD)[pw].d == -1){
								curr_d =  (*TD)[curr].d + (pw->P()-curr->P()).Norm();
								expansion.push_back(std::pair<typename MeshType::VertexPointer,typename MeshType::ScalarType>(pw,curr_d));
					}
					continue;
					}

		assert(  (*TD)[pw1].d != -1);
		assert( (curr!=pw) && (pw!=pw1) && (pw1 != curr));				
		assert(d_pw1!=-1.0);

		//************** calcolo della distanza di pw in base alle distanze note di pw1 e curr
		//************** sapendo che (curr,pw,pw1) e'una faccia della mesh
		//************** (vedi figura in file distance.gif)
		Point3<MeshType::ScalarType> w_c = pw->cP()- curr->cP();
		Point3<MeshType::ScalarType> w_w1 = pw->cP()- pw1->cP();
		Point3<MeshType::ScalarType> w1_c = pw1->cP()- curr->cP();

		ScalarType ew_c  = (w_c).Norm();
		ScalarType ew_w1 = (w_w1).Norm();
		ScalarType ec_w1 = (w1_c).Norm();
		ScalarType	alpha,alpha_, beta,beta_,theta,h,delta,s,a,b;

		alpha = acos((w_c*w1_c)/(ew_c*ec_w1));
		s = (d_curr + d_pw1+ec_w1)/2;
		a = s/ec_w1;
		b = a*s;
		alpha_ = 2*acos ( math::Min<ScalarType>(1.0,sqrt(  (b- a* d_pw1)/d_curr)));

		if ( alpha+alpha_ > M_PI){
			curr_d = d_curr + ew_c;		
			}else
			{
				beta_ = 2*acos ( math::Min<ScalarType>(1.0,sqrt(  (b- a* d_curr)/d_pw1)));
				beta  = acos((w_w1)*(-w1_c)/(ew_w1*ec_w1));

				if ( beta+beta_ > M_PI)
					curr_d = d_pw1  + ew_w1;
				else 
					{
					theta	= ScalarType(M_PI)-alpha-alpha_;
					delta	= cos(theta)* ew_c;
					h		= sin(theta)* ew_c;
					curr_d = sqrt( pow(h,2)+ pow(d_curr + delta,2));
					}
			}
		//**************************************************************************************
		toQueue = ( (*TD)[(pw)].d==-1);

		if(toQueue){// se non e'gia' in coda ce lo mette
			expansion.push_back(std::pair<typename MeshType::VertexPointer,typename MeshType::ScalarType>(pw,curr_d));
			}else
			{
				if(   (*TD)[(pw)].d > curr_d )
						(*TD)[(pw)].d = curr_d;
			}
				
		if(isLeaf){
			if(d_curr > max_distance){
				max_distance = d_curr;
				farthest = curr;
				}
			}


				}
		std::vector <std::pair<typename MeshType::VertexPointer,typename MeshType::ScalarType> > ::iterator i;
		for(i = expansion.begin(); i!= expansion.end(); ++i)
				{
					(*TD)[(*i).first].d = (*i).second;
					frontier.push_back(VertDist((*i).first,(*TD)[(*i).first].d));
					push_heap(frontier.begin(),frontier.end(),pred());
						} // end for
	}// end while

	// scrivi le distanze sul campo qualita' (nn: farlo parametrico)
	MeshType::VertexIterator vi;
	for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
		(*vi).Q() =  (*TD)[&(*vi)].d; 


	(*TD).Stop();

	delete TD;

	return farthest;

 }
	

public:
	/*
	Given a mesh and  a vector of pointers to vertices (sources), assigns the approximated geodesic
	distance from the cloasest source to all the mesh vertices and returns the pointer to the farthest.
	Note: update the field Q() of the vertices
	*/
 void FartestVertex( MeshType & m,
									std::vector<typename MeshType::VertexPointer> & fro, 
									typename MeshType::VertexPointer & farthest,		 
									ScalarType & distance){					

									std::vector<typename MeshType::VertexPointer>::iterator fi; 
									std::vector<VertDist>fr;

									for( fi  = fro.begin(); fi != fro.end() ; ++fi)
										fr.push_back(VertDist(*fi,-1));
									farthest = Visit(m,fr,distance,false); 
					  }
	/*
	Given a mesh and  a  pointers to a vertex-source (source), assigns the approximated geodesic
	distance from the vertex-source to all the mesh vertices and returns the pointer to the farthest
	Note: update the field Q() of the vertices
	*/
 void FartestVertex( MeshType & m, 
									typename MeshType::VertexPointer seed,
									typename MeshType::VertexPointer & farthest,		 
									ScalarType & distance){	
	std::vector<typename MeshType::VertexPointer>  fro;
	fro.push_back( seed );
	typename MeshType::VertexPointer v0;
	FartestVertex(m,fro,v0,distance);
	farthest = v0;
}

/* 
	Same as FartestPoint but the returned pointer is to a border vertex
	Note: update the field Q() of the vertices
*/
 void FartestBVertex(MeshType & m,
										std::vector<typename MeshType::VertexPointer> & fro,  
										typename MeshType::VertexPointer & farthest,	     
										ScalarType & distance){

	std::vector<typename MeshType::VertexPointer>::iterator fi; 
	std::vector<VertDist>fr;

	for( fi  = fro.begin(); fi != fro.end() ; ++fi)
		fr.push_back(VertDist(*fi,-1));
	farthest =  Visit(m,fr,distance,true); 
}
/* 
	Same as FartestPoint but the returned pointer is to a border vertex
	Note: update the field Q() of the vertices
*/
 void FartestBVertex( MeshType & m, 
									typename MeshType::VertexPointer seed,
									typename MeshType::VertexPointer & farthest,		 
									ScalarType & distance){	
	std::vector<typename MeshType::VertexPointer>  fro;
	fro.push_back( seed );
	typename MeshType::VertexPointer v0;
	FartestBVertex(m,fro,v0,distance);
	farthest = v0;
 }

/* 
	Assigns to each vertex of the mesh its distance to the closest vertex on the border
	Note: update the field Q() of the vertices
*/
 void DistanceFromBorder(	MeshType & m,
									  typename MeshType::VertexPointer & v0,	 
										ScalarType & distance			 
					){
	std::vector<typename MeshType::VertexPointer> fro;
	MeshType::VertexIterator vi;
	MeshType::VertexPointer farthest;
	for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
		if( (*vi).IsB())
			fro.push_back(&(*vi));
	FartestVertex(m,fro,farthest,distance);
}

 };
};// end namespace