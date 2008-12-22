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
#include <vcg/simplex/face/topology.h>
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
namespace tri{
template <class MeshType>
class Geo{

	public:

	typedef typename MeshType::VertexType VertexType;
	typedef typename MeshType::VertexIterator VertexIterator;
	typedef typename MeshType::VertexPointer VertexPointer;
	typedef typename MeshType::FaceType  FaceType;
	typedef typename MeshType::CoordType  CoordType;
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
	struct TempData{
		TempData(){}
		TempData(const ScalarType & d_){d=d_;visited=false;}
		ScalarType d;
		bool visited;
		};

	typedef SimpleTempData<std::vector<VertexType>, TempData >  TempDataType;
	//TempDataType  * TD;
	
	
	struct pred: public std::binary_function<VertDist,VertDist,bool>{
		pred(){};
			bool operator()(const VertDist& v0, const VertDist& v1) const
				{return (v0.d > v1.d);}
		};
	struct pred_addr: public std::binary_function<VertDist,VertDist,bool>{
		pred_addr(){};
			bool operator()(const VertDist& v0, const VertDist& v1) const
				{return (v0.v > v1.v);}
		};

	//************** calcolo della distanza di pw in base alle distanze note di pw1 e curr
	//************** sapendo che (curr,pw,pw1) e'una faccia della mesh
	//************** (vedi figura in file distance.gif)
	static ScalarType Distance(const VertexPointer &pw,
										   const VertexPointer &pw1,
										   const VertexPointer &curr,
										   const ScalarType &d_pw1,
										   const ScalarType &d_curr)
	{
		ScalarType curr_d=0;
		CoordType w_c = pw->cP()- curr->cP();
		CoordType w_w1 = pw->cP()- pw1->cP();
		CoordType w1_c = pw1->cP()- curr->cP();

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
		return (curr_d);
	}

	/*
	  starting from the seeds, it assign a distance value to each vertex. The distance of a vertex is its 
	  approximated geodesic distance to the closest seeds.
	  This is function is not meant to be called (although is not prevented). Instead, it is invoked by 
	  wrapping function.
	*/
 	static  VertexPointer Visit( 
		 MeshType & m,
		 std::vector<VertDist> & _inputfrontier,
		 ScalarType & max_distance,
		 bool farthestOnBorder = false
		 )
{
	bool isLeaf,toQueue;
	std::vector<VertDist> frontier;
	VertexIterator ii;
	std::list<VertexPointer> children;
	VertexPointer curr,farthest,pw1;	
	typename std::list<VertexPointer>::iterator is;
	std::deque<VertexPointer> leaves;
	std::vector<VertDist> _frontier;

	std::vector <std::pair<VertexPointer,ScalarType> > expansion;
	typename std::vector <VertDist >::iterator ifr;
	face::VFIterator<FaceType> x;int k;
	VertexPointer pw;

	//Requirements
	assert(m.HasVFTopology());
	assert(!_inputfrontier.empty());

	ScalarType unreached  = std::numeric_limits<ScalarType>::max();
	TempDataType  * TD;
	TD = new TempDataType(m.vert,unreached);

	bool singleSource = _inputfrontier.size() ==1;
	for(ifr = _inputfrontier.begin(); ifr != _inputfrontier.end(); ++ifr){
		(*TD)[(*ifr).v].d = 0.0;	
		(*ifr).d = 0.0;
		frontier.push_back(VertDist((*ifr).v,0.0));
		}
		// initialize Heap
		make_heap(frontier.begin(),frontier.end(),pred());	

	//	for(ifr = frontier.begin(); ifr != frontier.end(); ++ifr)
	//		printf("%d %f\n",(*ifr).v,(*ifr).d);


		ScalarType curr_d,d_curr = 0.0,d_heap;
		max_distance=0.0;
		typename std::vector<VertDist >:: iterator iv;

 		while(!frontier.empty())
			{ //printf("size: %d\n", frontier.size());
				 
				pop_heap(frontier.begin(),frontier.end(),pred());
				curr = (frontier.back()).v;
				d_heap = (frontier.back()).d;
				frontier.pop_back();

				float fff = (*TD)[curr].d;
				bool vis = (*TD)[curr].visited;

				assert((*TD)[curr].d <= d_heap); 

				if((*TD)[curr].d < d_heap )
					continue; 
				assert((*TD)[curr].d == d_heap); 

			//	printf("extracted %d %f\n",curr,fff);
				d_curr =  (*TD)[curr].d;

				(*TD)[curr].visited = true;
				isLeaf = (!farthestOnBorder || curr->IsB());

				face::VFIterator<FaceType> x;int k;

			
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

					if(singleSource){
						const ScalarType & d_pw1  =  (*TD)[pw1].d;

						if((*TD)[curr].d==0.0)// numerical 
							curr_d = (pw->P()-curr->P()).Norm();
						else
						if(d_pw1==0.0)// numerical 
							curr_d = (pw->P()-pw1->P()).Norm();
						else
						if(  (*TD)[pw1].d == unreached ){
						 		curr_d =  d_curr + (pw->P()-curr->P()).Norm();
						}
						else{

							ScalarType inter  = (curr->P() - pw1->P()).Norm();

							if (	(inter + d_curr < d_pw1 + 0.01 ) ||
									(inter + d_pw1  < d_curr + 0.01 ) ||
									(d_curr + d_pw1  < inter + 0.01 ))
								curr_d = d_curr + (pw->P()-pw1->P()).Norm();// triangular inequality
							else{
								//curr_d = d_pw1 + (pw->P()-pw1->P()).Norm();
								curr_d = Distance(pw,pw1,curr,d_pw1,d_curr);
			
							}
							
						}
					}else{
							curr_d = d_curr + (pw->P()-pw1->P()).Norm();
					}

					//printf("%f %f \n",
					//	curr_d,
					//	d_curr);
					// queue if the estimation is lower or if it is its first
					toQueue = ((*TD)[(pw)].d > curr_d) ;
				 
					if(toQueue ){
	//					printf("push: %d %f\n",pw,curr_d);
							(*TD)[(pw)].d = curr_d; 
						//	printf("from %f estim. %f\n",d_curr,curr_d);
							frontier.push_back(VertDist(pw,curr_d));
							push_heap(frontier.begin(),frontier.end(),pred());
							  
			 			} 
					if(isLeaf){
						if(d_curr > max_distance){
							max_distance = d_curr;
							farthest = curr;
							}
						}
					}
	}// end while

	// scrivi le distanze sul campo qualita' (nn: farlo parametrico)
	VertexIterator vi;
	for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
		(*vi).Q() =  (*TD)[&(*vi)].d; 

	delete TD;
 	return farthest;

 }
	

public:
	/*
	Given a mesh and  a vector of pointers to vertices (sources), assigns the approximated geodesic
	distance from the cloasest source to all the mesh vertices and returns the pointer to the farthest.
	Note: update the field Q() of the vertices
	*/
 static void FarthestVertex( MeshType & m,
									std::vector<VertexPointer> & fro, 
									VertexPointer & farthest,		 
									ScalarType & distance){					

									typename std::vector<VertexPointer>::iterator fi; 
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
 static void FarthestVertex( MeshType & m, 
									VertexPointer seed,
									VertexPointer & farthest,		 
									ScalarType & distance){	
	std::vector<VertexPointer>  fro;
	fro.push_back( seed );
	VertexPointer v0;
	FarthestVertex(m,fro,v0,distance);
	farthest = v0;
}

/* 
	Same as FarthestPoint but the returned pointer is to a border vertex
	Note: update the field Q() of the vertices
*/
 static void FarthestBVertex(MeshType & m,
										std::vector<VertexPointer> & fro,  
										VertexPointer & farthest,	     
										ScalarType & distance){

	typename std::vector<VertexPointer>::iterator fi; 
	std::vector<VertDist>fr;

	for( fi  = fro.begin(); fi != fro.end() ; ++fi)
		fr.push_back(VertDist(*fi,-1));
	farthest =  Visit(m,fr,distance,true); 
}
/* 
	Same as FarthestPoint but the returned pointer is to a border vertex
	Note: update the field Q() of the vertices
*/
 static void FarthestBVertex( MeshType & m, 
									VertexPointer seed,
									VertexPointer & farthest,		 
									ScalarType & distance){	
	std::vector<VertexPointer>  fro;
	fro.push_back( seed );
	VertexPointer v0;
	FarthestBVertex(m,fro,v0,distance);
	farthest = v0;
 }

/* 
	Assigns to each vertex of the mesh its distance to the closest vertex on the border
	Note: update the field Q() of the vertices
*/
 static void DistanceFromBorder(	MeshType & m,
										ScalarType & distance			 
					){
	std::vector<VertexPointer> fro;
	VertexIterator vi;
	VertexPointer farthest;
	for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
		if( (*vi).IsB())
			fro.push_back(&(*vi));
	FarthestVertex(m,fro,farthest,distance);
}

 };
};// end namespace tri
};// end namespace vcg
