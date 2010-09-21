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
#include <vcg/math/base.h>
#include <vcg/container/simple_temporary_data.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/complex/trimesh/update/quality.h>
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
struct EuclideanDistance{
	typedef typename MeshType::VertexType VertexType;
	typedef typename MeshType::ScalarType  ScalarType;

	EuclideanDistance(){}
	ScalarType operator()(const VertexType * v0, const VertexType * v1) const
	{return vcg::Distance(v0->cP(),v1->cP());}
};

template <class MeshType, class DistanceFunctor = EuclideanDistance<MeshType> >
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
		TempData(const ScalarType & d_){d=d_;source = NULL;}
		ScalarType d;
		VertexPointer source;//closest source

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

		ScalarType ew_c  = DistanceFunctor()(pw,curr);
		ScalarType ew_w1 = DistanceFunctor()(pw,pw1);
		ScalarType ec_w1 = DistanceFunctor()(pw1,curr);
		CoordType w_c =  (pw->cP()-curr->cP()).Normalize() * ew_c;
		CoordType w_w1 = (pw->cP() - pw1->cP()).Normalize() * ew_w1;
		CoordType w1_c =  (pw1->cP() - curr->cP()).Normalize() * ec_w1;

		ScalarType	alpha,alpha_, beta,beta_,theta,h,delta,s,a,b;

		alpha = acos((w_c.dot(w1_c))/(ew_c*ec_w1));
		s = (d_curr + d_pw1+ec_w1)/2;
		a = s/ec_w1;
		b = a*s;
    alpha_ = 2*acos ( std::min<ScalarType>(1.0,sqrt(  (b- a* d_pw1)/d_curr)));

		if ( alpha+alpha_ > M_PI){
			curr_d = d_curr + ew_c;		
			}else
			{
        beta_ = 2*acos ( std::min<ScalarType>(1.0,sqrt(  (b- a* d_curr)/d_pw1)));
				beta  = acos((w_w1).dot(-w1_c)/(ew_w1*ec_w1));

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
		 bool farthestOnBorder = false,
		 typename MeshType::template PerVertexAttributeHandle<VertexPointer> * sources = NULL
		 )
{
	bool isLeaf;
	std::vector<VertDist> frontier;
	VertexIterator ii;
	std::list<VertexPointer> children;
  VertexPointer curr,farthest=0,pw1;
	typename std::list<VertexPointer>::iterator is;
	std::deque<VertexPointer> leaves;
	std::vector<VertDist> _frontier;
	ScalarType unreached  = std::numeric_limits<ScalarType>::max();

	std::vector <std::pair<VertexPointer,ScalarType> > expansion;
	typename std::vector <VertDist >::iterator ifr;
  face::VFIterator<FaceType> x;
	VertexPointer pw;

	//Requirements
	assert(m.HasVFTopology());
	assert(!_inputfrontier.empty());

	TempDataType  * TD;
	TD = new TempDataType(m.vert,unreached);

	for(ifr = _inputfrontier.begin(); ifr != _inputfrontier.end(); ++ifr){
		(*TD)[(*ifr).v].d = 0.0;	
		(*ifr).d = 0.0;
		(*TD)[(*ifr).v].source  = (*ifr).v;
		frontier.push_back(VertDist((*ifr).v,0.0));
		}
		// initialize Heap
		make_heap(frontier.begin(),frontier.end(),pred());	

		ScalarType curr_d,d_curr = 0.0,d_heap;
		VertexPointer curr_s = NULL;
		max_distance=0.0;
		typename std::vector<VertDist >:: iterator iv;

 		while(!frontier.empty())
			{  
				pop_heap(frontier.begin(),frontier.end(),pred());
				curr = (frontier.back()).v;
				curr_s = (*TD)[curr].source;
				if(sources!=NULL) 
 					(*sources)[curr] = curr_s;
				d_heap = (frontier.back()).d;
				frontier.pop_back();

				assert((*TD)[curr].d <= d_heap); 
				assert(curr_s != NULL);
				if((*TD)[curr].d < d_heap )// a vertex whose distance has been improved after it was inserted in the queue
					continue; 
				assert((*TD)[curr].d == d_heap); 

				d_curr =  (*TD)[curr].d;

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

						const ScalarType & d_pw1  =  (*TD)[pw1].d;
						{

							const ScalarType inter  = DistanceFunctor()(curr,pw1);//(curr->P() - pw1->P()).Norm();
							const ScalarType tol = (inter + d_curr + d_pw1)*.0001f;

							if (	((*TD)[pw1].source != (*TD)[curr].source)||// not the same source
									(inter + d_curr < d_pw1  +tol   ) ||
									(inter + d_pw1  < d_curr +tol  ) ||
									(d_curr + d_pw1  < inter +tol  )   // triangular inequality
									 )  
								curr_d = d_curr + DistanceFunctor()(pw,curr);//(pw->P()-curr->P()).Norm();
							else
								curr_d = Distance(pw,pw1,curr,d_pw1,d_curr);
							
						}
				 
					if((*TD)[(pw)].d > curr_d){
							(*TD)[(pw)].d = curr_d; 
							(*TD)[pw].source = curr_s;
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
  assert(farthest);
 	return farthest;

 }
	

public:
	/*
	Given a mesh and  a vector of pointers to vertices (sources), assigns the approximated geodesic
	distance from the cloasest source to all the mesh vertices and returns the pointer to the farthest.
	Note: update the field Q() of the vertices
	*/
 static bool FarthestVertex( MeshType & m,
									std::vector<VertexPointer> & fro, 
									VertexPointer & farthest,		 
									ScalarType & distance,
									typename MeshType::template PerVertexAttributeHandle<VertexPointer> * sources = NULL){					

									typename std::vector<VertexPointer>::iterator fi; 
									std::vector<VertDist>fr;
									if(fro.empty())	return false;
										
									for( fi  = fro.begin(); fi != fro.end() ; ++fi)
										fr.push_back(VertDist(*fi,0.0));
									farthest = Visit(m,fr,distance,false,sources); 
									return true;
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
										ScalarType & distance,
										typename MeshType::template PerVertexAttributeHandle<VertexPointer> * sources = NULL
										){

	typename std::vector<VertexPointer>::iterator fi; 
	std::vector<VertDist>fr;

	for( fi  = fro.begin(); fi != fro.end() ; ++fi)
		fr.push_back(VertDist(*fi,-1));
	farthest =  Visit(m,fr,distance,true,sources); 
}
/* 
	Same as FarthestPoint but the returned pointer is to a border vertex
	Note: update the field Q() of the vertices
*/
 static void FarthestBVertex(	MeshType & m, 
								VertexPointer seed,
								VertexPointer & farthest,		 
								ScalarType & distance,
								typename MeshType::template PerVertexAttributeHandle<VertexPointer> * sources = NULL){	
	std::vector<VertexPointer>  fro;
	fro.push_back( seed );
	VertexPointer v0;
	FarthestBVertex(m,fro,v0,distance,sources);
	farthest = v0;
 }

/* 
	Assigns to each vertex of the mesh its distance to the closest vertex on the border
	Note: update the field Q() of the vertices
*/
 static bool DistanceFromBorder(	MeshType & m,
									ScalarType & distance,
									typename MeshType::template PerVertexAttributeHandle<VertexPointer> * sources = NULL
					){
	std::vector<VertexPointer> fro;
	VertexIterator vi;
	VertexPointer farthest;
	for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
		if( (*vi).IsB())
			fro.push_back(&(*vi));
	if(fro.empty()) return false;
	
	tri::UpdateQuality<CMeshO>::VertexConstant(m,0);
		
	return FarthestVertex(m,fro,farthest,distance,sources);
}

 };
};// end namespace tri
};// end namespace vcg
