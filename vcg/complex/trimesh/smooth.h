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
/****************************************************************************
  History
$Log: not supported by cvs2svn $

/****************************************************************************/


#ifndef __VCGLIB__SMOOTH
#define __VCGLIB__SMOOTH

#include <vcg/space/point3.h>
#include <vcg/space/line3.h>
#include <vcg/container/simple_temporary_data.h>

namespace vcg
{

template<class FLT> 
class ScaleLaplacianInfo 
{
public:
	Point3<FLT> PntSum;
	FLT LenSum;
};

// Scale dependent laplacian smoothing [fujimori 95]
// Nuova versione, l'idea e'quella di usare anche gli angoli delle facce per pesare lo spostamento.
// 
// in pratica si sposta solo lungo la componente che e' parallela alla normale al vertice 
// (che si suppone esserci!!)
 

// Non ha bisogno della topologia
// Non fa assunzioni sull'ordinamento delle facce, ma vuole che i border flag ci siano!
//
// 

template<class MESH_TYPE>
void ScaleLaplacianSmooth(MESH_TYPE &m, int step, typename  MESH_TYPE::ScalarType delta)
{
	SimpleTempData<typename MESH_TYPE::VertContainer, ScaleLaplacianInfo<typename MESH_TYPE::ScalarType> > TD(m.vert);
	ScaleLaplacianInfo<typename MESH_TYPE::ScalarType> lpz;
	lpz.PntSum=typename  MESH_TYPE::CoordType(0,0,0);
	lpz.LenSum=0;
	TD.Start(lpz);
	typename  MESH_TYPE::FaceIterator fi;
	for(int i=0;i<step;++i)
	{
		typename  MESH_TYPE::VertexIterator vi;
		for(vi=m.vert.begin();vi!=m.vert.end();++vi)
			 TD[*vi]=lpz;
		typename MESH_TYPE::ScalarType a[3];
		for(fi=m.face.begin();fi!=m.face.end();++fi)if(!(*fi).IsD())
		{
			typename  MESH_TYPE::CoordType  mp=((*fi).V(0)->P() + (*fi).V(1)->P() + (*fi).V(2)->P())/3.0;
			typename  MESH_TYPE::CoordType  e0=((*fi).V(0)->P() - (*fi).V(1)->P()).Normalize(); 
			typename  MESH_TYPE::CoordType  e1=((*fi).V(1)->P() - (*fi).V(2)->P()).Normalize();
			typename  MESH_TYPE::CoordType  e2=((*fi).V(2)->P() - (*fi).V(0)->P()).Normalize();

			a[0]=AngleN(-e0,e2);
			a[1]=AngleN(-e1,e0);
			a[2]=AngleN(-e2,e1);
			//assert(fabs(M_PI -a[0] -a[1] -a[2])<0.0000001);

			for(int j=0;j<3;++j){
						typename  MESH_TYPE::CoordType dir= (mp-(*fi).V(j)->P()).Normalize();
						TD[(*fi).V(j)].PntSum+=dir*a[j];
						TD[(*fi).V(j)].LenSum+=a[j];
			}
		}		
		for(vi=m.vert.begin();vi!=m.vert.end();++vi)
				if(!(*vi).IsD() && TD[*vi].LenSum>0 )
				 (*vi).P() = (*vi).P() +  (TD[*vi].PntSum/TD[*vi].LenSum ) * delta;
				
	}		 			
	TD.Stop();
};

// Scale dependent laplacian smoothing [fujimori 95]
// Non ha bisogno della topologia
// Non fa assunzioni sull'ordinamento delle facce, ma vuole che i border flag ci siano!
//
// 
template<class MESH_TYPE>
void ScaleLaplacianSmoothOld(MESH_TYPE &m, int step, typename  MESH_TYPE::ScalarType delta)
{
	SimpleTempData<typename MESH_TYPE::VertContainer, ScaleLaplacianInfo<typename MESH_TYPE::ScalarType> > TD(m.vert);
	ScaleLaplacianInfo<typename MESH_TYPE::ScalarType> lpz;
	lpz.PntSum=typename  MESH_TYPE::CoordType(0,0,0);
	lpz.LenSum=0;
	TD.Start(lpz);
	typename  MESH_TYPE::FaceIterator fi;
	for(int i=0;i<step;++i)
	{
			typename  MESH_TYPE::VertexIterator vi;
			for(vi=m.vert.begin();vi!=m.vert.end();++vi)
				 TD[*vi]=lpz;

			for(fi=m.face.begin();fi!=m.face.end();++fi)if(!(*fi).IsD())
				for(int j=0;j<3;++j)
					if(!(*fi).IsB(j)) {
						typename  MESH_TYPE::CoordType edge= (*fi).V1(j)->P() -(*fi).V(j)->P();
						typename MESH_TYPE::ScalarType len=Norm(edge);
						edge/=len;
						TD[(*fi).V(j)].PntSum+=edge;
						TD[(*fi).V1(j)].PntSum-=edge;
						TD[(*fi).V(j)].LenSum+=len;
						TD[(*fi).V1(j)].LenSum+=len;
					}
				
			for(fi=m.face.begin();fi!=m.face.end();++fi)if(!(*fi).IsD())
				for(int j=0;j<3;++j)
					// se l'edge j e' di bordo si riazzera tutto e si riparte
					if((*fi).IsB(j)) {
							TD[(*fi).V(j)].PntSum=typename  MESH_TYPE::CoordType(0,0,0);
							TD[(*fi).V1(j)].PntSum=typename  MESH_TYPE::CoordType(0,0,0);
							TD[(*fi).V(j)].LenSum=0;
							TD[(*fi).V1(j)].LenSum=0;
					}
				

			for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
					for(int j=0;j<3;++j)
						if((*fi).IsB(j)) 
						{ 
							typename  MESH_TYPE::CoordType edge= (*fi).V1(j)->P() -(*fi).V(j)->P();
							typename MESH_TYPE::ScalarType len=Norm(edge);
							edge/=len;
							TD[(*fi).V(j)].PntSum+=edge;
							TD[(*fi).V1(j)].PntSum-=edge;
							TD[(*fi).V(j)].LenSum+=len;
							TD[(*fi).V1(j)].LenSum+=len;
						}

			for(vi=m.vert.begin();vi!=m.vert.end();++vi)
				if(!(*vi).IsD() && TD[*vi].LenSum>0 )
				 (*vi).P() = (*vi).P() + (TD[*vi].PntSum/TD[*vi].LenSum)*delta;
	}		 			
	TD.Stop();
};


template<class FLT> 
class LaplacianInfo 
{
public:
	Point3<FLT> sum;
	FLT cnt;
};

template<class MESH_TYPE>
void LaplacianSmooth(MESH_TYPE &m, int step,bool SmoothSelected=false)
{
	SimpleTempData<typename MESH_TYPE::VertContainer,LaplacianInfo<typename MESH_TYPE::ScalarType> > TD(m.vert);
  LaplacianInfo<typename MESH_TYPE::ScalarType> lpz;
	lpz.sum=typename  MESH_TYPE::CoordType(0,0,0);
	lpz.cnt=0;
	TD.Start(lpz);
	for(int i=0;i<step;++i)
	{
		typename  MESH_TYPE::VertexIterator vi;
		for(vi=m.vert.begin();vi!=m.vert.end();++vi)
			 TD[*vi]=lpz;

		typename  MESH_TYPE::FaceIterator fi;
		for(fi=m.face.begin();fi!=m.face.end();++fi)
			if(!(*fi).IsD()) 
				for(int j=0;j<3;++j)
					if(!(*fi).IsB(j)) 
						{
							TD[(*fi).V(j)].sum+=(*fi).V1(j)->P();
							TD[(*fi).V1(j)].sum+=(*fi).V(j)->P();
							++TD[(*fi).V(j)].cnt;
							++TD[(*fi).V1(j)].cnt;
					}

			// si azzaera i dati per i vertici di bordo
			for(fi=m.face.begin();fi!=m.face.end();++fi)
				if(!(*fi).IsD()) 
					for(int j=0;j<3;++j)
						if((*fi).IsB(j))
							{
								TD[(*fi).V(j)]=lpz;
								TD[(*fi).V1(j)]=lpz;
							}

			// se l'edge j e' di bordo si deve mediare solo con gli adiacenti
			for(fi=m.face.begin();fi!=m.face.end();++fi)
				if(!(*fi).IsD()) 
					for(int j=0;j<3;++j)
						if((*fi).IsB(j)) 
							{
								TD[(*fi).V(j)].sum+=(*fi).V1(j)->P();
								TD[(*fi).V1(j)].sum+=(*fi).V(j)->P();
								++TD[(*fi).V(j)].cnt;
								++TD[(*fi).V1(j)].cnt;
						}
	
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if(!(*vi).IsD() && TD[*vi].cnt>0 )
			if(!SmoothSelected || (*vi).IsS())
					(*vi).P()=TD[*vi].sum/TD[*vi].cnt;
	}
	 	
	TD.Stop();
};


template<class FLT> 
class HCSmoothInfo 
{
public:
	Point3<FLT> dif;
	Point3<FLT> sum;
	int cnt;
};
template<class MESH_TYPE>
void HCSmooth(MESH_TYPE &m, int step)
{
	typename MESH_TYPE::ScalarType beta=0.5;
	SimpleTempData<typename MESH_TYPE::VertContainer,HCSmoothInfo<typename MESH_TYPE::ScalarType> > TD(m.vert);
  HCSmoothInfo<typename MESH_TYPE::ScalarType> lpz;
	lpz.sum=typename  MESH_TYPE::CoordType(0,0,0);
	lpz.dif=typename  MESH_TYPE::CoordType(0,0,0);
	lpz.cnt=0;
	TD.Start(lpz);
	// First Loop compute the laplacian
	typename  MESH_TYPE::FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
		{
			for(int j=0;j<3;++j)
			{
				TD[(*fi).V(j)].sum+=(*fi).V1(j)->P();
				TD[(*fi).V1(j)].sum+=(*fi).V(j)->P();
				++TD[(*fi).V(j)].cnt;
				++TD[(*fi).V1(j)].cnt;
				// se l'edge j e' di bordo si deve sommare due volte
				if((*fi).IsB(j)) 
				{ 
					TD[(*fi).V(j)].sum+=(*fi).V1(j)->P(); 
					TD[(*fi).V1(j)].sum+=(*fi).V(j)->P(); 
					++TD[(*fi).V(j)].cnt;
					++TD[(*fi).V1(j)].cnt;
				}
			}
		}
	typename  MESH_TYPE::VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		 TD[*vi].sum/=(float)TD[*vi].cnt;
	
	// Second Loop compute average difference
	for(fi=m.face.begin();fi!=m.face.end();++fi)
		{
			for(int j=0;j<3;++j)
			{
				TD[(*fi).V(j)].dif +=TD[(*fi).V1(j)].sum-(*fi).V1(j)->P();
				TD[(*fi).V1(j)].dif+=TD[(*fi).V(j)].sum-(*fi).V(j)->P();
				// se l'edge j e' di bordo si deve sommare due volte
				if((*fi).IsB(j)) 
				{ 
					TD[(*fi).V(j)].dif +=TD[(*fi).V1(j)].sum-(*fi).V1(j)->P();
					TD[(*fi).V1(j)].dif+=TD[(*fi).V(j)].sum-(*fi).V(j)->P();
				}
			}
		}

	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		{
		 TD[*vi].dif/=(float)TD[*vi].cnt;
		 (*vi).P()=TD[*vi].sum -((TD[*vi].sum-(*vi).P()*beta) + TD[*vi].dif)*(1.f-beta);
		}
		 	
	TD.Stop();
};


template<class FLT> 
class QualitySmoothInfo 
{
public:
	FLT sum;
	int cnt;
};

template<class MESH_TYPE>
void LaplacianSmoothQuality(MESH_TYPE &m, int step,bool SmoothSelected=false)
{ 
	SimpleTempData<typename MESH_TYPE::VertContainer,QualitySmoothInfo<typename MESH_TYPE::ScalarType> > TD(m.vert);
  QualitySmoothInfo<typename MESH_TYPE::ScalarType> lpz;
	lpz.sum=0;
	lpz.cnt=0;
	TD.Start(lpz);
	for(int i=0;i<step;++i)
	{
		typename  MESH_TYPE::VertexIterator vi;
		for(vi=m.vert.begin();vi!=m.vert.end();++vi)
			 TD[*vi]=lpz;

		typename  MESH_TYPE::FaceIterator fi;
		for(fi=m.face.begin();fi!=m.face.end();++fi)
			if(!(*fi).IsD()) 
				for(int j=0;j<3;++j)
					if(!(*fi).IsB(j)) 
						{
							TD[(*fi).V(j)].sum+=(*fi).V1(j)->Q();
							TD[(*fi).V1(j)].sum+=(*fi).V(j)->Q();
							++TD[(*fi).V(j)].cnt;
							++TD[(*fi).V1(j)].cnt;
					}

			// si azzaera i dati per i vertici di bordo
			for(fi=m.face.begin();fi!=m.face.end();++fi)
				if(!(*fi).IsD()) 
					for(int j=0;j<3;++j)
						if((*fi).IsB(j))
							{
								TD[(*fi).V(j)]=lpz;
								TD[(*fi).V1(j)]=lpz;
							}

			// se l'edge j e' di bordo si deve mediare solo con gli adiacenti
			for(fi=m.face.begin();fi!=m.face.end();++fi)
				if(!(*fi).IsD()) 
					for(int j=0;j<3;++j)
						if((*fi).IsB(j)) 
							{
								TD[(*fi).V(j)].sum+=(*fi).V1(j)->Q();
								TD[(*fi).V1(j)].sum+=(*fi).V(j)->Q();
								++TD[(*fi).V(j)].cnt;
								++TD[(*fi).V1(j)].cnt;
						}

	//typename  MESH_TYPE::VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if(!(*vi).IsD() && TD[*vi].cnt>0 )
			if(!SmoothSelected || (*vi).IsS())
					(*vi).Q()=TD[*vi].sum/TD[*vi].cnt;
	}
		 	
	TD.Stop();
};
template<class MESH_TYPE>
void LaplacianSmoothNormals(MESH_TYPE &m, int step,bool SmoothSelected=false)
{
	SimpleTempData<typename MESH_TYPE::VertContainer,LaplacianInfo<typename MESH_TYPE::ScalarType> > TD(m.vert);
  LaplacianInfo<typename MESH_TYPE::ScalarType> lpz;
	lpz.sum=typename  MESH_TYPE::CoordType(0,0,0);
	lpz.cnt=0;
	TD.Start(lpz);
	for(int i=0;i<step;++i)
	{
		typename  MESH_TYPE::VertexIterator vi;
		for(vi=m.vert.begin();vi!=m.vert.end();++vi)
			 TD[*vi]=lpz;

		typename  MESH_TYPE::FaceIterator fi;
		for(fi=m.face.begin();fi!=m.face.end();++fi)
			if(!(*fi).IsD()) 
				for(int j=0;j<3;++j)
					if(!(*fi).IsB(j)) 
						{
							TD[(*fi).V(j)].sum+=(*fi).V1(j)->N();
							TD[(*fi).V1(j)].sum+=(*fi).V(j)->N();
							++TD[(*fi).V(j)].cnt;
							++TD[(*fi).V1(j)].cnt;
					}

			// si azzaera i dati per i vertici di bordo
			for(fi=m.face.begin();fi!=m.face.end();++fi)
				if(!(*fi).IsD()) 
					for(int j=0;j<3;++j)
						if((*fi).IsB(j))
							{
								TD[(*fi).V(j)]=lpz;
								TD[(*fi).V1(j)]=lpz;
							}

			// se l'edge j e' di bordo si deve mediare solo con gli adiacenti
			for(fi=m.face.begin();fi!=m.face.end();++fi)
				if(!(*fi).IsD()) 
					for(int j=0;j<3;++j)
						if((*fi).IsB(j)) 
							{
								TD[(*fi).V(j)].sum+=(*fi).V1(j)->N();
								TD[(*fi).V1(j)].sum+=(*fi).V(j)->N();
								++TD[(*fi).V(j)].cnt;
								++TD[(*fi).V1(j)].cnt;
						}

	//typename  MESH_TYPE::VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if(!(*vi).IsD() && TD[*vi].cnt>0 )
			if(!SmoothSelected || (*vi).IsS())
					(*vi).N()=TD[*vi].sum/TD[*vi].cnt;
	}
		 	
	TD.Stop();
};

// Smooth solo lungo la direzione di vista
	// alpha e' compreso fra 0(no smoot) e 1 (tutto smoot)
  // Nota che se smootare il bordo puo far fare bandierine.
template<class MESH_TYPE>
void DepthSmooth(MESH_TYPE &m,
				 const typename  MESH_TYPE::CoordType & viewpoint,
				 const typename  MESH_TYPE::ScalarType alpha,
				 int step, bool SmoothBorder=false )
{
	typedef typename  MESH_TYPE::CoordType v_type;
	typedef typename MESH_TYPE::ScalarType    s_type;


	//const typename  MESH_TYPE::CoordType viewpoint;
	//const typename MESH_TYPE::ScalarType alpha;

	SimpleTempData<typename MESH_TYPE::VertContainer,LaplacianInfo<typename MESH_TYPE::ScalarType> > TD(m.vert);
	LaplacianInfo<typename MESH_TYPE::ScalarType> lpz;
	lpz.sum=typename  MESH_TYPE::CoordType(0,0,0);
	lpz.cnt=0;
	TD.Start(lpz);
	for(int i=0;i<step;++i)
	{
		typename  MESH_TYPE::VertexIterator vi;
		for(vi=m.vert.begin();vi!=m.vert.end();++vi)
			 TD[*vi]=lpz;

		typename  MESH_TYPE::FaceIterator fi;
		for(fi=m.face.begin();fi!=m.face.end();++fi)
			if(!(*fi).IsD()) 
				for(int j=0;j<3;++j)
					if(!(*fi).IsB(j)) 
						{
							TD[(*fi).V(j)].sum+=(*fi).V1(j)->Supervisor_P();
							TD[(*fi).V1(j)].sum+=(*fi).V(j)->Supervisor_P();
							++TD[(*fi).V(j)].cnt;
							++TD[(*fi).V1(j)].cnt;
					}

			// si azzaera i dati per i vertici di bordo
			for(fi=m.face.begin();fi!=m.face.end();++fi)
				if(!(*fi).IsD()) 
					for(int j=0;j<3;++j)
						if((*fi).IsB(j))
							{
								TD[(*fi).V(j)]=lpz;
								TD[(*fi).V1(j)]=lpz;
							}

			// se l'edge j e' di bordo si deve mediare solo con gli adiacenti
     if(SmoothBorder)
			for(fi=m.face.begin();fi!=m.face.end();++fi)
				if(!(*fi).IsD()) 
					for(int j=0;j<3;++j)
						if((*fi).IsB(j)) 
							{
								TD[(*fi).V(j)].sum+=(*fi).V1(j)->Supervisor_P();
								TD[(*fi).V1(j)].sum+=(*fi).V(j)->Supervisor_P();
								++TD[(*fi).V(j)].cnt;
								++TD[(*fi).V1(j)].cnt;
						}
	
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if(!(*vi).IsD() && TD[*vi].cnt>0 )
			{
				v_type np = TD[*vi].sum/TD[*vi].cnt;
				v_type d = (*vi).Supervisor_P() - viewpoint; d.Normalize();
				s_type s = d * ( np - (*vi).Supervisor_P() );
				(*vi).Supervisor_P() += d * (s*alpha);
			}
	}
	 	
	TD.Stop();
}

}		// End namespace vcg

#endif //  VCG_SMOOTH