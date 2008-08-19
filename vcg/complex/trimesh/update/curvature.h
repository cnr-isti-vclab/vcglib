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
Revision 1.8  2008/05/14 10:03:29  ganovelli
Point3f->Coordtype

Revision 1.7  2008/04/23 16:37:15  onnis
VertexCurvature method added.

Revision 1.6  2008/04/04 10:26:12  cignoni
Cleaned up names, now Kg() gives back Gaussian Curvature (k1*k2), while Kh() gives back Mean Curvature 1/2(k1+k2)

Revision 1.5  2008/03/25 11:00:56  ganovelli
fixed  bugs sign of principal direction and mean curvature value

Revision 1.4  2008/03/17 11:29:59  ganovelli
taubin and desbrun estimates added (-> see vcg/simplex/vertexplus/component.h [component_ocf.h|component_occ.h ]

Revision 1.3  2006/02/27 18:02:57  ponchio
Area -> doublearea/2

added some typename

Revision 1.2  2005/10/25 09:17:41  spinelli
correct IsBorder

Revision 1.1  2005/02/22 16:40:29  ganovelli
created. This version writes the gaussian curvature on the Q() member of
the vertex

/****************************************************************************/

#ifndef VCGLIB_UPDATE_CURVATURE_
#define VCGLIB_UPDATE_CURVATURE_

#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/math/base.h>
#include <vcg/math/matrix.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/simplex/face/jumping_pos.h>
#include <vcg/container/simple_temporary_data.h>
#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/point_sampling.h>
#include <vcg/complex/trimesh/append.h>
#include <vcg/complex/intersection.h>
#include <vcg/complex/trimesh/inertia.h>

#include <wrap/io_trimesh/export_PLY.h>

namespace vcg {
namespace tri {

/// \ingroup trimesh 

/// \headerfile curvature.h vcg/complex/trimesh/update/curvature.h

/// \brief Management, updating and computation of per-vertex and per-face normals.
/** 
This class is used to compute or update the normals that can be stored in the vertex or face component of a mesh.
*/

template <class MeshType>
class UpdateCurvature
{

public:
	typedef typename MeshType  MeshType;
	typedef typename MeshType::FaceType FaceType;
	typedef typename MeshType::FacePointer FacePointer;
	typedef typename MeshType::FaceIterator FaceIterator;
	typedef typename MeshType::VertexIterator VertexIterator;
	typedef typename MeshType::VertContainer VertContainer;
	typedef typename MeshType::VertexType VertexType;
	typedef typename MeshType::VertexPointer VertexPointer;
	typedef vcg::face::VFIterator<FaceType> VFIteratorType;
	typedef typename MeshType::CoordType CoordType;
	typedef typename CoordType::ScalarType ScalarType;

	
private:
		typedef struct AdjVertex {
			VertexType * vert;
			float doubleArea;
 			bool isBorder;
		};
		

public:
	/// \brief Compute principal direction and magniuto of curvature.
	
	/** 
	 Based on the paper  <a href="http://mesh.caltech.edu/taubin/publications/taubin-iccv95b.pdf">  <em> "Estimating the Tensor of Curvature of a Surface from a Polyhedral Approximation" </em> </a>
	*/
		static void PrincipalDirections(MeshType &m) {

			assert(m.HasVFTopology());

			vcg::tri::UpdateNormals<MeshType>::PerVertexNormalized(m);
			vcg::tri::UpdateFlags<MeshType>::VertexBorderFromFace(m);

			VertexIterator vi; 
			for (vi =m.vert.begin(); vi !=m.vert.end(); ++vi) {
				if ( ! (*vi).IsD() && (*vi).VFp() != NULL) {

					VertexType * central_vertex = &(*vi);

					std::vector<float> weights;
					std::vector<AdjVertex> vertices_dup,vertices;

					assert((*vi).VFp() != NULL);
					vcg::face::JumpingPos<FaceType> pos((*vi).VFp(), central_vertex);
					float totalDoubleAreaSize = 0.0f;

					FaceType * startf = pos.F();
					FaceType* tempF;
					int hh = 0;
					do 
					{ hh++;
						AdjVertex v;

						pos.FlipE();
						v.vert = pos.VFlip();			
						v.doubleArea = vcg::DoubleArea(*pos.F());
						vertices_dup.push_back(v);

						pos.FlipE();
						v.vert = pos.VFlip();			
						v.doubleArea = vcg::DoubleArea(*pos.F());
						vertices_dup.push_back(v);

						pos.NextFE();
						tempF = pos.F();
					} 
					while(tempF != startf);	

					AdjVertex v;
					for(int i = 1 ; i <= vertices_dup.size();  )
					{
							v.vert = vertices_dup[(i)%vertices_dup.size()].vert;
							v.doubleArea = vertices_dup[i%vertices_dup.size()].doubleArea ;
						if( vertices_dup[(i)%vertices_dup.size()].vert == vertices_dup[(i+1)%vertices_dup.size()].vert){
							v.doubleArea += vertices_dup[(i+1)%vertices_dup.size()].doubleArea;
							i+=2;
						}else 
							++i;

						totalDoubleAreaSize+=v.doubleArea;
						vertices.push_back(v);
					}

					for (int i = 0; i < vertices.size(); ++i)  
							weights.push_back(vertices[i].doubleArea / totalDoubleAreaSize);

					Matrix33<ScalarType> Tp;
					for (int i = 0; i < 3; ++i)
						Tp[i][i] = 1.0f - powf(central_vertex->cN()[i],2);
					Tp[0][1] = Tp[1][0] = -1.0f * (central_vertex->cN()[0] * central_vertex->cN()[1]);
					Tp[1][2] = Tp[2][1] = -1.0f * (central_vertex->cN()[1] * central_vertex->cN()[2]);
					Tp[0][2] = Tp[2][0] = -1.0f * (central_vertex->cN()[0] * central_vertex->cN()[2]);

					Matrix33<ScalarType> tempMatrix;
					Matrix33<ScalarType> M;
					M.SetZero();
					for (int i = 0; i < vertices.size(); ++i) {
						CoordType edge = (central_vertex->cP() - vertices[i].vert->cP());
						float curvature = (2.0f * (central_vertex->cN() * edge) ) / edge.SquaredNorm();
						CoordType T = (Tp*edge).Normalize()*(-1.0); // -1.0 useless, just to conform the paper
						tempMatrix.ExternalProduct(T,T);
						M += tempMatrix * weights[i] * curvature ;
					}

					CoordType W;
					CoordType e1(1.0f,0.0f,0.0f);
					if ((e1 - central_vertex->cN()).SquaredNorm() > (e1 + central_vertex->cN()).SquaredNorm())
						W = e1 - central_vertex->cN();
					else 
						W = e1 + central_vertex->cN();
					W.Normalize();

					Matrix33<ScalarType> Q;
					Q.SetIdentity();
					tempMatrix.ExternalProduct(W,W);
					Q -= tempMatrix * 2.0f;

					Matrix33<ScalarType> Qt(Q);
					Qt.Transpose();
					Matrix33<ScalarType> QtMQ = (Qt * M * Q);

					CoordType bl = Q.GetColumn(0);
					CoordType T1 = Q.GetColumn(1);
					CoordType T2 = Q.GetColumn(2);

					float s,c;
					// Gabriel Taubin hint and Valentino Fiorin impementation
					float qt21 = QtMQ[2][1];
					float qt12 = QtMQ[1][2];


					float alpha = QtMQ[1][1]-QtMQ[2][2];
					float beta  = QtMQ[2][1];

					float h[2];
					float delta = sqrtf(4.0f*powf(alpha, 2) +16.0f*powf(beta, 2));
					h[0] = (2.0f*alpha + delta) / (2.0f*beta);
					h[1] = (2.0f*alpha - delta) / (2.0f*beta);

					float t[2];
					float best_c, best_s;
					float min_error = std::numeric_limits<ScalarType>::infinity();
					for (int i=0; i<2; i++)
					{
						delta = sqrtf(powf(h[1], 2) + 4.0f);
						t[0] = (h[i]+delta) / 2.0f;
						t[1] = (h[i]-delta) / 2.0f;

						for (int j=0; j<2; j++)
						{
							float squared_t = powf(t[j], 2);
							float denominator = 1.0f + squared_t;
							s = (2.0f*t[j])		/ denominator;
							c = (1-squared_t) / denominator;

							float approximation = c*s*alpha + (powf(c, 2) - powf(s, 2))*beta;
							float angle_similarity = fabs(acosf(c)/asinf(s));
							float error = fabs(1.0f-angle_similarity)+fabs(approximation);
							if (error<min_error)
							{
								min_error = error;
								best_c = c;
								best_s = s;
							}
						}
					}
					c = best_c;
					s = best_s;

					vcg::Matrix33<ScalarType> minor22(QtMQ);
					// clean up
					minor22[0][0] = minor22[0][1] = minor22[0][2] = 0.0;
					minor22[0][0] = minor22[1][0] = minor22[2][0] = 0.0;

					vcg::Matrix33<ScalarType> S; S.SetIdentity();
					S[1][1] = S[2][2] = c;
					S[1][2] = s;
					S[2][1] = -1.0f * s;

					vcg::Matrix33<ScalarType>  St (S);
					St.Transpose();					

					vcg::Matrix33<ScalarType>  StMS(St * minor22 * S);

					float Principal_Curvature1 = (3.0f * StMS[1][1]) - StMS[2][2];
					float Principal_Curvature2 = (3.0f * StMS[2][2]) - StMS[1][1];

					CoordType Principal_Direction1 = T1 * c - T2 * s;
					CoordType Principal_Direction2 = T1 * s + T2 * c; 

					(*vi).PD1() = Principal_Direction1;
					(*vi).PD2() = Principal_Direction2;
					(*vi).K1() = Principal_Curvature1;
					(*vi).K2() = Principal_Curvature2;
					}
			}
		}
	 


  class AreaData
  {
  public:
    float A;
  };

	/** Curvature meseaure as described in the paper: 
	Robust principal curvatures on Multiple Scales, Yong-Liang Yang, Yu-Kun Lai, Shi-Min Hu Helmut Pottmann
	SGP 2004
	If pointVSfaceInt==true the covariance is computed by montecarlo sampling on the mesh (faster)
	If pointVSfaceInt==false the covariance is computed by (analytic)integration over the surface (slower)
	*/

	typedef vcg::GridStaticPtr	<typename MeshType::FaceType, typename MeshType::ScalarType >		MeshGridType;
	typedef vcg::GridStaticPtr	<typename MeshType::VertexType, typename MeshType::ScalarType >		PointsGridType;

		static void PrincipalDirectionsPCA(MeshType &m, ScalarType r, bool pointVSfaceInt = true) {
			std::vector<MeshType:: VertexType*> closests;
			std::vector<MeshType::ScalarType> distances;
			std::vector<MeshType::CoordType> points;
			VertexIterator vi;
			ScalarType area;
			MeshType tmpM;
			std::vector<typename MeshType::CoordType>::iterator ii;
			vcg::tri::TrivialSampler<MeshType> vs;

			MeshGridType mGrid;
			PointsGridType pGrid;

			// Fill the grid used
			if(pointVSfaceInt){ 
					area = Stat<MeshType>::ComputeMeshArea(m);
					vcg::tri::SurfaceSampling<MeshType,vcg::tri::TrivialSampler<MeshType> >::Montecarlo(m,vs,1000 * area / (2*M_PI*r*r )); 
					vi = vcg::tri::Allocator<MeshType>::AddVertices(tmpM,m.vert.size());
					for(int y  = 0; y <   m.vert.size(); ++y,++vi)  (*vi).P() =  m.vert[y].P(); 
					pGrid.Set(tmpM.vert.begin(),tmpM.vert.end());
				}	else{	mGrid.Set(m.face.begin(),m.face.end()); }

				for(vi  = m.vert.begin(); vi != m.vert.end(); ++vi){
						vcg::Matrix33<ScalarType> A,eigenvectors;
						vcg::Point3<ScalarType> bp,eigenvalues;
						int nrot;

						// sample the neighborhood
						if(pointVSfaceInt)
						{ 
							vcg::trimesh::GetInSphereVertex<
								MeshType,
								PointsGridType,std::vector<MeshType::VertexType*>,
								std::vector<MeshType::ScalarType>,
								std::vector<MeshType::CoordType> >(tmpM,pGrid,  (*vi).cP(),r ,closests,distances,points);

							A.Covariance(points,bp);
							A*=area*area/1000;
						}
					else{
						IntersectionBallMesh<MeshType,ScalarType>( m ,vcg::Sphere3<ScalarType>((*vi).cP(),r),tmpM );
						vcg::Point3<ScalarType> _bary;
						vcg::tri::Inertia<MeshType>::Covariance(tmpM,_bary,A);
					}

					Jacobi(A,  eigenvalues , eigenvectors, nrot); 

					// get the estimate of curvatures from eigenvalues and eigenvectors
					// find the 2 most tangent eigenvectors (by finding the one closest to the normal)
					int best = 0; ScalarType bestv = fabs( (*vi).cN() * eigenvectors.GetColumn(0).Normalize());
					for(int i  = 1 ; i < 3; ++i){
						ScalarType prod = fabs((*vi).cN() * eigenvectors.GetColumn(i).Normalize());
						if( prod > bestv){bestv = prod; best = i;}
					}

					(*vi).PD1()  = eigenvectors.GetColumn( (best+1)%3).Normalize();
					(*vi).PD2()  = eigenvectors.GetColumn( (best+2)%3).Normalize();

					// project them to the plane identified by the normal
					vcg::Matrix33<ScalarType> rot;
					ScalarType angle = acos((*vi).PD1()*(*vi).N());
					rot.SetRotateRad(  - (M_PI*0.5 - angle),(*vi).PD1()^(*vi).N());
					(*vi).PD1() = rot*(*vi).PD1();
					angle = acos((*vi).PD2()*(*vi).N());
					rot.SetRotateRad(  - (M_PI*0.5 - angle),(*vi).PD2()^(*vi).N());
					(*vi).PD2() = rot*(*vi).PD2();


					// copmutes the curvature values
					const ScalarType r5 = r*r*r*r*r;
					const ScalarType r6 = r*r5;
					(*vi).K1() = (2.0/5.0) * (4.0*M_PI*r5 + 15*eigenvalues[(best+2)%3]-45.0*eigenvalues[(best+1)%3])/(M_PI*r6);
					(*vi).K2() = (2.0/5.0) * (4.0*M_PI*r5 + 15*eigenvalues[(best+1)%3]-45.0*eigenvalues[(best+2)%3])/(M_PI*r6);
					if((*vi).K1() < (*vi).K2())	{	std::swap((*vi).K1(),(*vi).K2());
																				std::swap((*vi).PD1(),(*vi).PD2());
																			}
			}
			

		}
 /// \brief Computes the discrete gaussian curvature. 
 
/** For further details, please, refer to: \n

- <em> Discrete Differential-Geometry Operators for Triangulated 2-Manifolds Mark Meyer,
 Mathieu Desbrun, Peter Schroder, Alan H. Barr VisMath '02, Berlin </em>
*/   
	static void MeanAndGaussian(MeshType & m) 
    {
      float area0, area1, area2, angle0, angle1, angle2, e01, e12, e20;
			FaceIterator fi;
      VertexIterator vi;   
			typename MeshType::CoordType  e01v ,e12v ,e20v;

			SimpleTempData<VertContainer, AreaData> TDAreaPtr(m.vert);  
			SimpleTempData<VertContainer, typename MeshType::CoordType> TDContr(m.vert);  

 			vcg::tri::UpdateNormals<MeshType>::PerVertexNormalized(m);
     //Compute AreaMix in H (vale anche per K)
      for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi) if(!(*vi).IsD())
      {
        (TDAreaPtr)[*vi].A = 0.0;
				(TDContr)[*vi]  =typename MeshType::CoordType(0.0,0.0,0.0);
				(*vi).Kh() = 0.0;
        (*vi).Kg() = (float)(2.0 * M_PI);
      }

      for(fi=m.face.begin();fi!=m.face.end();++fi) if( !(*fi).IsD())
      {
        // angles
        angle0 = math::Abs(Angle(	(*fi).P(1)-(*fi).P(0),(*fi).P(2)-(*fi).P(0) ));
        angle1 = math::Abs(Angle(	(*fi).P(0)-(*fi).P(1),(*fi).P(2)-(*fi).P(1) ));
        angle2 = M_PI-(angle0+angle1);
        
        if((angle0 < M_PI/2) && (angle1 < M_PI/2) && (angle2 < M_PI/2))  // triangolo non ottuso
        { 
	        float e01 = SquaredDistance( (*fi).V(1)->cP() , (*fi).V(0)->cP() );
	        float e12 = SquaredDistance( (*fi).V(2)->cP() , (*fi).V(1)->cP() );
	        float e20 = SquaredDistance( (*fi).V(0)->cP() , (*fi).V(2)->cP() );
      	
          area0 = ( e20*(1.0/tan(angle1)) + e01*(1.0/tan(angle2)) ) / 8.0;
	        area1 = ( e01*(1.0/tan(angle2)) + e12*(1.0/tan(angle0)) ) / 8.0;
	        area2 = ( e12*(1.0/tan(angle0)) + e20*(1.0/tan(angle1)) ) / 8.0;
      	
	        (TDAreaPtr)[(*fi).V(0)].A  += area0;
	        (TDAreaPtr)[(*fi).V(1)].A  += area1;
	        (TDAreaPtr)[(*fi).V(2)].A  += area2;

	      }
        else // obtuse
	      { 
					(TDAreaPtr)[(*fi).V(0)].A += vcg::DoubleArea<typename MeshType::FaceType>((*fi)) / 6.0;
	        (TDAreaPtr)[(*fi).V(1)].A += vcg::DoubleArea<typename MeshType::FaceType>((*fi)) / 6.0;
	        (TDAreaPtr)[(*fi).V(2)].A += vcg::DoubleArea<typename MeshType::FaceType>((*fi)) / 6.0;      
	      }
      }   
     
      for(fi=m.face.begin();fi!=m.face.end();++fi) if( !(*fi).IsD() )
      {    
        angle0 = math::Abs(Angle(	(*fi).P(1)-(*fi).P(0),(*fi).P(2)-(*fi).P(0) ));
        angle1 = math::Abs(Angle(	(*fi).P(0)-(*fi).P(1),(*fi).P(2)-(*fi).P(1) ));
        angle2 = M_PI-(angle0+angle1);
        
        e01v = ( (*fi).V(1)->cP() - (*fi).V(0)->cP() ) ;
        e12v = ( (*fi).V(2)->cP() - (*fi).V(1)->cP() ) ;
        e20v = ( (*fi).V(0)->cP() - (*fi).V(2)->cP() ) ;
        
        TDContr[(*fi).V(0)] += ( e20v * (1.0/tan(angle1)) - e01v * (1.0/tan(angle2)) ) / 4.0;
	      TDContr[(*fi).V(1)] += ( e01v * (1.0/tan(angle2)) - e12v * (1.0/tan(angle0)) ) / 4.0;
	      TDContr[(*fi).V(2)] += ( e12v * (1.0/tan(angle0)) - e20v * (1.0/tan(angle1)) ) / 4.0;
          
        (*fi).V(0)->Kg() -= angle0;
        (*fi).V(1)->Kg() -= angle1;
        (*fi).V(2)->Kg() -= angle2;

        
        for(int i=0;i<3;i++)
		    {
			    if(vcg::face::IsBorder((*fi), i))
			    {
				    CoordType e1,e2;
				    vcg::face::Pos<FaceType> hp(&*fi, i, (*fi).V(i));
				    vcg::face::Pos<FaceType> hp1=hp;
				    
            hp1.FlipV();
    	      e1=hp1.v->cP() - hp.v->cP();
				    hp1.FlipV();
				    hp1.NextB();
				    e2=hp1.v->cP() - hp.v->cP();
            (*fi).V(i)->Kg() -= math::Abs(Angle(e1,e2));
			    }
	      }
      }
         
      for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi) if(!(*vi).IsD() /*&& !(*vi).IsB()*/)
      {
        if((TDAreaPtr)[*vi].A<=std::numeric_limits<ScalarType>::epsilon())
        {
          (*vi).Kh() = 0;
          (*vi).Kg() = 0;
        }
        else
        {
					(*vi).Kh()  = (((TDContr)[*vi]* (*vi).cN()>0)?1.0:-1.0)*((TDContr)[*vi] / (TDAreaPtr) [*vi].A).Norm();
          (*vi).Kg() /= (TDAreaPtr)[*vi].A;
        }
			}
    }
	
	
	/// \brief Update the mean and the gaussian curvature of a vertex.
	
	/**
	The function uses the VF adiacency to walk around the vertex. 
	\return It will return the voronoi area around the vertex.  If (norm == true) the mean and the gaussian curvature are normalized.
	 Based on the paper  <a href="http://www2.in.tu-clausthal.de/~hormann/papers/Dyn.2001.OTU.pdf">  <em> "Optimizing 3d triangulations using discrete curvature analysis" </em> </a>
	  */
	  
	static float VertexCurvature(VertexPointer v, bool norm = true)
	{
		// VFAdjacency required!
		assert(FaceType::HasVFAdjacency());
		assert(VertexType::HasVFAdjacency());
		
		VFIteratorType vfi(v);
		float A = 0;
		
		v->Kh() = 0;
		v->Kg() = 2 * M_PI;

		while (!vfi.End()) {
			if (!vfi.F()->IsD()) {
				FacePointer f = vfi.F();
				int i = vfi.I();
				VertexPointer v0 = f->V0(i), v1 = f->V1(i), v2 = f->V2(i);
				
				float ang0 = math::Abs(Angle(v1->P() - v0->P(), v2->P() - v0->P() ));
				float ang1 = math::Abs(Angle(v0->P() - v1->P(), v2->P() - v1->P() ));
				float ang2 = M_PI - ang0 - ang1;

				float s01 = SquaredDistance(v1->P(), v0->P());
				float s02 = SquaredDistance(v2->P(), v0->P());

				// voronoi cell of current vertex
				if (ang0 >= M_PI/2)
					A += (0.5f * DoubleArea(*f) - (s01 * tan(ang1) + s02 * tan(ang2)) / 8.0 );
				else if (ang1 >= M_PI/2)
					A += (s01 * tan(ang0)) / 8.0;
				else if (ang2 >= M_PI/2)
					A += (s02 * tan(ang0)) / 8.0;
				else  // non obctuse triangle
					A += ((s02 / tan(ang1)) + (s01 / tan(ang2))) / 8.0;
				
				// gaussian curvature update
				v->Kg() -= ang0;

				// mean curvature update
				ang1 = math::Abs(Angle(f->N(), v1->N()));
				ang2 = math::Abs(Angle(f->N(), v2->N()));
				v->Kh() += ( (math::Sqrt(s01) / 2.0) * ang1 + 
				             (math::Sqrt(s02) / 2.0) * ang2 );
			}
			
			++vfi;
		}
		
		v->Kh() /= 4.0f;
		
		if(norm) {
			if(A <= std::numeric_limits<float>::epsilon()) {
				v->Kh() = 0;
				v->Kg() = 0;
			}
			else {
				v->Kh() /= A;
				v->Kg() /= A;
			}
		}

		return A;
	}

	static void VertexCurvature(MeshType & m){

		for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
			VertexCurvature(&*vi,false);
	}

};


}
}
#endif
