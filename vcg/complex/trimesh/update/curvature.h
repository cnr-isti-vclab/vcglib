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

#include <vcg/math/base.h>
#include <vcg/math/matrix.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/simplex/face/jumping_pos.h>
#include <vcg/container/simple_temporary_data.h>

namespace vcg {
namespace tri {

/** \addtogroup trimesh */
/*@{*/

/// Management, updating and computation of per-vertex and per-face normals.
/// This class is used to compute or update the normals that can be stored in the vertex or face component of a mesh.
template <class MeshType>
class UpdateCurvature
{

public:
		typedef typename MeshType::FaceIterator FaceIterator;
		typedef typename MeshType::VertexIterator VertexIterator;
		typedef typename MeshType::VertContainer VertContainer;
		typedef typename MeshType::VertexType VertexType;
		typedef typename MeshType::FaceType FaceType;
		typedef typename MeshType::CoordType CoordType;		
		typedef typename CoordType::ScalarType ScalarType;

	
private:
		typedef struct AdjVertex {
			VertexType * vert;
			float doubleArea;
			bool isBorder;
		};
public:
	/*
	Compute principal direction and magniuto of curvature as describe in the paper:
	@InProceedings{bb33922,
  author =	"G. Taubin",
  title =	"Estimating the Tensor of Curvature of a Surface from a
		 Polyhedral Approximation",
  booktitle =	"International Conference on Computer Vision",
  year = 	"1995",
  pages =	"902--907",
  URL =  	"http://dx.doi.org/10.1109/ICCV.1995.466840",
  bibsource =	"http://www.visionbib.com/bibliography/describe440.html#TT32253",
	}
	*/
		static void PrincipalDirections(MeshType &m) {

			assert(m.HasVFTopology());

			vcg::tri::UpdateNormals<MeshType>::PerVertexNormalized(m);

			VertexIterator vi;
			for (vi =m.vert.begin(); vi !=m.vert.end(); ++vi) {
				if ( ! (*vi).IsD() && (*vi).VFp() != NULL) {

					VertexType * central_vertex = &(*vi);

					std::vector<float> weights;
					std::vector<AdjVertex> vertices;

					vcg::face::JumpingPos<FaceType> pos((*vi).VFp(), central_vertex);

					VertexType* firstV = pos.VFlip();
					VertexType* tempV;
					float totalDoubleAreaSize = 0.0f;

					if (((firstV->P()-central_vertex->P())^(pos.VFlip()->P()-central_vertex->P()))*central_vertex->N()<=0.0f)
					{
						pos.Set(central_vertex->VFp(), central_vertex);
						pos.FlipE();
						firstV = pos.VFlip();
					}	
					else pos.Set(central_vertex->VFp(), central_vertex);

					do 
					{
						pos.NextE();
						tempV = pos.VFlip();

						AdjVertex v;

						v.isBorder = pos.IsBorder();
						v.vert = tempV;			
						v.doubleArea = ((pos.F()->V(1)->P() - pos.F()->V(0)->P()) ^ (pos.F()->V(2)->P()- pos.F()->V(0)->P())).Norm();;
						totalDoubleAreaSize += v.doubleArea;

						vertices.push_back(v);						
					} 
					while(tempV != firstV);	

					for (int i = 0; i < vertices.size(); ++i) {
						if (vertices[i].isBorder) {
							weights.push_back(vertices[i].doubleArea / totalDoubleAreaSize);
						} else {
							weights.push_back(0.5f * (vertices[i].doubleArea + vertices[(i-1)%vertices.size()].doubleArea) / totalDoubleAreaSize);
						}
						assert(weights.back() < 1.0f);
					}

					Matrix33f Tp;
					for (int i = 0; i < 3; ++i)
						Tp[i][i] = 1.0f - powf(central_vertex->N()[i],2);
					Tp[0][1] = Tp[1][0] = -1.0f * (central_vertex->N()[0] * central_vertex->N()[1]);
					Tp[1][2] = Tp[2][1] = -1.0f * (central_vertex->N()[1] * central_vertex->N()[2]);
					Tp[0][2] = Tp[2][0] = -1.0f * (central_vertex->N()[0] * central_vertex->N()[2]);

					Matrix33f tempMatrix;
					Matrix33f M;
					M.SetZero();
					for (int i = 0; i < vertices.size(); ++i) {
						Point3f edge = (central_vertex->P() - vertices[i].vert->P());
						float curvature = (2.0f * (central_vertex->N() * edge) ) / edge.SquaredNorm();
						Point3f T = (Tp*edge).Normalize();
						tempMatrix.ExternalProduct(T,T);
						M += tempMatrix * weights[i] * curvature ;
					}

					Point3f W;
					Point3f e1(1.0f,0.0f,0.0f);
					if ((e1 - central_vertex->N()).SquaredNorm() > (e1 + central_vertex->N()).SquaredNorm())
						W = e1 - central_vertex->N();
					else 
						W = e1 + central_vertex->N();
					W.Normalize();

					Matrix33f Q;
					Q.SetIdentity();
					tempMatrix.ExternalProduct(W,W);
					Q -= tempMatrix * 2.0f;

					Matrix33f Qt(Q);
					Qt.Transpose();

					Matrix33f QtMQ = (Qt * M * Q);

					Point3f bl = Q.GetColumn(0);
					Point3f T1 = Q.GetColumn(1);
					Point3f T2 = Q.GetColumn(2);

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
					float min_error = std::numeric_limits<float>::infinity();
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

					vcg::ndim::MatrixMNf minor2x2 (2,2);
					vcg::ndim::MatrixMNf S (2,2);


					minor2x2[0][0] = QtMQ[1][1];
					minor2x2[0][1] = QtMQ[1][2];
					minor2x2[1][0] = QtMQ[2][1];
					minor2x2[1][1] = QtMQ[2][2];

					S[0][0] = S[1][1] = c;
					S[0][1] = s;
					S[1][0] = -1.0f * s;

					vcg::ndim::MatrixMNf St (S);
					St.Transpose();					

					vcg::ndim::MatrixMNf StMS(St * minor2x2 * S);

					float Principal_Curvature1 = (3.0f * StMS[0][0]) - StMS[1][1];
					float Principal_Curvature2 = (3.0f * StMS[1][1]) - StMS[0][0];

					Point3f Principal_Direction1 = T1 * c - T2 * s;
					Point3f Principal_Direction2 = T1 * s + T2 * c; 

					(*vi).PD1() = Principal_Direction1 ;
					(*vi).PD2() = Principal_Direction2 ;
					(*vi).K1() = -Principal_Curvature1;
					(*vi).K2() = -Principal_Curvature2;

				}
			}
		}
	 


  class AreaData
  {
  public:
    float A;
  };


 /** computes the discrete gaussian curvature as proposed in 
Discrete Differential-Geometry Operators for Triangulated 2-Manifolds Mark Meyer,
 Mathieu Desbrun, Peter Schroder, Alan H. Barr VisMath '02, Berlin
*/   
	static void MeanAndGaussian(MeshType & m) 
    {
      float area0, area1, area2, angle0, angle1, angle2, e01, e12, e20;
			FaceIterator fi;
      VertexIterator vi;   

			SimpleTempData<VertContainer, AreaData> TDAreaPtr(m.vert); TDAreaPtr.Start();

      //Calcola AreaMix in H (vale anche per K)
      for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi) if(!(*vi).IsD())
      {
        (TDAreaPtr)[*vi].A = 0;
        (*vi).H() = 0;
        (*vi).K() = (float)(2.0 * M_PI);
      }

      for(fi=m.face.begin();fi!=m.face.end();++fi) if( !(*fi).IsD())
      {
        // angles
        angle0 = math::Abs(Angle(	(*fi).P(1)-(*fi).P(0),(*fi).P(2)-(*fi).P(0) ));
        angle1 = math::Abs(Angle(	(*fi).P(0)-(*fi).P(1),(*fi).P(2)-(*fi).P(1) ));
        angle2 = M_PI-(angle0+angle1);
        
        if((angle0 < M_PI/2) && (angle1 < M_PI/2) && (angle2 < M_PI/2))  // triangolo non ottuso
        { 
	        float e01 = SquaredDistance( (*fi).V(1)->P() , (*fi).V(0)->P() );
	        float e12 = SquaredDistance( (*fi).V(2)->P() , (*fi).V(1)->P() );
	        float e20 = SquaredDistance( (*fi).V(0)->P() , (*fi).V(2)->P() );
      	
          area0 = ( e20*(1.0/tan(angle1)) + e01*(1.0/tan(angle2)) ) / 8.0;
	        area1 = ( e01*(1.0/tan(angle2)) + e12*(1.0/tan(angle0)) ) / 8.0;
	        area2 = ( e12*(1.0/tan(angle0)) + e20*(1.0/tan(angle1)) ) / 8.0;
      	
	        (TDAreaPtr)[(*fi).V(0)].A  += area0;
	        (TDAreaPtr)[(*fi).V(1)].A  += area1;
	        (TDAreaPtr)[(*fi).V(2)].A  += area2;

	      }
        else // triangolo ottuso
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
        
        e01 = ( (*fi).V(1)->P() - (*fi).V(0)->P() ) * (*fi).V(0)->N();
        e12 = ( (*fi).V(2)->P() - (*fi).V(1)->P() ) * (*fi).V(1)->N();
        e20 = ( (*fi).V(0)->P() - (*fi).V(2)->P() ) * (*fi).V(2)->N();
        
        area0 = ( e20 * (1.0/tan(angle1)) + e01 * (1.0/tan(angle2)) ) / 4.0;
	      area1 = ( e01 * (1.0/tan(angle2)) + e12 * (1.0/tan(angle0)) ) / 4.0;
	      area2 = ( e12 * (1.0/tan(angle0)) + e20 * (1.0/tan(angle1)) ) / 4.0;
          
        (*fi).V(0)->H()  += area0;
	      (*fi).V(1)->H()  += area1;
	      (*fi).V(2)->H()  += area2;

        (*fi).V(0)->K() -= angle0;
        (*fi).V(1)->K() -= angle1;
        (*fi).V(2)->K() -= angle2;

        
        for(int i=0;i<3;i++)
		    {
			    if(vcg::face::IsBorder((*fi), i))
			    {
				    CoordType e1,e2;
				    vcg::face::Pos<FaceType> hp(&*fi, i, (*fi).V(i));
				    vcg::face::Pos<FaceType> hp1=hp;
				    
            hp1.FlipV();
    	      e1=hp1.v->P() - hp.v->P();
				    hp1.FlipV();
				    hp1.NextB();
				    e2=hp1.v->P() - hp.v->P();
            (*fi).V(i)->K() -= math::Abs(Angle(e1,e2));
			    }
	      }
      }
         
      for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi) if(!(*vi).IsD() /*&& !(*vi).IsB()*/)
      {
        if((TDAreaPtr)[*vi].A<=std::numeric_limits<float>::epsilon())
        {
          (*vi).H() = 0;
          (*vi).K() = 0;
        }
        else
        {
          (*vi).H() /= (TDAreaPtr)[*vi].A;
          (*vi).K() /= (TDAreaPtr)[*vi].A;
        }
			}

			TDAreaPtr.Stop();

    }
    
};


}
}
#endif
