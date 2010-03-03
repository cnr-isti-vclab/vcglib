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
Revision 1.10  2006/01/22 10:06:23  cignoni
Corrected use of Area with the unambiguous DoubleArea

Revision 1.9  2005/09/28 19:35:06  m_di_benedetto
Added class PointDistanceFunctor.

Revision 1.8  2005/09/14 12:58:44  pietroni
changed min calls to Min<ScalarType> of math.h of vcglib

Revision 1.7  2005/09/14 09:58:32  pietroni
removed vcg::math::Min<ScalarType> definition generate warnings

Revision 1.6  2005/09/14 09:03:54  pietroni
added definition of vcg::math::Min<ScalarType> function

Revision 1.5  2005/02/02 16:44:34  pietroni
1 warning corrected added casting in const ScalarType EPSILON = ScalarType( 0.000001);

Revision 1.4  2005/01/28 12:00:33  cignoni
small gcc compiling issues for namespaces

Revision 1.3  2005/01/24 15:35:25  cignoni
Removed a 'using namespace'

Revision 1.2  2005/01/21 17:11:03  pietroni
changed Dist Function to PointDistance... the function is on vcg::face::PointDistance this file will contain all distance functions between a face and othe entities

Revision 1.1  2004/05/12 18:50:25  ganovelli
created


****************************************************************************/

#ifndef __VCGLIB_FACE_DISTANCE
#define __VCGLIB_FACE_DISTANCE

#include <vcg/math/base.h>
#include <vcg/space/point3.h>
#include <vcg/space/segment3.h>


namespace vcg {
	namespace face{
/*
   Point face distance
   trova il punto <p> sulla faccia piu' vicino a <q>, con possibilit� di 
   rejection veloce su se la distanza trovata � maggiore di <rejdist>

 Commenti del 12/11/02
 Funziona solo se la faccia e di quelle di tipo E (con edge e piano per faccia gia' calcolati)
 algoritmo:
	1) si calcola la proiezione <p> di q sul piano della faccia
	2) se la distanza punto piano e' > rejdist ritorna
	3) si lavora sul piano migliore e si cerca di capire se il punto sta dentro il triangolo:
	   a) prodotto vettore tra edge triangolo (v[i+1]-v[i]) e (p-v[i])
		 b) se il risultato e' negativo (gira in senso orario) allora il punto
		    sta fuori da quella parte e si fa la distanza punto segmento.
     c) se il risultato sempre positivo allora sta dentro il triangolo
	4) e si restituisce la distanza punto /piano gia` calcolata 

	Note sulla robustezza:
	il calcolo del prodotto vettore e` la cosa piu` delicata:
	possibili fallimenti quando a^b ~= 0
	1) doveva essere <= 0 e viene positivo (q era fuori o sulla linea dell'edge)
	   allora capita che si faccia la distanza punto piano anziche` la distanza punto seg
  2) doveva essere > 0 e viene <=0 (q era dentro il triangolo)

*/
	template <class FaceType>
	bool PointDistance(	const FaceType &f, 
							const vcg::Point3<typename FaceType::ScalarType> & q, 
							typename FaceType::ScalarType & dist, 
							vcg::Point3<typename FaceType::ScalarType> & p )
	{
		typedef typename FaceType::ScalarType ScalarType;
		
                const ScalarType EPS = ScalarType( 0.000001);

                //const ScalarType EPSILON = 0.00000001;
		ScalarType b,b0,b1,b2;
			// Calcolo distanza punto piano
		ScalarType d = Distance( f.cPlane(), q );
		if( d>dist || d<-dist )			// Risultato peggiore: niente di fatto
			return false;

			// Calcolo del punto sul piano
		// NOTA: aggiunto un '-d' in fondo Paolo C.
		Point3<ScalarType> t = f.cPlane().Direction();
		t[0] *= -d;
		t[1] *= -d;
		t[2] *= -d;
		p = q; p += t;
		    
		switch( f.Flags() & (FaceType::NORMX|FaceType::NORMY|FaceType::NORMZ) )
		{
		case FaceType::NORMX:
			b0 = f.cEdge(1)[1]*(p[2] - f.cP(1)[2]) - f.cEdge(1)[2]*(p[1] - f.cP(1)[1]);
			if(b0<=0)
			{
				b0 = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
				if(dist>b0) { dist = b0; return true; }
				else return false;
			}
			b1 = f.cEdge(2)[1]*(p[2] - f.cP(2)[2]) - f.cEdge(2)[2]*(p[1] - f.cP(2)[1]);
			if(b1<=0)
			{
				b1 = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
				if(dist>b1) { dist = b1; return true; }
				else return false;
			}
			b2 = f.cEdge(0)[1]*(p[2] - f.cP(0)[2]) - f.cEdge(0)[2]*(p[1] - f.cP(0)[1]);
			if(b2<=0)
			{
				b2 = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
				if(dist>b2) { dist = b2; return true; }
				else return false;
			}
			// sono tutti e tre > 0 quindi dovrebbe essere dentro;
			// per sicurezza se il piu' piccolo dei tre e' < epsilon (scalato rispetto all'area della faccia
			// per renderlo dimension independent.) allora si usa ancora la distanza punto 
			// segmento che e' piu robusta della punto piano, e si fa dalla parte a cui siamo piu' 
			// vicini (come prodotto vettore)
			// Nota: si potrebbe rendere un pochino piu' veloce sostituendo Area()
			// con il prodotto vettore dei due edge in 2d lungo il piano migliore.
                        if( (b=vcg::math::Min<ScalarType>(b0,vcg::math::Min<ScalarType>(b1,b2))) < EPS*DoubleArea(f))
      {
				ScalarType bt;
				if(b==b0) 	    bt = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
				else if(b==b1) 	bt = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
				else if(b==b2) 	bt = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
                                //printf("Warning area:%g %g %g %g thr:%g bt:%g\n",Area(), b0,b1,b2,EPS*Area(),bt);
				if(dist>bt) { dist = bt; return true; }
				else return false;
			}
			break;

		  case FaceType::NORMY:
			b0 = f.cEdge(1)[2]*(p[0] - f.cP(1)[0]) - f.cEdge(1)[0]*(p[2] - f.cP(1)[2]);
			if(b0<=0)
			{
				b0 = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
				if(dist>b0) { dist = b0; return true; }
				else return false;
			}
			b1 = f.cEdge(2)[2]*(p[0] - f.cP(2)[0]) - f.cEdge(2)[0]*(p[2] - f.cP(2)[2]);
			if(b1<=0)
			{
				b1 = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
				if(dist>b1) { dist = b1; return true; }
				else return false;
			}
			b2 = f.cEdge(0)[2]*(p[0] - f.cP(0)[0]) - f.cEdge(0)[0]*(p[2] - f.cP(0)[2]);
			if(b2<=0)
			{
				b2 = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
				if(dist>b2) { dist = b2; return true; }
				else return false;
			}
                        if( (b=vcg::math::Min<ScalarType>(b0,vcg::math::Min<ScalarType>(b1,b2))) < EPS*DoubleArea(f))
      {
				ScalarType bt;
				if(b==b0) 	    bt = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
				else if(b==b1) 	bt = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
				else if(b==b2) 	bt = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
				//printf("Warning area:%g %g %g %g thr:%g bt:%g\n",Area(), b0,b1,b2,EPSILON*Area(),bt);
				if(dist>bt) { dist = bt; return true; }
				else return false;
			}
			break;

		  case FaceType::NORMZ:
			b0 = f.cEdge(1)[0]*(p[1] - f.cP(1)[1]) - f.cEdge(1)[1]*(p[0] - f.cP(1)[0]);
			if(b0<=0)
			{
				b0 = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
				if(dist>b0) { dist = b0; return true; }
				else return false;
			}
			b1 = f.cEdge(2)[0]*(p[1] - f.cP(2)[1]) - f.cEdge(2)[1]*(p[0] - f.cP(2)[0]);
			if(b1<=0)
			{
				b1 = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
				if(dist>b1) { dist = b1; return true; }
				else return false;
			}
			b2 = f.cEdge(0)[0]*(p[1] - f.cP(0)[1]) - f.cEdge(0)[1]*(p[0] - f.cP(0)[0]);
			if(b2<=0)
			{
				b2 = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
				if(dist>b2) { dist = b2; return true; }
				else return false;
			}
                        if( (b=vcg::math::Min<ScalarType>(b0,vcg::math::Min<ScalarType>(b1,b2))) < EPS*DoubleArea(f))
      {
				ScalarType bt;
				if(b==b0) 	    bt = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
				else if(b==b1) 	bt = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
				else if(b==b2) 	bt = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
				//printf("Warning area:%g %g %g %g thr:%g bt:%g\n",Area(), b0,b1,b2,EPSILON*Area(),bt);
				
				if(dist>bt) { dist = bt; return true; }
				else return false;
			}
			break;

		}

		dist = ScalarType(fabs(d));
		//dist = Distance(p,q);
		return true;
	}

	template <class S>
	class PointDistanceFunctor {
	public:
		typedef S ScalarType;
		typedef Point3<ScalarType> QueryType;
		static inline const Point3<ScalarType> &  Pos(const QueryType & qt)  {return qt;}

		template <class FACETYPE, class SCALARTYPE>
		inline bool operator () (const FACETYPE & f, const Point3<SCALARTYPE> & p, SCALARTYPE & minDist, Point3<SCALARTYPE> & q) {
			const Point3<typename FACETYPE::ScalarType> fp = Point3<typename FACETYPE::ScalarType>::Construct(p);
			Point3<typename FACETYPE::ScalarType> fq;
			typename FACETYPE::ScalarType md = (typename FACETYPE::ScalarType)(minDist);
			const bool ret = PointDistance(f, fp, md, fq);
			minDist = (SCALARTYPE)(md);
			q = Point3<SCALARTYPE>::Construct(fq);
			return (ret);
		}
	};



	template <class S>
	class PointNormalDistanceFunctor {
	public:
		typedef typename S::ScalarType ScalarType;
		typedef S QueryType;
		static inline const Point3<ScalarType> &  Pos(const QueryType & qt)  {return qt.P();}


		static ScalarType & Alpha(){static ScalarType alpha = 1.0; return alpha;}
		static ScalarType & Beta (){static ScalarType beta  = 1.0; return beta;}
		static ScalarType & Gamma(){static ScalarType gamma = 1.0; return gamma;}
		static ScalarType & InterPoint(){static ScalarType interpoint= 1.0; return interpoint;} 


		template <class FACETYPE, class SCALARTYPE>
		inline bool operator () (const FACETYPE &f, const typename FACETYPE::VertexType &p, 
			SCALARTYPE & minDist,Point3<SCALARTYPE> & q)
		{		
			const Point3<typename FACETYPE::ScalarType> fp = Point3<typename FACETYPE::ScalarType>::Construct(p.cP());
			const Point3<typename FACETYPE::ScalarType> fn = Point3<typename FACETYPE::ScalarType>::Construct(p.cN());
			Point3<typename FACETYPE::ScalarType> fq;
			typename FACETYPE::ScalarType md = (typename FACETYPE::ScalarType)(minDist);
			const bool ret=PointDistance(f,fp,md,fq);

			SCALARTYPE  dev=InterPoint()*(pow((ScalarType)(1-f.cN().dot(fn)),(ScalarType)Beta())/(Gamma()*md+0.1));

			if (md+dev < minDist){
				minDist = (SCALARTYPE)(md+dev);
				q = Point3<SCALARTYPE>::Construct(fq);
				//q.N() = f.N();
				return (ret);
			}
			return false;
		}
	};
		
		/// BASIC VERSION of the Point-face distance that does not require the EdgePlane Additional data.
		/// Given a face and a point, returns the closest point of the face to p.
		/// it assumes that the face has Normalized Normal and on the flags stored the preferred orientation.
		// UpdateNormals::PerFaceNormalized(m)
		// UpdateFlags<>::FaceProjection(m);
		
		template <class FaceType>
			bool PointDistanceBase(
													const FaceType &f,																		/// the face to be tested
													const vcg::Point3<typename FaceType::ScalarType> & q, /// the point tested
													typename FaceType::ScalarType & dist,                 /// bailout distance. It must be initialized with the max admittable value. 
													vcg::Point3<typename FaceType::ScalarType> & p )      
		{
				typedef typename FaceType::ScalarType ScalarType;
				// remember that the macro NDEBUG is defined when you want to optimize a lot. 
				#ifndef NDEBUG
				static int staticCnt=0; // small piece of code that sometime check that face normals are really normalized
				if((staticCnt++%100)==0) 
            assert((f.cN().SquaredNorm() ==0) || (f.cN().SquaredNorm() > 0.9999 && f.cN().SquaredNorm()<1.0001)); // if you get this assert you have forgot to make a UpdateNormals::PerFaceNormalized(m)
        #endif

        if(f.cN()==Point3<ScalarType>(0,0,0)) // to correctly manage the case of degenerate triangles we consider them as segments.
        {
            Box3<ScalarType> bb;
            f.GetBBox(bb);
            Segment3<ScalarType> degenTri(bb.min,bb.max);
            Point3<ScalarType> closest= ClosestPoint( degenTri, q );
            ScalarType d = Distance(closest, q);
            if( d>dist || d<-dist )			// Risultato peggiore: niente di fatto
                      return false;
            dist=d;
            p=closest;
            return true;
        }

				Plane3<ScalarType> fPlane;
				fPlane.Init(f.cP(0),f.cN());
        const ScalarType EPS = ScalarType( 0.000001);
				ScalarType b,b0,b1,b2;
				// Calcolo distanza punto piano
				ScalarType d = Distance( fPlane, q );
				if( d>dist || d<-dist )			// Risultato peggiore: niente di fatto
					return false;
				
				// Calcolo del punto sul piano
				// NOTA: aggiunto un '-d' in fondo Paolo C.
				Point3<ScalarType> t = fPlane.Direction();
				t[0] *= -d;
				t[1] *= -d;
				t[2] *= -d;
				p = q; p += t;

				Point3<ScalarType> fEdge[3];				
				fEdge[0] = f.cP(1); fEdge[0] -= f.cP(0);
				fEdge[1] = f.cP(2); fEdge[1] -= f.cP(1);
				fEdge[2] = f.cP(0); fEdge[2] -= f.cP(2);
				
				/* 
				This piece of code is part of the EdgePlane initialization structure: note that the edges are scaled!. 
				 
				if(nx>ny && nx>nz) { f.Flags() |= FaceType::NORMX; d = 1/f.Plane().Direction()[0]; }
				else if(ny>nz)     { f.Flags() |= FaceType::NORMY; d = 1/f.Plane().Direction()[1]; }
				else               { f.Flags() |= FaceType::NORMZ; d = 1/f.Plane().Direction()[2]; }
				f.Edge(0)*=d; f.Edge(1)*=d;f.Edge(2)*=d;
				 
				So we must apply the same scaling according to the plane orientation, eg in the case of NORMX
				
				scaleFactor= 1/fPlane.Direction()[0];
				fEdge[0]*=d; fEdge[1]*=d;fEdge[2]*=d;
				*/
				
				ScalarType scaleFactor;

				switch( f.Flags() & (FaceType::NORMX|FaceType::NORMY|FaceType::NORMZ) )
				{
					case FaceType::NORMX:
						scaleFactor= 1/fPlane.Direction()[0];
						fEdge[0]*=scaleFactor; fEdge[1]*=scaleFactor; fEdge[2]*=scaleFactor;
						
						b0 = fEdge[1][1]*(p[2] - f.cP(1)[2]) - fEdge[1][2]*(p[1] - f.cP(1)[1]);
						if(b0<=0)
						{
							b0 = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
							if(dist>b0) { dist = b0; return true; }
							else return false;
						}
							b1 = fEdge[2][1]*(p[2] - f.cP(2)[2]) - fEdge[2][2]*(p[1] - f.cP(2)[1]);
						if(b1<=0)
						{
							b1 = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
							if(dist>b1) { dist = b1; return true; }
							else return false;
						}
							b2 = fEdge[0][1]*(p[2] - f.cP(0)[2]) - fEdge[0][2]*(p[1] - f.cP(0)[1]);
						if(b2<=0)
						{
							b2 = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
							if(dist>b2) { dist = b2; return true; }
							else return false;
						}
							// sono tutti e tre > 0 quindi dovrebbe essere dentro;
							// per sicurezza se il piu' piccolo dei tre e' < epsilon (scalato rispetto all'area della faccia
							// per renderlo dimension independent.) allora si usa ancora la distanza punto 
							// segmento che e' piu robusta della punto piano, e si fa dalla parte a cui siamo piu' 
							// vicini (come prodotto vettore)
							// Nota: si potrebbe rendere un pochino piu' veloce sostituendo Area()
							// con il prodotto vettore dei due edge in 2d lungo il piano migliore.
              if( (b=vcg::math::Min<ScalarType>(b0,vcg::math::Min<ScalarType>(b1,b2))) < EPS*DoubleArea(f))
							{
								ScalarType bt;
								if(b==b0) 	    bt = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
								else if(b==b1) 	bt = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
								else if(b==b2) 	bt = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
								//printf("Warning area:%g %g %g %g thr:%g bt:%g\n",Area(), b0,b1,b2,EPSILON*Area(),bt);
								if(dist>bt) { dist = bt; return true; }
								else return false;
							}
							break;
						
					case FaceType::NORMY:
						scaleFactor= 1/fPlane.Direction()[1];
						fEdge[0]*=scaleFactor; fEdge[1]*=scaleFactor; fEdge[2]*=scaleFactor;
						
						b0 = fEdge[1][2]*(p[0] - f.cP(1)[0]) - fEdge[1][0]*(p[2] - f.cP(1)[2]);
						if(b0<=0)
						{
							b0 = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
							if(dist>b0) { dist = b0; return true; }
							else return false;
						}
							b1 = fEdge[2][2]*(p[0] - f.cP(2)[0]) - fEdge[2][0]*(p[2] - f.cP(2)[2]);
						if(b1<=0)
						{
							b1 = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
							if(dist>b1) { dist = b1; return true; }
							else return false;
						}
							b2 = fEdge[0][2]*(p[0] - f.cP(0)[0]) - fEdge[0][0]*(p[2] - f.cP(0)[2]);
						if(b2<=0)
						{
							b2 = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
							if(dist>b2) { dist = b2; return true; }
							else return false;
						}
            if( (b=vcg::math::Min<ScalarType>(b0,vcg::math::Min<ScalarType>(b1,b2))) < EPS*DoubleArea(f))
							{
								ScalarType bt;
								if(b==b0) 	    bt = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
								else if(b==b1) 	bt = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
								else if(b==b2) 	bt = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
								//printf("Warning area:%g %g %g %g thr:%g bt:%g\n",Area(), b0,b1,b2,EPSILON*Area(),bt);
								if(dist>bt) { dist = bt; return true; }
								else return false;
							}
							break;
						
					case FaceType::NORMZ:
						scaleFactor= 1/fPlane.Direction()[2];
						fEdge[0]*=scaleFactor; fEdge[1]*=scaleFactor; fEdge[2]*=scaleFactor;
						
						b0 = fEdge[1][0]*(p[1] - f.cP(1)[1]) - fEdge[1][1]*(p[0] - f.cP(1)[0]);
						if(b0<=0)
						{
							b0 = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
							if(dist>b0) { dist = b0; return true; }
							else return false;
						}
							b1 = fEdge[2][0]*(p[1] - f.cP(2)[1]) - fEdge[2][1]*(p[0] - f.cP(2)[0]);
						if(b1<=0)
						{
							b1 = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
							if(dist>b1) { dist = b1; return true; }
							else return false;
						}
							b2 = fEdge[0][0]*(p[1] - f.cP(0)[1]) - fEdge[0][1]*(p[0] - f.cP(0)[0]);
						if(b2<=0)
						{
							b2 = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
							if(dist>b2) { dist = b2; return true; }
							else return false;
						}
            if( (b=vcg::math::Min<ScalarType>(b0,vcg::math::Min<ScalarType>(b1,b2))) < EPS*DoubleArea(f))
							{
								ScalarType bt;
								if(b==b0) 	    bt = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
								else if(b==b1) 	bt = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
								else if(b==b2) 	bt = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
								//printf("Warning area:%g %g %g %g thr:%g bt:%g\n",Area(), b0,b1,b2,EPSILON*Area(),bt);
								
								if(dist>bt) { dist = bt; return true; }
								else return false;
							}
							break;
						default: assert(0); // if you get this assert it means that you forgot to set the required UpdateFlags<MeshType>::FaceProjection(m);

				}
				
				dist = ScalarType(fabs(d));
				//dist = Distance(p,q);
				return true;
		}
		
	template <class S>
	class PointDistanceBaseFunctor {
public:
			typedef S ScalarType;
			typedef Point3<ScalarType> QueryType;

 		  static inline const Point3<ScalarType> & Pos(const Point3<ScalarType> & qt)  {return qt;}
			template <class FACETYPE, class SCALARTYPE>
			inline bool operator () (const FACETYPE & f, const Point3<SCALARTYPE> & p, SCALARTYPE & minDist, Point3<SCALARTYPE> & q) {
				const Point3<typename FACETYPE::ScalarType> fp = Point3<typename FACETYPE::ScalarType>::Construct(p);
				Point3<typename FACETYPE::ScalarType> fq;
				typename FACETYPE::ScalarType md = (typename FACETYPE::ScalarType)(minDist);
				const bool ret = PointDistanceBase(f, fp, md, fq);
				minDist = (SCALARTYPE)(md);
				q = Point3<SCALARTYPE>::Construct(fq);
				return (ret);
			}
		};
		
		
		
}	 // end namespace face
	
}	 // end namespace vcg


#endif

