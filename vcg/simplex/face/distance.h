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
Revision 1.2  2005/01/21 17:11:03  pietroni
changed Dist Function to PointDistance... the function is on vcg::face::PointDistance this file will contain all distance functions between a face and othe entities

Revision 1.1  2004/05/12 18:50:25  ganovelli
created


****************************************************************************/

#ifndef __VCGLIB_FACE_DISTANCE
#define __VCGLIB_FACE_DISTANCE


#include <vcg/space/point3.h>


namespace vcg {
	namespace face{
/*
   Point face distance
   trova il punto <p> sulla faccia piu' vicino a <q>, con possibilità di 
   rejection veloce su se la distanza trovata è maggiore di <rejdist>

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
							const Point3<typename FaceType::ScalarType> & q, 
							typename FaceType::ScalarType & dist, 
							Point3<typename FaceType::ScalarType> & p )
	{
		typedef typename FaceType::ScalarType ScalarType;

		//const ScalarType EPSILON = ScalarType( 0.000001);
		const ScalarType EPSILON = 0.00000001;
		ScalarType b,b0,b1,b2;
			// Calcolo distanza punto piano
		ScalarType d = Distance( f.plane, q );
		if( d>dist || d<-dist )			// Risultato peggiore: niente di fatto
			return false;

			// Calcolo del punto sul piano
		// NOTA: aggiunto un '-d' in fondo Paolo C.
		Point3<ScalarType> t = f.plane.Direction();
		t[0] *= -d;
		t[1] *= -d;
		t[2] *= -d;
		p = q; p += t;
		    
		switch( f.Flags() & (FaceType::NORMX|FaceType::NORMY|FaceType::NORMZ) )
		{
		case FaceType::NORMX:
			b0 = f.edge[1][1]*(p[2] - f.cP(1)[2]) - f.edge[1][2]*(p[1] - f.cP(1)[1]);
			if(b0<=0)
			{
				b0 = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
				if(dist>b0) { dist = b0; return true; }
				else return false;
			}
			b1 = f.edge[2][1]*(p[2] - f.cP(2)[2]) - f.edge[2][2]*(p[1] - f.cP(2)[1]);
			if(b1<=0)
			{
				b1 = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
				if(dist>b1) { dist = b1; return true; }
				else return false;
			}
			b2 = f.edge[0][1]*(p[2] - f.cP(0)[2]) - f.edge[0][2]*(p[1] - f.cP(0)[1]);
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
			if( (b=min(b0,min(b1,b2))) < EPSILON*Area(f)) 
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
			b0 = f.edge[1][2]*(p[0] - f.cP(1)[0]) - f.edge[1][0]*(p[2] - f.cP(1)[2]);
			if(b0<=0)
			{
				b0 = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
				if(dist>b0) { dist = b0; return true; }
				else return false;
			}
			b1 = f.edge[2][2]*(p[0] - f.cP(2)[0]) - f.edge[2][0]*(p[2] - f.cP(2)[2]);
			if(b1<=0)
			{
				b1 = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
				if(dist>b1) { dist = b1; return true; }
				else return false;
			}
			b2 = f.edge[0][2]*(p[0] - f.cP(0)[0]) - f.edge[0][0]*(p[2] - f.cP(0)[2]);
			if(b2<=0)
			{
				b2 = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
				if(dist>b2) { dist = b2; return true; }
				else return false;
			}
			if( (b=min(b0,min(b1,b2))) < EPSILON*Area(f)) 
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
			b0 = f.edge[1][0]*(p[1] - f.cP(1)[1]) - f.edge[1][1]*(p[0] - f.cP(1)[0]);
			if(b0<=0)
			{
				b0 = PSDist(q,f.V(1)->cP(),f.V(2)->cP(),p);
				if(dist>b0) { dist = b0; return true; }
				else return false;
			}
			b1 = f.edge[2][0]*(p[1] - f.cP(2)[1]) - f.edge[2][1]*(p[0] - f.cP(2)[0]);
			if(b1<=0)
			{
				b1 = PSDist(q,f.V(2)->cP(),f.V(0)->cP(),p);
				if(dist>b1) { dist = b1; return true; }
				else return false;
			}
			b2 = f.edge[0][0]*(p[1] - f.cP(0)[1]) - f.edge[0][1]*(p[0] - f.cP(0)[0]);
			if(b2<=0)
			{
				b2 = PSDist(q,f.V(0)->cP(),f.V(1)->cP(),p);
				if(dist>b2) { dist = b2; return true; }
				else return false;
			}
			if( (b=min(b0,min(b1,b2))) < EPSILON*Area(f)) 
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
	
}	 // end namespace face
	
}	 // end namespace vcg


#endif

