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
Revision 1.14  2006/11/07 15:13:56  zifnab1974
Necessary changes for compilation with gcc 3.4.6. Especially the hash function is a problem

Revision 1.13  2006/11/07 11:47:11  cignoni
gcc compiling issues

Revision 1.12  2006/11/07 07:56:43  cignoni
Added missing std::

Revision 1.11  2006/11/06 16:12:29  giec
Leipa ear now compute max dihedral angle.

Revision 1.10  2006/10/31 11:30:41  ganovelli
changed access throught iterator with static call to comply 2005 compiler

Revision 1.9  2006/10/20 07:44:45  cignoni
Added missing std::

Revision 1.8  2006/10/18 15:06:47  giec
New policy for compute quality in TrivialEar.
Bugfixed LeipaEar.
Added new algorithm "selfintersection" with test for self intersection.

Revision 1.7  2006/10/10 09:12:02  giec
Bugfix and added a new type of ear (Liepa like)

Revision 1.6  2006/10/09 10:07:07  giec
Optimized version of "EAR HOLE FILLING", the Ear is selected according to its dihedral angle.

Revision 1.5  2006/10/06 15:28:14  giec
first working implementationof "EAR HOLE FILLING".

Revision 1.4  2006/10/02 12:06:40  giec
BugFix

Revision 1.3  2006/09/27 15:33:32  giec
It close one simple hole . . .

Revision 1.2  2006/09/27 09:29:53  giec
Frist working release whit a few bugs.
It almost fills the hole ...

Revision 1.1  2006/09/25 09:17:44  cignoni
First Non working Version

****************************************************************************/
#ifndef __VCG_TRI_UPDATE_HOLE
#define __VCG_TRI_UPDATE_HOLE

#include <vcg/math/base.h>

/*
Questa Classe serve per gestire la non duplicazione degli edge durante la chiusura 
di un buco.
*/
namespace vcg {
	namespace tri {

		template<class MESH>
		class SimpleEdge
		{
		public:
			typename MESH::VertexType v[2];
			SimpleEdge()
			{}

			SimpleEdge(typename MESH::VertexType v0, typename MESH::VertexType v1)
			{
				if(v0.P().X() != v1.P().X() && 
					v0.P().Y() != v1.P().Y() &&
					v0.P().Z() != v1.P().Z())
				{v[0]=v1; v[1]=v0;}
				else {v[0]=v0; v[1]=v1;}
			}

			SimpleEdge(typename face::Pos<typename MESH::FaceType> &ep)		{
				//*this=SimpleEdge(*ep.VFlip(), *ep.v);
				typename MESH::VertexType v0 ,v1;
				v0 = *ep.VFlip(); 
				v1 = *ep.v;
				if(v0.P().X() != v1.P().X() && 
					v0.P().Y() != v1.P().Y() &&
					v0.P().Z() != v1.P().Z())
				{v[0]=v1; v[1]=v0;}
				else {v[0]=v0; v[1]=v1;}
			}


			bool operator < (const SimpleEdge & e) const
			{		
                 v[0] = e.v[0]; v[1]=e.v[1];
			}

			bool operator != (const SimpleEdge & e)
			{		
				if(v[0].P().X() != e.v[0].P().X() && 
					v[0].P().Y() != e.v[0].P().Y() &&
					v[0].P().Z() != e.v[0].P().Z())
					return true;
				else return false;

			}
		};

		template<class MESH>
		class HoleInfo
		{
		public: 
			HoleInfo(){}
			HoleInfo(face::Pos<typename MESH::FaceType> const &pHole, int  const pHoleSize, vcg::Box3<typename MESH::ScalarType> &pHoleBB)
			{
				p=pHole;	
				size=pHoleSize;
				bb=pHoleBB;
			}

			HoleInfo(face::Pos<typename MESH::FaceType> const &pHole, int  const pHoleSize, vcg::Box3<typename MESH::ScalarType> &pHoleBB, int FI)
			{
				p=pHole;	
				size=pHoleSize;
				bb=pHoleBB;
				faceindex = FI;
			}


			typename face::Pos<typename MESH::FaceType> p;
			int size;
			vcg::Box3<typename MESH::ScalarType>  bb;
			int faceindex;

			void Refresh(MESH &m)
			{
				p.f = (typename MESH::FacePointer)(faceindex + &(*(m.face.begin())));
			}

			bool operator <  (const  HoleInfo & hh) const {return size <  hh.size;}
			bool operator >  (const  HoleInfo & hh) const {return size >  hh.size;}
			bool operator == (const  HoleInfo & hh) const {return size == hh.size;}
			bool operator != (const  HoleInfo & hh) const {return size != hh.size;}
			bool operator >= (const  HoleInfo & hh) const {return size >= hh.size;}
			bool operator <= (const  HoleInfo & hh) const {return size <= hh.size;}

			typename MESH::ScalarType Perimeter()
			{
				typename MESH::ScalarType sum=0;
				face::Pos<typename MESH::FaceType> ip = p;
				do
				{
					sum+=Distance(ip.v->cP(),ip.VFlip()->cP());
					ip.NextB();
				}
				while (ip != p);
				return sum;				
			}


			int CollectEdges(std::vector< SimpleEdge<MESH> > &EV)
			{
				assert(p.IsBorder());
				EV.clear();
				int tsz=0;
				face::Pos<typename MESH::FaceType> ip=p;
				face::Pos<typename MESH::FaceType> tp;

				do
				{
					// Stesso codice della nextb
					do
					{
						ip.NextE();
						EV.push_back(SimpleEdge<MESH>(ip)); // l'edge che sto scorrendo
						tp=ip;
						tp.FlipV();tp.FlipE();
						EV.push_back(SimpleEdge<MESH>(tp)); // l'edge della faccia su cui sono e opposto al vertice su cui ruoto
						tp.FlipF(); tp.FlipE();
						EV.push_back(SimpleEdge<MESH>(tp));  // gli altri due edge della faccia opposta a questa
						tp.FlipE();
						EV.push_back(SimpleEdge<MESH>(tp));
					}
					while(!ip.f->IsB(ip.z));
					ip.FlipV();
					++tsz;
				}
				while (ip != p);
				assert(tsz==size);

				return EV.size();
			}
		};

		template<class MESH>
			void FindHole(MESH &m, face::Pos<typename MESH::FaceType> ep, HoleInfo<MESH> &h)
		{
			if(!ep.IsBorder()) return;

			int holesize = 0;

			Box3<typename MESH::ScalarType> hbox;
			if(ep.v->IsR()) hbox.Add(ep.v->cP());
			face::Pos<typename MESH::FaceType> init;
			init = ep;
			do
			{
				ep.NextB();
				ep.f->SetV();
				if(ep.v->IsR()) hbox.Add(ep.v->cP());
				++holesize;
			}
			while (ep != init);
			h=HoleInfo<MESH>(ep,holesize,hbox);
		}

		template<class MESH,class STL_CONTAINER_HOLES>
			void FindHole(MESH &m, STL_CONTAINER_HOLES & H)
		{
			typename MESH::FaceIterator pf;
			int holesize;
			for (pf=m.face.begin(); pf!=m.face.end(); ++pf)
				if( !(*pf).IsD() && (*pf).IsW() )
					(*pf).ClearS();

			face::Pos<typename MESH::FaceType> ep;
			for (pf=m.face.begin(); pf!=m.face.end(); ++pf)
			{
				if( !(*pf).IsD() && !(*pf).IsS() && (*pf).IsR() )
				{
					for(int j=0; j<3; ++j)
						if( (*pf).IsB(j) && !(*pf).IsS() && (*pf).IsR() )
						{
							(*pf).SetS();
							ep.Set(&*pf, j, (*pf).V(j));
							holesize = 0;

							Box3<typename MESH::ScalarType> hbox;
							if(ep.v->IsR()) hbox.Add(ep.v->cP());
							face::Pos<typename MESH::FaceType> init;
							init = ep;
							do
							{
								ep.NextB();
								ep.f->SetS();
								if(ep.v->IsR()) hbox.Add(ep.v->cP());
								++holesize;
							}
							while (ep != init);
							H.push_back(HoleInfo<MESH>(ep,holesize,hbox));
							break;
						}
				}
			}
		};


		/*
		Un ear e' identificato da due hedge pos.

		i vertici dell'ear sono
		e0.FlipV().v
		e0.v
		e1.v

		Vale che e1== e0.NextB();
		e che e1.FlipV() == e0;

		Situazioni ear non manifold, e degeneri (buco triangolare) 


		T XXXXXXXXXXXXX    A        /XXXXX        B      en/XXXXX
		/XXXXXXXXXXXXXXX            /XXXXXX                /XXXXXX
		XXXXXXep==en XXX     ep\   /en XXXX               /e1 XXXX
		XXXXXX ----/| XX   ------ ----/| XX       ------ ----/|XXX
		XXXXXX|   /e1 XX   XXXXXX|   /e1 XX       XXXXXX|  o/e0 XX 
		XXXXXX|  /XXXXXX   XXXXXX|  /XXXXXX       XXXXXX|  /XXXXXX 
		XXX e0|o/XXXXXXX   XXX e0|o/XXXXXXX       XXX ep| /XXXXXXX
		XXX  \|/XXXXXXXX   XXX  \|/XXXXXXXX       XXX  \|/XXXXXXXX
		XXXXXXXXXXXXXXXX   XXXXXXXXXXXXXXXX       XXXXXXXXXXXXXXXX   

		*/
		template<class MSH_TYPE> class TrivialEar
		{
		public:
			face::Pos<typename MSH_TYPE::FaceType> e0;	 
			face::Pos<typename MSH_TYPE::FaceType> e1;	 
			typedef typename MSH_TYPE::ScalarType ScalarType;
			ScalarType quality;
			ScalarType angle;
			TrivialEar(){}
			TrivialEar(const face::Pos<typename MSH_TYPE::FaceType> & ep)
			{
				e0=ep;
				assert(e0.IsBorder());
				e1=e0;
				e1.NextB();
				ComputeQuality();
				computeAngle();
			}

			void computeAngle()
			{
				Point3f p1 = e0.VFlip()->P() - e0.v->P();
				Point3f p2 = e1.v->P() - e0.v->P();

				ScalarType  w = p2.Norm()*p1.Norm();
				if(w==0) angle =90;
				ScalarType p = (p2*p1);
				p= p/w;
				p = acos(p);
				if(p < -1) p = -1;
				if(p > 1) p = 1;

				Point3f t = p2^p1;
				ScalarType n = t* e0.v->N();
				if(n<0)
				{
					p = 2.0 *(float)M_PI - p;
				}
				angle = p;
			}

			inline bool operator < ( const TrivialEar & c ) const { return quality <  c.quality; }

			bool IsNull(){return e0.IsNull() || e1.IsNull();}
			void SetNull(){e0.SetNull();e1.SetNull();}
			void ComputeQuality()
			{ 
				ScalarType ar;
				ar = ( (e0.VFlip()->P() - e0.v->P()) ^ ( e1.v->P() - e0.v->P()) ).Norm() ;
				ScalarType area = (ar);

				ScalarType l1 = Distance( e0.v->P(),e1.v->P());
				ScalarType l2 = Distance( e0.v->P(),e0.VFlip()->P());
				ScalarType l3 = Distance( e0.VFlip()->P(),e1.v->P());

				quality = area / ( (l1 *l1) + (l2 * l2) + (l3 * l3) );
			};
			bool IsUpToDate()	{return (e0.IsBorder() && e1.IsBorder());};

			bool IsConvex(){return (angle > (float)M_PI);}


			bool Degen()
			{
				face::Pos<typename MSH_TYPE::FaceType>	ep=e0; ep.FlipV(); ep.NextB(); ep.FlipV(); // he precedente a e0 
				face::Pos<typename MSH_TYPE::FaceType>	en=e1; en.NextB();												 // he successivo a e1

				// caso ear degenere per buco triangolare
				if(ep==en) return true;//provo a togliere sto controllo
				// Caso ear non manifold a
				if(ep.v==en.v)	return true;
				// Caso ear non manifold b
				if(ep.VFlip()==e1.v) return true;

				return false;
			}

			bool Close(TrivialEar &ne0, TrivialEar &ne1, typename MSH_TYPE::FaceType * f)
			{
				// simple topological check
				if(e0.f==e1.f) {
					printf("Avoided bad ear");
					return false;
				}

				//usato per generare una delle due nuove orecchie.
				face::Pos<typename MSH_TYPE::FaceType>	ep=e0; ep.FlipV(); ep.NextB(); ep.FlipV(); // he precedente a e0 
				face::Pos<typename MSH_TYPE::FaceType>	en=e1; en.NextB();												 // he successivo a e1

				(*f).V(0) = e0.VFlip();
				(*f).V(1) = e0.v;
				(*f).V(2) = e1.v;

				(*f).FFp(0) = e0.f;
				(*f).FFi(0) = e0.z;
				(*f).FFp(1) = e1.f;
				(*f).FFi(1) = e1.z;
				(*f).FFp(2) = f;
				(*f).FFi(2) = 2;

				e0.f->FFp(e0.z)=f;
				e0.f->FFi(e0.z)=0;	

				e1.f->FFp(e1.z)=f;
				e1.f->FFi(e1.z)=1;	

				// caso ear degenere per buco triangolare
				if(ep==en)
				{
					printf("Closing the last triangle");
					f->FFp(2)=en.f;
					f->FFi(2)=en.z;
					en.f->FFp(en.z)=f;
					en.f->FFi(en.z)=2;
					ne0.SetNull();
					ne1.SetNull();
				}
				// Caso ear non manifold a
				else if(ep.v==en.v)
				{
					printf("Ear Non manif A\n");
					face::Pos<typename MSH_TYPE::FaceType>	enold=en;
					en.NextB();
					f->FFp(2)=enold.f;
					f->FFi(2)=enold.z;
					enold.f->FFp(enold.z)=f;
					enold.f->FFi(enold.z)=2;
					ne0=TrivialEar(ep);
					ne1=TrivialEar(en);
				}
				// Caso ear non manifold b
				else if(ep.VFlip()==e1.v)
				{
					printf("Ear Non manif B\n");
					face::Pos<typename MSH_TYPE::FaceType>	epold=ep; 
					ep.FlipV(); ep.NextB(); ep.FlipV();
					f->FFp(2)=epold.f;
					f->FFi(2)=epold.z;
					epold.f->FFp(epold.z)=f;
					epold.f->FFi(epold.z)=2;
					ne0=TrivialEar(ep);
					ne1=TrivialEar(en);
				}
				else // caso standard
					// Now compute the new ears;
				{
					ne0=TrivialEar(ep);
					ne1=TrivialEar(face::Pos<typename MSH_TYPE::FaceType>(f,2,e1.v));
				}

				return true;
			}
		};

		//Ear with Leipa's quality policy
		template<class MSH_TYPE> class LeipaEar
		{
		public:
			face::Pos<typename MSH_TYPE::FaceType> e0;	 
			face::Pos<typename MSH_TYPE::FaceType> e1;	 
			typedef typename MSH_TYPE::ScalarType ScalarType;
			ScalarType angle;
			ScalarType dihedral;
			ScalarType area;
			LeipaEar(){}
			LeipaEar(const face::Pos<typename MSH_TYPE::FaceType> & ep)
			{
				e0=ep;
				assert(e0.IsBorder());
				e1=e0;
				e1.NextB();
				ComputeQuality();
				computeAngle();
			}

			void computeAngle()
			{
				Point3f p1 = e0.VFlip()->P() - e0.v->P();
				Point3f p2 = e1.v->P() - e0.v->P();

				ScalarType  w = p2.Norm()*p1.Norm();
				if(w==0) angle =90;
				ScalarType p = (p2*p1);
				p= p/w;
				p = acos(p);
				if(p < -1) p = -1;
				if(p > 1) p = 1;

				Point3f t = p2^p1;
				ScalarType n = t* e0.v->N();
				if(n<0)
				{
					p = 2.0 *(float)M_PI - p;
				}
				angle = p;
			}

			// Nota: minori invertiti		
			inline bool operator <  ( const LeipaEar & c ) const 
			{ 
				if(dihedral < c.dihedral)return true;
				else return ((dihedral == c.dihedral) && (area < c.area));
			}

			bool IsNull(){return e0.IsNull() || e1.IsNull();}
			void SetNull(){e0.SetNull();e1.SetNull();}
			void ComputeQuality()
			{ 
				//comute quality by (dihedral ancgle, area/sum(edge^2) )
				Point3f n1 = (e0.v->N() + e1.v->N() + e0.VFlip()->N() ) / 3;
				face::Pos<typename MSH_TYPE::FaceType> tmp = e1;
				tmp.FlipE();tmp.FlipV();
				Point3f n2=(e1.VFlip()->N() + e1.v->N() + tmp.v->N() ) / 3;
				tmp = e0;
				tmp.FlipE(); tmp.FlipV();
				Point3f n3=(e0.VFlip()->N() + e0.v->N() + tmp.v->N() ) / 3; 
				dihedral = std::max(Angle(n1,n2),Angle(n1,n3));

				ScalarType ar;
				ar = ( (e0.VFlip()->P() - e0.v->P()) ^ ( e1.v->P() - e0.v->P()) ).Norm() ;

				/*ScalarType l1 = Distance( e0.v->P(),e1.v->P());
				ScalarType l2 = Distance( e0.v->P(),e0.VFlip()->P());
				ScalarType l3 = Distance( e0.VFlip()->P(),e1.v->P());*/


				area = ar ;// ( (l1 *l1) + (l2 * l2) + (l3 * l3) );


			};//dovrebbe
			bool IsUpToDate()	{return (e0.IsBorder() && e1.IsBorder());};
			bool IsConvex(){return angle > (float)M_PI;}
			bool Degen()
			{
				face::Pos<typename MSH_TYPE::FaceType>	ep=e0; ep.FlipV(); ep.NextB(); ep.FlipV(); // he precedente a e0 
				face::Pos<typename MSH_TYPE::FaceType>	en=e1; en.NextB();												 // he successivo a e1

				// caso ear degenere per buco triangolare
				if(ep==en) return true;//provo a togliere sto controllo
				// Caso ear non manifold a
				if(ep.v==en.v)	return true;
				// Caso ear non manifold b
				if(ep.VFlip()==e1.v) return true;

				return false;
			}

			bool Close(LeipaEar &ne0, LeipaEar &ne1, typename MSH_TYPE::FaceType * f)
			{
				// simple topological check
				if(e0.f==e1.f) {
					printf("Avoided bad ear");
					return false;
				}

				//usato per generare una delle due nuove orecchie.
				face::Pos<typename MSH_TYPE::FaceType>	ep=e0; ep.FlipV(); ep.NextB(); ep.FlipV(); // he precedente a e0 
				face::Pos<typename MSH_TYPE::FaceType>	en=e1; en.NextB();												 // he successivo a e1

				(*f).V(0) = e0.VFlip();
				(*f).V(1) = e0.v;
				(*f).V(2) = e1.v;

				(*f).FFp(0) = e0.f;
				(*f).FFi(0) = e0.z;
				(*f).FFp(1) = e1.f;
				(*f).FFi(1) = e1.z;
				(*f).FFp(2) = f;
				(*f).FFi(2) = 2;

				e0.f->FFp(e0.z)=f;
				e0.f->FFi(e0.z)=0;	

				e1.f->FFp(e1.z)=f;
				e1.f->FFi(e1.z)=1;	

				// caso ear degenere per buco triangolare
				if(ep==en)
				{
					printf("Closing the last triangle");
					f->FFp(2)=en.f;
					f->FFi(2)=en.z;
					en.f->FFp(en.z)=f;
					en.f->FFi(en.z)=2;
					ne0.SetNull();
					ne1.SetNull();
				}
				// Caso ear non manifold a
				else if(ep.v==en.v)
				{
					printf("Ear Non manif A\n");
					face::Pos<typename MSH_TYPE::FaceType>	enold=en;
					en.NextB();
					f->FFp(2)=enold.f;
					f->FFi(2)=enold.z;
					enold.f->FFp(enold.z)=f;
					enold.f->FFi(enold.z)=2;
					ne0=LeipaEar(ep);
					ne1=LeipaEar(en);
				}
				// Caso ear non manifold b
				else if(ep.VFlip()==e1.v)
				{
					printf("Ear Non manif B\n");
					face::Pos<typename MSH_TYPE::FaceType>	epold=ep; 
					ep.FlipV(); ep.NextB(); ep.FlipV();
					f->FFp(2)=epold.f;
					f->FFi(2)=epold.z;
					epold.f->FFp(epold.z)=f;
					epold.f->FFi(epold.z)=2;
					ne0=LeipaEar(ep);
					ne1=LeipaEar(en);
				}
				else // caso standard
					// Now compute the new ears;
				{
					ne0=LeipaEar(ep);
					ne1=LeipaEar(face::Pos<typename MSH_TYPE::FaceType>(f,2,e1.v));
				}

				return true;
			}
		};
		//Ear for selfintersection algorithm
		template<class MSH_TYPE> class SelfIntersection
		{
		public:
			face::Pos<typename MSH_TYPE::FaceType> e0;	 
			face::Pos<typename MSH_TYPE::FaceType> e1;	 
			typedef typename MSH_TYPE::ScalarType ScalarType;
			ScalarType quality;
			ScalarType angle;
			SelfIntersection(){}
			SelfIntersection(const face::Pos<typename MSH_TYPE::FaceType> & ep)
			{
				e0=ep;
				assert(e0.IsBorder());
				e1=e0;
				e1.NextB();
				ComputeQuality();
				computeAngle();
			}

			inline bool operator <  ( const SelfIntersection & c ) const 
			{ 
				return (quality < c.quality);
			}
			void computeAngle()
			{
				Point3f p1 = e0.VFlip()->P() - e0.v->P();
				Point3f p2 = e1.v->P() - e0.v->P();

				ScalarType  w = p2.Norm()*p1.Norm();
				if(w==0) angle =90;
				ScalarType p = (p2*p1);
				p= p/w;
				p = acos(p);
				if(p < -1) p = -1;
				if(p > 1) p = 1;

				Point3f t = p2^p1;
				ScalarType n = t* e0.v->N();
				if(n<0)
				{
					p = 2.0 *(float)M_PI - p;
				}
				angle = p;
			}

			bool IsNull(){return e0.IsNull() || e1.IsNull();}
			void SetNull(){e0.SetNull();e1.SetNull();}
			void ComputeQuality()
			{ 
				ScalarType ar;
				ar = ( (e0.VFlip()->P() - e0.v->P()) ^ ( e1.v->P() - e0.v->P()) ).Norm() ;
				ScalarType area = (ar);

				ScalarType l1 = Distance( e0.v->P(),e1.v->P());
				ScalarType l2 = Distance( e0.v->P(),e0.VFlip()->P());
				ScalarType l3 = Distance( e0.VFlip()->P(),e1.v->P());

				quality = area / ( (l1 *l1) + (l2 * l2) + (l3 * l3) );
			};
			bool IsUpToDate()	{return (e0.IsBorder() && e1.IsBorder());};
			bool IsConvex(){ return (angle > (float)M_PI);}
			bool Degen()
			{
				face::Pos<typename MSH_TYPE::FaceType>	ep=e0; ep.FlipV(); ep.NextB(); ep.FlipV(); // he precedente a e0 
				face::Pos<typename MSH_TYPE::FaceType>	en=e1; en.NextB();												 // he successivo a e1

				// caso ear degenere per buco triangolare
				if(ep==en) return true;//provo a togliere sto controllo
				// Caso ear non manifold a
				if(ep.v==en.v)	return true;
				// Caso ear non manifold b
				if(ep.VFlip()==e1.v) return true;

				return false;
			}

			bool Close(SelfIntersection &ne0, SelfIntersection &ne1, typename MSH_TYPE::FaceType * f, std::vector<typename MSH_TYPE::FaceType> vf)
			{
				// simple topological check
				if(e0.f==e1.f) {
					printf("Avoided bad ear");
					return false;
				}

				face::Pos<typename MSH_TYPE::FaceType>	ep=e0; ep.FlipV(); ep.NextB(); ep.FlipV(); // he precedente a e0 
				face::Pos<typename MSH_TYPE::FaceType>	en=e1; en.NextB();												 // he successivo a e1
				//costruisco la faccia e poi testo, o copio o butto via.				
				(*f).V(0) = e0.VFlip();
				(*f).V(1) = e0.v;
				(*f).V(2) = e1.v;

				(*f).FFp(0) = e0.f;
				(*f).FFi(0) = e0.z;
				(*f).FFp(1) = e1.f;
				(*f).FFi(1) = e1.z;
				(*f).FFp(2) = f;
				(*f).FFi(2) = 2;

				int a1, a2;
				a1=e0.z;
				a2=e1.z;

				e0.f->FFp(e0.z)=f;
				e0.f->FFi(e0.z)=0;	

				e1.f->FFp(e1.z)=f;
				e1.f->FFi(e1.z)=1;
				typename std::vector<typename MSH_TYPE::FaceType>::iterator it;
				for(it = vf.begin();it!= vf.end();++it)
				{
					if(!it->IsD())
						if(		tri::Clean<MSH_TYPE>::TestIntersection(&(*f),&(*it)))
						{
							//rimetto a posto

							e0.f->FFp(e0.z)=e0.f;
							e0.f->FFi(e0.z)=a1;	

							e1.f->FFp(e1.z)=e1.f;
							e1.f->FFi(e1.z)=a2;

							return false;
						}
				}
				// caso ear degenere per buco triangolare
				if(ep==en)
				{
					printf("Closing the last triangle");
					f->FFp(2)=en.f;
					f->FFi(2)=en.z;
					en.f->FFp(en.z)=f;
					en.f->FFi(en.z)=2;
					ne0.SetNull();
					ne1.SetNull();
				}
				// Caso ear non manifold a
				else if(ep.v==en.v)
				{
					printf("Ear Non manif A\n");
					face::Pos<typename MSH_TYPE::FaceType>	enold=en;
					en.NextB();
					f->FFp(2)=enold.f;
					f->FFi(2)=enold.z;
					enold.f->FFp(enold.z)=f;
					enold.f->FFi(enold.z)=2;
					ne0=SelfIntersection(ep);
					ne1=SelfIntersection(en);
				}
				// Caso ear non manifold b
				else if(ep.VFlip()==e1.v)
				{
					printf("Ear Non manif B\n");
					face::Pos<typename MSH_TYPE::FaceType>	epold=ep; 
					ep.FlipV(); ep.NextB(); ep.FlipV();
					f->FFp(2)=epold.f;
					f->FFi(2)=epold.z;
					epold.f->FFp(epold.z)=f;
					epold.f->FFi(epold.z)=2;
					ne0=SelfIntersection(ep);
					ne1=SelfIntersection(en);
				}
				else // caso standard
					// Now compute the new ears;
				{
					ne0=SelfIntersection(ep);
					ne1=SelfIntersection(face::Pos<typename MSH_TYPE::FaceType>(f,2,e1.v));
				}
				return true;
			}
		};



		// Funzione principale per chiudier un buco in maniera topologicamente corretta.
		// Gestisce situazioni non manifold ragionevoli 
		// (tutte eccetto quelle piu' di 2 facce per 1 edge).
		// Controlla che non si generino nuove situazioni non manifold chiudendo orecchie
		// che sottendono un edge che gia'esiste.
		template<class MESH>
			tri::HoleInfo<MESH> getHoleInfo(MESH &m, face::Pos<typename MESH::FaceType> sp, 
			face::Pos<typename MESH::FaceType> fp, 
			int UBIT)
		{
			int holesize=0;

			Box3<typename MESH::ScalarType> hbox;
			hbox.Add(sp.v->cP());

			do
			{
				sp.f->SetUserBit(UBIT);
				hbox.Add(sp.v->cP());
				++holesize;
				sp.NextB();
				assert(sp.IsBorder());
			}while(sp != fp);

			int tmp = ((int)(sp.f - &(*(m.face.begin()))));
			return tri::HoleInfo<MESH>(sp,holesize,hbox, tmp );
		}

		template<class MESH,class EAR , class VECTOR_EAR>
			void refreshHole(MESH &m, VECTOR_EAR &ve, face::Pos<typename MESH::FaceType> &fp)
		{
			face::Pos<typename MESH::FaceType> ff = fp;

			do{
				ve.push_back(EAR(fp));
				fp.NextB();//semmai da provare a sostituire il codice della NextB();
				assert(fp.IsBorder());
			}while(fp!=ff);

		}

		template <class MESH, class EAR>
			void fillHoleEar(MESH &m, tri::HoleInfo<MESH> &h ,int UBIT)
		{
			//Aggiungo le facce e aggiorno il puntatore alla faccia!
			std::vector<typename MESH::FacePointer *> app;
			app.push_back( &h.p.f );
			typename MESH::FaceIterator f = tri::Allocator<MESH>::AddFaces(m, h.size-2, app);
			h.Refresh(m);	
			assert(h.p.IsBorder());//test fondamentale altrimenti qualcosa s'e' rotto!

			std::vector<EAR > H; //vettore di orecchie

			H.reserve(h.size);

			//prendo le informazioni sul buco
			refreshHole<MESH,EAR, std::vector<EAR> >(m,H,h.p);

			bool fitted = false;
			int cnt=h.size;
			typename MESH::FaceIterator tmp;

			make_heap(H.begin(), H.end());

			while( cnt > 2 && !H.empty()  && !fitted) //finche' il buco non e' chiuso o non ci sono piu' orecchie da analizzare
			{

				pop_heap(H.begin(), H.end());		

				EAR en0,en1;

				typename MESH::FaceIterator Fadd = f;

				if(H.back().IsUpToDate() && !H.back().IsConvex())	
				{
					if(H.back().Degen()){ 
						// Nota che nel caso di ear degeneri si DEVE permettere la creazione di un edge che gia'esiste
						printf("\n -> Evitata orecchia brutta!");
					}
					else 
					{

						if(H.back().Close(en0,en1,&*f))
						{
							if(!en0.IsNull()){
								H.push_back(en0);
								push_heap( H.begin(), H.end());
							}
							if(!en1.IsNull()){
								H.push_back(en1);
								push_heap( H.begin(), H.end());
							}
							--cnt;
							f->SetUserBit(UBIT);
							++f;
							fitted = true;
						}
					}

					if(cnt == 3 && !fitted)
					{//ultimo buco o unico buco
						if(H.back().Close(en0,en1,&*f))
						{
							--cnt;
							tmp = f;
							++f;
						}
					}
				}//is update()
				fitted = false;
				//non ho messo il triangolo quindi tolgo l'orecchio e continuo
				H.pop_back();
			}//fine del while principale

			while(f!=m.face.end())
			{
				(*f).SetD();
				++f;
				m.fn--;
			}

		}


		template<class MESH, class EAR>
			void holeFillingEar(MESH &m, int sizeHole,bool Selected = false)
		{
			typename MESH::FaceIterator fi;
			std::vector<tri::HoleInfo<MESH> > vinfo;
			int UBIT = MESH::FaceType::LastBitFlag();

			for(fi = m.face.begin(); fi!=m.face.end(); ++fi)
			{
				if(!(*fi).IsD())
				{
					if(Selected && !(*fi).IsS())
					{
						//se devo considerare solo i triangoli selezionati e 
						//quello che sto considerando non lo e' lo marchio e vado avanti
						(*fi).SetUserBit(UBIT);
					}
					else
					{
						if( !(*fi).IsUserBit(UBIT) )
						{
							(*fi).SetUserBit(UBIT);
							for(int j =0; j<3 ; ++j)
							{
								if( (*fi).IsB(j) )
								{//Trovato una faccia di bordo non ancora visitata.
									face::Pos<typename MESH::FaceType> sp(&*fi, j, (*fi).V(j));

									//	if(!(*fi).IsR())return;
									tri::HoleInfo<MESH> HI = getHoleInfo<MESH>(m,sp,sp,  UBIT);
									//ho recuperato l'inofrmazione su tutto il buco
									vinfo.push_back(HI);
								}
							}//for sugli edge del triangolo
						}//se e' gia stato visitato
					}//S & !S
				}//!IsD()
			}//for principale!!!

			typename std::vector<typename tri::HoleInfo<MESH> >::iterator ith;
			typename tri::HoleInfo<MESH> app;
			for(ith = vinfo.begin(); ith!= vinfo.end(); ++ith)
			{
				app=(tri::HoleInfo<MESH>)*ith;
				if(app.size < sizeHole){		
					fillHoleEar<MESH, EAR >(m, app,UBIT);
				}
			}

			for(fi = m.face.begin(); fi!=m.face.end(); ++fi)
			{
				if(!(*fi).IsD())
					(*fi).ClearUserBit(UBIT);
			}
		}

		/*
		FillHoleSelfIntersection
		*/
		template <class MESH, class EAR>
			void fillHoleInt(MESH &m, tri::HoleInfo<MESH> &h ,int UBIT, std::vector<typename MESH::FaceType > vf)
		{
			//Aggiungo le facce e aggiorno il puntatore alla faccia!
			std::vector<typename MESH::FacePointer *> app;
			app.push_back( &h.p.f );
			typename MESH::FaceIterator f = tri::Allocator<MESH>::AddFaces(m, h.size-2, app);
			h.Refresh(m);	
			assert(h.p.IsBorder());//test fondamentale altrimenti qualcosa s'e' rotto!
			std::vector<EAR > H; //vettore di orecchie
			H.reserve(h.size);

			//prendo le informazioni sul buco
			tri::refreshHole<MESH,EAR, std::vector<EAR> >(m,H,h.p);

			bool fitted = false;
			int cnt=h.size;
			typename MESH::FaceIterator tmp;

			make_heap(H.begin(), H.end());
			//finche' il buco non e' chiuso o non ci sono piu' orecchie da analizzare.
			while( cnt > 2 && !H.empty() ) 
			{
				pop_heap(H.begin(), H.end());		
				EAR en0,en1;
				typename MESH::FaceIterator Fadd = f;
				if(H.back().IsUpToDate() && !H.back().IsConvex())	
				{
					if(H.back().Degen()){ 
						// Nota che nel caso di ear degeneri si DEVE permettere la creazione di un edge che gia'esiste.
						printf("\n -> Evitata orecchia brutta!");
					}
					else 
					{
						if(H.back().Close(en0,en1,&*f,vf))
						{
							if(!en0.IsNull()){
								H.push_back(en0);
								push_heap( H.begin(), H.end());
							}
							if(!en1.IsNull()){
								H.push_back(en1);
								push_heap( H.begin(), H.end());
							}
							--cnt;
							f->SetUserBit(UBIT);
							vf.push_back(*f);
							++f;
							fitted = true;
						}
					}
					//ultimo buco o unico buco.
					if(cnt == 3 && !fitted)
					{
						if(H.back().Close(en0,en1,&*f,vf))
						{
							--cnt;
							tmp = f;
							vf.push_back(*f);
							++f;
						}
					}
				}//is update()
				fitted = false;
				//non ho messo il triangolo quindi tolgo l'orecchio e continuo.
				H.pop_back();
			}//fine del while principale.
			//tolgo le facce non utilizzate.
			while(f!=m.face.end())
			{
				(*f).SetD();
				++f;
				m.fn--;
			}

		}

		//hole filling selfintersection main algorithm
		template<class MESH, class EAR>
			void holeFillingIntersection(MESH &m, int sizeHole,bool Selected = false)
		{
			typename MESH::FaceIterator fi;
			std::vector<typename tri::HoleInfo<MESH> > vinfo;
			int UBIT = fi->LastBitFlag();

			for(fi = m.face.begin(); fi!=m.face.end(); ++fi)
			{
				if(!(*fi).IsD())
				{
					if(Selected && !(*fi).IsS())
					{
						//se devo considerare solo i triangoli selezionati e 
						//quello che sto considerando non lo e' lo marchio e vado avanti
						(*fi).SetUserBit(UBIT);
					}
					else
					{
						if( !(*fi).IsUserBit(UBIT) )
						{
							(*fi).SetUserBit(UBIT);
							for(int j =0; j<3 ; ++j)
							{
								if( (*fi).IsB(j) )
								{//Trovato una faccia di bordo non ancora visitata.
									face::Pos<typename MESH::FaceType> sp(&*fi, j, (*fi).V(j));

									//	if(!(*fi).IsR())return;
									tri::HoleInfo<MESH> HI = tri::getHoleInfo<MESH>(m,sp,sp,  UBIT);
									//ho recuperato l'inofrmazione su tutto il buco
									vinfo.push_back(HI);
								}
							}//for sugli edge del triangolo
						}//se e' gia stato visitato
					}//S & !S
				}//!IsD()
			}//for principale!!!

			std::vector<typename MESH::FaceType > vf;
			face::Pos<typename MESH::FaceType>sp;
			face::Pos<typename MESH::FaceType>ap;
			typename std::vector<tri::HoleInfo<MESH> >::iterator ith;
			tri::HoleInfo<MESH> app;
			for(ith = vinfo.begin(); ith!= vinfo.end(); ++ith)
			{
				app=(tri::HoleInfo<MESH>)*ith;
				if(app.size < sizeHole){		
					app.Refresh(m);
					//colleziono il ring intorno al buco per poi fare il test sul'intersezione
					sp = app.p;
					do
					{
						ap = sp;
						do
						{
							ap.FlipE();
							ap.FlipF();
							vf.push_back(*ap.f);
						}while(!ap.IsBorder());
						sp.NextB();

					}while(sp != app.p);	

					fillHoleInt<MESH, EAR >(m, app,UBIT,vf);
					vf.clear();
				}
			}
			for(fi = m.face.begin(); fi!=m.face.end(); ++fi)
			{
				if(!(*fi).IsD())
					(*fi).ClearUserBit(UBIT);
			}
		}
	} // end namespace
}
#endif
