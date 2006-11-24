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
Revision 1.16  2006/11/22 13:43:28  giec
Code refactory and added minimum weight triangolation.

Revision 1.15  2006/11/13 10:11:38  giec
Clear some useless code

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
#include <vcg/complex/trimesh/clean.h>
#include <vcg/space/point3.h>
#include <vector>

#define FLT_MAX     3.402823466e+38F        /* max float rappresentable */
/*
Questa Classe serve per gestire la non duplicazione degli edge durante la chiusura 
di un buco.
*/
namespace vcg {
	namespace tri {


		template<class MESH>
		class HoleInfo
		{
		public: 
			HoleInfo(){}
			HoleInfo(face::Pos<typename MESH::FaceType> const &pHole, int  const pHoleSize, Box3<typename MESH::ScalarType> &pHoleBB)
			{
				p=pHole;	
				size=pHoleSize;
				bb=pHoleBB;
			}

			HoleInfo(face::Pos<typename MESH::FaceType> const &pHole, int  const pHoleSize, Box3<typename MESH::ScalarType> &pHoleBB, int FI)
			{
				p=pHole;	
				size=pHoleSize;
				bb=pHoleBB;
				faceindex = FI;
			}


			typename face::Pos<typename MESH::FaceType> p;
			int size;
			Box3<typename MESH::ScalarType>  bb;
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

		};

		//function prototype
		template <class MESH>
		int GetHoleInfo(MESH &m,bool Selected ,std::vector<typename tri::HoleInfo<MESH> >& VHI);
		
		template<class MESH>
		void triangulate(std::vector<typename MESH::VertexPointer > &m,int i, int j, std::vector< std::vector<int> > vi, 
										 std::vector<face::Pos<typename MESH::FaceType> > vv);

		template <class MESH>
		void getBoundHole (face::Pos<typename MESH::FaceType> sp,std::vector<face::Pos<typename MESH::FaceType> >&ret);
		


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
			std::vector<typename MSH_TYPE::FaceType>* vf;
			TrivialEar(){}
			TrivialEar(const face::Pos<typename MSH_TYPE::FaceType> & ep)
			{
				e0=ep;
				assert(e0.IsBorder());
				e1=e0;
				e1.NextB();
				ComputeQuality();
				ComputeAngle();
			}

			void SetAdiacenseRing(std::vector<typename MSH_TYPE::FaceType>* ar){vf = ar;}

			void ComputeAngle()
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

			virtual inline bool operator < ( const TrivialEar & c ) const { return quality <  c.quality; }

			bool IsNull(){return e0.IsNull() || e1.IsNull();}
			void SetNull(){e0.SetNull();e1.SetNull();}
			virtual	void ComputeQuality()
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

			virtual bool Close(TrivialEar &ne0, TrivialEar &ne1, typename MSH_TYPE::FaceType * f)
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

		//Ear with FillHoleMinimumWeight's quality policy
		template<class MSH_TYPE> class MinimumWeightEar : public TrivialEar<MSH_TYPE>
		{
		public:
			typename MSH_TYPE::ScalarType dihedral;
			typename MSH_TYPE::ScalarType area;
			MinimumWeightEar(){}
			MinimumWeightEar(const face::Pos<typename MSH_TYPE::FaceType> & ep)
			{
				this->e0=ep;
				assert(this->e0.IsBorder());
				this->e1=this->e0;
				this->e1.NextB();
				this->ComputeQuality();
				this->ComputeAngle();
			}

			virtual inline bool operator <  ( const MinimumWeightEar & c ) const 
			{ 
				if(dihedral < c.dihedral)return true;
				else return ((dihedral == c.dihedral) && (area < c.area));
			}

			virtual void ComputeQuality()
			{ 
				//comute quality by (dihedral ancgle, area/sum(edge^2) )
				Point3f n1 = (this->e0.v->N() + this->e1.v->N() + this->e0.VFlip()->N() ) / 3;
				face::Pos<typename MSH_TYPE::FaceType> tmp = this->e1;
				tmp.FlipE();tmp.FlipV();
				Point3f n2=(this->e1.VFlip()->N() + this->e1.v->N() + tmp.v->N() ) / 3;
				tmp = this->e0;
				tmp.FlipE(); tmp.FlipV();
				Point3f n3=(this->e0.VFlip()->N() + this->e0.v->N() + tmp.v->N() ) / 3; 
				dihedral = std::max(Angle(n1,n2),Angle(n1,n3));

				typename MSH_TYPE::ScalarType ar;
				ar = ( (this->e0.VFlip()->P() - this->e0.v->P()) ^ ( this->e1.v->P() - this->e0.v->P()) ).Norm() ;

				area = ar ;
			}

		};
		//Ear for selfintersection algorithm
		template<class MSH_TYPE> class SelfIntersectionEar : public TrivialEar<MSH_TYPE>
		{
		public:

			SelfIntersectionEar(){}
			SelfIntersectionEar(const face::Pos<typename MSH_TYPE::FaceType> & ep)
			{
				this->e0=ep;
				assert(this->e0.IsBorder());
				this->e1=this->e0;
				this->e1.NextB();
				this->ComputeQuality();
				this->ComputeAngle();
			}

			virtual bool Close(SelfIntersectionEar &ne0, SelfIntersectionEar &ne1, typename MSH_TYPE::FaceType * f)
			{
				// simple topological check
				if(this->e0.f==this->e1.f) {
					printf("Avoided bad ear");
					return false;
				}

				face::Pos<typename MSH_TYPE::FaceType>	ep=this->e0; ep.FlipV(); ep.NextB(); ep.FlipV(); // he precedente a e0 
				face::Pos<typename MSH_TYPE::FaceType>	en=this->e1; en.NextB();												 // he successivo a e1
				//costruisco la faccia e poi testo, o copio o butto via.				
				(*f).V(0) = this->e0.VFlip();
				(*f).V(1) = this->e0.v;
				(*f).V(2) = this->e1.v;

				(*f).FFp(0) = this->e0.f;
				(*f).FFi(0) = this->e0.z;
				(*f).FFp(1) = this->e1.f;
				(*f).FFi(1) = this->e1.z;
				(*f).FFp(2) = f;
				(*f).FFi(2) = 2; 

				int a1, a2;
				a1= this->e0.z;
				a2= this->e1.z;

				this->e0.f->FFp(this->e0.z)=f;
				this->e0.f->FFi(this->e0.z)=0;	
				
				this->e1.f->FFp(this->e1.z)=f;
				this->e1.f->FFi(this->e1.z)=1;
				typename std::vector<typename MSH_TYPE::FaceType>::iterator it;
				for(it = (*	this->vf).begin();it!= (*	this->vf).end();++it)
				{
					if(!it->IsD())
						if(		tri::Clean<MSH_TYPE>::TestIntersection(&(*f),&(*it)))
						{
							this->e0.f->FFp(this->e0.z)= this->e0.f;
							this->e0.f->FFi(this->e0.z)=a1;	

							this->e1.f->FFp(this->e1.z)=this->e1.f;
							this->e1.f->FFi(this->e1.z)=a2;
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
					ne0=SelfIntersectionEar(ep);
					ne0.SetAdiacenseRing(this->vf);
					ne1=SelfIntersectionEar(en);
					ne1.SetAdiacenseRing(this->vf);
				}
				// Caso ear non manifold b
				else if(ep.VFlip()==this->e1.v)
				{
					printf("Ear Non manif B\n");
					face::Pos<typename MSH_TYPE::FaceType>	epold=ep; 
					ep.FlipV(); ep.NextB(); ep.FlipV();
					f->FFp(2)=epold.f;
					f->FFi(2)=epold.z;
					epold.f->FFp(epold.z)=f;
					epold.f->FFi(epold.z)=2;
					ne0=SelfIntersectionEar(ep);
					ne0.SetAdiacenseRing(this->vf);
					ne1=SelfIntersectionEar(en);
					ne1.SetAdiacenseRing(this->vf);
				}
				else// Now compute the new ears;
				{
					ne0=SelfIntersectionEar(ep);
					ne0.SetAdiacenseRing(this->vf);
					ne1=SelfIntersectionEar(face::Pos<typename MSH_TYPE::FaceType>(f,2,this->e1.v));
					ne1.SetAdiacenseRing(this->vf);
				}
				return true;
			}
		};

		// Funzione principale per chiudier un buco in maniera topologicamente corretta.
		// Gestisce situazioni non manifold ragionevoli 
		// (tutte eccetto quelle piu' di 2 facce per 1 edge).
		// Controlla che non si generino nuove situazioni non manifold chiudendo orecchie
		// che sottendono un edge che gia'esiste.

		template <class MESH, class EAR>
			void FillHoleEar(MESH &m, tri::HoleInfo<MESH> &h ,int UBIT, std::vector<typename MESH::FaceType > *vf =0)
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
			face::Pos<typename MESH::FaceType> ff = h.p;
			face::Pos<typename MESH::FaceType> fp = h.p;
			do{
				EAR app = EAR(fp);
				app.SetAdiacenseRing(vf);
				H.push_back( app );
				fp.NextB();//semmai da provare a sostituire il codice della NextB();
				assert(fp.IsBorder());
			}while(fp!=ff);

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
				if(H.back().IsUpToDate() && H.back().IsConvex())	
				{
					if(H.back().Degen()){ 
						// Nota che nel caso di ear degeneri si DEVE permettere la creazione di un edge che gia'esiste.
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
							if(vf != 0)	(*vf).push_back(*f);
							++f;
							fitted = true;
						}
					}
					//ultimo buco o unico buco.
					if(cnt == 3 && !fitted)
					{
						if(H.back().Close(en0,en1,&*f))
						{
							--cnt;
							tmp = f;
							if(vf != 0)(*vf).push_back(*f);
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

		template<class MESH, class EAR>
			void holeFillingEar(MESH &m, int sizeHole,bool Selected = false)
		{
			std::vector<typename tri::HoleInfo<MESH> > vinfo;
			int UBIT = GetHoleInfo<MESH>(m, Selected,vinfo);

			typename std::vector<typename tri::HoleInfo<MESH> >::iterator ith;
			typename tri::HoleInfo<MESH> app;
			for(ith = vinfo.begin(); ith!= vinfo.end(); ++ith)
			{
				app=(tri::HoleInfo<MESH>)*ith;
				if(app.size < sizeHole){		
					FillHoleEar<MESH, EAR >(m, app,UBIT);
				}
			}
			typename MESH::FaceIterator fi;
			for(fi = m.face.begin(); fi!=m.face.end(); ++fi)
			{
				if(!(*fi).IsD())
					(*fi).ClearUserBit(UBIT);
			}
		}

		template<class MESH, class EAR>
			void holeFillingIntersection(MESH &m, int sizeHole,bool Selected = false)
		{
			std::vector<typename tri::HoleInfo<MESH> > vinfo;
			int UBIT = GetHoleInfo<MESH>(m, Selected,vinfo);
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

					FillHoleEar<MESH, EAR >(m, app,UBIT,&vf);
					vf.clear();
				}
			}
			typename MESH::FaceIterator fi;
			for(fi = m.face.begin(); fi!=m.face.end(); ++fi)
			{
				if(!(*fi).IsD())
					(*fi).ClearUserBit(UBIT);
			}
		}

		template <class MESH>
			int GetHoleInfo(MESH &m,bool Selected ,std::vector<typename tri::HoleInfo<MESH> >& VHI)
		{
			typename MESH::FaceIterator fi;
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
									face::Pos<typename MESH::FaceType> fp=sp;
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

									//ho recuperato l'inofrmazione su tutto il buco
									VHI.push_back( tri::HoleInfo<MESH>(sp,holesize,hbox, tmp) );
								}
							}//for sugli edge del triangolo
						}//se e' gia stato visitato
					}//S & !S
				}//!IsD()
			}//for principale!!!
			return UBIT;
		}

		//Minimum Weight Algorithm
		class Weight
		{
		public:

			Weight() { ang = 180; ar = FLT_MAX ;}
			Weight( float An, float Ar ) { ang=An ; ar= Ar;}
			~Weight() {}

			float angle() const { return ang; }
			float area()  const { return ar; }

			Weight operator+( const Weight & other ) const {return Weight( std::max( angle(), other.angle() ), area() + other.area());}
			bool operator<( const Weight & rhs ) const {return ( angle() < rhs.angle() ||(angle() == rhs.angle() && area() < rhs.area()));	}

		private:
			float ang;
			float ar;
		};
		/*     
		 \ /        \/ 
		v1*---------*v4
		 / \       /
hole/   \     /
	 / 	   \   /
	/ear		\ /
-*---------*-
/ v3      v2\
		*/
		template <class MESH>
			float ComputeDihedralAngle(typename MESH::VertexPointer v1,typename MESH::VertexPointer v2,
			typename MESH::VertexPointer v3,typename MESH::VertexPointer v4)
		{
			typename MESH::CoordType  n1 = ((v1->P() - v2->P()) ^ (v3->P() - v1->P()) ).Normalize();
			typename MESH::CoordType	 n2 = ((v2->P() - v1->P()) ^ (v4->P() - v2->P()) ).Normalize();
			typename MESH::ScalarType t = (n1 * n2 )  ;
			return ( acos(t)* 180.0 / M_PI);
		}

		template<class MESH>
			bool existEdge(face::Pos<typename MESH::FaceType> pi,face::Pos<typename MESH::FaceType> pf)
		{
			face::Pos<typename MESH::FaceType> app = pi;
			face::Pos<typename MESH::FaceType> appF = pi;
			face::Pos<typename MESH::FaceType> tmp;
			assert(pi.IsBorder());

			appF.NextB();
			appF.FlipV();

			do
			{
				tmp = app;
				tmp.FlipV();
				if(tmp.v == pf.v)
					return true;
				app.FlipE();
				app.FlipF();

				if(app == pi)return false;
			}while(app != appF);
			return false;
		}

		template<class MESH>
			Weight computeWeight( int i, int j, int k,
			std::vector<face::Pos<typename MESH::FaceType> > pv,
			std::vector< std::vector< int > >  v)
		{
			face::Pos<typename MESH::FaceType> pi = pv[i];
			face::Pos<typename MESH::FaceType> pj = pv[j];
			face::Pos<typename MESH::FaceType> pk = pv[k];

			//test complex edge
			if(existEdge<MESH>(pi,pj) || existEdge<MESH>(pj,pk)|| existEdge<MESH>(pk,pi)	)	
			{
				return Weight();
			}
			// Return an infinite weight, if one of the neighboring patches
			// could not be created.
			if(v[i][j] == -1){return Weight();}
			if(v[j][k] == -1){return Weight();}

			//calcolo il massimo angolo diedrale, se esiste.
			float angle = 0.0f;
			face::Pos<typename MESH::FaceType> px;
			if(i + 1 == j)
			{
				px = pj; 
				px.FlipE(); px.FlipV();
				angle = std::max<float>(angle , ComputeDihedralAngle<MESH>(pi.v, pj.v, pk.v, px.v)	);
			}
			else
			{
				angle = std::max<float>( angle, ComputeDihedralAngle<MESH>(pi.v,pj.v, pk.v, pv[ v[i][j] ].v));
			}

			if(j + 1 == k)
			{
				px = pk; 
				px.FlipE(); px.FlipV();
				angle = std::max<float>(angle , ComputeDihedralAngle<MESH>(pj.v, pk.v, pi.v, px.v)	);
			}
			else
			{
				angle = std::max<float>( angle, ComputeDihedralAngle<MESH>(pj.v,pk.v, pi.v, pv[ v[j][k] ].v));
			}

			if( i == 0 && k == (int)v.size() - 1)
			{
				px = pi; 
				px.FlipE(); px.FlipV();
				angle = std::max<float>(angle , ComputeDihedralAngle<MESH>(pk.v, pi.v,  pj.v,px.v )	);
			}

			typename MESH::ScalarType area = ( (pj.v->P() - pi.v->P()) ^ (pk.v->P() - pi.v->P()) ).Norm() * 0.5;

			return Weight(angle, area);
		}

		template <class MESH>
			std::vector<typename MESH::VertexPointer > calculateMinimumWeightTriangulation(MESH &m, std::vector<face::Pos<typename MESH::FaceType> > vv )
		{
			std::vector< std::vector< Weight > > w; //matrice dei pesi minimali di ogni orecchio preso in conzideraione
			std::vector< std::vector< int    > > vi;//memorizza l'indice del terzo vertice del triangolo

			//hole size
			int nv = vv.size();

			w.clear();
			w.resize( nv, std::vector<Weight>( nv, Weight() ) );

			vi.resize( nv, std::vector<int>( nv, 0 ) );

			//inizializzo tutti i pesi possibili del buco
			for ( int i = 0; i < nv-1; ++i )
				w[i][i+1] = Weight( 0, 0 );

			//doppio ciclo for per calcolare di tutti i possibili triangoli i loro pesi.
			for ( int j = 2; j < nv; ++j )
			{
				for ( int i = 0; i + j < nv; ++i )
				{
					//per ogni triangolazione mi mantengo il minimo valore del peso tra i triangoli possibili
					Weight minval;

					//indice del vertice che da il peso minimo nella triangolazione corrente
					int minIndex = -1;

					//ciclo tra i vertici in mezzo a i due prefissati
					for ( int m = i + 1; m < i + j; ++m )
					{
						Weight a = w[i][m];
						Weight b = w[m][i+j];
						Weight newval =  a + b + computeWeight<MESH>( i, m, i+j, vv, vi);
						if ( newval < minval )
						{
							minval = newval;
							minIndex = m;
						}
					}
					w[i][i+j] = minval;
					vi[i][i+j] = minIndex;
				}
			}

			//Triangulate
			int i, j;
			i=0; j=nv-1;
			std::vector<typename MESH::VertexPointer > vf;

			vf.clear();

			triangulate<MESH>(vf, i, j, vi, vv);

			return vf;
		}

		template<class MESH>
			void triangulate(std::vector<typename MESH::VertexPointer > &m,int i, int j, std::vector< std::vector<int> > vi, 
			std::vector<face::Pos<typename MESH::FaceType> > vv)
		{
			if(i + 1 == j){return;}
			if(i==j)return;

			int k = vi[i][j];

			if(k == -1)	return;

			m.push_back(vv[i].v);
			m.push_back(vv[k].v);
			m.push_back(vv[j].v);

			triangulate<MESH>(m, i, k, vi, vv);
			triangulate<MESH>(m, k, j, vi, vv);
		}

		template <class MESH>
			void FillHoleMinimumWeight(MESH &m, bool Selected)
		{
			typename MESH::FaceIterator fi;
			std::vector<face::Pos<typename MESH::FaceType> > vvi;
			std::vector<typename MESH::FacePointer * > vfp;

			std::vector<typename tri::HoleInfo<MESH> > vinfo;
			typename std::vector<typename tri::HoleInfo<MESH> >::iterator VIT;
			int UBIT = GetHoleInfo<MESH>(m, Selected,vinfo);

			for(VIT = vinfo.begin(); VIT != vinfo.end();++VIT)
			{
				vvi.push_back(VIT->p);
			}

			typename std::vector<face::Pos<typename MESH::FaceType> >::iterator ith;
			typename std::vector<face::Pos<typename MESH::FaceType> >::iterator ithn;
			typename std::vector<typename MESH::VertexPointer >::iterator itf;

			std::vector<face::Pos<typename MESH::FaceType> > app;
			face::Pos<typename MESH::FaceType> ps;
			std::vector<typename MESH::FaceType > tr;
			std::vector<typename MESH::VertexPointer > vf;

			for(ith = vvi.begin(); ith!= vvi.end(); ++ith)
			{
				tr.clear();
				vf.clear();
				app.clear();
				vfp.clear();

				for(ithn = vvi.begin(); ithn!= vvi.end(); ++ithn)
					vfp.push_back(&(ithn->f));

				ps = *ith;
				getBoundHole<MESH>(ps,app);

				vf = calculateMinimumWeightTriangulation(m, app);

				if(vf.size() == 0)continue;//non e' stata trovata la triangolazione

				typename MESH::FaceIterator f = tri::Allocator<MESH>::AddFaces(m, app.size()-2, vfp);

				for(itf = vf.begin();itf != vf.end(); )
				{
					(*f).V(0) = (*itf++);
					(*f).V(1) = (*itf++);
					(*f).V(2) = (*itf++);
					++f;
				}
			}		
		}

		template <class MESH>
			void getBoundHole (face::Pos<typename MESH::FaceType> sp,std::vector<face::Pos<typename MESH::FaceType> >&ret)
		{
			face::Pos<typename MESH::FaceType> fp = sp;
			//take vertex around the hole
			do
			{
				assert(fp.IsBorder());
				ret.push_back(fp);
				fp.NextB();
			}while(sp != fp);
		}

	} // end namespace
}
#endif
