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
****************************************************************************/
#ifndef __VCG_TRI_UPDATE_HOLE
#define __VCG_TRI_UPDATE_HOLE

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
	MESH::VertexPointer v[2];
	SimpleEdge(MESH::VertexPointer v0, MESH::VertexPointer v1)
		{
			if(v0>v1) {v[0]=v1; v[1]=v0;}
			else {v[0]=v0; v[1]=v1;}
		}

	SimpleEdge(MESH::hedgepos_type &ep)		{
			*this=SimpleEdge(ep.VFlip(), ep.v);
	}
	
	bool operator < (const SimpleEdge & e) const
		{		return (v[0]!=e.v[0])?(v[0]<e.v[0]):(v[1]<e.v[1]); }
						
};

template<class MESH>
class HoleInfo
{
public: 
	HoleInfo(){}
	HoleInfo(MESH::hedgepos_type const &pHole, int  const pHoleSize, Box3<MESH::scalar_type> &pHoleBB)
		{
			p=pHole;	
			size=pHoleSize;
			bb=pHoleBB;
		}
	MESH::hedgepos_type p;
	int size;
	Box3<MESH::scalar_type>  bb;

	bool operator <  (const  HoleInfo & hh) const {return size <  hh.size;}
	bool operator >  (const  HoleInfo & hh) const {return size >  hh.size;}
	bool operator == (const  HoleInfo & hh) const {return size == hh.size;}
	bool operator != (const  HoleInfo & hh) const {return size != hh.size;}
	bool operator >= (const  HoleInfo & hh) const {return size >= hh.size;}
	bool operator <= (const  HoleInfo & hh) const {return size <= hh.size;}

	MESH::scalar_type Perimeter()
		{
			MESH::scalar_type sum=0;
			MESH::hedgepos_type ip = p;
						do
						{
							sum+=Distance(ip.v->cP(),ip.VFlip()->cP());
							ip.NextB();
						}
						while (ip != p);
						return sum;				
		}


	int CollectEdges(set< SimpleEdge<MESH> > &EV)
	{
		assert(p.IsBorder());
		EV.clear();
		int tsz=0;
		MESH::hedgepos_type ip=p;
		MESH::hedgepos_type tp;

		do
		{
			// Stesso codice della nextb
			do
			{
				ip.NextE();
				EV.insert(SimpleEdge<MESH>(ip)); // l'edge che sto scorrendo
				tp=ip;
				tp.FlipV();tp.FlipE();
				EV.insert(SimpleEdge<MESH>(tp)); // l'edge della faccia su cui sono e opposto al vertice su cui ruoto
				tp.FlipF(); tp.FlipE();
				EV.insert(SimpleEdge<MESH>(tp));  // gli altri due edge della faccia opposta a questa
				tp.FlipE();
				EV.insert(SimpleEdge<MESH>(tp));
			}
			while(!ip.f->IsBorder(ip.z));
			ip.FlipV();
			++tsz;
		}
		while (ip != p);
		assert(tsz==size);

		return EV.size();
	}
};

template<class MESH>
void FindHole(MESH &m, MESH::hedgepos_type ep, HoleInfo<MESH> &h)
{
	if(!ep.IsBorder()) return;

	int holesize = 0;
						
	Box3<MESH::scalar_type> hbox;
	if(ep.v->IsR()) hbox.Add(ep.v->cP());
	MESH::hedgepos_type init;
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
	MESH::face_iterator pf;
	int holesize;
	for (pf=m.face.begin(); pf!=m.face.end(); ++pf)
			if( !(*pf).IsD() && (*pf).IsW() )
				(*pf).ClearV();

		MESH::hedgepos_type ep;
		for (pf=m.face.begin(); pf!=m.face.end(); ++pf)
		{
			if( !(*pf).IsDeleted() && !(*pf).IsV() && (*pf).IsR() )
			{
				for(int j=0; j<3; ++j)
					if( (*pf).IsBorder(j) && !(*pf).IsV() && (*pf).IsR() )
					{
						(*pf).SetV();
						ep.Set(&*pf, j, (*pf).V(j));
						holesize = 0;
						
						Box3<MESH::scalar_type> hbox;
						if(ep.v->IsR()) hbox.Add(ep.v->cP());
						MESH::hedgepos_type init;
						init = ep;
						do
						{
							ep.NextB();
							ep.f->SetV();
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
	MSH_TYPE::hedgepos_type e0;	// 
	MSH_TYPE::hedgepos_type e1;	// 
	MSH_TYPE::scalar_type quality;
	TrivialEar(){}
	TrivialEar(const MSH_TYPE::hedgepos_type & ep)
	{
		e0=ep;
		assert(e0.IsBorder());
		e1=e0;
		e1.NextB();
		ComputeQuality();
	}

		// Nota: minori invertiti
	inline bool operator <  ( const TrivialEar & c ) const { return quality >  c.quality; }
	inline bool operator >  ( const TrivialEar & c ) const { return quality <  c.quality; }
	inline bool operator == ( const TrivialEar & c ) const { return quality == c.quality; }
	inline bool operator != ( const TrivialEar & c ) const { return quality != c.quality; }
	inline bool operator >= ( const TrivialEar & c ) const { return quality <= c.quality; }
	inline bool operator <= ( const TrivialEar & c ) const { return quality >= c.quality; }

	bool IsNull(){return e0.IsNull() || e1.IsNull();}
	void SetNull(){e0.SetNull();e1.SetNull();}
	void ComputeQuality(){ quality = Distance(e0.VFlip()->P(),e1.v->P());};
	bool IsUpToDate()	{return (e0.IsBorder() && e1.IsBorder());};

	bool Degen()
		{
			MSH_TYPE::hedgepos_type	ep=e0; ep.FlipV(); ep.NextB(); ep.FlipV(); // he precedente a e0 
			MSH_TYPE::hedgepos_type	en=e1; en.NextB();												 // he successivo a e1

			// caso ear degenere per buco triangolare
			if(ep==en) return true; 
			// Caso ear non manifold a
			if(ep.v==en.v)	return true;
			// Caso ear non manifold b
			if(ep.VFlip()==e1.v) return true;

			return false;
		}

	bool Close(TrivialEar &ne0, TrivialEar &ne1, MSH_TYPE::face_type* f)
	{
		// simple topological check
		if(e0.f==e1.f) {
			TRACE("Avoided bad ear");
			return false;
		}

		//usato per generare una delle due nuove orecchie.
		MSH_TYPE::hedgepos_type	ep=e0; ep.FlipV(); ep.NextB(); ep.FlipV(); // he precedente a e0 
		MSH_TYPE::hedgepos_type	en=e1; en.NextB();												 // he successivo a e1
		
		(*f).V(0) = e0.VFlip();
		(*f).V(1) = e0.v;
		(*f).V(2) = e1.v;

		(*f).F(0) = e0.f;
		(*f).Z(0) = e0.z;
		(*f).F(1) = e1.f;
		(*f).Z(1) = e1.z;
		(*f).F(2) = f;
		(*f).Z(2) = 2;

		e0.f->F(e0.z)=f;
		e0.f->Z(e0.z)=0;	
		
		e1.f->F(e1.z)=f;
		e1.f->Z(e1.z)=1;	

		// caso ear degenere per buco triangolare
		if(ep==en)
		{
			TRACE("Closing the last triangle");
			f->F(2)=en.f;
			f->Z(2)=en.z;
			en.f->F(en.z)=f;
			en.f->Z(en.z)=2;
			ne0.SetNull();
			ne1.SetNull();
		}
	  // Caso ear non manifold a
		else if(ep.v==en.v)
		{
			TRACE("Ear Non manif A\n");
			MSH_TYPE::hedgepos_type	enold=en;
			en.NextB();
			f->F(2)=enold.f;
			f->Z(2)=enold.z;
			enold.f->F(enold.z)=f;
			enold.f->Z(enold.z)=2;
			ne0=TrivialEar(ep);
			ne1=TrivialEar(en);
		}
		// Caso ear non manifold b
		else if(ep.VFlip()==e1.v)
		{
			TRACE("Ear Non manif B\n");
			MSH_TYPE::hedgepos_type	epold=ep; 
			ep.FlipV(); ep.NextB(); ep.FlipV();
			f->F(2)=epold.f;
			f->Z(2)=epold.z;
			epold.f->F(epold.z)=f;
			epold.f->Z(epold.z)=2;
			ne0=TrivialEar(ep);
			ne1=TrivialEar(en);
		}
		else // caso standard
		// Now compute the new ears;
		{
			ne0=TrivialEar(ep);
			ne1=TrivialEar(MSH_TYPE::hedgepos_type(f,2,e1.v));
		}

		return true;
	}
};




// Funzione principale per chiudier un buco in maniera topologicamente corretta.
// Gestisce situazioni non manifold ragionevoli 
// (tutte eccetto quelle piu' di 2 facce per 1 edge).
// Controlla che non si generino nuove situazioni non manifold chiudendo orecchie
// che sottendono un edge che gia'esiste.
//
// Attenzione: se per riaggiungere facce deve riallocare il vettore non funge!!!!
// 
template<class MESH, class EAR>
MESH::face_iterator CloseHole(MESH &m, HoleInfo<MESH> &h)
{

	set<SimpleEdge<MESH> > ES;  // vettore con tutti gli edge adiacenti al buco.
	h.CollectEdges(ES);
	vector<EAR> H;			// Heap delle ear da chiudere
	H.reserve(h.size);
	
	assert(h.p.IsBorder());
	MESH::hedgepos_type ep=h.p;
	do {
		H.push_back(EAR(ep));
		ep.NextB();
	}	while(ep!=h.p);

	make_heap(H.begin(),H.end());
	int cnt=h.size;
	EAR en0,en1;
	MESH::face_iterator f=m.AddFaces(h.size-2);
	
	MESH::face_iterator firstf=f;
  SimpleEdge<MESH> se(0,0);
	while(cnt>2 && !H.empty())
	{
		pop_heap(H.begin(),H.end());
		se=SimpleEdge<MESH>(H.back().e0.VFlip(), H.back().e1.v);
		if(H.back().IsUpToDate())	
		{
			if(!H.back().Degen() && ES.find(se)!=ES.end()){ 
				// Nota che nel caso di ear degeneri si DEVE permettere la creazione di un edge che gia'esiste
				TRACE("Evitata orecchia brutta!");
			}
			else if(H.back().Close(en0,en1,&*f))
			{
				ES.insert(se);
				if(!en0.IsNull()){
					H.push_back(en0);
					push_heap( H.begin(), H.end());
				}
				if(!en1.IsNull()){
					H.push_back(en1);
					push_heap( H.begin(), H.end());
				}
				--cnt;
				++f;
				//return firstf;///////////////dbug
			}
		}
		H.pop_back();
	}
	//Delete the unused faces (caused by non 1-manifold vertexes)
	while(f!=m.face.end())
		{
			(*f).SetD();
			++f;
			m.fn--;
		}
	return firstf;
};

/*
Trivial Ear con preferred Normal
*/
template<class MSH_TYPE> class TrivialEarN : public TrivialEar<MSH_TYPE>
{
	public:

		TrivialEarN(){}
	TrivialEarN(const MSH_TYPE::hedgepos_type & ep)
	{
		e0=ep;
		assert(e0.IsBorder());
		e1=e0;
		e1.NextB();
		ComputeQuality();
	}


	static MSH_TYPE::vectorial_type &PreferredNormal()
	{
			static MSH_TYPE::vectorial_type nn;
			return nn;
	}

	void ComputeQuality(){ 
		Point3d nn= -Normal(	e0.VFlip()->P(), e0.v->P(), e1.v->P());
		quality = Distance(e0.VFlip()->P(),e1.v->P());
		if(nn*PreferredNormal() < -0.1) 
			quality*=1000000;
		
	};

};

/* 2d Triangulation Code */ 
class Triangulate2D 
{

static double Area(const vector<Point2d> &contour)
{
  int n = contour.size();

  double A=0.0f;

  for(int p=n-1,q=0; q<n; p=q++) {
    A+= contour[p].x()*contour[q].y() - contour[q].x()*contour[p].y();
  }
  return A*0.5f;
}

   /*
     InsideTriangle decides if a point P is Inside of the triangle
     defined by A, B, C.
   */
static bool InsideTriangle(double Ax, double Ay,
                      double Bx, double By,
                      double Cx, double Cy,
                      double Px, double Py)

{
  double ax, ay, bx, by, cx, cy, apx, apy, bpx, bpy, cpx, cpy;
  double cCROSSap, bCROSScp, aCROSSbp;

  ax = Cx - Bx;  ay = Cy - By;
  bx = Ax - Cx;  by = Ay - Cy;
  cx = Bx - Ax;  cy = By - Ay;
  apx= Px - Ax;  apy= Py - Ay;
  bpx= Px - Bx;  bpy= Py - By;
  cpx= Px - Cx;  cpy= Py - Cy;

  aCROSSbp = ax*bpy - ay*bpx;
  cCROSSap = cx*apy - cy*apx;
  bCROSScp = bx*cpy - by*cpx;

  return ((aCROSSbp >= 0.0f) && (bCROSScp >= 0.0f) && (cCROSSap >= 0.0f));
};

static bool Snip(const vector<Point2d> &contour,int u,int v,int w,int n,int *V)
{
  int p;
  double Ax, Ay, Bx, By, Cx, Cy, Px, Py;
	const double epsilon =1e-2;

  Ax = contour[V[u]].x();
  Ay = contour[V[u]].y();

  Bx = contour[V[v]].x();
  By = contour[V[v]].y();

  Cx = contour[V[w]].x();
  Cy = contour[V[w]].y();

  if ( epsilon> (((Bx-Ax)*(Cy-Ay)) - ((By-Ay)*(Cx-Ax))) ) return false;

  for (p=0;p<n;p++)
  {
    if( (p == u) || (p == v) || (p == w) ) continue;
    Px = contour[V[p]].x();
    Py = contour[V[p]].y();
    if (InsideTriangle(Ax,Ay,Bx,By,Cx,Cy,Px,Py)) return false;
  }

  return true;
}
public:
static bool Process(const vector<Point2d> &contour,vector<int> &result)
{
  /* allocate and initialize list of Vertices in polygon */

  int n = contour.size();
	double area=Area(contour);
  if ( n < 3 ) return false;

  int *V = new int[n];

  /* we want a counter-clockwise polygon in V */

  if ( 0.0f < area )   for (int v=0; v<n; v++) V[v] = v;
	else{
			    for(int v=0; v<n; v++) V[v] = (n-1)-v;
					area=-area;
	}


  int nv = n;

  /*  remove nv-2 Vertices, creating 1 triangle every time */
  int count = 2*nv;   /* error detection */

	double CurrBest=Sqrt(area)/1000;

  for(int m=0, v=nv-1; nv>2; )
  {
		count--;
    /* if we loop, it is probably a non-simple polygon */
    if( count<0)
    {
			CurrBest*=1.3;
			count = 2*nv;
			
			if(CurrBest>Sqrt(area)*2)
				return false;
    }

    /* three consecutive vertices in current polygon, <u,v,w> */
    int u = v  ; if (nv <= u) u = 0;     /* previous */
    v = u+1; if (nv <= v) v = 0;     /* new v    */
    int w = v+1; if (nv <= w) w = 0;     /* next     */

		if(Distance(contour[u],contour[w]) < CurrBest)
				if ( Snip(contour,u,v,w,nv,V) )
				{
					int a,b,c,s,t;

					/* true names of the vertices */
					a = V[u]; b = V[v]; c = V[w];

					/* output Triangle */
					result.push_back( a );
					result.push_back( b );
					result.push_back( c );

					m++;

					/* remove v from remaining polygon */
					for(s=v,t=v+1;t<nv;s++,t++) V[s] = V[t];
					nv--;

					/* resest error detection counter */
					count = 2*nv;
				}
  }

  delete V;

  return true;
}

};
} // end namespace
#endif