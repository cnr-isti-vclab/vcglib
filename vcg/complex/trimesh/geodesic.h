/*#***************************************************************************
 * Geodesic.h                                                         o o    *
 *                                                                  o     o  *
 * Visual Computing Group                                           _  O  _  *
 * IEI Institute, CNUCE Institute, CNR Pisa                          \/)\/   *
 *                                                                  /\/|     *
 * Copyright(C) 1999 by Paolo Cignoni,                                 |     *
 * All rights reserved.                                                \     *
 *                                                                           *
 * Permission  to use, copy, modify, distribute  and sell this  software and *
 * its documentation for any purpose is hereby granted without fee, provided *
 * that  the above copyright notice appear  in all copies and that both that *
 * copyright   notice  and  this  permission  notice  appear  in  supporting *
 * documentation. the author makes  no representations about the suitability *
 * of this software for any purpose. It is provided  "as is" without express *
 * or implied warranty.                                                      *
 *                                                                           *
 ***************************************************************************#*/
/*#**************************************************************************
  History

 2002 Apr 01 First Working release copiata, ripulita e templatata dal plyvcg(pc)
	  Dic 27 Aggiunta classe Geo. First release (gano)
 2003 Jan 10 Aggiunto controllo IsD() in ComputeGeodesicQuality() (pc)
	  Feb 02 Cambiato l'algoritmo e aggiunta FindPath 
	  May 05 Corretti un po' di bug sui casi limite e templatata la BuildSP
			 Aggiunta SelectRegion
	  May 13 Variazioni alla SelectRegion per la selezione di triangoli
	         in casi limite (longobardi)
	  May 19 Corretto baco in BuildSP, aggiunto template, ripulito un po' (gano)
		July 3 *** aggiornati i tipi di meshPos alla nuova versione
					 aggiunta SelectRegionFF (gano)
		    16 Aggiunto namespace
		Nov 14 Corretto bug nella ComputeGeodesicQuality. (pc)
 *#**************************************************************************/

#ifndef __VCG_GEODESIC__
#define __VCG_GEODESIC__

namespace vcg {
	
template <class MESH>
class VQualityHeap
{
public:
	float q;
	MESH::vertex_pointer p;
	inline VQualityHeap( MESH::vertex_pointer np )
	{
		q = np->Q();
		p = np;
	}
		// Attenzione il minore e' maggiore
	inline bool operator <  ( const VQualityHeap & vq ) const { return q >  vq.q; }
	inline bool operator == ( const VQualityHeap & vq ) const { return q == vq.q; }
	inline bool operator >  ( const VQualityHeap & vq ) const { return q <  vq.q; }
	inline bool operator != ( const VQualityHeap & vq ) const { return q != vq.q; }
	inline bool operator <= ( const VQualityHeap & vq ) const { return q >= vq.q; }
	inline bool operator >= ( const VQualityHeap & vq ) const { return q <= vq.q; }
	inline bool is_valid() const { return q==p->Q(); }
};

// Calcola la qualita' come distanza geodesica dal bordo della mesh.
// Robusta funziona anche per mesh non manifold.
// La qualita' memorizzata indica la distanza assoluta dal bordo della mesh.
// Nota prima del 13/11/03 in alcuni casi rari SPT andava in loop perche' poteva capitare
// che per approx numeriche ben strane pw->Q() > pv->Q()+d ma durante la memorizzazione 
// della nuova distanza essa rimanesse uguale a prima. Patchato rimettendo i vertici nello 
// heap solo se migliorano la distanza di un epsilon == 1/100000 della mesh diag.

template <class MESH>
void ComputeGeodesicQuality(MESH &m, bool per_face )	// R1
{
	//Requirements
	assert(m.HasVFTopology());
	assert(m.HasPerVertexQuality());
  if(per_face) assert(m.HasPerFaceQuality());

	vector<VQualityHeap<MESH> > heap;
	MESH::vertex_iterator v;
	MESH::face_iterator   f;
	int j;

	m.VFTopology();									// Creazione adiacenza vertici
	m.ComputeBorderFlagVF();				// Marco gli edge di bordo

	for(v=m.vert.begin();v!=m.vert.end();++v)
		(*v).Q() = -1;
	for(f=m.face.begin();f!=m.face.end();++f)			// Inserisco nell'heap i v di bordo
		if(!(*f).IsD())
			for(j=0;j<3;++j)
				if( (*f).IsB(j) )
				{
					for(int k=0;k<2;++k)
					{
						MESH::vertex_pointer pv = (*f).V((j+k)%3);
						if( pv->Q()==-1 )
						{
							pv->Q() = 0;
							heap.push_back(VQualityHeap<MESH>(pv));
						}
					}
				}
	
 const MESH::scalar_type loc_eps=m.bbox.Diag()/MESH::scalar_type(100000);
 while( heap.size()!=0 )							// Shortest path tree
	{
		MESH::vertex_pointer pv;
		pop_heap(heap.begin(),heap.end());
		if( ! heap.back().is_valid() )
		{
			heap.pop_back();
			continue;
		}
		pv = heap.back().p;
		heap.pop_back();
	 MESH::vedgepos_type x;
		for( x.f = pv->Fp(), x.z = pv->Zp(); x.f!=0; x.NextF() )
		{
			for(int k=0;k<2;++k)
			{
				MESH::vertex_pointer pw;
				float d;
				if(k==0) pw = x.f->V1(x.z);
				else     pw = x.f->V2(x.z);
				d = Distance(pv->P(),pw->P());
				if( pw->Q()==-1 || pw->Q() > pv->Q()+d + loc_eps)
				{
					pw->Q() = pv->Q()+d;
					heap.push_back(VQualityHeap<MESH>(pw));
					push_heap(heap.begin(),heap.end());
				}
			}
		}
	}

	for(v=m.vert.begin();v!=m.vert.end();++v)
		if(v->Q()==-1)
			v->Q() = 0;

	if(per_face)
	{
		for(f=m.face.begin();f!=m.face.end();++f)
			(*f).Q() = ((*f).V(0)->Q() + (*f).V(1)->Q() + (*f).V(2)->Q())/3;
	}
}





// ************************************** classe GEo
// first release
// Modo d'uso
//		Geo<TYPE_MESH> geo(m); //dichiare e definisce
//			m.VFTopology();
//	    geo.FartestPoint(....);//dato un insieme di punti a distanza 0 trova
//							   // il punto piu'lontano sulla superficie (una versione ritorna anche il 
//							   // cammino)
//		
//		geo.FindPath(...);	// determina il cammino minimo tra due vertici	
//	m.ComputeVBorderFlag(); (SOLO SE SI USA MostInternal)
//		geo.MostInternal(....);// punto piu'lontano dal bordo: invoca FartestPoint passando
//							   // tutti i vertici di bordo come insieme di partenza.
//							   // NB: questa fa la stessa cosa di ComputeGeodesicQuality
//							   // ma non scrive sul campo qualita' del vertice
//	
//
// *****************************************

template <class TMESH>
class Geo{

	public:

		Geo( TMESH & m):TD(m.vert),TDP(m.vert),TDS(m.vert),TDF(m.face),M(m)
			{lessThan.pt = this;}
		~Geo(){}

	TMESH &M;	
	template <class TMESH>
	struct TempData{
		TempData(){}
		TempData(const double & d_){d=d_;visited=false;}
		double d;
		bool visited;
		};

	template <class TMESH>
	struct TempDataPath{
		TempDataPath(){}
		TempDataPath(const TMESH::vertex_pointer & parent_){parent = parent_;}
		TMESH::vertex_pointer  parent;		
		};

	template <class TMESH>
	struct TempDataSource{
		TempDataSource(){}
		TempDataSource(const TMESH::vertex_pointer & source_){source = source_;}
		TMESH::vertex_pointer  source;		
		};

	TempData_vect<TMESH::vertex_type, TempData<TMESH> >  TD;
	TempData_vect<TMESH::vertex_type, TempDataPath<TMESH> > TDP;
	TempData_vect<TMESH::vertex_type, TempDataSource<TMESH> > TDS;

	
	template<class VERTEX_POINTER, class PT> 
	struct pred: public binary_function<VERTEX_POINTER,VERTEX_POINTER,bool>{
			PT * pt;
			bool operator()(const VERTEX_POINTER& v0, const VERTEX_POINTER& v1) const
				{return (pt->TD[v0].d > pt->TD[v1].d);}
		};

	 pred<TMESH::vertex_pointer,Geo<TMESH> > lessThan;
	
	 template<class VERTEX_POINTER>
	 struct VPath{
		 VPath(VERTEX_POINTER v0,VERTEX_POINTER v1, double dist):pathLenght(dist){
			 sources[0] = v0;sources[1] = v1;
			 }
		VERTEX_POINTER sources[2];
		double pathLenght;
		 };

	 typedef typename VPath<TMESH::vertex_pointer> Path;

	 template <bool RETURN_PATH, bool CONNECT_TWO_SOURCES,bool FARTEST_ON_BORDER,bool SKIP_SELECTED>
	 TMESH::vertex_pointer BuildSP( 
		 deque<TMESH::vertex_pointer> &path,
		 vector<TMESH::vertex_pointer> & _frontier,
		 double & max_distance
		 )
		 {
		 TMESH::vertex_pointer curr=NULL,fartest=NULL,pw1;	
		 bool isLeaf;	
		 Path shortestPath(NULL,NULL,-1.0);
		 TMESH::vertex_pointer sources[2];
		 assert(!CONNECT_TWO_SOURCES || (_frontier.size()==2));

		 if(CONNECT_TWO_SOURCES){
			 if( _frontier.size() != 2) return (TMESH::vertex_pointer)0;
			 if( _frontier[0] == _frontier[1]) return (TMESH::vertex_pointer)0;
			 sources[0] = _frontier[0];
			 sources[1] = _frontier[1];
			 }
	
		 vector<TMESH::vertex_pointer> frontier;
		 frontier.clear();

		 //Requirements
		 assert(M.HasVFTopology());

		 if(M.vn==0) return NULL;

		 list<TMESH::vertex_pointer> children;
		 bool toQueue;

		 TD.Start(TempData<TMESH>(-1.0));

		 if(RETURN_PATH)
			 TDP.Start(TempDataPath<TMESH>(NULL));

		 list<TMESH::vertex_pointer>::iterator is;
		 deque<TMESH::vertex_pointer> leaves;
		 vector <pair<TMESH::vertex_pointer,TMESH::scalar_type> > expansion;

		 vector <TMESH::vertex_pointer >::const_iterator ifr;
		 TMESH::vedgepos_type x;int k;
		 TMESH::vertex_pointer pw;

		 for(ifr = _frontier.begin(); ifr != _frontier.end(); ++ifr){
			 TD[*ifr].visited= true;
			 TD[*ifr].d = 0.0;	
			 }
		if(CONNECT_TWO_SOURCES){
			TDS.Start(TempDataSource<TMESH>(NULL));
			for(ifr = _frontier.begin(); ifr != _frontier.end(); ++ifr)
				TDS[*ifr].source= (*ifr);
			}

 for(ifr = _frontier.begin(); ifr != _frontier.end(); ++ifr)
	 {	
	 if(RETURN_PATH)
		 TDP[*ifr].parent=NULL;
	 // determina la distanza dei vertici della fan

	 for( x.f = (*ifr)->Fp(), x.z = (*ifr)->Zp(); x.f!=0; x.NextF() )
		 for(k=0;k<2;++k)
			 {
			 if(k==0) pw = x.f->V1(x.z);
			 else     pw = x.f->V2(x.z);
	


			if(SKIP_SELECTED)
			if(pw->IsS()) 
				continue;

			 if(CONNECT_TWO_SOURCES)
					 {
					 if( (TDS[pw].source != NULL)&&(TDS[pw].source!= TDS[*ifr].source)){
						 double l = TD[pw].d + Distance((*ifr)->cP(),pw->cP());
						 if( (shortestPath.pathLenght== -1)||(l < shortestPath.pathLenght))
								 shortestPath = Path((*ifr),pw,l );
						 }
					 }

			 if(TD[pw].d ==-1){
				 TD[pw].d = Distance(pw->cP(),(*ifr)->cP());
				 frontier.push_back(pw);				

				 if(RETURN_PATH)
					 if(!TD[pw].visited)
						 TDP[pw].parent= (*ifr);

				if(CONNECT_TWO_SOURCES)
					TDS[pw].source  = (*ifr);

				 }
			 }
	 }

 if(CONNECT_TWO_SOURCES)
	if(shortestPath.pathLenght != -1)
		goto found;


 double curr_d,d_curr = 0.0;
 max_distance=0.0;

 make_heap(frontier.begin(),frontier.end(),lessThan);
 while(!frontier.empty()){
	 expansion.clear();
	 pop_heap(frontier.begin(),frontier.end(),lessThan);
	 curr = frontier.back();
	 frontier.pop_back();
	 d_curr = TD[curr].d;
	 TD[curr].visited=true;
	
	 if(CONNECT_TWO_SOURCES)
		 {
		 double otherDist;

		 if(shortestPath.pathLenght!= -1)
			 {
				if( TDS[curr].source == TDS[shortestPath.sources[0]].source)
					otherDist = TD[shortestPath.sources[1]].d;
				else
					otherDist = TD[shortestPath.sources[0]].d;

			 if (d_curr > shortestPath.pathLenght-otherDist)
				continue;
			 }
		 }
	 isLeaf = (!FARTEST_ON_BORDER || curr->IsB());

	 TMESH::vedgepos_type x;int k;

	 for( x.f = curr->Fp(), x.z = curr->Zp(); x.f!=0; x.NextF() )
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
		 const double & d_pw1  = TD[pw1].d;

		 if(SKIP_SELECTED)
			 if(pw->IsS()) 
				 continue;

		 if((TD[pw].visited &&
			(!CONNECT_TWO_SOURCES || (TDS[pw].source == TDS[curr].source))))
			 continue;

		 if(CONNECT_TWO_SOURCES){	
			 if( (TDS[pw].source!=TDS[curr].source) && (TDS[pw].source!=NULL)&& (TDS[curr].source!=NULL)){
					double l = TD[curr].d+ TD[pw].d + Distance(curr->cP(),pw->cP());
					if( (shortestPath.pathLenght == -1) || (shortestPath.pathLenght > l) )
					shortestPath = Path(curr,pw,l);
					continue;
				 }
		 if( (TDS[pw1].source!=TDS[curr].source) && (TDS[pw1].source!=NULL)&& (TDS[curr].source!=NULL)){
				double l = TD[curr].d+ TD[pw1].d + Distance(curr->cP(),pw1->cP());
				if( (shortestPath.pathLenght == -1) || (shortestPath.pathLenght > l) )
					shortestPath  =  Path(curr,pw1,l);

				continue;
				}
			 }
	 if((!TD[pw1].visited ) ||(d_pw1 == 0.0) )
		 continue;		

#ifdef _DEBUG
	 if(CONNECT_TWO_SOURCES){
		 assert(TDS[pw1].source  == TDS[curr].source);
		 assert((TDP[curr].parent == NULL) || TDS[curr].source == TDS[TDP[curr].parent].source);
		 }
		 assert( TD[pw1].d != -1);
		 assert( TD[pw1].d != 0.0);
		 assert( (curr!=pw) && (pw!=pw1) && (pw1 != curr));				
		 assert(d_pw1!=-1.0);
#endif

		 //************** calcolo della distanza di pw in base alle distanze note di pw1 e curr
		 //************** sapendo che (curr,pw,pw1) e'una faccia della mesh
		 //************** (vedi figura in file distance.gif)
		 Point3<TMESH::scalar_type> w_c = pw->cP()- curr->cP();
		 Point3<TMESH::scalar_type> w_w1 = pw->cP()- pw1->cP();
		 Point3<TMESH::scalar_type> w1_c = pw1->cP()- curr->cP();

		 double ew_c  = (w_c).Norm();
		 double ew_w1 = (w_w1).Norm();
		 double ec_w1 = (w1_c).Norm();
		 double	alpha,alpha_, beta,beta_,theta,h,delta,s,a,b;

		 alpha = acos((w_c*w1_c)/(ew_c*ec_w1));
		 s = (d_curr + d_pw1+ec_w1)/2;
		 a = s/ec_w1;
		 b = a*s;
		 alpha_ = 2*acos ( min(1.0,sqrt(  (b- a* d_pw1)/d_curr)));

		 if ( alpha+alpha_ > M_PI){
			 curr_d = d_curr + ew_c;		
			 }else
			 {
			 beta_ = 2*acos ( min(1.0,sqrt(  (b- a* d_curr)/d_pw1)));
			 beta  = acos((w_w1)*(-w1_c)/(ew_w1*ec_w1));

			 if ( beta+beta_ > M_PI){
				 curr_d = d_pw1  + ew_w1;
				 }
			 else 
				 {
				 theta	= M_PI-alpha-alpha_;
				 delta	= cos(theta)* ew_c;
				 h		= sin(theta)* ew_c;
				 curr_d = sqrt( pow(h,2)+ pow(d_curr + delta,2));
				 }
			 }
		 //**************************************************************************************

		 toQueue = (TD[(pw)].d==-1);

		 if(toQueue){// se non e'gia' in coda ce lo mette
			 if(RETURN_PATH)//aggiorna il puntatore al genitore se richiesto 
				 {
				 if(( TD[curr].d + ew_c) < (d_pw1  + ew_w1))
					 TDP[(pw)].parent = curr;
				 else
					 TDP[(pw)].parent = pw1;
				 assert(curr_d!=0.0);
				 }

			 if(CONNECT_TWO_SOURCES)
				 TDS[pw].source = TDS[TDP[pw].parent].source;

			 expansion.push_back(pair<TMESH::vertex_pointer,double>(pw,curr_d));
			 isLeaf =false;
			 }else{
				 if(  TD[pw].d > curr_d )
					 {
					 if(RETURN_PATH)
						 if(( d_curr + ew_c) < (d_pw1  + ew_w1))
							 TDP[(pw)].parent = curr;
						 else
							 TDP[(pw)].parent = pw1;

					TD[(pw)].d = curr_d;
					 if(CONNECT_TWO_SOURCES)
						 TDS[pw].source = TDS[curr].source;
				
						assert(!((TD[pw].visited && (TDS[pw].source != TDS[curr].source))));

					 isLeaf =false;
					 }
				 }
			 }

		 if(isLeaf){
			 if(d_curr > max_distance){
				 max_distance = d_curr;
				 fartest = curr;
				 }
			 }

	 vector <pair<TMESH::vertex_pointer,double> > ::iterator i;
	
	 for(i = expansion.begin(); i!= expansion.end(); ++i){
		 if(TD[(*i).first].d== -1.0){
			 TD[(*i).first].d = (*i).second;
			 frontier.push_back((*i).first);
			 push_heap(frontier.begin(),frontier.end(),lessThan);
			 }
		 }
	 }// end while: la frontiera e'vuota 


	 if(CONNECT_TWO_SOURCES)
		if(shortestPath.pathLenght == -1)// i path non si sono incontrati
			return NULL;

found:
// scrivi le distanze sul campo qualita' (nn: farlo parametrico)
 TMESH::vertex_iterator vi;
 for(vi = M.vert.begin(); vi != M.vert.end(); ++vi)
	 (*vi).Q() = TD[&(*vi)].d; 
 
 if(RETURN_PATH)
	 {
	 TMESH::vertex_pointer s; 
	 if(CONNECT_TWO_SOURCES)
		 {
			curr  = shortestPath.sources[0];
			pw  = shortestPath.sources[1];
		
			if(TDS[curr].source != sources[0]) swap(curr,pw);

			// curr e' il vertice estremo di una sorgente
			// pw quello dell'altra
			max_distance = 	TD[curr].d+ TD[pw].d + Distance(curr->cP(),pw->cP());	
			s = TDS[pw].source;
			while(pw!=NULL){// pw->Q() = 2.0 *max_distance;
				assert(TDS[pw].source == s);
				path.push_front(pw);
				pw = TDP[pw].parent;
				}	
			s=TDS[curr].source;
			}
			while(curr!=NULL){
			//curr->Q() = max_distance;		
					assert(!CONNECT_TWO_SOURCES || (TDS[curr].source == s));
				path.push_back(curr);
				curr = TDP[curr].parent;
				}		
	 }

 TD.Stop();	
 TDP.Stop();
 TDS.Stop();
 return  CONNECT_TWO_SOURCES? 0 : fartest;
}

// SelectRegion: dato un vertice esegue una visita di  superficie.
// Un ramo della visita si ferma quando incontra un vertice marcato S(elected)
// Serve per selecionare delle regioni delimitate da liste di vertici gia' definite
template <class STL_CONT_VERT,class STL_CONT_FACE, bool LIMITED, bool USE_FLAG_V>
bool SelectRegion(int iter,STL_CONT_VERT & frontier, STL_CONT_FACE & result ) {
	//result.clear();
	int i =0 ;
	TMESH::vertex_pointer curr=NULL;	
	if (!USE_FLAG_V)
		TD.Start(TempData<TMESH>(-1.0));

	//Requirements
	assert(M.HasVFTopology());

	if(M.vn==0) return false;


	bool toQueue;
	TMESH::EdgePos x;int k;
	TMESH::vertex_pointer pw;

	while((!frontier.empty())&&( LIMITED?(i<iter):true))
	{	if(LIMITED)
			i++;
		curr = frontier.back();	
		frontier.pop_back();
	
		if (USE_FLAG_V)
			curr->SetV();
		else
			TD[curr].visited = true;

		for( x.f = curr->Fp(), x.z = curr->Zp(); x.f!=0; x.NextF() )
			for(k=0;k<2;++k)
			{
			if(k==0) pw = x.f->V1(x.z);
			else     pw = x.f->V2(x.z);

			if(!(x.f)->IsS())  
				result.push_back(x.f);

			if(!pw->IsV() && !pw->IsS())
				frontier.push_back(pw);
			}
	}

	if (!USE_FLAG_V)
		TD.Stop();

	return frontier.empty();

}

// versione temporanea: data una faccia e un vertice seleziona la regione 
// procedendo per adiacenze faccia-faccia

struct TempDataFace{
		TempDataFace():front(0){}
		TempDataFace(const short int & _front):front(_front){}
		short int front;
		bool IsV(const int &i){
			return (front & (1<<i));
			}
		bool SetV(const int &i){
			return (front|= (1<<i));
			}
		};


TempData_vect<TMESH::face_type, TempDataFace >  TDF;

// versione rozza: non fa region growing, se si limita il numero di iterazioni
// puo'terminare con un insieme di facce a genus diverso da zero
template <class STL_CONT_VERT,class STL_CONT_FACE, bool LIMITED,bool USE_FLAG_V>
bool SelectRegionFF(int iter,STL_CONT_FACE & frontier, STL_CONT_FACE & result ) {
	//result.clear();
	int i =0, h=0 ;
	TDF.Start(TempDataFace(0));

	TMESH::face_pointer  curr=  NULL;	

	if (!USE_FLAG_V)
		TD.Start(TempData<TMESH>(-1.0));
	//Requirements
	assert(M.HasFFTopology());

	if(M.vn==0) return false;

	bool toQueue;
	TMESH::vedgepos_type x;int k;
	TMESH::vertex_pointer pw;

	while((!frontier.empty())&&( LIMITED?(i<iter):true))
	{	if(LIMITED)
			i++;
		curr = frontier.back();	
		frontier.pop_back();
	
		curr->SetV();
		for(h = 0; h<3;++h){
		TMESH::face_type* iii =(TMESH::face_type* )curr->F(h);
		TempDataFace & tdf = TDF[iii];
		if(	!tdf.IsV(curr->Z(h)) && (curr->F(h) != curr)){
			if(!(curr->F(h))->IsS())  {
					result.push_back((TMESH::face_type* )curr->F(h));
					frontier.push_back((TMESH::face_type* )curr->F(h));
				}
			tdf.SetV(curr->Z(h));
			}
			}
	}
			
	TDF.Stop();
	return frontier.empty();
	}


// multiple sorgenti -  no max iter
template <class STL_CONT_VERT,class STL_CONT_FACE, bool USE_FLAG_V>
bool SelectRegion(STL_CONT_VERT & frontier, STL_CONT_FACE & result){
		return SelectRegion<STL_CONT_VERT,STL_CONT_FACE, false,USE_FLAG_V>
			(STL_CONT_VERT & frontier, STL_CONT_FACE & result );
	}

// multiple sorgenti -  max iter
template <class STL_CONT_VERT,class STL_CONT_FACE, bool USE_FLAG_V>
bool SelectRegion(STL_CONT_VERT & frontier, STL_CONT_FACE & result, int iter){
		return SelectRegion<STL_CONT_VERT,STL_CONT_FACE, false,USE_FLAG_V>
			(iter,STL_CONT_VERT & frontier, STL_CONT_FACE & result);
	}

// singola sorgente -  max iter
template <class STL_CONT_VERT,class STL_CONT_FACE, bool USE_FLAG_V>
bool SelectRegion(const TMESH::vertex_pointer & v0, STL_CONT_FACE & result, int iter) {
	
	STL_CONT_VERT frontier;
	frontier.push_back(v0);
	return SelectRegion<STL_CONT_VERT,STL_CONT_FACE,true,USE_FLAG_V>( iter, frontier,result);
	}

// singola sorgente - no max iter
template <class STL_CONT_VERT,class STL_CONT_FACE, bool USE_FLAG_V>
bool SelectRegion(const TMESH::vertex_pointer & v0, STL_CONT_FACE & result) {
	
	STL_CONT_VERT frontier;
	frontier.push_back(v0);
	return SelectRegion<STL_CONT_VERT,STL_CONT_FACE,false,USE_FLAG_V>( 0, frontier,result);
	}


template<class TMESH>
struct cp{
	cp(){};
	cp(TMESH::vertex_pointer v0_,TMESH::vertex_pointer v1_):v0(v0_),v1(v1_){};
	TMESH::vertex_pointer v0,v1;
	const bool operator != (const cp<TMESH> & other)const {
		return ( ( v0 != other.v0) ||  ( v1 != other.v1));
	}

	const bool operator < (const cp<TMESH> & other)const {
		return ( ( v0 == other.v0)? ( v1 < other.v1):v0 < other.v0);
	}
	};

set<cp<TMESH> > connected;
bool AreConn(TMESH::vertex_pointer v0,TMESH::vertex_pointer v1){
	return (connected.find(cp<TMESH>(v0,v1))!=connected.end());
	}
void Conn(TMESH::vertex_pointer v0,TMESH::vertex_pointer v1){
	connected.insert(cp<TMESH>(v0,v1));
	connected.insert(cp<TMESH>(v1,v0));
	}


void Delaunay( 
		 deque<TMESH::vertex_pointer> &path,
		 vector<TMESH::vertex_pointer> & _frontier)
	{
		TMESH::vertex_pointer s;
		 bool isLeaf,found=false;
		 vector<TMESH::vertex_pointer> frontier;
		 frontier.clear();
		 //Requirements
		 assert(M.HasVFTopology());
		 TMESH::vertex_iterator ii;
		 list<TMESH::vertex_pointer> children;
		 TMESH::vertex_pointer curr = NULL,fartest,pw1;	
		 bool toQueue;

		 TD.Start(TempData<TMESH>(-1.0));
		 TDP.Start(TempDataPath<TMESH>(NULL));

		 list<TMESH::vertex_pointer>::iterator is;
		 deque<TMESH::vertex_pointer> leaves;
		 vector <pair<TMESH::vertex_pointer,TMESH::scalar_type> > expansion;

		 vector <TMESH::vertex_pointer >::const_iterator ifr;
		 TMESH::EdgePos x;int k;
		 TMESH::vertex_pointer pw;

		 for(ifr = _frontier.begin(); ifr != _frontier.end(); ++ifr){
			 TD[*ifr].visited= true;
			 TD[*ifr].d = 0.0;	
			 }

		TDS.Start(TempDataSource<TMESH>(NULL));
		for(ifr = _frontier.begin(); ifr != _frontier.end(); ++ifr)
			TDS[*ifr].source= (*ifr);

 for(ifr = _frontier.begin(); ifr != _frontier.end(); ++ifr)
	 {	

		 TDP[*ifr].parent=NULL;
	 // determina la distanza dei vertici della fan
	 for( x.f = (*ifr)->Fp(), x.z = (*ifr)->Zp(); x.f!=0; x.NextF() )
		 for(k=0;k<2;++k)
			 {
			 if(k==0) pw = x.f->V1(x.z);
			 else     pw = x.f->V2(x.z);

			 if(TD[pw].d ==-1){
				 TD[pw].d = Distance(pw->cP(),(*ifr)->cP());
				 frontier.push_back(pw);				

				if(TD[pw].visited)
					TDP[pw].parent= (*ifr);
				 TDS[pw].source  = (*ifr);
				 }
			 }
	 }

 make_heap(frontier.begin(),frontier.end(),lessThan);
 double curr_d,d_curr = 0.0;

 vector<TMESH::vertex_type >:: iterator iv;
 while(!frontier.empty()){

	 expansion.clear();

	 pop_heap(frontier.begin(),frontier.end(),lessThan);
	 curr = frontier.back();
	 frontier.pop_back();
	 d_curr = TD[curr].d;
	 TD[curr].visited=true;

	 isLeaf = true;

	 TMESH::EdgePos x;int k;

	 for( x.f = curr->Fp(), x.z = curr->Zp(); x.f!=0; x.NextF() )
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
		 const double & d_pw1  = TD[pw1].d;

		 if((!TD[pw1].visited ) ||(d_pw1 == 0.0) || 
			 (TD[pw].visited && (TDS[pw].source == TDS[curr].source))
				)
			 continue;

			if( (TDS[pw].source!=TDS[curr].source) && (TDS[pw].source!=NULL)&& (TDS[curr].source!=NULL)
				)
				if(AreConn(TDS[pw].source,TDS[curr].source)) continue;
				else
							found=true;

			if( (TDS[pw1].source!=TDS[curr].source) && (TDS[pw1].source!=NULL)&& (TDS[curr].source!=NULL))
				{
				if(AreConn(TDS[pw1].source,TDS[curr].source))
					continue;
				else
					{
					pw = pw1;
					found = true;
					}
				}
			if(found){
					s = TDS[pw].source;
					TMESH::vertex_pointer v0 = pw,v1=curr;
					while(v0!=NULL){ 
						//v0->Q() = 2.0;
						assert(TDS[v0].source == s);
						path.push_front(v0);
						v0 = TDP[v0].parent;
						}	

					s=TDS[v1].source;
					while(v1!=NULL){
					//v1->Q() = 10;		
							assert(TDS[v1].source == s);
						path.push_back(v1);
						v1 = TDP[v1].parent;
						}
				Conn(TDS[pw].source,TDS[curr].source);
				assert(AreConn(TDS[pw].source,TDS[curr].source));
				found = false;
				continue;
				}


		 assert(TDS[pw1].source  == TDS[curr].source);
		 assert((TDP[curr].parent == NULL) || TDS[curr].source == TDS[TDP[curr].parent].source);
		 assert( TD[pw1].d != -1);
		 assert( TD[pw1].d != 0.0);
		 assert( (curr!=pw) && (pw!=pw1) && (pw1 != curr));				
		 assert(d_pw1!=-1.0);

		 //************** calcolo della distanza di pw in base alle distanze note di pw1 e curr
		 //************** sapendo che (curr,pw,pw1) e'una faccia della mesh
		 //************** (vedi figura in file distance.gif)
		 Point3<TMESH::scalar_type> w_c = pw->cP()- curr->cP();
		 Point3<TMESH::scalar_type> w_w1 = pw->cP()- pw1->cP();
		 Point3<TMESH::scalar_type> w1_c = pw1->cP()- curr->cP();

		 double ew_c  = (w_c).Norm();
		 double ew_w1 = (w_w1).Norm();
		 double ec_w1 = (w1_c).Norm();
		 double	alpha,alpha_, beta,beta_,theta_c,theta,h,delta,s,a,b;

		 alpha = acos((w_c*w1_c)/(ew_c*ec_w1));
		 s = (d_curr + d_pw1+ec_w1)/2;
		 a = s/ec_w1;
		 b = a*s;
		 alpha_ = 2*acos ( min(1.0,sqrt(  (b- a* d_pw1)/d_curr)));

		 if ( alpha+alpha_ > M_PI){
			 curr_d = d_curr + ew_c;		
			 }else
			 {
			 beta_ = 2*acos ( min(1.0,sqrt(  (b- a* d_curr)/d_pw1)));
			 beta  = acos((w_w1)*(-w1_c)/(ew_w1*ec_w1));

			 if ( beta+beta_ > M_PI){
				 curr_d = d_pw1  + ew_w1;
				 }
			 else 
				 {
				 theta	= M_PI-alpha-alpha_;
				 delta	= cos(theta)* ew_c;
				 h		= sin(theta)* ew_c;
				 curr_d = sqrt( pow(h,2)+ pow(d_curr + delta,2));
				 }
			 }
		 //**************************************************************************************

		 toQueue = (TD[(pw)].d==-1);

		 if(toQueue){// se non e'gia' in coda ce lo mette
				 if(( TD[curr].d + ew_c) < (d_pw1  + ew_w1))
					 TDP[(pw)].parent = curr;
				 else
					 TDP[(pw)].parent = pw1;
				 assert(curr_d!=0.0);
				 TDS[pw].source = TDS[curr].source;

			 expansion.push_back(pair<TMESH::vertex_pointer,double>(pw,curr_d));
			 isLeaf =false;
			 }else{
				 if(  TD[(pw)].d > curr_d )
					 {
						if(( d_curr + ew_c) < (d_pw1  + ew_w1))
							TDP[(pw)].parent = curr;
						else
							TDP[(pw)].parent = pw1;

						TD[(pw)].d = curr_d;
						TDS[pw].source = TDS[curr].source;
					 isLeaf =false;
					 }
				 }

	 vector <pair<TMESH::vertex_pointer,double> > ::iterator i;
	// unique(expansion.begin(),expansion.end());
	 for(i = expansion.begin(); i!= expansion.end(); ++i){
		 if(TD[(*i).first].d== -1.0){
			 TD[(*i).first].d = (*i).second;
			 frontier.push_back((*i).first);
			 push_heap(frontier.begin(),frontier.end(),lessThan);
			 }
		 }
	 }
}// end while
// scrivi le distanze sul campo qualita' (nn: farlo parametrico)
/* 
TMESH::vertex_iterator vi;
 for(vi = M.vert.begin(); vi != M.vert.end(); ++vi)
	 (*vi).Q() = TD[&(*vi)].d; 
 */


	TD.Stop();	
	TDP.Stop();
	TDS.Stop;

}
public:
template <bool SKIP_SELECTED>
void FartestPoint(  vector<TMESH::vertex_pointer> & fro,//insieme di vertici da cui trovare le distanze
				  TMESH::vertex_pointer & fartest,		//punto piu'lontano
				  double & distance){					//distaza geodesica
				  deque<TMESH::vertex_pointer> _;
				  fartest = BuildSP<false,false,false,SKIP_SELECTED> (_,fro,distancee); 
					  }
template <bool SKIP_SELECTED>
void FartestPoint( vector<TMESH::vertex_pointer> & fro, //insieme di vertici da cui trovare le distanze
				  TMESH::vertex_pointer & fartest,	    //punto piu'lontano
				  double & distance,					//distaza geodesica
				  deque<TMESH::vertex_pointer> & path){// ritorna esplicitamente il cammino minimo
					vector<TMESH::vertex_pointer> _;
					fartest =  BuildSP<true,false,false,SKIP_SELECTED>(path,fro,distance); 
					  }
template <bool SKIP_SELECTED>
void FartestBPoint( vector<TMESH::vertex_pointer> & fro, //insieme di vertici da cui trovare le distanze
				  TMESH::vertex_pointer & fartest,	    //punto piu'lontano
				  double & distance,					//distaza geodesica
				  deque<TMESH::vertex_pointer> & path){// ritorna esplicitamente il cammino minimo
					vector<TMESH::vertex_pointer> _;
					fartest =  BuildSP<true,false,false,SKIP_SELECTED>(path,fro,distance); 
					  }
template <bool SKIP_SELECTED>
void FartestBPoint( vector<TMESH::vertex_pointer> & fro, //insieme di vertici da cui trovare le distanze
				  TMESH::vertex_pointer & fartest,	    //punto piu'lontano
				  double & distance){					//distaza geodesica
					deque<TMESH::vertex_pointer> _;
					fartest =  BuildSP<true,false,true,SKIP_SELECTED>(_,fro,distance); 
					  }

template <bool SKIP_SELECTED>
void FindPath(	 const TMESH::vertex_pointer & v0,	//vertice di partenza
				 const TMESH::vertex_pointer & v1,	//vertice cercato
				deque<TMESH::vertex_pointer> & path,
				double & dist)
					{
					vector<TMESH::vertex_pointer> fro;
					fro.push_back(v0);
					fro.push_back(v1);
					BuildSP<true,true,false,SKIP_SELECTED>(path,fro,dist);

				//	v0 = fro[0];
				//	v1 = fro[1];
					}

template <bool SKIP_SELECTED>
bool MostInternal(	TMESH::vertex_pointer & v0,	//ritorna il vertice piu'lontano da ogni punto sul bordo
					TMESH::vertex_pointer & v1, // ritorna il vertice di bordo piu'vicino a v0
					double & distance,			// distanza geodesica tra v0 e v1
					deque<TMESH::vertex_pointer> & path)// ritorna il path tra i due
					{
					vector<TMESH::vertex_pointer> border_vertices;
					
					TMESH::vertex_iterator ii;
					double maxDistance = 0.0;

					for(ii = M.vert.begin(); ii != M.vert.end();++ii)
						if((*ii).IsB())
							border_vertices.push_back(&(*ii));

					if(border_vertices.empty()) return false;

					vector<TMESH::vertex_pointer> _;
					
					BuildSP<true,false,false,SKIP_SELECTED>(path,border_vertices,distance);

					if(path.empty()) return false;	
					v0 = *path.begin();
					v1 = path.back();
					return true;
					}

template <bool SKIP_SELECTED>
bool MostInternal(	TMESH::vertex_pointer & v0,	//ritorna il vertice piu'lontano da ogni punto sul bordo
					TMESH::vertex_pointer & v1, // ritorna il vertice di bordo piu'vicino a v0
					double & distance// distanza geodesica tra v0 e v1
					)
					{
					vector<TMESH::vertex_pointer> border_vertices;
					deque<TMESH::vertex_pointer> path;
					TMESH::vertex_iterator ii;
					double maxDistance = 0.0;

					for(ii = M.vert.begin(); ii != M.vert.end();++ii)
						if((*ii).IsB())
							border_vertices.push_back(&(*ii));
					
					if(border_vertices.empty()) return false;
					BuildSP<true,false,false,SKIP_SELECTED>(path,border_vertices,distance);
					if(path.empty()) return false;	
					v0 = *path.begin();
					v1 = path.back();
					return true;		
					}

};
}
#endif