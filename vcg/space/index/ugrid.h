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

#ifndef __VCGLIB_UGRID
#define __VCGLIB_UGRID

namespace vcg {

	/** Static Uniform Grid
	    A spatial search structure for a accessing a container of objects. It is based on a uniform grid overlayed over a protion of space. 
      The grid partion the space into cells. Cells contains just pointers to the object that are stored elsewhere. 
      The set of object is meant to be static and pointer stable. 
      Useful for situation were many space related query are issued over the same dataset (ray tracing, measuring distances between meshes, re-detailing ecc.). Works well for distribution that ar reasonably uniform.
      How to use it:

     */
template < class ContainerType >
class UGrid
{
public:

	/** Internal class for keeping the first pointer of object.
	    Definizione Link dentro la griglia. Classe di supporto per UGrid.
	 */
		class Link
		{
		public:
				/// Costruttore di default
				inline Link(){};
				/// Costruttore con inizializzatori
				inline Link( ContainerType::iterator const nt, const int ni ){
							assert(ni>=0);
							t = nt;
							i = ni;
				};

				
				inline bool operator <  ( const Link & l ) const{ 	return i <   l.i; } 
				inline bool operator <= ( const Link & l ) const{ 	return i <=  l.i; }
				inline bool operator >  ( const Link & l ) const{ 	return i >   l.i; }
				inline bool operator >= ( const Link & l ) const{ 	return i >=  l.i; }
				inline bool operator == ( const Link & l ) const{ 	return i ==  l.i; }
				inline bool operator != ( const Link & l ) const{ 	return i !=  l.i; }
			
				inline ContainerType::iterator & Elem() {
					return t;
				}
				inline int & Index() {
					return i;
				}

		private:
				/// Puntatore all'elemento T
				ContainerType::iterator t;
				/// Indirizzo del voxel dentro la griglia
				int i;
				

		};//end class Link

		typedef ContainerType::value_type ObjType;
		typedef ObjType::ScalarType ScalarType;
		typedef Point3<ScalarType> Point3x;
		typedef Box3<ScalarType> Box3x;
		typedef Line3<ScalarType> Line3x;
		
    Box3x   bbox;
		Point3x dim;   /// Dimensione spaziale (lunghezza lati) del bbox
    Point3i siz;   /// Dimensioni griglia in celle
		Point3x voxel; /// Dimensioni di una cella
    
		/// Insieme di tutti i links
    vector<Link>   links;
		/// Griglia vera e propria
    vector<Link *> grid;

		/// Dato un punto, ritorna la cella che lo contiene
    inline Link ** Grid( const Point3d & p )
		{
			int x = int( (p[0]-bbox.min[0])/voxel[0] );
			int y = int( (p[1]-bbox.min[1])/voxel[1] );
			int z = int( (p[2]-bbox.min[2])/voxel[2] );

#ifndef NDEBUG
			if ( x<0 || x>=siz[0] || y<0 || y>=siz[1] || z<0 || z>=siz[2] )
				return NULL;
			else
#endif

			return grid.begin() + ( x+siz[0]*(y+siz[1]*z) );
		}
		/// Date le coordinate ritorna la cella
    inline Link ** Grid( const int x, const int y, const int z )
		{
#ifndef NDEBUG
			if ( x<0 || x>=siz[0] || y<0 || y>=siz[1] || z<0 || z>=siz[2] )
				assert(0);
				//return NULL;
			else
#endif
			assert(x+siz[0]*(y+siz[1]*z<grid.size()));
			return &*grid.begin() + ( x+siz[0]*(y+siz[1]*z) );
		}

	void Grid( const Point3d & p, Link * & first, Link * & last )
	{
		Link ** g = Grid(s);
		
		first = *g;
		last  = *(g+1);
	}

	void Grid( const int x, const int y, const int z, Link * & first, Link * & last )
	{
		Link ** g = Grid(x,y,z);
		
		first = *g;
		last  = *(g+1);
	}

		/// Setta il bounding box della griglia
    void SetBBox( const Box3x & b )
		{
			bbox = b;
			dim  = b.max - b.min;
		}

    void SetSafeBBox( const Box3x & b )
		{
			Box3x btmp=b;
			btmp.InflateFix(0.01);
			bbox = btmp;
			dim  = bbox.max - bbox.min;
		}

		/// Dato un punto 3d ritorna l'indice del box corrispondente
    inline void PToIP(const Point3x & p, Point3i &pi ) const
		{
			Point3x t = p - bbox.min;
			pi[0] = int( t[0]/voxel[0] );
			pi[1] = int( t[1]/voxel[1] );
			pi[2] = int( t[2]/voxel[2] );
		}
		/// Dato un box reale ritorna gli indici dei voxel compresi dentro un ibox
    void BoxToIBox( const Box3x & b, Box3i & ib ) const
		{
			PToIP(b.min,ib.min);
			PToIP(b.max,ib.max);
		}

	
#ifdef MONITOR_GRID
		void ShowRayCastingStats() {
			if (n_rays>0) {
			printf("\n%d rays (%f faces in %f voxel per ray)\n",
				n_rays,
				double(n_traversed_elem)/n_rays,
				double(n_traversed_voxel)/n_rays);
			printf("HIT: %d rays (%d%%) (%f faces in %f voxel per ray)\n",
				n_rays_hit,
				(n_rays_hit*100)/n_rays,
				double(n_traversed_elem_hit)/n_rays_hit,
				double(n_traversed_voxel_hit)/n_rays_hit);
			}
		}
#endif

	void ShowStats(FILE *fp)
	{
				// Conto le entry
		int nentry = 0;
		Hist H;
		H.SetRange(0,1000,1000);
		int pg;
		for(pg=0;pg<grid.size()-1;++pg)
			if( grid[pg]!=grid[pg+1] )
			{
				++nentry;
				H.Add(grid[pg+1]-grid[pg]);
			}

			fprintf(fp,"Uniform Grid: %d x %d x %d (%d voxels), %.1f%% full, %d links \nNon empty Cell Occupancy Distribution Avg: %f (%4.0f %4.0f %4.0f) \n",
			siz[0],siz[1],siz[2],grid.size()-1,
			double(nentry)*100.0/(grid.size()-1),links.size(),H.Avg(),H.Percentile(.25),H.Percentile(.5),H.Percentile(.75)
      
		);
	}
	void SlicedPPM( const char * filename )
	{
		xstring basename=filename;
		xstring name;
		int ix,iy,iz;
    
		for(iz=0;iz<siz[2];++iz)
		{
			name.format("%s%03i.ppm",filename,iz);
			FILE * fp = fopen(name,"wb");
			fprintf(fp,
				"P6\n"
				"%d %d\n"
				"255\n"
				,siz[0]
				,siz[1]
			);
			for(ix=0;ix<siz[0];++ix)
			{
				for(iy=0;iy<siz[1];++iy)
				{
					Link ** ii= Grid(ix,iy,iz);
					if( *ii != *(ii+1) )
					{
						unsigned char c = 255;
						fwrite(&c,1,1,fp);
						fwrite(&c,1,1,fp);
						fwrite(&c,1,1,fp);
					}
					else
					{
						unsigned char c = 0;
						fwrite(&c,1,1,fp);
						fwrite(&c,1,1,fp);
						fwrite(&c,1,1,fp);
					}
				}
			}
			fclose(fp);
		}
	}
	void PPM( const char * filename )
	{
		int ix,iy,iz;

		FILE * fp = fopen(filename,"wb");
		fprintf(fp,
		  "P6\n"
		  "%d %d\n"
		  "255\n"
		  ,(siz[1]+1)*siz[0]
		  ,siz[2]
		);
		for(iz=0;iz<siz[2];++iz)
		  for(ix=0;ix<siz[0];++ix)
		  {
			for(iy=0;iy<siz[1];++iy)
		{
		  int i = ix+siz[0]*(iy+siz[1]*iz);
		  if( grid[i]!=grid[i+1] )
		  {
			  unsigned char c = 255;
			  fwrite(&c,1,1,fp);
			  fwrite(&c,1,1,fp);
			  fwrite(&c,1,1,fp);
		  }
		  else
		  {
			  unsigned char c = 0;
			  fwrite(&c,1,1,fp);
			  fwrite(&c,1,1,fp);
			  fwrite(&c,1,1,fp);
		  }
			}
		{
			unsigned char c = 0;
			fwrite(&c,1,1,fp);
			fwrite(&c,1,1,fp);
			c = 255;
			fwrite(&c,1,1,fp);
		}
		  }
		fclose(fp);
	}


	/** Returns the closest posistion of a point p and its distance
		@param p a 3d point
		@return The closest element
	*/
	T * GetClosest( const Point3x & p, ScalarType & min_dist, Point3x  & res)
	{
		ScalarType dx = ( (p[0]-bbox.min[0])/voxel[0] );
		ScalarType dy = ( (p[1]-bbox.min[1])/voxel[1] );
		ScalarType dz = ( (p[2]-bbox.min[2])/voxel[2] );

		int ix = int( dx );
		int iy = int( dy );
		int iz = int( dz );

		double voxel_min=voxel[0];
		if (voxel_min<voxel[1]) voxel_min=voxel[1];
		if (voxel_min<voxel[2]) voxel_min=voxel[2];

		ScalarType radius=(dx-ScalarType(ix));
		if (radius>0.5)  radius=(1.0-radius);	radius*=voxel[0];

		/*ScalarType radio_if[6];
		radio_if[0]= (dx-double(ix)    )*voxel[0] ;
		radio_if[1]= (dy-double(iy)    )*voxel[1] ;
		radio_if[2]= (dz-double(iz)    )*voxel[2] ;
		radio_if[3]= (1.0+double(ix)-dx)*voxel[0] ;
		radio_if[4]= (1.0+double(iy)-dy)*voxel[1] ;
		radio_if[5]= (1.0+double(iz)-dz)*voxel[2] ;*/

		ScalarType 
		tmp=dy-ScalarType(iy); if (tmp>0.5) tmp=1.0-tmp; tmp*=voxel[1]; if (radius>tmp) radius=tmp;
		tmp=dz-ScalarType(iz); if (tmp>0.5) tmp=1.0-tmp; tmp*=voxel[2]; if (radius>tmp) radius=tmp;

		Point3x t_res;
		//ScalarType min_dist=1e10;
		T* winner=NULL;

		Link *first, *last, *l;
		if ((ix>=0) && (iy>=0) && (iz>=0) && (ix<siz[0]) && (iy<siz[1]) && (iz<siz[2])) {

			Grid( ix, iy, iz, first, last );
			for(l=first;l!=last;++l)
			{
				if (l->Elem()->Dist(p,min_dist,t_res)) {
					winner=&*(l->Elem());
					res=t_res;

				}
			};
		};

		//return winner;

		Point3i done_min=Point3i(ix,iy,iz), done_max=Point3i(ix,iy,iz);

		//printf(".");

		while (min_dist>radius) {
			//if (dy-ScalarType(iy))
			done_min[0]--; if (done_min[0]<0) done_min[0]=0;
			done_min[1]--; if (done_min[1]<0) done_min[1]=0;
			done_min[2]--; if (done_min[2]<0) done_min[2]=0;
			done_max[0]++; if (done_max[0]>=siz[0]-1) done_max[0]=siz[0]-1;
			done_max[1]++; if (done_max[1]>=siz[1]-1) done_max[1]=siz[1]-1;
			done_max[2]++; if (done_max[2]>=siz[2]-1) done_max[2]=siz[2]-1;
			radius+=voxel_min;
			//printf("+");
			for (ix=done_min[0]; ix<=done_max[0]; ix++) 
			for (iy=done_min[1]; iy<=done_max[1]; iy++) 
			for (iz=done_min[2]; iz<=done_max[2]; iz++) 
			{
				Grid( ix, iy, iz, first, last );
				for(l=first;l!=last;++l)
				{
					if (l->Elem()->Dist(p,min_dist,t_res)) {
						winner=&*(l->Elem());
						res=t_res;
					};
				};
			}
		};
		return winner;
	};

	T * DoRay( Line3x & ray, double & dist , double &p_a, double &p_b)
	{
		int ix,iy,iz;
		double gx,gy,gz;

		if(!bbox.IsIn(ray.orig))
		{
			//printf(".");

			Point3x ip;
			if(Intersection(bbox,ray,ip))
			{
				ix = int( (ip.x() - bbox.min.x())/voxel.x() );
				iy = int( (ip.y() - bbox.min.y())/voxel.y() );
				iz = int( (ip.z() - bbox.min.z())/voxel.z() );
				if(ix<0) ix = 0;
				if(iy<0) iy = 0;
				if(iz<0) iz = 0;
				if(ix>=siz[0]) ix = siz[0]-1;
				if(iy>=siz[1]) iy = siz[1]-1;
				if(iz>=siz[2]) iz = siz[2]-1;
			}
			else
				return 0;
		}
		else
		{
				/* Indici di voxel */
			ix = int( (ray.orig.x() - bbox.min.x())/voxel.x() );
			iy = int( (ray.orig.y() - bbox.min.y())/voxel.y() );
			iz = int( (ray.orig.z() - bbox.min.z())/voxel.z() );
		}

				/* Punti goal */
		if(ray.dire.x()>0.0) gx=double(ix+1)*voxel.x()+bbox.min.x();
		else                 gx=double(ix  )*voxel.x()+bbox.min.x();
		if(ray.dire.y()>0.0) gy=double(iy+1)*voxel.y()+bbox.min.y();
		else                 gy=double(iy  )*voxel.y()+bbox.min.y();
		if(ray.dire.z()>0.0) gz=double(iz+1)*voxel.z()+bbox.min.z();
		else                 gz=double( iz )*voxel.z()+bbox.min.z();

		const ScalarType  MAXFLOAT = 1e50;
		const ScalarType  EPSILON = 1e-50;

		/* Parametri della linea */
		double tx,ty,tz;
		if( fabs(ray.dire.x())>EPSILON ) tx = (gx-ray.orig.x())/ray.dire.x();
		else tx = MAXFLOAT;
		if( fabs(ray.dire.y())>EPSILON ) ty = (gy-ray.orig.y())/ray.dire.y();
		else ty = MAXFLOAT;
		if( fabs(ray.dire.z())>EPSILON ) tz = (gz-ray.orig.z())/ray.dire.z();
		else tz = MAXFLOAT;

		T  * bestf = NULL;

#ifdef MONITOR_GRID
		int n_traversed_voxel0=0,n_traversed_elem0=0; 
		n_rays++;
#endif
		for(;;)
		{
			Link *first, *last, *l;
			Grid( ix, iy, iz, first, last );
#ifdef MONITOR_GRID
			n_traversed_voxel++;n_traversed_voxel0++;
#endif
			for(l=first;l!=last;++l)
			{
#ifdef MONITOR_GRID
				n_traversed_elem++;n_traversed_elem0++;
#endif
				//return &(*(l->Elem()));
				if(l->Elem()->Intersect(ray,dist, p_a, p_b)) {
#ifdef MONITOR_GRID
					n_traversed_voxel_hit+=n_traversed_voxel0;
					n_traversed_elem_hit+=n_traversed_elem0;
					n_rays_hit++;
#endif

					return &(*(l->Elem()));
					//{
					//	bestf = &(*(l->Elem())); break;
					//}
				}
				//if(bestf) break;
			}

			if( tx<ty && tx<tz ){
				if(ray.dire.x()<0.0) { gx -= voxel.x(); --ix; if(ix<0      ) break; }
				else                 { gx += voxel.x(); ++ix; if(ix>=siz[0]) break; }
				tx = (gx-ray.orig.x())/ray.dire.x();
			} else if( ty<tz ){
				if(ray.dire.y()<0.0) { gy -= voxel.y(); --iy; if(iy<0      ) break; }
				else                 { gy += voxel.y(); ++iy; if(iy>=siz[1]) break; }
				ty = (gy-ray.orig.y())/ray.dire.y();
			} else {
				if(ray.dire.z()<0.0) { gz -= voxel.z(); --iz; if(iz<0      ) break; }
				else                 { gz += voxel.z(); ++iz; if(iz>=siz[2]) break; }
				tz = (gz-ray.orig.z())/ray.dire.z();
			}
		}
		return bestf;
	}
  
	
	/// Inserisce una mesh nella griglia. Nota: prima bisogna 
	/// chiamare SetBBox che setta dim in maniera corretta
	void Set( S & s )
	{
		Set(s,s.size());
	}


	/// Inserisce una mesh nella griglia. Nota: prima bisogna 
	/// chiamare SetBBox che setta dim in maniera corretta
	void Set( S & s,int _size )
	{
		Point3i _siz;

		BestDim( _size, dim, _siz );
		Set(s,_siz);
	}
	void Set(S & s, Point3i _siz)
	{	
		siz=_siz;
			// Calcola la dimensione della griglia
		voxel[0] = dim[0]/siz[0];
		voxel[1] = dim[1]/siz[1];
		voxel[2] = dim[2]/siz[2];

		// "Alloca" la griglia: +1 per la sentinella
		grid.resize( siz[0]*siz[1]*siz[2]+1 );

    		// Ciclo inserimento dei tetraedri: creazione link
		links.clear();
		S::iterator pt;
		for(pt=s.begin(); pt!=s.end(); ++pt)
		{
			Box3x bb;			// Boundig box del tetraedro corrente
				(*pt).GetBBox(bb);
				bb.Intersect(bbox);
				if(! bb.IsNull() )
				{
					Box3i ib;		// Boundig box in voxels
					BoxToIBox( bb,ib );
					/*ib.min-=Point3i(1,1,1);
					ib.max+=Point3i(1,1,1);
					if (ib.min[0]<0) ib.min[0]=0;
					if (ib.min[1]<0) ib.min[1]=0;
					if (ib.min[2]<0) ib.min[2]=0;
					if (ib.max[0]>siz[0]) ib.max[0]=siz[0];
					if (ib.max[1]>siz[1]) ib.max[1]=siz[1];
					if (ib.max[2]>siz[2]) ib.max[2]=siz[2];*/

					int x,y,z;
					for(z=ib.min[2];z<=ib.max[2];++z)
					{
						 int bz = z*siz[1];
						 for(y=ib.min[1];y<=ib.max[1];++y)
						 {
							 int by = (y+bz)*siz[0];
							 for(x=ib.min[0];x<=ib.max[0];++x)
								// Inserire calcolo cella corrente
								// if( pt->Intersect( ... )
								links.push_back( Link(pt,by+x) );
						 }
					}
				}
		}
		// Push della sentinella
		links.push_back( Link(NULL,grid.size()-1) );

		// Ordinamento dei links
		sort( links.begin(), links.end() );

		// Creazione puntatori ai links
		vector<Link>::iterator pl;
		int pg;
		pl = links.begin();
		for(pg=0;pg<grid.size();++pg)
		{
			assert(pl!=links.end());

			grid[pg] = &*pl;
			while( pg == pl->Index() )	// Trovato inizio
			{
				++pl;		// Ricerca prossimo blocco
				if(pl==links.end())
				break;
			}
		}

	}

		/** Calcolo dimensioni griglia.
		 Calcola la dimensione della griglia in funzione
			 della ratio del bounding box e del numero di elementi
	 */
		static void BestDim( const int elems, const Point3x & size, Point3i & dim )
		{
				const int mincells   = 27;		// Numero minimo di celle
				const double GFactor = 1.0;		// GridEntry = NumElem*GFactor
				double diag = size.Norm();		// Diagonale del box
				double eps  = diag*1e-4;		// Fattore di tolleranza

				assert(elems>0);
				assert(size[0]>=0.0);
				assert(size[1]>=0.0);
				assert(size[2]>=0.0);

				int ncell = int(elems*GFactor);	// Calcolo numero di voxel
				if(ncell<mincells)
				ncell = mincells;

				dim[0] = 1;
				dim[1] = 1;
				dim[2] = 1;

				if(size[0]>eps)
				{
					if(size[1]>eps)
					{
						if(size[2]>eps)
						{
							double k = pow(ncell/(size[0]*size[1]*size[2]),1.0/3.0);
							dim[0] = int(size[0] * k);
							dim[1] = int(size[1] * k);
							dim[2] = int(size[2] * k);
						} 
						else 
						{
							dim[0] = int(::sqrt(ncell*size[0]/size[1]));
							dim[1] = int(::sqrt(ncell*size[1]/size[0]));
						}
					}
					else
					{
						if(size[2]>eps)
						{
							dim[0] = int(::sqrt(ncell*size[0]/size[2]));
							dim[2] = int(::sqrt(ncell*size[2]/size[0]));
						}
						else
							dim[0] = int(ncell);
					}
				}
				else
				{
					if(size[1]>eps)
					{
						if(size[2]>eps)
						{
							dim[1] = int(::sqrt(ncell*size[1]/size[2]));
							dim[2] = int(::sqrt(ncell*size[2]/size[1]));
						}
						else
							dim[1] = int(ncell);
					}
					else if(size[2]>eps)
						dim[2] = int(ncell);
				}
		}


	int MemUsed()
	{
		return sizeof(UGrid)+ sizeof(Link)*links.size() + sizeof(Link *) * grid.size();
	}
}; //end class UGrid

} // end namespace

#endif
