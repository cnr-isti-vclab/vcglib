#ifndef __VCGLIB_SPATIAL_ITERATORS
#define __VCGLIB_SPATIAL_ITERATORS

#include <vector>
//#include <vcg/space/intersection3.h>
#include <vcg/space/point3.h>
#include <vcg/space/box3.h>
//#include <vcg/space/point4.h>
#include <vcg/space/ray3.h>
#include <vcg/math/base.h>
#include <algorithm>
#include <float.h>


namespace vcg{
	template <class Spatial_Idexing,class INTFUNCTOR,class TMARKER> 
	class RayIterator
	{	
	protected:
		typedef typename Spatial_Idexing::ObjType ObjType;
		typedef typename Spatial_Idexing::ScalarType ScalarType;
		typedef typename vcg::Point3<ScalarType>  CoordType;
		typedef typename vcg::Ray3<ScalarType> RayType;
		typedef typename Spatial_Idexing::CellIterator CellIterator;

		///control right bonding current cell index (only on initialization)
		void _ControlLimits()
		{
			for (int i=0;i<3;i++)
			{
				vcg::Point3i dim=Si.siz;
				if (CurrentCell.V(i)<0)
					CurrentCell.V(i) = 0;
				else 
					if (CurrentCell.V(i)>=dim.V(i))
						CurrentCell.V(i)=dim.V(i)-1;
			}
		}

		///find initial line parameters
		void _FindLinePar()
		{
			/* Punti goal */

			///da verificare se vanno oltre ai limiti
			vcg::Point3i ip;
			Si.PToIP(start,ip);
			Si.IPToP(ip,goal);
			for (int i=0;i<3;i++)
				if(r.Direction().V(i)>0.0) 
					goal.V(i)+=Si.voxel.V(i);

			ScalarType gx=goal.X();
			ScalarType gy=goal.Y();
			ScalarType gz=goal.Z();

			dist=(r.Origin()-goal).Norm(); 

      const float LocalMaxScalar = std::numeric_limits<float>::max();
			const float	EPSILON = 1e-50f;

			/* Parametri della linea */
			ScalarType tx,ty,tz;

			if(	fabs(r.Direction().X())>EPSILON	) 
				tx = (gx-r.Origin().X())/r.Direction().X();
			else 
				tx	=LocalMaxScalar;

			if(	fabs(r.Direction().Y())>EPSILON	)
				ty = (gy-r.Origin().Y())/r.Direction().Y();
			else 
				ty	=LocalMaxScalar;

			if(	fabs(r.Direction().Z())>EPSILON	) 
				tz = (gz-r.Origin().Z())/r.Direction().Z();
			else 
				tz	=LocalMaxScalar;

			t=CoordType(tx,ty,tz);
		}

		bool _controlEnd()
		{
			return  (((CurrentCell.X()<0)||(CurrentCell.Y()<0)||(CurrentCell.Z()<0))||
				((CurrentCell.X()>=Si.siz.X())||(CurrentCell.Y()>=Si.siz.Y())||(CurrentCell.Z()>=Si.siz.Z())));
		}

		void _NextCell()
		{
			assert(!end);
			if( t.X()<t.Y() && t.X()<t.Z() )
			{
				if(r.Direction().X()<0.0) 
				{goal.X() -= Si.voxel.X(); --CurrentCell.X();}
				else 
				{goal.X() += Si.voxel.X(); ++CurrentCell.X();}
				t.X() = (goal.X()-r.Origin().X())/r.Direction().X();
			} 
			else if( t.Y()<t.Z() ){
				if(r.Direction().Y()<0.0) 
				{goal.Y() -= Si.voxel.Y(); --CurrentCell.Y();}
				else  
				{goal.Y() += Si.voxel.Y(); ++CurrentCell.Y();}
				t.Y() = (goal.Y()-r.Origin().Y())/r.Direction().Y();
			} else {
				if(r.Direction().Z()<0.0) 
				{ goal.Z() -= Si.voxel.Z(); --CurrentCell.Z();}
				else
				{ goal.Z() += Si.voxel.Z(); ++CurrentCell.Z();}
				t.Z() = (goal.Z()-r.Origin().Z())/r.Direction().Z();
			}

			dist=(r.Origin()-goal).Norm();
			end=_controlEnd();
		}

	public:


		///contructor
		RayIterator(Spatial_Idexing &_Si,INTFUNCTOR _int_funct):Si(_Si),int_funct(_int_funct){
		};

		void SetMarker(TMARKER _tm)
		{
			tm=_tm;
		}

		void Init(const RayType _r)
		{
			r=_r;
			end=false;
			tm.UnMarkAll();
			Elems.clear();
			//CoordType ip;
			//control if intersect the bounding box of the mesh
			if(vcg::Intersection_Ray_Box<ScalarType>(Si.bbox,r,start))
			{
				Si.PToIP(start,CurrentCell);
				_ControlLimits();
				_FindLinePar();
				//go to first intersection
				while ((!End())&& Refresh())
					_NextCell();
			}
			else
				end=true;
		}

		bool End()
		{return end;}


		///refresh current cell intersection , return true if there are 
		///at lest 1 intersection
		bool Refresh()
		{
			//Elems.clear();

      typename Spatial_Idexing::CellIterator first,last,l;

			///take first, last iterators to elements in the cell
			Si.Grid(CurrentCell.X(),CurrentCell.Y(),CurrentCell.Z(),first,last);
			for(l=first;l!=last;++l)
			{
				ObjType* elem=&(*(*l));
				ScalarType t;
				CoordType Int;
				if((!tm.IsMarked(elem))&&(int_funct((**l),r,t)))
				{
					Int=r.Origin()+r.Direction()*t;
					Elems.push_back(Entry_Type(elem,t,Int));
					tm.Mark(elem);
				}
			}
			////then control if there are more than 1 element
			if (Elems.size()>1)
				std::sort(Elems.begin(),Elems.end());

			CurrentElem=Elems.end();
			if (Elems.size() > 0) {
				CurrentElem--;
			}

			return((Elems.size()==0)||(Dist()>dist));
		}

		void operator ++()
		{
			//if (CurrentElem!=Elems.end())
			if (Elems.size()>0)
			{
				CurrentElem--;
				//std::pop_heap<ElemIterator>(Elems.begin(),Elems.end());
				Elems.pop_back();
			}
			/*if (CurrentElem==Elems.end())
			{*/
			if (Dist()>dist)
			{
				if (!End())
				{
					_NextCell();
					while ((!End())&&Refresh())
						_NextCell();
				}
			}
		}

		ObjType &operator *(){return *((*CurrentElem).elem);}

		CoordType IntPoint()
		{return ((*CurrentElem).intersection);}

		ScalarType Dist()
		{
			if (Elems.size()>0)
				return ((*CurrentElem).dist);
			else 
				return ((ScalarType)FLT_MAX);
		}
		//{return ((*CurrentElem).dist);}

		///set the current spatial indexing structure used
		void SetIndexStructure(Spatial_Idexing &_Si)
		{Si=_Si;}



	protected:

		///structure that mantain for the current cell pre-calculated data 
		typedef struct Entry_Type
		{
		public:

			Entry_Type(ObjType* _elem,ScalarType _dist,CoordType _intersection)
			{
				elem=_elem;
				dist=_dist;
				intersection=_intersection;
			}
			inline bool operator <  ( const Entry_Type & l ) const{return (dist > l.dist); } 
			ObjType* elem;
			ScalarType dist;
			CoordType intersection;
		};

		RayType r;				  //ray to find intersections
		Spatial_Idexing &Si;	  //reference to spatial index algorithm
		bool end;				  //true if the scan is terminated
		INTFUNCTOR &int_funct;
		TMARKER tm;

		std::vector<Entry_Type> Elems;					//element loaded from curren cell
		typedef typename std::vector<Entry_Type>::iterator ElemIterator;
		ElemIterator CurrentElem;	//iterator to current element

		vcg::Point3i CurrentCell;						//current cell

		//used for raterization
		CoordType start;
		CoordType goal;	
		ScalarType dist;
		CoordType t;

	};


	template <class Spatial_Idexing,class DISTFUNCTOR,class TMARKER> 
	class ClosestIterator
	{	
		typedef typename Spatial_Idexing::ObjType ObjType;
		typedef typename Spatial_Idexing::ScalarType ScalarType;
		typedef typename vcg::Point3<ScalarType>  CoordType;
		typedef typename Spatial_Idexing::CellIterator CellIterator;



		///control the end of scanning
		bool  _EndGrid()
		{
			if ((explored.min==vcg::Point3i(0,0,0))&&(explored.max==Si.siz-vcg::Point3i(1,1,1)))
				end =true;
			return end;
		}

		void _UpdateRadius()
		{
			if (radius>=max_dist)
				end=true;

			radius+=step_size;
			//control bounds
			if (radius>max_dist)
				radius=max_dist;
		}

		///add cell to the curren set of explored cells
		bool _NextShell()
		{

			//then expand the box
			explored=to_explore;
			_UpdateRadius();
			Box3<ScalarType> b3d(p,radius);
			/*b3d.Intersect(Si.bbox);
			Si.BoxToIBox(b3d,to_explore);*/
			Si.BoxToIBox(b3d,to_explore);
			Box3i ibox(Point3i(0,0,0),Si.siz-Point3i(1,1,1));
			to_explore.Intersect(ibox);
			if (!to_explore.IsNull())
			{
				assert(!( to_explore.min.X()<0 || to_explore.max.X()>=Si.siz[0] ||
					to_explore.min.Y()<0 || to_explore.max.Y()>=Si.siz[1] ||  to_explore.min.Z()<0 
					|| to_explore.max.Z()>=Si.siz[2] ));
				return true;
			}
			return false;
		}



	public:

		///contructor
		ClosestIterator(Spatial_Idexing &_Si,DISTFUNCTOR _dist_funct):Si(_Si),dist_funct(_dist_funct){}

		///set the current spatial indexing structure used
		void SetIndexStructure(Spatial_Idexing &_Si)
		{Si=_Si;}

		void SetMarker(TMARKER _tm)
		{
			tm=_tm;
		}

		///initialize the Itarator
		void Init(CoordType _p,const ScalarType &_max_dist)
		{
			explored.SetNull();
			to_explore.SetNull();
			p=_p;
			max_dist=_max_dist;
			Elems.clear();
			end=false;
			tm.UnMarkAll();
			//step_size=Si.voxel.X();
			step_size=Si.voxel.Norm();
			radius=0;

			///inflate the bbox until find a valid bbox
			while ((!_NextShell())&&(!End()));

			if (!_EndGrid())
				Refresh();///load elements form currect cell

			///until don't find an element
			///that is inside the radius
			while ((!End())&&(Dist()>radius))
			{
				if ((_NextShell())&&(!_EndGrid()))
					Refresh();
			}

			//set to the last element ..the nearest
			CurrentElem=Elems.end();
			CurrentElem--;

		}

		//return true if the scan is complete
		bool End()
		{return end;}

		///refresh Object found	also considering current shere radius, 
		//and object comes from	previos	that are already in	the	stack
		void Refresh()
		{
			int	ix,iy,iz;
			for( iz = to_explore.min.Z();iz <=	to_explore.max.Z(); ++iz)
				for(iy =to_explore.min.Y(); iy	<=to_explore.max.Y(); ++iy)
					for(ix =to_explore.min.X(); ix	<= to_explore.max.X();++ix)
					{
						// this test is to avoid to re-process already analyzed cells.
						if((explored.IsNull())||
							(ix<explored.min[0] || ix>explored.max[0] ||  
							iy<explored.min[1] || iy>explored.max[1] ||
							iz<explored.min[2] || iz>explored.max[2] ))
						{
							typename Spatial_Idexing::CellIterator first,last,l;

							Si.Grid(ix,iy,iz,first,last);
							for(l=first;l!=last;++l)
							{
								ObjType	*elem=&(**l);
								if (!tm.IsMarked(elem))
								{

									CoordType nearest;
									ScalarType dist=max_dist;
									if (dist_funct((**l),p,dist,nearest))
										Elems.push_back(Entry_Type(elem,fabs(dist),nearest));
									tm.Mark(elem);
								}
							}
						}

					}

					std::sort(Elems.begin(),Elems.end());

					CurrentElem=Elems.end();
					CurrentElem--;
		}

		void operator ++()
		{
			/*if (Dist()<=radius)
			{
				CurrentElem--;
				Elems.pop_back();
			}

			while ((!End())&&(Dist()>radius))
				if (_NextShell()&&!_EndGrid())
					Refresh();*/

			if (Elems.size()>0)
			{
				CurrentElem--;
				Elems.pop_back();
			}
			while ((!End())&&(Dist()>radius))
				if (_NextShell()&&!_EndGrid())
					Refresh();
		}

		ObjType &operator *(){return *((*CurrentElem).elem);}

		//return distance of the element form the point if no element 
		//are in the vector then return max dinstance
		ScalarType Dist()
		{
			if (Elems.size()>0)
				return ((*CurrentElem).dist);
			else 
				return ((ScalarType)FLT_MAX);
		}

		CoordType NearestPoint()
		{return ((*CurrentElem).intersection);}

	protected:

		///structure that mantain for the current cell pre-calculated data 
		typedef struct Entry_Type
		{
		public:

			Entry_Type(ObjType* _elem,ScalarType _dist,CoordType _intersection)
			{
				elem=_elem;
				dist=_dist;
				intersection=_intersection;
			}

			inline bool operator <  ( const Entry_Type & l ) const{return (dist > l.dist); } 

			inline bool operator ==  ( const Entry_Type & l ) const{return (elem == l.elem); } 

			ObjType* elem;
			ScalarType dist;
			CoordType intersection;
		};

		CoordType p;				  //initial point
		Spatial_Idexing &Si;		  //reference to spatial index algorithm
		bool end;					  //true if the scan is terminated
		ScalarType max_dist;		  //max distance when the scan terminate
		vcg::Box3i explored;		  //current bounding box explored
		vcg::Box3i to_explore;		  //current bounding box explored
		ScalarType radius;			  //curret radius for sphere expansion
		ScalarType step_size;		  //radius step
		std::vector<Entry_Type> Elems; //element loaded from the current sphere

		DISTFUNCTOR &dist_funct;
		TMARKER tm;

		typedef typename std::vector<Entry_Type>::iterator ElemIterator;
		ElemIterator CurrentElem;	//iterator to current element

};
}

#endif
