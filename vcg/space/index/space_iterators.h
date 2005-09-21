#include <vector>
#include <vcg/space/box3.h>
#include <algorithm>
#include <float.h>
#include <vcg/space/intersection3.h>
#include <vcg/simplex/face/distance.h>//to remove
#ifndef _VCG_SPACE_ITERATORS
#define _VCG_SPACE_ITERATORS

namespace vcg{
template <class Spatial_Idexing,class INTFUNCTOR,class TMARKER> 
class RayIterator
{	
protected:
	typedef typename Spatial_Idexing::ObjType ObjType;
	typedef typename Spatial_Idexing::ScalarType ScalarType;
	typedef typename vcg::Point3<ScalarType>  CoordType;
	typedef typename vcg::Line3<ScalarType> RayType;
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

	const float	 MAXFLOAT =	 FLT_MAX;
	const float	 EPSILON = 1e-50f;

	/* Parametri della linea */
	ScalarType tx,ty,tz;

	if(	fabs(r.Direction().X())>EPSILON	) 
		tx = (gx-r.Origin().X())/r.Direction().X();
	else 
		tx	=MAXFLOAT;

	if(	fabs(r.Direction().Y())>EPSILON	)
		ty = (gy-r.Origin().Y())/r.Direction().Y();
	else 
		ty	=MAXFLOAT;

	if(	fabs(r.Direction().Z())>EPSILON	) 
		tz = (gz-r.Origin().Z())/r.Direction().Z();
	else 
		tz	=MAXFLOAT;

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

	void Init(RayType _r)
	{
		r=_r;
		end=false;
		tm.UnMarkAll();
		Elems.clear();
		//CoordType ip;
		//control if intersect the bounding box of the mesh
		if(vcg::Intersection<ScalarType>(Si.bbox,r,start))
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

		Spatial_Idexing::CellIterator first,last,l;

		///take first, last iterators to elements in the cell
		Si.Grid(CurrentCell.X(),CurrentCell.Y(),CurrentCell.Z(),first,last);
		for(l=first;l!=last;++l)
		{
			ObjType* elem=&(*(*l));
		/*	CoordType p0 =  elem->V(0)->P();
			CoordType p1 =  elem->V(1)->P();
			CoordType p2 =  elem->V(2)->P();*/

			//ScalarType dist=0;
			CoordType Int;
			if((!tm.IsMarked(elem))&&(int_funct((**l),r,Int)))
			{
				Elems.push_back(Entry_Type(elem,(r.Origin()-Int).Norm(),Int));
				//std::push_heap<ElemIterator>(Elems.begin(),Elems.end());
				tm.Mark(elem);
			}
			
			/////different types of intersections as required
			//if ((!_Unique)&&(_DoubleIntersect(r,p0,p1,p2,Int)))
			//	Elems.push_back(Entry_Type(elem,(r.Origin()-Int).Norm(),Int));
			//else
			//if ((_Unique)&&(_UniqueIntersection(r,p0,p1,p2,Int)))
			//	Elems.push_back(Entry_Type(elem,(r.Origin()-Int).Norm(),Int));
		}
		////then control if there are more than 1 element
		if (Elems.size()>1)
				std::sort(Elems.begin(),Elems.end());
		
		CurrentElem=Elems.end();
		CurrentElem--;

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


	///find the radius of curren sphere
	///considering more nearest cell to current radius
	void _FindSphereRadius()
	{
		//find diffence between the initial point 
		//and the cell narest to the sphere
		CoordType min;
		CoordType max;
		Si.IPToP(explored.min,min);
		Si.IPToP(explored.max,max);
		CoordType diff_min=p-min;
		CoordType diff_max=p-max;

		ScalarType diffx=std::min<ScalarType>(fabs(diff_min.X()),fabs(diff_max.X()));
		ScalarType diffy=std::min<ScalarType>(fabs(diff_min.Y()),fabs(diff_max.Y()));
		ScalarType diffz=std::min<ScalarType>(fabs(diff_min.Z()),fabs(diff_max.Z()));

		ScalarType diff=std::min<ScalarType>(diffx,std::min<ScalarType>(diffy,diffz));

		//radius_min=0;
		radius=diff;
		//radius+=diff;
		if (radius>max_dist)
			radius=max_dist;
	}
	
	///control right cell value of current bounding box 
	/// of explored cells
	void _ControlLimits()
	{
		vcg::Point3i dim=Si.siz;
		for (int i=0;i<3;i++)
		{
			if (explored.min.V(i)<-1)
				explored.min.V(i) = -1;
			if (explored.max.V(i)>Si.siz.V(i))
				explored.max.V(i) =Si.siz.V(i);
		}
	}

	bool _OutOfLimits(vcg::Point3i p)
	{
		for (int i=0;i<3;i++)
		 if ((p.V(i)==-1)||(p.V(i)>=Si.siz.V(i)))
			return true	;

		return false;
	}

	///control the end of scanning
	void  _ControlEnd()
	{
		if ((explored.min==vcg::Point3i(-1,-1,-1))&&(explored.max==Si.siz)&&(Elems.size()==0))
				end =true;
	}

	///add cell to the curren set of explored cells
	void _NextShell()
	{
			if (radius>=max_dist)
				end=true;

			radius+=voxel_min;
			//control bounds
			if (radius>max_dist)
				radius=max_dist;

			//expand the box
			explored.min-=vcg::Point3i(1,1,1);
			explored.max+=vcg::Point3i(1,1,1);

			//control right limits of the bound
			_ControlLimits();			
	}



public:

	///contructor
	ClosestIterator(Spatial_Idexing &_Si,DISTFUNCTOR _dist_funct):Si(_Si),dist_funct(_dist_funct){}

	///set the current spatial indexing structure used
	void SetIndexStructure(Spatial_Idexing &_Si)
	{Si=_Si;}

	///initialize the Itarator
	void Init(CoordType _p,const ScalarType &_max_dist)
	{
		//CoordType vox=Si.Voxel();
		CoordType vox=Si.voxel;

		voxel_min=std::min<ScalarType>(vox.V(0),std::min<ScalarType>(vox.V(1),vox.V(2)));
		p=_p;
		max_dist=_max_dist;
		
		//initialize the explored region
		//finding cell coordinate of initial point
		vcg::Point3i c;
		Si.PToIP(p,c);

		assert( c.X()>=0 && c.X()<Si.siz.X() && c.Y()>=0 && c.Y()<Si.siz.Y() && c.Z()>=0 && c.Z()<Si.siz.Z() );
			
		//explored=vcg::Box3i(c,c+vcg::Point3i(1,1,1));
		explored=vcg::Box3i(c,c);
		
		radius=0;
		//radius_min=0;

		Elems.clear();
		end=false;

		tm.UnMarkAll();

		_FindSphereRadius();
		Refresh();
		///until don't find an element
		///that is inside the radius
		while ((!End())&&(Dist()>radius))
		{
			if (radius>=max_dist)
				end=true;
			_NextShell();
			Refresh();
			_ControlEnd();
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
		int	x,y,z;
		for( z = explored.min.Z(); z <=	explored.max.Z(); ++z)
			for(y =	explored.min.Y(); y	<=explored.max.Y();	 ++y)
				for(x =	explored.min.X(); x	<= explored.max.X();)
				{
					/*vcg::Point3i CurrentCell=vcg::Point3i(x,y,z);*/
					Spatial_Idexing::CellIterator first,last,l;
					///take	first, last	iterators to elements in the cell
					if (!_OutOfLimits(vcg::Point3i(x,y,z)))
					{
						Si.Grid(x,y,z,first,last);
						for(l=first;l!=last;++l)
						{
							ObjType	*elem=&(**l);
							if (!tm.IsMarked(elem))
							{
								
								CoordType nearest;
								ScalarType dist;
								dist_funct((**l),p,dist,nearest);
								Elems.push_back(Entry_Type(elem,fabs(dist),nearest));
								tm.Mark(elem);
							}
						}
					}
					if(	( (	y == explored.min.Y()) || (	y == explored.max.Y()))	||
						( (	z == explored.min.Z()) || (	z == explored.max.Z()))	||
						( x	== explored.max.X()))
						++x;
					else 
						x=explored.max.X();
				}

				std::sort(Elems.begin(),Elems.end());
				//std::unique(Elems.begin(),Elems.end());
				
				CurrentElem=Elems.end();
				CurrentElem--;
	}

	void operator ++()
	{
		if (Elems.size()>0)
		{
			CurrentElem--;
			Elems.pop_back();
		}
		if (Dist()>radius)
		{
			_NextShell();
			Refresh();
			//continue to scan until finish the scanning or the
			//first element (the nearest for ordering of the structure) 
			//is at distance<radius
			while ((!End())&&(Dist()>radius))
			{
				if (radius>=max_dist)
					end=true;
				_NextShell();
				Refresh();
				_ControlEnd();
			}
		}
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
	ScalarType radius;			  //curret radius for sphere expansion
	//ScalarType radius_min;		  //curret radius of explored simplexes
	ScalarType voxel_min;		  //minimum value of the voxel
	std::vector<Entry_Type> Elems; //element loaded from the current sphere
	
	DISTFUNCTOR &dist_funct;
	TMARKER tm;

	typedef typename std::vector<Entry_Type>::iterator ElemIterator;
	ElemIterator CurrentElem;	//iterator to current element

};
}
#endif