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
Revision 1.3  2005/03/11 15:25:29  ganovelli
added ClosersIterator and other minor changes. Not compatible with the previous version.
Still other modifications to do (temporary commit)

Revision 1.2  2005/02/21 12:13:25  ganovelli
added vcg header



****************************************************************************/

#ifndef VCGLIB_SPATIAL_HASHING
#define VCGLIB_SPATIAL_HASHING

#define P0 73856093
#define P1 19349663
#define P2 83492791

#include <map>
#include <vector>
#include <algorithm>
#ifdef WIN32
	#include <hash_map>
	#define STDEXT stdext
#else
	#include <ext/hash_map>
	#define STDEXT __gnu_cxx
#endif

namespace vcg{
/** Spatial Hash Table
     Spatial Hashing as described in
		 "Optimized Spatial Hashing for Collision Detection of Deformable Objects", 
		 Matthias Teschner and Bruno Heidelberger and Matthias Muller and Danat Pomeranets and Markus Gross
     */
template <class ElemType>
class SpatialHashTable{

public:
	
	typedef typename ElemType::CoordType CoordType;
	typedef typename CoordType::ScalarType ScalarType;
		
	
	//element of a cell
	typedef typename std::pair<ElemType*,int> MapCellElem;

	//element stored in the hash table
	struct HElement
		{

			//iterator to the map element into the cell
			typedef typename std::map<ElemType*,int>::iterator IteMap;

			std::map<ElemType*,int> elem;
			//int flag;

			public:
			
			HElement()
			{
			//	flag=0;
			}

			HElement(ElemType* sim,const int  &_tempMark)
			{
				elem.insert(MapCellElem(sim,_tempMark));
			//	flag=0;
			}

			///return true if the element is in the cell
			bool IsIn(ElemType* sim)
			{
				int n=elem.count(sim);
				return (n==1);
			}
			
			int Size()
			{
				return (int)(elem.size());
			}
			
			///update or insert an element into a cell
			void Update(ElemType* sim, const int & _tempMark)
			{
				std::pair<IteMap, bool> res=elem.insert(MapCellElem(sim,_tempMark));
				//the element was already in the map structure so update the temporary mark
				if (res.second==false)
				{
					//update the temporary mark
					IteMap ite=res.first;
					(*ite).second=_tempMark;
				}
			}
			
			//return an array of all simplexes of the map that have a right timestamp or are not deleted
			void  Elems(const int & _tempMark,std::vector<ElemType*> & res)
			{

				for (IteMap ite=elem.begin();ite!=elem.end();ite++)
				{
					ElemType* sim=(*ite).first;
					int t=(*ite).second;
					if ((!sim->IsD())&&(t>=_tempMark))
						res.push_back(sim);
				}
			}
		}; // end struct HElement
	

		struct ClosersIterator{
			CoordType p;
			SpatialHashTable<ElemType> * sh;
			vcg::Point3i mincorner,maxcorner;
			ScalarType sq_radius;

			// current position
			vcg::Point3i curr_ic; // triple corresponding to the cell
			HElement * curr_c;		// current cell
			typename HElement::IteMap curr_i; // current iterator
			bool end;

			bool  Advance(){
				if(curr_ic[0] < maxcorner[0]) ++curr_ic[0]; 
				else{
						if(curr_ic[1] < maxcorner[1]) ++curr_ic[1]; 
						else{
								if(curr_ic[2] < maxcorner[2])	++curr_ic[2];
								else
									return false;
								curr_ic[1] = mincorner[1];
						}
					curr_ic[0] = mincorner[0];
				}
				curr_c = &(*(sh->hash_table.find(sh->Hash(curr_ic)))).second;
				return true;
			}
			void Init(SpatialHashTable<ElemType> * _sh, CoordType _p, const ScalarType  &_radius)
			{
					sh = _sh;
					p =_p;
					CoordType halfDiag(_radius,_radius,_radius);
					mincorner = sh->Cell(p-halfDiag);
					maxcorner = sh->Cell(p+halfDiag);
					curr_ic = mincorner;
					sq_radius = _radius * _radius;

					IteHtable iht = sh->hash_table.find(sh->Hash(curr_ic));

					// initialize the iterator to the first element
					bool isempty  = (iht == sh->hash_table.end());
					if(isempty)
						while( Advance() && (isempty=sh->IsEmptyCell(curr_ic)));

					if(!isempty){
						curr_c = &(*(sh->hash_table.find(sh->Hash(curr_ic)))).second;
						curr_i =  curr_c->elem.begin();
						end = false;
					}
					else
						end = true;

			}

			void operator ++()  {
				bool isempty = true;
				HElement::IteMap e = curr_c->elem.end();
				--e;
				if(curr_i != e) 
					++curr_i;
				else{
					while( Advance() && (isempty=sh->IsEmptyCell(curr_ic)));
					if(!isempty){
						curr_c = &(*(sh->hash_table.find(sh->Hash(curr_ic)))).second;
						curr_i =  curr_c->elem.begin();
					}
					else
						end = true;
				}
			}
			 ElemType * operator *(){
				 vcg::Point3d __ = (*curr_i).first->P();
						return (*curr_i).first;
			}

			 bool End(){
				 //bool __  = (curr_i == curr_c->elem.end());
				 //return ( (curr_ic == maxcorner) && (curr_i == curr_c->elem.end()) );
				 return end;
			 }
		}; // end struct CloserIterator


	//hash table definition
	typedef typename STDEXT::hash_map<int,HElement> Htable;
	//record of the hash table
	typedef typename std::pair<int,HElement> HRecord;
	//iterator to the hash table
	typedef typename Htable::iterator IteHtable;

	SpatialHashTable(){};
    ~SpatialHashTable(){};
	
	
	//ContSimplex & _simplex;
	int tempMark;
	Htable hash_table;
	
	int num;
	float l;
	

	CoordType min;
	CoordType max;


	void Init(CoordType _min,CoordType _max,ScalarType _l)
	{
		min=_min;
		max=_max;
		l=_l;
		CoordType d=max-min;
		//num = (int) floor(d.V(0)*d.V(1)*d.V(2)/l);
		num = (int) floor(100*d.V(0)*d.V(1)*d.V(2)/l);
		tempMark=0;
	}

	void InsertInCell(ElemType* s,Point3i cell)
	{
		int h=Hash(cell);
		//insert a cell if there isn't
		if (hash_table.count(h)==0)
				hash_table.insert(HRecord(h,HElement(s,tempMark)));
		//otherwise insert the element or update the temporary mark
		else
			{
				IteHtable HI=hash_table.find(h);
//				(*HI).second.flag|=_flag;
				(*HI).second.Update(s,tempMark);
			}
	}
	

	std::vector<Point3i> AddElem( ElemType* s)
	{
		std::vector<Point3i> box=BoxCells(s->BBox().min,s->BBox().max);
		for (std::vector<Point3i>::iterator bi=box.begin();bi<box.end();bi++)
			InsertInCell(s,*bi);
		return box;
	}

	template<class ContElemType>
		void AddElems(  ContElemType & elem_set)
	{
		typename ContElemType::iterator i; 
		for(i = elem_set.begin(); i!= elem_set.end(); ++i)
			AddElem(&(*i));
	}

	std::vector<Point3i> BoxCells(CoordType _min,CoordType _max)
	{
		std::vector<Point3i> ret;
		ret.clear();
		Point3i MinI=Cell(_min);
		Point3i MaxI=Cell(_max);
		int dimx=abs(MaxI.V(0)-MinI.V(0));
		int dimy=abs(MaxI.V(1)-MinI.V(1));
		int dimz=abs(MaxI.V(2)-MinI.V(2));

		for (int x=0;x<=dimx;x++)
			for (int y=0;y<=dimy;y++)
				for (int z=0;z<=dimz;z++)
				{
					Point3i cell=Point3i(MinI.V(0)+x,MinI.V(1)+y,MinI.V(2)+z);
					ret.push_back(cell);
				}
		return ret;
	}
	
//*********************************************************************
	template <class A> 
		bool usefirst(const A & a,const A & b)const {return a.first < b.first;}

	int ClosestK(const int& k,	
		ElemType* e, 
	std::vector<ElemType*>& res) {

		typedef std::pair<ScalarType,ElemType*> ElemDist;
		std::vector<ElemDist > neigh_dist;
		std::vector<ElemDist >::iterator ite_nd;
		std::vector<ElemType* > neigh;
		std::vector<ElemType*>::iterator i_neigh;
		typename ElemType::CoordType p = e->P();
		ScalarType radius,tmp,d;

		// set the radius as the distance to the closest face
		radius =	p[2]-floor(p[2]/l)*l;
		if(radius > l*0.5) radius = l -radius;
		tmp =	p[1]-floor(p[1]/l)*l;
		if(tmp > l*0.5) tmp = l -tmp;
		if(radius > tmp) tmp = radius;
		tmp =	p[0]-floor(p[0]/l)*l;
		if(tmp > l*0.5) tmp = l -tmp;
		if(radius > tmp) radius = tmp;

		int x,y,z;
		vcg::Point3i mincorner,maxcorner,c;
		c = Cell(p);
		mincorner = maxcorner = c;	
		neigh_dist.push_back(ElemDist(-1,e));
		ite_nd = neigh_dist.begin();
		
		while((int)res.size() < k)
		{

		 //run on the border
			for( z = mincorner[2]; z <= maxcorner[2]; ++z)
				for(	y = mincorner[1]; y <= maxcorner[1];  ++y)
					for(	x = mincorner[0]; x <= maxcorner[0];)
					{

						neigh.clear();
						getAtCell(vcg::Point3i(x,y,z),neigh);
								for(i_neigh = 	neigh.begin(); i_neigh != neigh.end(); ++i_neigh)
								{
										d = Distance(p,(*i_neigh)->P());
										if( (*i_neigh) != e) 
											neigh_dist.push_back(ElemDist(d,*i_neigh));
									}
						if(
						( ( y == mincorner[1]) || ( y == maxcorner[1])) ||
						( ( z == mincorner[2]) || ( z == maxcorner[2])) ||
						( x == maxcorner[0])
						)++x; else x=maxcorner[0];
				}
//		,usefirst<ElemDist> ---<std::vector<ElemDist >::iterator >
			ite_nd =neigh_dist.begin();
			std::advance(ite_nd,res.size());
			std::sort(ite_nd,neigh_dist.end());
			while ( ( (int)res.size() < k ) && (ite_nd != neigh_dist.end()))
			{
				if((*ite_nd).first < radius)
					res.push_back( (*ite_nd).second );
				++ite_nd;
			}

			mincorner -= vcg::Point3i(1,1,1);
			maxcorner += vcg::Point3i(1,1,1);
			radius+=l;

		}
				return 0;	
	}
//**********************************************************************


	// return the elem closer than radius
	int CloserThan(	typename ElemType::CoordType p, 
									typename ElemType::ScalarType radius, 
									std::vector<ElemType*> & closers){
			ClosersIterator cli;
			cli.Init(this,p,radius);
			while(!cli.End()){
				if ( (((*cli)->P() -p )*((*cli)->P() -p ) < radius*radius) &&
					(*cli.curr_i).second >= tempMark)
				closers.push_back(*cli);
				++cli;
			}
	return (int)closers.size();
									}

	std::vector<Point3i> Cells(ElemType *s)
	{
		return BoxCells(s,s->BBox().min,s->BBox().max);
	}

	inline Point3i MinCell()
	{
		return Cell(min);
	}

	inline Point3i MaxCell()
	{
		return Cell(max);
	}

	inline int numElemCell(Point3i _c)
	{
		int h=Hash(_c);
		if (hash_table.count(h)==0)
			return 0;
		else 
			{
				IteHtable Ih=hash_table.find(h);
				return ((*Ih).second.Size());
			}
	}

	inline bool IsEmptyCell(Point3i _c)
	{
		int h=Hash(_c);
		if (hash_table.count(h)==0)
			return true;
		else
			return false;
	}
	
	
	void Clear()
	{
		hash_table.clear();
	}

	//void	std::vector<ElemType> getAt(CoordType _p,std::vector<ElemType> & res)
	//{
	//	std::vector<ElemType> result;
	//	Point3i c=Cell(_p);
	//	return (getAtCell(c,res));
	//}

	void getAtCell(Point3i _c,std::vector<ElemType*> & res)
	{
		std::vector<ElemType> result;
		int h=Hash(_c);
		if (numElemCell(_c)!=0){
			IteHtable h_res=hash_table.find(h);
			((*h_res).second.Elems(tempMark,res));
		}
	}

	const Point3i Cell(const CoordType & p) const 
		{
			int x=(int)floor(p.V(0)/l);
			int y=(int)floor(p.V(1)/l);
			int z=(int)floor(p.V(2)/l);
			return Point3i(x,y,z);
		}

		// hashing
		const int Hash(Point3i p) const
		{
			vcg::Point3i dim(100,100,100);
			return ((p.V(0)*P0 ^ p.V(1)*P1 ^ p.V(2)*P2)%num);
//			return ( p[2]-min[2] )* dim[0]*dim[1] +
//						 ( p[1]-min[1] )* dim[1] +
//						 ( p[0]-min[0] );
		}
private:
}; // end class

}// end namespace

#undef P0
#undef P1
#undef P2

#endif
