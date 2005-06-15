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
Revision 1.6  2005/06/01 13:47:59  pietroni
resolved hash code conflicts

Revision 1.5  2005/03/15 09:50:44  ganovelli
there was a debug line, now removed

Revision 1.4  2005/03/14 15:11:18  ganovelli
ClosestK added and other minor changes

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


		//type of entries element of a cell
		typedef typename std::pair<ElemType*,int> EntryType;

		//This Class Identify the cell
		struct HElement
		{
		private:

			///the min and max point corresponding to the cell - used to inverse hashing
			CoordType min;//coordinate min of the cell
			CoordType max;//coordinate max of the cell
			Point3i cell_n;//cell number
			
		public:

			//elements 
			std::map<ElemType*,int> _entries;
			
			//iterator to the map element into the cell
			typedef typename std::map<ElemType*,int>::iterator IteMap;

			HElement()
			{}

			HElement(ElemType* sim,const int  &_tempMark,CoordType _min,CoordType _max,Point3i _cell)
			{
				_entries.insert(EntryType(sim,_tempMark));
				min=_min;
				max=_max;
				assert(min<max);
				cell_n=_cell;
			}

			///return true if the element is in the cell
			bool IsIn(ElemType* sim)
			{int n=elem.count(sim);
			return (n==1);
			}

			///return the number of elements stored in the cell
			int Size()
			{return (int)(_entries.size());}

			///update or insert an element into a cell
			void Update(ElemType* sim, const int & _tempMark)
			{
				IteMap I=_entries.find(sim);

				if (I!=_entries.end())//the entry exist in the cell
					(*I).second=_tempMark;
				else
					_entries.insert(_entries.begin(),EntryType(sim,_tempMark));

				//at the end update the temporary mark on the simplex
				sim->Mark()=_tempMark;
				//Assert();
			}

			///given an iterator to the instance of the entry in the cell
			///return true if the the entry is valid
			///(using temporary mark).
			bool IsUpdated(IteMap &I)
			{
				return ((*I).second >= (*I).first->Mark());
			}

			///given an simplex pointer
			///return true if the the entry corripondent to that 
			///simplex is valid or not
			///(using temporary mark).
			bool IsUpdated(ElemType* sim)
			{
				IteMap I=_entries.find(sim);
				if (I!=_entries.end())
					return(IsUpdated(I));
				else
					return false;
			}

			//add to the vector all simplexes of the map that have a right timestamp or are not deleted
			void  Elems(std::vector<ElemType*> & res)
			{
				for (IteMap ite=_entries.begin();ite!=_entries.end();ite++)
				{
					ElemType* sim=(*ite).first;
					if (IsUpdated(ite)&&(!sim->IsD()))
						res.push_back(sim);
				}
			}

			
			CoordType Min()
			{return min;}
			
			CoordType Max()
			{return max;}

			Point3i CellN()
			{return cell_n;}

			bool operator ==(const HElement &h)  
			{return (cell_n==h.CellN());}

			bool operator !=(const HElement &h)  
			{return ((cell_n!=h.CellN()));}
			
			void Assert()
			{
				for (IteMap ite=_entries.begin();ite!=_entries.end();ite++)
				{
					ElemType* sim=(*ite).first;
					if (IsUpdated(sim))
					{
						ScalarType Xs=sim->P().X();
						ScalarType Ys=sim->P().Y();
						ScalarType Zs=sim->P().Z();
						ScalarType Xm=Min().X()-0.2;
						ScalarType Ym=Min().Y()-0.2;
						ScalarType Zm=Min().Z()-0.2;
						ScalarType XM=Max().X()+0.2;
						ScalarType YM=Max().Y()+0.2;
						ScalarType ZM=Max().Z()+0.2;
						if ((Xs<Xm)||(Xs>XM)||(Ys<Ym)||(Ys>YM)||(Zs<Zm)||(Zs>ZM))
						{
							printf("---ERROR---\n");
							printf("Point=%f,%f,%f.\n",Xs,Ys,Zs);
							printf("cellMin=%f,%f,%f\n",Xm,Ym,Zm);
							printf("cellMax=%f,%f,%f\n",XM,YM,ZM);
					
						}
					}
				}
			}

		}; // end struct HElement

		///this Iterator returns all the elements that
		///are in a specified box.
		struct ClosersIterator{

		private:
			CoordType p;
			SpatialHashTable<ElemType> * sh;	///pointer to spatial hash table structure
			vcg::Point3i mincorner,maxcorner;	///corners of the box where the scannig is performed
			vcg::Point3i curr_ic;				/// triple corresponding to the cell coordinate
			HElement * curr_c;				    /// current cell pointer
			typename HElement::IteMap curr_i;   /// current iterator inside the cell
			bool end;							///return true if the scanning of the elements is terminated
			
			///advance the current coordinate of one step inside the space box
			///set and to true if scannig fo cells is complete
			void  Advance(){
				if(curr_ic[0] < maxcorner[0]) 
					++curr_ic[0]; 
				else{
					if(curr_ic[1] < maxcorner[1]) 
						++curr_ic[1]; 
					else{
						if(curr_ic[2] < maxcorner[2])	
							++curr_ic[2];
						else
							end=true;
						curr_ic[1] = mincorner[1];
					}
					curr_ic[0] = mincorner[0];
				}
				
			}

			//operator Next
			//go to next simplex without considering temporary mark
			void Next()  {
				SpatialHashTable<ElemType>::IteHtable I;

				HElement::IteMap e = curr_c->_entries.end();
				--e;

				///if the current index
				//is not at the end of element in the cell so advance of one step
				if(curr_i != e) 
					++curr_i;
				else{
					///Advance until find a cell that isn't empty or the scan is complete
					Advance();

					while((!End())&&(sh->numElemCell(curr_ic,I)==0))
					{Advance();}

					if (!End())
					{
							curr_c = &((*I).second);
							curr_i = curr_c->_entries.begin();
					}
				}
			}

		public:

			void AssertUpdated()
			{
				ElemType* sim=(*curr_i).first;
				ScalarType Xs=sim->P().X();
				ScalarType Ys=sim->P().Y();
				ScalarType Zs=sim->P().Z();
				ScalarType Xm=curr_c->Min().X()-0.2;
				ScalarType Ym=curr_c->Min().Y()-0.2;
				ScalarType Zm=curr_c->Min().Z()-0.2;
				ScalarType XM=curr_c->Max().X()+0.2;
				ScalarType YM=curr_c->Max().Y()+0.2;
				ScalarType ZM=curr_c->Max().Z()+0.2;
				if ((Xs<Xm)||(Xs>XM)||(Ys<Ym)||(Ys>YM)||(Zs<Zm)||(Zs>ZM))
				{
					printf("---ERROR---\n");
					printf("Point=%f,%f,%f.\n",Xs,Ys,Zs);
					printf("cellMin=%f,%f,%f\n",Xm,Ym,Zm);
					printf("cellMax=%f,%f,%f\n",XM,YM,ZM);
					printf("tempMark Simplex=%d\n",sim->Mark());
					printf("tempMark Global=%d\n",sh->tempMark);
					if (curr_c->IsUpdated(sim))
						printf("updated");
					else
						printf("disupdated");
				}
			}

			///Initialize the iterator, p is the center of the box and _edge is his size
			void Init(SpatialHashTable<ElemType> * _sh, CoordType _p, const ScalarType  &_edge)
			{
				end=false;
				sh = _sh;
				p =_p;
				SpatialHashTable<ElemType>::IteHtable I;

				///find the box 
				CoordType halfDiag(_edge,_edge,_edge);
				mincorner = sh->PointToCell(p-halfDiag);
				maxcorner = sh->PointToCell(p+halfDiag);

				///set the initial cell of the iterator as the one that stay on the
				/// left lower position
				curr_ic = mincorner;

				//if the fist position isn't empty
				if (sh->numElemCell(curr_ic,I)>0)
				{
					curr_c = &((*I).second);
					curr_i = curr_c->_entries.begin();
				}
				else
				{///advance until don't find an non empty cell
				while((!End())&&(sh->numElemCell(curr_ic,I)==0))
					Advance();
				}
				///then if is not finished find the first updated occorrency
				if (!End())
				{
					curr_c = &((*I).second);
					curr_i = curr_c->_entries.begin();
					while ((!curr_c->IsUpdated(curr_i))&&(!end))
						Next();
				}
			}

			/////operator ++ of the itearator
			//void Next()  {
			//	bool isempty = true;
			//	HElement::IteMap e = curr_c->_entries.end();
			//	--e;

			//	///if the current index
			//	//is not at the end of elemnt in the cell so advance of one step
			//	if(curr_i != e) 
			//		++curr_i;

			//	///if the elemants on the cell are terminated then switch to another cell
			//	if(curr_i == e)
			//	{
			//		///Advance until find a cell that isn't empty
			//		while( Advance() && (isempty=sh->IsEmptyCell(curr_ic)));
			//		if(!isempty){
			//			curr_c = &(*(sh->hash_table.find(sh->Hash(curr_ic)))).second;
			//			curr_i =  curr_c->_entries.begin();
			//		}
			//		else
			//			end = true;
			//	}
			//}
			
			

			///operator ++ of the itearator
			void operator ++()  
			{
				Next();
				while ((!curr_c->IsUpdated(curr_i))&&(!end))
					Next();
				
				assert(curr_c->IsUpdated(curr_i)||(end));
			}

			///dereferent operator
			ElemType * operator *(){
				return (*curr_i).first;
			}

			///return true if the scanning of elements is complete
			bool End(){
				return end;
			}

		}; // end struct CloserIterator


		//hash table definition
		//typedef typename STDEXT::hash_map<int,HElement> Htable;
		typedef typename STDEXT::hash_multimap<int,HElement> Htable;
		//hash table definition
		//typedef typename STDEXT::hash_map<vcg::Point3i,HElement> Htable;
		//record of the hash table
		typedef typename std::pair<int,HElement> HRecord;
		//iterator to the hash table
		typedef typename Htable::iterator IteHtable;

		SpatialHashTable(){};
		~SpatialHashTable(){};

		int tempMark;

	protected:

		//ContSimplex & _simplex;
		
		Htable hash_table;

		int num;
		float l;


		CoordType min;
		CoordType max;

		int conflicts;
		/*Point3i min;
		Point3i max;*/

		///insert a new cell
		void _InsertNewHentry(ElemType* s,Point3i cell)
		{
			int h=Hash(cell);
			CoordType _min;
			CoordType _max;
			_min=CellToPoint(cell);
			_max=_min+CoordType(l,l,l);
			hash_table.insert(HRecord(h,HElement(s,tempMark,_min,_max,cell)));
			//Assert();
			s->Mark()=tempMark;
		}
		
		bool _IsInHtable(Point3i cell,IteHtable &result)
		{
			int h=Hash(cell);
			int count=hash_table.count(h);
			if (count==0)///in this case there is no entry for that key
				return false;
			else
			{
				////std::pair<Htable::const_iterator, Htable::const_iterator> p =hash_table.equal_range(h);
				std::pair<IteHtable, IteHtable> p =hash_table.equal_range(h);
				IteHtable i = p.first;

				while((i != p.second)&&((*i).second.CellN()!=cell))++i;

				if (i==p.second)///the scan is terminated and we have not fuond the cell
				{
					conflicts++;
					return false;
				}
				else	///we have found the right cell
				{
					result=i;
					return true;
				}
			}
		}

		///insert an element in a specified cell if the cell doesn't exist than
		///create it.
		void _InsertInCell(ElemType* s,Point3i cell)
		{
			IteHtable I;
			if (!_IsInHtable(cell,I))
				_InsertNewHentry(s,cell);
			else///there is the entry specified by the iterator I so update only the temporary mark
				(*I).second.Update(s,tempMark);

			//Assert();
		}
		
		// hashing
		const int Hash(Point3i p) const
		{
			//vcg::Point3i dim(100,100,100);
			return ((p.V(0)*P0 ^ p.V(1)*P1 ^ p.V(2)*P2)%num);
		}

		///return the cells intersected by the Bounding box of the simplex
		virtual void BoxCells(CoordType _min,CoordType _max,std::vector<Point3i>& ret)
		{
			ret.clear();
			Point3i MinI=PointToCell(_min);
			Point3i MaxI=PointToCell(_max);
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
			assert(ret.size()!=0);
		}
		
		/*void getAtCell(Point3i _c,std::vector<ElemType*> & res)
					{
						std::vector<ElemType> result;
						int h=Hash(_c);
						if (numElemCell(_c)!=0){
							IteHtable h_res=hash_table.find(h);
							((*h_res).second.Elems(tempMark,res));
						}
					}*/
		

				Point3i PointToCell(CoordType p) 
					{
						int x=(int)floor((p.V(0)-(ScalarType)min.V(0))/(ScalarType)l);
						int y=(int)floor((p.V(1)-(ScalarType)min.V(1))/(ScalarType)l);
						int z=(int)floor((p.V(2)-(ScalarType)min.V(2))/(ScalarType)l);
						return (vcg::Point3i(x,y,z));
					}

				CoordType CellToPoint(Point3i c)
					{
						ScalarType x=(((ScalarType)c.V(0)+min.V(0))*(ScalarType)l);
						ScalarType y=(((ScalarType)c.V(1)+min.V(1))*(ScalarType)l);
						ScalarType z=(((ScalarType)c.V(2)+min.V(2))*(ScalarType)l);
						return (CoordType(x,y,z));
					}
	public:

		///initialize the structure HashSpace is one estimation about 
		///how many keys the system have to generate in order to obtain as less
		///conflicts as possible
		void Init(CoordType _min,CoordType _max,ScalarType _l,int HashSpace=1000)
		{
			l=_l;
			min=_min;
			max=_max;
			num=HashSpace;
			tempMark=0;
			conflicts=0;
		}	

		virtual std::vector<Point3i> AddElem( ElemType* s)
		{
			std::vector<Point3i> box;
			BoxCells(s->BBox().min,s->BBox().max,box);
			for (std::vector<Point3i>::iterator bi=box.begin();bi<box.end();bi++)
				_InsertInCell(s,*bi);

			return box;
		}
		
		/*void AddElem( ElemType* s)
		{
			std::vector<Point3i> box;
			BoxCells(s->BBox().min,s->BBox().max,box);
			for (std::vector<Point3i>::iterator bi=box.begin();bi<box.end();bi++)
				_InsertInCell(s,*bi);
		}*/

		template<class ContElemType>
			void AddElems(  ContElemType & elem_set)
		{
			typename ContElemType::iterator i; 
			for(i = elem_set.begin(); i!= elem_set.end(); ++i)
				AddElem(&(*i));
		}

		

		////*********************************************************************
		//template <class A> 
		//	bool usefirst(const A & a,const A & b)const {return a.first < b.first;}

			int ClosestK(const int& k,ElemType* e, std::vector<ElemType*>& res) 
			{
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
					c = PointToCell(p);
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

				///return the simplexes on a specified cell
				void getAtCell(Point3i _c,std::vector<ElemType*> & res)
					{
						IteHtable I;
						if (_IsInHtable(_c,I))//if there is the cell then
							(*I).second.Elems(res);
					}

				// return the elem closer than radius
				int CloserThan(	typename ElemType::CoordType p, 
					typename ElemType::ScalarType radius, 
					std::vector<ElemType*> & closers){
						ClosersIterator cli;
						cli.Init(this,p,radius);
						while(!cli.End()){
							if ( (((*cli)->P() -p )*((*cli)->P() -p ) < radius*radius))// &&(*cli.curr_i).second >= tempMark)
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
						return PointToCell(min);
					}

					inline Point3i MaxCell()
					{
						return PointToCell(max);
					}

					/*inline int numElemCell(Point3i _c)
					{
						int h=Hash(_c);
						if (hash_table.count(h)==0)
							return 0;
						else 
						{
							IteHtable Ih=hash_table.find(h);
							return ((*Ih).second.Size());
						}
					}*/
					
					///return the number of elemnts in the cell and the iterator to the cell
					///if the cell exist
					int numElemCell(Point3i _c,IteHtable &I)
					{
						if (_IsInHtable(_c,I))
							return ((*I).second.Size());
						else
							return 0;
					}
					
					/*inline bool IsEmptyCell(Point3i _c)
					{
						return(numElemCell(_c)==0);
					}*/

					/*inline bool IsEmptyCell(Point3i _c)
					{
						int h=Hash(_c);
						if (hash_table.count(h)==0)
							return true;
						else
							return false;
					}*/

					///return the number of cell created
					int CellNumber()
					{return (hash_table.size());}
					
					int Conflicts()
					{return conflicts;}

					void Clear()
					{
						hash_table.clear();
					}
					
					void UpdateTmark()
					{tempMark++;}
					
					///only debug
					void DrawCell(HElement &c)
					{
						glPushMatrix();
						glTranslate(c.Min()+vcg::Point3d(l/2,l/2,l/2));
						glutWireCube(l);
						glPopMatrix();
					}

					void Draw()
					{
						for (IteHtable I=hash_table.begin();I!=hash_table.end();I++)
							DrawCell((*I).second);
					}
					
					///only debug
					void Assert()
					{
						for (IteHtable I=hash_table.begin();I!=hash_table.end();I++)
							(((*I).second).Assert());
					}
	}; // end class

}// end namespace

#undef P0
#undef P1
#undef P2

#endif
