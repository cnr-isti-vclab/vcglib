#ifndef SPATIAL_HASHING
#define SPATIAL_HASHING

//#define P0 73856093
//#define P1 19349663
//#define P2 83492791

#define P0 73
#define P1 19
#define P2 83

#include <map>
#include <hash_map>
#include <vector>

template <class SimplexType>

class SpatialHashTable{

public:
	
	typedef typename SimplexType* SimplexPointer;
	typedef typename SimplexType::CoordType Point3x;
	typedef typename Point3x::ScalarType ScalarType;
		
	
	//iterator to the map element into the cell
	typedef typename std::map<SimplexType*,int>::iterator IteMap;
	//element of a cell
	typedef typename std::pair<SimplexType*,int> MapCellElem;

	//element stored in the hash table
	struct Helement
		{
		
			std::map<SimplexType*,int> Elem;
			//int flag;

			public:
			
			Helement()
			{
			//	flag=0;
			}

			Helement(SimplexType* sim,int _tempMark)
			{
				Elem.insert(MapCellElem(sim,_tempMark));
			//	flag=0;
			}

			///return true if the element is in the cell
			bool IsIn(SimplexType* sim)
			{
				int n=Elem.count(sim);
				return (n==1);
			}
			
			int Size()
			{
				return (Elem.size());
			}
			
			///update or insert an element into a cell
			void Update(SimplexType* sim,int _tempMark)
			{
				std::pair<IteMap, bool> res=Elem.insert(MapCellElem(sim,_tempMark));
				//the element was already in the map structure so update the temporary mark
				if (res.second==false)
				{
					//update the temporary mark
					IteMap ite=res.first;
					(*ite).second=_tempMark;
				}
			}
			
			//return an array of all simplexes of the map that have a right timastamp or are not deleted
			std::vector<SimplexType*> Simplexes(int _tempMark)
			{
				std::vector<SimplexType*> result;
				result.clear();
				for (IteMap ite=Elem.begin();ite!=Elem.end();ite++)
				{
					SimplexType* sim=(*ite).first;
					int t=(*ite).second;
					if ((!sim->IsD())&&(t>=_tempMark))
						result.push_back(sim);
				}
				return (result);
			}
		};
	
	//enum { 
	//	// First user bit
	//	USER0      = 0x0001			// First user bit
	//		};

	//static int &LastBitFlag()
	//	{
	//		static int b =USER0;
	//		return b;
	//	}

	///// allocate a bit among the flags that can be used by user.
	//static inline int NewBitFlag()
	//	{
	//		LastBitFlag()=LastBitFlag()<<1;
	//		return LastBitFlag();
	//	}

	//// de-allocate a bit among the flags that can be used by user.
	//static inline bool DeleteBitFlag(int bitval)
	//	{	
	//		if(LastBitFlag()==bitval) {
	//				LastBitFlag()= LastBitFlag()>>1;
	//				return true;
	//		}
	//		assert(0);
	//		return false;
	//	}

	//hash table definition
	typedef typename std::hash_map<int,Helement> Htable;
	//record of the hash table
	typedef typename std::pair<int,Helement> HRecord;
	//iterator to the hash table
	typedef typename Htable::iterator IteHtable;


	//SpatialHashTable(ContSimplex & r_):_simplex(r_){};
	SpatialHashTable(){};
    ~SpatialHashTable(){};
	
	
	//ContSimplex & _simplex;
	int TempMark;
	Htable hash_table;
	
	int num;
	float l;
	

	Point3x min;
	Point3x max;


	void Init(Point3x _min,Point3x _max,ScalarType _l)
	{
		min=_min;
		max=_max;
		l=_l;
		Point3x d=max-min;
		num=floor(d.V(0)*d.V(1)*d.V(2)/l);
		TempMark=0;
	}

	void InsertInCell(SimplexType* s,Point3i cell)
	{
		int h=Hash(cell);
		//insert a cell if there isn't
		if (hash_table.count(h)==0)
				hash_table.insert(HRecord(h,Helement(s,TempMark)));
		//otherwise insert the element or update the temporary mark
		else
			{
				IteHtable HI=hash_table.find(h);
//				(*HI).second.flag|=_flag;
				(*HI).second.Update(s,TempMark);
			}
	}
	

	std::vector<Point3i> addBox(SimplexType* s,Point3x _min,Point3x _max)
	{
		
		std::vector<Point3i> box=BoxCells(_min,_max);
		for (std::vector<Point3i>::iterator bi=box.begin();bi<box.end();bi++)
			InsertInCell(s,*bi);
		return box;
	}

	std::vector<Point3i> BoxCells(Point3x _min,Point3x _max)
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

	std::vector<Point3i> Cells(SimplexType *s)
	{
		Point3x min=Point3x((*s).V(0)->P());
		Point3x max=min;

		//find min max coordinate of bounding box
		for (int i=1;i<3;i++)
		{
			if (min.V(0)>(*s).V(i)->P().V(0))
				min.V(0)=(*s).V(i)->P().V(0);
			if (min.V(1)>(*s).V(i)->P().V(1))
				min.V(1)=(*s).V(i)->P().V(1);
			if (min.V(2)>(*s).V(i)->P().V(2))
				min.V(2)=(*s).V(i)->P().V(2);

			if (max.V(0)<(*s).V(i)->P().V(0))
				max.V(0)=(*s).V(i)->P().V(0);
			if (max.V(1)<(*s).V(i)->P().V(1))
				max.V(1)=(*s).V(i)->P().V(1);
			if (max.V(2)>(*s).V(i)->P().V(2))
				max.V(2)=(*s).V(i)->P().V(2);
		}

		return BoxCells(s,min,max);
	}

	std::vector<Point3i> addSimplex(SimplexType *s)
	{
		Point3x min=Point3x((*s).V(0)->P());
		Point3x max=min;

		//find min max coordinate of bounding box
		for (int i=1;i<3;i++)
		{
			if (min.V(0)>(*s).V(i)->P().V(0))
				min.V(0)=(*s).V(i)->P().V(0);
			if (min.V(1)>(*s).V(i)->P().V(1))
				min.V(1)=(*s).V(i)->P().V(1);
			if (min.V(2)>(*s).V(i)->P().V(2))
				min.V(2)=(*s).V(i)->P().V(2);

			if (max.V(0)<(*s).V(i)->P().V(0))
				max.V(0)=(*s).V(i)->P().V(0);
			if (max.V(1)<(*s).V(i)->P().V(1))
				max.V(1)=(*s).V(i)->P().V(1);
			if (max.V(2)>(*s).V(i)->P().V(2))
				max.V(2)=(*s).V(i)->P().V(2);
		}

		///now set the cell touched by the bounding box of the triangle
		return (addBox(s,min,max));
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

	///refresh all the alement in the container
	/*void RefreshElements()
	{
		TempMark++;
		for (SimplexIterator si=_simplex.begin();si<_simplex.end();++si)
		{
			if (!(*si).IsD())
				addSimplex(&*si);
		}
	}*/

	///erase all element of the hash table where the timestamp is old
	//void Clean()
	//{
	//	for (IteHtable Ih=hash_table.begin();Ih!=hash_table.end();Ih++)
	//		//if after clean there are no element erase the cell
	//		if ((*Ih).second.Clean(TempMark))
	//			hash_table.erase(Ih);
	//}

	/*void Clean()
	{
		hash_table.clear();
		RefreshElements();
	}*/
	
	void Clear()
	{
		hash_table.clear();
	}

	/*int FlagCell(Point3i p)
	{
		int h=Hash(p);
		if (hash_table.count(h)==0)
			return 0;
		else 
		{
			IteHtable Ih=hash_table.find(h);
			return ((*Ih).second.Flag());
		}
	}*/

	std::vector<SimplexPointer> getAt(Point3x _p)
	{
		std::vector<SimplexPointer> result;
		Point3i c=Cell(p);
		return (getAtCell(c));
	}

	std::vector<SimplexPointer> getAtCell(Point3i _c)
	{
		std::vector<SimplexPointer> result;
		int h=Hash(_c);
		if (numElemCell(_c)==0)
		{
			return result;
		}
		else
		{
			IteHtable res=hash_table.find(h);
			return ((*res).second.Simplexes(TempMark));
		}
	}


private:

		Point3i Cell(Point3x p)
		{
			int x=floor(p.V(0)/l);
			int y=floor(p.V(1)/l);
			int z=floor(p.V(2)/l);
			return Point3i(x,y,z);
		}

		int Hash(Point3i p)
		{
			return ((p.V(0)*P0 ^ p.V(1)*P1 ^ p.V(2)*P2)%num);
		}

};
#endif