#ifndef COLLISION_DETECTION
#define COLLISION_DETECTION

#include <set>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/space/intersection3.h>

template <class ContSimplex>
class Collision_Detector{

public:
	
		typedef typename ContSimplex::value_type SimplexType;
		typedef typename ContSimplex::value_type* SimplexPointer;
		typedef typename ContSimplex::iterator SimplexIterator;
		typedef typename SimplexType::CoordType Point3x;
		typedef typename Point3x::ScalarType ScalarType;

		typedef  SpatialHashTable<SimplexType> HashingTable;
		
		Collision_Detector(ContSimplex & r_):_simplex(r_){};
		~Collision_Detector(){};

		ContSimplex & _simplex;
		
		HashingTable *HTable;
		
		std::set<Point3i> vactive;

		int active;

	//control if two faces share an edge
	bool ShareEdge(SimplexType *f0,SimplexType *f1)
	{
		for (int i=0;i<3;i++)
				if (f0->FFp(i)==f1)
					return (true);

		return(false);
	}

	///initialize the box for collision detection and the dimension of a cell
	void Init(Point3x _min,Point3x _max,ScalarType _l)
	{
		HTable=new HashingTable();
		HTable->Init(_min,_max,_l);
	}

	//control if two faces share a vertex
	bool ShareVertex(SimplexType *f0,SimplexType *f1)
	{
		for (int i=0;i<3;i++)
			for (int j=0;j<3;j++)
				if (f0->V(i)==f1->V(j))
					return (true);

		return(false);
	}
	
	//test real intersection between faces
	bool TestRealIntersection(SimplexType *f0,SimplexType *f1)
	{
		
		if (((!f0->IsActive())&&(!f1->IsActive()))||(f0->IsD())||(f1->IsD()))
			return false;
		//no adiacent faces
		if ((f0!=f1)&& (!ShareEdge(f0,f1))
			&& (!ShareVertex(f0,f1)))
			return (vcg::Intersection<SimplexType>((*f0),(*f1)));
		return false;
	}

	///refresh all the elements of spatial hashing table 
	///this function must called sometimes
	void RefreshElements()
	{
		HTable->Clear();
		vactive.clear();///new
		for (SimplexIterator si=_simplex.begin();si<_simplex.end();++si)
		{
			if ((!(*si).IsD())&&(!(*si).IsActive()))
					HTable->addSimplex(&*si);
			///new now
			else
			{
				std::vector<Point3i> cells=HTable->addSimplex(&*si);
				for(std::vector<Point3i>::iterator it=cells.begin();it<cells.end();it++)
					vactive.insert(*it);
			}
			///end new now
		}

		//UpdateStep();	commented  now
	}
	
	/////put active cells on apposite structure
	//void UpdateStep()
	//{
	//	vactive.clear();
	//	for (SimplexIterator si=_simplex.begin();si<_simplex.end();++si)
	//	{
	//		if ((((!(*si).IsD()))&&(*si).IsActive()))
	//		{
	//			std::vector<Point3i> cells=HTable->addSimplex(&*si);
	//			for(std::vector<Point3i>::iterator it=cells.begin();it<cells.end();it++)
	//				vactive.insert(*it);
	//		}
	//	}
	//}


	///put active cells on apposite structure
	template <class Container_Type>
	void UpdateStep(Container_Type &simplex)
	{
		vactive.clear();
		for (Container_Type::iterator si=simplex.begin();si<simplex.end();++si)
		{
			if ((((!(*si).IsD()))&&(*si).IsActive()))
			{
				std::vector<Point3i> cells=HTable->addSimplex(&*si);
				for(std::vector<Point3i>::iterator it=cells.begin();it<cells.end();it++)
					vactive.insert(*it);
			}
		}
	}


	///control the real self intersection in the mesh and returns the elements that intersect with someone
	std::vector<SimplexType*> computeSelfIntersection()
	{
		std::vector<SimplexType*> ret;
		std::set<Point3i>::iterator act;
		for (act=vactive.begin();act!=vactive.end();act++)
				{
					Point3i p=*act;
					if (HTable->numElemCell(p)>=2)
					{
 						std::vector<SimplexType*> inCell=HTable->getAtCell(p);
						int nelem=inCell.size();
						if (nelem>=2)
						{
							//test combinations of elements
							for (int i=0;i<nelem-1;i++)
								for (int j=i+1;j<nelem;j++)
									if ((!inCell[i]->IsD())&&(!inCell[j]->IsD())&&(TestRealIntersection(inCell[i],inCell[j])))
									{
										ret.push_back(inCell[i]);
										ret.push_back(inCell[j]);
									}
						}
					}
				}
			return ret;
	}

};
#endif