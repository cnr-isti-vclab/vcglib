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
	typedef typename SimplexType::CoordType CoordType;
	typedef typename CoordType::ScalarType ScalarType;
	typedef typename vcg::Box3<ScalarType> Box3x;

	typedef  DynamicSpatialHashTable<SimplexType,float> HashingTable;

	Collision_Detector(ContSimplex & r_):_simplex(r_){};
	~Collision_Detector(){};

	ContSimplex & _simplex;

	HashingTable *HTable;

	std::set<Point3i> vactive;

	int active;

	//control if two faces share an edge
	bool ShareEdge(SimplexType *f0,SimplexType *f1)
	{
		assert((!f0->IsD())&&(!f1->IsD()));
		for (int i=0;i<3;i++)
			if (f0->FFp(i)==f1)
				return (true);

		return(false);
	}

	///initialize the box for collision detection and the dimension of a cell
	void Init(CoordType _min,CoordType _max,ScalarType _l)
	{
		HTable=new HashingTable();
		Box3x bb(_min,_max);
		CoordType d=((_max-_min)/_l);
		vcg::Point3i dim;
		dim.Import<ScalarType>(d);
		HTable->InitEmpty(bb,dim);
	}

	//control if two faces share a vertex
	bool ShareVertex(SimplexType *f0,SimplexType *f1)
	{
		assert((!f0->IsD())&&(!f1->IsD()));
		for (int i=0;i<3;i++)
			for (int j=0;j<3;j++)
				if (f0->V(i)==f1->V(j))
					return (true);

		return(false);
	}

	//test real intersection between faces
	bool TestRealIntersection(SimplexType *f0,SimplexType *f1)
	{
		assert((!f0->IsD())&&(!f1->IsD()));
		if ((!f0->IsActive())&&(!f1->IsActive()))
			return false;
		//no adiacent faces
		assert(f0!=f1);
		if ((f0!=f1)&& (!ShareEdge(f0,f1))&&!ShareVertex(f0,f1))
		{
			//vcg::Segment3<ScalarType> segm;
			//bool copl=false;
			return (vcg::Intersection<SimplexType>((*f0),(*f1)));//,copl,segm))
				//return ((copl)||(segm.Length()>0.001));
		}
		return false;
	}

	///refresh all the elements of spatial hashing table 
	void RefreshElements()
	{
		HTable->Clear();
		vactive.clear();

		HTable->tempMark=0;

		for (SimplexIterator si=_simplex.begin();si<_simplex.end();++si)
		{
			if (!(*si).IsD())
			{
				(*si).HMark()=0;
				vcg::Box3i cells=HTable->Add(&*si);
				if ((*si).IsActive())
				{
					vcg::Box3i cells=HTable->Add(&*si);
					for (int x=cells.min.X(); x<=cells.max.X();x++)
						for (int y=cells.min.Y(); y<=cells.max.Y();y++)
							for (int z=cells.min.Z(); z<=cells.max.Z();z++)
								vactive.insert(vcg::Point3i(x,y,z));
				}
			}
		}
	}




	///put active cells on apposite structure
	template <class Container_Type>
		void UpdateStep(Container_Type &simplex)
	{
		vactive.clear();
		HTable->UpdateTmark();
		for (Container_Type::iterator si=simplex.begin();si<simplex.end();++si)
		{
			if ((!(*si).IsD())&&((*si).IsActive()))
			{
				vcg::Box3i cells=HTable->Add(&*si);
				for (int x=cells.min.X();x<=cells.max.X();x++)
					for (int y=cells.min.Y();y<=cells.max.Y();y++)
						for (int z=cells.min.Z();z<=cells.max.Z();z++)
							vactive.insert(vcg::Point3i(x,y,x));
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
			HashingTable::IteHtable I;
			if (HTable->numElemCell(p,I)>=2)
			{
				std::vector<SimplexType*> inCell;
				inCell.clear();
				HTable->getInCellUpdated(p,inCell);
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