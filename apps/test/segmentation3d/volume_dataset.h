
#include <qimage.h>
#include <qdir.h>
#include <qcolor.h>
#include <stdio.h>
#include <list>
#include <algorithm>

#define LimX 512
#define LimY 512
#define LimZ 240

#define dimXCell 20
#define dimYCell 20
#define dimZCell 20

#define TLBx 30
#define TLBy 30
#define TLBz 30


template <class ScalarType>
class Volume_Dataset_Optimized{

public:	

	Volume_Dataset_Optimized(){pFile=0;};
	~Volume_Dataset_Optimized(){};

private:
	class Cell
	{		
	public:

		void Clear()
		{
			for (int i=0;i<dimZCell;i++)
				for (int j=0;j<dimYCell;j++)
					for (int k=0;k<dimXCell;k++)
						Data[i][j][k]=0;
		}
		
		ScalarType Data[dimXCell][dimYCell][dimZCell] ;

		///operatorn to perform sorting
		inline bool operator <(const Cell &c)
		{
			return (c.timestamp>timestamp);
		}

		///operatorn to perform sorting
		inline bool operator ==(const Cell &c)
		{
			return (c.timestamp==timestamp);
		}
		Point3i index;
		int timestamp;
	};

	

	typedef typename std::list<Cell> StackType;
	typedef typename StackType::iterator IteStack;
	
	///the class of grid of iterator to stack structure (corrispondence grid - stack)
	struct TLBelem
	{
	private:
		IteStack StackPoint;

	public:
		
		inline IteStack & I()
		{
			return StackPoint;
		}

		bool inMem;
	};

	FILE * pFile;
	StackType CurrStack;
	
	TLBelem TLB[TLBx][TLBy][TLBz];
	Cell buffer[TLBx][TLBy] ;

	int n_element;
	int max_element;
	int timestamp;

	int lx;
	int ly;
	int lz;
	int TLBdx;
	int TLBdy;
	int TLBdz;

	int timesort;

	void SortStack()
	{
		/*std::sort<Cell>(CurrStack.begin(),CurrStack.end());
		timestamp=0;
		for (IteStack i=CurrStack.begin();i!=CurrStack.end();i++)
			(*i).timestamp=0;*/
	}

	///allocate momory for the buffer
	void InitBuffer(Cell _buffer[TLBx][TLBy])
	{
		for (int j=0;j<TLBy;j++)
			for (int i=0;i<TLBx;i++)
			_buffer[i][j].Clear();
	}

	///set total size of the dataset
	void Resize(int _lx,int _ly,int _lz)
	{
		lx=_lx;
		ly=_ly;
		lz=_lz;
		TLBdx=ceil(((float)_lx)/((float)dimXCell));
		TLBdy=ceil(((float)_ly)/((float)dimYCell));
		TLBdz=ceil(((float)_lz)/((float)dimZCell));
	}

	///return a pair made by cell coordinate an coordinate relative to the cell
	std::pair<vcg::Point3i,vcg::Point3i> GetCoordinate(int x,int y, int z)
	{
		int xd=(int) x/dimXCell;
		int yd=(int) y/dimYCell;
		int zd=(int) z/dimZCell;

		int xx= x%dimXCell;
		int yy= y%dimYCell;
		int zz= z%dimZCell;
		
		return (std::pair<vcg::Point3i,vcg::Point3i> (Point3i(xd,yd,zd),Point3i(xx,yy,zz)));
	}

	///add an element to the structure
	void Add(Cell _buffer[TLBx][TLBy],int x, int y, int z, ScalarType value)
	{
		assert(value<256);
		std::pair<vcg::Point3i,vcg::Point3i> coords=GetCoordinate(x,y,z);
		Point3i cell=coords.first;
		Point3i inter=coords.second;		
		_buffer[cell.V(0)][cell.V(1)].Data[inter.V(0)][inter.V(1)][inter.V(2)]=value;
	}

	///store a cell in the file
	void SaveCell(Cell &cell)
	{
		int size_unit=sizeof( ScalarType );
		int dim=dimXCell*dimYCell*dimZCell;
		int numwritten = fwrite(cell.Data, size_unit,dim,pFile );
	}

	///save to file the current buffer
	void SwapBuffer(Cell _buffer[TLBx][TLBy])
	{
		for (int x=0;x<TLBdx;x++)
			for (int y=0;y<TLBdy;y++)
				SaveCell(_buffer[x][y]);
	}
	
	///load a cell from file
	Cell& loadFromFile(Point3i p)
	{
		if (pFile!=0)
		{
			//8 for the 2 initial integer that describe the data set
			/*long offset=2*sizeof( ScalarType );*/
			long sizeCell=dimXCell*dimYCell*dimZCell;
			long offset=0;
			int index=(p.V(2)*TLBdx*TLBdy)+(p.V(1)*TLBdx)+p.V(0);
			offset+=sizeCell*(index)*sizeof( ScalarType );
			fseek(pFile,offset,SEEK_SET);
			Cell c;
			fread(c.Data, sizeof( ScalarType ),dimXCell*dimYCell*dimZCell,pFile);
			//assert(!ControlCell(c));
			c.index=p;
			return c;
		}
	}

	///funtion used to control correct values of a cell use only for debbugging
	bool ControlCell(Cell &c)
	{
		for (int i=0;i<dimXCell;i++)
				for (int j=0;j<dimYCell;j++)
					for (int k=0;k<dimZCell;k++)
							if (c.Data[i][j][k]>=256)
								return true;
		return false;
	}

public:

	///build and save the strucure made fo blocks
	void Resample(char *path,char *_newFile)
	{
		pFile=fopen (_newFile,"w+b");

		//std::vector<Cell> buffer;
		/*Cell buffer[TLBx][TLBy] ;*/

		//load first one image to see dimensions
		QImage qI=QImage();
		QDir Qd=QDir(path);
		QString qformat;
		QString Path=QString(path);
		Qd.setNameFilter("*.jpg");
		int levels=Qd.count();

		Qd.setSorting(QDir::Name);
		QString PathFile=Path;

		if (PathFile.right(1)!="/")
			PathFile.append("/");

		PathFile.append(Qd[0]);
		bool b=qI.load(PathFile,qformat);

		Resize(qI.width(),qI.height(),levels);

		InitBuffer(buffer);
		
		for (int z=0;z<levels; z++)
		{
			PathFile=Path;
			PathFile.append(Qd[z]);
			qformat=qI.imageFormat(PathFile);
			bool gray=qI.allGray();
			b=qI.load(PathFile,qformat);

			//first time I set limits of scalar field
			for (int x=0;x<qI.width();x++)
				for (int y=0;y<qI.height();y++)
				{
					QRgb color=qI.pixel(x,y);
					if (gray)//all tree component are the same 
						Add(buffer,x,y,z,qRed (color)); 
					else ///otherwise transform
						Add(buffer,x,y,z,qGray (color));
				}
			
			//if they are the last i swap the buffer
			if (((z%dimZCell)==(dimZCell-1))||(z==levels-1))
			{
				SwapBuffer(buffer);
				InitBuffer(buffer);
			}
		}
		fclose(pFile);
		pFile=0;
	}

	///set initial size of buffer in terms of cell
	void Init(int size,char* archive,int timeSort=3000)
	{
		timestamp=0;
		n_element=0;
		max_element=size;
		CurrStack.clear();
		pFile=fopen (archive,"r+b");
		timestamp=0;
		for(int z=0;z<TLBz;z++)
				for(int x=0;x<TLBx;x++)
					for(int y=0;y<TLBy;y++)
						TLB[x][y][z].inMem=false;
	}

	///control if i must sort the buffer
	bool TimeSort()
	{
	static clock_t time;
	clock_t elapsedsecs=abs(time-clock());
	if (elapsedsecs>timesort)
	{
		time=clock();
		return true;
	}
	return false;
	}
	
	inline Point3i Min()//to change in case of difefrent space between sections
	{return Point3i(0,0,0);}

	inline Point3i Max()//to change in case of difefrent space between sections
	{return Point3i(lx,ly,lz);}

	///return the x dimension of the dataset
	inline int dimX()
	{return lx;}

	///return the y dimension of the dataset
	inline int dimY()
	{return ly;}

	///return the z dimension of the dataset
	inline int dimZ()
	{return lz;}


	///return the lenght of the diagonal
	inline float Diag()
	{
		Point3f diag=Point3f((float) _X,(float) _Y,(float) _Z);
		return (diag.Norm());
	}

	/////erase the element not used for long time
	//void EraseLastUsed()
	//{
	//	if (TimeSort())
	//		SortStack();

	//	int mint=TLB[0][0][0].timestamp;
	//	Point3i minElem=Point3i(0,0,0);

	//	for(int z=0;z<TLBdz;z++)
	//		for(int x=0;x<TLBdx;x++)
	//			for(int y=0;y<TLBdy;y++)
	//				{
	//					if ((TLB[x][y][z].timestamp<mint)&&(TLB[x][y][z].inMem==true))
	//					{
	//						mint=TLB[x][y][z].timestamp;
	//						minElem=Point3i(x,y,z);
	//					}
	//				}

	//	TLB[minElem.V(0)][minElem.V(1)][minElem.V(2)].inMem=false;
	//	IteStack ite=TLB[minElem.V(0)][minElem.V(1)][minElem.V(2)].I();
	//	CurrStack.erase(ite);
	//}

	
	///return the value of the element in position point3i(i0,i1,i2)
	ScalarType getAt(Point3i p)
	{
		assert ((p.V(0)<dimX())&&(p.V(1)<dimY())&&(p.V(2)<dimZ()));
		assert ((p.V(0)>=0)&&(p.V(1)>=0)&&(p.V(2)>=0));

		std::pair<vcg::Point3i,vcg::Point3i> co=GetCoordinate(p.V(0),p.V(1),p.V(2));
		Point3i cell=co.first;
		Point3i inter=co.second;
		TLBelem e=TLB[cell.V(0)][cell.V(1)][cell.V(2)];
		if (e.inMem)///the element is in the buffer
		{
			IteStack i=e.I();
			ScalarType ret=(*i).Data[inter.V(1)][inter.V(0)][inter.V(2)];
			timestamp++;
			(*i).timestamp=timestamp;
			return (ret);
		}
		else///element fault, must load from file ///
		{
			//insert new element in the TLB table
			Cell c=loadFromFile(cell);
			CurrStack.push_front(c);
			TLB[cell.V(0)][cell.V(1)][cell.V(2)].I()=CurrStack.begin();
			TLB[cell.V(0)][cell.V(1)][cell.V(2)].inMem=true;
			(*CurrStack.begin()).timestamp=timestamp;
			n_element++;
			///if the number of element is the meximum , then erase one with second chanche algorithm
			if (n_element>=max_element)
			{
				if (TimeSort())
					SortStack();
				CurrStack.pop_back();
				n_element--;
			}
			ScalarType ret=(*CurrStack.begin()).Data[inter.V(1)][inter.V(0)][inter.V(2)];
			return ret;
		}
	}

	
	

};


template <class ScalarType>
class Volume_Dataset{

public:	

	Volume_Dataset(){};
	~Volume_Dataset(){};

	ScalarType Data[LimX][LimY][LimZ] ;
	
	int lx;
	int ly;
	int lz;
	
	
	///set total size of the dataset
	void Resize(int _lx,int _ly,int _lz)
	{
		lx=_lx;
		ly=_ly;
		lz=_lz;
	}

public:
	
	///build and save the strucure made fo blocks
	void Load(char *path)
	{
		//load first one image to see dimensions
		QImage qI=QImage();
		QDir Qd=QDir(path);
		QString qformat;
		QString Path=QString(path);
		Qd.setNameFilter("*.jpg");
		int levels=Qd.count();

		Qd.setSorting(QDir::Name);
		QString PathFile=Path;

		if (PathFile.right(1)!="/")
			PathFile.append("/");

		PathFile.append(Qd[0]);
		bool b=qI.load(PathFile,qformat);

		Resize(qI.width(),qI.height(),levels);

		
		for (int z=0;z<levels; z++)
		{
			PathFile=Path;
			PathFile.append(Qd[z]);
			qformat=qI.imageFormat(PathFile);
			bool gray=qI.allGray();
			b=qI.load(PathFile,qformat);

			//first time I set limits of scalar field
			for (int x=0;x<qI.width();x++)
				for (int y=0;y<qI.height();y++)
				{
					QRgb color=qI.pixel(x,y);
					if (gray)//all tree component are the same 
						Data[x][y][z]=qRed (color); 
					else ///otherwise transform
						Data[x][y][z]=qGray (color);
				}	
		}
		
	}

	
	inline Point3i Min()//to change in case of difefrent space between sections
	{return Point3i(0,0,0);}

	inline Point3i Max()//to change in case of different space between sections
	{return Point3i(lx,ly,lz);}

	///return the x dimension of the dataset
	inline int dimX()
	{return lx;}

	///return the y dimension of the dataset
	inline int dimY()
	{return ly;}

	///return the z dimension of the dataset
	inline int dimZ()
	{return lz;}

	///return the lenght of the diagonal
	inline float Diag()
	{
		Point3f diag=Point3f((float) _X,(float) _Y,(float) _Z);
		return (diag.Norm());
	}

	///return the value of the element in position point3i(i0,i1,i2)
	ScalarType getAt(Point3i p)
	{
		assert ((p.V(0)<dimX())&&(p.V(1)<dimY())&&(p.V(2)<dimZ()));
		assert ((p.V(0)>=0)&&(p.V(1)>=0)&&(p.V(2)>=0));
		return (Data[p.V(0)][p.V(1)][p.V(2)]);
	}

};
