#ifndef __PARTIALCONT__
#define __PARTIALCONT__

#include<vector>

template <class STL_CONT, class ELEM> 
struct Partial_Container : STL_CONT
{
	
	
	typedef typename STL_CONT::iterator ite_father;

	typedef ELEM value_type;

public:
struct iterator{

		ite_father i;
		iterator (){}
		iterator (ite_father i_):i(i_){}

		ELEM &operator *(){return *(*i);}

		void operator ++()
		{
			///((i!=(STL_CONT::end()))&& da controllare la fine
			//while ((*i)->IsInvalid())++i;
			++i;
		}

		iterator operator =(const iterator & oth){
				i=oth.i;
				return *this;
			}

		bool operator ==(const iterator & oth){
			return (i==oth.i);
			}
		bool operator !=(const iterator & oth){
			return (i!=oth.i);
			}

		bool operator <(const iterator & oth){
			return (i<oth.i);
			}
		};
	
	Partial_Container(){};
	iterator  begin(){return iterator(STL_CONT::begin());}
	iterator  end(){return iterator(STL_CONT::end());}

};
#endif