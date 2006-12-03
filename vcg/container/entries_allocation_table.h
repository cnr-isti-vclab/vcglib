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
Revision 1.5  2005/07/06 15:28:11  ganovelli
aggiornamento di alcuni path

Revision 1.4  2004/04/05 13:53:37  ganovelli
Aggiunto typename

Revision 1.3  2004/03/31 22:36:44  ganovelli
First Working Release (with this comment)


/****************************************************************************/ 
#ifndef __VCGLIB_ENTRIES__
#define __VCGLIB_ENTRIES__

namespace vcg {

// EntryCATBase: base class for the entry of the allocation table
// templated over the container type
template <class STL_CONT>
struct EntryCATBase{

EntryCATBase(STL_CONT & _c):c(_c){};
typename typename typename STL_CONT::value_type * Start() const;
virtual bool Empty(){return true;};
const STL_CONT *  C();
virtual void Push_back(const int &){};

virtual void Reserve(const int & s){};
virtual void Resize(const int & s){};

const bool operator < (const EntryCATBase<STL_CONT> & other) const;

private:
	STL_CONT & c;
};

//EntryCAT: entry for the case of optional core types (matches with CAT)
template <class STL_CONT,class ATTR_TYPE >
struct EntryCAT: public EntryCATBase<STL_CONT>{
typedef ATTR_TYPE attr_type;
EntryCAT(STL_CONT & _c) : EntryCATBase<STL_CONT>(_c){};
std::vector<ATTR_TYPE> & Data(){return data;}
void Push_back(const int & n){ for(int i = 0; i < n ; ++i) data.push_back(ATTR_TYPE());}
virtual void Reserve(const int & s){data.reserve(s);};
virtual void Resize(const int & s){data.resize(s);};


private:
std::vector<ATTR_TYPE> data;
};

//----------------------EntryCAT: implementation ----------------------------------------
template <class STL_CONT>
const bool EntryCATBase<STL_CONT>:: operator < (const EntryCATBase<STL_CONT> & other) const{
	return (Start() < other.Start());
}

template <class STL_CONT>
 typename typename STL_CONT::value_type  * EntryCATBase<STL_CONT>::Start()const {
	return  c.Pointer2begin();
	}

template <class STL_CONT>
	const STL_CONT * EntryCATBase<STL_CONT>::C(){
	return &c;
	}



// -----------------------------------------------------------------------------------------
// WrapBase: used to implement a list of pointers to std::vector of different types
// Wrap: derived from WrapBase (to take the function and from std::vector)
struct WrapBase{
virtual void Push_back(const int & n)=0;
virtual void Reserve(const int & n)=0;
virtual void Resize(const int & n)=0;	
	};
// (note) double hineritance is not necessary, just handy
template <class ATTR_TYPE>
struct Wrap: public WrapBase,std::vector<ATTR_TYPE>{
	virtual void Push_back(const int & n){for (int i = 0 ; i < n;	++i) push_back(		ATTR_TYPE());}	
	virtual void Reserve(const int & n){reserve(n);}
	virtual void Resize(const int & n){resize(n);}	
	};
//-------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------
// EntryCATMulti: entry type for multiple user data
template <class STL_CONT>
struct EntryCATMulti: public EntryCATBase<STL_CONT>{

EntryCATMulti(STL_CONT & _c) : EntryCATBase<STL_CONT>(_c){};
std::list<WrapBase * > & Data(){return data;}
void push_back(const int & n ){
	std::list< void * >::iterator ite;
	for(ite = data.begin(); ite != data.end(); ++ite)
		(*ite)->Push_back(n);
	}
virtual bool Empty(){return data.empty();};

virtual void Reserve(const int & n){
																			std::list<WrapBase * >::iterator ite;
																			for(ite = data.begin(); ite != data.end(); ++ite)
																				(*ite)->Reserve(n);
																		};
virtual void Resize(const int & n){		
																			std::list<WrapBase * >::iterator ite;
																			for(ite = data.begin(); ite != data.end(); ++ite)
																				(*ite)->Resize(n);
																	};

private:
	std::list< WrapBase * > data;
};
//----------------------------------------------------------------------------------


//----------------------------------------------------------------------------------
// TempData implements a handle to one of the vector od data stored in EntryCATMulti
template <class STL_CONT, class ATTR_TYPE>
class TempData{
public:
	TempData(std::vector<ATTR_TYPE>  *d):item(d){};
		typedef typename ATTR_TYPE attr_type;

		std::vector<ATTR_TYPE>  * Item(){return item;};
		std::vector<ATTR_TYPE>  * item;
		ATTR_TYPE & operator [](typename typename STL_CONT::value_type * v)
			{
				int pos = CATEntry<STL_CONT, EntryCATMulti<STL_CONT> >::Ord(v);
				return (*item)[pos];
			}
	};
//----------------------------------------------------------------------------------


}; // end namespace vcg

#endif