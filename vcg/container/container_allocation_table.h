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
Revision 1.9  2006/12/03 18:01:01  ganovelli
versione compliant vs2005

Revision 1.8  2006/06/08 20:28:57  ganovelli
aggiunto qualche const sui parametri

Revision 1.7  2005/10/15 16:21:48  ganovelli
Working release (compilata solo su MSVC), vector_occ � migrato da component_opt

Revision 1.6  2005/07/07 13:33:51  ganovelli
some comment

Revision 1.5  2005/07/06 15:28:10  ganovelli
aggiornamento di alcuni path

Revision 1.4  2004/04/05 18:20:50  ganovelli
Aggiunto typename
Eliminata bug di ricorsione nell'istanzazione dei template

Revision 1.3  2004/03/31 22:36:44  ganovelli
First Working Release (with this comment)


****************************************************************************/
  
#ifndef __VCGLIB_CAT__
#define __VCGLIB_CAT__

#include <vector>
#include <list>
#include <algorithm>
#include <assert.h>
#include <vcg/container/entries_allocation_table.h>


namespace vcg {
	/*@{*/
/*!
 * CATBase is the abstract class for all the allocation tables. These table keep track of
 * where the traced vector (see traced_ector.h) are kept in memory.
 * The goal is to know (given a pointer to a memory location), which is the vector the pointed 
 * element is in
 * 
 */

template <typename STL_CONT>
class CATBase{
public:
typedef  typename STL_CONT::value_type ValueType;

virtual void Resort(ValueType*,ValueType*) =0;
virtual void Remove(const STL_CONT&) = 0;
virtual void AddDataElem(ValueType*,int)=0;
virtual void Resize(ValueType*,int)=0;

public:
// ID serves as a type trait. 
static int & Id(){
			static int id=0;
			return id;
			}	
};

/// CATEntry: first derivation templated on the type of entry
/// It implements all the methods to trace and access vector_occ's elements
template <typename STL_CONT, class ENTRY_TYPE>
class CATEntry: public CATBase<STL_CONT>{
public:
typedef  typename STL_CONT::value_type ValueType;
typedef  ENTRY_TYPE EntryType;

CATEntry(){if(Id()==0){
							Id() = CATBase<STL_CONT>::Id()+1;
							CATBase<STL_CONT>::Id() = Id();
							}
					}


static unsigned int Ord(const ValueType *);
static ENTRY_TYPE & GetEntry(typename STL_CONT::value_type*pt);

static	void  Insert( STL_CONT & c,bool cond=false );				// insert a vector to trace
virtual void	Remove(  const STL_CONT  &	c);								// remove the container c
static  void	RemoveIfEmpty(  const STL_CONT  &	c);					// remove the container c
static  void	Remove(  ValueType  *	v);										  // remove the container that contains v

virtual void Resort(	ValueType* old_start,			// resort the allocation table
											ValueType* new_start);			// after a container was moved
	
protected:

static std::list<ENTRY_TYPE >& AT(){							// tallocation table
static std::list<ENTRY_TYPE > allocation_table; 
		return allocation_table;
	}
static bool & UTD(){
	static bool upToDate;														// true if Lower() and  Upper() 
	return upToDate;																// are up to date
	}

static ValueType *& Lower() {
	static ValueType * lower;											// pointer to the first element
	return lower;																		// of the last container accessed
	}
static ValueType *& Upper() {
	static ValueType * upper;											// pointer to the first element
	return upper;																		// if the container next to the last accessed
}		

static typename std::list<ENTRY_TYPE>::iterator	 & Curr(){		// container that was last accessed
	static typename std::list<ENTRY_TYPE>::iterator currEntry;
	return currEntry;
}


static bool IsTheSameAsLast(const ValueType *pt);	// true if pt is in the  container
																							// that was accessed last
static void Update(const ValueType*);							// set Upper() e Lower() 
static typename std::list<ENTRY_TYPE>::iterator FindBase(const ValueType * pt);	
																							// find the container that contains pt (naive)
virtual  void  AddDataElem(typename STL_CONT::value_type * pt,int n);// add n element to the auxiliary data
virtual  void  Resize(typename STL_CONT::value_type * pt,int n);// resize the  auxiliary data

public:
static int & Id(){															// unique identifier of the istance
		static int id=0;														// (used as type trait)
		return id;
		}
};

// --------------------------- CATEntry: implementation --------------------
// derivazione fatta per i membri Occ (Optional Component Compact)
template <typename STL_CONT, class ENTRY_TYPE>
unsigned int CATEntry<STL_CONT,ENTRY_TYPE>::

Ord(const ValueType * pt)
{
	Update(pt);
	return (pt-Lower());
}


template <typename STL_CONT, class ENTRY_TYPE>
typename std::list<ENTRY_TYPE>::iterator CATEntry<STL_CONT,ENTRY_TYPE>::

FindBase(const ValueType * pt)
{
typename std::list<ENTRY_TYPE>::iterator ite,curr_base,_;
ite = AT().begin();
curr_base = AT().end();

for(;ite != AT().end();ite++)
	if( pt < (*ite).Start())
		return curr_base;
	else
		curr_base = ite;

return curr_base;
}


template <typename STL_CONT, class ENTRY_TYPE>
 bool CATEntry< STL_CONT, ENTRY_TYPE>::
 
IsTheSameAsLast(const ValueType * pt)
{
return ( UTD() && ( !(Lower()> pt)) && (pt < Upper()) );
}

template <typename STL_CONT, class ENTRY_TYPE>
void CATEntry< STL_CONT, ENTRY_TYPE>::

Update(const ValueType * pt)
{
if(!IsTheSameAsLast(pt)){
	typename std::list<ENTRY_TYPE>::iterator lower_ite,upper_ite;
	lower_ite = FindBase(pt);

	assert(	lower_ite!=AT().end());

	Lower() = (*lower_ite).Start();
	if( (*lower_ite).Start() == AT().back().Start())
		Upper() = (ValueType *) 0xffffffff;
	else
	{
		upper_ite = lower_ite;	++upper_ite;	
		Upper() = (*upper_ite).Start();
	}
	
	Curr() = lower_ite;
	UTD() = true;
	}
}

template <typename STL_CONT, class ENTRY_TYPE>
void CATEntry< STL_CONT,  ENTRY_TYPE>::
Resort(ValueType* old_start,ValueType* new_start)
{
AT().sort();
UTD() = false;
}

template <typename STL_CONT, class ENTRY_TYPE>
void CATEntry<STL_CONT, ENTRY_TYPE>::

Remove( const STL_CONT & c )
{
typename std::list<ENTRY_TYPE>::iterator ite;
for(ite = AT().begin(); ite != AT().end();  ++ite)
	if((*ite).C() == &c)
		{
			AT().erase(ite);
			break;
		}
UTD() = false;
}

template <typename STL_CONT, class ENTRY_TYPE>
void CATEntry<STL_CONT, ENTRY_TYPE>::

RemoveIfEmpty( const STL_CONT & c )
{
typename std::list<ENTRY_TYPE>::iterator ite;
for(ite = AT().begin(); ite != AT().end();  ++ite)
	if((*ite).C() == &c)
			if(!(*ite).Empty())
				AT().erase(ite);
UTD() = false;
}

template <typename STL_CONT, class ENTRY_TYPE>
void CATEntry<STL_CONT, ENTRY_TYPE>::

Remove(ValueType  *	pt)
{
	typename std::list<ENTRY_TYPE>::iterator lower_ite;
	lower_ite = FindBase(pt);
	AT().erase(lower_ite);
	UTD() = false;
		
}

template <typename STL_CONT, class ENTRY_TYPE>
void CATEntry<STL_CONT, ENTRY_TYPE>::

Insert( STL_CONT & c,bool cond )
{
ENTRY_TYPE entry(c);
typename std::list<ENTRY_TYPE>::iterator lower_ite,upper_ite;
upper_ite = FindBase( c.Pointer2begin());
bool isIn = (upper_ite != AT().end());
if(isIn){
	if((*upper_ite).C() != &c )
				++upper_ite; 
	else
		return;
	}
lower_ite = AT().insert(upper_ite,entry);
lower_ite->Reserve(c.capacity());
lower_ite->Resize(c.size());
UTD() = false;
}

template <typename STL_CONT, class ENTRY_TYPE>
ENTRY_TYPE & CATEntry<STL_CONT, ENTRY_TYPE>::
GetEntry(typename STL_CONT::value_type*pt){
Update(pt);
return *Curr();
}

template <typename STL_CONT, class ENTRY_TYPE>
void CATEntry<STL_CONT, ENTRY_TYPE>::

AddDataElem(typename STL_CONT::value_type * pt,int n)
{
Update(pt);
Curr()->Push_back(n);
}


template <typename STL_CONT, class ENTRY_TYPE>
void CATEntry<STL_CONT, ENTRY_TYPE>::

Resize(typename STL_CONT::value_type * pt,int n)
{
Update(pt);
Curr()->Resize(n);
}


//--------------------------------------------------------------------------------------------
template <typename STL_CONT,class ATTR_TYPE>
class CAT:public CATEntry<STL_CONT, EntryCAT<STL_CONT,ATTR_TYPE> >{
typedef typename STL_CONT::value_type ValueType;
typedef CATEntry<STL_CONT, EntryCAT<STL_CONT,ATTR_TYPE> > TT;
public:
static ATTR_TYPE & Get(const ValueType * pt);
static CAT<STL_CONT,ATTR_TYPE> * New();
static CAT<STL_CONT,ATTR_TYPE> *& Instance(){ static CAT<STL_CONT,ATTR_TYPE> *  instance=NULL; return instance;}
}; 
//---------------------- CAT: implementation---------------------------------------------------
template <typename STL_CONT, class ATTR_TYPE>
ATTR_TYPE & CAT<STL_CONT,ATTR_TYPE>::

Get(const ValueType * pt)
{
int ord = Ord(pt);
//int ord = pt-  &(*  ((*AT().begin()).C()->begin())); se AT() contiene un solo elemento funziona anche cos
return TT::Curr()->Data()[ord];
}

template <typename STL_CONT, class ATTR_TYPE>
CAT<STL_CONT,ATTR_TYPE> * CAT<STL_CONT,ATTR_TYPE>::

New(){
	if(Instance()==NULL) 
		{
		 Instance() =  new CAT<STL_CONT,ATTR_TYPE>();
		}
	return Instance();
	}


//---------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// TempData implements a handle to one of the vector od data stored in EntryCATMulti
template <class STL_CONT, class ATTR_TYPE>
class TempData{
public:
	TempData(std::vector<ATTR_TYPE>  *d):item(d){};
		typedef ATTR_TYPE attr_type;

		std::vector<ATTR_TYPE>  * Item(){return item;};
		std::vector<ATTR_TYPE>  * item;
		ATTR_TYPE & operator []( typename STL_CONT::value_type * v)
			{
				int pos = CATEntry<STL_CONT, EntryCATMulti<STL_CONT> >::Ord(v);
				return (*item)[pos];
			}
	};
//----------------------------------------------------------------------------------

};//end namespace vcg

#endif
