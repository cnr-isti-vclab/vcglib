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


#ifndef __VCGLIB_CAT__
#define __VCGLIB_CAT__

#include <vector>
#include <list>
#include <algorithm>
#include <assert.h>
#include <vcg/container/entries_allocation_table.h>

namespace vcg {

// CATBase: abstract base class for all the allocation tables
template <class STL_CONT>
class CATBase{
public:
typedef STL_CONT::value_type VALUE_TYPE;

virtual void Resort(VALUE_TYPE*,VALUE_TYPE*) =0;
virtual void Remove(const STL_CONT&) = 0;
virtual void AddDataElem(VALUE_TYPE*,int)=0;

public:
// ID serves as a type trait. 
static int & Id(){
			static int id=0;
			return id;
			}	
};

// CATEntry: first derivation templated on the type of entry
// It implements all the methods to trace and access TVector element
template <class STL_CONT, class ENTRY_TYPE>
class CATEntry: public CATBase<STL_CONT>{
public:
typedef STL_CONT::value_type VALUE_TYPE;
typedef ENTRY_TYPE EntryType;

CATEntry(){if(Id()==0){
							Id() = CATBase<STL_CONT>::Id()+1;
							CATBase<STL_CONT>::Id() = Id();
							}
					}


static unsigned int Ord(VALUE_TYPE *);
static ENTRY_TYPE & GetEntry(STL_CONT::value_type*pt);

static	void  Insert( STL_CONT & c,bool cond=false );				// insert a vector to trace
virtual void	Remove(  const STL_CONT  &	c);								// remove the container c
static  void	RemoveIfEmpty(  const STL_CONT  &	c);								// remove the container c
static  void	Remove(  VALUE_TYPE  *	v);										// remove the container that contains v

virtual void Resort(	VALUE_TYPE* old_start,			// resort the allocation table
											VALUE_TYPE* new_start);			// after a container was moved
	
protected:

static std::list<ENTRY_TYPE >& AT(){							// tallocation table
static std::list<ENTRY_TYPE > allocation_table; 
		return allocation_table;
	}
static bool & UTD(){
	static bool upToDate;														// true if Lower() and  Upper() 
	return upToDate;																// are up to date
	}

static VALUE_TYPE *& Lower() {
	static VALUE_TYPE * lower;											// pointer to the first element
	return lower;																		// of the last container accessed
	}
static VALUE_TYPE *& Upper() {
	static VALUE_TYPE * upper;											// pointer to the first element
	return upper;																		// if the container next to the last accessed
}		

static std::list<ENTRY_TYPE>::iterator	 & Curr(){		// container that was last accessed
	static std::list<ENTRY_TYPE>::iterator currEntry;
	return currEntry;
}


static bool IsTheSameAsLast(VALUE_TYPE *pt);	// true if pt is in the  container
																							// that was accessed last
static void Update(VALUE_TYPE*);							// set Upper() e Lower() 
static std::list<ENTRY_TYPE>::iterator FindBase(const VALUE_TYPE * pt);	
																							// find the container that contains pt (naive)
virtual  void  AddDataElem(STL_CONT::value_type * pt,int n);// add n element to the auxiliary data

public:
static int & Id(){															// unique identifier of the istance
		static int id=0;														// (used as type trait)
		return id;
		}
};

// --------------------------- CATEntry: implementation --------------------

template <class STL_CONT, class ENTRY_TYPE>
unsigned int CATEntry<STL_CONT,ENTRY_TYPE>::

Ord(VALUE_TYPE * pt)
{
	Update(pt);
	return (pt-Lower());
}


template <class STL_CONT, class ENTRY_TYPE>
std::list<ENTRY_TYPE>::iterator CATEntry<STL_CONT,ENTRY_TYPE>::

FindBase(const VALUE_TYPE * pt)
{
std::list<ENTRY_TYPE>::iterator ite,curr_base,_;
ite = AT().begin();
curr_base = AT().end();

for(;ite != AT().end();ite++)
	if( pt < (*ite).Start())
		return curr_base;
	else
		curr_base = ite;

return curr_base;
}


template <class STL_CONT, class ENTRY_TYPE>
 bool CATEntry< STL_CONT, ENTRY_TYPE>::
 
IsTheSameAsLast(VALUE_TYPE * pt)
{
return ( UTD() && ( !(Lower()> pt)) && (pt < Upper()) );
}

template <class STL_CONT, class ENTRY_TYPE>
void CATEntry< STL_CONT, ENTRY_TYPE>::

Update(VALUE_TYPE * pt)
{
if(!IsTheSameAsLast(pt)){
	std::list<ENTRY_TYPE>::iterator lower_ite;
	lower_ite = FindBase(pt);

	assert(	lower_ite!=AT().end());

	Lower() = (*lower_ite).Start();
	if( (*lower_ite).Start() == AT().back().Start())
		Upper() = (VALUE_TYPE *) 0xffffffff;
	else
	{
		lower_ite++;		
		Upper() = (*lower_ite).Start();
	}
	
	Curr() = lower_ite;
	UTD() = true;
	}
}

template <class STL_CONT, class ENTRY_TYPE>
void CATEntry< STL_CONT,  ENTRY_TYPE>::
Resort(VALUE_TYPE* old_start,VALUE_TYPE* new_start)
{
AT().sort();
UTD() = false;
}

template <class STL_CONT, class ENTRY_TYPE>
void CATEntry<STL_CONT, ENTRY_TYPE>::

Remove( const STL_CONT & c )
{
std::list<ENTRY_TYPE>::iterator ite;
for(ite = AT().begin(); ite != AT().end();  ++ite)
	if((*ite).C() == &c)
		{
			AT().erase(ite);
			break;
		}
UTD() = false;
}

template <class STL_CONT, class ENTRY_TYPE>
void CATEntry<STL_CONT, ENTRY_TYPE>::

RemoveIfEmpty( const STL_CONT & c )
{
std::list<ENTRY_TYPE>::iterator ite;
for(ite = AT().begin(); ite != AT().end();  ++ite)
	if((*ite).C() == &c)
			if(!(*ite).Empty())
				AT().erase(ite);
UTD() = false;
}

template <class STL_CONT, class ENTRY_TYPE>
void CATEntry<STL_CONT, ENTRY_TYPE>::

Remove(VALUE_TYPE  *	pt)
{
	std::list<ENTRY_TYPE>::iterator lower_ite;
	lower_ite = FindBase(pt);
	AT().erase(lower_ite);
	UTD() = false;
		
}

template <class STL_CONT, class ENTRY_TYPE>
void CATEntry<STL_CONT, ENTRY_TYPE>::

Insert( STL_CONT & c,bool cond )
{
ENTRY_TYPE entry(c);
std::list<ENTRY_TYPE>::iterator lower_ite,upper_ite;
upper_ite = FindBase(&*c.begin());
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

template <class STL_CONT, class ENTRY_TYPE>
ENTRY_TYPE & CATEntry<STL_CONT, ENTRY_TYPE>::
GetEntry(STL_CONT::value_type*pt){
Update(pt);
return *Curr();
}

template <class STL_CONT, class ENTRY_TYPE>
void CATEntry<STL_CONT, ENTRY_TYPE>::

AddDataElem(STL_CONT::value_type * pt,int n)
{
Update(pt);
Curr()->Push_back(n);
}

//--------------------------------------------------------------------------------------------
// CAT: derivation of CATEntry for the case where the temporary data is unique for each type.
// VERY IMPORTANT: there cannot be two vector of value with the same type of temporary datya
// This class is used to implement optional core data (NormalOpt, CoordOpt etc...)
template <class STL_CONT,class ATTR_TYPE>
class CAT:public CATEntry<STL_CONT, EntryCAT<STL_CONT,ATTR_TYPE> >{
public:
static ATTR_TYPE & Get(STL_CONT::value_type * pt);
}; 
//---------------------- CAT: implementation---------------------------------------------------
template <class STL_CONT, class ATTR_TYPE>
ATTR_TYPE & CAT<STL_CONT,ATTR_TYPE>::

Get(STL_CONT::value_type * pt)
{
int ord = Ord(pt);
return Curr()->Data()[ord];
}
//---------------------------------------------------------------------------------------------

};//end namespace vcg

#endif