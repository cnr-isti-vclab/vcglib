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
Revision 1.3  2004/03/31 22:36:44  ganovelli
First Working Release (with this comment)

/****************************************************************************/

#ifndef __VCGLIB_TRACED_VECTOR__
#define __VCGLIB_TRACED_VECTOR__ 


#include <vcg/container/container_allocation_table.h>
#include <vcg/container/entries_allocation_table.h>

#include <assert.h>

namespace vcg {

template <class VALUE_TYPE>
class TVector: public std::vector<VALUE_TYPE>{
	typedef typename TVector<VALUE_TYPE> ThisType;

public:
	TVector():std::vector<VALUE_TYPE>(){reserve(1);}
	~TVector();
	

	std::list < CATBase<ThisType>* > attributes;
	// override di tutte le funzioni che possono spostare 
	// l'allocazione in memoria del container
	void push_back(const VALUE_TYPE & v);
	void pop_back();
	void resize(const unsigned int & size);
	void reserve(const unsigned int & size);

	template <class ATTR_TYPE>
		void EnableAttribute(){
			CAT<ThisType,ATTR_TYPE> * cat = new CAT<ThisType,ATTR_TYPE>();
			cat->Insert(*this);
			attributes.push_back(cat);
			}

	template <class ATTR_TYPE>
		void DisableAttribute(){
				std::list < CATBase<ThisType> * >::iterator ia; 
				for(ia = attributes.begin(); ia != attributes.end(); ++ia)
					if((*ia)->Id() == CAT<ThisType,ATTR_TYPE>::Id())
						{
							(*ia)->Remove(*this);
							delete (*ia);
							attributes.erase(ia);
							break;
						}
				}

	template <class ATTR_TYPE>
		TempData<ThisType,ATTR_TYPE> NewTempData(){
			typedef typename CATEntry<ThisType,EntryCATMulti<ThisType> >::EntryType EntryTypeMulti;
			CATEntry<ThisType,EntryTypeMulti>::Insert(*this);
			EntryTypeMulti	entry = CATEntry<ThisType,EntryTypeMulti >::GetEntry(&*begin());
			entry.Data().push_back(new Wrap< ATTR_TYPE>);

			((Wrap<ATTR_TYPE>*)entry.Data().back())->reserve(capacity());
			((Wrap<ATTR_TYPE>*)entry.Data().back())->resize(size());

			return TempData<ThisType,ATTR_TYPE>((Wrap<ATTR_TYPE>*) entry.Data().back());
			}
			
	template <class ATTR_TYPE>
		void DeleteTempData(TempData<ThisType,ATTR_TYPE> & td){
			typedef typename CATEntry<ThisType,EntryCATMulti<ThisType> >::EntryType EntryTypeMulti;
			CATEntry<ThisType,EntryTypeMulti >::RemoveIfEmpty(*this);
			EntryTypeMulti
				entry = CATEntry<ThisType,EntryCATMulti<ThisType> >::GetEntry(&*begin());

			entry.Data().remove((Wrap<ATTR_TYPE>*)td.Item());
			delete ((Wrap<ATTR_TYPE>*)td.Item());
			}


private:	
	VALUE_TYPE * old_start;
	void Update();
};

template <class VALUE_TYPE>
void TVector<VALUE_TYPE>::push_back(const VALUE_TYPE & v){
	std::vector<VALUE_TYPE>::push_back(v);
	Update();	
	std::list < CATBase<ThisType> * >::iterator ia; 
	for(ia = attributes.begin(); ia != attributes.end(); ++ia)
		(*ia)->AddDataElem(&(*(this->begin())),1);

}
template <class VALUE_TYPE>
void TVector<VALUE_TYPE>::pop_back(){
	std::vector<VALUE_TYPE>::pop_back();
	Update();
}

template <class VALUE_TYPE>
void TVector<VALUE_TYPE>::resize(const unsigned int & size){
	std::vector<VALUE_TYPE>::resize(size);
	std::list < CATBase<ThisType> * >::iterator ia; 
	for(ia = attributes.begin(); ia != attributes.end(); ++ia)
		(*ia)->
	Update();
}

template <class VALUE_TYPE>
void TVector<VALUE_TYPE>::reserve(const unsigned int & size){
	std::vector<VALUE_TYPE>::reserve(size);
	Update();
}

template <class VALUE_TYPE>
	void TVector<VALUE_TYPE>::
		Update(){
		std::list < CATBase<ThisType> * >::iterator ia; 
		if(&(*begin()) != old_start)
			for(ia = attributes.begin(); ia != attributes.end(); ++ia)
				(*ia)->Resort(old_start,&(*begin()));

		old_start = &(*begin());
	}



template <class VALUE_TYPE>
TVector<VALUE_TYPE>::~TVector(){
		std::list < CATBase<ThisType> * >::iterator ia; 
		for(ia = attributes.begin(); ia != attributes.end(); ++ia)
			{	
				(*ia)->Remove(*this);
				delete *ia;
			}
		}

}; // end namespace
#endif
