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


#ifndef __VCGLIB_TRACED_VECTOR__
#define __VCGLIB_TRACED_VECTOR__


#include <vcg/container_allocation_table.h>
#include <assert.h>

namespace vcg {

template <class VALUE_TYPE>
class TVector: public std::vector<VALUE_TYPE>{
	typedef TVector<VALUE_TYPE> THIS_TYPE;
public:
	TVector():std::vector<VALUE_TYPE>(){reserve(1);}
	~TVector();
	

	std::list < CATBase<THIS_TYPE>* > attributes;
	// override di tutte le funzioni che possono spostare 
	// l'allocazione in memoria del container
	void push_back(const VALUE_TYPE & v);
	void pop_back();
	void resize(const unsigned int & size);
	void reserve(const unsigned int & size);

	template <class ATTR_TYPE>
		void EnableAttribute(){
			CAT<THIS_TYPE,ATTR_TYPE> * cat = new CAT<THIS_TYPE,ATTR_TYPE>();
			cat->Insert(*this);
			attributes.push_back(cat);
			}

	template <class ATTR_TYPE>
		void DisableAttribute(){
				std::list < CATBase<THIS_TYPE> * >::iterator ia; 
				for(ia = attributes.begin(); ia != attributes.end(); ++ia)
					if((*ia)->Id() == CAT<THIS_TYPE,ATTR_TYPE>::Id())
						{
							(*ia)->Remove(*this);
							delete (*ia);
							attributes.erase(ia);
							break;
						}
				}

	template <class ATTR_TYPE>
		TempData<THIS_TYPE,ATTR_TYPE> NewTempData(){
			//CAT<THIS_TYPE,EntryCATMulti>::Insert(*this)
			CATMulti<THIS_TYPE,EntryCATMulti<THIS_TYPE> >::EntryType 
				entry = CATMulti<THIS_TYPE,EntryCATMulti<THIS_TYPE> >::GetEntry(&*begin());
			entry.Data().push_back(new std::vector<ATTR_TYPE>);

			((std::vector<ATTR_TYPE>*)entry.Data().back())->reserve(capacity());
			((std::vector<ATTR_TYPE>*)entry.Data().back())->resize(size());

			return TempData<THIS_TYPE,ATTR_TYPE>((std::vector<ATTR_TYPE>*) entry.Data().back());
			}
			
	template <class ATTR_TYPE>
		void DeleteTempData(TempData<THIS_TYPE,ATTR_TYPE> & td){
			//CAT<THIS_TYPE,EntryCATMulti>::Insert(*this)
			CATMulti<THIS_TYPE,EntryCATMulti<THIS_TYPE> >::EntryType 
				entry = CATMulti<THIS_TYPE,EntryCATMulti<THIS_TYPE> >::GetEntry(&*begin());

			entry.Data().remove(td.Item());
			delete td.Item();
			}


private:	
	VALUE_TYPE * old_start;
	void Update();
};

template <class VALUE_TYPE>
void TVector<VALUE_TYPE>::push_back(const VALUE_TYPE & v){
	std::vector<VALUE_TYPE>::push_back(v);
	std::list < CATBase<THIS_TYPE> * >::iterator ia; 
	for(ia = attributes.begin(); ia != attributes.end(); ++ia)
		(*ia)->AddDataElem(&(*(this->begin())),1);
	Update();
}
template <class VALUE_TYPE>
void TVector<VALUE_TYPE>::pop_back(){
	std::vector<VALUE_TYPE>::pop_back();
	Update();
}

template <class VALUE_TYPE>
void TVector<VALUE_TYPE>::resize(const unsigned int & size){
	std::vector<VALUE_TYPE>::resize(size);
	std::list < CATBase<THIS_TYPE> * >::iterator ia; 
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
		std::list < CATBase<THIS_TYPE> * >::iterator ia; 
		if(&(*begin()) != old_start)
			for(ia = attributes.begin(); ia != attributes.end(); ++ia)
				(*ia)->Resort(old_start,&(*begin()));

		old_start = &(*begin());
	}



template <class VALUE_TYPE>
TVector<VALUE_TYPE>::~TVector(){
		std::list < CATBase<THIS_TYPE> * >::iterator ia; 
		for(ia = attributes.begin(); ia != attributes.end(); ++ia)
			{	
				(*ia)->Remove(*this);
				delete *ia;
			}
		}

}; // end namespace
#endif
