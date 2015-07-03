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

#ifndef __MEMORY_INFO_H
#define __MEMORY_INFO_H

#include <stdexcept>

namespace vcg 
{
	//WARNING: All the classes derived from MemoryInfo has been intended to be instantiated as a singleton in the host application 
	//(i.e. in every application using it just an instance of a class derived from MemoryInfo should be declared).

	class MemoryInfo
	{
	public:
		class MemoryInfoException : public std::exception
		{
		public:
			MemoryInfoException(const char* text)
				:std::exception(),_exctext(text){}

			~MemoryInfoException() throw() {}
			inline const char* what() const throw() {return _exctext;}
		private:
			const char* _exctext;
		};

		MemoryInfo(long long unsigned int originalmem)
			:_originaltotalmemory(originalmem),_currentfreememory(_originaltotalmemory)
		{       
		}

		virtual ~MemoryInfo() {}
		virtual void acquiredMemory(long long unsigned int mem) = 0;
		virtual long long unsigned int usedMemory() const = 0;
		virtual long long unsigned int currentFreeMemory() const = 0;
		virtual void releasedMemory(long long unsigned int mem = 0) = 0;
		virtual bool isAdditionalMemoryAvailable(long long unsigned int mem) = 0;

	protected:
		const long long unsigned int _originaltotalmemory;
		long long unsigned int _currentfreememory;
	};

	//WARNING: this is not a thread safe class. The object derived from MemoryInfo are intended to be used inside GLMeshAttributeFeeder as static variable in order to manage the available GPUMemory.
	//We strongly recommend you to define in your code a thread safe version of the class, defining mutexed access member functions. 
	//This class should be consider just as a basic example for the implementations of the required functionalities. 
	//It is safe to use it just when the user has only one mesh to pass to the GPU.

	class NotThreadSafeMemoryInfo : public MemoryInfo
	{
	public:
		NotThreadSafeMemoryInfo(long long unsigned int originalmem)
			:MemoryInfo(originalmem)
		{
		}

		~NotThreadSafeMemoryInfo() {}

		void acquiredMemory(long long unsigned int mem)
		{
			if (mem > _originaltotalmemory)
				throw MemoryInfo::MemoryInfoException("It has been requested more memory than the total one.\\n");
			else 
				if (mem > _currentfreememory)
					throw MemoryInfo::MemoryInfoException("It has been requested more memory than the free available one.\\n");
				else
					_currentfreememory -= mem;
		}

		long long unsigned int usedMemory() const
		{
			return _originaltotalmemory - _currentfreememory;
		}

		long long unsigned int currentFreeMemory() const
		{
			return _currentfreememory;
		}

		void releasedMemory(long long unsigned int mem = 0)
		{
			if (mem > _originaltotalmemory)
				throw MemoryInfo::MemoryInfoException("It has been released more memory than the total one. Something strange happened!\\n");
			else
				_currentfreememory += mem;
		}

		bool isAdditionalMemoryAvailable(long long unsigned int mem)
		{
			return (_currentfreememory >= mem);
		}
	};
};

#endif