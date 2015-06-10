/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
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

#include "ml_thread_safe_memory_info.h"

MLThreadSafeMemoryInfo::MLThreadSafeMemoryInfo( long long unsigned int originalmem ) 
	:vcg::NotThreadSafeMemoryInfo(originalmem),lock(QReadWriteLock::Recursive)
{

}

MLThreadSafeMemoryInfo::~MLThreadSafeMemoryInfo()
{
}

void MLThreadSafeMemoryInfo::acquiredMemory(long long unsigned int mem)
{
	lock.lockForWrite();
	vcg::NotThreadSafeMemoryInfo::acquiredMemory(mem);
	lock.unlock();
}

long long unsigned int MLThreadSafeMemoryInfo::usedMemory() const
{
	lock.lockForRead();
	long long unsigned int tmp = vcg::NotThreadSafeMemoryInfo::usedMemory();
	lock.unlock();
	return tmp;
}

long long unsigned int MLThreadSafeMemoryInfo::currentFreeMemory() const
{
	lock.lockForRead();
	long long unsigned int tmp = vcg::NotThreadSafeMemoryInfo::currentFreeMemory();
	lock.unlock();
	return tmp;
}

void MLThreadSafeMemoryInfo::releasedMemory(long long unsigned int mem)
{
	lock.lockForWrite();
	vcg::NotThreadSafeMemoryInfo::releasedMemory(mem);
	lock.unlock();
}

bool MLThreadSafeMemoryInfo::isAdditionalMemoryAvailable( long long unsigned int mem )
{
	lock.lockForRead();
	bool tmp = vcg::NotThreadSafeMemoryInfo::isAdditionalMemoryAvailable(mem);
	lock.unlock();
	return tmp;
}
