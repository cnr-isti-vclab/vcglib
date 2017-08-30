/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2017                                           \/)\/    *
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

#ifndef VCG__FOREACH_H
#define VCG__FOREACH_H

#ifndef __VCG_MESH
#error "This file should not be included alone. It is automatically included by complex.h"
#endif

namespace vcg {
namespace tri {
/** \addtogroup trimesh
@{
*/

template <class MeshType>
inline void ForEachFacePos(MeshType &m, std::function<void (typename face::Pos<typename MeshType::FaceType>  &)> action)
{
  typedef typename face::Pos<typename MeshType::FaceType> PosType;
  
    for(auto fi=m.face.begin();fi!=m.face.end();++fi)
        if(!(*fi).IsD())
        {
            for(int i=0;i<3;++i)
            {
                PosType pi(&*fi,i);
                action(pi);
            }
        }
}

template <class MeshType>
inline void ForEachFace(const MeshType &m, std::function<void (const typename MeshType::FaceType &)> action)
{
    for(auto fi=m.face.begin();fi!=m.face.end();++fi)
        if(!(*fi).IsD())
        {
          action(*fi);
        }
}

template <class MeshType>
inline void ForEachFace(MeshType &m, std::function<void (typename MeshType::FaceType &)> action)
{
    for(auto fi=m.face.begin();fi!=m.face.end();++fi)
        if(!(*fi).IsD())
        {
          action(*fi);
        }
}


template <class MeshType>
inline void forEachVertex(MeshType &m, std::function<void (typename MeshType::FaceType &)> action)
{
    for(auto fi=m.face.begin();fi!=m.face.end();++fi)
        if(!(*fi).IsD())
        {
          action(*fi);
        }
}



/** @} */ // end doxygen group trimesh
} // end namespace tri
} // end namespace vcg

#endif // VCG__FOREACH_H
