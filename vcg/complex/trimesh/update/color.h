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
Revision 1.6  2004/08/25 15:15:26  ganovelli
minor changes to comply gcc compiler (typename's and stuff)

Revision 1.5  2004/07/15 00:13:39  cignoni
Better doxigen documentation

Revision 1.4  2004/06/24 07:56:54  cignoni
now use std::numeric_limits instead of old max val()

Revision 1.3  2004/03/12 15:22:19  cignoni
Written some documentation and added to the trimes doxygen module

Revision 1.2  2004/03/10 00:48:06  cignoni
changed to the face::IsBorder() style

Revision 1.1  2004/03/05 10:59:24  cignoni
Changed name from plural to singular (normals->normal)


****************************************************************************/
#ifndef __VCG_TRI_UPDATE_COLOR
#define __VCG_TRI_UPDATE_COLOR
#include <limits>
namespace vcg {
namespace tri {
/** \addtogroup trimesh */
/*@{*/
/// Generation of per-vertex and per-face colors according to various strategy.
/// This class is used to compute per face or per vertex color with respect to for example Border (UpdateColor::VertexBorderFlag), Selection (UpdateColor::FaceSelected), Quality .

template <class UpdateMeshType>
class UpdateColor
{
public:
typedef  UpdateMeshType MeshType; 
typedef typename UpdateMeshType::VertexType     VertexType;
typedef typename UpdateMeshType::VertexPointer  VertexPointer;
typedef typename UpdateMeshType::VertexIterator VertexIterator;
typedef typename UpdateMeshType::FaceType       FaceType;
typedef typename UpdateMeshType::FacePointer    FacePointer;
typedef typename UpdateMeshType::FaceIterator   FaceIterator;

/** Color the vertex of the mesh that are on the border
It uses the information in the flags, and not any topology. So it just require that you have correctly computed (or loaded) the flags; See the 
**/
static void VertexBorderFlag( UpdateMeshType &m, Color4b vb=Color4b::Blue)
{
	typename UpdateMeshType::VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if(!(*vi).IsD())
			(*vi).C()=Color4b::White;

	typename UpdateMeshType::FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
		for(int j=0;j<3;++j)
			if((*fi).IsB(j)){
				(*fi).V(j)->C() = vb;
				(*fi).V1(j)->C() = vb;
				(*fi).C() = vb;
			}
}

static void FaceBF( UpdateMeshType &m, Color4b vn=Color4b::White, Color4b vb=Color4b::Blue, 
						Color4b vc=Color4b::Red, Color4b vs=Color4b::LightBlue)
{
	typename UpdateMeshType::FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
		if(!(*fi).IsD())
			(*fi).C() = vn;

	for(fi=m.face.begin();fi!=m.face.end();++fi)		
		if(!(*fi).IsD())
		{
			if((*fi).IsS())
				(*fi).C() = vs;
			else
			{
				for(int j=0;j<3;++j)
					if(*fi.IsManifold(j)){
						if((*fi).IsB(j)){
							(*fi).C() = vb;
							(*fi).C() = vb;
						}	
					}
					else
					{
						(*fi).C() = vc;
						(*fi).C() = vc;
					}
			}
		}
}


static int FaceSelected(UpdateMeshType &m, Color4b vs=Color4b::LightBlue)
{
	int cnt=0;
	typename UpdateMeshType::FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
		if(!(*fi).IsD()) 
			if((*fi).IsS()) { (*fi).C() = vs; ++cnt; }
			else (*fi).C() = Color4b::White;
	return cnt;
}

static void FaceColorStrip(UpdateMeshType &m, std::vector<FacePointer> &TStripF)
{
	 vcg::Color4b cc[7]={
								  vcg::Color4b::White ,
									vcg::Color4b::Red    ,
									vcg::Color4b::Green  ,
									vcg::Color4b::Blue   ,
									vcg::Color4b::Cyan   ,
									vcg::Color4b::Yellow ,
									vcg::Color4b::Magenta
	};
	int cnt=0;

	typename std::vector<FacePointer>::iterator fi;
	for(fi=TStripF.begin();fi!=TStripF.end();++fi)
		if(*fi) (**fi).C().ColorRamp(0,16,cnt);
		else cnt=(cnt+1)%16;
	//	if(*fi) (**fi).C()=cc[cnt];
  //			else cnt=(cnt+1)%7;
	
}


static int VertexSelected(UpdateMeshType &m, Color4b vs=Color4b::LightBlue)
{
	int cnt=0;
	typename UpdateMeshType::VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if(!(*vi).IsD()) 
			if((*vi).IsS()) {(*vi).C() = vs;  ++cnt; }
			else (*vi).C() = Color4b::White;
	
	return cnt;
}

static void VertexConstant(UpdateMeshType &m, Color4b c=Color4b::White)
{
	typename UpdateMeshType::VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) if(!(*vi).IsD()) 
		(*vi).C()=c;
}

static void FaceConstant(UpdateMeshType &m, Color4b c=Color4b::White)
{
	typename UpdateMeshType::FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)		
		(*fi).C()=c;
}

static void VertexBorderManifoldFlag(UpdateMeshType &m, Color4b vn=Color4b::White, Color4b vb=Color4b::Blue, Color4b vc=Color4b::Red)
{
	typename UpdateMeshType::VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) if(!(*vi).IsD()) 
		(*vi).C()=vn;

	typename UpdateMeshType::FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)		
		if(!(*fi).IsD()) 
		for(int j=0;j<3;++j)
			if((*fi).IsManifold(j)){
				if((*fi).IsB(j)){
				(*fi).V(j)->C()=vb;
				(*fi).V1(j)->C()=vb;
				}	
			}
			else
			{
				(*fi).V(j)->C()=vc;
				(*fi).V1(j)->C()=vc;
			}
}



static void FaceQuality(UpdateMeshType &m)
{
	// step 1: find the range
	typename UpdateMeshType::FaceIterator fi;
	float minq=m.face[0].Q(),
				maxq=m.face[0].Q();
	for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
		if(!(*fi).IsD()) 
		{
			minq=min(minq,(*fi).Q());
			maxq=max(maxq,(*fi).Q());
		}

	FaceQuality(m,minq,maxq);
}

static void FaceQuality(UpdateMeshType &m, float minq, float maxq)
{
	typename UpdateMeshType::FaceIterator fi;

	for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD()) 		
		(*fi).C().ColorRamp(minq,maxq,(*fi).Q());
}

static void VertexQuality(UpdateMeshType &m, float minq, float maxq)
{
	typename UpdateMeshType::VertexIterator vi;

	for(vi=m.vert.begin();vi!=m.vert.end();++vi)		
		if(!(*vi).IsD()) 
			(*vi).C().ColorRamp(minq,maxq,(*vi).Q());
}

static void VertexQuality(UpdateMeshType &m)
{
	// step 1: find the range
	typename UpdateMeshType::VertexIterator vi;
  float minq=std::numeric_limits<float>::max(),
				maxq=std::numeric_limits<float>::min();
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)		
		if(!(*vi).IsD()) 
		{	
			minq=vcg::math::Min(minq,(*vi).Q());
			maxq=vcg::math::Max(maxq,(*vi).Q());
		}
	VertexQuality(m,minq,maxq);
}

static void VertexQualityHistEq(UpdateMeshType &m)
{
	// step 1: find the range
	typename UpdateMeshType::VertexIterator vi;
	float minq=MaxVal(0.0f),
				maxq=-MaxVal(0.0f);
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)		
		if(!(*vi).IsD()) 
		{
			minq=min(minq,(*vi).Q());
			maxq=max(maxq,(*vi).Q());
		}
	// step 2; Get the distribution
	//	Hist H;
	//H.SetRange(minq,maxq,1024);
	//for(vi=m.vert.begin();vi!=m.vert.end();++vi)		
	//	if(!(*vi).IsD()) H.Add((*vi).Q());
	
	// VertexQuality(m,H.Percentile(.05f),H.Percentile(.95f));
}

};
/*@}*/
}// end namespace
}// end namespace
#endif
