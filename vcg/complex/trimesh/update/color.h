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
typedef UpdateMeshType MeshType; 
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;

static void VertexBorderFlag(MeshType &m, Color4b vb=Color4b::Blue)
{
	MeshType::VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if(!(*vi).IsD())
			(*vi).C()=Color4b::White;

	MeshType::FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
		for(int j=0;j<3;++j)
			if((*fi).IsB(j)){
				(*fi).V(j)->C() = vb;
				(*fi).V1(j)->C() = vb;
				(*fi).C() = vb;
			}
}

static void FaceBF(MeshType &m, Color4b vn=Color4b::White, Color4b vb=Color4b::Blue, 
						Color4b vc=Color4b::Red, Color4b vs=Color4b::LightBlue)
{
	MeshType::FaceIterator fi;
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
          if(face::IsManifold(*fi,j)){
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


static int FaceSelected(MeshType &m, Color4b vs=Color4b::LightBlue)
{
	int cnt=0;
	MeshType::FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
		if(!(*fi).IsD()) 
			if((*fi).IsS()) { (*fi).C() = vs; ++cnt; }
			else (*fi).C() = Color4b::White;
	return cnt;
}

static void FaceColorStrip(MeshType &m, std::vector<FacePointer> &TStripF)
{
	Color4b::Color4b cc[7]={Color4b::White ,Color4b::Red    ,Color4b::Green  ,Color4b::Blue   ,Color4b::Cyan   ,Color4b::Yellow ,Color4b::Magenta};
	int cnt=0;

	vector<FacePointer>::iterator fi;
	for(fi=TStripF.begin();fi!=TStripF.end();++fi)
		if(*fi) (**fi).C().ColorRamp(0,16,cnt);
		else cnt=(cnt+1)%16;
	//	if(*fi) (**fi).C()=cc[cnt];
  //			else cnt=(cnt+1)%7;
	
}


static int VertexSelected(MeshType &m, Color4b vs=Color4b::LightBlue)
{
	int cnt=0;
	MeshType::VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if(!(*vi).IsD()) 
			if((*vi).IsS()) {(*vi).C() = vs;  ++cnt; }
			else (*vi).C() = Color4b::White;
	
	return cnt;
}

static void VertexBorderManifoldFlag(MeshType &m, Color4b vn=Color4b::White, Color4b vb=Color4b::Blue, Color4b vc=Color4b::Red)
{
	MeshType::VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) if(!(*vi).IsD()) 
		(*vi).C()=vn;

	MeshType::FaceIterator fi;
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



static void FaceQuality(MeshType &m)
{
	// step 1: find the range
	MeshType::FaceIterator fi;
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

static void FaceQuality(MeshType &m, float minq, float maxq)
{
	MeshType::FaceIterator fi;

	for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD()) 		
		(*fi).C().ColorRamp(minq,maxq,(*fi).Q());
}

static void VertexQuality(MeshType &m, float minq, float maxq)
{
	MeshType::VertexIterator vi;

	for(vi=m.vert.begin();vi!=m.vert.end();++vi)		
		if(!(*vi).IsD()) 
			(*vi).C().ColorRamp(minq,maxq,(*vi).Q());
}

static void VertexQuality(MeshType &m)
{
	// step 1: find the range
	MeshType::VertexIterator vi;
  float minq=std::numeric_limits<float>::max(),
				maxq=std::numeric_limits<float>::min();
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)		
		if(!(*vi).IsD()) 
		{
			minq=min(minq,(*vi).Q());
			maxq=max(maxq,(*vi).Q());
		}
	VertexQuality(m,minq,maxq);
}

static void VertexQualityHistEq(MeshType &m)
{
	// step 1: find the range
	MeshType::VertexIterator vi;
	float minq=MaxVal(0.0f),
				maxq=-MaxVal(0.0f);
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)		
		if(!(*vi).IsD()) 
		{
			minq=min(minq,(*vi).Q());
			maxq=max(maxq,(*vi).Q());
		}
	// step 2; Get the distribution
	Hist H;
	H.SetRange(minq,maxq,1024);
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)		
		if(!(*vi).IsD()) H.Add((*vi).Q());
	
	VertexQuality(m,H.Percentile(.05f),H.Percentile(.95f));
}

};
/*@}*/
}// end namespace
}// end namespace
#endif