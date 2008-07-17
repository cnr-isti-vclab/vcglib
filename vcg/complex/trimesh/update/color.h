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
Revision 1.11  2007/05/04 16:34:31  ganovelli
changes to comply "plus" types

Revision 1.10  2006/05/21 07:00:01  cignoni
Removed not working Equalized color (use funcs in stat.h)

Revision 1.9  2006/03/01 10:29:55  ponchio
HACK: MaxVal(0.0f) not defined in vcg/math/base.h as it should be,
changing it to 1e36 (pretty close :P)

Revision 1.8  2005/12/19 16:47:42  cignoni
Better comment and a parameter more for UpdateColor::VertexBorderFlag

Revision 1.7  2005/08/08 10:28:13  ganovelli
added math:: namespace before min and max

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
#include <math.h>
#include <vcg/space/color4.h>
#include <vcg/math/histogram.h>

namespace vcg {
namespace tri {

/// \ingroup trimesh

/// \headerfile color.h vcg/complex/trimesh/update/color.h

/// \brief Generation of per-vertex and per-face colors according to various strategy.
/**
This class is used to compute per face or per vertex color with respect to for example Border (UpdateColor::VertexBorderFlag), Selection (UpdateColor::FaceSelected), Quality .
*/

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

/// \brief Color the vertexes of the mesh that are on the border

/**
It uses the information in the Vertex flags, and not any topology.
So it just require that you have correctly computed the flags;

   - vcg::tri::UpdateTopology<Mesh>::FaceFace(m.cm);
   - vcg::tri::UpdateFlags<Mesh>::FaceBorderFromFF(m.cm);
   - vcg::tri::UpdateFlags<Mesh>::VertexBorderFromFace (m.cm);
   - vcg::tri::UpdateColor<Mesh>::VertexBorderFlag(m.cm);

*/
static void VertexBorderFlag( UpdateMeshType &m, Color4b BorderColor=Color4b::Blue, Color4b InternalColor=Color4b::White, Color4b MixColor=Color4b::Cyan)
{
	Color4b BaseColor = Color4b::Green;

	VertexConstant(m,BaseColor);

	typename UpdateMeshType::FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
		for(int j=0;j<3;++j)
			if((*fi).IsB(j)){
				if( (*fi).V(j)->C() == BaseColor)     (*fi).V(j)->C() = BorderColor;
				if( (*fi).V(j)->C() == InternalColor) (*fi).V(j)->C() = MixColor;
				if( (*fi).V1(j)->C() == BaseColor)     (*fi).V1(j)->C() = BorderColor;
				if( (*fi).V1(j)->C() == InternalColor) (*fi).V1(j)->C() = MixColor;
			} else
			{
				if( (*fi).V(j)->C() == BaseColor)     (*fi).V(j)->C() = InternalColor;
				if( (*fi).V(j)->C() == BorderColor) (*fi).V(j)->C() = MixColor;
				if( (*fi).V1(j)->C() == BaseColor)     (*fi).V1(j)->C() = InternalColor;
				if( (*fi).V1(j)->C() == BorderColor) (*fi).V1(j)->C() = MixColor;
			}

}

/// This function colores the face of a mesh randomly.
/// The feature bit is used to color polygonal faces uniformly

static void MultiFaceRandom( UpdateMeshType &m)
{
	FaceIterator fi;
	Color4b BaseColor = Color4b::Black;
	FaceConstant(m,BaseColor);
	int id=0;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
			if(!(*fi).IsD())
			{
				id++;
				if((*fi).C() == BaseColor) (*fi).C() = Color4b::Scatter(50, id%50,.4f,.7f);
				for(int j=0;j<3;++j)
					if((*fi).IsF(j))
					{
						assert(!IsBorder((*fi),j));
						(*fi).FFp(j)->C()= (*fi).C();
					}
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
				maxq=-std::numeric_limits<float>::max();
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if(!(*vi).IsD())
		{
			minq=vcg::math::Min<float>(minq,(*vi).Q());
			maxq=vcg::math::Max<float>(maxq,(*vi).Q());
		}
	VertexQuality(m,minq,maxq);
}

//Fill the mesh with the selected color.
static int Filling(UpdateMeshType &m, Color4b c, const bool ProcessSelected=false)
{
  int counter=0;
  VertexIterator vi;
  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        (*vi).C() = c;
        ++counter;
      }
    }
  }
  return counter;
}

//Reduces the mesh to two colors according to a treshold.
static int Tresholding(UpdateMeshType &m, float treshold, Color4b c1 = Color4<unsigned char>::Black, Color4b c2 = Color4<unsigned char>::White, const bool ProcessSelected=false)
{
  int counter=0;
  VertexIterator vi;
  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        float value = ComputeLightness((*vi).C());

        if(value<=treshold) (*vi).C() = c1;
        else (*vi).C() = c2;
        ++counter;
      }
    }
  }
  return counter;
}

//Computes the luminance value for a specified color.
static float ComputeLightness(Color4b c)
{
  float min_rgb = math::Min((float)c[0],(float)c[1]);
  min_rgb = math::Min(min_rgb,(float)c[2]);
  float max_rgb = math::Max((float)c[0],(float)c[1]);
  max_rgb = math::Max(max_rgb,(float)c[2]);
  return (max_rgb + min_rgb)/2;
}

//Apply the brightness filter, with the given amount, to the mesh.
static int Brighting(UpdateMeshType &m, float amount, const bool ProcessSelected=false)
{
	int counter=0;
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
	{
		if(!(*vi).IsD()) //if it has not been deleted...
		{
			if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        (*vi).C() = Color4b(
                            math::Clamp(int((*vi).C()[0]+amount),0,255),
                            math::Clamp(int((*vi).C()[1]+amount),0,255),
                            math::Clamp(int((*vi).C()[2]+amount),0,255),
                            255);
        ++counter;
      }
    }
  }
	return counter;
}

static int Contrast(UpdateMeshType &m, float factor, const bool ProcessSelected=false)
{
  int counter=0;
  VertexIterator vi;
  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        (*vi).C() = ColorMul((*vi).C(),factor);
        ++counter;
      }
    }
  }
  return counter;
}

//Subtracts a middle value, multiplies the rgb components of the color for a factor,
//and adds the middle value back.This is used for contrast operation.
static Color4b ColorMul(Color4b c, float factor)
{
  return Color4b( ValueMul(c[0], factor), ValueMul(c[1], factor), ValueMul(c[2], factor), 1);
}

static int ValueMul(int value, float factor)
{
  return math::Clamp<int>((int)((value - 128)*factor + 128), 0, 255);
}

static int ContrastBrightness(UpdateMeshType &m, float factor, float amount, const bool ProcessSelected=false)
{
  int counter=0;
  VertexIterator vi;
  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        (*vi).C() = ColorMulAdd((*vi).C(),factor,amount);
        ++counter;
      }
    }
  }
  return counter;
}

//This is a composition of ColorMul() and ColorAdd(), used for Contrast&Brightness operations.
//The result is clamped just one time after all computations; this get a more accurate result.
static Color4b ColorMulAdd(Color4b c, float factor, float amount)
{
  return Color4b( ValueMulAdd(c[0], factor, amount), ValueMulAdd(c[1], factor, amount), ValueMulAdd(c[2], factor, amount), 1 );
}

static int ValueMulAdd(int value, float factor, float amount)
{
  return math::Clamp<int>((int)((value - 128)*factor + 128 + amount), 0, 255);
}

//Invert the rgb components of the color.
static int Invert(UpdateMeshType &m, const bool ProcessSelected=false)
{
  int counter=0;
  VertexIterator vi;
  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        (*vi).C() = ColorInvert((*vi).C());
        ++counter;
      }
    }
  }
  return counter;
}

//invert the given color
static Color4b ColorInvert(Color4b c)
{
  return Color4b( ValueInvert(c[0]), ValueInvert(c[1]), ValueInvert(c[2]), 1);
}

static int ValueInvert(int value)
{
    return 255-value;
}
//Apply the gamma correction filter, with the given gamma exponet, to the mesh.
static int Gamma(UpdateMeshType &m, float gamma, const bool ProcessSelected=false)
{
  int counter=0;

  VertexIterator vi;
  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        (*vi).C() = ColorPow((*vi).C(), 1/gamma);
        ++counter;
      }
    }
  }
  return counter;
}

//computes the gamma transformation on a given color, according to new_val = old_val^gamma
static Color4b ColorPow(Color4b c, float exponent)
{
  return vcg::Color4b(
                      math::Clamp((int)( ValuePow(float(c[0])/255, exponent)*255), 0, 255),
                      math::Clamp((int)( ValuePow(float(c[1])/255, exponent)*255), 0, 255),
                      math::Clamp((int)( ValuePow(float(c[2])/255, exponent)*255), 0, 255),
                      255);
}

static float ValuePow(float value, float exponent)
{
  return powf(value, exponent);
}

static int Levels(UpdateMeshType &m, float gamma, float in_min, float in_max, float out_min, float out_max, unsigned char rgbMask, const bool ProcessSelected=false)
{
  int counter=0;
  VertexIterator vi;
  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        (*vi).C() = ColorLevels((*vi).C(), gamma, in_min, in_max, out_min, out_max, rgbMask);
        ++counter;
      }
    }
  }
  return counter;
}

enum rgbChMask {ALL_CHANNELS = 7, RED_CHANNEL = 4, GREEN_CHANNEL = 2, BLUE_CHANNEL = 1, NO_CHANNELS = 0 };

static Color4b ColorLevels(Color4b c, float gamma, float in_min, float in_max, float out_min, float out_max, unsigned char rgbMask)
{
  unsigned char r = c[0], g = c[1], b = c[2];
  if(rgbMask & RED_CHANNEL) r = ValueLevels(c[0], gamma, in_min, in_max, out_min, out_max);
  if(rgbMask & GREEN_CHANNEL) g = ValueLevels(c[1], gamma, in_min, in_max, out_min, out_max);
  if(rgbMask & BLUE_CHANNEL) b = ValueLevels(c[2], gamma, in_min, in_max, out_min, out_max);
  return Color4b(r, g, b, 255);
}

static int ValueLevels(int value, float gamma, float in_min, float in_max, float out_min, float out_max)
{
  float fvalue = value/255.0f;
  // normalize
  fvalue = math::Clamp<float>(fvalue - in_min, 0.0f, 1.0f) / math::Clamp<float>(in_max - in_min, 1.0f/255.0f, 1.0f);
  // transform gamma
  fvalue = pow(fvalue,1/gamma);
  // rescale range
  fvalue = fvalue * (out_max - out_min) + out_min;
  //back in interval [0,255] and clamp
  return math::Clamp<int>((int)(fvalue * 255), 0, 255);
}

static int Colourisation(UpdateMeshType &m, Color4b c, float intensity, const bool ProcessSelected=false)
{
  int counter=0;
  VertexIterator vi;
  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        (*vi).C() = ColorApplyDiff((*vi).C(), c, intensity);
        ++counter;
      }
    }
  }
  return counter;
}

static Color4b ColorApplyDiff(Color4b old_color, Color4b new_color, float intensity)
{
  return Color4b( ValueApplyDiff(old_color[0],new_color[0],intensity), ValueApplyDiff(old_color[1],new_color[1],intensity), ValueApplyDiff(old_color[2], new_color[2],intensity), 255);
}

static int ValueApplyDiff(int old_value, int new_value, float intensity)
{
  return  math::Clamp<int>((int)(old_value + intensity * (new_value - old_value)), 0, 255);
}

static int Desaturation(UpdateMeshType &m, int method, const bool ProcessSelected=false)
{
  int counter=0;
  VertexIterator vi;
  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        (*vi).C() = ColorDesaturate((*vi).C(), method);
        ++counter;
      }
    }
  }
  return counter;
}

enum DesaturationMethods {M_LIGHTNESS = 0, M_MEAN = 1};

static Color4b ColorDesaturate(Color4b c, int method)
{
  switch(method){
    case M_LIGHTNESS:{
      int val = (int)ComputeLightness(c);
      return Color4b( val, val, val, 255);
    }
    case M_MEAN:{
      int val = (int)ComputeMeanLightness(c);
      return Color4b( val, val, val, 255);
    }
    default: assert(0);
  }
}

static float ComputeMeanLightness(Color4b c)
{
    return float(c[0]+c[1]+c[2])/3.0f;
}

static int Equalize(UpdateMeshType &m, unsigned int rgbMask, const bool ProcessSelected=false)
{
  Histogramf Hl, Hr, Hg, Hb;
  Hl.Clear(); Hr.Clear(); Hg.Clear(); Hb.Clear();
	Hl.SetRange(0, 255, 255); Hr.SetRange(0, 255, 255); Hg.SetRange(0, 255, 255); Hb.SetRange(0, 255, 255);

  int counter=0;
  VertexIterator vi;

  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        int v = (int)(ComputeLightness((*vi).C())+0.5);
        Hl.Add(v); Hr.Add((*vi).C()[0]); Hg.Add((*vi).C()[1]); Hb.Add((*vi).C()[2]);
      }
    }
  }

	int cdf_l[256], cdf_r[256], cdf_g[256], cdf_b[256];
	cdf_l[0] = Hl.BinCount(0); cdf_r[0] = Hr.BinCount(0); cdf_g[0] = Hg.BinCount(0); cdf_b[0] = Hb.BinCount(0);
	for(int i=1; i<256; i++){
    cdf_l[i] = Hl.BinCount(float(i)) + cdf_l[i-1];
    cdf_r[i] = Hr.BinCount(float(i)) + cdf_r[i-1];
    cdf_g[i] = Hg.BinCount(float(i)) + cdf_g[i-1];
    cdf_b[i] = Hb.BinCount(float(i)) + cdf_b[i-1];
  }

  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        (*vi).C()=ColorEqualize((*vi).C(), cdf_l, cdf_r, cdf_g, cdf_b, rgbMask);
        ++counter;
      }
    }
  }
  return counter;
}

static Color4b ColorEqualize(Color4b c, int cdf_l[256], int cdf_r[256], int cdf_g[256], int cdf_b[256], unsigned int rgbMask)
{
  unsigned char r = c[0], g = c[1], b = c[2];
  if(rgbMask == NO_CHANNELS){
    int v = ValueEqualize(cdf_l[(int)(ComputeLightness(c)+0.5f)], cdf_l[0], cdf_l[255]);
    return Color4b(v, v, v, 255);
  }
  if(rgbMask & RED_CHANNEL) r = ValueEqualize(cdf_r[c[0]], cdf_r[0], cdf_r[255]);
  if(rgbMask & GREEN_CHANNEL) g = ValueEqualize(cdf_g[c[1]], cdf_g[0], cdf_g[255]);
  if(rgbMask & BLUE_CHANNEL) b = ValueEqualize(cdf_b[c[2]], cdf_b[0], cdf_b[255]);
  return Color4b(r, g, b, 255);
}

static int ValueEqualize(int cdfValue, int cdfMin, int cdfMax)
{
  return int(float((cdfValue - cdfMin)/float(cdfMax - cdfMin)) * 255.0f);
}

};

}// end namespace
}// end namespace
#endif
