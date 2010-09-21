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
#include <vcg/complex/trimesh/stat.h>
#include <vcg/math/perlin_noise.h>
#include <vcg/math/random_generator.h>

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

	
	
	class ColorAvgInfo 
	{
	public:
		unsigned int r;
		unsigned int g;
		unsigned int b;
		unsigned int a;
		int cnt;
	};
	
	static void VertexFromFace( UpdateMeshType &m)
	{
		ColorAvgInfo csi;
		csi.r=0; csi.g=0; csi.b=0; csi.cnt=0;
		SimpleTempData<typename UpdateMeshType::VertContainer, ColorAvgInfo> TD(m.vert,csi);
		
		FaceIterator fi;
			for(fi=m.face.begin();fi!=m.face.end();++fi)
				if(!(*fi).IsD()) 
					for(int j=0;j<3;++j)
						{
							TD[(*fi).V(j)].r+=(*fi).C()[0];
							TD[(*fi).V(j)].g+=(*fi).C()[1];
							TD[(*fi).V(j)].b+=(*fi).C()[2];
							TD[(*fi).V(j)].a+=(*fi).C()[3];
							++TD[(*fi).V(j)].cnt;
						}
		
		VertexIterator vi;
		for(vi=m.vert.begin();vi!=m.vert.end();++vi)
				if(!(*vi).IsD() && TD[*vi].cnt>0 )
					{
						(*vi).C()[0] = TD[*vi].r / TD[*vi].cnt;
						(*vi).C()[1] = TD[*vi].g / TD[*vi].cnt;
						(*vi).C()[2] = TD[*vi].b / TD[*vi].cnt;
						(*vi).C()[3] = TD[*vi].a / TD[*vi].cnt;
					}
		}

	static void FaceFromVertex( UpdateMeshType &m)
	{
		FaceIterator fi;
		for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD()) 
		{
			Color4f avg = (Color4f::Construct((*fi).V(0)->C()) + 
										 Color4f::Construct((*fi).V(1)->C()) + 
										 Color4f::Construct((*fi).V(2)->C()) )/ 3.0;
			(*fi).C().Import(avg);
		}
	}
	
	
	
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
/// The faux bit is used to color polygonal faces uniformly
static void MultiFaceRandom( UpdateMeshType &m)
{
	FaceIterator fi;
	Color4b BaseColor = Color4b::Black;
	FaceConstant(m,BaseColor);
    int id_num=0;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
			if(!(*fi).IsD())
			{
                id_num++;
                if((*fi).C() == BaseColor) (*fi).C() = Color4b::Scatter(50, id_num%50,.4f,.7f);
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
            if(!(*fi).IsD()){
			if((*fi).IsS()) { (*fi).C() = vs; ++cnt; }
			else (*fi).C() = Color4b::White;
                    }
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
            if(!(*vi).IsD()){
			if((*vi).IsS()) {(*vi).C() = vs;  ++cnt; }
			else (*vi).C() = Color4b::White;
                    }
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

static void FaceQualityGray(UpdateMeshType &m, float minq, float maxq)
{
	typename UpdateMeshType::FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
		(*fi).C().SetGrayShade( ((*fi).Q()-minq)/(maxq-minq));
}

static void FaceQualityGray(UpdateMeshType &m)
{
	std::pair<float,float> minmax = Stat<UpdateMeshType>::ComputePerFaceQualityMinMax(m);
	FaceQualityGray(m,minmax.first,minmax.second);
}	

static void FaceQualityRamp(UpdateMeshType &m,bool selected=false)
{
  // step 1: find the range
	typename UpdateMeshType::FaceIterator fi;
  float minq=std::numeric_limits<float>::max(),
        maxq=-std::numeric_limits<float>::max();

  for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
  if(!selected || (*fi).IsS())
    {
    minq=std::min(minq,(*fi).Q());
    maxq=std::max(maxq,(*fi).Q());
		}

  FaceQualityRamp(m,minq,maxq,selected);
}

static void FaceQualityRamp(UpdateMeshType &m, float minq, float maxq,bool selected=false)
{
	typename UpdateMeshType::FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
    if(!selected || (*fi).IsS())
      (*fi).C().ColorRamp(minq,maxq,(*fi).Q());
}

static void VertexQualityRamp(UpdateMeshType &m, float minq, float maxq)
{
	typename UpdateMeshType::VertexIterator vi;

	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if(!(*vi).IsD())
			(*vi).C().ColorRamp(minq,maxq,(*vi).Q());
}

static void VertexQualityRamp(UpdateMeshType &m)
{
	std::pair<float,float> minmax = Stat<UpdateMeshType>::ComputePerVertexQualityMinMax(m);
	VertexQualityRamp(m,minmax.first,minmax.second);
}

static void VertexQualityGray(UpdateMeshType &m, const float minq, const float maxq)
{
	typename UpdateMeshType::VertexIterator vi;
	
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if(!(*vi).IsD())
			(*vi).C().SetGrayShade( ((*vi).Q()-minq)/(maxq-minq));
}

static void VertexQualityGray(UpdateMeshType &m)
{
	std::pair<float,float> minmax = Stat<UpdateMeshType>::ComputePerVertexQualityMinMax( m);
	VertexQualityGray(m,minmax.first,minmax.second);
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

//Reduces the mesh to two colors according to a threshold.
static int Thresholding(UpdateMeshType &m, float threshold, Color4b c1 = Color4<unsigned char>::Black, Color4b c2 = Color4<unsigned char>::White, const bool ProcessSelected=false)
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

        if(value<=threshold) (*vi).C() = c1;
        else (*vi).C() = c2;
        ++counter;
      }
    }
  }
  return counter;
}

//Computes the lightness value for a specified color. lightness = 0.5*(Max(R,G,B)+Min(R,G,B))
static float ComputeLightness(Color4b c)
{
  float min_rgb = (float)math::Min(c[0],c[1],c[2]);
  float max_rgb = (float)math::Max(c[0],c[1],c[2]);
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

//Apply Contrast filter to the mesh with the given contrast factor.
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

//Performs contrast operations on color, i.e expands or compress the histogram around
//the midpoint value.  NewValue = (OldValue - 128) ◊ factor + 128
static Color4b ColorMul(Color4b c, float factor)
{
  return Color4b( ValueMul(c[0], factor), ValueMul(c[1], factor), ValueMul(c[2], factor), 1);
}

static int ValueMul(int value, float factor)
{
  return math::Clamp<int>((int)((value - 128)*factor + 128), 0, 255);
}

//Apply  Brightness and Contrast  filter to the mesh, with the given contrast factor and brightness amount.
static int BrightnessContrast(UpdateMeshType &m, float brightness, float contrast, const bool ProcessSelected=false)
{
  int counter=0;
  VertexIterator vi;
  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        (*vi).C() = ColorBrightnessContrast((*vi).C(),brightness,contrast);
        ++counter;
      }
    }
  }
  return counter;
}

//Performs contrast and brightness operations on color, i.e NewValue = (OldValue - 128) ◊ contrast + 128 + amount
//The result is clamped just one time after all computations; this get a more accurate result.
// The formula used here is the one of GIMP.
static Color4b ColorBrightnessContrast(Color4b c, float brightness, float contrast)
{
  return Color4b( ValueBrightnessContrast(c[0], brightness, contrast),
									ValueBrightnessContrast(c[1], brightness, contrast),
									ValueBrightnessContrast(c[2], brightness, contrast), 1 );
}

static int ValueBrightnessContrast(unsigned char ivalue, float brightness, float contrast)
{
	float value = float(ivalue)/255.0f;
  if (brightness < 0.0)  value = value * ( 1.0 + brightness);
                    else value = value + ((1.0 - value) * brightness);
	value = (value - 0.5) * (tan ((contrast + 1) * M_PI/4) ) + 0.5;
	return math::Clamp<int>(255.0*value, 0, 255);
}

//Invert the colors of the mesh.
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

//Invert the rgb component of the color
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

//computes the standard gamma transformation on a given color, according to NewVal = OldVal^(1/gamma)
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

//useful bit masks for RGB channels, used for Levels filter.
enum rgbChMask {ALL_CHANNELS = 7, RED_CHANNEL = 4, GREEN_CHANNEL = 2, BLUE_CHANNEL = 1, NO_CHANNELS = 0 };

//Adjusts color levels of the mesh. Filter can be applied to all RGB channels or to each channel separately.
//in_min, gamma and in_max are respectively the black point, the gray point and the white point.
//out_min and out_max are the output level for black and white respectively.
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

//Performs levels transformation on each channel set to 1 in the rgbMask.
static Color4b ColorLevels(Color4b c, float gamma, float in_min, float in_max, float out_min, float out_max, unsigned char rgbMask)
{
  unsigned char r = c[0], g = c[1], b = c[2];
  if(rgbMask & RED_CHANNEL) r = ValueLevels(c[0], gamma, in_min, in_max, out_min, out_max);
  if(rgbMask & GREEN_CHANNEL) g = ValueLevels(c[1], gamma, in_min, in_max, out_min, out_max);
  if(rgbMask & BLUE_CHANNEL) b = ValueLevels(c[2], gamma, in_min, in_max, out_min, out_max);
  return Color4b(r, g, b, 255);
}

//Transform on levels
static int ValueLevels(int value, float gamma, float in_min, float in_max, float out_min, float out_max)
{
  float fvalue = value/255.0f;
  // normalize
  fvalue = math::Clamp<float>(fvalue - in_min, 0.0f, 1.0f) / math::Clamp<float>(in_max - in_min, 1.0f/255.0f, 1.0f);
  // transform gamma
  fvalue = powf(fvalue,1/gamma);
  // rescale range
  fvalue = fvalue * (out_max - out_min) + out_min;
  //back in interval [0,255] and clamp
  return math::Clamp<int>((int)(fvalue * 255), 0, 255);
}

//Colors the mesh. Color is blended to the mesh with the given intensity.
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

//Perform colourisation operation. For each channel C: newC = origC + intensity * (newC - origC)
static Color4b ColorApplyDiff(Color4b old_color, Color4b new_color, float intensity)
{
  return Color4b( ValueApplyDiff(old_color[0],new_color[0],intensity), ValueApplyDiff(old_color[1],new_color[1],intensity), ValueApplyDiff(old_color[2], new_color[2],intensity), 255);
}

static int ValueApplyDiff(int old_value, int new_value, float intensity)
{
  return  math::Clamp<int>((int)(old_value + intensity * (new_value - old_value)), 0, 255);
}

//An useful ENUM to hold all desaturation methods.
enum DesaturationMethods {M_LIGHTNESS = 0, M_LUMINOSITY = 1, M_AVERAGE = 2};

//Desaturates the mesh according the selected method. Method belongs to DesaturationMethods's ENUM.
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

//Desature the color. Ausiliary functions to calculate lightness/luminosity/average.
static Color4b ColorDesaturate(Color4b c, int method)
{
  switch(method){
    case M_LIGHTNESS:{
      int val = (int)ComputeLightness(c);
      return Color4b( val, val, val, 255);
    }
    case M_AVERAGE:{
      int val = (int)ComputeAvgLightness(c);
      return Color4b( val, val, val, 255);
    }
    case M_LUMINOSITY:{
      int val = (int)ComputeLuminosity(c);
      return Color4b( val, val, val, 255);
    }
    default: assert(0);
  }
}

//ausiliary function to compute average lightness. value = (R+G+B)/3
static float ComputeAvgLightness(Color4b c)
{
    return float(c[0]+c[1]+c[2])/3.0f;
}

//ausiliary function to compute luminosity. value = 0.21*R+0.71*G+0.7*B
static float ComputeLuminosity(Color4b c)
{
    return float(0.2126f*c[0]+0.7152f*c[1]+0.0722f*c[2]);
}

//Equalize the histogram of colors. It can equalize any combination of rgb channels or
//it can work on lightness.
static int Equalize(UpdateMeshType &m, unsigned int rgbMask, const bool ProcessSelected=false)
{
  //declares , resets and set up 4 histograms, for Red, Green, Blue and Lightness
  Histogramf Hl, Hr, Hg, Hb;
  Hl.Clear(); Hr.Clear(); Hg.Clear(); Hb.Clear();
	Hl.SetRange(0, 255, 255); Hr.SetRange(0, 255, 255); Hg.SetRange(0, 255, 255); Hb.SetRange(0, 255, 255);

  int counter=0;
  VertexIterator vi;

  //Scan the mesh to build the histograms
  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, put it in the histograms
      {
        float v = ComputeLightness((*vi).C())+0.5; //compute and round lightness value
        Hl.Add(v); Hr.Add((float)(*vi).C()[0]); Hg.Add((float)(*vi).C()[1]); Hb.Add((float)(*vi).C()[2]);
      }
    }
  }

  //for each histogram, compute the cumulative distribution function, and build a lookup table
	int cdf_l[256], cdf_r[256], cdf_g[256], cdf_b[256];
	cdf_l[0] = Hl.BinCount(0); cdf_r[0] = Hr.BinCount(0); cdf_g[0] = Hg.BinCount(0); cdf_b[0] = Hb.BinCount(0);
	for(int i=1; i<256; i++){
    cdf_l[i] = Hl.BinCount(float(i)) + cdf_l[i-1];
    cdf_r[i] = Hr.BinCount(float(i)) + cdf_r[i-1];
    cdf_g[i] = Hg.BinCount(float(i)) + cdf_g[i-1];
    cdf_b[i] = Hb.BinCount(float(i)) + cdf_b[i-1];
  }

  //this loop aaplies the transformation to colors
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

//Applies equalization to the components of the color according to rgbmask
static Color4b ColorEqualize(Color4b c, int cdf_l[256], int cdf_r[256], int cdf_g[256], int cdf_b[256], unsigned int rgbMask)
{
  unsigned char r = c[0], g = c[1], b = c[2];
  if(rgbMask == NO_CHANNELS) //in this case, equalization is done on lightness
  {
    int v = ValueEqualize(cdf_l[(int)(ComputeLightness(c)+0.5f)], cdf_l[0], cdf_l[255]);
    return Color4b(v, v, v, 255); //return the equalized gray color
  }
  if(rgbMask & RED_CHANNEL) r = ValueEqualize(cdf_r[c[0]], cdf_r[0], cdf_r[255]); //Equalizes red
  if(rgbMask & GREEN_CHANNEL) g = ValueEqualize(cdf_g[c[1]], cdf_g[0], cdf_g[255]); //Equalizes green
  if(rgbMask & BLUE_CHANNEL) b = ValueEqualize(cdf_b[c[2]], cdf_b[0], cdf_b[255]); //Equalizes blue
  return Color4b(r, g, b, 255); //return the equalized color
}

//Compute the equalized value
static int ValueEqualize(int cdfValue, int cdfMin, int cdfMax)
{
  return int(float((cdfValue - cdfMin)/float(cdfMax - cdfMin)) * 255.0f);
}

//applies the white balance filter. It may works with an auto regulation of white, or based on a user
//color that is supposed to be white.
static int WhiteBalance(UpdateMeshType &m, bool automatic, Color4b userColor, const bool ProcessSelected=false)
{
  Color4b unbalancedWhite;
  float lightness = 0;
  int counter=0;
  VertexIterator vi;

  if(!automatic) unbalancedWhite = userColor;  //no auto regolation required, user has provided a color.
  else  //else, we need to scan the mesh and pick its lighter color...
  {
    for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
    {
      if(!(*vi).IsD()) //if it has not been deleted...
      {
        if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected...
        {
          //the lighter color is selected with an incremental approach...
          float v = ComputeLightness((*vi).C());
          if( v > lightness){
            lightness = v;                //save lightness
            unbalancedWhite = (*vi).C();  //save the color
          }
        }
      }
    }
  }

  //in this loop the transformation is applied to the mesh
  for(vi=m.vert.begin();vi!=m.vert.end();++vi) //scan all the vertex...
  {
    if(!(*vi).IsD()) //if it has not been deleted...
    {
      if(!ProcessSelected || (*vi).IsS()) //if this vertex has been selected, do transormation
      {
        (*vi).C()=ColorWhiteBalance((*vi).C(),unbalancedWhite);
        ++counter;
      }
    }
  }
  return counter;
}

//Balnce the white of the color, applying a correction factor based on the unbalancedWhite color.
static Color4b ColorWhiteBalance(Color4b c, Color4b unbalancedWhite)
{
  //sanity check to avoid division by zero...
  if(unbalancedWhite[0]==0) unbalancedWhite[0]=1;
  if(unbalancedWhite[1]==0) unbalancedWhite[1]=1;
  if(unbalancedWhite[2]==0) unbalancedWhite[2]=1;

  return Color4b(
                 math::Clamp<int>((int)(c[0]*(255.0f/unbalancedWhite[0])), 0, 255),
                 math::Clamp<int>((int)(c[1]*(255.0f/unbalancedWhite[1])), 0, 255),
                 math::Clamp<int>((int)(c[2]*(255.0f/unbalancedWhite[2])), 0, 255),
                 255);
}

static void PerlinColor(MeshType& m, Box3f bbox, float freq, Point3i channelOffsets)
{
    typedef typename MeshType::ScalarType ScalarType;

    Point3<ScalarType> p;
    for(VertexIterator vi = m.vert.begin(); vi!=m.vert.end(); ++vi)
    {
        if(!(*vi).IsD()){
            p = bbox.GlobalToLocal(m.Tr * (*vi).P());   //actual vertex position scaled to bbox
            (*vi).C() = Color4b( int(255*math::Perlin::Noise(channelOffsets[0]+p[0]*freq,channelOffsets[0]+p[1]*freq,channelOffsets[0]+p[2]*freq)),
                                 int(255*math::Perlin::Noise(channelOffsets[1]+p[0]*freq,channelOffsets[1]+p[1]*freq,channelOffsets[1]+p[2]*freq)),
                                 int(255*math::Perlin::Noise(channelOffsets[2]+p[0]*freq,channelOffsets[2]+p[1]*freq,channelOffsets[2]+p[2]*freq)),
                                 255 );
        }
    }
}

static void ColorNoise(MeshType& m, int noiseBits)
{
    if(noiseBits>8) noiseBits = 8;
    if(noiseBits<1) return;

    math::SubtractiveRingRNG randomGen = math::SubtractiveRingRNG(time(NULL));
    for(VertexIterator vi = m.vert.begin(); vi!=m.vert.end(); ++vi)
    {
        if(!(*vi).IsD()){
            (*vi).C()[0] = math::Clamp<int>((*vi).C()[0] + randomGen.generate(int(2*pow(2.0f,noiseBits))) - int(pow(2.0f,noiseBits)),0,255);
            (*vi).C()[1] = math::Clamp<int>((*vi).C()[1] + randomGen.generate(int(2*pow(2.0f,noiseBits))) - int(pow(2.0f,noiseBits)),0,255);
            (*vi).C()[2] = math::Clamp<int>((*vi).C()[2] + randomGen.generate(int(2*pow(2.0f,noiseBits))) - int(pow(2.0f,noiseBits)),0,255);
        }
    }
}

};

}// end namespace
}// end namespace
#endif
