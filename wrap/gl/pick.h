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

****************************************************************************/

#include <vector>
#include <algorithm>

namespace vcg{

template <class TETRA_MESH_TYPE>
class GLPickTetra
{
	typedef typename TETRA_MESH_TYPE::TetraIterator TetraIterator;
	typedef typename TETRA_MESH_TYPE::TetraPointer  TetraPointer;
	typedef typename TETRA_MESH_TYPE::VertexType  VertexType;

public:
static bool PickNearestTetra(int x, int y,TETRA_MESH_TYPE &m, TetraIterator &ti, int width=4, int height=4)
{
 std::vector<TetraPointer> result;
 int val=PickTetra(x,y,m,result,width,height);
 if(val)
 {
  ti=result[0];
	return true;
 }
 ti=0;
 return false; 
}

static int PickTetra(int x, int y, TETRA_MESH_TYPE &m, std::vector<TetraPointer> &result, int width=4, int height=4)
{
	result.clear();
	long hits;	
  int sz=m.tetra.size()*5;
	unsigned int *selectBuf =new unsigned int[sz];
	//  static unsigned int selectBuf[16384];
  glSelectBuffer(sz, selectBuf);
  glRenderMode(GL_SELECT);
  glInitNames();

  /* Because LoadName() won't work with no names on the stack */
  glPushName(-1);
	double mp[16];
	
  int viewport[4];
  glGetIntegerv(GL_VIEWPORT,viewport);
	glMatrixMode(GL_PROJECTION);
	glGetDoublev(GL_PROJECTION_MATRIX ,mp);
	glPushMatrix();
  glLoadIdentity();
  //gluPickMatrix(x, viewport[3]-y, 4, 4, viewport);
  gluPickMatrix(x, y, width, height, viewport);
	glMultMatrixd(mp);

	glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  int tetracnt=0; 
	TetraIterator ti;
	for(ti=m.tetra.begin();ti!=m.tetra.end();++ti)
	{
		if(!(*ti).IsD())
		{
			glLoadName(tetracnt);
			glBegin(GL_TRIANGLES);
			for (int face=0;face<4;face++)
			{
				VertexType *v0=ti->V(Tetra::VofF(face,0));
				VertexType *v1=ti->V(Tetra::VofF(face,1));
				VertexType *v2=ti->V(Tetra::VofF(face,2));
				glVertex(v0->P());
				glVertex(v1->P());
				glVertex(v2->P());
			}
			glEnd();
			tetracnt++;
		}
		
	}

  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
	glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  hits = glRenderMode(GL_RENDER);
	//xstring buf;
	//if (hits <= 0)     return 0;
	std::vector< std::pair<double,unsigned int> > H;
	for(int ii=0;ii<hits;ii++){
		//TRACE("%ui %ui %ui %ui\n",selectBuf[ii*4],selectBuf[ii*4+1],selectBuf[ii*4+2],selectBuf[ii*4+3]);
		H.push_back( std::pair<double,unsigned int>(selectBuf[ii*4+1]/4294967295.0,selectBuf[ii*4+3]));
	}
	std::sort(H.begin(),H.end());
//  if(H.size()>0) TRACE("\n Closest is %i\n",H[0].second);
  result.resize(H.size());
	for(ii=0;ii<hits;ii++){
		TetraIterator ti=m.tetra.begin();
		advance(ti ,H[ii].second);
		result[ii]=&*ti;
	}
	
	delete [] selectBuf;
  return result.size();
}

};
}
