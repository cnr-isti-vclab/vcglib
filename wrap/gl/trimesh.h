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
Revision 1.7  2005/02/26 12:45:23  spinelli
ripristinata la modalita' di render bbox....

Revision 1.6  2005/01/12 14:39:41  tommyfranken
constructor name was wrong (old class name)

Revision 1.5  2004/12/15 18:45:07  tommyfranken
*** empty log message ***

Revision 1.4  2004/09/30 01:40:39  ponchio
<gl/glew.h>  --> <GL/glew.h>

Revision 1.3  2004/07/13 11:25:57  pietroni
changed order of initial include ( it had problems with extension of openGL)

Revision 1.2  2004/07/12 15:57:33  ganovelli
first draft: it includes glew !




****************************************************************************/

#ifndef __VCG_GLTRIMESH
#define __VCG_GLTRIMESH

#include <queue>
#include <vector>

#include <GL/glew.h>
#include <wrap/gl/space.h>
#include <vcg/space/Color4.h>

namespace vcg {


// classe base di glwrap usata solo per poter usare i vari drawmode, normalmode senza dover 
// specificare tutto il tipo (a volte lunghissimo) 
// della particolare classe glwrap usata.
class GLW
{
public:
	enum DrawMode  {DMNone, DMBox, DMPoints, DMWire, DMHidden, DMFlat, DMSmooth, DMFlatWire, DMRadar, DMLast} ;
	enum NormalMode{NMNone, NMPerVert, NMPerFace, NMPerWedge, NMLast};
	enum ColorMode {CMNone, CMPerMesh, CMPerFace, CMPerVert, CMLast};
	enum TextureMode   {TMNone, TMPerVert, TMPerWedge, TMPerWedgeMulti};
	enum Hint {
		HNUseTriStrip		  = 0x0001,				// ha bisogno che ci sia la fftopology gia calcolata!
//		HNUseEdgeStrip		  = 0x0002,			// 
		HNUseDisplayList	  = 0x0004, 
		HNCacheDisplayList	  = 0x0008,		// Each mode has its dl;
		HNLazyDisplayList	  = 0x0010,			// Display list are generated only when requested 
		HNIsTwoManifold		  = 0x0020,			// There is no need to make DetachComplex before . 
		HNUsePerWedgeNormal	  = 0x0040,		// 
		HNHasFFTopology       = 0x0080,		// E' l'utente che si preoccupa di tenere aggiornata la topologia ff
		HNHasVFTopology       = 0x0100,		// E' l'utente che si preoccupa di tenere aggiornata la topologia vf
		HNHasVertNormal       = 0x0200,		// E' l'utente che si preoccupa di tenere aggiornata le normali per faccia
		HNHasFaceNormal       = 0x0400,		// E' l'utente che si preoccupa di tenere aggiornata le normali per vertice
		HNUseVArray           = 0x0800, 
		HNUseLazyEdgeStrip	  = 0x1000,		// Edge Strip are generated only when requested
		HNUseVBO              = 0x2000		// Use Vertex Buffer Object
	};

	enum Change {
		CHVertex		= 0x01,
		CHNormal		= 0x02,
		CHColor			= 0x04,
		CHFace			= 0x08,
		CHFaceNormal= 0x10,
		CHAll       =  0xff
	};
	enum HintParami {
		HNPDisplayListSize =0		
	};
	enum HintParamf {
		HNPCreaseAngle =0,	// crease angle in radians
		HNPZTwist = 1		// Z offset used in Flatwire and hiddenline modality
	};

	template<class MESH_TYPE>
	class VertToSplit
	{
	public:
		typename MESH_TYPE::face_base_pointer f;
		char z;
		char edge;
		bool newp;
		typename MESH_TYPE::vertex_pointer v;
	};

	// GL Array Elemet
	class GLAElem {
	public : 
		int glmode;
		int len;
		int start;
	};

	
};

template <class MESH_TYPE,  bool partial = false , class FACE_POINTER_CONTAINER = std::vector<MESH_TYPE::FacePointer> > 
class GlTrimesh : public GLW
	{
	public:
		
	MESH_TYPE *m;
	GlTrimesh(){
		m=0;
		dl=0xffffffff;
		h=HNUseLazyEdgeStrip;
		cdm=DMNone;
		ccm=CMNone;
		cnm=NMNone;
    
		SetHintParamf(HNPCreaseAngle,float(M_PI/5));
		SetHintParamf(HNPZTwist,0.00005f);

	}
		typedef MESH_TYPE mesh_type;
		FACE_POINTER_CONTAINER face_pointers;


	unsigned int b[3];
	int h;      // the current hints
	// The parameters of hints
  int   HNParami[8];
	float HNParamf[8];
	void SetHintParami(const HintParami hip, const int value)
		{
			HNParamI[hip]=value;
		}
	int GetHintParami(const HintParami hip) const
		{
			return HNParamI[hip];
		}
	void SetHintParamf(const HintParamf hip, const float value)
		{
			HNParamf[hip]=value;
		}
	float GetHintParamf(const HintParamf hip) const
		{
			return HNParamf[hip];
		}

void SetHint(Hint hn) 
{
	h |= hn;
}
void ClearHint(Hint hn) 
{
	h&=(~hn);
}



	unsigned int dl;  
	std::vector<unsigned int> indices;

	DrawMode cdm; // Current DrawMode
	NormalMode cnm; // Current NormalMode
	ColorMode ccm; // Current ColorMode

void Update(Change c=CHAll)
{
	if(m==0) return;
		if(h&HNUseVArray){

		MESH_TYPE::FaceIterator fi;
		indices.clear();
		for(fi = m->face.begin(); fi != m->face.end(); ++fi)
		{
			indices.push_back((unsigned int)((*fi).V(0) - &(*m->vert.begin())));
			indices.push_back((unsigned int)((*fi).V(1) - &(*m->vert.begin())));
			indices.push_back((unsigned int)((*fi).V(2) - &(*m->vert.begin())));
		}
		
		if(h&HNUseVBO){
		if(!glIsBuffer(b[1]))
				glGenBuffers(2,b);
		glBindBuffer(GL_ARRAY_BUFFER,b[0]);   
		glBufferData(GL_ARRAY_BUFFER_ARB, m->vn * sizeof(MESH_TYPE::VertexType), 
								(char *)&(m->vert[0].P()), GL_STATIC_DRAW_ARB); 
		
		glBindBuffer(GL_ARRAY_BUFFER,b[1]);   
		glBufferData(GL_ARRAY_BUFFER_ARB, m->vn * sizeof(MESH_TYPE::VertexType), 
								(char *)&(m->vert[0].N()), GL_STATIC_DRAW_ARB); 

		}
		glVertexPointer(3,GL_FLOAT,sizeof(MESH_TYPE::VertexType),0);
		glNormalPointer(GL_FLOAT,sizeof(MESH_TYPE::VertexType),0);
	}

	//int C=c;
	//if((C&CHVertex) || (C&CHFace)) {
	//	ComputeBBox(*m);
	//	if(!(h&HNHasFaceNormal)) m->ComputeFaceNormal();
	//	if(!(h&HNHasVertNormal)) m->ComputeVertexNormal();
	//	C= (C | CHFaceNormal);
	//}
	//if((C&CHFace) && (h&HNUseEdgeStrip)) 		 ComputeEdges();
	//if((C&CHFace) && (h&HNUseLazyEdgeStrip)) ClearEdges();
	//if(MESH_TYPE::HasFFTopology())
	//	if((C&CHFace) && (h&HNUseTriStrip)) 		{
	//			if(!(h&HNHasFFTopology)) m->FFTopology();
	//			ComputeTriStrip();
	//		}
	//if((C&CHFaceNormal) && (h&HNUsePerWedgeNormal))		{
	//	  if(!(h&HNHasVFTopology)) m->VFTopology();
	//		CreaseWN(*m,MESH_TYPE::scalar_type(GetHintParamf(HNPCreaseAngle)));
	//}
	//if(C!=0) { // force the recomputation of display list 
	//	cdm=DMNone;
	//	ccm=CMNone;
	//	cnm=NMNone;
	//}
	//if((h&HNUseVArray) && (h&HNUseTriStrip)) 
	//	{
	//	 ConvertTriStrip<MESH_TYPE>(*m,TStrip,TStripF,TStripVED,TStripVEI);
	//	}
}

void Draw(DrawMode dm ,ColorMode cm, TextureMode tm)
{
	switch(dm)
	{
		case	DMNone    : Draw<DMNone    >(cm,tm); break;
		case	DMBox     : Draw<DMBox     >(cm,tm); break;
		case	DMPoints  : Draw<DMPoints  >(cm,tm); break;
		case	DMWire    : Draw<DMWire    >(cm,tm); break;
		case	DMHidden  : Draw<DMHidden  >(cm,tm); break;
		case	DMFlat    : Draw<DMFlat    >(cm,tm); break;
		case	DMSmooth  : Draw<DMSmooth  >(cm,tm); break;
		case	DMFlatWire: Draw<DMFlatWire>(cm,tm); break;
		default : break;
	}
}

template< DrawMode dm >
void Draw(ColorMode cm, TextureMode tm)
{
	switch(cm)
	{
		case	CMNone    : Draw<dm,CMNone   >(tm); break;
		case	CMPerMesh : Draw<dm,CMPerMesh>(tm); break;
		case	CMPerFace : Draw<dm,CMPerFace>(tm); break;
		case	CMPerVert : Draw<dm,CMPerVert>(tm); break;
		default : break;
	}
}

template< DrawMode dm, ColorMode cm >
void Draw(TextureMode tm)
{
	switch(tm)
	{
		case	TMNone          : Draw<dm,cm,TMNone          >(); break;
		case	TMPerVert       : Draw<dm,cm,TMPerVert       >(); break;
		case	TMPerWedge      : Draw<dm,cm,TMPerWedge      >(); break;
		case	TMPerWedgeMulti : Draw<dm,cm,TMPerWedgeMulti >(); break;
		default : break;
	}
}



template< DrawMode dm, ColorMode cm, TextureMode tm>
void Draw()
{
	if(!m) return;
	if((h & HNUseDisplayList)){
				if (cdm==dm && ccm==cm){
						glCallList(dl);
						return;
				}
				else {
					if(dl==0xffffffff) dl=glGenLists(1);
					glNewList(dl,GL_COMPILE);
				}
	}

	glPushMatrix();
	switch(dm)
		{
			case DMNone		  : break;
			case DMBox		  : DrawBBox(cm);break;
			case DMPoints   : DrawPoints<NMPerVert,cm>();break;
			case DMHidden		:	DrawHidden();break;
			case DMFlat			:	DrawFill<NMPerFace,cm,tm>();break;
			case DMFlatWire :	DrawFlatWire<NMPerFace,cm,tm>();break;
			case DMRadar		:	DrawRadar<NMPerFace,cm>();break;
			case DMWire		  :	DrawWire<NMPerVert,cm>();break;
			case DMSmooth   : DrawFill<NMPerVert,cm,tm>();break;
			default : break;
		}
	glPopMatrix();

	if((h & HNUseDisplayList)){
		cdm=dm;
		ccm=cm;
		glEndList();
		glCallList(dl);
	}
}


/*********************************************************************************************/
/*********************************************************************************************/

template <NormalMode nm, ColorMode cm, TextureMode tm>
void DrawFill()
{
	FACE_POINTER_CONTAINER::iterator fp;

	MESH_TYPE::FaceIterator fi;
	

	std::vector<MESH_TYPE::FaceType*>::iterator fip;
	short curtexname=-1;

	if(cm == CMPerMesh) glColor(m->C());
 
	if(h&HNUseVArray) 
	{
		if( (nm==NMPerVert) && ((cm==CMNone) || (cm==CMPerMesh)))
		{
			
			glEnableClientState (GL_NORMAL_ARRAY);
			glEnableClientState (GL_VERTEX_ARRAY);

			if(h&HNUseVBO){
				glBindBuffer(GL_ARRAY_BUFFER,b[1]);   
				glNormalPointer(GL_FLOAT,sizeof(MESH_TYPE::VertexType),0);
				glBindBuffer(GL_ARRAY_BUFFER,b[0]);   
				glVertexPointer(3,GL_FLOAT,sizeof(MESH_TYPE::VertexType),0);
			}
			else
			{
				glNormalPointer(GL_FLOAT,sizeof(MESH_TYPE::VertexType),&(m->vert.begin()->N()[0]));
				glVertexPointer(3,GL_FLOAT,sizeof(MESH_TYPE::VertexType),&(m->vert.begin()->P()[0])); 
			}
		

			glDrawElements(GL_TRIANGLES ,m->fn*3,GL_UNSIGNED_INT, &(*indices.begin()) );
			glDisableClientState (GL_VERTEX_ARRAY);						
			glDisableClientState (GL_NORMAL_ARRAY  );

			return;
		}
	}
	else

 	if(h&HNUseTriStrip) 
	{
		//if( (nm==NMPerVert) && ((cm==CMNone) || (cm==CMPerMesh)))
		//	if(h&HNUseVArray){
		//		glEnableClientState (GL_NORMAL_ARRAY  );
		//		glNormalPointer(GL_FLOAT,sizeof(MESH_TYPE::VertexType),&(m->vert[0].cN()));
		//		glEnableClientState (GL_VERTEX_ARRAY);
		//		glVertexPointer(3,GL_FLOAT,sizeof(MESH_TYPE::VertexType),&(m->vert[0].cP()));
		//		std::vector<GLAElem>::iterator vi;
		//		for(vi=TStripVED.begin();vi!=TStripVED.end();++vi)
		//					glDrawElements(vi->glmode ,vi->len,GL_UNSIGNED_SHORT,&TStripVEI[vi->start] );
		//					
		//		glDisableClientState (GL_NORMAL_ARRAY  );
		//		glDisableClientState (GL_VERTEX_ARRAY);
		//		return;
		//	}

		//std::vector< MESH_TYPE::VertexType *>::iterator vi;
		//glBegin(GL_TRIANGLE_STRIP);
		//if(nm == NMPerFace) fip=TStripF.begin();

		//for(vi=TStrip.begin();vi!=TStrip.end(); ++vi){
		//	if((*vi)){
		//		if(nm==NMPerVert) glNormal((*vi)->cN());
		//		if(nm==NMPerFace) glNormal((*fip)->cN());
		//		glVertex((*vi)->P());
		//		}
		//	else
		//		{
		//			glEnd();
		//			glBegin(GL_TRIANGLE_STRIP);
		//		}
		//	if(nm == NMPerFace) ++fip;
		//	}
		//glEnd();

	}
 	else
	 	{
 
			glBegin(GL_TRIANGLES);
			if(partial) 
				fp = face_pointers.begin();
			else
				fi = m->face.begin();

			while( (partial)?(fp!=face_pointers.end()):(fi!=m->face.end()))
			{
 				MESH_TYPE::FaceType & f = (partial)?(*(*fp)): *fi;

				if(!f.IsD()){
//ATLTRACE("Drawing face\n");
				//if(tm==TMPerWedgeMulti)
				//	if(f.WT(0).n() != curtexname)
				//		{
				//		curtexname=(*fi).WT(0).n();
				//		glEnd();
				//		glBindTexture(GL_TEXTURE_2D,TMId(curtexname));
				//		glBegin(GL_TRIANGLES);
				//		}
				if(nm == NMPerFace) glNormal(f.cN());
				if(cm == CMPerFace) glColor(f.C());
	
				
				if(nm==NMPerVert) glNormal(f.V(0)->cN());
				if(nm==NMPerWedge)glNormal(f.WN(0));
				if(cm==CMPerVert) glColor(f.V(0)->C());
//				if(tm==TMPerVert) glTexCoord(f.V(0)->T());
				if( (tm==TMPerWedge)|| (tm==TMPerWedgeMulti)) glTexCoord(f.WT(0).t(0));
				glVertex(f.V(0)->P());

				if(nm==NMPerVert) glNormal(f.V(1)->cN());
				if(nm==NMPerWedge)glNormal(f.WN(1));
				if(cm == CMPerVert) glColor(f.V(1)->C());
//				if(tm==TMPerVert) glTexCoord(f.V(1)->T());
				if( (tm==TMPerWedge)|| (tm==TMPerWedgeMulti)) glTexCoord(f.WT(1).t(0));
				glVertex(f.V(1)->P());

				if(nm==NMPerVert) glNormal(f.V(2)->cN());
				if(nm==NMPerWedge)glNormal(f.WN(2));
				if(cm == CMPerVert) glColor(f.V(2)->C());
	//			if(tm==TMPerVert) glTexCoord(f.V(2)->T());
				if( (tm==TMPerWedge)|| (tm==TMPerWedgeMulti)) glTexCoord(f.WT(2).t(0));
				glVertex(f.V(2)->P());
				}
				if(partial) ++fp; else ++fi;
			}
			glEnd();
			
		}

}

template<NormalMode nm, ColorMode cm>
void DrawPoints()
{
	MESH_TYPE::VertexIterator vi;
	glBegin(GL_POINTS);
	if(cm==CMPerMesh) glColor(m->C());
		
	for(vi=m->vert.begin();vi!=m->vert.end();++vi)if(!(*vi).IsD())
	{
			if(nm==NMPerVert) glNormal((*vi).cN());			
			if(cm==CMPerVert) glColor((*vi).C());			
			glVertex((*vi).P());			
	}
	glEnd();
}	


void DrawHidden()
{
	const float ZTWIST=HNParamf[HNPZTwist];
	glDepthRange(ZTWIST,1.0f);
	glDisable(GL_LIGHTING);
	glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
	DrawFill<NMNone,CMNone,TMNone>();
	glEnable(GL_LIGHTING);
	glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
	
	glDepthRange(0.0f,1.0f-ZTWIST);
	DrawWire<NMPerVert,CMNone>();
	glDepthRange(0,1.0f);
}

template <NormalMode nm, ColorMode cm, TextureMode tm>
void DrawFlatWire()
{
	const float ZTWIST=HNParamf[HNPZTwist];
	glDepthRange(ZTWIST,1.0f);
	DrawFill<nm,cm,tm>();
	glDepthRange(0.0f,1.0f-ZTWIST);
	glPushAttrib(GL_CURRENT_BIT);
	glColor3f(.3f,.3f,.3f);
	DrawWire<NMPerVert,CMNone>();
	glPopAttrib();
	glDepthRange(0,1.0f);
}	

template <NormalMode nm, ColorMode cm>
void DrawRadar()
{
		const float ZTWIST=HNParamf[HNPZTwist];
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDepthMask(0);
	glDepthRange(ZTWIST,1.0f);

	if (cm == CMNone) 
		glColor4f(0.2f, 1.0f, 0.4f, 0.2f);
//	DrawFill<nm,cm,TMNone>();
	Draw<DMFlat,CMNone,TMNone>();

	glDepthMask(1);
	glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
//	DrawFill<nm,cm,TMNone>();
	Draw<DMFlat,CMNone,TMNone>();

	glDepthRange(0.0f,1.0f-ZTWIST);
	glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
	glColor4f(0.1f, 1.0f, 0.2f, 0.6f);
	Draw<DMWire,CMNone,TMNone>();
	glDisable(GL_BLEND);
	glDepthRange(0,1.0f);

}



#ifdef GL_TEXTURE0_ARB
// Multitexturing nel caso voglia usare due texture unit.
void DrawTexture_NPV_TPW2()
{
	unsigned int texname=(*(m->face.begin())).WT(0).n(0);
	glBindTexture(GL_TEXTURE_2D,TMId(texname));
	MESH_TYPE::FaceIterator fi;
	glBegin(GL_TRIANGLES);
	for(fi=m->face.begin();fi!=m->face.end();++fi)if(!(*fi).IsD()){
		  if(texname!=(*fi).WT(0).n(0))	{
				texname=(*fi).WT(0).n(0);
				glEnd();
				glBindTexture(GL_TEXTURE_2D,TMId(texname));
				glBegin(GL_TRIANGLES);
			}
			glMultiTexCoordARB(GL_TEXTURE0_ARB, (*fi).WT(0).t(0));
			glMultiTexCoordARB(GL_TEXTURE1_ARB, (*fi).WT(0).t(0));
			glNormal((*fi).V(0)->N());
			glVertex((*fi).V(0)->P());
			
			glMultiTexCoordARB(GL_TEXTURE0_ARB, (*fi).WT(1).t(0));
			glMultiTexCoordARB(GL_TEXTURE1_ARB, (*fi).WT(1).t(0));
			glNormal((*fi).V(1)->N());
			glVertex((*fi).V(1)->P());
			
			glMultiTexCoordARB(GL_TEXTURE0_ARB, (*fi).WT(2).t(0));
			glMultiTexCoordARB(GL_TEXTURE1_ARB, (*fi).WT(2).t(0));
			glNormal((*fi).V(2)->N());
			glVertex((*fi).V(2)->P());
	}
	glEnd();
}

#endif


int MemUsed()
{
	int tot=sizeof(GLWrap);
	tot+=sizeof(edge_type)*edge.size();
	tot+=sizeof(MESH_TYPE::VertexType *) * EStrip.size();
	tot+=sizeof(MESH_TYPE::VertexType *) * TStrip.size();
	tot+=sizeof(MESH_TYPE::FaceType *)   * TStripF.size();
	return tot;
}

private:

template <NormalMode nm, ColorMode cm>
void DrawWire()
{
	//if(!(h & (HNUseEdgeStrip | HNUseLazyEdgeStrip) ) )
	//	{
			glPushAttrib(GL_POLYGON_BIT);
			glPolygonMode(GL_FRONT_AND_BACK ,GL_LINE);
			DrawFill<nm,cm,TMNone>();
			glPopAttrib();
	//	}
	//else 
	//	{
//			if(!HasEdges()) ComputeEdges();
			
			//if(cm==CMPerMesh)	glColor(m->C());
			//std::vector< MESH_TYPE::VertexType *>::iterator vi;
			//glBegin(GL_LINE_STRIP);
			//for(vi=EStrip.begin();vi!=EStrip.end(); ++vi){
			//	if((*vi)){
			//			glNormal((*vi)->N());
			//			glVertex((*vi)->P());
			//		}
			//	else
			//		{
			//			glEnd();
			//			glBegin(GL_LINE_STRIP);
			//		}
			//}
			//glEnd();
	//	}			
}

void DrawBBox(ColorMode cm)
{
	if(cm==CMPerMesh) glColor(m->C());
	glBoxWire(m->bbox);
}


};// end class

/*
Crease Angle
Assume che:
la mesh abbia la topologia ff
la mesh non abbia complex (o se li aveva fossero stati detached)
Abbia le normali per faccia normalizzate!!


Prende una mesh e duplica tutti gli edge le cui normali nelle facce incidenti formano un angolo maggiore 
di <angle> (espresso in rad).
foreach face
 foreach unvisited vert vi
   scan the star of triangles around vi duplicating vi each time we encounter a crease angle.

the new (and old) vertexes are put in a std::vector that is swapped with the original one at the end.
*/
// uncomment one of the following line to enable the Verbose Trace for Crease
#define VCTRACE (void)0
//#define VCTRACE TRACE

template<class MESH_TYPE>
void Crease(MESH_TYPE &m, typename MESH_TYPE::scalar_type angleRad)
{
	assert(m.HasFFTopology());
	MESH_TYPE::scalar_type cosangle=Cos(angleRad);

	std::vector<GLW::VertToSplit<MESH_TYPE> > SPL;
	std::vector<MESH_TYPE::VertexType> newvert;
	newvert.reserve(m.fn*3);
	// indica se un il vertice z della faccia e' stato processato 
	enum {VISITED_0= MESH_TYPE::FaceType::USER0,      
				VISITED_1= MESH_TYPE::FaceType::USER0<<1,
				VISITED_2= MESH_TYPE::FaceType::USER0<<2} ;
	int vis[3]={VISITED_0,VISITED_1,VISITED_2}; 
						
	int _t2=clock();
	MESH_TYPE::FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
		if(!(*fi).IsD())	(*fi).Supervisor_Flags()&= (~(VISITED_0 | VISITED_1 | VISITED_2));
	
	for(fi=m.face.begin();fi!=m.face.end();++fi)
		if(!(*fi).IsD())
		 for(int j=0;j<3;++j)
			 if(!((*fi).Supervisor_Flags() & (vis[j])))
			 {
				 //VCTRACE("Face %i Spinning around vertex %i\n",fi-m.face.begin(), (*fi).V(j)-m.vert.begin());
				//(*fi).Supervisor_Flags() |= vis[j];
				MESH_TYPE::hedgepos_type he(&*fi,j,(*fi).V(j));
				MESH_TYPE::hedgepos_type she=he;
				MESH_TYPE::face_base_pointer nextf;
				GLW::VertToSplit<MESH_TYPE>  spl;
				spl.newp=false;
				spl.edge=-1;
					
				//Primo giro per trovare un bordo da cui partire
				do {
					he.FlipF();
					he.FlipE();
					if(he.IsBorder()) break;
				}	while(he!=she);
       if(he==she) // non c'e'bordi allora si cerca un crease
			 {
				 do {
						he.FlipF();
						he.FlipE();
						nextf=he.f->F(he.z);
						MESH_TYPE::scalar_type ps=nextf->N()*he.f->N();
						if(ps<cosangle)   break;		
						int vz=0;
						if(he.v == he.f->V(he.z)) vz=he.z;
						if(he.v == he.f->V((he.z+1)%3)) vz=(he.z+1)%3;
						assert((he.f->Supervisor_Flags() & vis[vz] )==0);
					}	while(he!=she);
				}
				he.FlipE(); 
				
				she=he;
				newvert.push_back(*(*fi).V(j));
				MESH_TYPE::vertex_pointer curvert=&newvert.back();
//				VCTRACE("Starting from face %i edge %i vert %i \n",he.f-m.face.begin(), he.z, he.v-m.vert.begin());

				// Secondo giro in cui si riempie il vettore SPL con tutte le info per fare i nuovi vertici
				do{
					//TRACE("     -- spinning face %i edge %i vert %i\n",he.f-m.face.begin(), he.z, he.v-m.vert.begin());
					spl.v=curvert;
					spl.f=he.f;
					spl.z=-1;
					if(he.v == he.f->V(he.z)) spl.z=he.z;
					if(he.v == he.f->V((he.z+1)%3)) spl.z=(he.z+1)%3;
					assert(spl.z>=0);
					//VCTRACE("     -- spinning face vert %i Adding spl face %i vert %i\n",\
					//	he.v-m.vert.begin(),	spl.f-m.face.begin(), spl.z );
					assert((spl.f->Supervisor_Flags() & vis[spl.z] )==0);
					spl.f->Supervisor_Flags() |= vis[spl.z];
					SPL.push_back(spl);
					spl.newp=false;
					spl.edge=-1;
					if(he.IsBorder()) break; 
					nextf=he.f->F(he.z);
					if(nextf==she.f) break;
					MESH_TYPE::scalar_type ps=nextf->N()*he.f->N();
					if(ps<cosangle){
//									VCTRACE("splitting faces %i-%i edge %i vert %i\n",nextf-m.face.begin(),he.f-m.face.begin(), he.z, he.v-m.vert.begin());
									newvert.push_back(*(he.v));
									curvert=&newvert.back();
									spl.newp=true;
									//spl.edge=he.z;
					}
					he.FlipF();
					if(spl.newp) spl.edge=he.z;
					he.FlipE();
					
				}while(he!=she);
			 }
	assert(SPL.size()==m.fn*3);

	std::vector<GLW::VertToSplit<MESH_TYPE> >::iterator vsi;
	for(vsi=SPL.begin();vsi!=SPL.end();++vsi)
		{
			(*vsi).f->V((*vsi).z)=(*vsi).v;
			if((*vsi).newp){
				assert((*vsi).edge>=0 && (*vsi).edge<3);
					if(!(*vsi).f->IsBorder( (*vsi).edge) )      
						(*vsi).f->Detach((*vsi).edge);
					
			}
		}

	m.vert.math::Swap(newvert);
	m.vn=m.vert.size();
}

/*
	Secondo tipo di crease angle. ha bisogno del per wedge normal
	e delle adiacence per vertice faccia gia fatte;
	Assume che le normali per faccia siano gia'state fatte (se ci sono)
 */

template<class MESH_TYPE>
void CreaseWN(MESH_TYPE &m, typename MESH_TYPE::scalar_type angle)
{
	if(!(MESH_TYPE::FaceType::OBJ_TYPE & MESH_TYPE::FaceType::OBJ_TYPE_WN) )
		{
			assert(0); // You needs a mesh with faces having per wedge normals
			return;
		}

	MESH_TYPE::scalar_type cosangle=Cos(angle);
	
	MESH_TYPE::FaceIterator fi;

	// Clear the per wedge normals
	for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())	
		{	
			(*fi).WN(0)=MESH_TYPE::vectorial_type(0,0,0);
			(*fi).WN(1)=MESH_TYPE::vectorial_type(0,0,0);
			(*fi).WN(2)=MESH_TYPE::vectorial_type(0,0,0);
		}

	MESH_TYPE::FaceType::vectorial_type nn;
		
	for(fi=m.face.begin();fi!=m.face.end();++fi)		 if(!(*fi).IsD())	
		{
			nn=(*fi).cN();
			for(int i=0;i<3;++i)
			 {
		 		VEdgePosB<MESH_TYPE::FaceType::face_base> x;
						for(x.f = (*fi).V(i)->Fp(), x.z = (*fi).V(i)->Zp(); x.f!=0; x.NextF() )	{
							assert(x.f->V(x.z)==(*fi).V(i));
							if(x.f->cN()*nn > cosangle)		x.f->WN(x.z)+=nn;
						}
				}
		}
	
}
} // end namespace 

 #endif
