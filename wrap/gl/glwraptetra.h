#ifndef __GLWRAPTETRA__
#define __GLWRAPTETRA__

#include<GL/GL.h>
#include <vcg/space/color4.h>
#include <vcg/space/Tetra3.h>
#include <wrap/gl/space.h>

namespace vcg {
class GLW {
public:
  enum DrawMode  {DMNone, DMSmallTetra,DMFlat,DMWire, DMHidden,DMTransparent,DMFlatWire} ;
	enum NormalMode{NMFlat,NMSmooth, NMUser, NMPerMesh};
	enum ColorMode {CMNone, CMPerMesh,CMUser,CMPerTetraF,CMPerVertexF,CMPerVertex};
	enum Hint {HShrinkFactor};
};

template <typename CONT_TETRA>
class GLWrapTetra:public GLW{
public:
  
	typedef typename CONT_TETRA::value_type TetraType; 
  typedef typename TetraType::VertexType VertexType;
  typedef typename VertexType::ScalarType ScalarType;
  typedef typename VertexType::CoordType Point3x;

	GLWrapTetra(CONT_TETRA & _t):tetra(_t){}

	CONT_TETRA	& tetra;	

	private:
	double shrink_factor;

	public:

		void SetHint(Hint h, double value){
			switch(h){
				case HShrinkFactor: shrink_factor = value; break;
				}
			}

  typedef Color4b (*color_func_vertex)(VertexType&v);
	color_func_vertex  color_vertex;

	typedef Color4b (*color_func_tetra)(TetraType&v);
	color_func_tetra  color_tetra;
	

	template <DrawMode dm,NormalMode nm,ColorMode cm >
    void	Draw(){
			switch (dm){
				case DMNone: break;
        case DMSmallTetra:	_DrawSmallTetra<cm>();break;
        case DMFlat:_DrawSurface<dm,nm,cm>();break;	
        case DMWire:_DrawSurface<dm,nm,cm>();break;
				case DMHidden:_DrawSurface<dm,nm,cm>();break;
				case DMFlatWire:_DrawFlatWire<nm,cm>(); break;
        case DMTransparent:break;
				}
			}

private:
template <ColorMode cm >
 void _DrawSmallTetra(){
		Point3x p[4],br;
		CONT_TETRA::iterator it;
    glPushAttrib(0xffffffff);
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
    glPolygonMode(GL_FRONT,GL_FILL);
		glBegin(GL_TRIANGLES);
		for( it = tetra.begin(); it != tetra.end(); ++it)
			if(!(*it).IsD()){
					_DrawSmallTetra<cm>(*it);
			}
		glEnd();
    glPopAttrib();
		}

template <NormalMode nm,ColorMode cm >
void 	_DrawFlatWire(){
		glPushAttrib(0xffff);
		glEnable(GL_DEPTH);
		glDepthRange(0.001,1.0);
		Draw<DMFlat,nm,cm>();
		glDisable(GL_LIGHTING);
		glColor3f(0.0,0.0,0.0);
		glDepthRange(0.0,0.999);
		Draw<DMWire,nm,cm>();
		glPopAttrib();
}


template <DrawMode dm,NormalMode nm,ColorMode cm >
void _DrawSurface(){
		CONT_TETRA::iterator it;

		glPushAttrib(0xffffffff);
		
		if((dm == DMWire)||(dm ==DMHidden))
    {
      glDisable(GL_LIGHTING);
      glDisable(GL_NORMALIZE);
			glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    }
		else
    {
      glEnable(GL_LIGHTING);
      glEnable(GL_NORMALIZE);
			glPolygonMode(GL_FRONT,GL_FILL);
    }

		glBegin(GL_TRIANGLES);
		for( it = tetra.begin(); it != tetra.end(); ++it)
        _DrawTetra<dm,nm,cm>((*it));
	  glEnd();
	  glPopAttrib();
}

template <DrawMode dm,NormalMode nm,ColorMode cm >
void _DrawTetra(TetraType &t)
{
  if(!(t.IsD()))
      {
       _ChooseColorTetra<cm>(t);
       for(int i = 0; i < 4; ++i){
         if (dm == DMWire)
           _DrawFace<cm>(t,i);
         else
         {
           if (t.IsBorderF(i))
           {
              if(nm==NMSmooth)
					      _DrawFaceSmooth<cm>(t,i);
              else
					    if(nm==NMFlat)
						    _DrawFace<cm>(t,i);
           }
         }
       }
      }
}

template <ColorMode cm >
void _ChooseColorTetra(TetraType &t)
{
  if (cm==CMNone)
      glColor3d(0.8,0.8,0.8);
  else
  if(cm == CMPerTetraF)
      {
				 Color4b c;
				 c = color_tetra(t);
				 GLint ic[4]; ic[0] = c[0];ic[1] = c[1];ic[2] = c[2];ic[3] = c[3];
				 glMaterialiv(GL_FRONT,GL_DIFFUSE ,ic);
			}
}

template <ColorMode cm >
void _ChooseColorVertex(VertexType &v)
{
  if (cm!=CMNone)
  {
  if(cm == CMPerVertexF)
      {
				 Color4b c;
				 c = color_vertex(v);
				 GLint ic[4]; ic[0] = c[0];ic[1] = c[1];ic[2] = c[2];ic[3] = c[3];
				 glMaterialiv(GL_FRONT,GL_DIFFUSE ,ic);
			}
  else
  if(cm == CMPerVertex)
						glColor3f(v.C()[0],v.C()[1],v.C()[2]);
  }
}

template <ColorMode cm >
void _DrawFaceSmooth(TetraType &t,int face)
{

  VertexType *v0=t.V(Tetra::VofF(face,0));
  VertexType *v1=t.V(Tetra::VofF(face,1));
  VertexType *v2=t.V(Tetra::VofF(face,2));
  _ChooseColorVertex<cm>(*v0);
  glNormal(v0->N());
  glVertex(v0->P());
  _ChooseColorVertex<cm>(*v1);
  glNormal(v1->N());
  glVertex(v1->P());
   _ChooseColorVertex<cm>(*v2);
  glNormal(v2->N());
  glVertex(v2->P());
}

template <ColorMode cm >
void _DrawFace(TetraType &t,int face)
{
  glNormal(t.N(face));
  VertexType *v0=t.V(Tetra::VofF(face,0));
  VertexType *v1=t.V(Tetra::VofF(face,1));
  VertexType *v2=t.V(Tetra::VofF(face,2));
  _ChooseColorVertex<cm>(*v0);
  glVertex(v0->P());
  _ChooseColorVertex<cm>(*v1);
  glVertex(v1->P());
  _ChooseColorVertex<cm>(*v2);
  glVertex(v2->P());
}

template <ColorMode cm >
void _DrawSmallTetra(TetraType &t)
{
  Tetra3<ScalarType> T=Tetra3<ScalarType>();
  T.P0(0)=t.V(0)->cP();
  T.P1(0)=t.V(1)->cP();
  T.P2(0)=t.V(2)->cP();
  T.P3(0)=t.V(3)->cP();
  Point3x p[4], br;
  br=T.ComputeBarycenter();
	for(int i = 0; i < 4; ++i)
						p[i] = t.V(i)->P()* shrink_factor + br *(1- shrink_factor);
  _ChooseColorTetra<cm>(t);
  for(int i = 0; i < 4; ++i)
    {
				glNormal(t.N(i));
        VertexType *v0=t.V(Tetra::VofF(i,0));
        VertexType *v1=t.V(Tetra::VofF(i,1));
        VertexType *v2=t.V(Tetra::VofF(i,2));
        _ChooseColorVertex<cm>(*v0);
        glVertex(p[Tetra::VofF(i,0)]);
        _ChooseColorVertex<cm>(*v1);
				glVertex(p[Tetra::VofF(i,1)]);
        _ChooseColorVertex<cm>(*v2);
				glVertex(p[Tetra::VofF(i,2)]);
			}
}


};

}
#endif