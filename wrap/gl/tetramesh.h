/****************************************************************************
* MeshLab                                                           o o     *
* An extendible mesh processor                                    o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005, 2009                                          \/)\/    *
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

#ifndef __GLWRAPTETRA__
#define __GLWRAPTETRA__
#include <GL/glew.h>
#include <GL/GL.h>
#include <vcg/space/color4.h>
#include <vcg/space/tetra3.h>
#include <wrap/gui/view.h>
#include <wrap/gl/space.h>
#include <wrap/gl/math.h>

namespace vcg {

namespace tetra {

class GLW {
public:
    enum DrawMode  {DMNone, DMSmallTetra,DMFlat,DMWire, DMHidden,DMTransparent,DMFlatWire} ;
    enum NormalMode{NMFlat,NMSmooth, NMUser, NMPerMesh};
    enum ColorMode {CMNone, CMPerMesh,CMUser,CMPerTetra,CMPerVertexF,CMPerVertex};
    enum Hint {HShrinkFactor};
};

template <typename MeshType>
class GlTetramesh:public GLW{


public:

    typedef typename MeshType::TetraType    TetraType;
    typedef typename TetraType::VertexType  VertexType;
    typedef typename VertexType::ScalarType ScalarType;
    typedef typename VertexType::CoordType  CoordType;

    //subclass for clipping planes
    class ClipPlane
    {
    private:

        CoordType D, D0;
        ScalarType dist;

        GLdouble eqn[4];


    public:
        bool active;

        ClipPlane (){active=false;}

        ~ClipPlane (){}

        ClipPlane(CoordType & p0, CoordType & p1, CoordType & p2)
        {
            CoordType N = ((p1-p0)^(p2-p0)).Normalize();

            D  = N;
            D0 = D;

            dist = vcg::Norm((p0 + p1 + p2) / 3.f);
        }

        //set normal of the clipping plane
        void SetD(CoordType d)
        {
            D  = d;
            D0 = d;
        }
        //set the point of the clipping plane
        void SetDist(ScalarType d)
        {
            dist = d;
        }
        bool IsClipped(CoordType p)
        {
            return D.V(0) * p.X() + D.V(1) * p.Y() + D.V(2) * p.Z() - dist > 0;
        }

        void GlClip()
        {
            if (active){
                GLdouble d=-(D.V(0)*P.V(0)+D.V(1)*P.V(1)+D.V(2)*P.V(2));
                eqn[0]=-D.V(0);
                eqn[1]=-D.V(1);
                eqn[2]=-D.V(2);
                eqn[3]=-d;
                glClipPlane(GL_CLIP_PLANE0, eqn);
                glEnable(GL_CLIP_PLANE0);
            }
        }

        void GlDraw()
        {
            const ScalarType w = 50;
            glColor4f(1., 1., 0., 0.3);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glDisable(GL_LIGHTING);

            glBegin(GL_LINES);
            glVertex3f(0, 0, 0);
            glVertex3f(dist, 0, 0);
            glEnd();
            glBegin(GL_TRIANGLES);
            glVertex3f(dist, -w,  w);
            glVertex3f(dist,  w,  w);
            glVertex3f(dist, -w, -w);
            glVertex3f(dist,  w,  w);
            glVertex3f(dist,  w, -w);
            glVertex3f(dist, -w, -w);
            glEnd();

            glDisable(GL_BLEND);
            glEnable(GL_LIGHTING);
        }

        void Transform(vcg::Matrix44<float> & tr)
        {
            D = (tr * D0).Normalize();
        }

        void offsetDist(ScalarType off)
        {
            dist = (off < -dist) ? 0.001f : dist + off;
        }

        bool IsActive()
        {
            return active;
        }

        bool switchActive()
        {
            return active ^= true;
        }
    };

    GlTetramesh(MeshType * m) : _m(m){}
    GlTetramesh( )  {}

    MeshType * _m;
    ClipPlane section;

private:
    ScalarType shrink_factor = 0.98f;


public:

    void SetHint(Hint h, double value){
        switch(h){
        case HShrinkFactor: shrink_factor = value; break;
        }
    }

    void AddClipSection(CoordType p0, CoordType p1, CoordType p2)
    {
        section=ClipPlane(p0,p1,p2);
        section.active=true;
    }

    void ClearClipSection()
    {
        section.active=false;
    }

    template <DrawMode dm,NormalMode nm,ColorMode cm >
    void Draw(){
        switch (dm){
        case DMNone: break;
        case DMSmallTetra: _DrawSmallTetra<cm>();break;
        case DMFlat:       _DrawSurface<dm,nm,cm>();break;
        case DMWire:       _DrawSurface<dm,nm,cm>();break;
        case DMHidden:     _DrawSurface<dm,nm,cm>();break;
        case DMFlatWire:   _DrawFlatWire<nm,cm>(); break;
        case DMTransparent: break;
        }
    }

private:
    template <ColorMode cm >
    void _DrawSmallTetra(){

        glEnable(GL_LIGHT0);
        glEnable(GL_LIGHTING);
        ForEachTetra(*_m, [&] (TetraType & t) {
            if (!t.IsD())
            {
                if (!t.IsS()) //draw as normal
                    _DrawSmallTetra<cm>(t);
                else          //draw in selected mode
                    _DrawSelectedTetra(t);
            }
        });
    }

    template <NormalMode nm,ColorMode cm >
    void _DrawFlatWire(){
        glPushAttrib(0xffff);
        glEnable(GL_COLOR_MATERIAL);
        glEnable(GL_DEPTH);
        glDepthRange(0.001,1.0);
        Draw<DMFlat,nm,cm>();
        glDisable(GL_LIGHTING);
        glColor3f(0.0,0.0,0.0);
        glDepthRange(0.0,0.999);
        Draw<DMHidden,nm,cm>();
        glPopAttrib();
    }


    template <DrawMode dm,NormalMode nm,ColorMode cm >
    void _DrawSurface(){


        glPushAttrib(0xffff);
        glEnable(GL_COLOR_MATERIAL);
        if((dm == DMWire)||(dm ==DMHidden))
        {
            glDisable(GL_LIGHTING);
            glDisable(GL_NORMALIZE);
            glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        }
        else
        {
            glEnable(GL_LIGHT0);
            glEnable(GL_LIGHTING);
            glEnable(GL_NORMALIZE);
            glPolygonMode(GL_FRONT,GL_FILL);
        }

        ForEachTetra(*_m, [&] (TetraType & t) {
            _DrawTetra<dm,nm,cm>(t);
        });

        glPopAttrib();
    }


    void _DrawSelectedTetra(TetraType &t)
    {
        glPushMatrix();
        glPushAttrib(0xffff);
        glDisable(GL_CLIP_PLANE0);
        glDisable(GL_BLEND);
        glDisable(GL_LIGHTING);
        glDisable(GL_NORMALIZE);
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

        glColor3d(1,0,0);

        glBegin(GL_TRIANGLES);
        for (int face=0;face<4;face++)
        {
            glVertex(t.V(Tetra::VofF(face,0))->P());
            glVertex(t.V(Tetra::VofF(face,1))->P());
            glVertex(t.V(Tetra::VofF(face,2))->P());
        }
        glEnd();

        //end drawing
        glPopAttrib();
        glPopMatrix();
    }

    template <DrawMode dm,NormalMode nm,ColorMode cm >
    void _DrawTetra(TetraType &t)
    {
        if((!t.IsD())&&(!t.IsS()))
        {
            if ((dm!=DMWire)&&(dm!=DMHidden))
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
        else
            if((!t.IsD())&&(t.IsS()))
                _DrawSelectedTetra(t);
    }

    template <ColorMode cm >
    void _ChooseColorTetra(TetraType &t)
    {
        if (cm==CMNone)
        {
            if (t.IsS())
                glColor3d(1,0,0);
            else
                glColor3d(0.8f,0.8f,0.8f);
        }
        if (cm == CMPerTetra)
            vcg::glColor(t.C());
        //        else
        //            if(cm == CMPerTetraF)
        //            {
        //                Color4b c;
        //                c = color_tetra(t);
        //                GLint ic[4]; ic[0] = c[0];ic[1] = c[1];ic[2] = c[2];ic[3] = c[3];
        //                glMaterialiv(GL_FRONT,GL_DIFFUSE ,ic);
        //            }
    }

    template <ColorMode cm >
    void _ChooseColorVertex(VertexType &v)
    {
        if (cm!=CMNone)
        {
            //            if(cm == CMPerVertexF)
            //            {
            //                Color4b c;
            //                c = color_vertex(v);
            //                GLint ic[4]; ic[0] = c[0];ic[1] = c[1];ic[2] = c[2];ic[3] = c[3];
            //                glMaterialiv(GL_FRONT,GL_DIFFUSE ,ic);
            //            }
            //            else
            if(cm == CMPerVertex)
                vcg::glColor(v.C());
        }
    }

    template <ColorMode cm >
    void _DrawFaceSmooth(TetraType &t,int face)
    {

        VertexType *v0=t.V(Tetra::VofF(face,0));
        VertexType *v1=t.V(Tetra::VofF(face,1));
        VertexType *v2=t.V(Tetra::VofF(face,2));

        glBegin(GL_TRIANGLES);
        _ChooseColorVertex<cm>(*v0);
        glNormal(v0->N());
        glVertex(v0->P());
        _ChooseColorVertex<cm>(*v1);
        glNormal(v1->N());
        glVertex(v1->P());
        _ChooseColorVertex<cm>(*v2);
        glNormal(v2->N());
        glVertex(v2->P());
        glEnd();
    }

    template < ColorMode cm >
    void _DrawFace(TetraType &t,int face)
    {
        glBegin(GL_TRIANGLES);
        VertexType *v0=t.V(Tetra::VofF(face,0));
        VertexType *v1=t.V(Tetra::VofF(face,1));
        VertexType *v2=t.V(Tetra::VofF(face,2));
        glNormal(vcg::Normal(v0->P(), v1->P(), v2->P()).normalized());
        _ChooseColorVertex<cm>(*v0);
        glVertex(v0->P());
        _ChooseColorVertex<cm>(*v1);
        glVertex(v1->P());
        _ChooseColorVertex<cm>(*v2);
        glVertex(v2->P());
        glEnd();
    }

    template < ColorMode cm >
    void _DrawSmallTetra(TetraType &t)
    {
        CoordType p[4], br;
        br = Tetra::Barycenter(t);

        if (section.active)
        {
            if (section.IsClipped(br))
                return;
            bool border = false;
            bool clipBorder = false;
            for (int i = 0; i < 4; ++i)
            {
                border = border || t.IsB(i);

                CoordType br1 = Tetra::Barycenter(*t.TTp(i));
                clipBorder = clipBorder || section.IsClipped(br1);
            }

            if (!border && !clipBorder)
                return;
        }

        for(int i = 0; i < 4; ++i)
            //            p[i] = t.V(i)->P();
            p[i] = t.V(i)->P() * shrink_factor + br * (1 - shrink_factor);


        _ChooseColorTetra<cm>(t);

        glBegin(GL_TRIANGLES);
        for(int i = 0; i < 4; ++i)
        {
            VertexType *v0=t.V(Tetra::VofF(i,0));
            VertexType *v1=t.V(Tetra::VofF(i,1));
            VertexType *v2=t.V(Tetra::VofF(i,2));


            glNormal(vcg::Normal(v0->P(), v1->P(), v2->P()).normalized());

            _ChooseColorVertex<cm>(*v0);
            glVertex(p[Tetra::VofF(i,0)]);
            _ChooseColorVertex<cm>(*v1);
            glVertex(p[Tetra::VofF(i,1)]);
            _ChooseColorVertex<cm>(*v2);
            glVertex(p[Tetra::VofF(i,2)]);
        }
        glEnd();
    }
};

} // end namespace tetra
} // end nemaspace tri
#endif
