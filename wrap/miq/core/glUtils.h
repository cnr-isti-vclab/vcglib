#ifndef MIQ_GL_UTILS
#define MIQ_GL_UTILS
//#include <wrap/gl/space.h>
#include "vertex_indexing.h"

class Miq_Gl_Utils
{
public:

   template <class VertexIndexingType>
   static void GLDrawVertexIndexing(VertexIndexingType &VI)
   {
       typedef typename VertexIndexingType::ScalarType ScalarType;

       glPushAttrib(GL_ALL_ATTRIB_BITS);
       glEnable(GL_COLOR_MATERIAL);
       glDisable(GL_LIGHTING);
       glDepthRange(0,0.999);
       typename VertexIndexingType::ScalarType size=5;
       glPointSize(size);
       vcg::glColor(vcg::Color4b(0,255,0,255));
       glBegin(GL_POINTS);
       for (unsigned int i=0;i<VI.duplicated.size();i++)
       {
           assert(!VI.duplicated[i]->IsD());
           vcg::glVertex(VI.duplicated[i]->P());
       }
       glEnd();
       glPopAttrib();
   }

   template <class MeshType>
   static void DrawFlippedFacesIfSelected(MeshType &Tmesh)
   {
       glPushAttrib(GL_ALL_ATTRIB_BITS);
       glDisable(GL_LIGHTING);
       glLineWidth(1.5);
       glDepthRange(0,0.998);
       vcg::glColor(vcg::Color4b(255,0,0,255));
       glBegin(GL_LINES);
       for (unsigned int i=0;i<Tmesh.face.size();i++)
       {
           if (!Tmesh.face[i].IsS())continue;
           for (int k=0;k<3;k++)
           {
               vcg::glVertex(Tmesh.face[i].P0(k));
               vcg::glVertex(Tmesh.face[i].P1(k));
           }
       }
       glEnd();
       glPopAttrib();
   }

   template <class QuadrangulatorType>
   static void QuadGLDrawIntegerVertices(QuadrangulatorType &Quadr)
   {
       glPushAttrib(GL_ALL_ATTRIB_BITS);
       glDisable(GL_LIGHTING);
       glEnable(GL_COLOR_MATERIAL);
       glDisable(GL_TEXTURE_2D);
       glPointSize(8);

       glDepthRange(0,0.997);
       /*glColor3d(1,0,0);*/
       glBegin(GL_POINTS);
       for (int i=0;i<Quadr.IntegerVertex.size();i++)
       {
           typename QuadrangulatorType::TriVertexType* v=Quadr.IntegerVertex[i];
           typename QuadrangulatorType::CoordType pos=v->P();
           if (v->IsV())
               glColor3d(1,0,0);
           else
               glColor3d(1,1,0);
           glVertex(pos);
       }
       glEnd();

       glPopAttrib();
   }

   template <class QuadrangulatorType>
   static void GLDrawIntegerLines(QuadrangulatorType &Quadr)
   {
       glPushAttrib(GL_ALL_ATTRIB_BITS);
       glDisable(GL_LIGHTING);
       glEnable(GL_COLOR_MATERIAL);
       glDisable(GL_TEXTURE_2D);
       glLineWidth(2);

       glColor3d(0,1,0);
       glDepthRange(0,0.998);

       for (int i=0;i<Quadr.IntegerLines.size();i++)
       {
           typename QuadrangulatorType::TriFaceType *f=Quadr.IntegerLines[i].first;
           int edge=Quadr.IntegerLines[i].second;
           typename QuadrangulatorType::TriVertexType* v0=f->V0(edge);
           typename QuadrangulatorType::TriVertexType* v1=f->V1(edge);
           glBegin(GL_LINES);
           glVertex(v0->P());
           glVertex(v1->P());
           glEnd();
       }

       glPopAttrib();
   }

   template <class QuadrangulatorType>
   static void GLDrawPolygons(QuadrangulatorType &Quadr)
   {
       glPushAttrib(GL_ALL_ATTRIB_BITS);
       glEnable(GL_LIGHTING);
       glEnable(GL_COLOR_MATERIAL);
       glDisable(GL_TEXTURE_2D);
       glColor3d(0.7,0.8,0.9);
       //glFrontFace(GL_CW);
       glDepthRange(0,0.998);
       for (unsigned int i=0;i<Quadr.polygons.size();i++)
       {
           glBegin(GL_POLYGON);
           for (unsigned int j=0;j<Quadr.polygons[i].size();j++)
           {
               typename QuadrangulatorType::TriVertexType* v=Quadr.polygons[i][j];
               glNormal(v->N());
               glVertex(v->P());
           }
           glEnd();
       }

       glDepthRange(0,0.997);
       glDisable(GL_LIGHTING);
       glEnable(GL_COLOR_MATERIAL);
       glColor3d(0,0,0);
       for (unsigned int i=0;i<Quadr.polygons.size();i++)
       {
           glBegin(GL_LINE_LOOP);
           for (unsigned int j=0;j<Quadr.polygons[i].size();j++)
           {
               typename QuadrangulatorType::TriVertexType* v=Quadr.polygons[i][j];
               glVertex(v->P());
           }
           glEnd();
       }

       glPopAttrib();
   }

   template <class PolyMesh>
   static void GLDrawPolygonalMesh(PolyMesh &polymesh)
   {
       glPushAttrib(GL_ALL_ATTRIB_BITS);
       glEnable(GL_LIGHTING);
       glEnable(GL_COLOR_MATERIAL);
       glDisable(GL_TEXTURE_2D);
       glColor3d(0.7,0.8,0.9);
       //glFrontFace(GL_CW);
       glDepthRange(0,0.998);
       for (unsigned int i=0;i<polymesh.face.size();i++)
       {
           glBegin(GL_POLYGON);
           for (int j=0;j<polymesh.face[i].VN();j++)
           {
               typename PolyMesh::VertexType* v=polymesh.face[i].V(j);
               glNormal(v->N());
               glVertex(v->P());
           }
           glEnd();
       }

       glDepthRange(0,0.997);
       glDisable(GL_LIGHTING);
       glEnable(GL_COLOR_MATERIAL);
       glColor3d(0,0,0);
       for (unsigned int i=0;i<polymesh.face.size();i++)
       {
           glBegin(GL_LINE_LOOP);
           for (int j=0;j<polymesh.face[i].VN();j++)
           {
               typename PolyMesh::VertexType* v=polymesh.face[i].V(j);
               glNormal(v->N());
               glVertex(v->P());
           }
           glEnd();
       }
       glPopAttrib();
   }
};
#endif
