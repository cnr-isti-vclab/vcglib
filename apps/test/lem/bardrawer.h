#ifndef __GLWRAPBAR__
#define __GLWRAPBAR__

#include<GL/GL.h>
#include <vcg/space/color4.h>
#include <vcg/physics/methods/lem/lem.h>
#include <wrap/gl/space.h>

namespace vcg {

template < class STL_BAR_CONT >
class GLWrapBar{
public:
	/// The bar container
  typedef STL_BAR_CONT BarContainer;
  /// The bar type 
  typedef typename STL_BAR_CONT::value_type BarType;	
  /// The type of bar iterator
  typedef typename STL_BAR_CONT::iterator BarIterator;
  ///the type of coordinates
  typedef typename BarType::CoordType CoordType;
  ///the type of scalar
  typedef typename CoordType::ScalarType ScalarType;

 /* typedef typename MESH_TYPE MeshType;
  typedef typename MeshType::FaceType FaceType;
  typedef typename MeshType::VertexType VertexType;
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::ScalarType ScalarType;
  typedef typename vcg::LemSolver<MeshType,FaceType> LemSolver; 
  typedef typename LemSolver::BarType BarType;*/

	GLWrapBar(BarContainer & _b):Bars(_b){}

	BarContainer & Bars;	
	
	public:


    void Draw()
    {

        BarIterator Bi;
        glLineWidth(3.f);
        glPushAttrib(GL_CURRENT_BIT|GL_ENABLE_BIT );
        glDisable(GL_NORMALIZE);
        glDisable(GL_LIGHTING);
        
        
       for (Bi=Bars.begin();Bi<Bars.end();Bi++)
       {  
         CoordType direction=CoordType(0,0,0);
		 ScalarType verse=1.f;
		 //invert verse of axis bar
		 if (Bi->D>2)
			verse=-1.f;

		 direction.V(Bi->D%3)=verse;

		 if (Bi->D==0)
			 glColor3d(1,1,1);
		 else
		 if (Bi->D==1)
			 glColor3d(1,0,0);
		 else
		 if (Bi->D==2)
			 glColor3d(0,1,0);
		 else
		 if (Bi->D==3)
			 glColor3d(0,0,1);
		 else
		 if (Bi->D==4)
			 glColor3d(0,1,1);
		 else
			 glColor3d(1,0,1);
		
		 if (Bi->IsTouched())
			glColor3d(0,0,0);

         glBegin(GL_LINE_STRIP);
            vcg::glVertex(Bi->P);
            vcg::glVertex(Bi->P+(direction*Bi->L));
			// vcg::glVertex(Bi->P+direction);
          glEnd();
		/*glBegin(GL_LINE_STRIP);
            vcg::glVertex(Bi->V0->P());
            vcg::glVertex(Bi->V1->P());
          glEnd();*/
	   }
       glPopAttrib();       
      }

	template<class MESH_TYPE>
    void DrawMesh(MESH_TYPE *m)
    {
       MESH_TYPE::FaceIterator Fi;
       glPushAttrib(GL_COLOR_BUFFER_BIT);
           
       glColor4d(0.8,0.8,0.8,0.9);
       for (Fi=m->face.begin();Fi<m->face.end();Fi++)
       {
	    glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);
	    glEnable(GL_LIGHTING);
        glEnable(GL_NORMALIZE);
        glBegin(GL_TRIANGLES);
          glNormal(Fi->NormalizedNormal());
          glVertex(Fi->V(0)->P());
          glVertex(Fi->V(1)->P());
          glVertex(Fi->V(2)->P());
        glEnd();

		glDisable(GL_BLEND);
		glDisable(GL_LIGHTING);
        glDisable(GL_NORMALIZE);
		glColor3d(0,0,0);
		glBegin(GL_LINE_LOOP);
          glVertex(Fi->V(0)->P());
          glVertex(Fi->V(1)->P());
          glVertex(Fi->V(2)->P());
        glEnd();
       }
       glPopAttrib(); 
    }
};

}
#endif