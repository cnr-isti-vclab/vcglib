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
Revision 1.2  2004/07/11 22:13:30  cignoni
Added GPL comments


****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <wrap/gl/space.h>


#include <wrap/callback.h>
#include <wrap/gui/trackball.h>
#include <vcg/simplex/vertex/with/vcvn.h>
#include <vcg/simplex/vertex/with/vcvn.h>
#include <vcg/simplex/face/with/fcfn.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/trimesh/base.h>
#include<wrap/io_trimesh/export_ply.h>
#include<wrap/io_trimesh/import_ply.h>
#include<vcg/complex/trimesh/update/normal.h>
#include<vcg/complex/trimesh/update/bounding.h>

#include "visshader.h"
using namespace vcg;
using namespace std;




// Vertex, Face, Mesh and Grid definitions.
class MyEdge;
class AFace;
class AVertex   : public VertexVCVN< float ,MyEdge,AFace > {};
class AFace     : public FaceFCFN< AVertex,MyEdge,AFace > {};
class AMesh     : public tri::TriMesh< vector<AVertex>, vector<AFace> > {};

///////// Global ////////

int SampleNum=32;
int WindowRes=800;

bool SwapFlag=false;

float lopass=0,hipass=1,Gamma=1;
bool LightFlag=true;
bool ColorFlag=true;

Trackball Q;

int ScreenH,ScreenW;
float ViewAngle=45;

class TimeOracle 
{
public:
	time_t start;
	time_t cur;

	char buf[128];

	char const *TimeToEndStr(double perc)
	{
		time(&cur);
		double diff=difftime(cur,start);
    diff= diff/perc - diff;

		int hh=diff/3600;
		int mm=(diff-hh*3600)/60;
		int ss=diff-hh*3600-mm*60;
		sprintf(buf,"%02i:%02i:%02i",hh,mm,ss);
		return buf;
	}

	void Start(){time(&start);};
};

bool cb(const char *buf)
{
	printf(buf);
	return true;
}

// prototypes
void SaveTexturedGround();

void Draw(AMesh &mm)
{
  AMesh::FaceIterator fi;
  glBegin(GL_TRIANGLES);
  for(fi=mm.face.begin();fi!=mm.face.end();++fi)
  {
    
    glNormal((*fi).V(0)->N()); glColor((*fi).V(0)->C());  glVertex((*fi).V(0)->P());
    glNormal((*fi).V(1)->N()); glColor((*fi).V(1)->C());  glVertex((*fi).V(1)->P());
    glNormal((*fi).V(2)->N()); glColor((*fi).V(2)->C());  glVertex((*fi).V(2)->P());
  }
  glEnd();
}

AMesh m;
VertexVisShader<AMesh> Vis(m);

string OutNameMsh;

	
/*  Called when the window is first opened and whenever 
 *  the window is reconfigured (moved or resized).
 */
void  ViewReshape(GLsizei w, GLsizei h)
{
	ScreenW=w; ScreenH=h;
	glMatrixMode (GL_PROJECTION);   
	glLoadIdentity (); 
	gluPerspective(ViewAngle,(float)w/(float)h,.1,10000);
	glViewport (0, 0, (GLsizei) w, (GLsizei) h); 
}

void  ViewDisplay (void)
{
  glMatrixMode (GL_PROJECTION);   
  glLoadIdentity (); 
  gluPerspective(ViewAngle,1,.1,10);
  glMatrixMode (GL_MODELVIEW);    
  glLoadIdentity ();  
  glTranslatef(0,0,-4);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  Q.GetView();
  Q.Apply();
  Q.Draw();
  float d = 2.0/m.bbox.Diag();
  glScalef(d, d, d);
  glTranslate(-m.bbox.Center());
  if(LightFlag) glEnable(GL_LIGHTING);
          else glDisable(GL_LIGHTING);
  if(ColorFlag) glEnable(GL_COLOR_MATERIAL);
          else glDisable(GL_COLOR_MATERIAL);
  Draw(m);
  glutSwapBuffers();
}
	
void ViewSpecialKey(int , int , int )
{
  glutPostRedisplay();
}
void Toggle(bool &flag) {flag = !flag;}
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
void ViewKey(unsigned char key, int , int )
{
	Point3f dir;
  switch (key) {
  case 27: exit(0);   	break;
	case 'l' :
		lopass=lopass+.05; printf("Lo %f, Hi %f Gamma %f\n",lopass,hipass,Gamma); 
		Vis.MapVisibility(Gamma,lopass,hipass);
		break;
	case 'L' :
		lopass=lopass-.05; printf("Lo %f, Hi %f Gamma %f\n",lopass,hipass,Gamma); 
		Vis.MapVisibility(Gamma,lopass,hipass);
		break;
	case 'h' :
		hipass=hipass-.05; printf("Lo %f, Hi %f Gamma %f\n",lopass,hipass,Gamma); 
		Vis.MapVisibility(Gamma,lopass,hipass);
		break;
	case 'H' :
		hipass=hipass+.05; printf("Lo %f, Hi %f Gamma %f\n",lopass,hipass,Gamma); 
		Vis.MapVisibility(Gamma,lopass,hipass);
		break;
	case 'g' :
		Gamma=Gamma-.05; printf("Lo %f, Hi %f Gamma %f\n",lopass,hipass,Gamma); 
		Vis.MapVisibility(Gamma,lopass,hipass);
		break;
	case 'G' :
		Gamma=Gamma+.05; printf("Lo %f, Hi %f Gamma %f\n",lopass,hipass,Gamma); 
		Vis.MapVisibility(Gamma,lopass,hipass);
		break;
	case 'c' : 
		Vis.ComputeUniform(SampleNum,cb); 
		Vis.MapVisibility(Gamma,lopass,hipass);
		break;
  case ' ' : {
    Point3f dir = Q.camera.ViewPoint();
    printf("ViewPoint %f %f %f\n",dir[0],dir[1],dir[2]);
    dir.Normalize();
    dir=Inverse(Q.track.Matrix())*dir;
    printf("ViewPoint %f %f %f\n",dir[0],dir[1],dir[2]);
    dir.Normalize();
		Vis.ComputeSingle(dir,cb); 
    Vis.MapVisibility(Gamma,lopass,hipass); }
		break;
	case 's' :
		Vis.SmoothVisibility();
		Vis.MapVisibility(Gamma,lopass,hipass); 
		break;
  case 'S' :
    { 
      vcg::tri::io::PlyInfo p; 
      p.mask|=vcg::ply::PLYMask::PM_VERTCOLOR /* | vcg::ply::PLYMask::PM_VERTQUALITY*/ ;
      tri::io::ExporterPLY<AMesh>::Save(m,OutNameMsh.c_str(),false,p);
    }
		break;
  case 'a' : LightFlag = !LightFlag; printf("Toggled Light\n"); break;
  case 'A' : ColorFlag = !ColorFlag; printf("Toggled Color\n"); break;
	}
	glutPostRedisplay(); ;
} 
void ViewMenu(int val)
{
 	 ViewKey(val, 0, 0); 
}
/*********************************************************************/
// TrackBall Functions
/*********************************************************************/

int PressedButton; // What is the button actually pressed? 
int KeyMod;
int GW,GH; // Grandezza della finestra
int B[3]={0,0,0}; // Variabile globale che tiene lo stato dei tre bottoni;


void ViewMouse(int button, int state, int x, int y)
{
  static int glut_buttons=0;

  int m_mask = 0;
  KeyMod=glutGetModifiers();
  if(GLUT_ACTIVE_SHIFT & KeyMod)		m_mask |=  Trackball::KEY_SHIFT;
  if(GLUT_ACTIVE_ALT & KeyMod)			m_mask |=  Trackball::KEY_ALT;
  if(GLUT_ACTIVE_CTRL & KeyMod)			m_mask |=  Trackball::KEY_CTRL;

  if(state == GLUT_DOWN) {
    glut_buttons |= (1<<button);
    Q.MouseDown(x, ScreenH-y, glut_buttons | m_mask);
  } else {
    //glut_buttons &= ~(1<<button);
    glut_buttons = 0;
    m_mask = 0;
    Q.MouseUp(x, ScreenH-y, glut_buttons);
  }
}

void ViewMouseMotion(int x, int y)
{
	Q.MouseMove(x,ScreenH-y);
	glutPostRedisplay();
}

void SetLight()
{
  GLfloat light_ambient0[] = {0.0, 0.0, 0.0, 1.0};  
  GLfloat light_diffuse0[] = {1.0, 1.0, 1.0, 1.0};  
  GLfloat light_position0[] = {0.0, 10.0, 300.0, 0.0};
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
	glEnable(GL_LIGHT0);
}
 
void  ViewInit (void) {
  SetLight();
	Q.Reset();
  Q.radius= 1;
  glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glClearColor (0.8, 0.8, 0.8, 0.0);
	glEnable(GL_NORMALIZE);
	glEnable(GL_LIGHTING);
	
//	glEnable(GL_BLEND);
	glShadeModel(GL_SMOOTH);
//  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);     
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT,GL_DIFFUSE);
	//glColorMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE);
	//glColorMaterial(GL_FRONT,GL_AMBIENT);
  glMateriali(GL_FRONT,GL_SHININESS,0);
  float spec[4]={0,0,0,1};
  glMaterialfv(GL_FRONT,GL_SPECULAR,spec);
  glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	
}



int main(int argc, char** argv)
{
	if(argc<2) {
		printf(
			"shadevis 1.0 \n"__DATE__"\n"
			"Copyright 2003-2004 Visual Computing Lab I.S.T.I. C.N.R.\n"
			"Paolo Cignoni (cignoni@isti.cnr.it)\n\n"
			"Usage: shadevis file.ply [options]\n"
			"Options:\n"
			"     -w#      WindowResolution (default 600)\n"
			"     -n#      Sample Directions (default 32)\n"
			"     -z#      z offset (default 1e-4)\n"
			"     -da #    Cone Direction Angle in degree (default 45)\n"
			"     -dv # # # Cone Direction vector (default 0 0 1)\n"			
			"     -c       Set IsClosed Flag\n"
			"     -f       Flip normal of the model\n"
					 );

		return 1;
	}

	
  int i=1;  
	while(i<argc 	&& (argv[i][0]=='-'))
		{
				switch(argv[i][1])
			{
				case 'n'  : SampleNum = atoi(argv[i]+2); break;
				case 'f'  : SwapFlag=false; break;
        case 'w'  : WindowRes= atoi(argv[i]+2); printf("Set WindowRes to %i\n",WindowRes ); break;
        case 's'  : Vis.SplitNum= atoi(argv[i]+2); printf("Set SplitNum to %i\n",Vis.SplitNum ); break;
        case 'z'  : Vis.ZTWIST = atof(argv[i]+2); printf("Set ZTWIST to %f\n",Vis.ZTWIST ); break;
        default: {printf("Error unable to parse option '%s'\n",argv[i]); exit(0);}
			}

			++i;
		}
	
  
  string basename = argv[i];
	if(!(basename.substr(basename.length()-4)==".ply"))	{
				printf("Error: Unknown file extension %s\n",basename.c_str());
				return 1;
	}
	
  // loading original mesh
  int ret=tri::io::ImporterPLY<AMesh>::Open(m,argv[i]);
	if(ret) {printf("Error unable to open mesh %s\n",argv[i]);exit(-1);}
  tri::UpdateNormals<AMesh>::PerVertexNormalized(m);
  tri::UpdateBounding<AMesh>::Box(m);

  printf("Mesh bbox (%f %f %f)-(%f %f %f)\n\n",m.bbox.min[0],m.bbox.min[1],m.bbox.min[2],m.bbox.max[0],m.bbox.max[1],m.bbox.max[2]);
  OutNameMsh=(string(argv[i]).substr(0,strlen(argv[i])-4));
	OutNameMsh+="_vis.ply";
	
	printf("Mesh       Output filename %s\n",OutNameMsh.c_str());

	printf("Mesh %iv %if bbox Diag %g\n",m.vn,m.fn,m.bbox.Diag());
	
	glutInit(&argc, argv);

  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(WindowRes, WindowRes);
  glutInitWindowPosition (10,10);
  glutCreateWindow ("shadevis - Visual Computing Lab - vcg.isti.cnr.it ");
  glutDisplayFunc(ViewDisplay);
	glutReshapeFunc(ViewReshape); 
	glutKeyboardFunc(ViewKey);	
	glutSpecialFunc(ViewSpecialKey);	
	glutMouseFunc(ViewMouse);
	glutMotionFunc(ViewMouseMotion);

	ViewInit();	
  glewInit();	
	glutMainLoop();

	return(0);
}
