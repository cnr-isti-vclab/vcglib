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
Revision 1.8  2004/09/30 00:27:42  ponchio
Lot of changes. Backup.

Revision 1.7  2004/09/28 10:26:35  ponchio
Backup

Revision 1.6  2004/09/17 15:25:09  ponchio
First working (hopefully) release.

Revision 1.5  2004/09/16 14:25:16  ponchio
Backup. (lot of changes).

Revision 1.4  2004/08/27 00:38:34  ponchio
Minor changes.

Revision 1.3  2004/07/15 14:32:49  ponchio
Debug.

Revision 1.2  2004/07/05 15:49:39  ponchio
Windows (DevCpp, mingw) port.

Revision 1.1  2004/07/04 15:30:00  ponchio
Changed directory structure.

Revision 1.2  2004/07/04 15:16:01  ponchio
*** empty log message ***

Revision 1.1  2004/07/02 17:41:57  ponchio
Created.

Revision 1.3  2004/07/02 13:03:34  ponchio
*** empty log message ***

Revision 1.2  2004/07/01 21:33:46  ponchio
Added remap reading.

Revision 1.1  2004/06/23 00:10:38  ponchio
Created


****************************************************************************/

#ifdef WIN32
#include "getopt.h"
#else
#include <unistd.h>
#endif

#include <iostream>
using namespace std;

#include <SDL/SDL.h>

#include <GL/gl.h>
#include <GL/glu.h>

#include <apps/nexus/nexusmt.h>
#include <wrap/gui/trackball.h>


using namespace vcg;
using namespace nxs;

bool fullscreen = false;
int width =1024;
int height = 768;

//TrackHand hand;

SDL_Surface *screen = NULL;

bool init() {
  
  if(SDL_Init(SDL_INIT_VIDEO) != 0) {
    return false;
  }

  const SDL_VideoInfo *info = SDL_GetVideoInfo();
  int bpp = info->vfmt->BitsPerPixel;

  SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 5);
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

  int flags = SDL_OPENGL;
  if(fullscreen) 
    flags |= SDL_FULLSCREEN;

  screen = SDL_SetVideoMode(width, height, bpp, flags);
  if(!screen) {
    return false;
  }
  
  SDL_WM_SetIcon(SDL_LoadBMP("inspector.bmp"), NULL);
  SDL_WM_SetCaption(" Inspector", "Inspector");


  glDisable(GL_DITHER);
  glShadeModel(GL_SMOOTH);
  glHint( GL_FOG_HINT, GL_NICEST );
  glEnable(GL_DEPTH_TEST);
  glDepthFunc( GL_LEQUAL );
  glDisable(GL_LIGHTING); 

  glEnableClientState(GL_VERTEX_ARRAY);
  return true;
}




int main(int argc, char *argv[]) {
  enum Mode { SCREEN, GEO };
  int level = 0;
  int apatch = -1;
  float error = 1;
  Mode mode = SCREEN;

  Trackball track;
  int option;

  while((option = getopt(argc, argv, "l:p:g:s:")) != EOF) {
    switch(option) {
    case 'l': level = atoi(optarg); break;
    case 'p': apatch = atoi(optarg); break;
    case 'g': mode = GEO; error = (float)atof(optarg); break;
    case 's': mode = SCREEN; error = (float)atof(optarg); break;
    default: cerr << "Unknown option: " << (char)option << endl;
      return -1;
    }
  }

  if(optind != argc - 1) {
    cerr << "Usage: " << argv[0] << " <nexus file> [options]\n"
	 << " -l <n>: show level n\n"
	 << " -p <n>: show patch n\n"
	 << " -g <e>: extract at geometry error e\n"
	 << " -s <e>: extract at screen error e\n\n";
    return -1;
  }      

  NexusMt nexus;
  if(!nexus.Load(argv[optind])) {
    cerr << "Could not load nexus file: " << argv[1] << endl;
    return -1;
  }
  Sphere3f sphere = nexus.sphere;

  if(!init()) {
    cerr << "Could not init SDL window\n";
    return -1;
  }
  
  //  FrustumPolicy frustum_policy;
  
  
  bool rotate = true;
  bool show_borders = true;
  bool show_colors = true;
  bool show_normals = true;
  glClearColor(0, 0, 0, 0); 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  int quit = 0;
  SDL_Event         event;
  int x, y;
  float alpha = 0;
  while( !quit ) {                    
    while( SDL_WaitEvent( &event ) ){                        
      switch( event.type ) {
      case SDL_QUIT:  quit = 1; break;      
      case SDL_KEYDOWN:                                        
	switch(event.key.keysym.sym) {
	case SDLK_q: exit(0); break;
	case SDLK_b: show_borders = !show_borders;break;
	case SDLK_c: 
	  show_colors = !show_colors;

	  break;
	case SDLK_n: 
	  show_normals = !show_normals;
	  break;

	case SDLK_r:
	case SDLK_SPACE: rotate = !rotate; break;
	  
	case SDLK_MINUS: error *= 0.9; 
	  cerr << "error: " << error << endl; break;
	  
	case SDLK_EQUALS:
	case SDLK_PLUS: error *= 1.1; 
	  cerr << "error: " << error << endl; break;
	}
	//quit = 1;           
	//error++;
	//if(error == 5) error = 0;
	//render.setMaxError(error/10.0);
	break;
      case SDL_MOUSEBUTTONDOWN:       
	x = event.button.x;
	y = height - event.button.y;          
	if(event.button.button == SDL_BUTTON_WHEELUP) {
	  track.MouseWheel(1);
	} else if(event.button.button == SDL_BUTTON_WHEELDOWN) {
	  track.MouseWheel(-1);
	} else 
	  track.MouseDown(x, y, 1);
	//          hand.buttonDown(x, y, 1);
        break;
      case SDL_MOUSEBUTTONUP:          
	x = event.button.x;
	y = height - event.button.y;      
	track.MouseUp(x, y, 1);    
	//          hand.buttonUp(x, y);
	break;
      case SDL_MOUSEMOTION: 
	while(SDL_PeepEvents(&event, 1, SDL_GETEVENT, SDL_MOUSEMOTIONMASK));
	x = event.motion.x;
	y = height - event.motion.y;
	track.MouseMove(x, y);
	//          hand.mouseMove(x, y);        
	break;  
      default: break;
      }
  
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluPerspective(40, 1, 0.1, 100);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      gluLookAt(0,0,5,   0,0,0,   0,1,0);    


      track.GetView();
      track.Apply();

      float scale = 2/sphere.Radius();
      //    glRotatef(alpha, 0, 1, 0);
      //    if(rotate)
      //      alpha++;
      //    if(alpha > 360) alpha = 0;
      glScalef(scale, scale, scale);       
      Point3f center = sphere.Center();
      glTranslatef(-center[0], -center[1], -center[2]);

      nexus.SetMode(NexusMt::DEBUG);
      nexus.SetPolicy(NexusMt::FRUSTUM, error);
      nexus.SetComponent(NexusMt::COLOR, show_colors);
      nexus.SetComponent(NexusMt::NORMAL, show_normals);
   
      nexus.Render();

      /*      vector<unsigned int> cells;
      if(apatch != -1) {
	cells.push_back(apatch);
      } else if(mode == GEO) {
	nexus.ExtractFixed(cells, error);
      } else if(mode == SCREEN) {
	frustum_policy.error = error;
	frustum_policy.GetView();
	nexus.Extract(cells, &frustum_policy);
      } else {
	for(int i = 0; i < nexus.index.size(); i++) {
	  if(nexus.index[i].error == 0)
	    cells.push_back(i);
	}
      }
    
      glColor3f(1, 1, 1);

      for(unsigned int i = 0; i < cells.size(); i++) {
	unsigned int cell = cells[i];
	Patch patch = nexus.GetPatch(cell);

	if(show_color) {
	  unsigned int val = cell + 1;
	  glColor3ub(((val * 27)%128) + 128, 
		     ((val * 37)%128) + 128, 
		     ((val * 87)%128) + 128);
	}
      
	glBegin(GL_TRIANGLES);
	unsigned short *f = patch.FaceBegin();      
	for(unsigned int j = 0; j < patch.nf*3; j+= 3) {
	  Point3f &p1 = patch.Vert(f[j]);
	  Point3f &p2 = patch.Vert(f[j+1]);
	  Point3f &p3 = patch.Vert(f[j+2]);
	  Point3f n = ((p2 - p1) ^ (p3 - p1));
	
	  if(!show_normals) {
	    glNormal3f(n[0], n[1], n[2]);
	    glVertex3f(p1[0], p1[1], p1[2]);
	    glVertex3f(p2[0], p2[1], p2[2]);
	    glVertex3f(p3[0], p3[1], p3[2]);
	  } else {
	    short *n1 = patch.Norm16(f[j]);
	    short *n2 = patch.Norm16(f[j+1]);
	    short *n3 = patch.Norm16(f[j+2]);
	    glNormal3s(n1[0], n1[1], n1[2]);
	    glVertex3f(p1[0], p1[1], p1[2]);
	    glNormal3s(n2[0], n2[1], n2[2]);
	    glVertex3f(p2[0], p2[1], p2[2]);
	    glNormal3s(n3[0], n3[1], n3[2]);
	    glVertex3f(p3[0], p3[1], p3[2]);
	  }

	}
	glEnd();
      }
      if(show_borders) {
	for(unsigned int i = 0; i < cells.size(); i++) {
	  unsigned int cell = cells[i];
	  Patch patch = nexus.GetPatch(cell);
	  //drawing borders
	  glColor3f(1, 1, 1);
	
	  Border border = nexus.GetBorder(cell);
	  glPointSize(4);
	  glDisable(GL_LIGHTING);
	  glDisable(GL_DEPTH_TEST);
	  glBegin(GL_POINTS);
	  for(unsigned int k = 0; k < border.Size(); k++) {
	    if(border[k].IsNull()) continue;
	    Point3f &p = patch.Vert(border[k].start_vert);
	    glVertex3f(p[0], p[1], p[2]);
	  }
	  glEnd();
	  glEnable(GL_DEPTH_TEST);
	  glEnable(GL_LIGHTING);
	}
      }
*/
    
      SDL_GL_SwapBuffers();
    }
  }

  // Clean up

  SDL_Quit();
  return -1;
}


