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
Revision 1.23  2004/12/01 18:46:21  ponchio
Microchanges.

Revision 1.22  2004/11/28 04:16:19  ponchio
*** empty log message ***

Revision 1.21  2004/11/28 01:23:26  ponchio
Fixing borders... let's hope.

Revision 1.20  2004/11/18 18:30:14  ponchio
Using baricenters... lotsa changes.

Revision 1.19  2004/10/30 20:17:03  ponchio
Fixed big patches problem.

Revision 1.18  2004/10/21 13:40:16  ponchio
Debugging.

Revision 1.17  2004/10/21 12:22:21  ponchio
Small changes.

Revision 1.16  2004/10/19 01:23:02  ponchio
Daily backup (fragment...)

Revision 1.15  2004/10/15 16:45:27  ponchio
Vbo added.

Revision 1.14  2004/10/14 13:52:02  ponchio
Small changes.

Revision 1.13  2004/10/14 13:41:34  ponchio
Added statistics.

Revision 1.12  2004/10/09 14:46:47  ponchio
Windows porting small changes.

Revision 1.11  2004/10/04 16:49:54  ponchio
Daily backup. Preparing for compression.

Revision 1.10  2004/10/01 16:54:57  ponchio
Daily backup.

Revision 1.9  2004/09/30 23:56:33  ponchio
Backup (added strips and normals)

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

#include <apps/nexus/nexusmt.h>

#include <iostream>
using namespace std;

#include <SDL/SDL.h>

#ifdef WIN32
#include <windows.h>
#endif

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <wrap/gui/trackball.h>
#include "watch.h"


using namespace vcg;
using namespace nxs;

bool fullscreen = false;
int width =1024;
int height = 768;

void gl_print(float x, float y, char *str);

SDL_Surface *screen = NULL;

bool init(const string &str) {
  
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
  SDL_WM_SetCaption(str.c_str(), str.c_str());


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
  int level = 0;
  int apatch = -1;
  float error = 4;

  Trackball track;
//  int option;

  if(argc != 2) {
    cerr << "Usage: " << argv[0] << " <nexus file>\n";    return -1;
  }      

  NexusMt nexus;
  if(!nexus.Load(argv[1])) {
    cerr << "Could not load nexus file: " << argv[1] << endl;
    return -1;
  }
  Sphere3f sphere = nexus.sphere;

  if(!init(argv[1])) {
    cerr << "Could not init SDL window\n";
    return -1;
  }
  
  //  FrustumPolicy frustum_policy;

  cerr << "Commands: \n"
    " q: quit\n"
    " s: screen error extraction\n"
    " g: geometry error extraction\n"
    " p: draw points\n"
    " t: show statistics\n"
    " r: toggle realtime mode (TODO)\n"
    " b: increase memory buffer\n"
    " B: decrease memory buffer\n"
    " d: debug mode (show patches colored)\n"
    " f: flas shading mode\n"
    " m: smooth mode\n"
    " c: show colors\n"
    " n: show normals\n"
    " u: rotate model\n"
    " -: decrease error\n"
    " +: increase error (= too)\n";
  
  Watch watch;
  
  bool rotate = false;
  bool show_borders = false;
  bool show_colors = true;
  bool show_normals = true;
  bool show_statistics = true;
  bool extract = true;
  
  NexusMt::MetricKind metric;
  NexusMt::Mode mode = NexusMt::SMOOTH;
  unsigned int ram_size = 640000;

  nexus.SetError(error);
  nexus.SetExtractionSize(ram_size);   
  nexus.SetMetric(NexusMt::FRUSTUM);    
  if(!nexus.InitGL()) {
    cerr << "Could not init glew.\n";
  }
  
  glClearColor(0, 0, 0, 0); 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  int quit = 0;
  SDL_Event event;
  int x, y;
  float alpha = 0;
  bool redraw = false;
  float fps = 0;
  float tframe = 0;
  
  while( !quit ) {   
    bool first = true;
    SDL_WaitEvent(&event);
    while( first || SDL_PollEvent( &event ) ){                        
      first = false;
      switch( event.type ) {
      case SDL_QUIT:  quit = 1; break;      
      case SDL_KEYDOWN:                                        
	switch(event.key.keysym.sym) {
	case SDLK_RCTRL:
	case SDLK_LCTRL: 
	  track.ButtonDown(Trackball::KEY_CTRL); break;
	case SDLK_q: exit(0); break;
	case SDLK_b: show_borders = !show_borders; break;
	case SDLK_e: extract = !extract; break;
	case SDLK_c: show_colors = !show_colors; break;
	case SDLK_n: show_normals = !show_normals; break;
	  //	case SDLK_9: nexus.patches->ram_size *= 0.8f; break;
	  //	case SDLK_0: nexus.patches->ram_size *= 1.2f; break;

  case SDLK_LEFT: 
    ram_size *= 0.7; 
    nexus.SetExtractionSize(ram_size);   
    cerr << "Max extraction ram size: " << ram_size << endl; break;
  case SDLK_RIGHT: 
    ram_size *= 1.5; 
    nexus.SetExtractionSize(ram_size);   
    cerr << "Max extraction ram size: " << ram_size << endl; break;
  
	case SDLK_s: metric = NexusMt::FRUSTUM; break;
	case SDLK_p: mode = NexusMt::POINTS; nexus.SetMode(mode); break;
	case SDLK_d: mode = NexusMt::PATCHES; nexus.SetMode(mode); break;
	case SDLK_f: mode = NexusMt::FLAT; nexus.SetMode(mode); break;
	case SDLK_m: mode = NexusMt::SMOOTH; nexus.SetMode(mode); break;

	case SDLK_r:
	case SDLK_SPACE: rotate = !rotate; break;
	  
	case SDLK_MINUS: 
    error *= 0.9f;
    nexus.SetError(error);
	  cerr << "Error: " << error << endl; break;
	  
	case SDLK_EQUALS:
	case SDLK_PLUS: 
    error *= 1.1f; 
    nexus.SetError(error);
	  cerr << "Error: " << error << endl; break;
	}
	break;
      case SDL_KEYUP: 
	switch(event.key.keysym.sym) {
	case SDLK_RCTRL:
	case SDLK_LCTRL:
	  track.ButtonUp(Trackball::KEY_CTRL); break;
	}
	break;
      case SDL_MOUSEBUTTONDOWN:   
	x = event.button.x;
	y = height - event.button.y;          
#ifdef SDL_BUTTON_WHEELUP
	if(event.button.button == SDL_BUTTON_WHEELUP) 
	  track.MouseWheel(1);
	else if(event.button.button == SDL_BUTTON_WHEELDOWN) 
	  track.MouseWheel(-1);
	else 
 #endif
    if(event.button.button == SDL_BUTTON_LEFT)
	  track.MouseDown(x, y, Trackball::BUTTON_LEFT);
	else if(event.button.button == SDL_BUTTON_RIGHT)
	  track.MouseDown(x, y, Trackball::BUTTON_RIGHT);
        break;
      case SDL_MOUSEBUTTONUP:          
	x = event.button.x;
	y = height - event.button.y; 
	if(event.button.button == SDL_BUTTON_LEFT)
	  track.MouseUp(x, y, Trackball::BUTTON_LEFT);
	else if(event.button.button == SDL_BUTTON_RIGHT)
	  track.MouseUp(x, y, Trackball::BUTTON_RIGHT);     
	break;
      case SDL_MOUSEMOTION: 
	while(SDL_PeepEvents(&event, 1, SDL_GETEVENT, SDL_MOUSEMOTIONMASK));
	x = event.motion.x;
	y = height - event.motion.y;
	track.MouseMove(x, y);
	break;  
      case SDL_VIDEOEXPOSE:
      default: break;
      }
      redraw = true;
    }

    if(!redraw) continue;
    redraw = false;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, 1, 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,5,   0,0,0,   0,1,0);    
    
    glRotatef(alpha, 0, 1, 0);
    if(rotate) {
      alpha++;
      if(alpha > 360) alpha = 0;
      SDL_Event redraw;
      redraw.type = SDL_VIDEOEXPOSE;
      SDL_PushEvent(&redraw);
    }
    
    
    track.GetView();
    track.Apply();
    
    float scale = 2/sphere.Radius();
    
    glScalef(scale, scale, scale);       
    Point3f center = sphere.Center();
    glTranslatef(-center[0], -center[1], -center[2]);
    
    Point3f &p = nexus.sphere.Center();
    float r = nexus.sphere.Radius();

    glColor3f(0.8f, 0.8f, 0.8f);
        
     
    //nexus.SetPolicy(policy, error);
    nexus.SetComponent(NexusMt::COLOR, show_colors);
    nexus.SetComponent(NexusMt::NORMAL, show_normals);

    static vector<unsigned int> cells;    
    watch.Start();
    if(extract) {
      //      nexus.patches.Flush();
      
      nexus.metric->GetView();
      //      nexus.policy.Init();
      nexus.tri_total = 0;
      nexus.tri_rendered = 0;
      nexus.Extract(cells);      
    } 
    nexus.Draw(cells);

    /*    if(show_borders) {
      for(unsigned int i = 0; i < cells.size(); i++) {
	Border border = nexus.GetBorder(cells[i]);
	Patch &patch = nexus.GetPatch(cells[i]);
	glPointSize(4);
	glColor3f(1.0f, 1.0f, 1.0f);
	glBegin(GL_POINTS);
	for(unsigned int b = 0; b < border.Size(); b++) {
	  Link &link = border[b];
	  Point3f &p = patch.Vert(link.start_vert);
	  glVertex3f(p[0], p[1], p[2]);
	}
	glEnd();
	glPointSize(1);
      }
      }*/

    //cerr Do some reporting:
    if(show_statistics) {
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glLoadIdentity();
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();
      gluOrtho2D(0, 1, 0, 1);


      glDisable(GL_DEPTH_TEST);
      glDisable(GL_LIGHTING);
      char buffer[1024];
      glColor3f(1.0f, 1.0f, 1.0f);

      /*      sprintf(buffer, "Ram size : %.3fMb (max)   %.3fMb (cur)", 
	      nexus.patches->ram_size * nexus.chunk_size/(float)(1<<20), 
	      nexus.patches->ram_used * nexus.chunk_size/(float)(1<<20));
	      gl_print(0.03, 0.12, buffer);*/

      sprintf(buffer, "Extr size: %.3fMb(max)   %.3fMb(cur)",
	      nexus.extraction_max * nexus.chunk_size/(float)(1<<20), 
	      nexus.extraction_used * nexus.chunk_size/(float)(1<<20));
      gl_print(0.03, 0.09, buffer);
      
      sprintf(buffer, "Vbo size : %.3fMb(max)   %.3fMb(cur)",
	      nexus.patches.vbo_max * nexus.chunk_size/(float)(1<<20), 
	      nexus.patches.vbo_used * nexus.chunk_size/(float)(1<<20));
      gl_print(0.03, 0.06, buffer);

      sprintf(buffer, "Triangles: %.2fK (tot)   %.2fK (vis)    "
                      "%.3f time    %.2f FPS",
	      nexus.tri_total/(float)(1<<10),
	      nexus.tri_rendered/(float)(1<<10),
	      tframe, 1/tframe);
      gl_print(0.03, 0.03, buffer);


      /* cerr << "Ram flushed: " << nexus.patches.ram_flushed << endl;
      	cerr << "Ram readed: " << nexus.patches.ram_readed << endl;*/
      nexus.patches.ram_flushed = 0;
      nexus.patches.ram_readed = 0;
      glEnable(GL_DEPTH_TEST);
      glEnable(GL_LIGHTING);
      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
      glPopMatrix();
    }
    
    SDL_GL_SwapBuffers();
    tframe = watch.Time();
    
  }

  // Clean up

  SDL_Quit();
  return -1;
}

void gl_print(float x, float y, char *str) {
  glRasterPos2f(x, y);
  int len = (int)strlen(str);
  for(int i = 0; i < len; i++) 
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]);
}


