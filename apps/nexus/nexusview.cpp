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
  int option;

  if(argc != 2) {
    cerr << "Usage: " << argv[0] << " <nexus file>\n";
    return -1;
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
    " d: debug mode (show patches colored)\n"
    " m: smooth mode\n"
    " c: show colors\n"
    " n: show normals\n"
    " r: rotate model\n"
    " -: decrease error\n"
    " +: increase error (= too)\n";
  
  
  bool rotate = false;
  bool show_borders = true;
  bool show_colors = true;
  bool show_normals = true;
  NexusMt::Mode mode = NexusMt::SMOOTH;
  NexusMt::PolicyKind policy = NexusMt::FRUSTUM;
  
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
	case SDLK_c: show_colors = !show_colors; break;
	case SDLK_n: show_normals = !show_normals; break;

	case SDLK_s: policy = NexusMt::FRUSTUM; break;
	case SDLK_p: mode = NexusMt::POINTS; break;
	case SDLK_d: mode = NexusMt::DEBUG; break;
	case SDLK_m: mode = NexusMt::SMOOTH; break;

	case SDLK_r:
	case SDLK_SPACE: rotate = !rotate; break;
	  
	case SDLK_MINUS: error *= 0.9; 
	  cerr << "error: " << error << endl; break;
	  
	case SDLK_EQUALS:
	case SDLK_PLUS: error *= 1.1; 
	  cerr << "error: " << error << endl; break;
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
	if(event.button.button == SDL_BUTTON_WHEELUP) 
	  track.MouseWheel(1);
	else if(event.button.button == SDL_BUTTON_WHEELDOWN) 
	  track.MouseWheel(-1);
	else if(event.button.button == SDL_BUTTON_LEFT)
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

    glColor3f(0.8, 0.8, 0.8);
    nexus.SetMode(mode);
    nexus.SetPolicy(policy, error);
    nexus.SetComponent(NexusMt::COLOR, show_colors);
    nexus.SetComponent(NexusMt::NORMAL, show_normals);
    
    nexus.Render();
    
    SDL_GL_SwapBuffers();
  }

  // Clean up

  SDL_Quit();
  return -1;
}


