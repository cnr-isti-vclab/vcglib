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
Revision 1.3  2004/09/17 15:25:09  ponchio
First working (hopefully) release.

Revision 1.2  2004/07/05 15:49:39  ponchio
Windows (DevCpp, mingw) port.

Revision 1.1  2004/07/04 15:30:00  ponchio
Changed directory structure.

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

#ifdef WIN32
#include <windows.h>
#endif

#include <GL/gl.h>
#include <GL/glu.h>

#include <apps/nexus/crude.h>
#include <apps/nexus/vert_remap.h>
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

  if(argc < 2) {
    cerr << "Usage: " << argv[0] << " <crude file>\n";
    return -1;
  }
  
  Crude crude;
  if(!crude.Load(argv[1])) {
    cerr << "Could not load crude file: " << argv[1] << endl;
    return -1;
  }
  Box3f box = crude.GetBox();

  bool vremap = false;
  bool fremap = false;
  
  VFile<unsigned int> face_remap;
  if(face_remap.Load(argv[1] + string(".rmf"))) {
    cerr << "Found face remap.\n";
    fremap = true;
  } else {
    cerr << "Face remap not found.\n";
  }
  
  VertRemap vert_remap;
  if(vert_remap.Load(argv[1])) {
    cerr << "Found vert remap.\n";
    vremap = true;
  }
  

  if(!init()) {
    cerr << "Could not init SDL window\n";
    return -1;
  }

  Trackball track;

  glClearColor(0, 0, 0, 0);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
 
  bool redraw = false;
  bool show_normals = true;
  int quit = 0;
  SDL_Event event;
  int x, y;
  float alpha = 0;
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
	}
	break;
      case SDL_MOUSEBUTTONDOWN:       
	x = event.button.x;
	y = event.button.y;          
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
	y = event.button.y;          
	if(event.button.button == SDL_BUTTON_LEFT)
	  track.MouseUp(x, y, Trackball::BUTTON_LEFT);
	else if(event.button.button == SDL_BUTTON_RIGHT)
	  track.MouseUp(x, y, Trackball::BUTTON_RIGHT);     
	break;
      case SDL_MOUSEMOTION: 
	while(SDL_PeepEvents(&event, 1, SDL_GETEVENT, SDL_MOUSEMOTIONMASK));
	x = event.motion.x;
	y = event.motion.y;
	track.MouseMove(x, y);
	break;  
      default: break;
      }
      redraw = true;
    }
    
    if(!redraw) continue;
    redraw = false;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, 1, 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,3,   0,0,0,   0,1,0);    

    track.GetView();
    track.Apply();

    float scale = 3/box.Diag();
    glScalef(0.4, 0.4, 0.4);       
    //    glRotatef(alpha, 0, 1, 0);
    //    alpha++;
    //    if(alpha > 360) alpha = 0;
    glScalef(scale, scale, scale);       
    Point3f center = box.Center();
    glTranslatef(-center[0], -center[1], -center[2]);
//    render.render();
    glColor3f(0, 1, 0);
   
    glBegin(GL_TRIANGLES);

    for(unsigned int i = 0;i < crude.Faces(); i++) {
      Crude::Face &face = crude.GetFace(i);
      if(fremap) {
	unsigned int val = face_remap[i];
	glColor3ub((val * 27)%255, (val * 37)%255, (val * 87)%255);
      }
      Point3f &p0 = crude.GetVertex(face[0]);
      Point3f &p1 = crude.GetVertex(face[1]);
      Point3f &p2 = crude.GetVertex(face[2]);
	
      if(show_normals) {
	Point3f n = ((p1 - p0) ^ (p2 - p0));
	glNormal3f(n[0], n[1], n[2]);
      }
      glVertex3f(p0[0], p0[1], p0[2]);      
      glVertex3f(p1[0], p1[1], p1[2]);
      glVertex3f(p2[0], p2[1], p2[2]);

    }
    glEnd();
    
    SDL_GL_SwapBuffers();
  }

        // Clean up

  SDL_Quit();
  return -1;
}


