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

#include <apps/nexus/nexus.h>

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
  char file[64];
  if(argc < 2) {
    cerr << "Usage: " << argv[0] << " <nexus file>\n";
    return -1;
  }
  
  Nexus nexus;
  if(!nexus.Load(argv[1])) {
    cerr << "Could not load nexus file: " << argv[1] << endl;
    return -1;
  }
  Sphere3f sphere = nexus.sphere;

  if(!init()) {
    cerr << "Could not init SDL window\n";
    return -1;
  }

  bool rotate = true;
  glClearColor(0, 0, 0, 0); 
  int quit = 0;
  SDL_Event         event;
  int x, y;
  float alpha = 0;
  while( !quit ) {                
    while( SDL_PollEvent( &event ) ){                        
      switch( event.type ) {
        case SDL_QUIT:  quit = 1; break;      
        case SDL_KEYDOWN:                                        
          switch(event.key.keysym.sym) {
	  case SDLK_q: exit(0); break;
	  case SDLK_SPACE: rotate = !rotate;
          }
          //quit = 1;           
          //error++;
          //if(error == 5) error = 0;
          //render.setMaxError(error/10.0);
          break;
        case SDL_MOUSEBUTTONDOWN:       
          x = event.button.x;
          y = event.button.y;          
	  //          hand.buttonDown(x, y, 1);
        break;
        case SDL_MOUSEBUTTONUP:          
          x = event.button.x;
          y = event.button.y;          
	  //          hand.buttonUp(x, y);
          break;
        case SDL_MOUSEMOTION: 
          while(SDL_PeepEvents(&event, 1, SDL_GETEVENT, SDL_MOUSEMOTIONMASK));
          x = event.motion.x;
          y = event.motion.y;
	  //          hand.mouseMove(x, y);        
          break;  
      default: break;
      }
    }
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, 1, 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,3,   0,0,0,   0,1,0);    


    float scale = 2/sphere.Radius();
    glRotatef(alpha, 0, 1, 0);
    if(rotate)
      alpha++;
    if(alpha > 360) alpha = 0;
    glScalef(scale, scale, scale);       
    Point3f center = sphere.Center();
    glTranslatef(-center[0], -center[1], -center[2]);
   


    for(unsigned int i = 7; i < nexus.index.size(); i++) {
      Patch patch = nexus.GetPatch(i);
      
      unsigned int val = i + 1;
      glColor3ub((val * 27)%255, (val * 37)%255, (val * 87)%255);

      glBegin(GL_TRIANGLES);
      unsigned short *f = patch.FaceBegin();      
      for(unsigned int j = 0; j < patch.FaceSize() * 3; j++) {
	Point3f &p = patch.Vert(f[j]);
	glVertex3f(p[0], p[1], p[2]);
      }
      glEnd();
      //      glColor3ub(((val * 27)%255)/2, ((val * 37)%255)/2, ((val * 87)%255)/2);
      glColor3f(1, 1, 1);
      //drawing borders
      Border border = nexus.GetBorder(i);
      glPointSize(4);
      glBegin(GL_POINTS);
      for(unsigned int k = 0; k < border.Size(); k++) {
	Point3f &p = patch.Vert(border[k].start_vert);
	glVertex3f(p[0], p[1], p[2]);
      }
      glEnd();
    }

    
    SDL_GL_SwapBuffers();
  }

        // Clean up

  SDL_Quit();
  return -1;
}


