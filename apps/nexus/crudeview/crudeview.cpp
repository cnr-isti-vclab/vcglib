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

****************************************************************************/

#include <iostream>
using namespace std;

#include <SDL/SDL.h>

#include <GL/gl.h>
#include <GL/glu.h>

#include <wrap/nexus/crude.h>
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
    cerr << "Usage: " << argv[0] << " <crude file>\n";
    return -1;
  }
  
  Crude crude;
  if(!crude.Load(argv[1])) {
    cerr << "Could not load crude file: " << argv[1] << endl;
    return -1;
  }
  Box3f box = crude.GetBox();

  if(!init()) {
    cerr << "Could not init SDL window\n";
    return -1;
  }

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
	  //    hand.glTransform();
	  //    hand.glDraw(ColorUB(255, 100, 100, 200), 3);
    
    float scale = 3/box.Diag();
    //glScalef(0.4, 0.4, 0.4);       
    glRotatef(alpha, 0, 1, 0);
    alpha++;
    if(alpha > 360) alpha = 0;
    glScalef(scale, scale, scale);       
    Point3f center = box.Center();
    glTranslatef(-center[0], -center[1], -center[2]);
//    render.render();
    glColor3f(0, 1, 0);
   
    glBegin(GL_TRIANGLES);

    for(unsigned int i = 0;i < crude.Faces(); i++) {
      Crude::Face &face = crude.GetFace(i);
      for(int k = 0; k < 3; k++) {
	Point3f &p = crude.GetVertex(face[k]);
	glVertex3f(p[0], p[1], p[2]);
      }
    }
    glEnd();
    
    SDL_GL_SwapBuffers();
  }

        // Clean up

  SDL_Quit();
  return -1;
}


