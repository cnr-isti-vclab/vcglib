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
  int level = 0;
  int apatch = -1;
  float geo_error = -1;
  float screen_error = -1;
  int option;

  while((option = getopt(argc, argv, "l:p:g:s:")) != EOF) {
    switch(option) {
    case 'l': level = atoi(optarg); break;
    case 'p': apatch = atoi(optarg); break;
    case 'g': geo_error = (float)atof(optarg); break;
    case 's': screen_error = (float)atof(optarg); break;
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
  nexus.LoadHistory();
  Sphere3f sphere = nexus.sphere;

  if(!init()) {
    cerr << "Could not init SDL window\n";
    return -1;
  }

  bool rotate = true;
  glClearColor(0, 0, 0, 0); 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
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
    gluPerspective(40, 1, 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,6,   0,0,0,   0,1,0);    


    float scale = 2/sphere.Radius();
    glRotatef(alpha, 0, 1, 0);
    if(rotate)
      alpha++;
    if(alpha > 360) alpha = 0;
    glScalef(scale, scale, scale);       
    Point3f center = sphere.Center();
    glTranslatef(-center[0], -center[1], -center[2]);
   
    vector<unsigned int> cells;
    if(apatch != -1) {
      cells.push_back(apatch);
    } else if(geo_error != -1) {
      nexus.ExtractFixed(cells, geo_error);
    } else {
      for(int i = 0; i < nexus.index.size(); i++) {
	if(nexus.index[i].error == 0)
	  cells.push_back(i);
      }
    }
    

    for(unsigned int i = 0; i < cells.size(); i++) {
      unsigned int cell = cells[i];
      Patch patch = nexus.GetPatch(cell);
      
      unsigned int val = cell + 1;
      glColor3ub(((val * 27)%128) + 128, 
		 ((val * 37)%128) + 128, 
		 ((val * 87)%128) + 128);

      glBegin(GL_TRIANGLES);
      unsigned short *f = patch.FaceBegin();      
      for(unsigned int j = 0; j < patch.FaceSize()*3; j+= 3) {
	Point3f &p1 = patch.Vert(f[j]);
	Point3f &p2 = patch.Vert(f[j+1]);
	Point3f &p3 = patch.Vert(f[j+2]);
	Point3f n = ((p2 - p1) ^ (p3 - p1));
	
	glNormal3f(n[0], n[1], n[2]);
	glVertex3f(p1[0], p1[1], p1[2]);
	glVertex3f(p2[0], p2[1], p2[2]);
	glVertex3f(p3[0], p3[1], p3[2]);
      }
      glEnd();

      //drawing borders
      glColor3f(1, 1, 1);

      Border border = nexus.GetBorder(cell);
      glPointSize(2);
      glDisable(GL_LIGHTING);
      glBegin(GL_POINTS);
      for(unsigned int k = 0; k < border.Size(); k++) {
	if(border[k].IsNull()) continue;
	Point3f &p = patch.Vert(border[k].start_vert);
	glVertex3f(p[0], p[1], p[2]);
      }
      glEnd();
      glEnable(GL_LIGHTING);
    }

    
    SDL_GL_SwapBuffers();
  }

        // Clean up

  SDL_Quit();
  return -1;
}


