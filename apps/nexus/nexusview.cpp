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
Revision 1.38  2005/02/16 15:52:09  ponchio
qualche opzione in piu' , tolti i grafici

Revision 1.37  2005/02/15 15:55:36  ponchio
aggiunta delle sphere

Revision 1.36  2005/02/14 17:11:08  ponchio
aggiunta delle sphere

Revision 1.35  2005/02/14 14:49:09  ponchio
*** empty log message ***

Revision 1.34  2005/02/14 14:21:24  ponchio
Preload disabled at startap (-p)

Revision 1.33  2005/02/10 09:18:20  ponchio
Statistics.

Revision 1.32  2005/02/03 12:35:01  ponchio
Patch cache -> heap

Revision 1.31  2005/02/01 16:42:30  ponchio
Trigger

Revision 1.30  2005/01/21 17:09:13  ponchio
Porting and debug.

Revision 1.29  2005/01/17 17:35:47  ponchio
Small changes and adding realtime extr.

Revision 1.28  2005/01/14 15:25:29  ponchio
Revolution.

Revision 1.27  2004/12/15 16:37:55  ponchio
Optimizing realtime vis.

Revision 1.26  2004/12/15 13:50:32  ponchio
Optimizing realtime vis.

Revision 1.25  2004/12/15 08:46:16  ponchio
Optimizing realtime vis.

Revision 1.24  2004/12/13 00:44:48  ponchio
Lotsa changes...

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
#ifdef WIN32
#include <wrap/system/getopt.h>
#else
#include <unistd.h>
#endif

#include <deque>
#include <iostream>
#include <fstream>

using namespace std;

#include <SDL/SDL.h>

//this include MUST precede GL includes.
#include <apps/nexus/nexusmt.h>

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

  if(argc < 2) {
    cerr << "Usage: " << argv[0] << " <nexus file> [options]\n";    
    cerr << "-e <error>: set initial target error\n"
	 << "-m <ram>: max ram used\n"
         << "-x <ram>: max extraction size\n"
	 << "-r <ram>: max draw size\n"
	 << "-d <ram>: max disk read per frame\n"
	 << "-p      : no preload\n"
    << "-o namefile: ouput stats";
    return -1;
  }      

  NexusMt nexus;
  if(!nexus.Load(argv[1])) {
    cerr << "Could not load nexus file: " << argv[1] << endl;
    return -1;
  }

  Sphere3f sphere = nexus.sphere;
  Extraction extraction;
  DrawContest contest;
  Stats stats;

 

  bool rotate = false;
  bool show_borders = false;
  bool show_colors = true;
  bool do_render = true;
  bool show_normals = true;
  bool show_statistics = true;
  bool extract = true;
  bool realtime = true;
  bool preload = true;
  bool step = true;
  bool output_stats = false;
  char output_filename[100];
  char window_name [100];
  sprintf(window_name,"%s", argv[1]);
  int option;
  while((option = getopt(argc, argv, "e:m:x:r:d:o:w:h:p:f")) != EOF) {
    switch(option) {
    case 'e': extraction.target_error = atof(optarg); break;
    case 'm': nexus.MaxRam() = atoi(optarg); break;
    case 'x': extraction.extr_max = atoi(optarg); break;
    case 'r': extraction.draw_max = atoi(optarg); break;
    case 'd': extraction.disk_max = atoi(optarg); break;
    case 'o': output_stats = true; sprintf(output_filename,"%s",optarg); break;
    case 'w': width =  atoi(optarg); break;
    case 'h': height =  atoi(optarg); break;
    case 'p': preload = false; nexus.SetPreload(preload); break;
    case 'f': fullscreen = true; break;
    default:
      cerr << "Unknow option.\n"; break;
    }
  }

   if(!init(window_name)) {
    cerr << "Could not init SDL window\n";
    return -1;
  }


  //  FrustumPolicy frustum_policy;

  cerr << "Commands: \n"
    " q          : quit\n"
    " t          : toggle statistics\n"
    " right arrow: increase memory buffer\n"
    " left arrow : decrease memory buffer\n"
    " page up    : increase disk space\n"
    " page down  : increase disk space\n"
    " 0          : decrease extraction size\n"
    " 1          : increase extraction size\n"
    " s: toggle preload\n"
    

    " d: debug mode (show patches colored)\n"
    " f: flat shading mode\n"
    " m: smooth mode\n"
    " p: draw points\n"
    " h: draw bounding spheres\n"
    " w: disable glcalls\n"

    " c: show colors\n"
    " n: show normals\n"
    " r: rotate model\n"
    " -: decrease error\n"
    " +: increase error (= too)\n";
  
  Watch watch;
  
  if(!nexus.InitGL()) {
    cerr << "Could not init glew.\n";
    return -1;
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
  float fps = 0;
  unsigned int nave = 5;
  unsigned int offset = 0;
  deque<float> tframe;
  deque<float> terror;
  unsigned int tlen = 256;

  bool keepdrawing = true;
  
  watch.Start();
  while( !quit ) {           
    unsigned int anything = SDL_PollEvent(&event);      
    if(!anything && !keepdrawing) {
      SDL_WaitEvent(&event);
      anything = true;
    }
    if(anything) {        
      switch( event.type ) {
      case SDL_QUIT:  quit = 1; break;      
      case SDL_KEYDOWN:                                        
	switch(event.key.keysym.sym) {
	case SDLK_RCTRL:
	case SDLK_LCTRL: track.ButtonDown(Trackball::KEY_CTRL); break;
	case SDLK_q: exit(0); break;	
	case SDLK_k: keepdrawing = !keepdrawing; break;
	case SDLK_e: extract = !extract; break;
	case SDLK_c: show_colors = !show_colors; break;
	case SDLK_n: show_normals = !show_normals; break;
	case SDLK_w: do_render = !do_render; break;
	    
	case SDLK_LEFT:     nexus.MaxRam() *= 0.8; break;
	case SDLK_RIGHT:    nexus.MaxRam() *= 1.3; break;
	case SDLK_UP:       extraction.draw_max *= 1.3; break;
	case SDLK_DOWN:     extraction.draw_max *= 0.8; break;
	case SDLK_PAGEUP:   extraction.disk_max *= 1.3; break;
	case SDLK_PAGEDOWN: extraction.disk_max *= 0.8; break;
	case SDLK_0:        extraction.extr_max *= 1.3; break;
	case SDLK_9:        extraction.extr_max *= 0.8; break;
  
	  //  case SDLK_s: metric = NexusMt::FRUSTUM; break;
	case SDLK_p: contest.mode = DrawContest::POINTS; break;
	case SDLK_d: contest.mode = DrawContest::PATCHES; break;
	case SDLK_f: contest.mode = DrawContest::FLAT; break;
	case SDLK_m: contest.mode = DrawContest::SMOOTH; break;
	case SDLK_h: if(contest.attrs&DrawContest::SPHERES)
                  contest.attrs &=~DrawContest::SPHERES;  
                else
                  contest.attrs |=DrawContest::SPHERES;
                break;
	    
	case SDLK_o: realtime = !realtime; break;
	case SDLK_s: preload = !preload; nexus.SetPreload(preload); break;
	case SDLK_t: show_statistics = !show_statistics; break;
	case SDLK_r:
	case SDLK_SPACE: rotate = !rotate; break;
	    
	case SDLK_MINUS: 
	  extraction.target_error *= 0.9f;
	  cerr << "Error: " << extraction.target_error << endl; 
	  break;
	    
	case SDLK_EQUALS:
	case SDLK_PLUS: 
	  extraction.target_error *= 1.1f; 
	  cerr << "Error: " << extraction.target_error << endl; 
	  break;
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
    }
                                                                               

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, width/(float)height, 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,5,   0,0,0,   0,1,0);    
    
    glViewport(0,0,width,height);
    glRotatef(alpha, 0, 1, 0);
    if(rotate) {
      alpha++;
      if(alpha > 360) alpha = 0;
      if(!keepdrawing) {
        SDL_Event redraw;
        redraw.type = SDL_VIDEOEXPOSE;
        SDL_PushEvent(&redraw);
      }
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
        
    if(extract) {
      extraction.frustum.GetView();
      extraction.metric->GetView();
      if(!realtime) {
	extraction.Extract(&nexus);
      } else {
	extraction.Update(&nexus);
      }
    }
   if(do_render)
    nexus.Render(extraction, contest, &stats);
   else
    stats.Init();

    /*    if(show_borders) {
	  for(unsigned int i = 0; i < cells.size(); i++) {
	  Border &border = nexus.GetBorder(cells[i]);
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

    tframe.push_front(watch.Time());
    if(tframe.size() > tlen) tframe.pop_back();

    terror.push_front(extraction.max_error);
    if(terror.size() > tlen) terror.pop_back();

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
      if(false){
      glColor4f(0.6f, 0.6f, 0.6f, 0.5f);

      glBegin(GL_LINE_STRIP);
      for(unsigned int i = 0; i < tframe.size() -1; i++) {
	double diff = (tframe[i] - tframe[i+1]);
	//glVertex2f(i/1024.0f,0);
	glVertex2f(i/1024.0f,2*diff);
      }
      glEnd();

      glColor4f(0.0f, 0.6f, 0.2f, 0.5f);

      glBegin(GL_LINE_STRIP);
      for(unsigned int i = 0; i < terror.size() -1; i++) {
//	glVertex2f(i/1024.0f,0);
	glVertex2f(i/1024.0f,terror[i]/300);
      }
      glEnd();
              }
      glColor3f(1.0f, 1.0f, 1.0f);

      sprintf(buffer, "Ram size : %.2f / %.2f Mb", 
	      nexus.ram_used * nexus.chunk_size/(float)(1<<20),
	      nexus.ram_max * nexus.chunk_size/(float)(1<<20));
      gl_print(0.03, 0.15, buffer);

      sprintf(buffer, "Extr size: %.2f / %.2f Mb",
      	      extraction.extr_used * nexus.chunk_size/(float)(1<<20), 
      	      extraction.extr_max * nexus.chunk_size/(float)(1<<20));
      gl_print(0.03, 0.12, buffer);

      sprintf(buffer, "Draw size: %.2f / %.2f Mb",
      	      extraction.draw_used * nexus.chunk_size/(float)(1<<20), 
      	      extraction.draw_max * nexus.chunk_size/(float)(1<<20));
      gl_print(0.03, 0.09, buffer);
      
      sprintf(buffer, "Disk size: %.2f  / %.2f Mb",
      	      extraction.disk_max * nexus.chunk_size/(float)(1<<20), 
      	      extraction.disk_used * nexus.chunk_size/(float)(1<<20));
      gl_print(0.03, 0.06, buffer);

      sprintf(buffer, "%.2f KTri  %.2f FPS %.0f M/s",
      	      stats.ktri/(float)(1<<10),
      	      stats.fps, stats.fps * stats.ktri/(float)(1<<20));
      gl_print(0.03, 0.03, buffer);

      
      glEnable(GL_DEPTH_TEST);
      glEnable(GL_LIGHTING);
      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
      glPopMatrix();
    
      if(output_stats){
      // statistics: output on file

      static Stats statsAcc;
      static      float ram_used ,float extr_used, float draw_used  ,float disk_used;
       static  std::ofstream outf(output_filename);
       static bool first=true;
       if(first) {
              outf<< "ktri\t fps\t ram \t extr \t draw \t disk \n"
              << "       \t        \t" << nexus.ram_max * nexus.chunk_size/(float)(1<<20) << "\t"
              << extraction.extr_max * nexus.chunk_size/(float)(1<<20) << "\t"
              <<  extraction.draw_max * nexus.chunk_size/(float)(1<<20)<<"\t"
              <<  extraction.disk_max * nexus.chunk_size/(float)(1<<20)<< "\n";
              first = false;
       }
       statsAcc.count++ ;
      if((statsAcc.count%30)==0) {
        outf
          <<       (statsAcc.ktri/(float)statsAcc.count)/(float)(1<<10)      << "\t" 
          <<       (statsAcc.fps/(float)statsAcc.count)           << "\t" 
         // <<       (statsAcc.kdisk/(float)statsAcc.count)       << "\t" 
          <<        ram_used  /(float)statsAcc.count  * nexus.chunk_size/(float)(1<<20) << "\t" 
          <<        extr_used/(float)statsAcc.count     * nexus.chunk_size/(float)(1<<20) << "\t" 
          <<        draw_used/(float)statsAcc.count   * nexus.chunk_size/(float)(1<<20) << "\t" 
          <<        disk_used/(float)statsAcc.count   * nexus.chunk_size/(float)(1<<20) << "\t" 
          <<        "\n";
        statsAcc.Init();
        statsAcc.count=0;
        statsAcc.fps=0;
        ram_used = extr_used= draw_used  = disk_used=0.0;
      }
      else{
        statsAcc.fps+=stats.fps;
        statsAcc.kdisk+=stats.kdisk;
        statsAcc.ktri+=stats.ktri;

        ram_used +=nexus.ram_used;
        extr_used+=extraction.extr_used;
        draw_used+=extraction.draw_used;
        disk_used+=extraction.disk_used;
      }
        
      }
    }
    
    SDL_GL_SwapBuffers();
  }

  SDL_Quit();
  return -1;
}

void gl_print(float x, float y, char *str) {
  glRasterPos2f(x, y);
  int len = (int)strlen(str);
  for(int i = 0; i < len; i++) 
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]);
}


