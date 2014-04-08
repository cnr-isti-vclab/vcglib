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

#ifndef __PICK______H
#define __PICK______H

#include <vector>
#include <algorithm>
namespace vcg{

template <class MESH_TYPE>
class GLPickTri
{
  typedef typename MESH_TYPE::FaceIterator FaceIterator;
  typedef typename MESH_TYPE::VertexIterator VertexIterator;
  typedef typename MESH_TYPE::FacePointer  FacePointer;
  typedef typename MESH_TYPE::VertexPointer  VertexPointer;
  typedef typename MESH_TYPE::VertexType  VertexType;

private:

  static Point3f Proj(const Eigen::Matrix4f &M, const float * viewport, const Point3f &p)
  {
    const float vx=viewport[0];
    const float vy=viewport[1];
    const float vw2=viewport[2]/2.0f;
    const float vh2=viewport[3]/2.0f;
    Eigen::Vector4f vp(p[0],p[1],p[2],1.0f);
    Eigen::Vector4f vpp = M*vp;
    Eigen::Vector4f ndc = vpp/vpp[3];

    Point3f sc(
        vw2*ndc[0] + vx+vw2,
        vh2*ndc[1] + vy+vh2,
        ndc[2]
        );

    return sc;
  }

  static void FillProjectedVector(MESH_TYPE &m, std::vector<Point3f> &pVec, const Eigen::Matrix4f &M, const float * viewportF)
  {
    pVec.resize(m.vert.size());
    for(size_t i=0;i<m.vert.size();++i) if(!m.vert[i].IsD())
    {
      pVec[i] = Proj(M, viewportF,m.vert[i].P());
    }
  }

  static void glGetMatrixAndViewport(Eigen::Matrix4f &M, float *viewportF)
  {
    Eigen::Matrix4d mp,mm;

    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    for(int i=0;i<4;++i) viewportF[i]=viewport[i];

    glGetDoublev(GL_PROJECTION_MATRIX, mp.data());
    glGetDoublev(GL_MODELVIEW_MATRIX,  mm.data());

    M = (mp*mm).cast<float>();
  }

  // compute a bbox in Device Coordinate (with the z without the near far normalization and ranged in -1 1)
  static Box3f ComputeDCBox(int x, int y, int width, int height)
  {
    Box3f bb;
    bb.SetNull();
    bb.Add(Point3f(x-width/2.0f,y-height/2.0f,-1.0f));
    bb.Add(Point3f(x+width/2.0f,y+height/2.0f, 1.0f));
    return bb;
  }

public:

  static bool PickClosestFace(int x, int y, MESH_TYPE &m, FacePointer &fp,int width=4, int height=4)
  {
    Eigen::Matrix4f M;
    float viewportF[4];
    glGetMatrixAndViewport(M,viewportF);
    Box3f reg=ComputeDCBox(x,y,width,height);

    float bzmin = std::numeric_limits<float>::max();
    fp=0;
    for(size_t i=0;i<m.face.size();++i) if(!m.face[i].IsD())
    {
      Point3f bary = vcg::Barycenter(m.face[i]);
      Point3f bz = Proj(M, viewportF,bary);

      if(bz[2]<bzmin && reg.IsIn(bz))
      {
        bzmin=bz[2];
        fp = &m.face[i];
      }
    }
    return fp!=0;
  }

  static bool PickClosestVert(int x, int y, MESH_TYPE &m, VertexPointer &vp,int width=4, int height=4)
  {
    Eigen::Matrix4f M;
    float viewportF[4];
    glGetMatrixAndViewport(M,viewportF);
    float bzmin = std::numeric_limits<float>::max();
    vp=0;

    Box3f reg=ComputeDCBox(x,y,width,height);

    for(size_t i=0;i<m.vert.size();++i) if(!m.vert[i].IsD())
    {
      Point3f bz = Proj(M, viewportF,m.vert[i].P());
      if(bz[2]<bzmin && reg.IsIn(bz))
      {
        bzmin=bz[2];
        vp = &m.vert[i];
      }
    }
    return vp!=0;
  }

  static int PickVert(int x, int y, MESH_TYPE &m, std::vector<VertexPointer> &result, int width=4, int height=4)
   {
     result.clear();
     static Eigen::Matrix4f lastM;
     static MESH_TYPE *lastm=0;
     static std::vector<Point3f> pVec;

     Eigen::Matrix4f M;
     float viewportF[4];
     glGetMatrixAndViewport(M,viewportF);

     Box3f reg =ComputeDCBox(x,y,width,height);

     if(M!=lastM || &m != lastm)
     {
       FillProjectedVector(m,pVec,M,viewportF);
       lastM = M;
       lastm = &m;
     }

     for(size_t i=0;i<m.vert.size();++i) if(!m.vert[i].IsD())
     {
       if(reg.IsIn(pVec[i]))
         result.push_back(&m.vert[i]);
     }
     return result.size();
   }

  static int PickFace(int x, int y, MESH_TYPE &m, std::vector<FacePointer> &result, int width=4, int height=4)
  {
    static Eigen::Matrix4f lastM;
    static MESH_TYPE *lastm=0;
    static std::vector<Point3f> pVec;

    float viewportF[4];
    Eigen::Matrix4f M;
    glGetMatrixAndViewport(M,viewportF);
    result.clear();
    Box3f reg;
    reg.Add(Point3f(x-width/2.0f,y-height/2.0f,-1.0f));
    reg.Add(Point3f(x+width/2.0f,y+height/2.0f,1.0f));

    if(M!=lastM || &m != lastm)
    {
      FillProjectedVector(m,pVec,M,viewportF);
      lastM = M;
      lastm = &m;
    }

    for(size_t i=0;i<m.face.size();++i) if(!m.face[i].IsD())
    {
      const Point3f &p0 = pVec[tri::Index(m,m.face[i].V(0))];
      const Point3f &p1 = pVec[tri::Index(m,m.face[i].V(1))];
      const Point3f &p2 = pVec[tri::Index(m,m.face[i].V(2))];
      if( (p0[2]>-1.0f) && (p1[2]>-1.0f) && (p2[2]>-1.0f) && IntersectionTriangleBox(reg,p0,p1,p2))
        result.push_back(&m.face[i]);
    }
    return result.size();
  }

  // Same of above but it also assumes that you want only visible faces.
  // Visibility is computed according to the current depth buffer.
  static int PickVisibleFace(int x, int y, MESH_TYPE &m, std::vector<FacePointer> &resultZ, int width=4, int height=4)
  {
    float vp[4];
    Eigen::Matrix4f M;
    glGetMatrixAndViewport(M,vp);

    int screenW = (int)(vp[2]-vp[0]);
    int screenH = (int)(vp[3]-vp[1]);

    GLfloat *buffer = new GLfloat[screenW*screenH];
    glReadPixels(vp[0],vp[1],vp[2],vp[3],GL_DEPTH_COMPONENT,GL_FLOAT,buffer);

    std::vector<FacePointer> result;
    PickFace(x,y,m,result,width,height);
    float LocalEpsilon = 0.001f;
    for(size_t i =0;i<result.size();++i)
    {
      Point3f p = Proj(M,vp,Barycenter(*(result[i])));
      if(p[0] >=0 && p[0]<screenW && p[1] >=0 && p[1]<screenH)
      {
        float bufZ = buffer[int(p[0])+int(p[1])*screenW];
        //qDebug("face %i txyz (%f %f %f)  bufz %f",i,tx,ty,tz,bufZ);
        if(bufZ + LocalEpsilon >= (p[2]+1.0)/2.0f)
          resultZ.push_back(result[i]);
      }
    }

    delete [] buffer;
    return resultZ.size();
  }



 // Same of above but it also assumes that you want only visible faces.
 // Visibility is computed according to the current depth buffer.
    static int OldPickFaceVisible(int x, int y, MESH_TYPE &m, std::vector<FacePointer> &resultZ, int width=4, int height=4, bool sorted=true)
    {
        // First step

        double mm[16];
        double mp[16];
        GLint vp[4];
        glGetIntegerv(GL_VIEWPORT,vp);
        glGetDoublev(GL_MODELVIEW_MATRIX ,mm);
        glGetDoublev(GL_PROJECTION_MATRIX ,mp);
        int screenW = 		vp[2]-vp[0];
        int screenH = 		vp[3]-vp[1];

        GLfloat *buffer = new GLfloat[screenW*screenH];
        glReadPixels(vp[0],vp[1],vp[2],vp[3],GL_DEPTH_COMPONENT,GL_FLOAT,buffer);

        std::vector<FacePointer> result;
        OldPickFace(x,y,m,result,width,height,sorted);
        float LocalEpsilon = 0.001f;
    for(size_t i =0;i<result.size();++i)
        {
            typename VertexType::CoordType v=Barycenter(*(result[i]));
            GLdouble tx,ty,tz;
            gluProject(v.X(),v.Y(),v.Z(), mm,mp,vp, &tx,&ty,&tz);
            if(tx >=0 && tx<screenW && ty >=0 && ty<screenH)
            {
                    float bufZ = buffer[int(tx)+int(ty)*screenW];
                    //qDebug("face %i txyz (%f %f %f)  bufz %f",i,tx,ty,tz,bufZ);
                    if(bufZ + LocalEpsilon >= tz)
                                resultZ.push_back(result[i]);
            }
        }

        delete [] buffer;
        return resultZ.size();
    }


    static int OldPickFace(int x, int y, MESH_TYPE &m, std::vector<FacePointer> &result, int width=4, int height=4,bool sorted=true)
      {
          result.clear();
      if(width==0 ||height==0) return 0;
          long hits;
          int sz=m.face.size()*5;
          GLuint *selectBuf =new GLuint[sz];
          //  static unsigned int selectBuf[16384];
          glSelectBuffer(sz, selectBuf);
          glRenderMode(GL_SELECT);
          glInitNames();

          /* Because LoadName() won't work with no names on the stack */
          glPushName(-1);
          double mp[16];

          GLint viewport[4];
          glGetIntegerv(GL_VIEWPORT,viewport);
          glMatrixMode(GL_PROJECTION);
          glGetDoublev(GL_PROJECTION_MATRIX ,mp);
          glPushMatrix();
          glLoadIdentity();
          //gluPickMatrix(x, viewport[3]-y, 4, 4, viewport);
          gluPickMatrix(x, y, width, height, viewport);
          glMultMatrixd(mp);

          glMatrixMode(GL_MODELVIEW);
          glPushMatrix();
          int fcnt=0;
          FaceIterator fi;
          for(fi=m.face.begin();fi!=m.face.end();++fi)
          {
              if(!(*fi).IsD())
              {
                  glLoadName(fcnt);
                  glBegin(GL_TRIANGLES);
                      glVertex( (*fi).V(0)->P() );
                      glVertex( (*fi).V(1)->P() );
                      glVertex( (*fi).V(2)->P() );
                  glEnd();
              }
        fcnt++; // the counter should advance even for deleted faces!
          }

          glPopMatrix();
          glMatrixMode(GL_PROJECTION);
          glPopMatrix();
          glMatrixMode(GL_MODELVIEW);
          hits = glRenderMode(GL_RENDER);
          //xstring buf;
          //if (hits <= 0)     return 0;
          std::vector< std::pair<double,unsigned int> > H;
          for(long ii=0;ii<hits;ii++){
              //TRACE("%ui %ui %ui %ui\n",selectBuf[ii*4],selectBuf[ii*4+1],selectBuf[ii*4+2],selectBuf[ii*4+3]);
              H.push_back( std::pair<double,unsigned int>(selectBuf[ii*4+1]/4294967295.0,selectBuf[ii*4+3]));
          }
          if(sorted)
        std::sort(H.begin(),H.end());
          //  if(H.size()>0) TRACE("\n Closest is %i\n",H[0].second);
          result.resize(H.size());
          for(long ii=0;ii<hits;ii++){
              FaceIterator fi=m.face.begin();
              advance(fi ,H[ii].second);
              result[ii]=&*fi;
          }

          delete [] selectBuf;
          return result.size();
      }

    static int OldPickVert(int x, int y, MESH_TYPE &m, std::vector<VertexPointer> &result, int width=4, int height=4,bool sorted=true)
    {
        result.clear();
        if(width==0 ||height==0) return 0;
        long hits;
        int sz=m.vert.size()*5;
        GLuint *selectBuf =new GLuint[sz];
        glSelectBuffer(sz, selectBuf);
        glRenderMode(GL_SELECT);
        glInitNames();

        /* Because LoadName() won't work with no names on the stack */
        glPushName(-1);
        double mp[16];

        GLint viewport[4];
        glGetIntegerv(GL_VIEWPORT,viewport);
        glMatrixMode(GL_PROJECTION);
        glGetDoublev(GL_PROJECTION_MATRIX ,mp);
        glPushMatrix();
        glLoadIdentity();
        gluPickMatrix(x, y, width, height, viewport);
        glMultMatrixd(mp);

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        int vcnt=0;
        VertexIterator vi;
        for(vi=m.vert.begin();vi!=m.vert.end();++vi)
        {
          if(!(*vi).IsD())
          {
            glLoadName(vcnt);
            glBegin(GL_POINTS);
              glVertex( (*vi).P() );
            glEnd();
          }
          vcnt++; // the counter should advance even for deleted faces!
        }

        glPopMatrix();
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        hits = glRenderMode(GL_RENDER);
        std::vector< std::pair<double,unsigned int> > H;
        for(long ii=0;ii<hits;ii++){
          H.push_back( std::pair<double,unsigned int>(selectBuf[ii*4+1]/4294967295.0,selectBuf[ii*4+3]));
        }
        if(sorted)
          std::sort(H.begin(),H.end());
        result.resize(H.size());
        for(long ii=0;ii<hits;ii++){
          VertexIterator vi=m.vert.begin();
          advance(vi ,H[ii].second);
          result[ii]=&*vi;
        }

        delete [] selectBuf;
        return result.size();
    }

};

}

#endif
