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

#ifndef __VCG_MESH_VISIBILITY
#define __VCG_MESH_VISIBILITY
#include <bitset>
#include <vcg/math/matrix44.h>
#include "simplepic.h"
namespace vcg {

template <class ScalarType>
void GenNormal(int vn, std::vector<Point3<ScalarType > > &NN)
{
  typedef Point3<ScalarType> Point3x;
    NN.clear();
    while(NN.size()<vn)
    {
      Point3x pp(float(rand())/RAND_MAX,float(rand())/RAND_MAX,float(rand())/RAND_MAX);
      pp=pp*2.0-Point3x(1,1,1);
      if(pp.SquaredNorm()<1)
      {
        Normalize(pp);
        NN.push_back(pp);
      }
    }
  }
  // Base Class che definisce le varie interfaccie;
template <class MESH_TYPE, int MAXVIS=2048> class VisShader 
{
	public :
		enum {VisMax=MAXVIS};
	VisShader(MESH_TYPE &me):m(me)
	{
			CullFlag= false;
			IsClosed = false;
			ZTWIST=1e-3;
      SplitNum=1;
			
			CameraViewing=false;
	}

		typedef Point3<typename MESH_TYPE::ScalarType> Point3x;

    typedef typename MESH_TYPE::CoordType CoordType;
    typedef typename MESH_TYPE::ScalarType ScalarType;
		typedef typename MESH_TYPE::VertexType  VertexType;
    typedef typename MESH_TYPE::VertexPointer  VertexPointer;
    typedef typename MESH_TYPE::VertexIterator VertexIterator;
    typedef typename MESH_TYPE::FaceIterator   FaceIterator;
    typedef typename MESH_TYPE::FaceType   FaceType;
	typedef Matrix44<ScalarType> Matrix44x;
	typedef Box3<ScalarType> Box3x;

// The Basic Data the mesh and its wrapper;
	MESH_TYPE &m;

  std::vector<MESH_TYPE *> OMV;  // Occluder Mesh Vector;

// la visibilita' e' in float, per ogni entita' 
// 1 significa che e' totalmente visibile per una data direzione.

	std::vector<float> VV;
	std::vector< Point3x > VN;				 // Vettore delle normali che ho usato per calcolare la mask e i float in W;

	 // User defined parameters and flags
	 bool IsClosed;
	 float ZTWIST;
	 bool CullFlag;   // Enable the frustum culling. Useful when the splitting value is larger than 2 
	 int SplitNum;
	 protected:
	 bool CameraViewing;
	 //Camera<ScalarType> Cam;
	 public: 

/********************************************************/
// Generic functions with Specialized code for every subclass
	 virtual void MapVisibility(float Gamma=1, float LowPass=0, float HighPass=1,bool FalseColor=false)=0;
   //virtual void ApplyLightingEnvironment(std::vector<float> &W, float Gamma);
	 
	 virtual int GLAccumPixel(	std::vector<int> &PixSeen)=0;
	 
	 virtual bool ReadVisibility(const char * /*filename*/){assert( 0); return false;}
	 virtual bool WriteVisibility(const char * /*filename*/){assert( 0); return false;}

/********************************************************/
// Generic functions with same code for every subclass

	void Clear() {		fill(VV.begin(),VV.end(),0); }

	void InitGL() 
		{
		  glPushAttrib(GL_COLOR_BUFFER_BIT );
      ::glClearColor (1.0, 1.0, 1.0, 0.0);
			glMatrixMode (GL_PROJECTION);   			
			glPushMatrix();
			glMatrixMode (GL_MODELVIEW);    
			glPushMatrix();
		}

	void RestoreGL()
		{
			glMatrixMode (GL_PROJECTION);   			
			glPopMatrix();
			glMatrixMode (GL_MODELVIEW);    
			glPopMatrix();
			glPopAttrib();
		}


/*
 Funzione principale di conversione in visibilita'
 Dati i due vettori PixSeen e PixNotSeen che indicano per ogni entita' (vertice o faccia) 
 quanti sono, rispettivamente,  i pixel visibili e occlusi,
 questa funzione calcola un valore float per ogni entita' che indica quanto  e' visibile lungo una data direzione camera 
 == 1 significa completamente visibile
 == 0 significa completamente occluso.

*/
	void AddPixelCount(std::vector<float> &_VV, const std::vector<int> &PixSeen)
		{
		assert(_VV.size()==PixSeen.size());
			for(int i=0;i<PixSeen.size();++i)
				if(PixSeen[i]>0) _VV[i]+= 1;
		}
 

 //void SetVisibilityMask(std::vector< std::bitset<MAXVIS> > &_VM, const std::vector<int> &PixSeen, const int dir)
	//	{
	//	assert(_VM.size()==PixSeen.size());
	//		for(int i=0;i<PixSeen.size();++i)
	//			if(PixSeen[i]>0) _VM[i][dir]=true;
	//	}

/*******************************
Funzioni ad alto livello che computano le Visibility Mask per varie distribuzioni di direzioni


*******************************/

// Funzione Generica chiamata da tutte le seguenti 
void Compute( CallBack *cb)
{
  //cb(buf.format("Start to compute %i dir\n",VN.size()));
	InitGL();
	VV.resize(m.vert.size());
	vector<int> PixSeen(VV.size(),0);
	
	for(int i=0;i<VN.size();++i)
		{
      int t0=clock(); 
			fill(PixSeen.begin(),PixSeen.end(),0);
      int added=SplittedRendering(VN[i], PixSeen,cb);	
		  AddPixelCount(VV,PixSeen);
      int t1=clock();
      printf("ComputeSingleDir %i on %i : %i msec\n",i,VN.size(),t1-t0);
		}
	
  RestoreGL();
}

void ComputeCone(int nn, Point3x &dir, ScalarType ConeAngleDeg, CallBack *cb)
{
	string buf;
	ScalarType ConeAngleRad=ToRad(ConeAngleDeg);
	ScalarType SolidAngleSter = (1.0 - Cos(ConeAngleRad))*2*M_PI;
	int Frac=(4.0*M_PI)/SolidAngleSter;
	printf("ComputeCone for an angle of %f , solidAngle =%f pi asked %i normals, forecasted we need to ask %i normals\n",
				ConeAngleDeg,SolidAngleSter/M_PI,nn,nn*Frac);

	VN.clear();
	vector<Point3x> nvt;
	GenNormal(nn*Frac,nvt);
	ScalarType CosConeAngle=Cos(ConeAngleRad);
	for(int i=0;i<nvt.size();++i)
		if(dir*nvt[i]>=CosConeAngle) VN.push_back(nvt[i]);
 
	printf("Asked %i normal, got %i normals\n",nn,VN.size());
  Compute(cb); 
}

void ComputeHalf(int nn, Point3x &dir, CallBack *cb)
{
	string buf;
	
	VN.clear();
	vector<Point3x> nvt;
	GenNormal(nn*2,nvt);
	for(int i=0;i<nvt.size();++i)
		if(dir*nvt[i]>0) VN.push_back(nvt[i]);
 
	printf("Asked %i normal, got %i normals\n",nn,VN.size());
  Compute(cb); 
}

void ComputeUniform(int nn, CallBack *cb)
{
	VN.clear();
	GenNormal(nn,VN);
	char buf[256];
	sprintf(buf,"Asked %i normal, got %i normals\n",nn,VN.size());
  cb(buf);
  Compute(cb); 
}

void ComputeSingle(Point3x &dir, CallBack *cb)
{
	VN.clear();
	VN.push_back(dir);
	printf("Computing one direction (%f %f %f)\n",dir[0],dir[1],dir[2]);
  Compute(cb); 
}

/**********************************************************/

int SplittedRendering(Point3x &ViewDir, std::vector<int> &PixSeen, CallBack *cb=DummyCallBack)
{
  int tt=0;
  int i,j;
  for(i=0;i<SplitNum;++i)
    for(j=0;j<SplitNum;++j){
        SetupOrthoViewMatrix(ViewDir, i,j,SplitNum);
 	      tt+=GLAccumPixel(PixSeen);
    }
    return tt;
}

// Compute a rotation matrix that bring Axis parallel to Z. 
void GenMatrix(Matrix44d &a, Point3d Axis, double angle)
{
	const double eps=1e-5;
	Point3d RotAx   = Axis ^ Point3d(0,0,1);
  double RotAngle = Angle(Axis,Point3d(0,0,1));

  if(math::Abs(RotAx.Norm())<eps) { // in questo caso Axis e' collineare con l'asse z
			RotAx=Point3d(0,1,0);
      double RotAngle = Angle(Axis,Point3d(0,1,0));
		}
  //printf("Rotating around (%5.3f %5.3f %5.3f) %5.3f\n",RotAx[0],RotAx[1],RotAx[2],RotAngle);
  a.SetRotate(-RotAngle,RotAx);
	//Matrix44d rr;
  //rr.SetRotate(-angle, Point3d(0,0,1));
	//a=rr*a;
}


// Genera la matrice di proj e model nel caso di un rendering ortogonale.
// subx e suby indicano la sottoparte che si vuole
void SetupOrthoViewMatrix(Point3x &ViewDir, int subx, int suby,int LocSplit)
{
	glMatrixMode (GL_PROJECTION);   			
	glLoadIdentity (); 
  float dlt=2.0f/LocSplit;

  glOrtho(-1+subx*dlt, -1+(subx+1)*dlt, -1+suby*dlt, -1+(suby+1)*dlt,-2,2);
	glMatrixMode (GL_MODELVIEW);    
	glLoadIdentity ();  
	Matrix44d rot;
	Point3d qq; qq.Import(ViewDir);
	GenMatrix(rot,qq,0);
	glMultMatrixd((const double *)&rot);
	double d=2.0/m.bbox.Diag();
  glScalef(d,d,d);
	glTranslate(-m.bbox.Center());
}

void ComputeSingleDirection(Point3x BaseDir, std::vector<int> &PixSeen, CallBack *cb=DummyCallBack)
{
	int t0=clock();
	string buf;
 
	int added=SplittedRendering(BaseDir, PixSeen,cb);	
  int t1=clock();
	printf("ComputeSingleDir %i msec\n",t1-t0);
}

void ComputeAverageVisibilityDirection()
{
	int i,j;
	VD.resize(VM.size());
	for(j=0;j<VM.size();++j)
		{
			Point3x &nn=VD[j];
			nn=Point3x(0,0,0);
			bitset<VisMax> &msk=VM[j];
				for(i=0;i<VN.size();++i)
				    if(msk[i]) nn+=VN[i];
		}
		for(j=0;j<VM.size();++j)
			 VD[j].Normalize();
		
}

// calcola un LightingEnvironment direzionale, cioe'un vettore di pesi per l'insieme di normali 
// corrente tale che 
// mette a 1 tutti i vettori che sono entro un angolo DegAngle1 
// a 0 tutti quelli oltre DegAngle2 e
// sfuma linearmente nel mezzo. 
void DirectionalLightingEnvironment(std::vector<float> &LE, Point3x dir, ScalarType DegAngle1, ScalarType DegAngle2)
{
	LE.clear();
	LE.resize(VN.size(),0);
	int i;
	for(i=0;i<VN.size();++i)
		{
			ScalarType a=ToDeg(Angle(dir,VN[i]));
			if(a<DegAngle1) { LE[i]=1; continue; }
			if(a>DegAngle2) { LE[i]=0; continue; }
			LE[i] = 1.0-(a-DegAngle1)/(DegAngle2-DegAngle1);

		}
	// last step normalize the weights;
	ScalarType sum=0;
	for(i=0;i<VN.size();++i)
		sum+=LE[i];
	for(i=0;i<VN.size();++i)
		LE[i]/=sum;
}


};
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

template <class MESH_TYPE> class VertexVisShader : public VisShader<MESH_TYPE>
{
	public :

	// Function Members
	VertexVisShader(MESH_TYPE &me):VisShader<MESH_TYPE>(me)
	{
     // la mesh DEVE avere colore per vertice
			if(! m.HasPerVertexColor()) assert(0);
	}

	void Init()  {		VV.resize(m.vert.size()); AssignColorId(); }
//	Vis::VisMode Id() {return Vis::VMPerVert;};
	void Compute(int nn);
	
	void AddPixelCount(std::vector<float> &_VV, std::vector<int> &PixSeen, std::vector<int> &PixNotSeen )
		{
			for(int i=0;i<m.vert.size();++i)
				_VV[i]+= float(PixSeen[i])/float(PixSeen[i]-PixNotSeen[i]);
		}

void DrawFill(MESH_TYPE &mm)
{
  glBegin(GL_TRIANGLES);
  MESH_TYPE::FaceIterator fi;
  for(fi=mm.face.begin();fi!=mm.face.end();++fi)
  {
    glVertex((*fi).V(0)->P());
    glVertex((*fi).V(1)->P());
    glVertex((*fi).V(2)->P());
  }
  glEnd();
}

/***************************************************************************/

//VertexVisibility
// Funzione Principale restituisce per ogni entita' quanti px si vedono o no.
int GLAccumPixel(	std::vector<int> &PixSeen)
{
	SimplePic<float> snapZ;
	SimplePic<Color4b> snapC;

  glClearColor(Color4b::White);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_POLYGON_BIT );
	glDisable(GL_LIGHTING);	
	glDepthRange(0.0f,1.0f);
	glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
	glDepthMask(GL_TRUE);
	glDrawBuffer(GL_BACK);
	glReadBuffer(GL_BACK);
	
  /////** Si disegnano le front face  **/////
  glDepthRange(2.0*ZTWIST,1.0f);
//	glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glColor(Color4b::Red);
	DrawFill(m);
  snapC.OpenGLSnap();
	
 //	if(!IsClosed) // sono necessarie due passate in piu!
	//{
	//	/////** Si disegnano le back face shiftate in avanti di 2 epsilon **/////
	//	glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
	//	glColor(Color4b::Black);
	//	glCullFace(GL_FRONT);
	//	glDepthRange(0.0f,1.0f-2.0*ZTWIST);
	//  DrawFill(m);
	//	//glw.DrawFill<GLW::NMNone, GLW::CMNone,GLW::TMNone>();  
	//}
	//if(OMV.size()>0)
	//{
	//	/////** Si disegnano le back face shiftate in avanti di 2 epsilon **/////
	//	glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
	//	glColor(Color4b::Black);
	//	glDisable(GL_CULL_FACE);
	//	glDepthRange(0.0f,1.0f-2.0*ZTWIST);
	//	for(unsigned int i=0;i<OMV.size();++i)
	//    DrawFill(*(OMV[i]));
	//		//OGV[i]->DrawFill<GLW::NMNone, GLW::CMNone,GLW::TMNone>();  
	//}
	int cnt=0;
	snapZ.OpenGLSnap(GL_DEPTH_COMPONENT);
	snapC.OpenGLSnap();
	double MM[16];
	glGetDoublev(GL_MODELVIEW_MATRIX,MM);
	double MP[16];
  glGetDoublev(GL_PROJECTION_MATRIX,MP);
	int VP[4];
	glGetIntegerv(GL_VIEWPORT,VP);
	double tx,ty,tz;
	for(int i=0;i<m.vert.size();++i)
	{
		gluProject(m.vert[i].P()[0],m.vert[i].P()[1],m.vert[i].P()[2],
			MM,MP,VP,
			&tx,&ty,&tz);
    if(tx>=0 && tx<snapZ.sx && ty>=0 && ty<snapZ.sy)
    {
		    int txi=floor(tx),tyi=floor(ty);
		    float sd=snapZ.Pix(tx,ty);
    		
		    int col = min( min(snapC.Pix(txi+0,tyi+0)[0],snapC.Pix(txi+1,tyi+0)[0]),
						          min(snapC.Pix(txi+0,tyi+1)[0],snapC.Pix(txi+1,tyi+1)[0]));
		    if(col!=0)
		    if(tz<sd+ZTWIST) {
			    PixSeen[i]++;
			    cnt++;
		    }
	    }
  }
  	glPopAttrib();
printf("Seen %i vertexes on %i\n",cnt,m.vert.size());
return cnt;
}

void SmoothVisibility()
{
	MESH_TYPE::FaceIterator fi;
	vector<float> VV2;
	vector<int> VC(VV.size(),1);
	VV2=VV;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
		for(int i=0;i<3;++i)
		{
			VV2[(*fi).V(i)-&*m.vert.begin()] += VV[(*fi).V1(i)-&*m.vert.begin()];
			++VC[(*fi).V(i)-&*m.vert.begin()];
		}

	for(unsigned int i=0;i<VV2.size();++i)
		VV[i]=VV2[i]/VC[i];
}

void MapVisibility(float Gamma=1, float LowPass=0, float HighPass=1, bool FalseColor = false)
{
	float minv=*min_element(VV.begin(),VV.end());
	float maxv=*max_element(VV.begin(),VV.end());
	printf("Visibility Range %f %f\n", minv,maxv);

			MESH_TYPE::VertexIterator vi;
			for(vi=m.vert.begin();vi!=m.vert.end();++vi){
				float gval=(VV[vi-m.vert.begin()]-minv)/(maxv-minv);
				if(gval<LowPass) gval=LowPass;
				if(gval>HighPass) gval=HighPass;
				(*vi).C().SetGrayShade(pow((gval-LowPass)/(HighPass-LowPass),Gamma));
			}
}

//void ApplyLightingEnvironment(std::vector<float> &W, float Gamma=1)
//	{
//		assert(W.size()==VN.size());
//		MESH_TYPE::VertexIterator vi;
//	
//		for(vi=m.vert.begin();vi!=m.vert.end();++vi)
//		{
//		float gray=0;
//		bitset<VisMax> &msk=VM[vi-m.vert.begin()];
//			for(int i=0;i<VN.size();++i)
//				if(msk[i]) gray+=W[i];
//			
//			(*vi).C().SetGrayShade(gray);
//		}
//	}

};



}
#endif // __VCG_MESH_VISIBILITY