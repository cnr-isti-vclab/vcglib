#ifndef SEGMENTATOR
#define SEGMENTATOR

//#include <vcg/simplex/vertex/with/afvn.h>
//#include <vcg/simplex/face/with/afav.h>

#include <vcg/simplex/vertex/with/vn.h>
#include <vcg/simplex/face/with/af.h>

#include <sim/particle/with/basic_physics.h>
#include <sim/methods/mass_spring/triangle.h>


#include <vcg/complex/trimesh/base.h>
#include <vcg/complex/trimesh/allocate.h>
#include <vcg/complex/trimesh/update/topology.h>
#include <vcg/complex/trimesh/update/normal.h>

#include <vcg/complex/trimesh/refine.h>
#include <vcg/complex/trimesh/platonic.h>
#include <volume_dataset.h>

//#include <vcg/simplex/face/pos.h>


#include <vcg/space/point3.h>
#include <vcg/space/box3.h>

#include <sim/pde_integrator.h>
#include <partial_container.h>
#include <vector>
#include <time.h>
#include <math.h>
#include <collision_detection.h>
#include <vcg/complex/trimesh/smooth.h>

class Segmentator{
	
public:

struct DummyEdge;
struct DummyTetra;
struct MyFace;

struct MyVertex: public ParticleBasic<vcg::VertexVNf<DummyEdge,MyFace,DummyTetra> >
{
public:

	bool blocked;//optimize after with vertex flags
	bool stopped;

	MyVertex()
	{
		blocked=false;
		stopped=false;
		Acc()=Point3f(0,0,0);
		Vel()=Point3f(0,0,0);
		//neeed call of the super class
	}

	void UpdateAcceleration()
	{
		//if ((!IsBlocked(this))&&(!IsStopped(this)))
		if ((!blocked)&&(!stopped))
		{
			Acc()=(IntForce()+ExtForce())/Mass();
		}
		else
		{
			Acc()=Point3f(0,0,0);
			Vel()=Point3f(0,0,0);
		}
	}
	
	void Reset()
	{
		IntForce()=Point3f(0.f,0.f,0.f);
	}



};

///this class implements the deformable triangle in a mass spring system
struct MyFace : public TriangleMassSpring< vcg::FaceAF<MyVertex,DummyEdge,MyFace> >
{
public:
	bool intersected;
	float kdihedral;

	MyFace()
	{
		intersected=false;
	}

	void Init ( double k, double  mass,float k_dihedral )
	{ 
		__super::Init(k,mass);
		kdihedral=k_dihedral;
	}

	bool IsActive()
	{
		return(!(((V(0)->blocked)||(V(0)->stopped))&&
			  ((V(1)->blocked)||(V(1)->stopped))&&
			  ((V(2)->blocked)||(V(2)->stopped))));
	}

	bool IsBlocked()
	{
	return((V(0)->blocked)&&(V(1)->blocked)&&(V(2)->blocked));
	}
	
		
	double DiedralAngle(int edge)
	{
		MyFace *fopp=FFp(edge);
		CoordType norm1=NormalizedNormal();
		CoordType norm2=fopp->NormalizedNormal();
		return (NormalizedNormal()*fopp->NormalizedNormal());
	}

	///update of the internal forces using the dihedral angle
	bool Update ( void )	
	{
		for (int i=0;i<3;i++)
		{
			MyFace *fopp=FFp(i);
			MyFace *myAddr=fopp->FFp(FFi(i));
			if ((fopp!=0)||(fopp<myAddr))//test do not duplicate updates per edge
			{
				//normal and area based diadedral angle calcolus
				CoordType DirEdge=(V(i)->P()-V((i+1)%3)->P()).Normalize();
				fopp=FFp(i);
				CoordType Ver=(NormalizedNormal()^fopp->NormalizedNormal()).Normalize();
				ScalarType diaedral=DiedralAngle(i);

				if ((Ver*DirEdge)<=0)///convex
				{
					ScalarType Force=(((-diaedral)+1.f)*kdihedral);
					V((i+2)%3)->IntForce()+=NormalizedNormal()*(Force);
					V(i)->IntForce()-=NormalizedNormal()*(Force)/2.f;
					V((i+1)%3)->IntForce()-=NormalizedNormal()*(Force)/2.f;
				}
				else	///non-convex
				{
					ScalarType Force=(((-diaedral)+1.f)*kdihedral);
					V((i+2)%3)->IntForce()-=NormalizedNormal()*(Force);
					V(i)->IntForce()+=NormalizedNormal()*(Force)/2.f;
					V((i+1)%3)->IntForce()+=NormalizedNormal()*(Force)/2.f;
				}
			}
		}
		return(__super::Update());
	}
};

struct MyTriMesh: public vcg::tri::TriMesh<std::vector<MyVertex>,std::vector<MyFace> >{};

typedef  Partial_Container<std::vector<MyVertex*>,MyVertex> Part_VertexContainer;
typedef  Partial_Container<std::vector<MyFace*>,MyFace> Part_FaceContainer;
typedef  PDEIntegrator<Part_FaceContainer,Part_VertexContainer,float> myIntegrator;
typedef	 Collision_Detector<std::vector<MyFace> > Collision;


public:
	Point3f scale;
	
	//VolumetricDataset<int> d;
	MyTriMesh m;
	Part_FaceContainer P_Faces;
	Part_VertexContainer P_Vertex;
	Part_VertexContainer V_Stopped;
	myIntegrator *TrINT;
	
	MyTriMesh::CoordType InitialBarycenter;
	
	float mass;
	float k_elanst;
	float k_dihedral;
	float edge_size;
	int tolerance;
	int gray_init;
	float time_stamp;
	float edge_precision;

	bool end_loop;
	bool refined;

	clock_t interval_reinit;
	clock_t interval_selfcollision;

	Collision *CollDet;

	//Volume_Dataset_Optimized<short> V;
	Volume_Dataset <short> V;

	vcg::Box3<float>  bbox;

	char *inDir;
	char *outDir;

	//attention static members
	/*int BlockFlag;
	int StoppedFlag;*/

	Segmentator()
	{
		CollDet=new Collision(m.face);
	}

	~Segmentator()
	{
	}

private:

	/////return integer coordinete in volumetric dataset
	//Point3i MapToDataset(MyTriMesh::CoordType p)
	//{
	//	MyTriMesh::ScalarType x=((MyTriMesh::ScalarType)p.V(0));
	//	MyTriMesh::ScalarType y=((MyTriMesh::ScalarType)p.V(1));
	//	MyTriMesh::ScalarType z=((MyTriMesh::ScalarType)p.V(2));
	//	return Point3i((int)p.V(0),(int)p.V(1),(int)p.V(2));
	//}
	//

	///map to space coordinate from dataset coordinates
	MyTriMesh::CoordType MapToSpace(Point3i p)
	{
		MyTriMesh::ScalarType x=((MyTriMesh::ScalarType)p.V(0));
		MyTriMesh::ScalarType y=((MyTriMesh::ScalarType)p.V(1));
		MyTriMesh::ScalarType z=((MyTriMesh::ScalarType)p.V(2));
		return (MyTriMesh::CoordType(x,y,z));
	}

	///return integer coordinete in volumetric dataset
	float getColor(MyTriMesh::CoordType p)
	{
		

		float lx=(p.V(0)-(int)p.V(0))*scale.V(0);//da rivedere bene per lo scale
		float ly=(p.V(1)-(int)p.V(1))*scale.V(1);//da rivedere bene per lo scale
		float lz=(p.V(2)-(int)p.V(2))*scale.V(2);//da rivedere bene per lo scale

		p=Scale(p);

		Point3i base=Point3i((int)p.V(0),(int)p.V(1),(int)p.V(2));
				
		float v[8];
		Point3i px;
		for (int i=0;i<8;i++)
		{
			px=base+Point3i((i%2),(i/4),((i/2)%2));
			v[i]=(float)V.getAt(px);
		}

		float color=lx*v[1]+v[0]+lz*v[0]*lx-v[0]*lx-ly*v[0]+ly*v[2]-lz*v[0]+ly*v[0]*lx-ly*lx*v[1]-ly*v[2]*lx
		-lz*ly*v[0]*lx+lz*ly*lx*v[1]+lz*ly*v[2]*lx-lz*ly*lx*v[3]+lz*ly*lx*v[4]-lz*ly*lx*v[5]-lz*ly*
		v[6]*lx+lz*ly*lx*v[7]+ly*lx*v[3]-lz*lx*v[1]+lz*ly*v[0]-lz*ly*v[2]-lz*lx*v[4]+lz*lx*v[5]-lz*ly*
		v[4]+lz*ly*v[6]+lz*v[4];

		return color;
	}
	
	///maximixe the gradient of the movement
	MyTriMesh::CoordType Gradient(MyTriMesh::CoordType p,float h=0.01f)
	{
		float value=getColor(p);
		MyTriMesh::CoordType h0=MyTriMesh::CoordType(h,0,0);
		MyTriMesh::CoordType h1=MyTriMesh::CoordType(0,h,0);
		MyTriMesh::CoordType h2=MyTriMesh::CoordType(0,0,h);
		float dx=(getColor(p+h0)-value)/h;

		/*dx=v[1]+lz*v[0]-v[0]*lx+ly*v[0]-ly*v[1]-ly*v[2]
		-lz*ly*v[0]+lz*ly*v[1]+lz*ly*v[2]-lz*ly*v[3]+lz*ly*v[4]-lz*ly*v[5]-lz*ly*
		v[6]+lz*ly*v[7]+ly*v[3]-lz*v[1]-lz*v[4]+lz*v[5]-lz*ly*v[4]+lz*ly*v[6]+lz*v[4];
		
		dy=-v[0]+v[2]+v[0]*lx-lx*v[1]-v[2]*lx
		-lz*v[0]*lx+lz*lx*v[1]+lz*v[2]*lx-lz*lx*v[3]+lz*lx*v[4]-lz*lx*v[5]-lz*v[6]*lx+lz*ly*lx*v[7]+ly*lx*v[3]-lz*lx*v[1]+lz*ly*v[0]-lz*ly*v[2]-lz*lx*v[4]+lz*lx*v[5]-lz*ly*
		v[4]+lz*ly*v[6]+lz*v[4];*/

		float dy=(getColor(p+h1)-value)/h;
		float dz=(getColor(p+h2)-value)/h;
		MyTriMesh::CoordType ret=MyTriMesh::CoordType(dx,dy,dz);
		return (ret);
	}

	///scale the coordinates of a point
	MyTriMesh::CoordType UnScale(MyTriMesh::CoordType p)
	{
		MyTriMesh::ScalarType x=(p.V(0))*scale.V(0);
		MyTriMesh::ScalarType y=(p.V(1))*scale.V(1);
		MyTriMesh::ScalarType z=(p.V(2))*scale.V(2);
		return (MyTriMesh::CoordType(x,y,z));
	}

	///scale the coordinates of a point
	MyTriMesh::CoordType Scale(MyTriMesh::CoordType p)
	{
		MyTriMesh::ScalarType x=(p.V(0))/scale.V(0);
		MyTriMesh::ScalarType y=(p.V(1))/scale.V(1);
		MyTriMesh::ScalarType z=(p.V(2))/scale.V(2);
		return (MyTriMesh::CoordType(x,y,z));
	}

	///return true if a coordinate is out of limits
	bool OutOfLimits(MyTriMesh::CoordType p)
	{
		/*Point3f max=Scale(MapToSpace(V.Max()));
		Point3f min=Scale(MapToSpace(V.Min()));
		Point3f test=(Scale(p));*/
		Point3f test=p;
		Point3f max=UnScale(MapToSpace(V.Max()));
		Point3f min=UnScale(MapToSpace(V.Min()));
		
		for (int i=0;i<3;i++)
		{
			if(((test.V(i)>=max.V(i))||(test.V(i)<=min.V(i))))
				return true;
		}
		return false;
	}
	
	
	bool IsBlocked(MyVertex *v)
	{
		//return ((v->Flags()& BlockFlag!=0));
		return (v->blocked);
	}

	void SetBlocked(MyVertex *v)
	{
		//v->Flags()|= BlockFlag;
		v->blocked=true;
		//v->SetS();//for refine
	}

	void SetBlockedFace(MyFace *f)
	{
		SetBlocked(f->V(0));
		SetBlocked(f->V(1));
		SetBlocked(f->V(2));
	}

	void SetIntersectedFace(MyFace *f)
	{
		f->intersected=true;
	}

	bool IsStopped(MyVertex *v)
	{
		//return ((v->Flags()& StoppedFlag!=0));
		return (v->stopped);
	}

	void SetStopped(MyVertex *v)
	{
		//v->Flags()|= StoppedFlag;
		v->stopped=true;
		V_Stopped.push_back(v);
	}

	void ClearStopped(MyVertex *v)
	{
		//v->Flags()&= ~StoppedFlag;
		v->stopped=false;
	}

///re-set physical pararmeters on the mesh
void InitPhysParam(float k_elanst,float mass,float k_dihedral)
{
	for (unsigned int i=0;i<m.face.size();i++)
		{
			m.face[i].Init(k_elanst,mass,k_dihedral);
		}
}

///set the initial mesh of deformable object
void InitMesh(MyTriMesh &m)
	{
		m.Clear();
		
		vcg::tri::Icosahedron<MyTriMesh>(m);

		vcg::tri::UpdateTopology<MyTriMesh>::FaceFace(m);
		
	/*	P_Vertex.clear();
		P_Faces.clear();*/
		
		for (unsigned int i=0;i<m.vert.size();i++)
		{
			m.vert[i].P()+=InitialBarycenter;
			m.vert[i].P()=UnScale(m.vert[i].P());
		//	P_Vertex.push_back(&m.vert[i]);
		}
		
		
		vcg::tri::UpdateNormals<MyTriMesh>::PerVertexNormalized(m);

	}


///return true if the gray level of the vertex v differ from graylevel less than tolerance
bool InTolerance(MyTriMesh::VertexType *v)
{
	return (abs(getColor(v->P())-gray_init)<tolerance);
}

///add to the vertex v a containing force basing on diffence from tolerance 
MyTriMesh::CoordType ContainingForce(MyTriMesh::VertexType *v)
{
	//float dinstance=fabs((FindGrayMedia(v))-gray_init);
	float dinstance=fabs((getColor(v->P()))-(float)gray_init);
	assert(dinstance<=tolerance);
	MyTriMesh::CoordType ret=(-v->N()*((dinstance)/(float)tolerance));
	return (ret);
}

///find the gradient factor
MyTriMesh::CoordType GradientFactor(MyTriMesh::VertexType *v)
{
	MyTriMesh::CoordType value=Gradient(v->P());
	/*float d0=getColor(v->P()+value);
	float d1=getColor(v->P()-value);
	if ((fabs(d0-(float)gray_init))>(fabs(d1-(float)gray_init)))
		return (-value);
	else */
	return (value*(gray_init-getColor(v->P())));
}

///add the external forces to the deformable mesh
void AddExtForces()
{
	Part_VertexContainer::iterator vi;

	end_loop=true;
	for (vi=P_Vertex.begin();vi<P_Vertex.end();++vi)
	{
		
		if (!(*vi).IsD())
		{
			if (OutOfLimits((*vi).P()))
				//vi->blocked=true;
				SetBlocked(&*vi);
			if ((!IsBlocked(&*vi))&&(!IsStopped(&*vi)))
			{
				end_loop=false;
				if (!InTolerance(&*vi))
				{
					SetBlocked(&*vi);
					(*vi).ExtForce()=MyTriMesh::CoordType(0,0,0);
				}
				else
				{
					MyTriMesh::CoordType Inflating=(*vi).N();
					MyTriMesh::CoordType Containing0=ContainingForce(&*vi);	
					//MyTriMesh::CoordType Containing1=GradientFactor(&*vi);
					Containing0*=0.75;
					if (Containing0.Norm()>1)
						Containing0.Normalize();
					(*vi).ExtForce()=Inflating+Containing0;/*+Containing1+Containing0*/;
				}
			}
			else
				(*vi).ExtForce()=MyTriMesh::CoordType(0,0,0);
		}
	}
}

///reinit the partial integration vectors that describe active vertices
void Reinit_PVectors()
{
	V_Stopped.clear();
    P_Vertex.clear();
	MyTriMesh::VertexIterator vi;

	for (vi=m.vert.begin();vi<m.vert.end();vi++)
	{
		if ((!vi->IsD())&&(!vi->blocked))
			P_Vertex.push_back(&(*vi));	
		if ((!vi->IsD())&&((*vi).stopped)&&(!vi->blocked))
			V_Stopped.push_back(&(*vi));		
	}

	P_Faces.clear();
	MyTriMesh::FaceIterator fi;
	for (fi=m.face.begin();fi<m.face.end();fi++)
	{
		//if ((!fi->IsBlocked()))
		if ((!fi->IsD())&&(!fi->IsBlocked()))
			P_Faces.push_back(&(*fi));
	}
}

///erase the stopped entities from the partial containers
void Refresh_PVectors()
{
	Part_FaceContainer P_FacesAux;
	Part_VertexContainer P_VertexAux;
	P_FacesAux.clear();
	P_VertexAux.clear();
	
	int i=0;
	for (i=0;i<P_Vertex.size();i++)
	{
		if (!P_Vertex[i]->blocked)
			P_VertexAux.push_back(P_Vertex[i]);
	}

	for (i=0;i<P_Faces.size();i++)
	{
		if (!P_Faces[i]->IsBlocked())
			P_FacesAux.push_back(P_Faces[i]);
	}

	P_Faces.clear();
	P_Vertex.clear();

	P_Faces=P_FacesAux;
	P_Vertex=P_VertexAux;
}

///add the new elements on partial vectors when allocate space for new vertices
void AddNewElements(MyTriMesh::VertexIterator vi,MyTriMesh::FaceIterator fi)
{
	while (vi!=m.vert.end())
	{
		if (!(*vi).IsD())
			P_Vertex.push_back(&(*vi));
		vi++;
	}
	while (fi!=m.face.end())
	{
		if (!(*fi).IsD())
			P_Faces.push_back(&(*fi));
		fi++;
	}
}

///verify and eventually stop the vertices of the mesh
void VerifyForces()
{
	float proj;
	Part_VertexContainer::iterator vi;
	for (vi=P_Vertex.begin();vi<P_Vertex.end();++vi)
	{
	 if (!IsStopped(&*vi))
	 {
		MyTriMesh::CoordType accn=(*vi).Acc();
		proj=accn*(*vi).N();
		if  ((proj)<=0)
			SetStopped(&*vi);
	 }
	}
}

bool TimeReinit()
{
	static clock_t time=0;
	clock_t elapsedsecs=abs(time-clock());
	if (elapsedsecs>interval_reinit)
	{
		time=clock();
		return true;
	}
	return false;
}

bool TimeSelfIntersection()
{
	static clock_t time=0;
	clock_t elapsedsecs=abs(time-clock());
	if (elapsedsecs>interval_selfcollision)
	{
		time=clock();
		return true;
	}
	return false;
}

///refine the mesh and re-update eventually 
void RefineStep(float _edge_size)
{
	MyTriMesh::VertexIterator vinit=m.vert.begin();
	MyTriMesh::FaceIterator finit=m.face.begin();
	MyTriMesh::VertexIterator vend=m.vert.end();
	MyTriMesh::FaceIterator fend=m.face.end();

	refined=vcg::Refine(m,MidPoint<MyTriMesh>(),_edge_size);
	
	if (refined)
	{
		MyTriMesh::VertexIterator vinit2=m.vert.begin();
		MyTriMesh::FaceIterator finit2=m.face.begin();

		if ((vinit2!=vinit)||(finit2!=finit))
			Reinit_PVectors();
		else
			AddNewElements(vend,fend);

		vcg::tri::UpdateNormals<MyTriMesh>::PerVertexNormalized(m);
		CollDet->RefreshElements();
	}
}

///reset vertex position and unblock them
void ReinitPhysicMesh()
{
	Part_FaceContainer::iterator pfi;

	for (pfi=P_Faces.begin();pfi<P_Faces.end();++pfi)
		(*pfi).Init(k_elanst,mass,k_dihedral);

	/*for (MyTriMesh::VertexIterator vi=m.vert.begin();vi<m.vert.end();vi++)
		ClearStopped(&*vi);*/
}

///clear the stopped vertex
void ClearStopped()
{
//for (MyTriMesh::VertexIterator vi=m.vert.begin();vi<m.vert.end();vi++)
//		ClearStopped(&*vi);

	Part_VertexContainer::iterator vi;
	for (vi=V_Stopped.begin();vi<V_Stopped.end();++vi)
	{
		ClearStopped(&*vi);
	}
	V_Stopped.clear();
}

///do one step of controls for self collision detetction
void CollisionDetection()
{
	CollDet->UpdateStep();
	std::vector<MyFace*> coll=CollDet->computeSelfIntersection();
	for (std::vector<MyFace*>::iterator it=coll.begin();it<coll.end();it++)
	{
		SetBlockedFace(*it);
		SetIntersectedFace(*it);
	}
}

public:

///set the initial barycenter where the triangle mesh start to expand
//void SetInitialBarycenter(MyTriMesh::CoordType b)
//{
//	InitialBarycenter=b;
//	gray_init=getColor(b);
//	/*InitMesh(m);*/
//}

///set the input output directory of images
void LoadFromDir(char *in, char *out)
{
	inDir=in;
	outDir=out;
	//caso optimized
	/*V.Resample(inDir,outDir);

	V.Init(1000,outDir);*/
	V.Load(inDir);
	bbox=vcg::Box3<float>(MapToSpace((V.Min())),(MapToSpace(V.Max())));
}

///set parameters for segmentation
void SetSegmentParameters(int color,int tol,float Mass=0.5f,float K_elanst=0.2f,float Dihedral=0.2f,float Time_stamp=0.2f,
					  float Edge_precision=4.f,Point3f ScaleFactor=Point3f(1.f,1.f,1.f),
					  clock_t _interval=1000,clock_t _interval2=250)
{
	mass=Mass;
	k_elanst=K_elanst;
	tolerance=tol;

	interval_reinit=_interval;
	interval_selfcollision=_interval2;

	edge_size=16.f;
	edge_precision=Edge_precision;
	time_stamp=Time_stamp;
	k_dihedral=Dihedral;
	scale=ScaleFactor;
}

///init the segmentation of the mesh
void InitSegmentation(MyTriMesh::CoordType b)
{
	InitialBarycenter=b;
	gray_init=getColor(b);

	TrINT= new myIntegrator(P_Faces,P_Vertex);
	
	TrINT->SetPDESolver(PDESolvers::EULER_METHOD);
	
	////caso optimized
	///*V.Resample(inDir,outDir);

	//V.Init(1000,outDir);*/
	//V.Load(inDir);
	//
	/*bbox=vcg::Box3<float>(MapToSpace(V.Min()),MapToSpace(V.Max()));*/
	
	InitMesh(m);
	//init the mesh with new 
	Reinit_PVectors();
	ReinitPhysicMesh();
	
	CollDet->Init(bbox.min,bbox.max,5.f);


}

///return the bounding box of the mesh
vcg::Box3<float> BBox()
{
	return (bbox);
}

///one step of moving for the deformable object
void Step(float t,float _edge_size)
{	
	if (m.face.size()!=0)
	{
		AddExtForces();
		TrINT->Step(t);
		VerifyForces();
		Refresh_PVectors();
		if (end_loop)
		{
			RefineStep(_edge_size);
			ReinitPhysicMesh();
			ClearStopped();
		}
		if (TimeSelfIntersection())
			CollisionDetection();
	}
}

void Smooth()
{
	ScaleLaplacianSmooth<MyTriMesh>(m,2,0.5);
}

void AutoStep()
{
	refined=false;
	Step(time_stamp,edge_size);
	//test on 80% of the vertex blocked
	if ((((float)P_Vertex.size()/(float)m.vn)<0.2)&&(end_loop)&&(!refined)&&(edge_size>edge_precision))
	{
				edge_size/=2.f;
				if (edge_size<edge_precision)
					edge_size=edge_precision;
				time_stamp/=2.f;
	}
}

};
#endif