namespace vcg{
template <class MeshType>
class GLField
{
	typedef typename MeshType::FaceType FaceType;
	typedef typename MeshType::VertexType VertexType;
	typedef typename MeshType::ScalarType ScalarType;

	///draw the cross field of a given face
	static void GLDrawField(MeshType &mesh,
							const FaceType &f,
							ScalarType &size)
	{
		CoordType center=(f.P0(0)+f.P0(1)+f.P0(2))/3;
		CoordType normal=f.cN();
		CoordType dir[4];
		vcg::tri::CrossField<MeshType>::CrossVector(mesh,f,dir);
		/*ScalarType ImportVal=ImportanceField(f);
		glLineWidth(20.f*ImportVal);*/
		glLineWidth(1);
		vcg::Color4b c;
		/*c.ColorRamp(0,1,float (ImportVal) );*/
		vcg::glColor(vcg::Color4b(0,0,0,255));

		glBegin(GL_LINES);
		for (int i=0;i<4;i++)
		{
			glVertex(center);
			glVertex(center+dir[i]*size);
		}
		glEnd();
	}


	//void GLDrawField(const VertexType &v,
	//	ScalarType &size)
	//{
	//	//if ((rand()%5)!=0)return;
	//	CoordType center=v.cP();
	//	CoordType normal=v.cN();
	//	CoordType dir[4];
	//	CoordType dir0=v.cPD1();
	//	ScalarType w0=v.cK1();
	//	ScalarType w1=v.cK2();
	//	vcg::tri::CrossField<MeshType>::CrossVector(dir0,normal,dir);
	//	ScalarType ImportVal=ImportanceField(v);
	//	//glLineWidth(10.f*ImportVal);
	//	vcg::Color4b c;
	//	//c.ColorRamp(0,1,float (ImportVal) );
	//	//vcg::glColor(c);
	//	//vcg::glColor(vcg::Color4b(100,100,100,255));
	//	if (rand()%2==0) std::swap(w0, w1);
	//	float s;
	//	if (w0>w1)
	//	{
	//		vcg::glColor(vcg::Color4b(50,50,50,255));
	//		glLineWidth(2.0);
	//		s=2;
	//	}
	//	else
	//	{
	//		vcg::glColor(vcg::Color4b(150,150,150,255));
	//		glLineWidth(1);
	//		s=1;
	//	}
	//	glBegin(GL_LINES);
	//	glVertex(center);
	//	glVertex(center+dir[0]*size*s);
	//	glVertex(center);
	//	glVertex(center+dir[2]*size*s);
	//	glEnd();
	//	if (w0<=w1)
	//	{
	//		vcg::glColor(vcg::Color4b(50,50,50,255));
	//		glLineWidth(2.0);
	//		s=2;
	//	}
	//	else
	//	{
	//		vcg::glColor(vcg::Color4b(150,150,150,255));
	//		glLineWidth(1);
	//		s=1;
	//	}
	//	glBegin(GL_LINES);
	//	glVertex(center);
	//	glVertex(center+dir[1]*size*s);
	//	glVertex(center);
	//	glVertex(center+dir[3]*size*s);
	//	glEnd();
	//}
public:

	///singular vertices should be selected
	static void GLDrawSingularities(const MeshType &mesh)
	{
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_LIGHTING);
		glDepthRange(0,0.999);
		MyScalarType size=10;
		glPointSize(size);
		glBegin(GL_POINTS);
		for (int i=0;i<mymesh.vert.size();i++)
		{
			if (mymesh.vert[i].IsD())continue;
			if (!mymesh.vert[i].IsS())continue;
			int mmatch;
			bool IsSing=vcg::tri::CrossField<MeshType>::IsSingular(mymesh,mymesh.vert[i],mmatch);
			if (!IsSing)continue;
			assert(IsSing);
			assert(mmatch!=0);
			/*vcg::glColor(vcg::Color4b(255,0,0,255));*/
			if (mmatch==1)vcg::glColor(vcg::Color4b(0,0,255,255));
			else
			if (mmatch==2)vcg::glColor(vcg::Color4b(255,0,0,255));
			else
			if (mmatch==3)vcg::glColor(vcg::Color4b(0,255,255,255));
			
			vcg::glVertex(mymesh.vert[i].P());
			
		}
		glEnd();
		glPopAttrib();
	}
	static void GLDrawFaceField(const MeshType &mesh)
	{
		srand(12345);
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_LIGHTING);
		MyScalarType size=mymesh.bbox.Diag()/100.0;
		vcg::Color4b c=vcg::Color4b(255,0,0,255);
		for (int i=0;i<mymesh.face.size();i++)
		{
			if (mymesh.face[i].IsD())continue;
			GLDrawField(mymesh,mymesh.face[i],size);
		}

		glPopAttrib();
	}
};
}