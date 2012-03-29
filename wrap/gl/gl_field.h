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
		vcg::tri::CrossField<MeshType>::CrossVector(f,dir);
		glLineWidth(1);
		vcg::Color4b c;
		vcg::glColor(vcg::Color4b(0,0,0,255));

		glBegin(GL_LINES);
		for (int i=0;i<4;i++)
		{
			glVertex(center);
			glVertex(center+dir[i]*size);
		}
		glEnd();
	}


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
			bool IsSing=vcg::tri::CrossField<MeshType>::IsSingular(mymesh.vert[i],mmatch);
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