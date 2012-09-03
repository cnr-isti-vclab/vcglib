namespace vcg{
template <class MeshType>
class GLField
{
	typedef typename MeshType::FaceType FaceType;
	typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
	typedef typename MeshType::ScalarType ScalarType;
	
	static void GLDrawField(CoordType dir[4],
							CoordType center,
							ScalarType &size)
	{

        glLineWidth(3);
        vcg::glColor(vcg::Color4b(0,0,255,255));
        glBegin(GL_LINES);
            glVertex(center);
            glVertex(center+dir[0]*size);
        glEnd();

        glLineWidth(1);
        vcg::glColor(vcg::Color4b(0,0,0,255));

		glBegin(GL_LINES);
        for (int i=1;i<4;i++)
		{            
            glVertex(center);
            glVertex(center+dir[i]*size);
		}
		glEnd();
	}


	///draw the cross field of a given face
	static void GLDrawFaceField(const MeshType &mesh,
							const FaceType &f,
							ScalarType &size)
	{
		CoordType center=(f.P0(0)+f.P0(1)+f.P0(2))/3;
		CoordType normal=f.cN();
		CoordType dir[4];
		vcg::tri::CrossField<MeshType>::CrossVector(f,dir);
		GLDrawField(dir,center,size);
	}
	
    static void GLDrawFaceSeams(const FaceType &f)
    {
        glLineWidth(3);

        glBegin(GL_LINES);
        for (int i=0;i<3;i++)
        {
            if (!f.IsSeam(i))continue;

            glVertex(f.V0(i)->P());
            glVertex(f.V1(i)->P());
        }
        glEnd();
    }

	static void GLDrawVertField(const MeshType &mesh,
								const VertexType &v,
								ScalarType &size)
	{
		CoordType center=v.cP();
		CoordType normal=v.cN();
		CoordType dir[4];
		vcg::tri::CrossField<MeshType>::CrossVector(v,dir);
		GLDrawField(dir,center,size);
	}

public:

	///singular vertices should be selected
	static void GLDrawSingularities(const MeshType &mesh)
	{
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_LIGHTING);
		glDepthRange(0,0.999);
        ScalarType size=10;
		glPointSize(size);
		glBegin(GL_POINTS);
		for (int i=0;i<mesh.vert.size();i++)
		{
			if (mesh.vert[i].IsD())continue;
			if (!mesh.vert[i].IsS())continue;
			int mmatch;
			bool IsSing=vcg::tri::CrossField<MeshType>::IsSingular(mesh.vert[i],mmatch);
			if (!IsSing)continue;
			assert(IsSing);
			assert(mmatch!=0);
			/*vcg::glColor(vcg::Color4b(255,0,0,255));*/
			if (mmatch==1)vcg::glColor(vcg::Color4b(0,0,255,255));
			else
			if (mmatch==2)vcg::glColor(vcg::Color4b(255,0,0,255));
			else
			if (mmatch==3)vcg::glColor(vcg::Color4b(0,255,255,255));
			
			vcg::glVertex(mesh.vert[i].P());
			
		}
		glEnd();
		glPopAttrib();
	}

	static void GLDrawFaceField(const MeshType &mesh)
	{
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_LIGHTING);
        ScalarType size=mesh.bbox.Diag()/400.0;
        for (int i=0;i<mesh.face.size();i++)
		{
			if (mesh.face[i].IsD())continue;
			//if (!mesh.face[i].leading)continue;
			GLDrawFaceField(mesh,mesh.face[i],size);
		}
		glPopAttrib();
	}

	static void GLDrawVertField(const MeshType &mesh)
	{
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_LIGHTING);
        ScalarType size=mesh.bbox.Diag()/100.0;
        for (int i=0;i<mesh.vert.size();i++)
		{
			if (mesh.vert[i].IsD())continue;
			//if (!mesh.face[i].leading)continue;
			GLDrawVertField(mesh,mesh.vert[i],size);
		}
		glPopAttrib();
	}

    static void GLDrawSeams(const MeshType &mesh)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glEnable(GL_COLOR_MATERIAL);
        glDisable(GL_LIGHTING);
        vcg::glColor(vcg::Color4b(255,0,0,255));
        glDepthRange(0,0.999);
        for (int i=0;i<mesh.face.size();i++)
        {
            if (mesh.face[i].IsD())continue;
            GLDrawFaceSeams(mesh.face[i]);
        }
        glPopAttrib();
    }
};

}
