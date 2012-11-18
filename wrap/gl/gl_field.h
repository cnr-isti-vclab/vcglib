#ifndef GL_FIELD
#define GL_FIELD

#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>

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

        glLineWidth(2);
        vcg::glColor(vcg::Color4b(0,0,255,255));
        glBegin(GL_LINES);
            glVertex(center);
            glVertex(center+dir[0]*size);
        glEnd();

        glLineWidth(2);
        vcg::glColor(vcg::Color4b(0,255,0,255));
        glBegin(GL_LINES);
            glVertex(center);
            glVertex(center+dir[1]*size);
        glEnd();
        /*glLineWidth(1);
        vcg::glColor(vcg::Color4b(0,0,0,255));

		glBegin(GL_LINES);
        for (int i=1;i<4;i++)
		{            
            glVertex(center);
            glVertex(center+dir[i]*size);
		}
        glEnd();*/
	}


	///draw the cross field of a given face
    static void GLDrawFaceField(const FaceType &f,
							ScalarType &size)
	{
        CoordType center=(f.cP0(0)+f.cP0(1)+f.cP0(2))/3;
		CoordType normal=f.cN();
		CoordType dir[4];
		vcg::tri::CrossField<MeshType>::CrossVector(f,dir);
		GLDrawField(dir,center,size);
	}
	
//    static void GLDrawFaceSeams(const FaceType &f,
//                                vcg::Point3<bool> seams,
//                                vcg::Color4b seamCol[3])
//    {
//        glLineWidth(2);

//        glBegin(GL_LINES);
//        for (int i=0;i<3;i++)
//        {
//            if (!seams[i])continue;
//            vcg::glColor(seamCol[i]);
//            glVertex(f.V0(i)->P());
//            glVertex(f.V1(i)->P());
//        }
//        glEnd();
//    }

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

//	///singular vertices should be selected
//    static void GLDrawSingularities(MeshType &mesh)
//	{
//        bool hasSingular = vcg::tri::HasPerVertexAttribute(mesh,std::string("Singular"));
//        bool hasSingularDegree = vcg::tri::HasPerVertexAttribute(mesh,std::string("SingularityDegree"));

//        if (!hasSingular)return;

//        typename MeshType::template PerVertexAttributeHandle<bool> Handle_Singular;
//        typename MeshType::template PerVertexAttributeHandle<int> Handle_SingularDegree;

//        Handle_Singular=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<bool>(mesh,std::string("Singular"));

//        Handle_SingularDegree=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("SingularityDegree"));

//		glPushAttrib(GL_ALL_ATTRIB_BITS);
//		glEnable(GL_COLOR_MATERIAL);
//		glDisable(GL_LIGHTING);
//		glDepthRange(0,0.999);
//        ScalarType size=10;
//		glPointSize(size);
//		glBegin(GL_POINTS);
//        for (unsigned int i=0;i<mesh.vert.size();i++)
//		{
//			if (mesh.vert[i].IsD())continue;
//            if (!Handle_Singular[i])continue;
//            int mmatch=3;
//            if (hasSingularDegree)
//                mmatch=Handle_SingularDegree[i];


//            if (mmatch==1)vcg::glColor(vcg::Color4b(0,0,255,255));
//            else
//            if (mmatch==2)vcg::glColor(vcg::Color4b(255,0,0,255));
//            else
//            if (mmatch==3)vcg::glColor(vcg::Color4b(0,255,255,255));

//            vcg::glVertex(mesh.vert[i].P());
			
//		}
//		glEnd();
//		glPopAttrib();
//	}

	static void GLDrawFaceField(const MeshType &mesh)
	{
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_LIGHTING);
        ScalarType size=mesh.bbox.Diag()/400.0;
        for (unsigned int i=0;i<mesh.face.size();i++)
		{
			if (mesh.face[i].IsD())continue;
			//if (!mesh.face[i].leading)continue;
            GLDrawFaceField(mesh.face[i],size);
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

//    static void GLDrawSeams(MeshType &mesh)
//    {
//        bool hasSeam = vcg::tri::HasPerFaceAttribute(mesh,std::string("Seams"));
//        if(!hasSeam)return;
//        bool HasSeamIndex=vcg::tri::HasPerFaceAttribute(mesh,std::string("SeamsIndex"));

//        typedef typename MeshType::template PerFaceAttributeHandle<vcg::Point3<bool> > SeamsHandleType;
//        typedef typename MeshType::template PerFaceAttributeHandle<vcg::Point3i > SeamsIndexHandleType;

//        typedef typename vcg::tri::Allocator<MeshType> SeamsAllocator;

//        SeamsHandleType Handle_Seam;
//        Handle_Seam=SeamsAllocator::template GetPerFaceAttribute<vcg::Point3<bool> >(mesh,std::string("Seams"));

//        SeamsIndexHandleType Handle_SeamIndex;
//        if (HasSeamIndex)
//        Handle_SeamIndex=SeamsAllocator::template GetPerFaceAttribute<vcg::Point3i >(mesh,std::string("SeamsIndex"));

//        glPushAttrib(GL_ALL_ATTRIB_BITS);
//        glEnable(GL_COLOR_MATERIAL);
//        glDisable(GL_LIGHTING);

//        glDepthRange(0,0.999);
//        for (unsigned int i=0;i<mesh.face.size();i++)
//        {
//            if (mesh.face[i].IsD())continue;
//            vcg::Point3<bool> seams=Handle_Seam[i];
//            vcg::Color4b seamCol[3];
//            for (int j=0;j<3;j++)
//            {
//                seamCol[j]=vcg::Color4b(0,255,0,255);
//                if (HasSeamIndex)
//                {
//                    int index=Handle_SeamIndex[i][j];
//                    //assert(index>0);
//                    if (index>=0)
//                        seamCol[j]=vcg::Color4b::Scatter(100,index);
//                }
//            }

//            GLDrawFaceSeams(mesh.face[i],seams,seamCol);

//        }
//        glPopAttrib();
//    }
};

}

#endif
