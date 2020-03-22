#ifndef GL_FIELD
#define GL_FIELD

#include <wrap/gl/space.h>
#include <wrap/gl/math.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
#include <vcg/complex/allocate.h>

namespace vcg{
template <class MeshType>
class GLField
{
	typedef typename MeshType::FaceType FaceType;
	typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
	typedef typename MeshType::ScalarType ScalarType;
	
public:

	static void GLDrawField(CoordType dir[4],
                            const CoordType &center,
                            const ScalarType &size,
                            const ScalarType &Width0,
                            const ScalarType &Width1,
                            const vcg::Color4b &Color0,
                            const vcg::Color4b &Color1,
                            bool oneside,
                            bool onlyPD1)
    {
        CoordType dirN[4];
        for (size_t i=0;i<4;i++)
        {
            dirN[i]=dir[i];
            dirN[i].Normalize();
        }

        ScalarType size1=size;
        if (oneside)size1=0;

        glLineWidth(Width0);
        vcg::glColor(Color0);
        glBegin(GL_LINES);
            glVertex(center+dirN[0]*size);
            glVertex(center+dirN[2]*size1);
        glEnd();

        if (onlyPD1)return;

        glLineWidth(Width1);
        vcg::glColor(Color1);
        glBegin(GL_LINES);
            glVertex(center+dirN[1]*size);
            glVertex(center+dirN[3]*size1);
        glEnd();

	}

//    ///draw the cross field of a given face in a given position
//    static void GLDrawSingleFaceField(const FaceType &f,
//                                     CoordType pos,
//                                     ScalarType &size,
//                                     bool onlyPD1,
//                                     bool oneside)
//    {
//        CoordType center=pos;
//        CoordType normal=f.cN();
//        CoordType dir[4];
//        vcg::tri::CrossField<MeshType>::CrossVector(f,dir);
//        GLDrawField(dir,center,size,onlyPD1,oneside);
//    }

	///draw the cross field of a given face
    static void GLDrawSingleFaceField(const FaceType &f,
                                      const ScalarType &size,
                                      const bool oneside,
                                      const bool onlyPD1,
                                      const ScalarType maxN,
                                      const ScalarType minN,
                                      const bool UseK)
	{
        assert(maxN>=minN);
        CoordType center=(f.cP(0)+f.cP(1)+f.cP(2))/3;
        //CoordType normal=f.cN();
		CoordType dir[4];
		vcg::tri::CrossField<MeshType>::CrossVector(f,dir);

        if (maxN<=0)
            GLDrawField(dir,center,size,2,2,vcg::Color4b(0,0,0,255),vcg::Color4b(0,0,0,255),oneside,onlyPD1);
        else
        {
            ScalarType Norm0=dir[0].Norm();
            ScalarType Norm1=dir[1].Norm();
            if (UseK)
            {
                Norm0=f.cK1();
                Norm1=f.cK2();
            }
            ScalarType MaxW=6;
            ScalarType MinW=0.5;
            ScalarType IntervW=MaxW-MinW;
            if (Norm0>maxN)Norm0=maxN;
            if (Norm1>maxN)Norm1=maxN;
            if (Norm0<minN)Norm0=minN;
            if (Norm1<minN)Norm1=minN;

            vcg::Color4b Col0,Col1;
            ScalarType W0,W1;
            if (!UseK)
            {
               Col0=vcg::Color4b::ColorRamp(minN,maxN,Norm0);
               Col1=vcg::Color4b::ColorRamp(minN,maxN,Norm1);
               W0=(Norm0/(maxN-minN))*IntervW+MinW;
               W1=(Norm1/(maxN-minN))*IntervW+MinW;
            }
            else
            {
               ScalarType MaxAbs=std::max(fabs(minN),fabs(maxN));
               Col0=vcg::Color4b::ColorRamp(-MaxAbs,MaxAbs,Norm0);
               Col1=vcg::Color4b::ColorRamp(-MaxAbs,MaxAbs,Norm1);
//               if (Norm0<0)
//               {
//                   assert(minN<0);
//                   //PUT green on ZERO
//                   Col0=vcg::Color4b::ColorRamp(minN,fabs(minN),Norm0);

//               }else
//               {
//                   //PUT green on ZERO
//                   Col0=vcg::Color4b::ColorRamp(-maxN,maxN,Norm0);
//               }

//               if (Norm1<0)
//               {
//                   assert(minN<0);
//                   //PUT green on ZERO
//                   Col1=vcg::Color4b::ColorRamp(minN,fabs(minN),Norm1);
//               }else
//               {
//                   //PUT green on ZERO
//                   Col1=vcg::Color4b::ColorRamp(-maxN,maxN,Norm1);
//               }
               W0=(fabs(Norm0)/std::max(fabs(maxN),fabs(minN)))*IntervW+MinW;
               W1=(fabs(Norm1)/std::max(fabs(maxN),fabs(minN)))*IntervW+MinW;

            }
            GLDrawField(dir,center,size,W0,W1,Col0,Col1,oneside,onlyPD1);
        }
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

    static void GLDrawVertField(const VertexType &v,
                                ScalarType &size)
    {
        CoordType center=v.cP();
        CoordType normal=v.cN();
        CoordType dir[4];
        vcg::tri::CrossField<MeshType>::CrossVector(v,dir);
        GLDrawField(dir,center,size,2,2,vcg::Color4b(0,0,0,255),vcg::Color4b(0,0,0,255),false,false);
    }


    static void GLDrawFaceField(const MeshType &mesh,
                                bool onlyPD1,
                                bool oneside,
                                ScalarType GlobalScale=0.002,
                                const ScalarType maxN=0,
                                const ScalarType minN=0,
                                bool UseK=false)
	{

		glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDepthRange(0.0,0.999);
		glEnable(GL_COLOR_MATERIAL);
        glDisable(GL_LIGHTING);
        glDisable(GL_BLEND);
        ScalarType size=mesh.bbox.Diag()*GlobalScale;
        for (unsigned int i=0;i<mesh.face.size();i++)
		{
            if (mesh.face[i].IsD())continue;
            GLDrawSingleFaceField(mesh.face[i],size,oneside,onlyPD1,maxN,minN,UseK);
		}
		glPopAttrib();
	}

    static void GLDrawVertField(const MeshType &mesh,ScalarType sizeF=0.01)
	{
		glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDepthRange(0.0,0.9999);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_LIGHTING);
        glDisable(GL_BLEND);
        ScalarType size=mesh.bbox.Diag()*sizeF;
        for (int i=0;i<mesh.vert.size();i++)
		{
            if (mesh.vert[i].IsD())continue;
            GLDrawVertField(mesh.vert[i],size);
		}
		glPopAttrib();
	}

    static void GLDrawSingularity(MeshType &mesh)
    {
        // query if an attribute is present or not
       bool hasSingular = vcg::tri::HasPerVertexAttribute(mesh,std::string("Singular"));
       bool hasSingularIndex = vcg::tri::HasPerVertexAttribute(mesh,std::string("SingularIndex"));

       if (!hasSingular)return;
       if(!hasSingularIndex)return;

       typename MeshType::template PerVertexAttributeHandle<bool> Handle_Singular;
       Handle_Singular=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<bool>(mesh,std::string("Singular"));
       typename MeshType::template PerVertexAttributeHandle<int> Handle_SingularIndex;
       Handle_SingularIndex =vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("SingularIndex"));

       glPushAttrib(GL_ALL_ATTRIB_BITS);

       glDepthRange(0.0,0.9999);
       glEnable(GL_COLOR_MATERIAL);
       glDisable(GL_LIGHTING);
       glDisable(GL_BLEND);
       glPointSize(20);
       glBegin(GL_POINTS);
       for (size_t i=0;i<mesh.vert.size();i++)
       {
           if (mesh.vert[i].IsD())continue;
           if (!Handle_Singular[i])continue;


           int SingIndex=Handle_SingularIndex[i];

           vcg::Color4b colSing;

           switch (SingIndex)
           {
             case 1:colSing=vcg::Color4b(0,0,255,255);      break;
             case 2:colSing=vcg::Color4b(0,255,0,255);    break;
             case 3:colSing=vcg::Color4b(255,0,0,255);      break;
             case 4:colSing=vcg::Color4b(255,255,0,255);      break;
             default:colSing=vcg::Color4b(255,0,255,255);
           }


           vcg::glColor(colSing);
           vcg::glVertex(mesh.vert[i].P());
       }
       glEnd();
       glPopAttrib();
    }
};

}

#endif
