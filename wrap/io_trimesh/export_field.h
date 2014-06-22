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
#ifndef __VCGLIB_EXPORTERFIELD
#define __VCGLIB_EXPORTERFIELD

namespace vcg {
namespace tri {
namespace io {

/** 
This class encapsulate a filter for saving field formats
*/
template <class MeshType>
class ExporterFIELD
{
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;

public:
    
    ///load a field on the mesh, it could be a vfield file (per vertex)
    ///or an ffield file (per face)
    static void SaveFaceFIELD(MeshType &mesh,
                              const char *path)
    {
        
        FILE *f = fopen(path,"wt");
        //if (!f)return false;
//            char word[512]; word[0]=0;
//            fscanf(f,"%s",word);
//            char c=0;
//            if (word[0]=='#') {
//                // skip comment line
//                while (fscanf(f,"%c",&c)!=EOF) if (c=='\n') break;
//            }
//            else
//            {
//                return false;
//            }
            int nf = mesh.fn;//-1;
            fprintf(f,"# frame generated with VCG \n");
            fprintf(f,"target frame \n");
            fprintf(f,"%d\n",nf);
//            if (fscanf(f,"%d",&nnv)!=1)
//            {
//                while (fscanf(f,"%c",&c)!=EOF) if (c=='\n') break; // skip
//                fscanf(f,"%d",&nnv);
//            }
//            int targetnum=mesh.fn;
//            if (per_vertex)
//                targetnum=mesh.vn;
//            if (nnv != (int)targetnum)
//            {
//                //if (errorMsg) sprintf(errorMsg,"Wrong element number. Found: %d. Expected: %d.",nnv,mesh->vn);
//                return false;
//            }
//            
//            if( per_vertex && !HasPerVertexCurvatureDir(mesh)) throw vcg::MissingComponentException("PerVertexCurvatureDir");
//            if(!per_vertex && !HasPerFaceCurvatureDir(mesh))   throw vcg::MissingComponentException("PerFaceCurvatureDir");
            if (!HasPerFaceCurvatureDir(mesh))
                throw vcg::MissingComponentException("PerFaceCurvatureDir");
            
            fprintf(f,"k1	 k2	 k1v_x	 k1v_y	 k1v_z	 k2v_x	 k2v_y	 k2v_z\n");
//            while (fscanf(f,"%c",&c)!=EOF) if (c=='\n') break; // skip
//            // skip strange string line
//            while (fscanf(f,"%c",&c)!=EOF) if (c=='\n') break;
            for (int i=0; i<nf; i++){
                vcg::Point3<float> u;
                u.Import(mesh.face[i].PD1());
                vcg::Point3<float> v;
                v.Import(mesh.face[i].PD2());
                
                fprintf(f,"1 1 %f %f %f %f %f %f\n",
                        (u.X()),(u.Y()),(u.Z()),
                        (v.X()),(v.Y()),(v.Z()));
//                if (fscanf(f,
//                           "%f %f %f %f %f %f %f %f",
//                           &a,&b,
//                           &(v.X()),&(v.Y()),&(v.Z()),
//                           &(u.X()),&(u.Y()),&(u.Z())
//                           )!=8) {
//                    //if (errorMsg) sprintf(errorMsg,"Format error reading vertex n. %d",i);
//                    return false;
//                }
//                
//                u.Normalize();
//                v.Normalize();
//                
//                if (per_vertex)
//                {
//                    mesh.vert[i].PD1().Import(u);
//                    mesh.vert[i].PD2().Import(v);
//                }
//                else
//                {
//                    mesh.face[i].PD1().Import(u);
//                    mesh.face[i].PD2().Import(v);
//                }
            }
        fclose(f);
    }

    
    ///Save a 4 rosy format file as used by
    ///Interactive Visualization of Rotational Symmetry Fields on Surfaces
    ///Jonathan Palacios and Eugene Zhang
    static void Save4ROSY(MeshType &mesh,
                        const char *path)
    {
        FILE *f = fopen(path,"wt");
        fprintf(f,"%d\n",mesh.vn);
        fprintf(f,"4\n");
        for (unsigned int i=0;i<mesh.vert.size();i++)
        {
            float dirX=(float)mesh.vert[i].PD1().X();
            float dirY=(float)mesh.vert[i].PD1().Y();
            float dirZ=(float)mesh.vert[i].PD1().Z();
            fprintf(f,"%f %f %f \n",dirX,dirY,dirZ);

        }
        fclose(f);
    }

    ///Save a 4 rosy format file as used by
    ///Interactive Visualization of Rotational Symmetry Fields on Surfaces
    ///Jonathan Palacios and Eugene Zhang
    static void Save4ROSYFace(MeshType &mesh,
                        const char *path)
    {
        FILE *f = fopen(path,"wt");
        fprintf(f,"%d\n",mesh.vn);
        fprintf(f,"4\n");
        for (unsigned int i=0;i<mesh.face.size();i++)
        {
            float dirX=(float)mesh.face[i].PD1().X();
            float dirY=(float)mesh.face[i].PD1().Y();
            float dirZ=(float)mesh.face[i].PD1().Z();
            fprintf(f,"%f %f %f \n",dirX,dirY,dirZ);

        }
        fclose(f);
    }

}; // end class



} // end namespace tri
} // end namespace io
} // end namespace vcg

#endif

