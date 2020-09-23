/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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
Revision 1.1  2004/06/03 13:16:32  ganovelli
created

****************************************************************************/
#ifndef __VCGLIB_TETRAIMPORTER_PLY
#define __VCGLIB_TETRAIMPORTER_PLY

#include<wrap/callback.h>
#include<wrap/ply/plylib.h>
#include<wrap/io_tetramesh/io_mask.h>
#include<wrap/io_tetramesh/io_ply.h>

namespace vcg {
namespace tetra {
namespace io {

template <class TYPE>
int PlyType ()  { return 0;}

template <> inline int PlyType <float >()  { return ply::T_FLOAT; }
template <> inline int PlyType <double>()  { return ply::T_DOUBLE; }
template <> inline int PlyType <int   >()  { return ply::T_INT; } 
template <> inline int PlyType <short >()  { return ply::T_SHORT; }
template <> inline int PlyType <unsigned char >()  { return ply::T_UCHAR; }

/** 
This class encapsulate a filter for opening ply meshes.
The ply file format is quite extensible...
*/
template <class OpenMeshType>
class ImporterPLY
{
public:

    typedef ::vcg::ply::PropDescriptor PropDescriptor ;
    typedef typename OpenMeshType::VertexPointer VertexPointer;
    typedef typename OpenMeshType::ScalarType ScalarType;
    typedef typename OpenMeshType::VertexType VertexType;
    typedef typename OpenMeshType::TetraType TetraType;
    typedef typename OpenMeshType::VertexIterator VertexIterator;
    typedef typename OpenMeshType::TetraIterator TetraIterator;

    //template <class T> int PlyType () {	assert(0);  return 0;}

#define MAX_USER_DATA 256
    // Struttura ausiliaria per la lettura del file ply
    struct LoadPly_TetraAux
    {
        unsigned char size;
        int v[512];
        int flags;
        float q;
        float texcoord[32];
        unsigned char ntexcoord;
        int texcoordind;
        float colors[32];
        unsigned char ncolors;

        unsigned char r;
        unsigned char g;
        unsigned char b;
        unsigned char a;

        unsigned char data[MAX_USER_DATA];
    };

    //struct LoadPly_TristripAux
    //{
    //	int size;
    //	int *v;
    //	unsigned char data[MAX_USER_DATA];
    //};

    // Struttura ausiliaria per la lettura del file ply
    template<class S>
    struct LoadPly_VertAux
    {
        S p[3];
        int flags;
        float q;
        unsigned char r;
        unsigned char g;
        unsigned char b;
        unsigned char a;
        unsigned char data[MAX_USER_DATA];
    };

    // Struttura ausiliaria caricamento camera
    //struct LoadPly_Camera
    //{
    //	float view_px;
    //	float view_py;
    //	float view_pz;
    //	float x_axisx;
    //	float x_axisy;
    //	float x_axisz;
    //	float y_axisx;
    //	float y_axisy;
    //	float y_axisz;
    //	float z_axisx;
    //	float z_axisy;
    //	float z_axisz;
    //	float focal;
    //	float scalex;
    //	float scaley;
    //	float centerx;
    //	float centery;
    //	int   viewportx;
    //	int   viewporty;
    //	float k1;
    //	float k2;
    //	float k3;
    //	float k4;
    //};

    static const  PropDescriptor &VertDesc(int i)
    {
        const static PropDescriptor pv[13]={
            /*00*/ {"vertex", "x",         ply::T_FLOAT, PlyType<ScalarType>(),          offsetof(LoadPly_VertAux<ScalarType>,p),0,0,0,0,0},
            /*01*/ {"vertex", "y",         ply::T_FLOAT, PlyType<ScalarType>(),          offsetof(LoadPly_VertAux<ScalarType>,p) +     sizeof(ScalarType),0,0,0,0,0},
            /*02*/ {"vertex", "z",         ply::T_FLOAT, PlyType<ScalarType>(),          offsetof(LoadPly_VertAux<ScalarType>,p) + 2 * sizeof(ScalarType),0,0,0,0,0},
            /*03*/ {"vertex", "flags",     ply::T_INT,   ply::T_INT,            offsetof(LoadPly_VertAux<ScalarType>,flags),0,0,0,0,0},
            /*04*/ {"vertex", "quality",   ply::T_FLOAT, ply::T_FLOAT,          offsetof(LoadPly_VertAux<ScalarType>,q),0,0,0,0,0},
            /*05*/ {"vertex", "red"  ,     ply::T_UCHAR, ply::T_UCHAR,          offsetof(LoadPly_VertAux<ScalarType>,r),0,0,0,0,0},
            /*06*/ {"vertex", "green",     ply::T_UCHAR, ply::T_UCHAR,          offsetof(LoadPly_VertAux<ScalarType>,g),0,0,0,0,0},
            /*07*/ {"vertex", "blue" ,     ply::T_UCHAR, ply::T_UCHAR,          offsetof(LoadPly_VertAux<ScalarType>,b),0,0,0,0,0},
            /*08*/ {"vertex", "alpha" ,    ply::T_UCHAR, ply::T_UCHAR,         offsetof(LoadPly_VertAux<ScalarType>,a),0,0,0,0,0},
            // DOUBLE
            /*09*/ {"vertex", "x",         ply::T_DOUBLE, PlyType<ScalarType>(),offsetof(LoadPly_VertAux<ScalarType>,p),0,0,0,0,0  ,0},
            /*10*/ {"vertex", "y",         ply::T_DOUBLE, PlyType<ScalarType>(),offsetof(LoadPly_VertAux<ScalarType>,p) + sizeof(ScalarType)  ,0,0,0,0,0  ,0},
            /*11*/ {"vertex", "z",         ply::T_DOUBLE, PlyType<ScalarType>(),offsetof(LoadPly_VertAux<ScalarType>,p) + 2*sizeof(ScalarType),0,0,0,0,0  ,0},
            /*12*/ {"vertex", "quality",   ply::T_DOUBLE, PlyType<ScalarType>(),offsetof(LoadPly_VertAux<ScalarType>,q),0,0,0,0,0  ,0}
        };
        return pv[i];
    }


    static const  PropDescriptor &TetraDesc(int i)
    {
        const static 	PropDescriptor qf[10]=
        {
            {"tetra", "vertex_indices", ply::T_INT,   ply::T_INT,   offsetof(LoadPly_TetraAux,v),		     1,0,ply::T_UCHAR,ply::T_UCHAR,offsetof(LoadPly_TetraAux,size) },
            {"tetra", "flags",          ply::T_INT,   ply::T_INT,   offsetof(LoadPly_TetraAux,flags),     0,0,0,0,0},
            {"tetra", "quality",        ply::T_FLOAT, ply::T_FLOAT, offsetof(LoadPly_TetraAux,q),         0,0,0,0,0},
            {"tetra", "texcoord",       ply::T_FLOAT, ply::T_FLOAT, offsetof(LoadPly_TetraAux,texcoord),    1,0,ply::T_UCHAR,ply::T_UCHAR,offsetof(LoadPly_TetraAux,ntexcoord) },
            {"tetra", "color",          ply::T_FLOAT, ply::T_FLOAT, offsetof(LoadPly_TetraAux,colors),     1,0,ply::T_UCHAR,ply::T_UCHAR,offsetof(LoadPly_TetraAux,ncolors) },
            {"tetra", "texnumber",      ply::T_INT,   ply::T_INT,   offsetof(LoadPly_TetraAux,texcoordind), 0,0,0,0,0},
            {"tetra", "red"  ,          ply::T_UCHAR, ply::T_UCHAR, offsetof(LoadPly_TetraAux,r),         0,0,0,0,0},
            {"tetra", "green",          ply::T_UCHAR, ply::T_UCHAR, offsetof(LoadPly_TetraAux,g),         0,0,0,0,0},
            {"tetra", "blue" ,          ply::T_UCHAR, ply::T_UCHAR, offsetof(LoadPly_TetraAux,b),         0,0,0,0,0},
            {"tetra", "alpha" ,         ply::T_UCHAR, ply::T_UCHAR, offsetof(LoadPly_TetraAux,a),         0,0,0,0,0},

        };
        return qf[i];
    }



    /// Standard call for reading a mesh
    static int Open( OpenMeshType &m, const char * filename, CallBackPos *cb=0)
    {
        PlyInfo pi;
        pi.cb=cb;
        return Open(m, filename, pi);
    }

    /// Read a mesh and store in loadmask the loaded field
    static int Open( OpenMeshType &m, const char * filename, int & loadmask, CallBackPos *cb =0)
    {
        PlyInfo pi;
        pi.mask=loadmask;
        return Open(m, filename,pi);
        loadmask=pi.mask;
    }


    /// read a mesh with all the possible option specified in the PlyInfo obj.
    static int Open( OpenMeshType &m, const char * filename, PlyInfo & pi )
    {
        assert(filename!=0);
        vector<VertexPointer> index;

        LoadPly_TetraAux fa;
        LoadPly_VertAux<ScalarType> va;

        pi.mask = 0;

        va.flags = 42;

        pi.status = ::vcg::ply::E_NOERROR;

        // The main descriptor of the ply file
        vcg::ply::PlyFile pf;

        // Open the file and parse the header
        std::cerr << "Parsing the header.." << std::endl;
        if( pf.Open(filename,vcg::ply::PlyFile::MODE_READ)== -1 )
        {
            pi.status = pf.GetError();
            return -1;
        }
        std::cerr << "Parsed the header.." << std::endl;
        pi.header = pf.GetHeader();

        // Descrittori dati standard (vertex coord e faces)
        if( pf.AddToRead(VertDesc(0))  == -1 && pf.AddToRead(VertDesc(9))  == -1 ) { pi.status = PlyInfo::E_NO_VERTEX; return pi.status; }
        if( pf.AddToRead(VertDesc(1))  == -1 && pf.AddToRead(VertDesc(10)) == -1 ) { pi.status = PlyInfo::E_NO_VERTEX; return pi.status; }
        if( pf.AddToRead(VertDesc(2))  == -1 && pf.AddToRead(VertDesc(11)) == -1 ) { pi.status = PlyInfo::E_NO_VERTEX; return pi.status; }
        if( pf.AddToRead(TetraDesc(0)) == -1 ) { pi.status = PlyInfo::E_NO_VERTEX; return pi.status; }

        // Descrittori facoltativi dei flags
        if( pf.AddToRead(VertDesc(3)) != -1 )
            pi.mask |= vcg::tetra::io::Mask::IOM_VERTFLAGS;
        if( VertexType::HasQuality() )
        {
            if( pf.AddToRead(VertDesc(4)) != -1 ||
                pf.AddToRead(VertDesc(9)) != -1 )
                pi.mask |= vcg::tetra::io::Mask::IOM_VERTQUALITY;
        }
        if( VertexType::HasColor() )
        {
            if( pf.AddToRead(VertDesc(5))!=-1 )
            {
                pf.AddToRead(VertDesc(6));
                pf.AddToRead(VertDesc(7));
                pf.AddToRead(VertDesc(8));

                pi.mask |= vcg::tetra::io::Mask::IOM_VERTCOLOR;
            }
        }
        // se ci sono i flag per vertice ci devono essere anche i flag per faccia
        if( pf.AddToRead(TetraDesc(1))!=-1 )
            pi.mask |= vcg::tetra::io::Mask::IOM_TETRAFLAGS;
        if( TetraType::HasQuality())
        {
            if( pf.AddToRead(TetraDesc(2))!=-1 )
                pi.mask |= vcg::tetra::io::Mask::IOM_TETRAQUALITY;
        }
        if( TetraType::HasColor() )
        {
            if( pf.AddToRead(TetraDesc(6))!=-1 )
            {
                pf.AddToRead(TetraDesc(7));
                pf.AddToRead(TetraDesc(8));
                pf.AddToRead(TetraDesc(9));
                pi.mask |= vcg::tetra::io::Mask::IOM_TETRACOLOR;
            }
        }
        // Descrittori definiti dall'utente,
        vector<PropDescriptor> VPV(pi.vdn); // property descriptor relative al tipo LoadPly_VertexAux
        vector<PropDescriptor> FPV(pi.fdn); // property descriptor relative al tipo LoadPly_FaceAux

        if(pi.vdn>0){
            // Compute the total size needed to load additional per vertex data.
            size_t totsz=0;
            for(int i=0;i<pi.vdn;i++){
                VPV[i] = pi.VertexData[i];
                VPV[i].offset1=offsetof(LoadPly_VertAux<ScalarType>,data)+totsz;
                totsz+=pi.VertexData[i].memtypesize();
                if( pf.AddToRead(VPV[i])==-1 ) { pi.status = pf.GetError(); return -1; }
            }
            if(totsz > MAX_USER_DATA)
            {
                pi.status = vcg::ply::E_BADTYPE;
                return -1;
            }
        }
        if(pi.fdn>0){
            size_t totsz=0;
            for(int i=0;i<pi.fdn;i++){
                FPV[i] = pi.TetraData[i];
                FPV[i].offset1=offsetof(LoadPly_TetraAux,data)+totsz;
                totsz+=pi.TetraData[i].memtypesize();
                if( pf.AddToRead(FPV[i])==-1 ) { pi.status = pf.GetError(); return -1; }
            }
            if(totsz > MAX_USER_DATA)
            {
                pi.status = vcg::ply::E_BADTYPE;
                return -1;
            }
        }
        /**************************************************************/
        /* Main Reading Loop */
        /**************************************************************/
        m.Clear();
        for(int i=0;i<int(pf.elements.size());i++)
        {
            int n = pf.ElemNumber(i);

            if( !strcmp( pf.ElemName(i),"vertex" ) )
            {
                int j;

                pf.SetCurElement(i);
                VertexIterator vi=Allocator<OpenMeshType>::AddVertices(m,n);

                for(j=0;j<n;++j)
                {
                    if(pi.cb && (j%1000)==0) pi.cb(j*50/n,"Vertex Loading");
                    if( pf.Read( (void *)&(va) )==-1 )
                    {
                        pi.status = PlyInfo::E_SHORTFILE;
                        return -1;
                    }

                    (*vi).P()[0] = va.p[0];
                    (*vi).P()[1] = va.p[1];
                    (*vi).P()[2] = va.p[2];

                    if( pi.mask & vcg::tetra::io::Mask::IOM_VERTFLAGS )
                        (*vi).Flags() = va.flags;

                    if( pi.mask & vcg::tetra::io::Mask::IOM_VERTQUALITY )
                        (*vi).Q() = va.q;
                    if( pi.mask & vcg::tetra::io::Mask::IOM_VERTCOLOR )
                    {

                        (*vi).C()[0] = va.r;
                        (*vi).C()[1] = va.g;
                        (*vi).C()[2] = va.b;
                        (*vi).C()[3] = va.a;
                    }


                    for(int k=0;k<pi.vdn;k++)
                        memcpy((char *)(&*vi) + pi.VertexData[k].offset1,
                               (char *)(&va) + VPV[k].offset1,
                               VPV[k].memtypesize());
                    ++vi;
                }

                index.resize(n);
                for(j=0,vi=m.vert.begin();j<n;++j,++vi)
                    index[j] = &*vi;
            }
            else if( !strcmp( pf.ElemName(i),"tetra") )/************************************************************/
            {
                int j;
                int k;
                TetraIterator ti = Allocator<OpenMeshType>::AddTetras(m, n);
                pf.SetCurElement(i);

                for(j=0;j<n;++j)
                {


                    if(pi.cb && (j%1000)==0) pi.cb(50+j*50/n,"Tetra Loading");
                    if( pf.Read(&fa)==-1 )
                    {
                        pi.status = PlyInfo::E_SHORTFILE;
                        return -1;
                    }
                    if(fa.size!=4)
                    {
                        pi.status = PlyInfo::E_NO_3VERTINFACE;
                        return -1;
                    }

                    for(k=0;k<4;++k)
                    {
                        if( fa.v[k]<0 || fa.v[k]>=m.vn )
                        {
                            pi.status = PlyInfo::E_BAD_VERT_INDEX;
                            return -1;
                        }
                        (*ti).V(k) = index[ fa.v[k] ];
                    }

                    if( pi.mask & vcg::tetra::io::Mask::IOM_TETRAFLAGS )
                    {
                        (*ti).Flags() = fa.flags;
                    }

                    if( pi.mask & vcg::tetra::io::Mask::IOM_TETRAQUALITY )
                    {
                        (*ti).Q() = fa.q;
                    }

                    if( pi.mask & vcg::tetra::io::Mask::IOM_TETRACOLOR )
                    {
                        (*ti).C()[0] = fa.r;
                        (*ti).C()[1] = fa.g;
                        (*ti).C()[2] = fa.b;
                        (*ti).C()[3] = fa.a;
                    }
                    ++ti;
                }

                for(k=0;k<pi.fdn;k++)
                    memcpy((char *)(&(*ti)) + pi.TetraData[k].offset1,
                           (char *)(&fa) + FPV[k].offset1,
                           FPV[k].memtypesize());
            }
            else
            {
                // Skippaggio elementi non gestiti
                int n = pf.ElemNumber(i);
                pf.SetCurElement(i);

                for(int j=0;j<n;j++)
                {
                    if( pf.Read(0)==-1)
                    {
                        pi.status = PlyInfo::E_SHORTFILE;
                        return -1;
                    }
                }
            }
        }

        // // Parsing texture names
        //textures.clear();
        //normalmaps.clear();

        //for(int co=0;co<int(pf.comments.size());++co)
        //{
        //	const char * TFILE = "TextureFile";
        //	const char * NFILE = "TextureNormalFile";
        //	const char * c = pf.comments[co];
        //	char buf[256];
        //	int i,j,n;

        //	if( !strncmp(c,TFILE,strlen(TFILE)) )
        //	{
        //		strcpy(buf,c+strlen(TFILE)+1);
        //		n = strlen(buf);
        //		for(i=j=0;i<n;i++)
        //			if( buf[i]!=' ' && buf[i]!='\t' && buf[i]>32 && buf[i]<125 )	buf[j++] = buf[i];
        //
        //		buf[j] = 0;
        //		char buf2[255];
        //		__interpret_texture_name( buf,filename,buf2 );
        //		textures.push_back( xstring(buf2) );
        //	}
        //	if( !strncmp(c,NFILE,strlen(NFILE)) )
        //	{
        //		strcpy(buf,c+strlen(NFILE)+1);
        //		n = strlen(buf);
        //		for(i=j=0;i<n;i++)
        //			if( buf[i]!=' ' && buf[i]!='\t' && buf[i]>32 && buf[i]<125 )	buf[j++] = buf[i];
        //
        //		buf[j] = 0;
        //		char buf2[255];
        //		__interpret_texture_name( buf,filename,buf2 );
        //		normalmaps.push_back( xstring(buf2) );
        //	}
        //}

        // vn and fn should be correct but if someone wrongly saved some deleted elements they can be wrong.
        m.vn = 0;
        VertexIterator vi;
        for(vi=m.vert.begin();vi!=m.vert.end();++vi)
            if( ! (*vi).IsD() )
                ++m.vn;

        m.tn = 0;
        TetraIterator fi;
        for(fi=m.tetra.begin();fi!=m.tetra.end();++fi)
            if( ! (*fi).IsD() )
                ++m.tn;

        return 0;
    }


//    // Caricamento camera da un ply
//    int LoadCamera(const char * filename)
//    {
//        vcg::ply::PlyFile pf;
//        if( pf.Open(filename,vcg::ply::PlyFile::MODE_READ)==-1 )
//        {
//            pi.status = pf.GetError();
//            return -1;
//        }


//        bool found = true;
//        int i;
//        for(i=0;i<23;++i)
//        {
//            if( pf.AddToRead(CameraDesc(i))==-1 )
//            {
//                found = false;
//                break;
//            }
//        }

//        if(!found)
//            return -1;

//        for(i=0;i<int(pf.elements.size());i++)
//        {
//            int n = pf.ElemNumber(i);

//            if( !strcmp( pf.ElemName(i),"camera" ) )
//            {
//                pf.SetCurElement(i);

//                //LoadPly_Camera ca;

//                for(int j=0;j<n;++j)
//                {
//                    if( pf.Read( (void *)&(ca) )==-1 )
//                    {
//                        pi.status = PlyInfo::E_SHORTFILE;
//                        return -1;
//                    }
//                    camera.valid     = true;
//                    camera.view_p[0] = ca.view_px;
//                    camera.view_p[1] = ca.view_py;
//                    camera.view_p[2] = ca.view_pz;
//                    camera.x_axis[0] = ca.x_axisx;
//                    camera.x_axis[1] = ca.x_axisy;
//                    camera.x_axis[2] = ca.x_axisz;
//                    camera.y_axis[0] = ca.y_axisx;
//                    camera.y_axis[1] = ca.y_axisy;
//                    camera.y_axis[2] = ca.y_axisz;
//                    camera.z_axis[0] = ca.z_axisx;
//                    camera.z_axis[1] = ca.z_axisy;
//                    camera.z_axis[2] = ca.z_axisz;
//                    camera.f         = ca.focal;
//                    camera.s[0]      = ca.scalex;
//                    camera.s[1]      = ca.scaley;
//                    camera.c[0]      = ca.centerx;
//                    camera.c[1]      = ca.centery;
//                    camera.viewport[0] = ca.viewportx;
//                    camera.viewport[1] = ca.viewporty;
//                    camera.k[0]      = ca.k1;
//                    camera.k[1]      = ca.k2;
//                    camera.k[2]      = ca.k3;
//                    camera.k[3]      = ca.k4;
//                }
//                break;
//            }
//        }

//        return 0;
//    }


    static bool LoadMask(const char * filename, int &mask)
    {
        PlyInfo pi;
        return LoadMask(filename, mask, pi);
    }

    static bool LoadMask(const char * filename, int &mask, PlyInfo &pi)
    {
        mask=0;
        vcg::ply::PlyFile pf;
        if( pf.Open(filename,vcg::ply::PlyFile::MODE_READ)==-1 )
        {
            pi.status = pf.GetError();
            return false;
        }

        if( pf.AddToRead(VertDesc(0))!=-1 &&
            pf.AddToRead(VertDesc(1))!=-1 &&
            pf.AddToRead(VertDesc(2))!=-1 )
            mask |= vcg::tetra::io::Mask::IOM_VERTCOORD;

        if( pf.AddToRead(VertDesc(3))!=-1 )
            mask |= vcg::tetra::io::Mask::IOM_VERTFLAGS;

        if( pf.AddToRead(VertDesc(4))!=-1 )
            mask |= vcg::tetra::io::Mask::IOM_VERTQUALITY;

        if( ( pf.AddToRead(VertDesc(5))!=-1 ) &&
            ( pf.AddToRead(VertDesc(6))!=-1 ) &&
            ( pf.AddToRead(VertDesc(7))!=-1 ) &&
            ( pf.AddToRead(VertDesc(8))!=-1 ))
            mask |= vcg::tetra::io::Mask::IOM_VERTCOLOR;

        if( pf.AddToRead(TetraDesc(0))!=-1 )
            mask |= vcg::tetra::io::Mask::IOM_TETRAINDEX;

        if( pf.AddToRead(TetraDesc(1))!=-1 )
            mask |= vcg::tetra::io::Mask::IOM_TETRAFLAGS;

        if( pf.AddToRead(TetraDesc(2))!=-1 )
            mask |= vcg::tetra::io::Mask::IOM_TETRAQUALITY;

        if( ( pf.AddToRead(TetraDesc(6))!=-1 ) &&
            ( pf.AddToRead(TetraDesc(7))!=-1 ) &&
            ( pf.AddToRead(TetraDesc(8))!=-1 ) &&
            ( pf.AddToRead(TetraDesc(9))!=-1 ))
            mask |= vcg::tetra::io::Mask::IOM_TETRACOLOR;


        return true;
    }


}; // end class



} // end namespace tri
} // end namespace io
} // end namespace vcg

#endif
