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

#ifndef __VOLUME_H__
#define __VOLUME_H__

#include "voxel.h"
#include <vcg/space/index/grid_static_ptr.h>

namespace vcg {

// forward definition
template < class VOL >
class VolumeIterator;

//******************************************
//******************************************
//typedef Voxel<float> Voxelf;

const char *SFormat( const char * f, ... )
    {
        static char buf[4096];
        va_list marker;
        va_start( marker, f );
        vsprintf(buf,f,marker);
        va_end( marker );
        return buf;
    }


template<class VOX_TYPE, class SCALAR_TYPE=float>
class Volume {
public:
  typedef SCALAR_TYPE scalar;
  typedef Point3<scalar> Point3x;
  typedef Box3<scalar> Box3x;

  typedef VOX_TYPE voxel_type;

    static int BLOCKSIDE() { return 8;}
    Volume(){
        SetDefaultParam();
    }
    // The actual data
    // They are contained in a block vector.
    std::vector<  std::vector<VOX_TYPE>  > rv;
    Box3x   bbox;

        _int64 AskedCells;
    Point3x dim;  /// Spatial dimension (side length) of the bbox
    Point3i sz;   /// Grid size as the number of cells per side
    Point3i ssz;  /// Sub-block dimensions under consideration as the number of cells per side
    Point3i rsz;  /// Macro grid dimensions of the blocks into which the volume is divided (each block is BLOCKSIDE () ^ 3 cells)
    Point3i asz;  /// Macro dimensions of the block grid relative to the sub-block in question (the one really allocated!)
    Point3x voxel; /// Size of a cell


    int WN,WP; // How many voxels I have to widen to set the manhattan distance in neg and pos
    int DeltaVoxelSafe; // How many voxes I have to expand to be sure in making only a subpart.

  const Point3i ISize() { return sz; }
private :
    // Various offsets and distances pre-calculated once and for all
    Point3f nnf[26],nni[26];
    float len[26],slen[26];

    /// Subpart management
    Point3i div,pos;
public:
    Box3i	  SubPart;                 // Subpart of the volume to be officially considered
    Box3x     SubBox;                  // BBox of the subpart of the volume to be considered in absolute coordinates
    Box3i	  SubPartSafe;             // as above but increased for safety.
    Box3x     SubBoxSafe;

 FILE *LogFP;
 bool Verbose; // if true print a lot of extra info on logfp;

    void SetDefaultParam(){
         WN=0;
         WP=1;
         //WN=-2;//
         //WP=3;
         DeltaVoxelSafe=BLOCKSIDE();
         Verbose=true;
         LogFP=stderr;
        }

    void Init(const Volume &VV)
    {
        SetDefaultParam();
        WN=VV.WN;
        WP=VV.WP;
        DeltaVoxelSafe=VV.DeltaVoxelSafe;
    Init(VV.AskedCells,VV.bbox,VV.div,VV.pos);
    }

        void Init(__int64 cells, Box3x bb, Point3i _div=Point3i(1,1,1), Point3i _pos=Point3i(0,0,0))
    {
        Point3d voxdim;voxdim.Import(bb.max-bb.min);
        AskedCells=cells;
        vcg::BestDim<double>( cells, voxdim, sz );
        bbox=bb;
/*
        printf("grid of ~%i kcells: %d x %d x %d \n",int(cells/1000),sz[0],sz[1],sz[2]);
        printf("grid voxel size of %f %f %f\n",voxdim[0]/sz[0],voxdim[1]/sz[1],voxdim[2]/sz[2]);
*/
        // the box must be multiples of BLOCKSIDE()
        sz=((sz/BLOCKSIDE())+Point3i(1,1,1))*BLOCKSIDE();


        rsz=sz/BLOCKSIDE();
        if(sz!=rsz*BLOCKSIDE()) {
            assert(0); // the box must be multiples of BLOCKSIDE()
            exit(-1);
        }

        dim=bbox.max-bbox.min;
        voxel[0]=dim[0]/sz[0];
        voxel[1]=dim[1]/sz[1];
        voxel[2]=dim[2]/sz[2];

        SetSubPart(_div,_pos);
        ssz=SubPartSafe.max-SubPartSafe.min;
        asz=ssz/BLOCKSIDE() + Point3i(1,1,1);
        rv.clear();
        rv.resize(asz[0]*asz[1]*asz[2]);
        for(size_t i=0;i<rv.size();++i)
            rv[i].resize(0,VOX_TYPE::Zero());
        SetDim(bb);
    }

private:

    // Always to be called AFTER the resize...
    void SetDim(const Box3x & /*bb*/)
    {

    // Setup the precomputed offsets and offset normals
        int cnt=0,x,y,z;
        for(z=-1;z<=1;++z){
         for(y=-1;y<=1;++y){
            for(x=-1;x<=1;++x)
                    if(x!=0 || y!=0 || z!=0)
                    {
                        nnf[cnt]=Point3f(x,y,z);
                        len[cnt]=nnf[cnt].Norm();
                        slen[cnt]=nnf[cnt].SquaredNorm();
                        nnf[cnt].Normalize();
                        nni[cnt]=Point3i(x,y,z);
                        cnt++;
                    }
                }
        }
    }

/*
Parameters
div indicates the number of blocks to be made along the various axes (always> = 1)
pos indicates the coord of the subbox to be taken into consideration (always> = 0 && <xdiv, ydiv, zdiv)
*/

void SetSubPart(Point3i _div, Point3i _pos)
{
    // Parameter correctness check.
    short k;
    for(k=0;k<3;++k)
        {
            assert(_div[k]>0);
            if(_div[k]==0){
                printf("Error in subbox definition:\n the subdiv settings must be greater than 0; it was %i %i %i\n",_div[0],_div[1],_div[2]);
                exit(-1);
            }
            if(_pos[k]<0 || _pos[k]>=_div[k]){
                printf("Error in subbox definition:\n the Position of the subbox must be between (0,0,0) and (%i,%i,%i); it was %i %i %i\n",_div[0],_div[1],_div[2],_pos[0],_pos[1],_pos[2]);
                exit(-1);
            }
        }

    div=_div;
    pos=_pos;

    // Setting the subpart under analisys
    for(k=0;k<3;++k)
        {
            SubPart.min[k]= pos[k]*sz[k]/div[k];
            SubPart.max[k]=(pos[k]+1)*sz[k]/div[k];
            SubBox.min[k]= bbox.min[k]+SubPart.min[k]*voxel[k];
            SubBox.max[k]= bbox.min[k]+SubPart.max[k]*voxel[k];
        }

    // Setting the Safe Subpart under analisys
    SubPartSafe=SubPart;
    for(k=0;k<3;++k)
        {
            SubPartSafe.min[k] -= DeltaVoxelSafe;
            SubPartSafe.max[k] += DeltaVoxelSafe;

            if( SubPartSafe.min[k]< 0     ) SubPartSafe.min[k] = 0;
            if( SubPartSafe.max[k]> sz[k] ) SubPartSafe.max[k] = sz[k];
            SubBoxSafe.min[k]= bbox.min[k]+SubPartSafe.min[k]*voxel[k];
            SubBoxSafe.max[k]= bbox.min[k]+SubPartSafe.max[k]*voxel[k];
        }
/*
        printf("  Computing only subvolume: (%d x %d x %d)= %dk cells  \n"
                 "                             %d,%d,%d -> %d,%d,%d\n"
                        ,SubPart.DimX(),SubPart.DimY(),SubPart.DimZ(),(int)(((__int64)SubPart.DimX()*(__int64)SubPart.DimY()*(__int64)SubPart.DimZ())/1000)
            ,SubPart.min[0]			,SubPart.min[1]			,SubPart.min[2]
            ,SubPart.max[0]			,SubPart.max[1]			,SubPart.max[2]		);
*/
}

public:

    // Sa
    /*bool Write(string filename, const float &minv, const float &maxv )
        {
          FILE *fp;
            if(div!=Point3i(1,1,1)) {
                    string subvoltag;
                    GetSubVolumeTag(subvoltag);
                    filename+=subvoltag;
            }
            string datname=filename;
            string rawname=filename;
            datname+=".dat";
            rawname+=".raw";

          fp=fopen(datname,"w");

            fprintf(fp,"ObjectFileName: %s\n",rawname.c_str());
            fprintf(fp,"TaggedFileName: ---\n");
            fprintf(fp,"Resolution:     %i %i %i\n",SubPart.max[0]-SubPart.min[0],SubPart.max[1]-SubPart.min[1],SubPart.max[2]-SubPart.min[2]);
            fprintf(fp,"SliceThickness: %f %f %f\n",voxel[2]/voxel[0],voxel[2]/voxel[1],voxel[2]/voxel[2]);
            fprintf(fp,"Format:         UCHAR\n");
            fprintf(fp,"NbrTags:        0\n");
            fprintf(fp,"ObjectType:     TEXTURE_VOLUME_OBJECT\n");
            fprintf(fp,"ObjectModel:    RGBA\n");
            fprintf(fp,"GridType:       EQUIDISTANT\n");

            fclose(fp);
      fp=fopen(rawname,"wb");
         if(!fp)
         {
             printf("Error: unable ro open output volume file '%s'\n",filename);
             return false;
         }

         int i,j,k;
         for(k=SubPart.min[2];k<SubPart.max[2];++k)
             for(j=SubPart.min[1];j<SubPart.max[1];++j)
                 for(i=SubPart.min[0];i<SubPart.max[0];++i)
                    {
                        float fv=V(i,j,k).V();
                      fv=(fv-minv)/(maxv-minv);
                        if(fv<0) fv=0;
                        else if(fv>1) fv=1;
                        fv=((fv*2.0f)-1.0f)*127;
                        char fs= (char) fv;
                        fwrite(&fs,sizeof(char),1,fp);
                    }
        fclose(fp);
        return true;
        }*/
    void AbsPos(Point3i pi, Point3x &p)
        {
            p[0]=bbox.min[0]+pi[0]*voxel[0];
            p[1]=bbox.min[1]+pi[1]*voxel[1];
            p[2]=bbox.min[2]+pi[2]*voxel[2];
        }


    void GetSubVolumeTag(std::string &subtag)
    {
    char buf[32];
        if     (div[0]<=  10 && div[1]<=  10 && div[2]<=  10 ) sprintf(buf,"_%01d%01d%01d",pos[0],pos[1],pos[2]);
        else if(div[0]<= 100 && div[1]<= 100 && div[2]<= 100 ) sprintf(buf,"_%02d%02d%02d",pos[0],pos[1],pos[2]);
                                                         else  sprintf(buf,"_%03d%03d%03d",pos[0],pos[1],pos[2]);
        subtag=buf;
    }

    /*
     * Compute the offset <lpos> inside the subblock <rpos> of voxel (x,y,z).    
     * return true if the subblock is allocated.
     */
    bool Pos(const int &_x,const int &_y,const int &_z, int & rpos,int &lpos) const
    {
        int x=_x-SubPartSafe.min[0],y=_y-SubPartSafe.min[1],z=_z-SubPartSafe.min[2];

        assert(_x>=SubPartSafe.min[0] && _x<SubPartSafe.max[0] &&
               _y>=SubPartSafe.min[1] && _y<SubPartSafe.max[1] &&
               _z>=SubPartSafe.min[2] && _z<SubPartSafe.max[2]);

    //	assert(x>=0 && x<sz[0] && y>=0 && y<sz[1] && z>=0 && z<sz[2]);

        int rx=x/BLOCKSIDE(),ry=y/BLOCKSIDE(),rz=z/BLOCKSIDE();
        assert(rx>=0 && rx<asz[0] && ry>=0 && ry<asz[1] && rz>=0 && rz<asz[2]);
        rpos = rz*asz[0]*asz[1]+ry*asz[0]+rx;
        assert(rpos < int(rv.size()));
        int lx = x%BLOCKSIDE(),ly = y%BLOCKSIDE(),lz = z % BLOCKSIDE();
        lpos = lz*BLOCKSIDE()*BLOCKSIDE()+ly*BLOCKSIDE()+lx;
        if(rv[rpos].empty()) return false;
        return true;
     }

    /* 
    Reverse function of the previous one
    Given two positions rpos and lpos returns absolute x, y, z
    */
    bool IPos(int &x,int &y,int &z, const int & rpos, const int &lpos) const
    {
        assert (rpos>=0 && lpos  >=0);

        int rz =   rpos / (asz[0]*asz[1]),remainder =  rpos % (asz[0]*asz[1]);
        int ry = ( remainder ) / asz[0] ;
        int rx =   remainder % asz[0];

        assert(rx>=0 && rx<asz[0] && ry>=0 && ry<asz[1] && rz>=0 && rz<asz[2]);

        int lz =   lpos / (BLOCKSIDE()*BLOCKSIDE()),lemaindel =  lpos % (BLOCKSIDE()*BLOCKSIDE());
        int ly = ( lemaindel ) / BLOCKSIDE();
        int lx =   lemaindel % BLOCKSIDE();

        x = rx*BLOCKSIDE()+lx;
        y = ry*BLOCKSIDE()+ly;
        z = rz*BLOCKSIDE()+lz;

        x+=SubPartSafe.min[0];
        y+=SubPartSafe.min[1];
        z+=SubPartSafe.min[2];

        assert(x>=0 && x<sz[0] && y>=0 && y<sz[1] && z>=0 && z<sz[2]);
        //int trpos,tlpos;
        //assert(rv[rpos].size()>0);
        //Pos(x,y,z,trpos,tlpos);
        //assert(trpos==rpos && tlpos == lpos);
        return true;
     }

    void Alloc(int rpos, const VOX_TYPE &zeroval)
    {
        rv[rpos].resize(BLOCKSIDE()*BLOCKSIDE()*BLOCKSIDE(),zeroval);
    }
    /************************************/
    // Data access functions
  bool ValidCell(const Point3i &p1, const Point3i &p2) const
  {
     if(!cV(p1.X(),p1.Y(),p1.Z()).B() ) return false;
     if(!cV(p2.X(),p1.Y(),p1.Z()).B() ) return false;
     if(!cV(p1.X(),p2.Y(),p1.Z()).B() ) return false;
     if(!cV(p2.X(),p2.Y(),p1.Z()).B() ) return false;
     if(!cV(p1.X(),p1.Y(),p2.Z()).B() ) return false;
     if(!cV(p2.X(),p1.Y(),p2.Z()).B() ) return false;
     if(!cV(p1.X(),p2.Y(),p2.Z()).B() ) return false;
     if(!cV(p2.X(),p2.Y(),p2.Z()).B() ) return false;
    return true;
  }

  float Val(const int &x,const int &y,const int &z) const {
    if(!cV(x,y,z).B())  return 1000;
      return cV(x,y,z).V();
    //else return numeric_limits<float>::quiet_NaN( );
  }

    VOX_TYPE &V(const int &x,const int &y,const int &z) {
        int rpos,lpos;
        if(!Pos(x,y,z,rpos,lpos)) Alloc(rpos,VOX_TYPE::Zero());
        return rv[rpos][lpos];
    }

    const VOX_TYPE &cV(const int &x,const int &y,const int &z) const
    {
        int rpos,lpos;
        if(!Pos(x,y,z,rpos,lpos)) return VOX_TYPE::Zero();
        else return rv[rpos][lpos];
    }
    const VOX_TYPE &V(const int &x,const int &y,const int &z) const
    {
        int rpos,lpos;
        if(!Pos(x,y,z,rpos,lpos)) return VOX_TYPE::Zero();
        else return rv[rpos][lpos];
    }
    /************************************/
    void Fill(VOX_TYPE (__cdecl *func)(const Point3i &p) )
{
    int x,y,z;
    for(z=0;z<sz[2];++z)
        for(y=0;y<sz[1];++y)
            for(x=0;x<sz[0];++x)
                        {
                                Point3i p(x,y,z);
                                V(x,y,z)=func(p);
                        }
}

void Fill(VOX_TYPE const p)
{
    int x,y,z;
    for(z=0;z<sz[2];++z)
        for(y=0;y<sz[1];++y)
            for(x=0;x<sz[0];++x)
                {
                                V(x,y,z)=p;
                }
}


// Copies a smoothed version of the incoming volume to the current volume S.
// the parameter is used to specify the range of field values close to zero that should not be averaged!
// this is because if you smooth even on zero you also smooth where it is well aligned

void CopySmooth( Volume<VOX_TYPE> &S, scalar SafeZone=1, scalar SafeQuality=0)
{
    if(sz!=S.sz)
        {
            printf("Error");
            exit(-1);
        }
 int lcnt=0;
 VolumeIterator< Volume > vi(S);
 vi.Restart();
 vi.FirstNotEmpty();
 // const Voxelf *VC;
 while(vi.IsValid())
     // scan the input volume, for each non-empty voxel of the volume
     // in input it calculates the average with the neighbors
    {
        if((*vi).B())
        {
            int x,y,z;
            IPos(x,y,z,vi.rpos,vi.lpos);
            if(Bound1(x,y,z))
                {
                  VOX_TYPE &VC =  V(x,y,z);
                    for(short int i=0;i<26;++i)
                    {
                        VOX_TYPE &VV= S.V(x+nni[i][0],y+nni[i][1],z+nni[i][2]);
                        if(VV.B()) VC+=VV;
                    }
                    lcnt++;

                /*
                    Voxelf &VV=V(x,y,z);
                    //Voxelf &VV=rv[vi.rpos][vi.lpos];
                    for(int i=0;i<26;++i)
                    {
                        VC = &(S.V(x+nni[i][0],y+nni[i][1],z+nni[i][2]));
                        if(VC->b)
                        {
                            VV+=*VC;
                            lcnt++;
                        }
                    }*/
                }
        }
        vi.Next();
        if(vi.IsValid()) vi.FirstNotEmpty();
        //if((lcnt%100)==0) vi.Dump();

    }
 // Step 2,
 // after calculating the average,

 VolumeIterator< Volume > svi(*this);
 svi.Restart();
 svi.FirstNotEmpty();
 int smoothcnt=0,preservedcnt=0,blendedcnt=0;
 
 const float FieldBorder = 1; // where the transition between the safe and smoothed zone ends
 const float EndFBorderZone = SafeZone+FieldBorder;
 const float EndQBorderZone = SafeQuality*1.5;
 const float QBorder = EndQBorderZone-SafeQuality; // where the transition between the safe and smoothed zone ends
 while(svi.IsValid())
    {
        if((*svi).Cnt()>0)
        {
            VOX_TYPE &sv=S.rv[svi.rpos][svi.lpos];
            (*svi).Normalize(1); // contains the averaged value
            float SafeThr = fabs(sv.V());

            // If the quality is low or if we are distant we always smoothes
            // if we are close to zero and with good quality, we must be careful
            if(SafeThr<EndFBorderZone && sv.Q() > EndQBorderZone)
            {		// if the current voxel had a value <safezone AND quality> SafeQuality
                    // then copy the original value of S
                    if((SafeThr <= SafeZone) && sv.Q() > SafeQuality )
                        {
                            (*svi)=sv;
                            (*svi).SetB(true);
                            ++preservedcnt;
                        }
                        else
                        {	// We are in the transition zone either by field or by quality
                            float blendq= std::max(0.0f,std::min(1.0f,(EndQBorderZone-sv.Q())/QBorder));
                            float blendf= std::max(0.0f,std::min(1.0f,(EndFBorderZone-SafeThr)/FieldBorder));
                            float BlendFactor = 1.0-std::max(blendf,blendq); // how much of the original voxel <sv> you take;
                            (*svi).Blend(sv,BlendFactor);
                            ++blendedcnt;
                        }
            }
            ++smoothcnt;
        }
        svi.Next();
        if(svi.IsValid()) svi.FirstNotEmpty();
    }

 if(Verbose) fprintf(LogFP,"CopySmooth %i voxels, %i preserved, %i blended\n",smoothcnt,preservedcnt,blendedcnt);
}

void Merge(Volume<VOX_TYPE> &S)
{
 VolumeIterator< Volume > svi(S);
 svi.Restart();
 svi.FirstNotEmpty();
 int loccnt=0;

 while(svi.IsValid())
    {
     if((*svi).B())
         {
            int x,y,z;
            IPos(x,y,z,svi.rpos,svi.lpos);
            if(cV(x,y,z).B())	V(x,y,z).Merge( (*svi));
                    else {
                        V(x,y,z).Set((*svi));
                        V(x,y,z).SetB(true);
                    }
            ++loccnt;
         }
    svi.Next();
  if(svi.IsValid()) svi.FirstNotEmpty();
    }

 printf("Merge2 %i voxels\n",loccnt);

}


void Interize( Point3x & vert ) const // OK
{
    for(short j=0;j<3;++j)
    {
        assert(vert[j]>=bbox.min[j]);
        assert(vert[j]<=bbox.max[j]);
        vert[j] = (vert[j] - bbox.min[j]) * sz[j] / (bbox.max[j] - bbox.min[j]);
    }
}


void DeInterize( Point3x & vert ) const	// OK
{
    for(short j=0;j<3;++j)
        vert[j] = vert[j] * (bbox.max[j] - bbox.min[j]) / sz[j] + bbox.min[j];
}

bool SplatVert( const Point3x & v0, double quality, const Point3x & nn, Color4b c)
{
    Box3i ibox;

    assert(math::Abs(SquaredNorm(nn) - 1.0) < 0.0001); // Just a safety check that the vertex normals are NORMALIZED!
    ibox.min=Point3i(floor(v0[0]),floor(v0[1]),floor(v0[2]));
    ibox.max=Point3i( ceil(v0[0]), ceil(v0[1]), ceil(v0[2]));
    ibox.Intersect(SubPartSafe);

    ibox.max[0] = std::min(SubPartSafe.max[0]-1,ibox.max[0]);
    ibox.max[1] = std::min(SubPartSafe.max[1]-1,ibox.max[1]);
    ibox.max[2] = std::min(SubPartSafe.max[2]-1,ibox.max[2]);


     // Skip faces not colliding current subvolume.
    if(ibox.IsEmpty())
        {
            // point outside the box do nothing
            return false;
        }

    Point3x iV, deltaIV;

    // Now scan the eight voxel surrounding the splat
    // and fill them with the distance from the plane
    for(iV[0]=ibox.min[0]; iV[0]<=ibox.max[0]; ++iV[0])
        for(iV[1]=ibox.min[1]; iV[1]<=ibox.max[1]; ++iV[1])
            for(iV[2]=ibox.min[2]; iV[2]<=ibox.max[2]; ++iV[2])
                {
                    deltaIV = v0-iV;
                    scalar offset = deltaIV * nn;
                    VOX_TYPE &VV=V(iV[0],iV[1],iV[2]);
                    VV=VOX_TYPE(offset,nn,quality,c);
                }
        return true;
}

template <const int CoordZ>
        void RasterFace(const int sx, const int ex, const int sy, const int ey,
                        scalar dist, const Point3x &norm, scalar quality,
                        const Point3x &v0,  const Point3x &v1,  const Point3x &v2,
                        const Point3x &d10, const Point3x &d21, const Point3x &d02)
{
  const scalar EPS     = scalar(1e-12);
  const int crd0 = CoordZ,crd1 = (CoordZ+1)%3,crd2 = (CoordZ+2)%3;
  assert(fabs(norm[crd0])+0.001 > fabs(norm[crd1]));
  assert(fabs(norm[crd0])+0.001 > fabs(norm[crd2]));
  scalar x,y;
    for(x=sx;x<=ex;++x)
        for(y=sy;y<=ey;++y)
        {
            scalar n0 = (x-v0[crd1])*d10[crd2] - (y-v0[crd2])*d10[crd1];
            scalar n1 = (x-v1[crd1])*d21[crd2] - (y-v1[crd2])*d21[crd1];
            scalar n2 = (x-v2[crd1])*d02[crd2] - (y-v2[crd2])*d02[crd1];

            if( (n0>-EPS && n1>-EPS && n2>-EPS) ||
                (n0< EPS && n1< EPS && n2< EPS ) )
            {
                scalar iz = ( dist - x*norm[crd1] - y*norm[crd2] ) / norm[crd0];
                //assert(iz>=fbox.min[2] && iz<=fbox.max[2]);
                AddIntercept<CoordZ>(x,y,iz, quality, norm );
            }
        }
}

// The face is known to have an intercept on the z-dir axis of coord xy at position z;
// then set the corresponding distance in the 2 vertices before and 2 after.
template<int CoordZ>
        void AddIntercept( const int x, const int y, const scalar z, const scalar q, const Point3f &n )
{
    scalar esgn = (n[CoordZ] > 0 ? -1 : 1);
    int  zint = floor(z);
    scalar dist=z-zint;  // always positive and between zero and one

    for(int k=WN;k<=WP;k++)
    {
        if(zint+k >= SubPartSafe.min[CoordZ] && zint+k < SubPartSafe.max[CoordZ])
        {
            VOX_TYPE *VV;
            if(CoordZ==2) VV=&V(x,y,zint+k);
            if(CoordZ==1) VV=&V(y,zint+k,x);
            if(CoordZ==0) VV=&V(zint+k,x,y);
            scalar nvv= esgn*( k-dist);
            if(!VV->B() || fabs(VV->V()) > fabs(nvv))		{
                *VV=VOX_TYPE(nvv,n,q);
            }
        }
    }
}

// assumes that the points of the incoming face have been interized
bool ScanFace2( const Point3x & v0, const Point3x & v1, const Point3x & v2,
                       scalar quality, const Point3x & norm)//, const int name )	// OK
{

    Box3x fbox;		// Bounding Box of the face (double)
    Box3i ibox;		// Bounding Box of the face (int)
    int sx,sy,sz;	// Bounding Box whole
    int ex,ey,ez;

    // Bbox calculation of the face
    fbox.Set(v0);
    fbox.Add(v1);
    fbox.Add(v2);

    // BBox integer (note that the cast to int does truncation (it only works because v0, v1, v2 are positive)

    ibox.min[0] =sx = floor(fbox.min[0]); if( ((scalar)sx)!=fbox.min[0] ) ++sx; // necessary if the point is close to .9999
    ibox.min[1] =sy = floor(fbox.min[1]); if( ((scalar)sy)!=fbox.min[1] ) ++sy;
    ibox.min[2] =sz = floor(fbox.min[2]); if( ((scalar)sz)!=fbox.min[2] ) ++sz;
    ibox.max[0] =ex = floor(fbox.max[0]);
    ibox.max[1] =ey = floor(fbox.max[1]);
    ibox.max[2] =ez = floor(fbox.max[2]);
     // Skip faces not colliding current subvolume.
    if(!ibox.Collide(SubPartSafe)) return false;

    Point3x d10 = v1 - v0;
    Point3x d21 = v2 - v1;
    Point3x d02 = v0 - v2;

    assert(norm.Norm() >0.999f && norm.Norm()<1.001f);
//    assert(0);
    scalar  dist = norm * v0;


        /**** Rasterization bbox ****/

    // Clamping of the rasterization values to the current subbox
    sx = std::max(SubPartSafe.min[0],sx); ex = std::min(SubPartSafe.max[0]-1,ex);
    sy = std::max(SubPartSafe.min[1],sy); ey = std::min(SubPartSafe.max[1]-1,ey);
    sz = std::max(SubPartSafe.min[2],sz); ez = std::min(SubPartSafe.max[2]-1,ez);

    if(fabs(norm[0]) > fabs(norm[1]) && fabs(norm[0])>fabs(norm[2])) RasterFace<0>(sy,ey,sz,ez,dist,norm,quality,v0,v1,v2,d10,d21,d02);
    else if(fabs(norm[1]) > fabs(norm[0]) && fabs(norm[1])>fabs(norm[2])) RasterFace<1>(sz,ez,sx,ex,dist,norm,quality,v0,v1,v2,d10,d21,d02);
    else RasterFace<2>(sx,ex,sy,ey,dist,norm,quality,v0,v1,v2,d10,d21,d02);


return true;
}

// assumes that the points of the incoming face have been interized
bool ScanFace( const Point3x & v0, const Point3x & v1, const Point3x & v2,
                       double quality, const Point3x & nn)//, const int name )	// OK
{
    const scalar EPS     = scalar(1e-12);
//	const scalar EPS_INT = scalar(1e-20);
    const scalar EPS_INT = .3f;//scalar(1e-20);

#ifndef _NDEBUG
    if(quality==0.0)
    {
        printf("Zero quality face are not allowed!\n");
        exit(-1);
    }
#endif

//	++nfaces;

    Box3x fbox;		// Bounding Box of the face (double)
    Box3i ibox;		// Bounding Box of the face (int)
    int sx,sy,sz;	// Bounding Box whole
    int ex,ey,ez;

        // Face bbox calculation
    fbox.Set(v0);
    fbox.Add(v1);
    fbox.Add(v2);

    // Integer BBox (note that the cast to int does truncation (it only works because v0, v1, v2 are positive)

    ibox.min[0] =sx = floor(fbox.min[0]); if( ((scalar)sx)!=fbox.min[0] ) ++sx; // necessary if the point is approx a .9999
    ibox.min[1] =sy = floor(fbox.min[1]); if( ((scalar)sy)!=fbox.min[1] ) ++sy;
    ibox.min[2] =sz = floor(fbox.min[2]); if( ((scalar)sz)!=fbox.min[2] ) ++sz;
    ibox.max[0] =ex = floor(fbox.max[0]);
    ibox.max[1] =ey = floor(fbox.max[1]);
    ibox.max[2] =ez = floor(fbox.max[2]);
     // Skip faces not colliding current subvolume.
    if(!ibox.Collide(SubPartSafe)) return false;

        /**** Data for intersection control ****/

        // Versors of the edges of the face

    Point3x d10 = v1 - v0;
    Point3x d21 = v2 - v1;
    Point3x d02 = v0 - v2;

        // Normal to face plane and origin distance

    Point3x norm = d10 ^ d21;
    norm.Normalize();
    double  dist = norm * v0;

        /**** Rasterization bbox ****/

    int x,y,z;

    // Clamping of the rasterization values to the current subbox
    sx = std::max(SubPartSafe.min[0],sx); ex = std::min(SubPartSafe.max[0]-1,ex);
    sy = std::max(SubPartSafe.min[1],sy); ey = std::min(SubPartSafe.max[1]-1,ey);
    sz = std::max(SubPartSafe.min[2],sz); ez = std::min(SubPartSafe.max[2]-1,ez);

        // Rasterization xy

    if(fabs(norm[2])>EPS_INT)
    for(x=sx;x<=ex;++x)
        for(y=sy;y<=ey;++y)
        {
            double n0 = ((double)x-v0[0])*d10[1] - ((double)y-v0[1])*d10[0];
            double n1 = ((double)x-v1[0])*d21[1] - ((double)y-v1[1])*d21[0];
            double n2 = ((double)x-v2[0])*d02[1] - ((double)y-v2[1])*d02[0];

            if( (n0>-EPS && n1>-EPS && n2>-EPS) ||
                (n0< EPS && n1< EPS && n2< EPS ))
            {
                double iz = ( dist - double(x)*norm[0] - double(y)*norm[1] ) / norm[2];
                //assert(iz>=fbox.min[2] && iz<=fbox.max[2]);
                AddXYInt(x,y,iz,-norm[2], quality, nn );
            }
        }

        // Rasterization xz

    if(fabs(norm[1])>EPS_INT)
    for(x=sx;x<=ex;++x)
        for(z=sz;z<=ez;++z)
        {
            double n0 = ((double)x-v0[0])*d10[2] - ((double)z-v0[2])*d10[0];
            double n1 = ((double)x-v1[0])*d21[2] - ((double)z-v1[2])*d21[0];
            double n2 = ((double)x-v2[0])*d02[2] - ((double)z-v2[2])*d02[0];

            if( (n0>-EPS && n1>-EPS && n2>-EPS) ||
                (n0< EPS && n1< EPS && n2< EPS ))
            {
                double iy = ( dist - double(x)*norm[0] - double(z)*norm[2] ) / norm[1];
                //assert(iy>=fbox.min[1] && iy<=fbox.max[1]);
                AddXZInt(x,z,iy,-norm[1], quality,nn  );
            }
        }

            // Rasterization yz

    if(fabs(norm[0])>EPS_INT)
    for(y=sy;y<=ey;++y)
        for(z=sz;z<=ez;++z)
        {
            double n0 = ((double)y-v0[1])*d10[2] - ((double)z-v0[2])*d10[1];
            double n1 = ((double)y-v1[1])*d21[2] - ((double)z-v1[2])*d21[1];
            double n2 = ((double)y-v2[1])*d02[2] - ((double)z-v2[2])*d02[1];

            if( (n0>-EPS && n1>-EPS && n2>-EPS) ||
                (n0< EPS && n1< EPS && n2< EPS ) )
            {
                double ix = ( dist - double(y)*norm[1] - double(z)*norm[2] ) / norm[0];
                //assert(ix>=fbox.min[0] && ix<=fbox.max[0]);
                AddYZInt(y,z,ix,-norm[0], quality, nn );
            }
        }
        return true;
}
// The face is known to have an intercept on the z-dir axis of coord xy at position z;
// then set the corresponding distance in the 2 vertices before and 2 after.

void AddXYInt( const int x, const int y, const double z, const double sgn, const double q, const Point3f &n )
{ double esgn = (sgn<0 ? -1 : 1);//*max(fabs(sgn),0.001);
    double dist=z-floor(z);  // always positive and between zero and one
    int  zint = floor(z);
    for(int k=WN;k<=WP;k++)
        if(zint+k >= SubPartSafe.min[2] && zint+k < SubPartSafe.max[2])
        {
            VOX_TYPE &VV=V(x,y,zint+k);
            double nvv= esgn*( k-dist);
            if(!VV.B() || fabs(VV.V()) > fabs(nvv))		{
                    VV=VOX_TYPE(nvv,n,q);
            }
        }
}
void AddYZInt( const int y, const int z, const double x, const double sgn, const double q, const Point3f &n  )
{ double esgn = (sgn<0 ? -1 : 1);//*max(fabs(sgn),0.001);
    double dist=x-floor(x);  // always positive and between zero and one
    int  xint = int(floor(x));
    for(int k=WN;k<=WP;k++)
        if(xint+k >= SubPartSafe.min[0] && xint+k < SubPartSafe.max[0])
        {
            VOX_TYPE &VV=V(xint+k,y,z);
            double nvv= esgn*( k-dist);
            if(!VV.B() || fabs(VV.V()) > fabs(nvv)) {
                    VV=VOX_TYPE(nvv,n,q);
            }
        }
}
void AddXZInt( const int x, const int z, const double y, const double sgn, const double q, const Point3f &n  )
{ double esgn = (sgn<0 ? -1 : 1);//*max(fabs(sgn),0.001);
    double dist=y-scalar(floor(y));  // always positive and between zero and one
    int  yint = floor(y);
    for(int k=WN;k<=WP;k++)
        if(yint+k >= SubPartSafe.min[1] && yint+k < SubPartSafe.max[1])
        {
            VOX_TYPE &VV=V(x,yint+k,z);
            double nvv= esgn*( k-dist);
            if(!VV.B() || fabs(VV.V()) > fabs(nvv))	{
                    VV=VOX_TYPE(nvv,n,q);
            }
        }
}

 void Dump(FILE *fp)
 {
     fprintf(fp,"Volume Info:\n");
     fprintf(fp,"  BBbox %7.4f %7.4f %7.4f - %7.4f %7.4f %7.4f:\n",bbox.min[0],bbox.min[1],bbox.min[2],bbox.max[0],bbox.max[1],bbox.max[2]);
     fprintf(fp,"  Size in voxels    %7i %7i %7i (tot: %7.3f M):\n",sz[0],sz[1],sz[2],(double(sz[0]*sz[1])/1000000.0)*sz[2]);
     fprintf(fp,"  Voxel dimension   %7.4f %7.4f %7.4f \n",voxel[0],voxel[1],voxel[2]);

     fprintf(fp,"  Size in MacroCell %7i %7i %7i (tot: %7.3f M):\n",rsz[0],rsz[1],rsz[2],double(rsz[0]*rsz[1]*rsz[2])/1000000.0);
     fprintf(fp," Memory Info: \n   Voxel Size %8li b Virtually needed mem %8i Mb\n",
                                        sizeof(VOX_TYPE),int(sizeof(VOX_TYPE)*(_int64)(sz[0])*(_int64)(sz[1])*(_int64)(sz[2])/(1024*1024)));
   if(div!=Point3i(1,1,1))
         {
            fprintf(fp,"  Subdivided in      %6i %6i %6i  (tot: %12i block):\n",div[0],div[1],div[2],div[0]*div[1]*div[2]);
            fprintf(fp,"  Computing subblock %6i %6i %6i :\n",pos[0],pos[1],pos[2]);
            fprintf(fp,"                %6i %6i %6i - %6i %6i %6i :\n",SubPart.min[0],SubPart.min[1],SubPart.min[2],SubPart.max[0],SubPart.max[1],SubPart.max[2]);
            fprintf(fp,"        Safe    %6i %6i %6i - %6i %6i %6i :\n",SubPartSafe.min[0],SubPartSafe.min[1],SubPartSafe.min[2],SubPartSafe.max[0],SubPartSafe.max[1],SubPartSafe.max[2]);

         }
     fprintf(fp,"\n");
 }

    int Allocated()
    {int cnt=0;
        for(size_t i=0;i<rv.size();++i)
            if(!rv[i].empty()) cnt++;
            return cnt;
    }

bool Bound1(const int x, const int y, const int z)
{
    return	(x>SubPartSafe.min[0] && x < SubPartSafe.max[0]-1 ) &&
                    (y>SubPartSafe.min[1] && y < SubPartSafe.max[1]-1 ) &&
                    (z>SubPartSafe.min[2] && z < SubPartSafe.max[2]-1 ) ;
}

/*
Notes on the expansion algorithm:

It fills the empty voxels

If the volume is initially filled with the values of the surface intercepts
in the 2 vertices immediately adjacent to the edge intersected by the surface
you have to expand this component in a sensible way.

Note that it is important that not the whole field is filled with approx. Hausdorf distance:
the field must be "cut" on the edges of the surface to avoid mess. Levoy fills the field
only along the direction of the scanner, I instead fill along the normal to the surface. In this way
the problem that the expansion is linked to the initial acquisition is avoided

*/


void Expand(scalar AngleThrRad)
{
 int i;
 VolumeIterator< Volume > vi(*this);

 float CosThr=math::Cos(AngleThrRad);
// printf("Expand2 angle %f, %f\n",AngleThrRad,CosThr);
 int loccnt=0;

 vi.Restart();
 vi.FirstNotEmpty();
 while(vi.IsValid())
 {
     if((*vi).B()) // expands only voxels with "valid" values
        {
            int x,y,z;
            IPos(x,y,z,vi.rpos,vi.lpos);
            Point3f n=(*vi).N();
            VOX_TYPE vtmp =  (*vi);
            if(Bound1(x,y,z))
            for(i=0;i<26;++i)
                        {
                            float angle = -(nnf[i]*n);    // cos angle between surface normal and expansion direction
                            if( fabs(angle)> CosThr )
                                        {
                                            //bbfloat tt=(*vi).V();
                                            vtmp.SetV((*vi).V()+len[i]*angle);   // the new distance is the distance from the plane passing through the point to which VV refers;
                                            VOX_TYPE &VV= V(x+nni[i][0],y+nni[i][1],z+nni[i][2]);
                                            if(!VV.B()){
                                                VV+=vtmp;
                                              loccnt++;
                                            }
                                        }
                        }
     }
    vi.Next();
  if(vi.IsValid()) vi.FirstNotEmpty();
 }
 printf("Expand  %8i ",loccnt);
 Normalize(1);
}


// filter the empty holes of a volume;
// iterates all the filled voxels and adds this value to the adjacent empty ones;
// check all empty voxels and if at least thr kept values have been added.
// basically the higher the value of thr the fewer holes you fill.

void Refill(const int thr,float maxdistance = std::numeric_limits<float>::max() )
{
 int lcnt=0;
 VolumeIterator< Volume > vi(*this);
 vi.Restart();
 vi.FirstNotEmpty();

 while(vi.IsValid())
    {
        if((*vi).B())
        {
            int x,y,z;
            IPos(x,y,z,vi.rpos,vi.lpos);
            if(Bound1(x,y,z))
                {
                    for(short i=0;i<26;++i)
                    {
                        VOX_TYPE &VC= V(x+nni[i][0],y+nni[i][1],z+nni[i][2]);
                        if(!VC.B()){
                            if(VC.Cnt()==0) lcnt++;
                            VC+=(*vi);
                        }
                    }
                }
        }

        vi.Next();

        if(vi.IsValid()) vi.FirstNotEmpty();


    }
 printf("ReFill  %8i ",lcnt);
 Normalize(thr,maxdistance);
}

/*
Traverse the volume and change the field value to create an offset surface that connects
well with the original surface
the parameter specifies where the offset surface should pass

The idea is to make another zero of the distance field at the specified threshold

The thing must be smooth
so I choose a function that has 2 zeros (in zero and in thr) and
*/
void Offset(const float thr)
{
 int lcnt=0;
 VolumeIterator< Volume > vi(*this);
 vi.Restart();
 vi.FirstNotEmpty();
 float thr2=thr/2.0;
 while(vi.IsValid())
    {
        if((*vi).B())
        {
            float vv=(*vi).V();
            if(thr<0) if(vv<thr2) vv=thr-vv;
            if(thr>0) if(vv>thr2) vv=thr-vv;

                (*vi).SetV(vv);
        }

        vi.Next();

        if(vi.IsValid()) vi.FirstNotEmpty();


    }
 printf("ReFill  %8i ",lcnt);
 Normalize(thr);
}

// takes a volume and sets field b of all voxels that have a significant value to true.
// thr indicates the minimum number of values that must be added on the voxel;
// usually equals 1;
int  Normalize(int thr, float maxdistance=std::numeric_limits<float>::max() )
{
 VolumeIterator< Volume > vi(*this);
 vi.Restart();
 vi.FirstNotEmpty();
 int loccnt=0;
 while(vi.IsValid())
    {
     if(!(*vi).B())
        {
          if((*vi).Normalize(thr))
                    ++loccnt;
      if(math::Abs((*vi).V())>maxdistance) *vi=VOX_TYPE::Zero();
         }
    vi.Next();
  if(vi.IsValid()) vi.FirstNotEmpty();
    }
 printf("Normalize(%i) %8i voxels\n",thr,loccnt);
 return loccnt;
}



// Save
void SlicedPPMQ( const char * filename,const char *tag,int SliceNum)
    {
        std::string basename=filename,name;
        int ix,iy,iz;
        Color4b Tab[100];
        for(short ii=1;ii<100;++ii)
            Tab[ii].SetColorRamp(0,100,ii);
        Tab[0]=Color4b::Gray;

    int ZStep=sz[2]/(SliceNum+1);
        for(iz=ZStep;iz<sz[2];iz+=ZStep)
        if(iz>=SubPartSafe.min[2] && iz<SubPartSafe.max[2])
        {
            name=SFormat("%s%03i%s_q.ppm",filename,iz,tag);
            FILE * fp = fopen(name.c_str(),"wb");
            fprintf(fp,
                "P6\n"
                "%d %d\n"
                "255\n"
                ,sz[1]
                ,sz[0]
            );
            unsigned char rgb[3];
            for(ix=0;ix<sz[0];++ix)
            {
                for(iy=0;iy<sz[1];++iy)
                {
                    if(	ix>=SubPartSafe.min[0] && ix<SubPartSafe.max[0] &&
                            iy>=SubPartSafe.min[1] && iy<SubPartSafe.max[1])
                        {
                            if(!V(ix,iy,iz).B())	{
                                rgb[0]=rgb[1]=rgb[2]=64;
                            }
                            else
                            {
                                float vv=V(ix,iy,iz).Q();
                                int qi=std::min(V(ix,iy,iz).Q()*100.0f,100.0f);

                                if( vv>0)		{
                                    rgb[0]=Tab[qi][0];
                                    rgb[1]=Tab[qi][1];
                                    rgb[2]=Tab[qi][2];
                                }
                                else if(vv<0)
                                {
                                    rgb[0]=128;
                                    rgb[1]=255+32*vv;
                                    rgb[2]=0;//V(ix,iy,iz).Q()*256;
                                }
                                else  	{
                                    rgb[0]=255;	rgb[1]=255; rgb[2]=255;
                                }
                            }
                    }
                    else{
                        rgb[0]=rgb[1]=rgb[2]=64;
                    }
                    fwrite(rgb,3,1,fp);
                }
            }
            fclose(fp);
        }
    }

void SlicedPPM( const char * filename,const char *tag,int SliceNum=1)
    {
        std::string basename=filename,name;
        int ix,iy,iz;
        int ZStep=sz[2]/(SliceNum+1);
        for(iz=ZStep;iz<sz[2];iz+=ZStep)
        if(iz>=SubPartSafe.min[2] && iz<SubPartSafe.max[2])
        {
            name=SFormat("%s_%03i_%s.ppm",filename,iz,tag);
      printf("Saving slice '%s'",name.c_str());
            FILE * fp = fopen(name.c_str(),"wb");
            if(!fp) return;
            fprintf(fp,
                "P6\n"
                "%d %d\n"
                "255\n"
                ,sz[1]
                ,sz[0]
            );
            unsigned char rgb[3];
            for(ix=0;ix<sz[0];++ix)
            {
                for(iy=0;iy<sz[1];++iy)
                {
                    if(	ix>=SubPartSafe.min[0] && ix<SubPartSafe.max[0] &&
                            iy>=SubPartSafe.min[1] && iy<SubPartSafe.max[1])
                        {
                            if(!V(ix,iy,iz).B())	{
                                rgb[0]=rgb[1]=rgb[2]=64;
                            }
                            else
                            {
                                float vv=V(ix,iy,iz).V();
                                if( vv>0)		{
                                    rgb[0]=255-32*vv;
                                    rgb[1]=128;
                                    rgb[2]=0;//V(ix,iy,iz).Q()*256;
                                }
                                else if(vv<0)
                                {
                                    rgb[0]=128;
                                    rgb[1]=255+32*vv;
                                    rgb[2]=0;//V(ix,iy,iz).Q()*256;
                                }
                                else  	{
                                    rgb[0]=255;	rgb[1]=255; rgb[2]=255;
                                }
                            }
                    }
                    else{
                        rgb[0]=rgb[1]=rgb[2]=64;
                    }
                    fwrite(rgb,3,1,fp);
                }
            }
            fclose(fp);
        }
    }
template < class VertexPointerType >
void GetXIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, float /*thr*/)
{
  float f1 = Val(p1.X(), p1.Y(), p1.Z());
  float f2 = Val(p2.X(), p2.Y(), p2.Z());
  float u = (float) f1/(f1-f2);
  v->P().X() = (float) p1.X()*(1-u) + u*p2.X();
  v->P().Y() = (float) p1.Y();
  v->P().Z() = (float) p1.Z();
  v->Q()=cV(p1.X(), p1.Y(), p1.Z()).Q();
  v->C()=cV(p1.X(), p1.Y(), p1.Z()).C4b();
}

template < class VertexPointerType >
void GetYIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, float /*thr*/)
{
  float f1 = Val(p1.X(), p1.Y(), p1.Z());
  float f2 = Val(p2.X(), p2.Y(), p2.Z());
  float u = (float) f1/(f1-f2);
  v->P().X() = (float) p1.X();
  v->P().Y() = (float) p1.Y()*(1-u) + u*p2.Y();
  v->P().Z() = (float) p1.Z();
  v->Q()=cV(p1.X(), p1.Y(), p1.Z()).Q();
  v->C()=cV(p1.X(), p1.Y(), p1.Z()).C4b();
}

template < class VertexPointerType>
void GetZIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, float /*thr*/)
{
  float f1 = Val(p1.X(), p1.Y(), p1.Z());
  float f2 = Val(p2.X(), p2.Y(), p2.Z());
  float u = (float) f1/(f1-f2);
  v->P().X() = (float) p1.X();
  v->P().Y() = (float) p1.Y();
  v->P().Z() = (float) p1.Z()*(1-u) + u*p2.Z();
  v->Q()=cV(p1.X(), p1.Y(), p1.Z()).Q();
  v->C()=cV(p1.X(), p1.Y(), p1.Z()).C4b();
}
};



template < class VOL >
class VolumeIterator
{
 public:
     VOL &V;
   //vector<VOL::voxel_type> vi;
     VolumeIterator(VOL &_VV):V(_VV) {}

     //Point3i curPos;
    int rpos;
    int lpos;
    void Restart(){rpos=0;lpos=0;}
 private:
 public:


    void Set(const Point3i &p)
        {
         //curPos=p;
         V.Pos(p[0],p[1],p[2],rpos,lpos);
        }
  bool FirstNotEmpty()
    {

        //Dump();
        typename std::vector<std::vector<typename VOL::voxel_type> >::iterator rvi=V.rv.begin()+rpos;
        do
        {
            if((*rvi).empty())
            {
                while(rvi!=V.rv.end() && (*rvi).empty()) ++rvi;
                if(rvi==V.rv.end())
                {
                    rpos=-1;
                    return false;
                }
                rpos= rvi-V.rv.begin();
                lpos=0;
            }
            typename std::vector<typename VOL::voxel_type>::iterator lvi= (*rvi).begin()+lpos;
            // a voxel is non-empty if it has b! = 0;
            while(lvi!=(*rvi).end() && !((*lvi).B() || (*lvi).Cnt()>0)) {
                ++lvi;
            }
            if(lvi!=(*rvi).end())
            {
                lpos= lvi-(*rvi).begin();
                //V.IPos(p[0],p[1],p[2],rpos,lpos);
                //Dump();
                return true;
            }
            else lpos=0;
            ++rvi;
            rpos= rvi-V.rv.begin();

        } while (rvi!=V.rv.end());
        rpos=-1;
        return false;
    }

    typename VOL::voxel_type &operator *()
        {
          assert(rpos>=0 && lpos >=0);
            return V.rv[rpos][lpos];
        }
    bool Next()
    {
        assert(IsValid());
        if(lpos< VOL::BLOCKSIDE() * VOL::BLOCKSIDE() * VOL::BLOCKSIDE() -1)
        {
            ++lpos;
            //V.IPos(p[0],p[1],p[2],rpos,lpos);
            return true;
        }
        if(rpos < int(V.rv.size()-1))
        {
            lpos=0;
            ++rpos;
            //V.IPos(p[0],p[1],p[2],rpos,lpos);
            return true;
        }
        rpos=-1;
        lpos=-1;
        return false;
    }

    bool IsValid() {
        return rpos>=0;
    }
  void Dump()
        {
        int x,y,z;
        V.IPos(x,y,z,rpos,lpos);
        printf("Iterator r %4i l %4i (%3i %3i %3i)\n",rpos,lpos,x,y,z);
    }
};
}
#endif
