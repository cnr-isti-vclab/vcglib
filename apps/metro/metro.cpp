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
Revision 1.6  2004/07/15 00:15:16  cignoni
inflate -> offset

Revision 1.5  2004/06/24 09:08:31  cignoni
Official Release of Metro 4.00

Revision 1.4  2004/05/14 13:53:12  ganovelli
GPL  added

****************************************************************************/

// -----------------------------------------------------------------------------------------------

// standard libraries
#include <time.h>
#include <locale>

// project definitions.
#include "defs.h"
#include "sampling.h"
#include "mesh_type.h"
#include <vcg/complex/trimesh/update/edges.h>
#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/math/histogram.h>
#include <vcg/complex/trimesh/clean.h>
//#include <wrap/io_trimesh/import_smf.h>
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/import_stl.h>
#include <wrap/io_trimesh/export_ply.h>


// -----------------------------------------------------------------------------------------------


////////////////// Command line Flags and parameters 
bool NumberOfSamples                = false;
bool SamplesPerAreaUnit             = false;
bool CleaningFlag=false;
// -----------------------------------------------------------------------------------------------

void Usage()
{
  printf("\nUsage:  "\
                                        "metro file1 file2  [opt]\n"\
                                        "where opt can be:\n"\
                                        "-v         disable vertex sampling\n"\
                                        "-e         disable edge sampling\n"\
                                        "-f         disable face sampling\n"\
                                        "-u         ignore unreferred vertices\n"\
                                        "-Sx        set the face sampling mode\n"\
                                        "           where x can be:\n"\
                                        "            -S0  montecarlo sampling\n"\
                                        "            -S1  subdivision sampling\n"\
                                        "            -S2  similar triangles sampling (Default)\n"\
                                        "-n#        set the required number of samples (overrides -A)\n"\
                                        "-a#        set the required number of samples per area unit (overrides -N)\n"\
                                        "-c         save error as vertex colour and quality"\
                                        "-C # #     Set the min/max values used for color mapping"\
                                        "-L         Remove duplicated and unreferenced vertices before processing"\
                                        "\n"
                                        "Default options are to sample vertexes, edge and faces, to take a number of sample that is \n"
                                        );
  exit(-1);
}


// simple aux function that compute the name for the file containing the stored computations
string SaveFileName(const string &filename)
{
 int pos=filename.find_last_of('.',filename.length());
 string fileout=filename.substr(0,pos)+"_metro.ply";
 return fileout;
}


// simple aux function that returns true if a given file has a given extesnion
bool FileExtension(string filename,  string extension)
{
  locale loc1 ;
  use_facet<ctype<char> > ( loc1 ).tolower(&*filename.begin(),&*filename.end());
  use_facet<ctype<char> > ( loc1 ).tolower(&*extension.begin(),&*extension.end());
  string end=filename.substr(filename.length()-extension.length(),extension.length());
  return end==extension;
}

// Open Mesh
void OpenMesh(const char *filename, CMesh &m)
{
  int err;
  if(FileExtension(filename,"ply"))
  {
    err = tri::io::ImporterPLY<CMesh>::Open(m,filename);
    if(err) {
      printf("Error in reading %s: '%s'\n",filename,tri::io::ImporterPLY<CMesh>::ErrorMsg(err));
      exit(-1);
    }
    printf("read mesh `%s'\n", filename);		  
  }
  else if(FileExtension(filename,"stl"))
  {
    err = tri::io::ImporterSTL<CMesh>::Open(m,filename);
    if(err) {
      printf("Error in reading %s: '%s'\n",filename,tri::io::ImporterSTL<CMesh>::ErrorMsg(err));
      exit(-1);
    }
    printf("read mesh `%s'\n", filename);		  
  }
  else {
    printf("Unknown file format for mesh '%s'\n",filename);
    exit(-1);
  }
  if(CleaningFlag){
  int dup = tri::Clean<CMesh>::RemoveDuplicateVertex(m);
  int unref =  tri::Clean<CMesh>::RemoveUnreferencedVertex(m);
  printf("Removed %i duplicate and %i unreferenced vertices from mesh %s\n",dup,unref,filename);
  }
}


int main(int argc, char**argv)
{
    CMesh                 S1, S2;
    float                ColorMin=0, ColorMax=0;
    double                dist1_max, dist2_max;
    unsigned long         n_samples_target, elapsed_time;
    double								n_samples_per_area_unit;
    int                   flags;
 
    // print program info
    printf("-------------------------------\n"
           "         Metro V.4.0 \n"
           "     http://vcg.isti.cnr.it\n"
           "   release date: "__DATE__"\n"
           "-------------------------------\n\n");

    if(argc <= 2)    Usage();
    // default parameters
    flags = SamplingFlags::VERTEX_SAMPLING |
          SamplingFlags::EDGE_SAMPLING |
          SamplingFlags::FACE_SAMPLING |
          SamplingFlags::SIMILAR_SAMPLING;

    // parse command line.
	  for(int i=3; i < argc;)
    {
      if(argv[i][0]=='-')
        switch(argv[i][1])
      { 
        case 'v' : flags &= ~SamplingFlags::VERTEX_SAMPLING; break;
        case 'e' : flags &= ~SamplingFlags::EDGE_SAMPLING; break;
        case 'f' : flags &= ~SamplingFlags::FACE_SAMPLING; break;
        case 'u' : flags |= SamplingFlags::INCLUDE_UNREFERENCED_VERTICES; break;
        case 's'   :
          switch(argv[i][2])
          {
            case 0:  flags = (flags | SamplingFlags::MONTECARLO_SAMPLING  ) & (~ SamplingFlags::NO_SAMPLING );break;
            case 1:  flags = (flags | SamplingFlags::SUBDIVISION_SAMPLING ) & (~ SamplingFlags::NO_SAMPLING );break;
            case 2:  flags = (flags | SamplingFlags::SIMILAR_SAMPLING     ) & (~ SamplingFlags::NO_SAMPLING );break;
            default  :  printf(MSG_ERR_INVALID_OPTION, argv[i]);
              exit(0);
          }
          break;
        case 'n':  NumberOfSamples       = true;     n_samples_target        = (unsigned long) atoi(&(argv[i][2]));          break;
        case 'a':  SamplesPerAreaUnit    = true;     n_samples_per_area_unit = (unsigned long) atoi(&(argv[i][2])); break;
        case 'c':  flags |= SamplingFlags::SAVE_ERROR;   break;
        case 'L':  CleaningFlag=true; break;
        case 'C':  ColorMin=atof(argv[i+1]); ColorMax=atof(argv[i+2]); i+=2; break;
        default  :  printf(MSG_ERR_INVALID_OPTION, argv[i]);
          exit(0);
      }
      i++;
    }
 
    // load input meshes.
    OpenMesh(argv[1],S1);
    OpenMesh(argv[2],S2);
    string S1NewName=SaveFileName(argv[1]);
    string S2NewName=SaveFileName(argv[2]);

   
    if(!NumberOfSamples && !SamplesPerAreaUnit)
    {
        NumberOfSamples = true;
        n_samples_target = 10 * max(S1.fn,S2.fn);// take 10 samples per face
    }

    // compute face information
		tri::UpdateEdges<CMesh>::Set(S1);
		tri::UpdateEdges<CMesh>::Set(S2);

	// set bounding boxes for S1 and S2
		tri::UpdateBounding<CMesh>::Box(S1);
		tri::UpdateBounding<CMesh>::Box(S2);

    // set Bounding Box.
    Box3d    bbox, tmp_bbox_M1=S1.bbox, tmp_bbox_M2=S2.bbox;
    bbox.Add(S1.bbox);
    bbox.Add(S2.bbox);
		bbox.Offset(bbox.Diag()*0.02);
	  S1.bbox = bbox;
	  S2.bbox = bbox;
    
    // initialize time info.
    int t0=clock();

    Sampling<CMesh> ForwardSampling(S1,S2);
    Sampling<CMesh> BackwardSampling(S2,S1);

    // print mesh info.
    printf("Mesh info:\n");
    printf(" M1: '%s'\n\t%vertices  %7i\n\tfaces     %7i\n\tarea      %12.4f\n", argv[1], S1.vn, S1.fn, ForwardSampling.GetArea());
    printf("\tbbox (%7.4f %7.4f %7.4f)-(%7.4f %7.4f %7.4f)\n", tmp_bbox_M1.min[0], tmp_bbox_M1.min[1], tmp_bbox_M1.min[2], tmp_bbox_M1.max[0], tmp_bbox_M1.max[1], tmp_bbox_M1.max[2]);
    printf("\tbbox diagonal %f\n", (float)tmp_bbox_M1.Diag());
    printf(" M2: '%s'\n\t%vertices  %7i\n\tfaces     %7i\n\tarea      %12.4f\n", argv[2], S2.vn, S2.fn, BackwardSampling.GetArea());
    printf("\tbbox (%7.4f %7.4f %7.4f)-(%7.4f %7.4f %7.4f)\n", tmp_bbox_M2.min[0], tmp_bbox_M2.min[1], tmp_bbox_M2.min[2], tmp_bbox_M2.max[0], tmp_bbox_M2.max[1], tmp_bbox_M2.max[2]);
    printf("\tbbox diagonal %f\n", (float)tmp_bbox_M2.Diag());

    // Forward distance.
    printf("\nForward distance (M1 -> M2):\n");
    ForwardSampling.SetFlags(flags);
    if(NumberOfSamples)
    {
        ForwardSampling.SetSamplesTarget(n_samples_target);
        n_samples_per_area_unit = ForwardSampling.GetNSamplesPerAreaUnit();
    }
    else
    {
        ForwardSampling.SetSamplesPerAreaUnit(n_samples_per_area_unit);
        n_samples_target = ForwardSampling.GetNSamplesTarget();
    }
    printf("target # samples      : %u\ntarget # samples/area : %f\n", n_samples_target, n_samples_per_area_unit);
    ForwardSampling.Hausdorff();
    dist1_max  = ForwardSampling.GetDistMax();
    printf("\ndistance:\n  max  : %f (%f  with respect to bounding box diagonal)\n", (float)dist1_max, (float)dist1_max/bbox.Diag());
    printf("  mean : %f\n", ForwardSampling.GetDistMean());
    printf("  RMS  : %f\n", ForwardSampling.GetDistRMS());
    printf("# vertex samples %9d\n", ForwardSampling.GetNVertexSamples());
    printf("# edge samples   %9d\n", ForwardSampling.GetNEdgeSamples());
    printf("# area samples   %9d\n", ForwardSampling.GetNAreaSamples());
    printf("# total samples  %9d\n", ForwardSampling.GetNSamples());
    printf("# samples per area unit: %f\n\n", ForwardSampling.GetNSamplesPerAreaUnit());

    // Backward distance.
    printf("\nBackward distance (M2 -> M1):\n");
    BackwardSampling.SetFlags(flags);
    if(NumberOfSamples)
    {
        BackwardSampling.SetSamplesTarget(n_samples_target);
        n_samples_per_area_unit = BackwardSampling.GetNSamplesPerAreaUnit();
    }
    else
    {
        BackwardSampling.SetSamplesPerAreaUnit(n_samples_per_area_unit);
        n_samples_target = BackwardSampling.GetNSamplesTarget();
    }
    printf("target # samples      : %u\ntarget # samples/area : %f\n", n_samples_target, n_samples_per_area_unit);
    BackwardSampling.Hausdorff();
    dist2_max  = BackwardSampling.GetDistMax();
    printf("\ndistance:\n  max  : %f (%f  with respect to bounding box diagonal)\n", (float)dist1_max, (float)dist1_max/bbox.Diag());
    printf("mean : %f\n", BackwardSampling.GetDistMean());
    printf("RMS  : %f\n", BackwardSampling.GetDistRMS());
    printf("# vertex samples %9d\n", BackwardSampling.GetNVertexSamples());
    printf("# edge samples   %9d\n", BackwardSampling.GetNEdgeSamples());
    printf("# area samples   %9d\n", BackwardSampling.GetNAreaSamples());
    printf("# total samples  %9d\n", BackwardSampling.GetNSamples());
    printf("# samples per area unit: %f\n\n", BackwardSampling.GetNSamplesPerAreaUnit());

    // compute time info.
    elapsed_time = clock() - t0;
    int n_total_sample=ForwardSampling.GetNSamples()+BackwardSampling.GetNSamples();
    double mesh_dist_max  = max(dist1_max , dist2_max);
    
    printf("\nHausdorff distance: %f (%f  with respect to bounding box diagonal)\n",(float)mesh_dist_max,(float)mesh_dist_max/bbox.Diag());
    printf("  Computation time  : %d ms\n", (int)elapsed_time);
    printf("  # samples/second  : %f\n\n", (float)n_total_sample/((float)elapsed_time/1000.0));

    // save error files.
    if(flags & SamplingFlags::SAVE_ERROR)
    {
      vcg::tri::io::PlyInfo p; 
      p.mask|=vcg::ply::PLYMask::PM_VERTCOLOR|vcg::ply::PLYMask::PM_VERTQUALITY;
      if(ColorMax!=0 || ColorMin != 0){
        vcg::tri::UpdateColor<CMesh>::VertexQuality(S1,ColorMin,ColorMax);
        vcg::tri::UpdateColor<CMesh>::VertexQuality(S2,ColorMin,ColorMax);
      }
      tri::io::ExporterPLY<CMesh>::Save( S1,S1NewName.c_str(),true,p);
      tri::io::ExporterPLY<CMesh>::Save( S2,S2NewName.c_str(),true,p);
    }

   return 0;
}

// -----------------------------------------------------------------------------------------------
