// -----------------------------------------------------------------------------------------------

// standard libraries
#include <time.h>

// project definitions.
#include "defs.h"
#include "sampling.h"
#include "mesh_type.h"
#include <vcg/complex/trimesh/update/edges.h>
#include <vcg/complex/trimesh/update/bounding.h>
//#include <wrap/io_trimesh/import_smf.h>
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>



// -----------------------------------------------------------------------------------------------


////////////////// Command line Flags and parameters 
bool ComputeHistFlag                = false;
bool VertexSampleFlag               = true;
bool EdgeSampleFlag                 = true;
bool FaceSampleFlag                 = true;
bool MontecarloSamplingFlag         = false;
bool SubdivisionSamplingFlag        = false;
bool SimilarTrianglesSamplingFlag   = false;
bool NumberOfSamples                = false;
bool SamplesPerAreaUnit             = false;
bool SaveErrorDisplacement          = false;
bool SaveErrorAsColour              = false;
bool IgnoreUnreferred							= true;
// -----------------------------------------------------------------------------------------------


inline char* GetExtension(char* filename)
{
		size_t i;
    for( i=strlen(filename)-1; i >= 0; i--)
        if(filename[i] == '.')
            break;
    if(i > 0)
        return &(filename[i+1]);
    else
        return NULL;
}


int main(int argc, char**argv)
{
    CMesh                 S1, S2;
    double                dist1_max, dist1_mean, dist1_RMS, volume_1;
    double                dist2_max, dist2_mean, dist2_RMS, volume_2;
    double                mesh_dist_max;
    unsigned long         n_samples_target, n_samples_output, elapsed_time;
    double								n_samples_per_area_unit;
    int                   flags, flags_fwd, flags_back, n_samples_area, n_samples_edge, n_samples_vertex, err;
    char                 *fmt, *hist_filename, *new_mesh_filename, *new_mesh_filename_2;
    char                  fname_1[] = STR_NEW_MESH_FILENAME_DEFAULT, fname_2[] = STR_NEW_MESH_FILENAME_DEFAULT_2;
    //FILE                 *fd;

    // print program info
    printf("-------------------------------\n"
           "             Metro\n"
           "   release date: "__DATE__"\n"
           "-------------------------------\n\n");

    // load input meshes.
    if(argc <= 2)
    {
        printf(MSG_ERR_N_ARGS);
        exit(-1);
    }

    // load mesh M1.
    if(!(fmt = GetExtension(argv[1])))
    {
        printf(MSG_ERR_UNKNOWN_FORMAT, fmt);
        exit(-1);
    }
    if(!strcmp("ply", fmt))
			{
		printf("reading the mesh `%s'...", argv[1]);		  
		err = tri::io::ImporterPLY<CMesh>::Open(S1,argv[1]);
	}
 /*   else
        if(!strcmp(file_ext_smf, fmt))
		{
			printf("reading the mesh `%s'...", argv[1]);
			err = tri::io::ImporterSMF::Open(S1,argv[1]);
		}*/
        else
        {
            printf(MSG_ERR_UNKNOWN_FORMAT, fmt);
            exit(-1);
        }
    if(err < 0)
    {
		printf("\n");
        printf(MSG_ERR_MESH_LOAD);
        printf(" error number %d",err);
        exit(-1);
    }
	else
		printf("done\n");

    // load mesh M2.
    if(!(fmt = GetExtension(argv[2])))
    {
        printf(MSG_ERR_UNKNOWN_FORMAT, fmt);
        exit(-1);
    }
    if(!strcmp("ply", fmt))
	{
		printf("reading the mesh `%s'...", argv[2]);
       err = tri::io::ImporterPLY<CMesh>::Open(S2,argv[2]);
	}
    /*else
        if(!strcmp(file_ext_smf, fmt))
     	{
			printf("reading the mesh `%s'...", argv[2]);
		    err = S2.Load_Smf(argv[2]);
		}*/
        else
        {
            printf(MSG_ERR_UNKNOWN_FORMAT, fmt);
            exit(-1);
        }
    if(err < 0)
    {
        printf(MSG_ERR_MESH_LOAD);
        exit(-1);
    }
	else
		printf("done\n");

    // parse command line.
	for(int i=3; i < argc;)
	{
		if(argv[i][0]=='-')
			switch(argv[i][1])
			{ 
				case CMD_LINE_ARG_HIST          :	ComputeHistFlag     = true;     hist_filename = &(argv[i][2]);
                                                    if(hist_filename[0] == '\0')
                                                        strcpy(hist_filename, STR_HIST_FILENAME_DEFAULT);
                                                    break;
				case CMD_LINE_ARG_VERTEX_SAMPLE :	VertexSampleFlag    = false;    break;
				case CMR_LINE_ARG_IGNORE_UNREF  : IgnoreUnreferred    = true;     break;
				case CMD_LINE_ARG_EDGE_SAMPLE   :	EdgeSampleFlag      = false;    break;
				case CMD_LINE_ARG_FACE_SAMPLE   :	FaceSampleFlag      = false;    break;
                case CMD_LINE_ARG_SAMPLE_TYPE   :
                    switch(argv[i][2])
                    {
                        case CMD_LINE_ARG_MONTECARLO_SAMPLING        :  MontecarloSamplingFlag       = true;     break;
                        case CMD_LINE_ARG_SUBDIVISION_SAMPLING       :  SubdivisionSamplingFlag      = true;     break;
                        case CMD_LINE_ARG_SIMILAR_TRIANGLES_SAMPLING :  SimilarTrianglesSamplingFlag = true;     break;
                        default  :  printf(MSG_ERR_INVALID_OPTION, argv[i]);
                                    exit(0);
                            }
                            break;
                case CMD_LINE_ARG_N_SAMPLES             :  NumberOfSamples       = true;     n_samples_target        = (unsigned long) atoi(&(argv[i][2]));          break;
                case CMD_LINE_ARG_SAMPLES_PER_AREA_UNIT :  SamplesPerAreaUnit    = true;     n_samples_per_area_unit = (unsigned long) atoi(&(argv[i][2])); break;
                case CMD_LINE_ARG_SAVE_DISPLACEMENT     :  SaveErrorDisplacement = true;     new_mesh_filename = &(argv[i][2]);
                                                           if(new_mesh_filename[0] == '\0')
                                                               new_mesh_filename = fname_1;
                                                           break;
                case CMD_LINE_ARG_SAVE_ERROR_AS_COLOUR  :  SaveErrorAsColour     = true;     new_mesh_filename_2 = &(argv[i][2]);
                                                           if(new_mesh_filename_2[0] == '\0')
                                                               new_mesh_filename_2 = fname_2;
                                                           break;
				default  :  printf(MSG_ERR_INVALID_OPTION, argv[i]);
                            exit(0);
			}
	  i++;
	}

    // set sampling scheme
    int sampling_method = MontecarloSamplingFlag + SubdivisionSamplingFlag + SimilarTrianglesSamplingFlag;
    // defaults
    if(!sampling_method)
        SimilarTrianglesSamplingFlag = true;
    if(sampling_method > 1)
    {
        printf("Cannot choose more than one sampling method. Similar Triangles sampling assumed.\n");
        SimilarTrianglesSamplingFlag = true;
        MontecarloSamplingFlag       = false;
        SubdivisionSamplingFlag      = false;
    }
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
		bbox.InflateFix(0.02);
	S1.bbox = bbox;
	S2.bbox = bbox;

    // set flags.
	flags = 0;
   if(IgnoreUnreferred)
        flags |= SamplingFlags::INCLUDE_UNREFERENCED_VERTICES;
    if(ComputeHistFlag)
        flags |= SamplingFlags::HIST;
    if(VertexSampleFlag)
        flags |= SamplingFlags::VERTEX_SAMPLING;
    if(EdgeSampleFlag)
        flags |= SamplingFlags::EDGE_SAMPLING;
    if(FaceSampleFlag)
        flags |= SamplingFlags::FACE_SAMPLING;
    if(MontecarloSamplingFlag)
        flags |= SamplingFlags::MONTECARLO_SAMPLING;
    if(SubdivisionSamplingFlag)
        flags |= SamplingFlags::SUBDIVISION_SAMPLING;
    if(SimilarTrianglesSamplingFlag)
        flags |= SamplingFlags::SIMILAR_TRIANGLES_SAMPLING;
    flags_fwd  = flags;
    flags_back = flags;
    if(SaveErrorDisplacement)
    {
        if(S1.vn >= S2.vn)
            flags_fwd  |= SamplingFlags::SAVE_ERROR_DISPLACEMENT;
        else
            flags_back |= SamplingFlags::SAVE_ERROR_DISPLACEMENT;
    }

    if(SaveErrorAsColour)
    {
        if(S1.vn >= S2.vn)
            flags_fwd  |= SamplingFlags::SAVE_ERROR_AS_COLOUR;
        else
            flags_back |= SamplingFlags::SAVE_ERROR_AS_COLOUR;
    }
    
    // initialize time info.
    int t0=clock();

    // print mesh info.
    Sampling<CMesh> ForwardSampling(S1,S2);
    Sampling<CMesh> BackwardSampling(S2,S1);

    printf("Mesh info:\n");
    printf(" M1: '%s'\n\t%vertices  %7i\n\tfaces     %7i\n\tarea      %12.4f\n", argv[1], S1.vn, S1.fn, ForwardSampling.GetArea());
    printf("\tbbox (%7.4f %7.4f %7.4f)-(%7.4f %7.4f %7.4f)\n", tmp_bbox_M1.min[0], tmp_bbox_M1.min[1], tmp_bbox_M1.min[2], tmp_bbox_M1.max[0], tmp_bbox_M1.max[1], tmp_bbox_M1.max[2]);
    printf("\tbbox diagonal %f\n", (float)tmp_bbox_M1.Diag());
    printf(" M2: '%s'\n\t%vertices  %7i\n\tfaces     %7i\n\tarea      %12.4f\n", argv[2], S2.vn, S2.fn, BackwardSampling.GetArea());
    printf("\tbbox (%7.4f %7.4f %7.4f)-(%7.4f %7.4f %7.4f)\n", tmp_bbox_M2.min[0], tmp_bbox_M2.min[1], tmp_bbox_M2.min[2], tmp_bbox_M2.max[0], tmp_bbox_M2.max[1], tmp_bbox_M2.max[2]);
    printf("\tbbox diagonal %f\n", (float)tmp_bbox_M2.Diag());

    // Forward distance.
    printf("\nForward distance (M1 -> M2):\n");
    ForwardSampling.SetFlags(flags_fwd);
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
    dist1_mean = ForwardSampling.GetDistMean();
    dist1_RMS  = ForwardSampling.GetDistRMS();
    volume_1   = ForwardSampling.GetDistVolume();
    n_samples_output = ForwardSampling.GetNSamples();
    n_samples_area   = ForwardSampling.GetNAreaSamples();
    n_samples_edge   = ForwardSampling.GetNEdgeSamples();
    n_samples_vertex = ForwardSampling.GetNVertexSamples();
    printf("\ndistance:\n  max  : %f (%f  with respect to bounding box diagonal)\n  mean : %f\n  RMS  : %f\n", (float)dist1_max, (float)dist1_max/bbox.Diag(), (float)dist1_mean, (float)dist1_RMS);
    if(VertexSampleFlag)
        printf("# vertex samples %d\n", n_samples_vertex);
    if(EdgeSampleFlag)
        printf("# edge samples   %d\n", n_samples_edge);
    printf("# area samples   %d\n# total samples  %d\nsamples per area unit: %f\n\n", n_samples_area, n_samples_output, ForwardSampling.GetNSamplesPerAreaUnit());

    // Backward distance.
    printf("\nBackward distance (M2 -> M1):\n");
    BackwardSampling.SetFlags(flags_back);
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
    dist2_mean = BackwardSampling.GetDistMean();
    dist2_RMS  = BackwardSampling.GetDistRMS();
    volume_2   = BackwardSampling.GetDistVolume();
    n_samples_output = BackwardSampling.GetNSamples();
    n_samples_area   = BackwardSampling.GetNAreaSamples();
    n_samples_edge   = BackwardSampling.GetNEdgeSamples();
    n_samples_vertex = BackwardSampling.GetNVertexSamples();
    printf("\ndistance:\n  max  : %f (%f  with respect to bounding box diagonal)\n  mean : %f\n  RMS  : %f\n", (float)dist2_max, (float)dist2_max/bbox.Diag(), (float)dist2_mean, (float)dist2_RMS);
    if(VertexSampleFlag)
        printf("# vertex samples %d\n", n_samples_vertex);
    if(EdgeSampleFlag)
        printf("# edge samples   %d\n", n_samples_edge);
    printf("# area samples   %d\n# total samples  %d\nsamples per area unit: %f\n\n", n_samples_area, n_samples_output, BackwardSampling.GetNSamplesPerAreaUnit());

    // compute time info.
    elapsed_time = clock() - t0;

    // save error distribution histogram
    /*if(ComputeHistFlag)
    {
        const Hist  &hist1 = ForwardSampling.GetHist();
        const Hist  &hist2 = BackwardSampling.GetHist();
        if(!(fd = fopen(hist_filename, "w")))
        {
            printf(MSG_ERR_FILE_OPEN);
            exit(-1);
        }
        vector<int>::const_iterator    ii;
        vector<float>::const_iterator  fi;

        fprintf(fd, "error distribution histogram (forward distance)\n\n");
        for(ii=hist1.H.begin(), fi=hist1.R.begin(); ii != hist1.H.end(); ++fi,ii++)
            fprintf(fd, "%6.4f\t%d\n", *fi, *ii);
            
        fprintf(fd, "\n\nerror distribution histogram (backward distance)\n");
        for(ii=hist2.H.begin(), fi=hist2.R.begin(); ii != hist2.H.end(); ++fi,ii++)
            fprintf(fd, "%6.4f\t%d\n", *fi, *ii);
        
        fclose(fd);
    }*/

    // max distance.
    mesh_dist_max  = max(dist1_max , dist2_max);
    printf("\nHausdorff distance: %f (%f  with respect to bounding box diagonal)\nComputation time  : %d ms\n# samples/second  : %f\n\n", (float)mesh_dist_max, (float)mesh_dist_max/bbox.Diag(), (int)elapsed_time, (float)n_samples_output/(float)elapsed_time*2000.0F);

    // save error files.
    if((flags_fwd & SamplingFlags::SAVE_ERROR_DISPLACEMENT) && (flags_fwd & SamplingFlags::SAVE_ERROR_AS_COLOUR))
 //       if(strcmp(new_mesh_filename, new_mesh_filename_2))
        {
				vcg::tri::io::PlyInfo p; 
				p.mask|=vcg::tri::io::PLYMask::PM_VERTCOLOR|vcg::tri::io::PLYMask::PM_VERTQUALITY;

						tri::io::ExporterPLY<CMesh>::Save( S1,new_mesh_filename,true,p);
	          exit(0);
        }
    if((flags_back & SamplingFlags::SAVE_ERROR_DISPLACEMENT) && (flags_back & SamplingFlags::SAVE_ERROR_AS_COLOUR))
 //       if(strcmp(new_mesh_filename, new_mesh_filename_2))
        {
						vcg::tri::io::PlyInfo p;
					p.mask|=vcg::tri::io::PLYMask::PM_VERTCOLOR|vcg::tri::io::PLYMask::PM_VERTQUALITY;
						tri::io::ExporterPLY<CMesh>::Save( S2,new_mesh_filename,true,p);
						exit(0);
        }

    //if(flags_fwd & SamplingFlags::SAVE_ERROR_DISPLACEMENT)
    //    S1.SavePly(new_mesh_filename, CMesh::SM_ALL & (CMesh::SM_ALL ^ CMesh::SM_VERTCOLOR));
    //else
    //    if(flags_back & SamplingFlags::SAVE_ERROR_DISPLACEMENT)
    //        S2.SavePly(new_mesh_filename, CMesh::SM_ALL & (CMesh::SM_ALL ^ CMesh::SM_VERTCOLOR));
    //if(flags_fwd & SamplingFlags::SAVE_ERROR_AS_COLOUR)
    //    S1.SavePly(new_mesh_filename_2, CMesh::SM_ALL & (CMesh::SM_ALL ^ CMesh::SM_VERTQUALITY));
    //else
    //    if(flags_back & SamplingFlags::SAVE_ERROR_AS_COLOUR)
    //        S2.SavePly(new_mesh_filename_2, CMesh::SM_ALL & (CMesh::SM_ALL ^ CMesh::SM_VERTQUALITY));
return 0;
}

// -----------------------------------------------------------------------------------------------
