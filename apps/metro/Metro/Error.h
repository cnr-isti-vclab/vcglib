#ifndef __METRO__VCG__HAUS
#define __METRO__VCG__HAUS

#include "sampling.h"

template <class MESH_TYPE>
double Distance(MESH_TYPE & S1,MESH_TYPE & S2){
	  int   flags = 0x18;
		int n_samples_target = 1.5 * __max(S1.fn,S2.fn);

		double dist12_max,dist21_max;
		// compute face information
		tri::UpdateEdges::Set(S1);
		tri::UpdateEdges::Set(S2);

		// set bounding boxes for S1 and S2
		tri::UpdateBounding::Box(S1);
		tri::UpdateBounding::Box(S2);

    // set Bounding Box.
		Box3<typename MESH_TYPE::ScalarType>    bbox, tmp_bbox_M1=S1.bbox, tmp_bbox_M2=S2.bbox;
    bbox.Add(S1.bbox);
    bbox.Add(S2.bbox);
		bbox.InflateFix(0.1);
		S1.bbox = bbox;
		S2.bbox = bbox;	
		
		Sampling<MESH_TYPE> ForwardSampling(S1,S2);
		Sampling<MESH_TYPE> BackwardSampling(S2,S1);

		ForwardSampling.SetFlags(flags);
		BackwardSampling.SetFlags(flags);

		ForwardSampling.SetSamplesTarget(n_samples_target);
		BackwardSampling.SetSamplesTarget(n_samples_target);

		ForwardSampling.Hausdorff();
    dist12_max  = ForwardSampling.GetDistMax();

		BackwardSampling.Hausdorff();
    dist21_max  = BackwardSampling.GetDistMax();

		return math::Max(dist12_max,dist12_max);

	}

#endif


