#ifndef NXS_BUILD_H
#define NXS_BUILD_H

#include "vert_remap.h"
#include "crude.h"
#include "nexus.h"

#include <vector>

namespace nxs {

struct RemapLink {
  unsigned int rel_vert;
  unsigned int patch;
  unsigned int abs_vert;
};

void RemapVertices(Crude &crude,
		   VertRemap &vert_remap,
		   VFile<unsigned int> &face_remap,	 
		   std::vector<unsigned int> &patch_verts);

void NexusAllocate(Crude &crude,
		   nxs::Nexus &nexus,
		   VFile<unsigned int> &face_remap,
		   std::vector<unsigned int> &patch_faces,
		   std::vector<unsigned int> &patch_verts);

void NexusFill(Crude &crude,
	       Nexus &nexus,
	       VertRemap &vert_remap,
	       VFile<RemapLink> &border_remap);

void NexusFixBorder(Nexus &nexus, 
		    VFile<RemapLink> &border_remap);

}

#endif
