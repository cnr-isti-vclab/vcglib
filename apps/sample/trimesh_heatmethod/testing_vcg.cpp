#include "trimesh_heatmethod.h"

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include <iostream>
#include <vector>
#include <cmath>

int main(int argc, char const *argv[])
{
    // VCG
    // load mesh
    MyMesh m;
    if (vcg::tri::io::ImporterOFF<MyMesh>::Open(m, "icosahedron.off") != 0)
    {
        printf("Error reading file  %s\n", argv[1]);
        exit(0);
    }
    vcg::tri::UpdateTopology<MyMesh>::VertexFace(m);

    // NOTE: storing quality as double will make the PLY file unreadable
    vcg::tri::io::ExporterPLY<MyMesh>::Save(m, "distance_mesh.ply"); // vcg::tri::io::Mask::IOM_VERTCOLOR | vcg::tri::io::Mask::IOM_VERTQUALITY
    return 0;
}
