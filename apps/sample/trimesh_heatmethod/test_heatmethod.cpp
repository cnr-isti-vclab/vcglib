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
    
    std::cout << "Initial conditions..." << std::endl;
    Eigen::VectorXd initialConditions(m.VN());
    int random_source = rand() % m.VN();
    for (int i = 0; i < m.VN(); ++i){
        if (i == random_source)
            initialConditions(i) = 1;
        else
            initialConditions(i) = 0;
    }
    printVectorXd(initialConditions);
    Eigen::VectorXd distance = computeHeatMethodGeodesicVerbose(m, initialConditions);

    for (int i = 0; i < m.VN(); ++i){
        m.vert[i].Q() = distance(i);
    }
    m.vert[random_source].C() = vcg::Color4b::Red;

    // NOTE: storing quality as double will make the PLY file unreadable
    vcg::tri::io::ExporterPLY<MyMesh>::Save(m, "distance_mesh.ply", vcg::tri::io::Mask::IOM_VERTQUALITY); // vcg::tri::io::Mask::IOM_VERTCOLOR | vcg::tri::io::Mask::IOM_VERTQUALITY
    return 0;
}
