#include "trimesh_heatmethod.h"

int main(int argc, char const *argv[])
{
    if (argc != 2) {
        printf("Usage: %s <input_mesh>\n", argv[0]);
        exit(0);
    }
    // load mesh
    MyMesh m;
    if (vcg::tri::io::Importer<MyMesh>::Open(m, argv[1]) != 0)
    {
        printf("Error reading file  %s\n", argv[1]);
        exit(0);
    }
    // set initial conditions of the system
    Eigen::VectorXd initialConditions(m.VN());
    // srand(time(NULL));
    int random_source = rand() % m.VN();
    for (int i = 0; i < m.VN(); ++i){
        if (i == random_source)
            initialConditions(i) = 1;
        else
            initialConditions(i) = 0;
    }
    // compute geodesic
    Eigen::VectorXd distance = computeHeatMethodGeodesic(m, initialConditions);
    // save geodesic and print to stdout
    for (int i = 0; i < m.VN(); ++i){
        m.vert[i].Q() = distance(i);
    }
    printVectorXd(distance);
    // set source vertex for reference (will be removed later)
    m.vert[random_source].C() = vcg::Color4b::Red;
    // save mesh
    // NOTE: storing quality as double will make the PLY file unreadable
    vcg::tri::io::ExporterPLY<MyMesh>::Save(m, "distance_mesh.ply", vcg::tri::io::Mask::IOM_VERTCOLOR | vcg::tri::io::Mask::IOM_VERTQUALITY); 
    return 0;
}
