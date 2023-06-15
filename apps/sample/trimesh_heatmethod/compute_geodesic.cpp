#include "trimesh_heatmethod.h"

inline void printUsage(char const * prog){
    printf("Usage: %s [--save-csv] <input_mesh>\n", prog);
    printf("Options:\n");
    printf("\t--save-csv\tSaves intermediate steps as csv files. Files are saved in ./output_csv.\n");
}

inline std::tuple<char const*, bool> parseArgs(int argc, char const *argv[]){
    if (argc < 2 || argc > 3) {
        printUsage(argv[0]);
        exit(0);
    }
    if (argc == 2){
        return std::make_tuple(argv[1], false);
    }
    else if (argc == 3){
        if (std::strcmp(argv[1], "--save-csv")){
            printUsage(argv[0]);
            std::cerr << "Invalid argument: " << argv[1] << std::endl;
            exit(0);
        }
        return std::make_tuple(argv[2], true);
    }
    else{
        printUsage(argv[0]);
        exit(0);
    }
}


int main(int argc, char const *argv[])
{
    std::tuple<char const*, bool> options = parseArgs(argc, argv);
    char const * path; bool save_csv;
    std::tie (path, save_csv) = options;

    // load mesh
    MyMesh m;
    if (vcg::tri::io::Importer<MyMesh>::Open(m, path) != 0)
    {
        printf("Error reading file  %s\n", path);
        exit(0);
    }

    // set initial conditions of the system
    Eigen::VectorXd initialConditions(m.VN());
    int random_source = rand() % m.VN();
    for (int i = 0; i < m.VN(); ++i){
        if (i == random_source)
            initialConditions(i) = 1;
        else
            initialConditions(i) = 0;
    }
    std::cout << "Source vertex: " << random_source << std::endl;
    std::cout << toEigen(m.vert[random_source].P()) << std::endl;

    // compute geodesic
    SaveMeshMask saveMask = SaveMeshMask::NONE;
    if (save_csv)
        saveMask = SaveMeshMask::CSV_ALL;
    Eigen::VectorXd distance = computeHeatMethodGeodesic(m, initialConditions, 1, saveMask);

    // save geodesic and print to stdout
    for (int i = 0; i < m.VN(); ++i){
        m.vert[i].Q() = distance(i);
    }
    m.vert[random_source].C() = vcg::Color4b::Red;

    // save mesh (NOTE: storing quality as double will make the PLY file unreadable)
    vcg::tri::io::ExporterPLY<MyMesh>::Save(m, "distance_mesh.ply", vcg::tri::io::Mask::IOM_VERTCOLOR | vcg::tri::io::Mask::IOM_VERTQUALITY); 
    return 0;
}
