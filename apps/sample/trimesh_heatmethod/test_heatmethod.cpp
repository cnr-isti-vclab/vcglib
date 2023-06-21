#include "trimesh_heatmethod.h"

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>


inline void printSparseMatrix(Eigen::SparseMatrix<double> &mat){
    for (int k=0; k < mat.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
            std::cout << "(" << it.row() << "," << it.col() << ") = " << it.value() << std::endl;
        }
    }
}

inline void printVectorXd(Eigen::VectorXd &vec){
    for (int i = 0; i < vec.size(); ++i){
        std::cout << i << " : " << vec[i] << std::endl;
    }
}

inline void printVectorX3d(Eigen::MatrixX3d &vec){
    for (int i = 0; i < vec.rows(); ++i){
        std::cout << vec.row(i) << std::endl;
    }
}

inline void printCompareSparseMatrices(Eigen::SparseMatrix<double> &mat1, Eigen::SparseMatrix<double> &mat2){
    for (int k=0; k < mat1.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat1,k); it; ++it){
            std::cout << "(" << it.row() << "," << it.col() << ") = " << it.value() << " | " << mat2.coeff(it.row(), it.col()) << std::endl;
        }
    }
}


inline Eigen::VectorXd computeHeatMethodGeodesicVerbose(CMeshO &mesh, const Eigen::VectorXd &init_cond, double m = 1){
    vcg::tri::UpdateTopology<CMeshO>::VertexFace(mesh);
    vcg::tri::UpdateTopology<CMeshO>::FaceFace(mesh);
    vcg::tri::UpdateNormal<CMeshO>::PerFaceNormalized(mesh);

    std::cout << "Computing Mass..." << std::endl;
    Eigen::SparseMatrix<double> mass(mesh.VN(), mesh.VN());
    buildMassMatrix(mesh, mass);
    // printSparseMatrix(mass);

    std::cout << "Computing Cotan..." << std::endl;
    Eigen::SparseMatrix<double> cotanOperator(mesh.VN(), mesh.VN());
    buildCotanLowerTriMatrix(mesh, cotanOperator);
    // printSparseMatrix(cotanOperator);

    std::cout << "Computing Edge Length..." << std::endl;
    double avg_edge_len = computeAverageEdgeLength(mesh);
    std::cout << "Average Edge: " << avg_edge_len << std::endl;
    double timestep = m * avg_edge_len * avg_edge_len;
    std::cout << "Timestep: " << timestep << std::endl;
    Eigen::SparseMatrix<double> system1(mesh.VN(), mesh.VN());
    system1 = mass - timestep * cotanOperator;
    // printSparseMatrix(system1);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    std::cout << "Cholesky Factorization 1..." << std::endl;
    solver.compute(system1);
    if(solver.info() != Eigen::Success) {
        std::cout << "ERROR: Cholesky Factorization 1 failed" << std::endl; 
    }
    Eigen::VectorXd heatflow = solver.solve(init_cond);
    if(solver.info() != Eigen::Success) {
        std::cout << "ERROR: Solving System 1 failed" << std::endl; 
    }
    // printVectorXd(heatflow);

    std::cout << "Computing Gradient..." << std::endl;
    Eigen::MatrixX3d heatGradient = computeVertexGradient(mesh, heatflow);
    // printVectorX3d(heatGradient);

    std::cout << "Normalizing Gradient..." << std::endl;
    Eigen::MatrixX3d normalizedVectorField = normalizeVectorField(-heatGradient);
    // printVectorX3d(normalizedVectorField);
    
    std::cout << "Computing Divergence..." << std::endl;
    Eigen::VectorXd divergence = computeVertexDivergence(mesh, normalizedVectorField);
    // printVectorXd(divergence);

    Eigen::SparseMatrix<double> system2(mesh.VN(), mesh.VN());
    system2 = cotanOperator; //+ 1e-6 * Eigen::Matrix<double,-1,-1>::Identity(mesh.VN(), mesh.VN());
    // printSparseMatrix(system2);

    std::cout << "Cholesky Factorization 2..." << std::endl;
    solver.compute(system2);
    if(solver.info() != Eigen::Success) {
        std::cout << "ERROR: Cholesky Factorization 2 failed" << std::endl; 
    }
    Eigen::VectorXd geodesicDistance = solver.solve(divergence);
    if(solver.info() != Eigen::Success) {
        std::cout << "ERROR: Solving System 2 failed" << std::endl; 
    }
    // shift to impose boundary conditions (dist(d) = 0 \forall d \in init_cond)
    geodesicDistance.array() -= geodesicDistance.minCoeff();
    printVectorXd(geodesicDistance);

    return geodesicDistance;
}


int main(int argc, char const *argv[])
{
    // load mesh
    CMeshO m;
    if (vcg::tri::io::ImporterOFF<CMeshO>::Open(m, "models/icosahedron.off") != 0)
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
    std::cout << "Source point id: " << random_source << std::endl;
    std::cout << toEigen(m.vert[random_source].P()) << std::endl;

    Eigen::VectorXd distance = computeHeatMethodGeodesicVerbose(m, initialConditions);
    computeHeatMethodGeodesic(m, initialConditions);

    for (int i = 0; i < m.VN(); ++i){
        m.vert[i].Q() = distance(i);
    }

    // NOTE: storing quality as double will make the PLY file unreadable
    vcg::tri::io::ExporterPLY<CMeshO>::Save(m, "distance_mesh.ply", vcg::tri::io::Mask::IOM_VERTQUALITY); // vcg::tri::io::Mask::IOM_VERTCOLOR | vcg::tri::io::Mask::IOM_VERTQUALITY
    return 0;
}
