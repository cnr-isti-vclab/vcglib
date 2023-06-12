#include "trimesh_heatmethod.h"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <iostream>
#include <vector>


void print_matrix(Eigen::MatrixX3d a){
    for (int i=0; i<a.rows(); ++i){
        for (int j=0; j<a.cols(); ++j){
            std::cout << a(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

int main(int argc, char const *argv[])
{
    // EIGEN
    Eigen::Vector3d v0(1, 3, 4);
    Eigen::Vector3d v1(2, 1, 2);

    // dot product
    std::cout << v0.dot(v1) << std::endl;
    // cross product
    std::cout << v0.cross(v1) << std::endl;
    // norm of cross product
    std::cout << v0.cross(v1).norm() << std::endl;
    // cotan
    std::cout << v0.dot(v1) / (v0.cross(v1)).norm() << std::endl;

    Eigen::MatrixX3d heatGradientField(10, 3);
    for (int i=0; i<10 ; ++i){
        heatGradientField.row(i) = Eigen::Vector3d(i, 2*i, 3*i);
    }
    
    print_matrix(heatGradientField);
    
    return 0;
}
