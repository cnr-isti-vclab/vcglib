/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2023                                           \/)\/    *
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
#ifndef __VCG_GEODESIC_HEAT
#define __VCG_GEODESIC_HEAT
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/stat.h>

#include <vector>
#include <memory>

namespace vcg{
namespace tri{

template <class MeshType>
class GeodesicHeat{
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::VertexIterator VertexIterator;
    typedef typename MeshType::VertexPointer VertexPointer;
    typedef typename MeshType::FacePointer FacePointer;
    typedef typename MeshType::FaceIterator FaceIterator;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::ScalarType ScalarType;

public:
    /** @brief Fills the matrix with dual cell area of each vertex.
     * Additionally saves the area of faces in the face quality field
     *
     * @param mesh the mesh
     * @param mass a sparse matrix to store values into
    */
    static void buildMassMatrix(MeshType &mesh, Eigen::SparseMatrix<double> &mass){
        mass.resize(mesh.VN(), mesh.VN());
        mass.reserve(Eigen::VectorXi::Constant(mesh.VN(),1));
        vcg::tri::UpdateQuality<MeshType>::FaceArea(mesh);

        // compute area of the dual cell for each vertex
        for (int i = 0; i < mesh.VN(); ++i){
            VertexType *vp = &mesh.vert[i];

            std::vector<FacePointer> faces;
            std::vector<int> indices;
            vcg::face::VFStarVF<FaceType>(vp, faces, indices);

            double area = 0;
            for (int j = 0; j < faces.size(); ++j)
            {
                area += faces[j]->Q();
            }
            area /= 3;
            mass.insert(i, i) = area;
        }

        mass.makeCompressed();
    }

    /** @brief Fills the matrix with cotan weights. Only lower triangular values are saved.
     *
     * @param mesh the mesh
     * @param mass a sparse matrix to store values into
    */
    static void buildCotanLowerTriMatrix(MeshType &mesh, Eigen::SparseMatrix<double> &cotanOperator){
        cotanOperator.resize(mesh.VN(), mesh.VN());

        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        tripletList.reserve(3*mesh.VN() + 3*mesh.FN());

        // compute cotan weights
        for (int i = 0; i < mesh.FN(); ++i){
            FacePointer fi = &mesh.face[i];

            VertexPointer v0 = fi->V(0);
            VertexPointer v1 = fi->V(1);
            VertexPointer v2 = fi->V(2);

            vcg::Point3f p0 = v0->P();
            vcg::Point3f p1 = v1->P();
            vcg::Point3f p2 = v2->P();

            Eigen::Vector3d e0 = toEigen(p2 - p1);
            Eigen::Vector3d e1 = toEigen(p0 - p2);
            Eigen::Vector3d e2 = toEigen(p1 - p0);

            // first edge is inverted to get correct orientation
            double alpha0 = cotan(-e1, e2) / 2;
            double alpha1 = cotan(-e2, e0) / 2;
            double alpha2 = cotan(-e0, e1) / 2;

            int i0 = vcg::tri::Index(mesh,v0);
            int i1 = vcg::tri::Index(mesh,v1);
            int i2 = vcg::tri::Index(mesh,v2);

            // save only lower triangular part
            if (i0 > i1)
                tripletList.push_back(T(i0,i1,alpha2));
            else
                tripletList.push_back(T(i1,i0,alpha2));
            if (i0 > i2)
                tripletList.push_back(T(i0,i2,alpha1));
            else
                tripletList.push_back(T(i2,i0,alpha1));
            if (i1 > i2)
                tripletList.push_back(T(i1,i2,alpha0));
            else
                tripletList.push_back(T(i2,i1,alpha0));

            tripletList.push_back(T(i0,i0,-(alpha1 + alpha2)));
            tripletList.push_back(T(i1,i1,-(alpha0 + alpha2)));
            tripletList.push_back(T(i2,i2,-(alpha0 + alpha1)));
        }

        cotanOperator.setFromTriplets(tripletList.begin(), tripletList.end());
        cotanOperator.makeCompressed();
    }

    /** @brief given a mesh returns the average weight length
     * @param mesh the mesh
     * @return double the average edge length
    */
    static double computeAverageEdgeLength(MeshType &mesh){
        return vcg::tri::Stat<MeshType>::ComputeFaceEdgeLengthAverage(mesh);
    }

    /** @brief Computes the gradient of a (vertex sampled) scalar field with respect to faces
     * @param mesh the mesh
     * @param scalarField a scalar field of size (VN)
     * @return Eigen::MatrixX3d a vector field of with size (FN, 3)
     *
     * Requires both face normals and quality. The latter containing face areas.
    */
    static Eigen::MatrixX3d computeFaceGradient(MeshType &mesh, const Eigen::VectorXd &scalarField){
        Eigen::MatrixX3d faceGradientField(mesh.FN(), 3);
        // compute gradient of the scalar function at each face
        for (int i = 0; i < mesh.FN(); ++i){
            FacePointer fp = &mesh.face[i];

            VertexPointer v0 = fp->V(0);
            VertexPointer v1 = fp->V(1);
            VertexPointer v2 = fp->V(2);

            vcg::Point3f p0 = v0->P();
            vcg::Point3f p1 = v1->P();
            vcg::Point3f p2 = v2->P();

            int i0 = vcg::tri::Index(mesh,v0);
            int i1 = vcg::tri::Index(mesh,v1);
            int i2 = vcg::tri::Index(mesh,v2);

            // normal unit vector
            Eigen::Vector3d n = toEigen(fp->N());
            n /= n.norm();
            // face area
            double faceArea = fp->Q();
            // edge unit vectors (counter-clockwise)
            Eigen::Vector3d e0 = toEigen(p2 - p1);
            Eigen::Vector3d e1 = toEigen(p0 - p2);
            Eigen::Vector3d e2 = toEigen(p1 - p0);
            // gradient unit vectors
            Eigen::Vector3d g0 = n.cross(e0); //v0 grad
            Eigen::Vector3d g1 = n.cross(e1); //v1 grad
            Eigen::Vector3d g2 = n.cross(e2); //v2 grad

            // add vertex gradient contributions
            Eigen::Vector3d tri_grad = (g0 * scalarField(i0) + g1 * scalarField(i1) + g2 * scalarField(i2)) / (2 * faceArea);

            faceGradientField.row(i) = tri_grad;
        }
        return faceGradientField;
    }

    /** @brief Computes the per-row normalization of a given vector field
     * @param field a vector field
     * @return Eigen::MatrixX3d the per-row normalized vector field
    */
    static Eigen::MatrixX3d normalizeVectorField(const Eigen::MatrixX3d &field){
        Eigen::MatrixX3d normalizedField(field.rows(), 3);
        // normalize vector field at each vertex
        for (int i = 0; i < field.rows(); ++i){
            Eigen::Vector3d v = field.row(i);
            normalizedField.row(i) = v / v.norm();
        }
        return normalizedField;
    }

    /** @brief Computes the per-vertex divergence of a (normalized) per-face gradient (vector field)
     * @param mesh the mesh
     * @param field a vector field of size (FN, 3)
     * @return Eigen::VectorXd the per-vertex divergence
    */
    static Eigen::VectorXd computeVertexDivergence(MeshType &mesh, const Eigen::MatrixX3d &field){
        Eigen::VectorXd divergence(mesh.VN());
        divergence.setZero();

        // compute divergence of vector field at each vertex
        for (int i = 0; i < mesh.VN(); ++i){
            VertexPointer vp = &mesh.vert[i];

            std::vector<FacePointer> faces;
            std::vector<int> indices;
            vcg::face::VFStarVF<FaceType>(vp, faces, indices);
            for (int j = 0; j < faces.size(); ++j)
            {
                FacePointer fp = faces[j];
                int index = indices[j];

                vcg::Point3f p0 = fp->V(0)->P();
                vcg::Point3f p1 = fp->V(1)->P();
                vcg::Point3f p2 = fp->V(2)->P();

                // edge vectors
                Eigen::Vector3d el, er, eo; //left, right, opposite to vp
                if (index == 0){
                    el = toEigen(p2 - p0); //-e1
                    er = toEigen(p1 - p0); //e2
                    eo = toEigen(p2 - p1); //e0
                } else if (index == 1){
                    el = toEigen(p0 - p1); //-e2
                    er = toEigen(p2 - p1); //e0
                    eo = toEigen(p0 - p2); //e1
                } else if (index == 2){
                    el = toEigen(p1 - p2); //-e0
                    er = toEigen(p0 - p2); //e1
                    eo = toEigen(p1 - p0); //e2
                }
                // compute left and right cotangents
                double cotl = cotan(-el, -eo);
                double cotr = cotan(-er, eo);
                // add divergence contribution of given face
                Eigen::Vector3d x = field.row(vcg::tri::Index(mesh,fp));
                divergence(i) += (cotl * er.dot(x) + cotr * el.dot(x)) / 2;
            }
        }
        return divergence;
    }

    /**
     * @brief Computes an approximated geodesic distance using the heat method.
     *
     * @param mesh the mesh
     * @param sources a vector of source points
     * @param m a parameter regulating the timestep taken in the backward Euler method (default = 1)
     * @return true in case the computation was successful
     *
     * Approximated geodesic distance is stored in the vector quality field.
     * If the underlying factorization fails or the linear system cannot be solved false is returned.
     *
     * THIS METHOD IS NOT APPROPRIATE FOR MULTIPLE CALLS ON DIFFERENT SOURCE POINTS.
     * If multiple calls are required its best to first build the factorizations with
     * GeodesicHeat::BuildCache(mesh, m) then call ComputeFromCache(mesh, source, cache)
     * this will avoid recomputing the underlying factorizations for each set of source points.
     *
     */
    static bool Compute(MeshType &mesh, const std::vector<VertexPointer> sources, float m = 1){
        vcg::tri::RequireVFAdjacency<MeshType>(mesh);
        vcg::tri::RequireFFAdjacency<MeshType>(mesh);
        vcg::tri::RequirePerFaceNormal<MeshType>(mesh);
        vcg::tri::RequirePerFaceQuality<MeshType>(mesh);
        vcg::tri::RequirePerVertexQuality<MeshType>(mesh);

        // initialization of required components
        vcg::tri::Allocator<MeshType>::CompactEveryVector(mesh);
        vcg::tri::UpdateTopology<MeshType>::VertexFace(mesh);
        vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(mesh);

        if(sources.empty())	return false;

        Eigen::VectorXd sourcePoints(mesh.VN());
        sourcePoints.setZero();
        for(VertexPointer vp : sources){
            sourcePoints(vcg::tri::Index(mesh, vp)) = 1;
        }

        Eigen::SparseMatrix<double> massMatrix;
        Eigen::SparseMatrix<double> cotanMatrix;
        double averageEdgeLength;

        buildMassMatrix(mesh, massMatrix);
        buildCotanLowerTriMatrix(mesh, cotanMatrix);
        averageEdgeLength = computeAverageEdgeLength(mesh);

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver1;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver2;

        // core of the heat method
        double timestep = m * averageEdgeLength * averageEdgeLength;
        solver1.compute(massMatrix - timestep * cotanMatrix);
        if (solver1.info() != Eigen::Success) return false;
        solver2.compute(cotanMatrix);
        if (solver2.info() != Eigen::Success) return false;

        Eigen::VectorXd heatflow = solver1.solve(sourcePoints); // (VN)
        if (solver1.info() != Eigen::Success) return false;
        Eigen::MatrixX3d heatGradient = computeFaceGradient(mesh, heatflow); // (FN, 3)
        Eigen::MatrixX3d unitVectorField = normalizeVectorField(-heatGradient); // (FN, 3)
        Eigen::VectorXd divergence = computeVertexDivergence(mesh, unitVectorField); // (VN)
        Eigen::VectorXd geodesicDistance = solver2.solve(divergence); // (VN)
        if (solver2.info() != Eigen::Success) return false;

        // shift to impose dist(source) = 0
        geodesicDistance.array() -= geodesicDistance.minCoeff();

        // save geodesic in vertex quaity field
        for(int i=0; i < mesh.VN(); i++){
            mesh.vert[i].Q() = geodesicDistance(i);
        }
        return true;
    }

    // used to cache factorizations
    typedef typename std::pair<
        std::shared_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>,
        std::shared_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>
    > GeodesicHeatCache;

    /**
     * @brief Precomputes matrix factorizations required by the heat method
     *
     * @param mesh the mesh
     * @param m a parameter regulating the timestep taken in the backward Euler method (default = 1)
     * @return GeodesicHeatCache object containing the matrix factorizations
     *
     * This method returns the cache to use in ComputeFromCache.
     * Note that when the mesh changes this cache should be rebuilt.
     *
     * If the factorization fails no error is thrown (errors can be caught during ComputeFromCache).
     */
    static GeodesicHeatCache BuildCache(MeshType &mesh, float m = 1){
        vcg::tri::RequireVFAdjacency<MeshType>(mesh);
        vcg::tri::RequireFFAdjacency<MeshType>(mesh);
        vcg::tri::RequirePerFaceNormal<MeshType>(mesh);
        vcg::tri::RequirePerFaceQuality<MeshType>(mesh);
        vcg::tri::RequirePerVertexQuality<MeshType>(mesh);

        // initialization of required components
        vcg::tri::Allocator<MeshType>::CompactEveryVector(mesh);
        vcg::tri::UpdateTopology<MeshType>::VertexFace(mesh);
        vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(mesh);

        GeodesicHeatCache cache = std::make_pair(
            std::make_shared<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>(),
            std::make_shared<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>()
        );

        Eigen::SparseMatrix<double> massMatrix;
        Eigen::SparseMatrix<double> cotanMatrix;
        double averageEdgeLength;

        buildMassMatrix(mesh, massMatrix);
        buildCotanLowerTriMatrix(mesh, cotanMatrix);
        averageEdgeLength = computeAverageEdgeLength(mesh);

        // compute factorizations
        double timestep = m * averageEdgeLength * averageEdgeLength;
        std::get<0>(cache)->compute(massMatrix - timestep * cotanMatrix);
        std::get<1>(cache)->compute(cotanMatrix);

        return cache;
    }

    /**
     * @brief Computes an approximated geodesic distance using the heat method.
     *
     * @param mesh the mesh
     * @param sources a vector of source points
     * @param cache an (up to date) GeodesicHeatCache of the mesh
     * @return true if computation was successful
     *
     * We assume that the face quality field has not changed after BuildCache was called.
     *
     * Approximated geodesic distance is stored in the vector quality field.
     */
    static bool ComputeFromCache(MeshType &mesh, const std::vector<VertexPointer> sources, GeodesicHeatCache &cache){
        if(sources.empty())	return false;

        Eigen::VectorXd sourcePoints(mesh.VN());
        sourcePoints.setZero();
        for(VertexPointer vp : sources){
            sourcePoints(vcg::tri::Index(mesh, vp)) = 1;
        }

        Eigen::VectorXd heatflow = std::get<0>(cache)->solve(sourcePoints); // (VN)
        if (std::get<0>(cache)->info() != Eigen::Success) return false;
        Eigen::MatrixX3d heatGradient = computeFaceGradient(mesh, heatflow); // (FN, 3)
        Eigen::MatrixX3d unitVectorField = normalizeVectorField(-heatGradient); // (FN, 3)
        Eigen::VectorXd divergence = computeVertexDivergence(mesh, unitVectorField); // (VN)
        Eigen::VectorXd geodesicDistance = std::get<1>(cache)->solve(divergence); // (VN)
        if (std::get<1>(cache)->info() != Eigen::Success) return false;

        // shift to impose dist(source) = 0
        geodesicDistance.array() -= geodesicDistance.minCoeff();

        // save geodesic in vertex quality field
        for(int i=0; i < mesh.VN(); i++){
            mesh.vert[i].Q() = geodesicDistance(i);
        }
        return true;
    }

    private:
    static inline Eigen::Vector3d toEigen(const vcg::Point3f& p){
        return Eigen::Vector3d(p.X(), p.Y(), p.Z());
    }
    static inline double cotan(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1){
        // cos(theta) / sin(theta)
        return v0.dot(v1) / v0.cross(v1).norm();
    }
};

}
}

#endif
