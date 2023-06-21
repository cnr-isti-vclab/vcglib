#ifndef TRIMESH_HEAT_METHOD
#define TRIMESH_HEAT_METHOD

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/quality.h>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>

namespace {
    class MyEdge;
    class MyFace;
    class MyVertex;

    struct MyUsedTypes : public vcg::UsedTypes<
        vcg::Use<MyVertex>::AsVertexType, 
        vcg::Use<MyEdge>::AsEdgeType, 
        vcg::Use<MyFace>::AsFaceType
    >{};

    class MyVertex : public vcg::Vertex<
        MyUsedTypes,
        vcg::vertex::Coord3f,
        vcg::vertex::VFAdj,
        vcg::vertex::Color4b, //
        vcg::vertex::Qualityf,
        vcg::vertex::BitFlags // needed for PLY export
    >{};
    class MyEdge : public vcg::Edge<
        MyUsedTypes
    >{};
    class MyFace : public vcg::Face<
        MyUsedTypes,
        vcg::face::VFAdj,
        vcg::face::FFAdj,
        vcg::face::VertexRef,
        vcg::face::Normal3f,
        vcg::face::Qualityf,
        vcg::face::Color4b //
    >{};

    inline Eigen::Vector3d toEigen(const vcg::Point3f& p)
    {
        return Eigen::Vector3d(p.X(), p.Y(), p.Z());
    };


    inline double cotan(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1)
    {
        // cos(theta) / sin(theta)
        return v0.dot(v1) / v0.cross(v1).norm();
    };
}

class CMeshO : public vcg::tri::TriMesh<
    std::vector<MyVertex>, 
    std::vector<MyFace>,
    std::vector<MyEdge>
>{};


// foward declarations
inline void buildMassMatrix(CMeshO &mesh, Eigen::SparseMatrix<double> &mass);
inline void buildCotanLowerTriMatrix(CMeshO &mesh, Eigen::SparseMatrix<double> &cotanOperator);
inline double computeAverageEdgeLength(CMeshO &mesh);
inline Eigen::MatrixX3d computeVertexGradient(CMeshO &mesh, const Eigen::VectorXd &heat);
inline Eigen::VectorXd computeVertexDivergence(CMeshO &mesh, const Eigen::MatrixX3d &field);
inline Eigen::MatrixX3d normalizeVectorField(const Eigen::MatrixX3d &field);
enum SaveMeshMask;
inline void saveVertexScalarFieldCSV(CMeshO &mesh, const Eigen::VectorXd scalarField, const char *fname);
inline void saveFaceVectorFieldCSV(CMeshO &mesh, const Eigen::MatrixX3d vectorField, const char *fname);


inline void buildMassMatrix(CMeshO &mesh, Eigen::SparseMatrix<double> &mass){
    // compute area of all faces
    for (CMeshO::FaceIterator fi = mesh.face.begin(); fi != mesh.face.end(); ++fi)
    {
        vcg::Point3f p0 = fi->V(0)->P();
        vcg::Point3f p1 = fi->V(1)->P();
        vcg::Point3f p2 = fi->V(2)->P();
        double e0 = toEigen(p1 - p0).norm();
        double e1 = toEigen(p2 - p0).norm();
        double e2 = toEigen(p2 - p1).norm();
        double s = (e0 + e1 + e2) / 2;
        double area = std::sqrt(s * (s - e0) * (s - e1) * (s - e2));
        // we store face areas in quality field to avoid a hash table
        // this will also be useful for the gradient computation
        fi->Q() = area;
    }
    // compute area of the dual cell for each vertex
    for (int i = 0; i < mesh.VN(); ++i){
        CMeshO::VertexType *vp = &mesh.vert[i];

        std::vector<CMeshO::FaceType*> faces;
        std::vector<int> indices;
        vcg::face::VFStarVF<CMeshO::FaceType>(vp, faces, indices);
        
        double area = 0;
        for (int j = 0; j < faces.size(); ++j)
        {
            area += faces[j]->Q();
        }
        area /= 3;
        mass.coeffRef(i, i) = area;
    }
}


inline void buildCotanLowerTriMatrix(CMeshO &mesh, Eigen::SparseMatrix<double> &cotanOperator){
    // initialize a hashtable from vertex pointers to ids
    std::unordered_map<CMeshO::VertexType*, int> vertex_ids;
    for (int i = 0; i < mesh.VN(); ++i){
        vertex_ids[&mesh.vert[i]] = i;
    }

    // compute cotan weights
    for (CMeshO::FaceIterator fi = mesh.face.begin(); fi != mesh.face.end(); ++fi){
        vcg::Point3f p0 = fi->V(0)->P();
        vcg::Point3f p1 = fi->V(1)->P();
        vcg::Point3f p2 = fi->V(2)->P();

        Eigen::Vector3d e0 = toEigen(p2 - p1);
        Eigen::Vector3d e1 = toEigen(p0 - p2);
        Eigen::Vector3d e2 = toEigen(p1 - p0);

        // first edge is inverted to get correct orientation
        double alpha0 = cotan(-e1, e2) / 2;
        double alpha1 = cotan(-e2, e0) / 2;
        double alpha2 = cotan(-e0, e1) / 2;

        int i0 = vertex_ids[fi->V(0)];
        int i1 = vertex_ids[fi->V(1)];
        int i2 = vertex_ids[fi->V(2)];

        // save only lower triangular part
        if (i0 > i1)
            cotanOperator.coeffRef(i0, i1) += alpha2;
        else
            cotanOperator.coeffRef(i1, i0) += alpha2;
        if (i0 > i2)
            cotanOperator.coeffRef(i0, i2) += alpha1;
        else
            cotanOperator.coeffRef(i2, i0) += alpha1;
        if (i1 > i2)
            cotanOperator.coeffRef(i1, i2) += alpha0;
        else
            cotanOperator.coeffRef(i2, i1) += alpha0;
        
        cotanOperator.coeffRef(i0, i0) -= (alpha1 + alpha2);
        cotanOperator.coeffRef(i1, i1) -= (alpha0 + alpha2);
        cotanOperator.coeffRef(i2, i2) -= (alpha0 + alpha1);
    }
}


inline double computeAverageEdgeLength(CMeshO &mesh){
    // compute total length of all edges
    double total_length = 0;
    for (CMeshO::FaceIterator fi = mesh.face.begin(); fi != mesh.face.end(); ++fi)
    {
        vcg::Point3f p0 = fi->V(0)->P();
        vcg::Point3f p1 = fi->V(1)->P();
        vcg::Point3f p2 = fi->V(2)->P();
        double e2 = toEigen(p1 - p0).norm();
        double e1 = toEigen(p2 - p0).norm();
        double e0 = toEigen(p2 - p1).norm();
        total_length += (e0 + e1 + e2) / 2;
    }
    return total_length / (3./2. * mesh.FN());
}


inline Eigen::MatrixX3d computeVertexGradient(CMeshO &mesh, const Eigen::VectorXd &heat){
    Eigen::MatrixX3d heatGradientField(mesh.FN(), 3);
    // initialize a hashtable from vertex pointers to ids
    std::unordered_map<CMeshO::VertexType*, int> vertex_ids;
    for (int i = 0; i < mesh.VN(); ++i){
        vertex_ids[&mesh.vert[i]] = i;
    }
    // compute gradient of heat function at each vertex
    for (int i = 0; i < mesh.FN(); ++i){
        CMeshO::FaceType *fp = &mesh.face[i];

        vcg::Point3f p0 = fp->V(0)->P();
        vcg::Point3f p1 = fp->V(1)->P();
        vcg::Point3f p2 = fp->V(2)->P();

        // normal unit vector
        Eigen::Vector3d n = toEigen(fp->N());
        n /= n.norm();
        // face area
        double faceArea = fp->Q();
        // edge unit vectors (counter-clockwise)
        Eigen::Vector3d e0 = toEigen(p2 - p1);
        e0 /= e0.norm();
        Eigen::Vector3d e1 = toEigen(p0 - p2);
        e1 /= e1.norm();
        Eigen::Vector3d e2 = toEigen(p1 - p0);
        e2 /= e2.norm();
        // gradient unit vectors
        Eigen::Vector3d g0 = n.cross(e0); //v0 grad
        Eigen::Vector3d g1 = n.cross(e1); //v1 grad
        Eigen::Vector3d g2 = n.cross(e2); //v2 grad

        // add vertex gradient contributions
        Eigen::Vector3d tri_grad = (
            g0 * heat(vertex_ids[fp->V(0)]) + 
            g1 * heat(vertex_ids[fp->V(1)]) + 
            g2 * heat(vertex_ids[fp->V(2)])
        ) / (2 * faceArea);

        heatGradientField.row(i) = tri_grad;
    }
    return heatGradientField;
}


inline Eigen::MatrixX3d normalizeVectorField(const Eigen::MatrixX3d &field){
    Eigen::MatrixX3d normalizedField(field.rows(), 3);
    normalizedField.setZero();
    // normalize vector field at each vertex
    for (int i = 0; i < field.rows(); ++i){
        Eigen::Vector3d v = field.row(i);
        normalizedField.row(i) = v / v.norm();
    }
    return normalizedField;
}


inline Eigen::VectorXd computeVertexDivergence(CMeshO &mesh, const Eigen::MatrixX3d &field){
    Eigen::VectorXd divergence(mesh.VN());
    divergence.setZero();
    // initialize a hashtable from face pointers to ids
    std::unordered_map<CMeshO::FaceType*, int> face_ids;
    for (int i = 0; i < mesh.FN(); ++i){
        face_ids[&mesh.face[i]] = i;
    }

    // compute divergence of vector field at each vertex
    for (int i = 0; i < mesh.VN(); ++i){
        CMeshO::VertexType *vp = &mesh.vert[i];

        std::vector<CMeshO::FaceType*> faces;
        std::vector<int> indices;
        vcg::face::VFStarVF<CMeshO::FaceType>(vp, faces, indices);
        for (int j = 0; j < faces.size(); ++j)
        {
            CMeshO::FaceType *fp = faces[j];
            int index = indices[j];
            vcg::Point3f p0 = fp->V(0)->P();
            vcg::Point3f p1 = fp->V(1)->P();
            vcg::Point3f p2 = fp->V(2)->P();
            // (ORDERING) edge vectors
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
            // normalize edge vectors after cotangent computation
            el /= el.norm();
            er /= er.norm();
            // add divergence contribution of given face
            Eigen::Vector3d x = field.row(face_ids[fp]);
            divergence(i) += (cotl * er.dot(x) + cotr * el.dot(x)) / 2;
        }
    }
    return divergence;
}


enum SaveMask : unsigned char
{
	NONE                 = 0x00,

    CSV_HEAT_FLOW         = 0x01,
    CSV_HEAT_GRADIENT     = 0x02,
    CSV_UNIT_VECTOR_FIELD = 0x04,
    CSV_DIVERGENCE        = 0x08,
    CSV_DISTANCE          = 0x10,

    ALL                  = 0x1F,
};

inline SaveMask operator|(SaveMask lhs, SaveMask rhs) {
    return static_cast<SaveMask>(
        static_cast<std::underlying_type<SaveMask>::type>(lhs) |
        static_cast<std::underlying_type<SaveMask>::type>(rhs)
    );
}
inline SaveMask operator&(SaveMask lhs, SaveMask rhs) {
    return static_cast<SaveMask>(
        static_cast<std::underlying_type<SaveMask>::type>(lhs) &
        static_cast<std::underlying_type<SaveMask>::type>(rhs)
    );
}


// DEBUGGING FUNCTIONS

inline void saveVertexScalarFieldCSV(CMeshO &mesh, const Eigen::VectorXd scalarField, const char *fname){
    // save vertex positions with scalar field values into CSV file
    std::ofstream file(fname);
    file << "x,y,z,s" << std::endl;
    for (int i = 0; i < mesh.VN(); ++i){
        Eigen::Vector3d point = toEigen(mesh.vert[i].P());
        file << point(0) << "," << point(1) << "," << point(2) << "," << scalarField(i) << std::endl;
    }
}

inline void saveFaceVectorFieldCSV(CMeshO &mesh, const Eigen::MatrixX3d vectorField, const char *fname){
    // save face barycenter positions with vector field values into CSV file
    std::ofstream file(fname);
    file << "px,py,pz,vx,vy,vz" << std::endl;
    for (int i = 0; i < mesh.FN(); ++i){
        auto face = mesh.face[i];
        Eigen::Vector3d bar = toEigen(face.P(0) + face.P(1) + face.P(2)) / 3;
        Eigen::Vector3d vec = vectorField.row(i);
        file << bar(0) << "," << bar(1) << "," << bar(2) << "," << vec(0) << "," << vec(1) << "," << vec(2) << std::endl;
    }
}

inline Eigen::VectorXd computeHeatMethodGeodesic(
    CMeshO &mesh, 
    const Eigen::VectorXd &init_cond, 
    double m = 1.0, 
    SaveMask save = SaveMask::NONE
){
    vcg::tri::UpdateTopology<CMeshO>::VertexFace(mesh);
    vcg::tri::UpdateTopology<CMeshO>::FaceFace(mesh);
    vcg::tri::UpdateNormal<CMeshO>::PerFaceNormalized(mesh);

    Eigen::SparseMatrix<double> mass(mesh.VN(), mesh.VN());
    buildMassMatrix(mesh, mass);

    Eigen::SparseMatrix<double> cotanOperator(mesh.VN(), mesh.VN());
    buildCotanLowerTriMatrix(mesh, cotanOperator);

    double avg_edge_len = computeAverageEdgeLength(mesh);
    double timestep = m * avg_edge_len * avg_edge_len;
    Eigen::SparseMatrix<double> system1(mesh.VN(), mesh.VN());
    system1 = mass - timestep * cotanOperator;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    solver.compute(system1);
    if(solver.info() != Eigen::Success) {
        std::cerr << "ERROR: Cholesky Factorization 1 failed" << std::endl; 
    }
    Eigen::VectorXd heatflow = solver.solve(init_cond); // (VN)
    if(solver.info() != Eigen::Success) {
        std::cerr << "ERROR: Solving System 1 failed" << std::endl; 
    }
    if(save & SaveMask::CSV_HEAT_FLOW) 
        saveVertexScalarFieldCSV(mesh, heatflow, "output_csv/1_heatflow.csv");

    Eigen::MatrixX3d heatGradient = computeVertexGradient(mesh, heatflow); // (FN, 3)
    if(save & SaveMask::CSV_HEAT_GRADIENT) 
        saveFaceVectorFieldCSV(mesh, heatGradient, "output_csv/2_heatGradient.csv");

    Eigen::MatrixX3d normalizedVectorField = normalizeVectorField(-heatGradient); // (FN, 3)
    if(save & SaveMask::CSV_UNIT_VECTOR_FIELD)
        saveFaceVectorFieldCSV(mesh, normalizedVectorField, "output_csv/3_normalizedVectorField.csv");
    
    Eigen::VectorXd divergence = computeVertexDivergence(mesh, normalizedVectorField); // (VN)
    if(save & SaveMask::CSV_DIVERGENCE) 
        saveVertexScalarFieldCSV(mesh, divergence, "output_csv/4_divergence.csv");

    Eigen::SparseMatrix<double> system2(mesh.VN(), mesh.VN());
    system2 = cotanOperator;
    
    solver.compute(system2);
    if(solver.info() != Eigen::Success) {
        std::cerr << "ERROR: Cholesky Factorization 2 failed" << std::endl; 
    }
    Eigen::VectorXd geodesicDistance = solver.solve(divergence);
    if(solver.info() != Eigen::Success) {
        std::cerr << "ERROR: Solving System 2 failed" << std::endl; 
    }
    // shift to impose boundary conditions (dist(d) = 0 \forall d \in init_cond)
    geodesicDistance.array() -= geodesicDistance.minCoeff();

    if(save & SaveMask::CSV_DISTANCE) 
        saveVertexScalarFieldCSV(mesh, geodesicDistance, "output_csv/5_distance.csv");

    return geodesicDistance;
}


#endif