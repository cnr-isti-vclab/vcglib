#ifndef TRIMESH_HEAT_METHOD
#define TRIMESH_HEAT_METHOD

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/quality.h>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <vector>
#include <cmath>
#include <unordered_map>

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
    vcg::vertex::Color4b, 
    vcg::vertex::Qualityd,
    vcg::vertex::BitFlags // needed for PLY export
>{};
class MyEdge : public vcg::Edge<
    MyUsedTypes,
    vcg::edge::VertexRef //
>{};
class MyFace : public vcg::Face<
    MyUsedTypes,
    vcg::face::VFAdj,
    vcg::face::FFAdj,
    vcg::face::VertexRef,
    vcg::face::Qualityd
>{};
class MyMesh : public vcg::tri::TriMesh<
    std::vector<MyVertex>, 
    std::vector<MyFace>,
    std::vector<MyEdge>
>{};


inline Eigen::Vector3d toEigen(const vcg::Point3f& p)
{
    return Eigen::Vector3d(p.X(), p.Y(), p.Z());
};


inline double cotan(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1)
{
    return v0.dot(v1) / v0.cross(v1).norm();
};


inline void buildMassMatrix(MyMesh &mesh, Eigen::SparseMatrix<double> &mass){
    // compute area of all faces
    for (MyMesh::FaceIterator fi = mesh.face.begin(); fi != mesh.face.end(); ++fi)
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
        MyVertex *vp = &mesh.vert[i];

        std::vector<MyFace*> faces;
        std::vector<int> indices;
        vcg::face::VFStarVF<MyFace>(vp, faces, indices);
        
        double area = 0;
        for (int j = 0; j < faces.size(); ++j)
        {
            area += faces[j]->Q();
        }
        area /= 3;
        mass.coeffRef(i, i) = area;
    }
}


inline void buildCotanMatrix(MyMesh &mesh, Eigen::SparseMatrix<double> &cotanOperator){
    // initialize a hashtable from vertex pointers to ids
    std::unordered_map<MyVertex*, int> vertex_ids;
    for (int i = 0; i < mesh.VN(); ++i){
        vertex_ids[&mesh.vert[i]] = i;
    }

    // iterate over all vertices to fill cotan matrix
    for (int i = 0; i < mesh.VN(); ++i){
        MyVertex *vp = &mesh.vert[i];
        MyFace *fp = vp->VFp();
        vcg::face::Pos<MyFace> pos(fp, vp);
        vcg::face::Pos<MyFace> start(fp, vp);
        // iterate over all incident edges of vp
        do {
            // get the vertex opposite to vp
            pos.FlipV();
            MyVertex *vo = pos.V();
            // move to left vertex
            pos.FlipE();pos.FlipV();
            MyVertex *vl = pos.V();
            // move back then to right vertex
            pos.FlipV();pos.FlipE(); // back to vo
            pos.FlipF();pos.FlipE();pos.FlipV();
            MyVertex *vr = pos.V();
            pos.FlipV();pos.FlipE();pos.FlipF();pos.FlipV(); // back to vp

            // compute cotan of left edges and right edges
            Eigen::Vector3d elf = toEigen(vo->P() - vl->P()); // far left edge
            Eigen::Vector3d eln = toEigen(vp->P() - vl->P()); // near left edge
            Eigen::Vector3d erf = toEigen(vp->P() - vr->P()); // far right edge
            Eigen::Vector3d ern = toEigen(vo->P() - vr->P()); // near right edge

            double cotan_l = cotan(elf, eln);
            double cotan_r = cotan(ern, erf);

            // add to the matrix
            cotanOperator.coeffRef(vertex_ids[vp], vertex_ids[vo]) = (cotan_l + cotan_r)/2;

            // move to the next edge
            pos.FlipF();pos.FlipE();
        } while (pos != start);
    }

    // compute diagonal entries
    for (int i = 0; i < mesh.VN(); ++i){
        cotanOperator.coeffRef(i, i) = -cotanOperator.row(i).sum();
    }
    
}


inline double computeAverageEdgeLength(MyMesh &mesh){
    // for simplicity we use edge vertex references
    double avg_edge_len = 0;
    for (MyMesh::EdgeIterator ei = mesh.edge.begin(); ei != mesh.edge.end(); ++ei)
    {
        vcg::Point3f p0 = ei->V(0)->P();
        vcg::Point3f p1 = ei->V(1)->P();
        avg_edge_len += toEigen(p1 - p0).norm();
    }
    avg_edge_len /= mesh.edge.size();
    return avg_edge_len;
}


inline Eigen::MatrixX3d coputeVertexGradient(MyMesh &mesh, const Eigen::VectorXd &heat){
    Eigen::MatrixX3d heatGradientField(mesh.VN(), 3);
    heatGradientField.setZero();
    // compute gradient of heat function at each vertex
    for (int i = 0; i < mesh.VN(); ++i){
        MyVertex *vp = &mesh.vert[i];

        std::vector<MyFace*> faces;
        std::vector<int> indices;
        vcg::face::VFStarVF<MyFace>(vp, faces, indices);
        for (int j = 0; j < faces.size(); ++j)
        {
            MyFace *fp = faces[j];
            int index = indices[j];
            vcg::Point3f p0 = fp->V(0)->P();
            vcg::Point3f p1 = fp->V(1)->P();
            vcg::Point3f p2 = fp->V(2)->P();
            // edge unit vector
            Eigen::Vector3d e;
            // assuming counter-clockwise ordering
            if (index == 0){
                e = toEigen(p2 - p1); //e0
            } else if (index == 1){
                e = toEigen(p0 - p2); //e1
            } else if (index == 2){ 
                e = toEigen(p1 - p0); //e2
            }
            e /= e.norm();
            // normal unit vector
            Eigen::Vector3d n = toEigen(fp->N());
            n /= n.norm();
            // gradient unit vector
            Eigen::Vector3d g = n.cross(e);
            // add gradient contribution of given face
            double faceArea = fp->Q();
            heatGradientField.row(i) += g * (heat[i] / (2 * faceArea));
        }
    }
    return heatGradientField;
}

inline Eigen::MatrixX3d normalizeVectorField(MyMesh &mesh, const Eigen::MatrixX3d &field){
    Eigen::MatrixX3d normalizedField(mesh.VN(), 3);
    normalizedField.setZero();
    // normalize vector field at each vertex
    for (int i = 0; i < mesh.VN(); ++i){
        Eigen::Vector3d v = field.row(i);
        normalizedField.row(i) = v / v.norm();
    }
    return normalizedField;
}


inline Eigen::VectorXd computeVertexDivergence(MyMesh &mesh, const Eigen::MatrixX3d &field){
    Eigen::VectorXd divergence(mesh.VN());
    divergence.setZero();
    // compute divergence of vector field at each vertex
    for (int i = 0; i < mesh.VN(); ++i){
        MyVertex *vp = &mesh.vert[i];

        std::vector<MyFace*> faces;
        std::vector<int> indices;
        vcg::face::VFStarVF<MyFace>(vp, faces, indices);
        for (int j = 0; j < faces.size(); ++j)
        {
            MyFace *fp = faces[j];
            int index = indices[j];
            vcg::Point3f p0 = fp->V(0)->P();
            vcg::Point3f p1 = fp->V(1)->P();
            vcg::Point3f p2 = fp->V(2)->P();
            // edge vectors
            Eigen::Vector3d el, er, eo; //left, right, opposite
            if (index == 0){
                el = toEigen(p2 - p0); //e1
                er = toEigen(p1 - p0); //e2
                eo = toEigen(p1 - p2); //[+-] e0
            } else if (index == 1){
                el = toEigen(p0 - p1); //e2
                er = toEigen(p2 - p1); //e0
                eo = toEigen(p0 - p2); //[+-] e1
            } else if (index == 2){ 
                el = toEigen(p1 - p2); //e0
                er = toEigen(p0 - p2); //e1
                eo = toEigen(p0 - p1); //[+-] e2
            }
            // compute left and right cotangents
            double cotl = cotan(el, eo);
            double cotr = cotan(er, eo);
            // normalize edge vectors after cotangent computation
            el /= el.norm();
            er /= er.norm();
            // add divergence contribution of given face
            divergence(i) += (cotl * er.dot(field.row(i)) + cotr * el.dot(field.row(i))) / 2;
        }
    }
    return divergence;
}


inline Eigen::VectorXd computeHeatMethodGeodesic(MyMesh &mesh, const Eigen::VectorXd &init_cond, double m = 1){
    vcg::tri::UpdateTopology<MyMesh>::VertexFace(mesh);
    vcg::tri::UpdateNormal<MyMesh>::PerFaceNormalized(mesh);

    Eigen::SparseMatrix<double> mass(mesh.VN(), mesh.VN());
    buildMassMatrix(mesh, mass);
    Eigen::SparseMatrix<double> cotanOperator(mesh.VN(), mesh.VN());
    buildCotanMatrix(mesh, cotanOperator);

    double avg_edge_len = computeAverageEdgeLength(mesh);
    double timestep = m * avg_edge_len * avg_edge_len;
    Eigen::SparseMatrix<double> system(mesh.VN(), mesh.VN());
    system = mass - timestep * cotanOperator;

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> cholesky1(system);
    Eigen::VectorXd heatflow = cholesky1.solve(init_cond);

    Eigen::MatrixX3d heatGradient = coputeVertexGradient(mesh, heatflow);

    Eigen::MatrixX3d normalizedVectorField = normalizeVectorField(mesh, -heatGradient);
    
    Eigen::VectorXd divergence = computeVertexDivergence(mesh, normalizedVectorField);

    // maybe precondition the matrix by adding a small multiple of the identity matrix
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> cholesky2(cotanOperator);
    Eigen::VectorXd geodesicDistance = cholesky2.solve(divergence);

    return geodesicDistance;
}


#endif