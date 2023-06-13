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
    vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);

    std::unordered_map<MyVertex*, int> vertex_ids;
    for (int i = 0; i < m.VN(); ++i){
        vertex_ids[&m.vert[i]] = i;
    }

    MyVertex *vp = &m.vert[std::rand() % m.VN()];
    MyFace *fp = vp->VFp();
    vcg::face::Pos<MyFace> pos(fp, vp);
    vcg::face::Pos<MyFace> start(fp, vp);

    // iterate over all vertices to fill cotan matrix
    int total_triplets = 0;
    for (int i = 0; i < m.VN(); ++i){
        MyVertex *vp = &m.vert[i];
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
            std::cout << "[ " << vertex_ids[vp] << ", " << vertex_ids[vo] << " ] : " << (cotan_l + cotan_r)/2 << std::endl;
            ++total_triplets;

            // move to the next edge
            pos.FlipF();pos.FlipE();
        } while (pos != start);
    }

    // print mesh info
    std::cout << "# Verts: " << m.VN() << std::endl;
    std::cout << "# Edges: " << m.VN() + m.FN() - 2 << std::endl;
    std::cout << "Total triplets: " << total_triplets << std::endl;

    return 0;
}
