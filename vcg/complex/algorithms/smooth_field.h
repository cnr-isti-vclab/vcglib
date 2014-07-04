#ifndef SMOOTHER_FIELD_H
#define SMOOTHER_FIELD_H

#include <vector>
#include <list>
#include <utility>
#include "mesh_to_matrix.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/curvature_fitting.h>

//#include "NRosyField.h"

using namespace std;
#define Delta 10e-6

template < typename TriMeshType >
class ImplicitSmoother
{
    // Temporary variable for the field
    Eigen::VectorXd angles;

    // Hard constraints
    Eigen::VectorXd hard;
    std::vector<bool> isHard;

    // Soft constraints
    Eigen::VectorXd soft;
    Eigen::VectorXd wSoft;
    double          softAlpha;

    // Face Topology
    Eigen::MatrixXi TT, TTi;

    // Edge Topology
    Eigen::MatrixXi EV, FE, EF;
    std::vector<bool> isBorderEdge;

    // Per Edge information
    // Angle between two reference frames
    Eigen::VectorXd k;

    // Jumps
    Eigen::VectorXi p;
    std::vector<bool> pFixed;

    // Mesh
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    // Normals per face
    Eigen::MatrixXd N;

    // Singularity index
    Eigen::VectorXd singularityIndex;

    // Reference frame per triangle
    std::vector<Eigen::MatrixXd> TPs;

    // System stuff
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    Eigen::VectorXi tag_t;
    Eigen::VectorXi tag_p;

    // define types
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename TriMeshType::ScalarType ScalarType;

    TriMeshType &mesh;

    void computek()
    {
        // For every non-border edge
        for (unsigned eid=0; eid<EF.rows(); ++eid)
        {
            if (!isBorderEdge[eid])
            {
                int fid0 = EF(eid,0);
                int fid1 = EF(eid,1);

                Eigen::Vector3d N0 = N.row(fid0);
                Eigen::Vector3d N1 = N.row(fid1);

                // find common edge on triangle 0 and 1
                int fid0_vc = -1;
                int fid1_vc = -1;
                for (unsigned i=0;i<3;++i)
                {
                    if (EV(eid,0) == F(fid0,i))
                        fid0_vc = i;
                    if (EV(eid,1) == F(fid1,i))
                        fid1_vc = i;
                }
                assert(fid0_vc != -1);
                assert(fid1_vc != -1);

                Eigen::Vector3d common_edge = V.row(F(fid0,(fid0_vc+1)%3)) - V.row(F(fid0,fid0_vc));
                common_edge.normalize();

                // Map the two triangles in a new space where the common edge is the x axis and the N0 the z axis
                Eigen::MatrixXd P(3,3);
                Eigen::VectorXd o = V.row(F(fid0,fid0_vc));

                Eigen::VectorXd tmp = -N0.cross(common_edge);
                P << common_edge, tmp, N0;
                P.transposeInPlace();


                Eigen::MatrixXd V0(3,3);
                V0.row(0) = V.row(F(fid0,0)).transpose() -o;
                V0.row(1) = V.row(F(fid0,1)).transpose() -o;
                V0.row(2) = V.row(F(fid0,2)).transpose() -o;

                V0 = (P*V0.transpose()).transpose();

                assert(V0(0,2) < Delta);
                assert(V0(1,2) < Delta);
                assert(V0(2,2) < Delta);

                Eigen::MatrixXd V1(3,3);

                V1.row(0) = V.row(F(fid1,0)).transpose() -o;
                V1.row(1) = V.row(F(fid1,1)).transpose() -o;
                V1.row(2) = V.row(F(fid1,2)).transpose() -o;

                V1 = (P*V1.transpose()).transpose();

                assert(V1(fid1_vc,2) < Delta);
                assert(V1((fid1_vc+1)%3,2) < Delta);

                // compute rotation R such that R * N1 = N0
                // i.e. map both triangles to the same plane
                double alpha = -atan2(V1((fid1_vc+2)%3,2),V1((fid1_vc+2)%3,1));

                Eigen::MatrixXd R(3,3);
                R << 1,          0,            0,
                        0, cos(alpha), -sin(alpha) ,
                        0, sin(alpha),  cos(alpha);
                V1 = (R*V1.transpose()).transpose();

                assert(V1(0,2) < Delta);
                assert(V1(1,2) < Delta);
                assert(V1(2,2) < Delta);

                // measure the angle between the reference frames
                // k_ij is the angle between the triangle on the left and the one on the right
                Eigen::VectorXd ref0 = V0.row(1) - V0.row(0);
                Eigen::VectorXd ref1 = V1.row(1) - V1.row(0);

                ref0.normalize();
                ref1.normalize();

                double ktemp = atan2(ref1(1),ref1(0)) - atan2(ref0(1),ref0(0));

                // just to be sure, rotate ref0 using angle ktemp...
                Eigen::MatrixXd R2(2,2);
                R2 << cos(ktemp), -sin(ktemp), sin(ktemp), cos(ktemp);

                tmp = R2*ref0.head<2>();

                assert(tmp(0) - ref1(0) < 10^10);
                assert(tmp(1) - ref1(1) < 10^10);

                k[eid] = ktemp;
            }
        }
    }

    void resetConstraints()
    {
        isHard.resize(F.rows());
        for(unsigned i=0; i<F.rows(); ++i)
            isHard[i] = false;
        hard   = Eigen::VectorXd::Zero(F.rows());

        wSoft  = Eigen::VectorXd::Zero(F.rows());
        soft   = Eigen::VectorXd::Zero(F.rows());
    }

    void initializeSmoother()
    {

        // Generate topological relations
        vcg::MeshToMatrix<TriMeshType>::GetTriMeshData(mesh,F,V);
        vcg::MeshToMatrix<TriMeshType>::GetTriFFAdjacency(mesh,TT,TTi);
        vcg::MeshToMatrix<TriMeshType>::GetTriEdgeAdjacency(mesh,EV,FE,EF);

        // Flag border edges
        isBorderEdge.resize(EV.rows());
        for(unsigned i=0; i<EV.rows(); ++i)
            isBorderEdge[i] = (EF(i,0) == -1) || ((EF(i,1) == -1));

        // Generate normals per face
        //igl::per_face_normals(V, F, N);
        Eigen::MatrixXd NV;
        vcg::MeshToMatrix<TriMeshType>::GetNormalData(mesh,NV,N);

        // Generate reference frames
        for(unsigned fid=0; fid<F.rows(); ++fid)
        {
            // First edge
            Eigen::Vector3d e1 = V.row(F(fid,1)) - V.row(F(fid,0));
            e1.normalize();
            Eigen::Vector3d e2 = N.row(fid);
            e2 = e2.cross(e1);
            e2.normalize();

            Eigen::MatrixXd TP(2,3);
            TP << e1.transpose(), e2.transpose();
            TPs.push_back(TP);
        }

        // Alloc internal variables
        angles = Eigen::VectorXd::Zero(F.rows());
        p = Eigen::VectorXi::Zero(EV.rows());
        pFixed.resize(EV.rows());
        k = Eigen::VectorXd::Zero(EV.rows());
        singularityIndex = Eigen::VectorXd::Zero(V.rows());

        // Reset the constraints
        resetConstraints();

        // Compute k, differences between reference frames
        computek();

        softAlpha = 0.5;
    }

    double convert3DtoLocal(unsigned fid, const Eigen::Vector3d& v)
    {
        // Project onto the tangent plane
        Eigen::Vector2d vp = TPs[fid] * v;

        // Convert to angle
        return atan2(vp(1),vp(0));
    }

    Eigen::Vector3d convertLocalto3D(unsigned fid, double a)
    {
        Eigen::Vector2d vp(cos(a),sin(a));
        return vp.transpose() * TPs[fid];
    }

    void reduceSpace()
    {
        // All variables are free in the beginning
        for(unsigned i=0; i<EV.rows(); ++i)
            pFixed[i] = false;

        vector<Eigen::VectorXd> debug;

        // debug
        //  MatrixXd B(F.rows(),3);
        //  for(unsigned i=0; i<F.rows(); ++i)
        //    B.row(i) = 1./3. * (V.row(F(i,0)) + V.row(F(i,1)) + V.row(F(i,2)));

        vector<bool> visited(EV.rows());
        for(unsigned i=0; i<EV.rows(); ++i)
            visited[i] = false;

        vector<bool> starting(EV.rows());
        for(unsigned i=0; i<EV.rows(); ++i)
            starting[i] = false;

        queue<int> q;
        for(unsigned i=0; i<F.rows(); ++i)
            if (isHard[i] || wSoft[i] != 0)
            {
                q.push(i);
                starting[i] = true;
            }

        // Reduce the search space (see MI paper)
        while (!q.empty())
        {
            int c = q.front();
            q.pop();

            visited[c] = true;
            for(int i=0; i<3; ++i)
            {
                int eid = FE(c,i);
                int fid = TT(c,i);

                // skip borders
                if (fid != -1)
                {
                    assert((EF(eid,0) == c && EF(eid,1) == fid) || (EF(eid,1) == c && EF(eid,0) == fid));
                    // for every neighbouring face
                    if (!visited[fid] && !starting[fid])
                    {
                        pFixed[eid] = true;
                        p[eid] = 0;
                        visited[fid] = true;
                        q.push(fid);

                    }
                }
                else
                {
                    // fix borders
                    pFixed[eid] = true;
                    p[eid] = 0;
                }
            }

        }

        // Force matchings between fixed faces
        for(unsigned i=0; i<F.rows();++i)
        {
            if (isHard[i])
            {
                for(unsigned int j=0; j<3; ++j)
                {
                    int fid = TT(i,j);
                    if ((fid!=-1) && (isHard[fid]))
                    {
                        // i and fid are adjacent and fixed
                        int eid = FE(i,j);
                        int fid0 = EF(eid,0);
                        int fid1 = EF(eid,1);

                        pFixed[eid] = true;
                        p[eid] = roundl(2.0/M_PI*(hard(fid1) - hard(fid0) - k(eid)));
                    }
                }
            }
        }

        //  std::ofstream s("./debug.txt");
        //  for(unsigned i=0; i<debug.size(); i += 2)
        //    s << debug[i].transpose() << " " << debug[i+1].transpose() << endl;
        //  s.close();

    }

    void prepareSystemMatrix(const int N)
    {
        double Nd = N;

        // Minimize the MIQ energy
        // Energy on edge ij is
        //     (t_i - t_j + kij + pij*(2*pi/N))^2
        // Partial derivatives:
        //   t_i: 2     ( t_i - t_j + kij + pij*(2*pi/N)) = 0
        //   t_j: 2     (-t_i + t_j - kij - pij*(2*pi/N)) = 0
        //   pij: 4pi/N ( t_i - t_j + kij + pij*(2*pi/N)) = 0
        //
        //          t_i      t_j         pij       kij
        // t_i [     2       -2           4pi/N      2    ]
        // t_j [    -2        2          -4pi/N     -2    ]
        // pij [   4pi/N   -4pi/N    2*(2pi/N)^2   4pi/N  ]

        // Count and tag the variables
        tag_t = Eigen::VectorXi::Constant(F.rows(),-1);
        vector<int> id_t;
        int count = 0;
        for(unsigned i=0; i<F.rows(); ++i)
            if (!isHard[i])
            {
                tag_t(i) = count++;
                id_t.push_back(i);
            }

        unsigned count_t = id_t.size();

        tag_p = Eigen::VectorXi::Constant(EF.rows(),-1);
        vector<int> id_p;
        for(unsigned i=0; i<EF.rows(); ++i)
        {
            if (!pFixed[i])
            {
                // if it is not fixed then it is a variable
                tag_p(i) = count++;
            }

            // if it is not a border edge,
            if (!isBorderEdge[i])
            {
                // and it is not between two fixed faces
                if (!(isHard[EF(i,0)] && isHard[EF(i,1)]))
                {
                    // then it participates in the energy!
                    id_p.push_back(i);
                }
            }
        }

        unsigned count_p = count - count_t;
        // System sizes: A (count_t + count_p) x (count_t + count_p)
        //               b (count_t + count_p)

        b = Eigen::VectorXd::Zero(count_t + count_p);

        std::vector<Eigen::Triplet<double> > T;
        T.reserve(3 * 4 * count_p);

        for(unsigned r=0; r<id_p.size(); ++r)
        {
            int eid = id_p[r];
            int i = EF(eid,0);
            int j = EF(eid,1);
            bool isFixed_i = isHard[i];
            bool isFixed_j = isHard[j];
            bool isFixed_p = pFixed[eid];
            int row;
            // (i)-th row: t_i [     2       -2           4pi/N      2    ]
            if (!isFixed_i)
            {
                row = tag_t[i];
                if (isFixed_i) b(row) += -2               * hard[i]; else T.push_back(Eigen::Triplet<double>(row,tag_t[i]  , 2             ));
                if (isFixed_j) b(row) +=  2               * hard[j]; else T.push_back(Eigen::Triplet<double>(row,tag_t[j]  ,-2             ));
                if (isFixed_p) b(row) += -((4 * M_PI)/Nd) * p[eid] ; else T.push_back(Eigen::Triplet<double>(row,tag_p[eid],((4 * M_PI)/Nd)));
                b(row) += -2 * k[eid];
                assert(hard[i] == hard[i]);
                assert(hard[j] == hard[j]);
                assert(p[eid] == p[eid]);
                assert(k[eid] == k[eid]);
                assert(b(row) == b(row));
            }
            // (j)+1 -th row: t_j [    -2        2          -4pi/N     -2    ]
            if (!isFixed_j)
            {
                row = tag_t[j];
                if (isFixed_i) b(row) += 2               * hard[i]; else T.push_back(Eigen::Triplet<double>(row,tag_t[i]  , -2             ));
                if (isFixed_j) b(row) += -2              * hard[j]; else T.push_back(Eigen::Triplet<double>(row,tag_t[j] ,  2              ));
                if (isFixed_p) b(row) += ((4 * M_PI)/Nd) * p[eid] ; else T.push_back(Eigen::Triplet<double>(row,tag_p[eid],-((4 * M_PI)/Nd)));
                b(row) += 2 * k[eid];
                assert(k[eid] == k[eid]);
                assert(b(row) == b(row));
            }
            // (r*3)+2 -th row: pij [   4pi/N   -4pi/N    2*(2pi/N)^2   4pi/N  ]
            if (!isFixed_p)
            {
                row = tag_p[eid];
                if (isFixed_i) b(row) += -(4 * M_PI)/Nd              * hard[i]; else T.push_back(Eigen::Triplet<double>(row,tag_t[i] ,   (4 * M_PI)/Nd             ));
                if (isFixed_j) b(row) +=  (4 * M_PI)/Nd              * hard[j]; else T.push_back(Eigen::Triplet<double>(row,tag_t[j] ,  -(4 * M_PI)/Nd             ));
                if (isFixed_p) b(row) += -(2 * pow(((2*M_PI)/Nd),2)) * p[eid] ;  else T.push_back(Eigen::Triplet<double>(row,tag_p[eid],  (2 * pow(((2*M_PI)/Nd),2))));
                b(row) += - (4 * M_PI)/Nd * k[eid];
                assert(k[eid] == k[eid]);
                assert(b(row) == b(row));
            }

        }

        A = Eigen::SparseMatrix<double>(count_t + count_p, count_t + count_p);
        A.setFromTriplets(T.begin(), T.end());

        // Soft constraints
        bool addSoft = false;

        for(unsigned i=0; i<wSoft.size();++i)
            if (wSoft[i] != 0)
                addSoft = true;

        if (addSoft)
        {
            cerr << " Adding soft here: " << endl;
            cerr << " softAplha: " << softAlpha << endl;
            Eigen::VectorXd bSoft = Eigen::VectorXd::Zero(count_t + count_p);

            std::vector<Eigen::Triplet<double> > TSoft;
            TSoft.reserve(2 * count_p);

            for(unsigned i=0; i<F.rows(); ++i)
            {
                int varid = tag_t[i];
                if (varid != -1) // if it is a variable in the system
                {
                    TSoft.push_back(Eigen::Triplet<double>(varid,varid,wSoft[i]));
                    bSoft[varid] += wSoft[i] * soft[i];
                }
            }
            Eigen::SparseMatrix<double> ASoft(count_t + count_p, count_t + count_p);
            ASoft.setFromTriplets(TSoft.begin(), TSoft.end());

            //    ofstream s("./As.txt");
            //    for(unsigned i=0; i<TSoft.size(); ++i)
            //      s << TSoft[i].row() << " " << TSoft[i].col() << " " << TSoft[i].value() << endl;
            //    s.close();

            //    ofstream s2("./bs.txt");
            //    for(unsigned i=0; i<bSoft.rows(); ++i)
            //      s2 << bSoft(i) << endl;
            //    s2.close();

            // Stupid Eigen bug
            Eigen::SparseMatrix<double> Atmp (count_t + count_p, count_t + count_p);
            Eigen::SparseMatrix<double> Atmp2(count_t + count_p, count_t + count_p);
            Eigen::SparseMatrix<double> Atmp3(count_t + count_p, count_t + count_p);

            // Merge the two part of the energy
            Atmp = (1.0 - softAlpha)*A;
            Atmp2 = softAlpha * ASoft;
            Atmp3 = Atmp+Atmp2;

            A = Atmp3;
            b = b*(1.0 - softAlpha) + bSoft * softAlpha;
        }

        //  ofstream s("./A.txt");
        //  for (int k=0; k<A.outerSize(); ++k)
        //    for (SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
        //    {
        //      s << it.row() << " " << it.col() << " " << it.value() << endl;
        //    }
        //  s.close();
        //
        //  ofstream s2("./b.txt");
        //  for(unsigned i=0; i<b.rows(); ++i)
        //    s2 << b(i) << endl;
        //  s2.close();
    }

#ifdef USECOMISO
    void solveRoundings()
    {
        unsigned n = A.rows();

        gmm::col_matrix< gmm::wsvector< double > > gmm_A;
        std::vector<double> gmm_b;
        std::vector<int> ids_to_round;
        std::vector<double> x;

        gmm_A.resize(n,n);
        gmm_b.resize(n);
        x.resize(n);

        // Copy A
        for (int k=0; k<A.outerSize(); ++k)
            for (SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
            {
                gmm_A(it.row(),it.col()) += it.value();
            }

        // Copy b
        for(unsigned i=0; i<n;++i)
            gmm_b[i] = b[i];

        // Set variables to round
        ids_to_round.clear();
        for(unsigned i=0; i<tag_p.size();++i)
            if(tag_p[i] != -1)
                ids_to_round.push_back(tag_p[i]);

        // Empty constraints
        gmm::row_matrix< gmm::wsvector< double > > gmm_C(0, n);

        COMISO::ConstrainedSolver cs;
        //print_miso_settings(cs.misolver());
        cs.solve(gmm_C, gmm_A, x, gmm_b, ids_to_round, 0.0, false, true);

        // Copy the result back
        for(unsigned i=0; i<F.rows(); ++i)
            if (tag_t[i] != -1)
                angles[i] = x[tag_t[i]];
            else
                angles[i] = hard[i];

        for(unsigned i=0; i<EF.rows(); ++i)
            if(tag_p[i]  != -1)
                p[i] = roundl(x[tag_p[i]]);

    }
#endif

    void solveNoRoundings()
    {
        // Solve the linear system
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
        solver.compute(A);
        Eigen::VectorXd x = solver.solve(b);

        // Copy the result back
        for(unsigned i=0; i<F.rows(); ++i)
            if (tag_t[i] != -1)
                angles[i] = x(tag_t[i]);
            else
                angles[i] = hard[i];

        for(unsigned i=0; i<EF.rows(); ++i)
            if(tag_p[i]  != -1)
                p[i] = roundl(x[tag_p[i]]);
    }

    void roundAndFix()
    {
        for(unsigned i=0; i<p.rows(); ++i)
            pFixed[i] = true;
    }

    void findCones(int N)
    {
        // Compute I0, see http://www.graphics.rwth-aachen.de/media/papers/bommes_zimmer_2009_siggraph_011.pdf for details

        Eigen::VectorXd I0 = Eigen::VectorXd::Zero(V.rows());

        // first the k
        for (unsigned i=0; i < EV.rows(); ++i)
        {
            if (!isBorderEdge[i])
            {
                I0(EV(i,0)) -= k(i);
                I0(EV(i,1)) += k(i);
            }
        }

        // then the A
        Eigen::VectorXd A = angleDefect();

        I0 = I0 + A;

        // normalize
        I0 = I0 / (2*M_PI);

        // round to integer (remove numerical noise)
        for (unsigned i=0; i < I0.size(); ++i)
            I0(i) = round(I0(i));

        // compute I
        Eigen::VectorXd I = I0;

        for (unsigned i=0; i < EV.rows(); ++i)
        {
            if (!isBorderEdge[i])
            {
                I(EV(i,0)) -= double(p(i))/double(N);
                I(EV(i,1)) += double(p(i))/double(N);
            }
        }

        // Clear the vertices on the edges
        for (unsigned i=0; i < EV.rows(); ++i)
        {
            if (isBorderEdge[i])
            {
                I0(EV(i,0)) = 0;
                I0(EV(i,1)) = 0;
                I(EV(i,0)) = 0;
                I(EV(i,1)) = 0;
                A(EV(i,0)) = 0;
                A(EV(i,1)) = 0;
            }
        }

        singularityIndex = I;
    }

    Eigen::VectorXd angleDefect()
    {
      Eigen::VectorXd A = Eigen::VectorXd::Constant(V.rows(),-2*M_PI);

      for (unsigned i=0; i < F.rows(); ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          Eigen::VectorXd a = V.row(F(i,(j+1)%3)) - V.row(F(i,j));
          Eigen::VectorXd b = V.row(F(i,(j+2)%3)) - V.row(F(i,j));
          double t = a.transpose()*b;
          t /= (a.norm() * b.norm());
          A(F(i,j)) += acos(t);
        }
      }

      return A;
    }

public:


    ImplicitSmoother(TriMeshType &_mesh):mesh(_mesh)
    {
        initializeSmoother();
    }

    void setSoftAlpha(double alpha)
    {
        assert(alpha >= 0 && alpha < 1);
        softAlpha = alpha;
    }

    void setConstraintSoft(const int fid, const double w, const CoordType &v)
    {
        // create eigen vector
        Eigen::Vector3d c=vcg::MeshToMatrix<TriMeshType>::VectorFromCoord(v);

        //        // copy coordinates
        //        for (int i = 0; i < 3; i++)
        //            c(i) = v[i];

        //        // set smoother soft constraint
        //        smoother->setConstraintSoft(fid, w, c);

        wSoft(fid) = w;
        soft(fid) = convert3DtoLocal(fid, c);
    }


    void setConstraintHard(const int fid, const CoordType &v)
    {
        Eigen::Vector3d c=vcg::MeshToMatrix<TriMeshType>::VectorFromCoord(v);

        isHard[fid] = true;
        hard(fid) = convert3DtoLocal(fid, c);
    }



    void solve(const int N)
    {
        // Reduce the search space by fixing matchings
        reduceSpace();

        // Build the system
        prepareSystemMatrix(N);

#ifdef USECOMISO
        // Solve with integer roundings
        solveRoundings();
#else
        // Solve with no roundings
        solveNoRoundings();

        // Round all p and fix them
        roundAndFix();

        // Build the system
        prepareSystemMatrix(N);

        // Solve with no roundings (they are all fixed)
        solveNoRoundings();
#endif

        // Find the cones
        findCones(N);
    }


    void getFieldPerFace(vector<CoordType> &fields)
    {
        //assert(smoother);

        // get fields
        //Eigen::MatrixXd fs = smoother->getFieldPerFace();

        Eigen::MatrixXd fs(F.rows(),3);
        for(unsigned i=0; i<F.rows(); ++i)
            fs.row(i) = convertLocalto3D(i, angles(i));

        // reset output vector of fields
        fields.clear();

        // resize output vector of fields
        fields.resize(fs.rows());

        // copy fields in output vector
        for (int i = 0; i < fs.rows(); i++)
            for (int j = 0; j < 3; j++)
                fields[i][j] = fs(i,j);
    }


    void getFFieldPerFace(vector< pair<CoordType,CoordType> > &ffields)
    {
        // get fields
        //Eigen::MatrixXd ffs = smoother->getFFieldPerFace();

        Eigen::MatrixXd ffs(F.rows(),6);
        for(unsigned i=0; i<F.rows(); ++i)
        {
            Eigen::Vector3d v1 = convertLocalto3D(i, angles(i));
            Eigen::Vector3d n = N.row(i);
            Eigen::Vector3d v2 = n.cross(v1);
            v1.normalize();
            v2.normalize();

            ffs.block(i,0,1,3) = v1.transpose();
            ffs.block(i,3,1,3) = v2.transpose();
        }
        //return result;

        // reset output vector of fields
        ffields.clear();

        // resize output vector of fields
        ffields.resize(ffs.rows());

        // copy fields in output vector
        for (int i = 0; i < ffs.rows(); i++)
            for (int j = 0; j < 3; j++)
            {
                ffields[i].first[j] = ffs(i,j);
                ffields[i].second[j] = ffs(i,3+j);
            }
    }

    void getSingularityIndexPerVertex(vector<ScalarType> &sings)
    {
        // get fields
        //Eigen::VectorXd s = smoother->getSingularityIndexPerVertex();

        // reset output vector of singularities
        sings.clear();

        // resize output vector of singularities
        sings.resize(singularityIndex.rows());

        // copy fields in output vector
        for (int i = 0; i < singularityIndex.rows(); i++)
            sings[i] = singularityIndex(i);
    }


    void getSingularitiesIndexPerVertexList(list<int> &sIndexes, const ScalarType t = ScalarType(0))
    {

        // get singularities vector
        vector<ScalarType> sings;
        getSingularityIndexPerVertex(sings);

        // reset output list of indexes
        sIndexes.clear();

        // search for indexes with singularty greater than or equal to t
        for (unsigned long i = 0; i < sings.size(); i++)
            if (sings[i] >= t)
                sIndexes.push_back(i);
    }

};


template <class MeshType>
class FieldSmoother
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType CoordType;


    static void InitQualityByAnisotropyDir(MeshType &mesh)
    {
        std::vector<ScalarType> QVal;
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            ScalarType N1=fabs(mesh.vert[i].K1());
            ScalarType N2=fabs(mesh.vert[i].K2());

            ScalarType NMax=std::max(N1,N2);
            //ScalarType NMin=std::min(N1,N2);

            ScalarType CurvAni=NMax;//fabs((NMax-NMin)/(NMax+NMin));
            QVal.push_back(CurvAni);
            mesh.vert[i].Q()=CurvAni;
        }
        std::sort(QVal.begin(),QVal.end());
        int percUp=int(floor(QVal.size()*0.95+0.5));
        int percDown=int(floor(QVal.size()*0.05+0.5));
        ScalarType trimUp=QVal[percUp];
        ScalarType trimDown=QVal[percDown];
        vcg::tri::UpdateQuality<MeshType>::VertexClamp(mesh,trimDown,trimUp);
        vcg::tri::UpdateQuality<MeshType>::FaceFromVertex(mesh);
        vcg::tri::UpdateColor<MeshType>::PerFaceQualityGray(mesh,trimUp,trimDown);
    }

    static void SetEdgeDirection(FaceType *f,int edge)
    {
        CoordType dir=f->P0(edge)-f->P1(edge);
        dir.Normalize();
        ScalarType prod1=fabs(dir*f->PD1());
        ScalarType prod2=fabs(dir*f->PD2());
        if (prod1>prod2)
        {
            f->PD1()=dir;
            f->PD2()=f->N()^dir;
        }else
        {
            f->PD2()=dir;
            f->PD1()=f->N()^dir;
        }
    }

    static void AddSharpEdgesConstraints(MeshType & mesh,
                                         const ScalarType &thr=0.2)
    {
        for (size_t i=0;i<mesh.face.size();i++)
            for (int j=0;j<mesh.face[i].VN();j++)
            {
                FaceType *f0=&mesh.face[i];
                FaceType *f1=f0->FFp(j);
                if (f0==f1)continue;
                CoordType N0=f0->N();
                CoordType N1=f1->N();
                if ((N0*N1)>thr)continue;
                SetEdgeDirection(f0,j);
                f0->SetS();
            }
    }

    static void AddBorderConstraints(MeshType & mesh)
    {
        for (size_t i=0;i<mesh.face.size();i++)
            for (int j=0;j<mesh.face[i].VN();j++)
            {
                FaceType *f0=&mesh.face[i];
                FaceType *f1=f0->FFp(j);
                assert(f1!=NULL);
                if (f0!=f1)continue;
                SetEdgeDirection(f0,j);
                f0->SetS();
            }
    }

    static void AddCurvatureConstraints(MeshType & mesh,const ScalarType &thr=0.5)
    {
        for (size_t i=0;i<mesh.face.size();i++)
            if (mesh.face[i].Q()>thr)mesh.face[i].SetS();
    }

public:

    struct SmoothParam
    {
        int Ndir;

        ScalarType alpha_soft;

        bool soft_weight;
        bool align_borders;
        bool align_sharp;
        bool hard_curvature;

        ScalarType sharp_thr;
        ScalarType curv_thr;

        SmoothParam()
        {
            Ndir=4;
            alpha_soft=0.0;
            soft_weight=false;

            align_borders=false;
            align_sharp=false;
            hard_curvature=false;

            sharp_thr=0.1;
            curv_thr=0.8;
        }

    };

    static void InitByCurvature(MeshType & mesh)
    {
        vcg::tri::UpdateCurvatureFitting<MeshType>::computeCurvature(mesh);
        vcg::tri::CrossField<MeshType>::SetFaceCrossVectorFromVert(mesh);
        InitQualityByAnisotropyDir(mesh);
    }

    static void GloballyOrient(MeshType &mesh)
    {
        vcg::tri::CrossField<MeshType>::MakeDirectionFaceCoherent(mesh,true);
    }

    static void SmoothDirections(MeshType &mesh,
                                 SmoothParam SParam=SmoothParam())
    {

        //calculathe the maximim weight if needed
        ScalarType MaxQ=-1;
        if (SParam.soft_weight)
        {
            std::pair<ScalarType,ScalarType> MinMax=vcg::tri::Stat<MeshType>::ComputePerFaceQualityMinMax(mesh);
            MaxQ=MinMax.second;
        }

        assert(SParam.alpha_soft>=0);

        ImplicitSmoother<MeshType> SmoothW(mesh);
        //SmoothW.initializeSmoother();

        //if alpha==0 then do not consider initial constraints and set to a value by default
        SmoothW.setSoftAlpha(0.001);

        if (SParam.alpha_soft>0)
            SmoothW.setSoftAlpha(SParam.alpha_soft);

        //set hard constraints if needed
        vcg::tri::UpdateFlags<MeshType>::FaceClearS(mesh);

        if (SParam.align_borders)
            AddBorderConstraints(mesh);
        if (SParam.align_sharp)
            AddSharpEdgesConstraints(mesh,SParam.sharp_thr);
        if (SParam.hard_curvature)
            AddCurvatureConstraints(mesh,SParam.curv_thr);

        //then set hard constraints
        int num_fixed=0;
        for (size_t i=0;i<mesh.face.size();i++)
        {
            CoordType dir=mesh.face[i].PD1();
            dir.Normalize();
            if (mesh.face[i].IsS())
            {
                SmoothW.setConstraintHard(i,dir);
                num_fixed++;
            }
        }

        //check if alpha >0
        if (SParam.alpha_soft>0)
        {
            for (size_t i=0;i<mesh.face.size();i++)
            {

                ScalarType W=1;
                if (SParam.soft_weight)
                    ScalarType W=mesh.face[i].Q()/MaxQ;

                CoordType dir=mesh.face[i].PD1();
                dir.Normalize();
                SmoothW.setConstraintSoft(i,W,dir);
            }
        }

        //fix one face by default if none has been fixed and no soft constraints eather
        if ((num_fixed==0)&&(SParam.alpha_soft==0))
        {
            CoordType dirN=mesh.face[0].N();
            dirN.Normalize();
            CoordType dir1=CoordType(1,0,0);
            if (fabs(dir1*dirN)>0.9)
                dir1=CoordType(0,1,0);
            if (fabs(dir1*dirN)>0.9)
                dir1=CoordType(0,0,1);

            dir1=dirN^dir1;
            dir1.Normalize();

            SmoothW.setConstraintHard(0,dir1);
        }

        //solve
        SmoothW.solve(SParam.Ndir);

        ///assign back to vertices
        std::vector<CoordType> Field;
        SmoothW.getFieldPerFace(Field);
        for (size_t i=0;i<mesh.face.size();i++)
        {
            CoordType dir1=Field[i];
            dir1.Normalize();
            CoordType dir2=mesh.face[i].N()^dir1;
            dir2.Normalize();

            mesh.face[i].PD1()=dir1;
            mesh.face[i].PD2()=dir2;
        }
    }

};

#endif // SMOOTHER_FIELD_H
