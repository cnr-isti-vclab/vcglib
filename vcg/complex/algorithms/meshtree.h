#ifndef VCGLIB_MESHTREE_H
#define VCGLIB_MESHTREE_H

#include <vcg/complex/algorithms/align_pair.h>
#include <vcg/complex/algorithms/align_global.h>
#include <vcg/complex/algorithms/occupancy_grid.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace vcg {

    template<class MeshType, class ScalarType>
    class MeshTree {

    public:

        class MeshNode {

        public:
            bool glued;
            MeshType *m;

            explicit MeshNode(MeshType *_m) : m{_m}, glued{false} {}

            vcg::Matrix44<ScalarType> &tr() {
                return m->cm.Tr;
            }

            const vcg::Box3<ScalarType> &bbox() const {
                return m->cm.bbox;
            }

            int Id() {
                return m->id();
            }
        };

        class Param {
        public:
            int OGSize = 50000;
            float arcThreshold = 0.3f;
            float recalcThreshold = 0.1f;
        };

        std::map<int, MeshNode*> nodeMap;
        std::vector<vcg::AlignPair::Result> resultList;

        vcg::OccupancyGrid<CMeshO, ScalarType> OG{};
        vcg::CallBackPos* cb = vcg::DummyCallBackPos;

        MeshTree() = default;

        ~MeshTree() { clear(); }

        MeshType *MM(unsigned int i) {
            return nodeMap[i]->m;
        }

        void clear() {

            for (auto& ni : nodeMap) {
                delete ni.second;
            }

            nodeMap.clear();
            resultList.clear();
        }

        void deleteResult(MeshTree::MeshNode *mp) {

            auto li = std::begin(resultList);
            while (li != resultList.end()) {

                if (li->MovName == mp->Id() || li->FixName == mp->Id()) {
                    li = resultList.erase(li);
                }
                else {
                    ++li;
                }
            }
        }

        vcg::AlignPair::Result* findResult(int id1, int id2) {

            for (auto& li : resultList) {
                if ((li.MovName == id1 && li.FixName == id2) || (li.MovName == id2 && li.FixName == id1) ) {
                    return &li;
                }
            }

            return nullptr;
        }

        MeshTree::MeshNode *find(int id) {

            MeshTree::MeshNode *mp = nodeMap[id];

            if (mp == nullptr || mp->Id() != id) {
                assert("You are trying to find a non existent mesh" == nullptr);
            }

            return mp;
        }

        MeshTree::MeshNode *find(MeshType *m) {

            for (auto& ni : nodeMap) {
                if (ni.second->m == m) return ni.second;
            }

            assert("You are trying to find a non existent mesh" == nullptr);
            return nullptr;
        }

        int gluedNum() {

            int count = 0;

            for (auto& ni : nodeMap) {
                if (ni.second->glued) ++count;
            }

            return count;
        }

        void Process(vcg::AlignPair::Param &ap, MeshTree::Param &mtp) {

            char buf[1024];
            std::sprintf(buf, "Starting Processing of %i glued meshes out of %zu meshes\n", gluedNum(), nodeMap.size());
            cb(0, buf);

            /******* Occupancy Grid Computation *************/
            std::memset(buf, '\0', 1024);
            std::sprintf(buf, "Computing Overlaps %i glued meshes...\n", gluedNum());
            cb(0, buf);

            OG.Init(static_cast<int>(nodeMap.size()), vcg::Box3<ScalarType>::Construct(gluedBBox()), mtp.OGSize);

            for (auto& ni : nodeMap) {
                MeshTree::MeshNode *mn = ni.second;
                if (mn->glued) {
                    OG.AddMesh(mn->m->cm, vcg::Matrix44<ScalarType>::Construct(mn->tr()), mn->Id());
                }
            }

            OG.Compute();
            OG.Dump(stdout);
            // Note: the s and t of the OG translate into fix and mov, respectively.

            /*************** The long loop of arc computing **************/

            // count existing arcs within current error threshold
            float percentileThr = 0;
            if (!resultList.empty()) {

                vcg::Distribution<float> H;
                for (auto& li : resultList) {
                    H.Add(li.err);
                }

                percentileThr = H.Percentile(1.0f - mtp.recalcThreshold);
            }

            std::size_t totalArcNum = 0;
            int preservedArcNum = 0, recalcArcNum = 0;

            while(totalArcNum < OG.SVA.size() && OG.SVA[totalArcNum].norm_area > mtp.arcThreshold)
            {
                AlignPair::Result *curResult = findResult(OG.SVA[totalArcNum].s, OG.SVA[totalArcNum].t);
                if (curResult) {
                    if (curResult->err < percentileThr) {
                        ++preservedArcNum;
                    }
                    else {
                        ++recalcArcNum;
                    }
                }
                else {
                    resultList.push_back(AlignPair::Result());
                    resultList.back().FixName = OG.SVA[totalArcNum].s;
                    resultList.back().MovName = OG.SVA[totalArcNum].t;
                    resultList.back().err = std::numeric_limits<double>::max();
                }
                ++totalArcNum;
            }

            //if there are no arcs at all complain and return
            if (totalArcNum == 0) {
                std::memset(buf, '\0', 1024);
                std::sprintf(buf, "\n Failure. There are no overlapping meshes?\n No candidate alignment arcs. Nothing Done.\n");
                cb(0, buf);
                return;
            }

            int num_max_thread = 1;
#ifdef _OPENMP
            if (totalArcNum > 32) num_max_thread = omp_get_max_threads();
#endif
            std::memset(buf, '\0', 1024);
            std::sprintf(buf, "Arc with good overlap %6zu (on  %6zu)\n", totalArcNum, OG.SVA.size());
            cb(0, buf);

            std::memset(buf, '\0', 1024);
            std::sprintf(buf, " %6i preserved %i Recalc \n", preservedArcNum, recalcArcNum);
            cb(0, buf);

            bool hasValidAlign = false;

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_max_thread)

            // on windows, omp does not support unsigned types for indices on cycles
            for (int i = 0 ;i < static_cast<int>(totalArcNum); ++i)  {

                std::fprintf(stdout,"%4i -> %4i Area:%5i NormArea:%5.3f\n",OG.SVA[i].s,OG.SVA[i].t,OG.SVA[i].area,OG.SVA[i].norm_area);
                AlignPair::Result *curResult = findResult(OG.SVA[i].s,OG.SVA[i].t);

                // // missing arc and arc with great error must be recomputed.
                if (curResult->err >= percentileThr) {

                    ProcessArc(OG.SVA[i].s, OG.SVA[i].t, *curResult, ap);
                    curResult->area = OG.SVA[i].norm_area;

                    if (curResult->isValid()) {
                        hasValidAlign = true;
                        std::pair<double, double> dd = curResult->computeAvgErr();
#pragma omp critical

                        std::memset(buf, '\0', 1024);
                        std::sprintf(buf, "(%3i/%3zu) %2i -> %2i Aligned AvgErr dd=%f -> dd=%f \n", i+1,totalArcNum,OG.SVA[i].s,OG.SVA[i].t,dd.first,dd.second);
                        cb(0, buf);
                    }
                    else {
#pragma omp critical
                        std::memset(buf, '\0', 1024);
                        std::sprintf(buf, "(%3i/%3zu) %2i -> %2i Failed Alignment of one arc %s\n",i+1,totalArcNum,OG.SVA[i].s,OG.SVA[i].t,vcg::AlignPair::errorMsg(curResult->status));
                        cb(0, buf);
                    }
                }
            }

            //if there are no valid arcs complain and return
            if (!hasValidAlign) {
                std::memset(buf, '\0', 1024);
                std::sprintf(buf, "\n Failure. No successful arc among candidate Alignment arcs. Nothing Done.\n");
                cb(0, buf);
                return;
            }

            vcg::Distribution<float> H; // stat for printing
            for (auto& li : resultList) {
                if (li.isValid()) H.Add(li.err);
            }

            std::memset(buf, '\0', 1024);
            std::sprintf(buf, "Completed Mesh-Mesh Alignment: Avg Err %5.3f; Median %5.3f; 90%% %5.3f\n", H.Avg(), H.Percentile(0.5f), H.Percentile(0.9f));
            cb(0, buf);

            ProcessGlobal(ap);
        }

        void ProcessGlobal(vcg::AlignPair::Param &ap) {

            char buff[1024];
            std::memset(buff, '\0', 1024);

            /************** Preparing Matrices for global alignment *************/
            std::vector<int> GluedIdVec;
            std::vector<vcg::Matrix44d> GluedTrVec;

            std::map<int, std::string> names;

            for (auto& ni : nodeMap) {

                MeshTree::MeshNode *mn = ni.second;
                if (mn->glued) {
                    GluedIdVec.push_back(mn->Id());
                    GluedTrVec.push_back(vcg::Matrix44d::Construct(mn->tr()));
                    names[mn->Id()] = qUtf8Printable(mn->m->label());
                }
            }

            vcg::AlignGlobal AG;
            std::vector<vcg::AlignPair::Result *> ResVecPtr;
            for (auto& li : resultList) {
                if (li.isValid()) {
                    ResVecPtr.push_back(&li);
                }
            }

            AG.BuildGraph(ResVecPtr, GluedTrVec, GluedIdVec);

            float StartGlobErr = 0.001f;
            while (!AG.GlobalAlign(names, StartGlobErr, 100, ap.MatchMode == vcg::AlignPair::Param::MMRigid, stdout, cb)) {
                StartGlobErr *= 2;
                AG.BuildGraph(ResVecPtr,GluedTrVec, GluedIdVec);
            }

            std::vector<vcg::Matrix44d> GluedTrVecOut(GluedTrVec.size());
            AG.GetMatrixVector(GluedTrVecOut,GluedIdVec);

            // Now get back the results!
            for (std::size_t ii = 0; ii < GluedTrVecOut.size(); ++ii) {
                MM(GluedIdVec[ii])->cm.Tr.Import(GluedTrVecOut[ii]);
            }

            std::sprintf(buff, "Completed Global Alignment (error bound %6.4f)\n", StartGlobErr);
            cb(0, buff);
        }

        void ProcessArc(int fixId, int movId, vcg::AlignPair::Result &result, vcg::AlignPair::Param ap) {

            // l'allineatore globale cambia le varie matrici di posizione di base delle mesh
            // per questo motivo si aspetta i punti nel sistema di riferimento locale della mesh fix
            // Si fanno tutti i conti rispetto al sistema di riferimento locale della mesh fix
            vcg::Matrix44d FixM = vcg::Matrix44d::Construct(find(fixId)->tr());
            vcg::Matrix44d MovM = vcg::Matrix44d::Construct(find(movId)->tr());
            vcg::Matrix44d MovToFix = Inverse(FixM) * MovM;

            ProcessArc(fixId, movId, MovToFix, result, ap);
        }

        void ProcessArc(int fixId, int movId, vcg::Matrix44d &MovM, vcg::AlignPair::Result &result, vcg::AlignPair::Param ap) {

            vcg::AlignPair::A2Mesh Fix;
            vcg::AlignPair aa;

            // 1) Convert fixed mesh and put it into the grid.
            MM(fixId)->updateDataMask(MeshType::MeshModel::MM_FACEMARK);
            aa.convertMesh<CMeshO>(MM(fixId)->cm,Fix);

            vcg::AlignPair::A2Grid UG;
            vcg::AlignPair::A2GridVert VG;

            if (MM(fixId)->cm.fn==0 || ap.UseVertexOnly) {
                Fix.initVert(vcg::Matrix44d::Identity());
                vcg::AlignPair::InitFixVert(&Fix,ap,VG);
            }
            else {
                Fix.init(vcg::Matrix44d::Identity());
                vcg::AlignPair::initFix(&Fix, ap, UG);
            }

            // 2) Convert the second mesh and sample a <ap.SampleNum> points on it.
            MM(movId)->updateDataMask(MeshType::MeshModel::MM_FACEMARK);
            std::vector<vcg::AlignPair::A2Vertex> tmpmv;
            aa.convertVertex(MM(movId)->cm.vert,tmpmv);
            aa.sampleMovVert(tmpmv, ap.SampleNum, ap.SampleMode);

            aa.mov=&tmpmv;
            aa.fix=&Fix;
            aa.ap = ap;

            // Perform the ICP algorithm
            aa.align(MovM,UG,VG,result);

            result.FixName=fixId;
            result.MovName=movId;
        }

        inline vcg::Box3<ScalarType> bbox() {

            vcg::Box3<ScalarType> FullBBox;
            for (auto& ni : nodeMap) {
                FullBBox.Add(vcg::Matrix44d::Construct(ni.second->tr()), ni.second->bbox());
            }

            return FullBBox;
        }

        inline vcg::Box3<ScalarType> gluedBBox() {

            vcg::Box3<ScalarType> FullBBox;

            for (auto& ni : nodeMap) {
                if (ni.second->glued) {
                    FullBBox.Add(vcg::Matrix44<ScalarType>::Construct(ni.second->tr()), ni.second->bbox());
                }
            }

            return FullBBox;
        }

    };

}

#endif //VCGLIB_MESHTREE_H
