#include <bitset>

// #include <wrap/ply/plystuff.h>
#include <wrap/io_trimesh/import.h>

#include <vcg/complex/algorithms/align_pair.h>
#include <vcg/space/index/grid_static_obj.h>

#ifndef VCGLIB_OCCUPANCY_GRID_H
#define VCGLIB_OCCUPANCY_GRID_H

#define OG_MAX_MCB_SIZE 2048
#define OG_MESH_INFO_MAX_STAT 64

namespace vcg {
    template<class MeshType, class ScalarType>
    class OccupancyGrid {

    public:

        /**
         * Class to keep for each voxel the id of the mesh passing through it.
         * based on bitset
         */
        class MeshCounter {

        private:
            std::bitset<OG_MAX_MCB_SIZE> cnt;

        public:

            static constexpr int MaxVal() {
                return OG_MAX_MCB_SIZE;
            }

            bool Empty() const {
                return cnt.none();
            }

            void Clear() {
                cnt.reset();
            }

            bool IsSet(size_t i) const {
                return cnt.test(i);
            }

            void Set(size_t i) {
                cnt.set(i);
            }

            void UnSet(size_t i) {
                cnt.reset(i);
            }

            size_t Count() const {
                return cnt.count();
            }

            /**
             * Return a vector with all the id of the meshes
             */
            void Pack(std::vector<int> &v) const {

                v.clear();

                for (size_t i = 0; i < MeshCounter::MaxVal(); ++i) {
                    if (cnt.test(i)) {
                        v.push_back(i);
                    }
                }
            }

            bool operator < (const MeshCounter &c) const {

                if (cnt == c.cnt) return false;

                std::size_t ii = 0;

                while (ii < MeshCounter::MaxVal()){
                    if (cnt[ii] != c.cnt[ii]) {
                        return cnt[ii] < c.cnt[ii];
                    }
                    ++ii;
                }
                return false;
            }
        };

        /**
         * Class for collecting cumulative information about each mesh in the OG.
         * This info are collected in the Compute() by scanning the OG after we filled it with all the meshes.
         */
        class OGMeshInfo {

        public:

            int id {-1};        // the id of the mesh
            int area {0};      // number of voxels in the OG touched by this mesh
            int coverage {0};  // quanto e' ricoperta da altre mesh eccetto se stessa (eg: se ho due mesh di 1000 con overlap al 30% la covrg e' 300)

            bool used = false;

            std::vector<int> densityDistribution; // Distribution of the of the density of the voxels touched by this mesh:
            // densityDistribution[i] says how many voxel (among the ones coverd by this mesh)
            // are covered by <i> othermeshes. Sum(densityDistribution) == area;
            // if densityDistribution[1] > 0 means that this mesh is the unique to cover some portion of the space.

            void Init(int _id) {
                id=_id;
            }

            bool operator < (OGMeshInfo &o) const {
                return area < o.area;
            }

            static constexpr int MaxStat() {
                return OG_MESH_INFO_MAX_STAT;
            }
        };

        /**
         * Classe con informazioni su un arco plausibile
         */
        class OGArcInfo {
        public:

            enum sort {
                AREA,
                NORM_AREA,
                DEGREE
            };

            int s, t; // source and target (come indici nel gruppo corrente)
            int area;  //
            float norm_area;

            OGArcInfo(int _s, int _t, int _area, float _n) : s{_s}, t{_t}, area{_area}, norm_area{_n} {}
            OGArcInfo(int _s, int _t, int _area) : s{_s}, t{_t}, area{_area} {}

            bool operator <  (const OGArcInfo &p) const {
                return norm_area <  p.norm_area;
            }
        };

        GridStaticObj<MeshCounter, float> G;

        int mn;
        int TotalArea;
        /**
         * Maximum number of meshes that cross a cell
         */
        int MaxCount;

        /**
         * SortedVisual Arcs
         */
        std::vector<OGArcInfo>  SVA;  // SortedVirtual Arcs;
        /**
         * High level information for each mesh. Mapped by mesh id
         */
        std::map<int, OGMeshInfo> VM;

        bool Init(int _mn, Box3<ScalarType> bb, int size) {

            // the number of meshes (including all the unused ones; eg it is the range of the possible id)
            mn = _mn;
            if (mn > MeshCounter::MaxVal()) return false;

            MeshCounter MC;

            MC.Clear();
            G.Create(bb,size,MC);
            VM.clear();

            return true;
        }

        void Add(const char *MeshName, Matrix44<ScalarType> &Tr, int id) {

            AlignPair::A2Mesh M;

            vcg::tri::io::Importer<AlignPair::A2Mesh>::Open(M, MeshName);
            vcg::tri::Clean<AlignPair::A2Mesh>::RemoveUnreferencedVertex(M);

            AddMesh(M,Tr,id);
        }

        void AddMeshes(std::vector<std::string> &names, std::vector<Matrix44<ScalarType>> &trv,int size) {

            Box3<ScalarType> bb, totalbb;

            bb.SetNull();
            totalbb.SetNull();

            std::fprintf(stdout, "OG::AddMesh:Scanning BBoxes\n");

            for (std::size_t i = 0; i < names.size(); ++i) {
                // vcg::ply::ScanBBox(names[i].c_str(), bb, true);
                totalbb.Add(trv[i], bb);
            }

            Init(names.size(),totalbb,size);

            for (std::size_t i = 0; i < names.size(); ++i) {
                std::fprintf(stdout, "OG::AddMesh:Adding Mesh %i '%s'\n", i, names[i].c_str());
                Add(names[i].c_str(), trv[i], i);
            }
        }

        void AddMesh(MeshType &mesh, const Matrix44<ScalarType> &Tr, int ind) {

            Matrix44f Trf;
            Trf.Import(Tr);

            for (auto vi = std::begin(mesh.vert); vi != std::end(mesh.vert); ++vi) {

                if (!(*vi).IsD()) {
                    G.Grid(Trf * Point3f::Construct((*vi).P())).Set(ind);
                }
            }

            VM[ind].Init(ind);
            VM[ind].used = true;
        }

        void RemoveMesh(int id) {

            MeshCounter *GridEnd = G.grid + G.size();

            for (MeshCounter* ig = G.grid; ig != GridEnd; ++ig) {
                ig->UnSet(id);
            }
        }

        /**
         * This function is called after we have <added> all the mesh to the OG
         * to collect the information about the interferences between the various meshes.
         */
        void Compute() {

            // Analisi della griglia
            // Si deve trovare l'insieme degli archi piu'plausibili
            // un arco ha "senso" in una cella se entrambe le mesh compaiono in quell'arco
            // Si considera tutti gli archi possibili e si conta in quante celle ha senso un arco

            std::vector<int> VA; // virtual arcs
            VA.resize(mn * mn, 0);

            std::map<std::pair<int, int>, int> VAMap;

            // First Loop:
            // Scan the grid and update possible arc count
            for (int i = 0; i < G.siz[0]; ++i) {

                for (int j = 0; j < G.siz[1]; ++j) {

                    for (int k = 0; k < G.siz[2]; ++k) {

                        std::vector<int> vv;
                        G.Grid(i, j, k).Pack(vv);
                        std::size_t meshInCell = vv.size();

                        for (std::size_t ii = 0; ii < vv.size(); ++ii) {

                            OccupancyGrid::OGMeshInfo &omi_ii = VM[vv[ii]];

                            ++omi_ii.area; // compute mesh area
                            if (meshInCell > omi_ii.densityDistribution.size()) {
                                omi_ii.densityDistribution.resize(meshInCell);
                            }

                            ++(omi_ii.densityDistribution[meshInCell - 1]);
                        }

                        for (std::size_t ii = 0; ii < vv.size(); ++ii) {
                            for (std::size_t jj = ii + 1; jj < vv.size(); ++jj) {
                                // count intersections of all mesh pairs
                                ++VAMap[std::make_pair(vv[ii], vv[jj])];
                            }
                        }
                    }
                }
            }

            // Find all the arcs, e.g. all the pair of meshes
            SVA.clear();
            for (auto vi = std::begin(VAMap); vi != std::end(VAMap); ++vi) {
                if (vi->second > 0) {
                    int m_s = vi->first.first;
                    int m_t = vi->first.second;
                    int area = vi->second;
                    SVA.push_back( OGArcInfo (m_s,m_t,area,float(area)/float(std::min(VM[m_s].area,VM[m_t].area)) ));
                }
            }

            // Compute Mesh Coverage
            for (std::size_t i = 0; i < SVA.size(); ++i) {
                VM[SVA[i].s].coverage += SVA[i].area;
                VM[SVA[i].t].coverage += SVA[i].area;
            }

            std::sort(std::begin(SVA), std::end(SVA));
            std::reverse(std::begin(SVA), std::end(SVA));
        }

        void ComputeUsefulMesh(FILE *elfp) {

            int mcnt = 0;
            std::vector<int> UpdArea(mn);
            std::vector<int> UpdCovg(mn);

            for (std::size_t m = 0; m < mn; ++m) {

                if (VM[m].used && VM[m].area > 0) {
                    mcnt++;
                    UpdCovg[m]=VM[m].coverage;
                    UpdArea[m]=VM[m].area;
                }
            }

            int sz = G.size();
            if (elfp) {
                std::fprintf(elfp, "\n\nComputing Usefulness of Meshes of %i(on %i) meshes\n Og with %i / %i fill ratio %i max mesh per cell\n\n", mcnt, mn, TotalArea, sz, MaxCount);
                std::fprintf(elfp, "\n");
            }

            int CumArea = 0;

            for (std::size_t m = 0; m < mn-1; ++m) {

                int best = max_element(std::begin(UpdArea), std::end(UpdArea)) - std::begin(UpdArea);

                CumArea += UpdArea[best];
                if (UpdCovg[best] < 0) break;

                // se era una mesh fuori del working group si salta tutto.
                if (VM[best].area == 0) continue;

                if (elfp) {
                    fprintf(elfp, "%3i %3i %7i (%7i) %7i %5.2f %7i(%7i)\n",
                            m, best, UpdArea[best], VM[best].area, TotalArea - CumArea,
                            100.0 - 100 * float(CumArea) / TotalArea, UpdCovg[best], VM[best].coverage);
                }

                UpdArea[best] = -1;
                UpdCovg[best] = -1;

                for (std::size_t i = 0; i < sz; ++i) {

                    MeshCounter &mc = G.grid[i];

                    if (mc.IsSet(best))	{

                        mc.UnSet(best);

                        for (std::size_t j = 0; j < mn; ++j) {
                            if (mc.IsSet(j)) {
                                --UpdArea[j];
                                UpdCovg[j]-=mc.Count();
                            }
                        }

                        mc.Clear();
                    }
                }
            }
        }

        void Dump(FILE *fp) {

            std::fprintf(fp, "Occupancy Grid\n");
            std::fprintf(fp, "grid of ~%i kcells: %d x %d x %d\n", G.size(), G.siz[0], G.siz[1], G.siz[2]);
            std::fprintf(fp, "grid voxel size of %f %f %f\n", G.voxel[0], G.voxel[1], G.voxel[2]);

            std::fprintf(fp,"Computed %lu arcs for %i meshes\n", SVA.size(), mn);

            for (std::size_t i=0;i<VM.size();++i) {

                if (VM[i].used) {

                    std::fprintf(fp, "mesh %3lu area %6i covg %7i (%5.2f%%) DensDistr:", i, VM[i].area, VM[i].coverage, float(VM[i].coverage)/float(VM[i].area));

                    for (std::size_t j = 0; j < std::min(static_cast<std::size_t>(8), VM[i].densityDistribution.size()); ++j) {
                        std::fprintf(fp," %3i ", VM[i].densityDistribution[j]);
                    }

                    std::fprintf(fp,"\n");
                }
                else {
                    std::fprintf(fp, "mesh %3lu ---- UNUSED\n", i);
                }
            }

            std::fprintf(fp, "Computed %lu Arcs :\n", SVA.size());

            for (std::size_t i = 0; i < SVA.size() && SVA[i].norm_area > .1; ++i) {
                std::fprintf(fp, "%4i -> %4i Area:%5i NormArea:%5.3f\n", SVA[i].s, SVA[i].t, SVA[i].area, SVA[i].norm_area);
            }

            std::fprintf(fp, "End OG Dump\n");
        }

        void ComputeTotalArea() {

            using uint = unsigned int;

            int ccnt = 0;
            MaxCount = 0;

            int sz = G.size();

            for (int i = 0; i < sz; ++i)

                if (!G.grid[i].Empty()) {

                    ccnt++;
                    if (G.grid[i].Count() > static_cast<uint>(MaxCount)) {
                        MaxCount = G.grid[i].Count();
                    }
                }

            TotalArea = ccnt;
        }

    };
}

#endif //VCGLIB_OCCUPANCY_GRID_H
