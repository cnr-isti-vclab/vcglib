#include <list>
#include <queue>

#include <wrap/callback.h>

#ifndef MESHLAB_ALIGNGLOBAL_H
#define MESHLAB_ALIGNGLOBAL_H

namespace vcg {
    class AlignGlobal {
    public:

        /**
         * Forward declaration for the `VirtAlign` class.
         */
        class Node;

        /**
         * Allineamento virtuale tra due mesh (estratto da un alignresult).
         * Nota Importante: la trasformazione e i punti qui memorizzati si intendono al netto delle trasf di base delle due mesh in gioco.
         * Quindi se qualcuno sposta una mesh le pos dei punti sono ancora valide ma non la trasf da applicarvi.
         */
        class VirtAlign
        {
        public:

            AlignGlobal::Node *Fix, *Mov; // allineamento tra i e j
            std::vector<vcg::Point3d> FixP; // punti su Fix
            std::vector<vcg::Point3d> MovP; // punti su Mov
            std::vector<vcg::Point3d> FixN; // Normali su Fix
            std::vector<vcg::Point3d> MovN; // Normali su Mov
            vcg::Matrix44d M2F; //la matrice da applicare ai punti di Mov per ottenere quelli su Fix
            vcg::Matrix44d F2M; //la matrice da applicare ai punti di Fix per ottenere quelli su Mov
            /*
                Nel caso semplificato che le mesh avessero come trasf di base l'identita' deve valere:

                         N2A(N).Apply(   P(N)) ~= AdjP(N)
                         A2N(N).Apply(AdjP(N)) ~=    P(N)

                In generale un nodo N qualsiasi dell'VirtAlign vale che:

                    N2A(N).Apply(       N->M.Apply(   P(N)) ) ~= AdjN(N)->M.Apply( AdjP(N) );
                    A2M(N).Apply( AdjN(N)->M.Apply(AdjP(N)) ) ~=       N->M.Apply(    P(N) );

                in cui il ~= significa uguale al netto dell'errore di allineamento.

                Per ottenere i virtualmate relativi ad un nodo n:
            */

            inline vcg::Matrix44d       &N2A(AlignGlobal::Node *n) {if(n==Fix) return F2M; else return M2F;}
            inline vcg::Matrix44d       &A2N(AlignGlobal::Node *n) {if(n==Fix) return M2F; else return F2M;}

            inline std::vector<vcg::Point3d>    &P(AlignGlobal::Node *n) {if(n==Fix) return FixP; else return MovP;}
            inline std::vector<vcg::Point3d>    &N(AlignGlobal::Node *n) {if(n==Fix) return FixN; else return MovN;}

            inline std::vector<vcg::Point3d> &AdjP(AlignGlobal::Node *n) {if(n==Fix) return MovP; else return FixP;}
            inline std::vector<vcg::Point3d> &AdjN(AlignGlobal::Node *n) {if(n==Fix) return MovN; else return FixN;}

            AlignGlobal::Node *Adj(Node *n) const {

                assert(n == Fix || n == Mov);
                if (n == Fix) return Mov;

                return Fix;
            }

            bool Check() const {

                if (FixP.size() != MovP.size()) return false;

                Point3d mp, fp;

                double md = 0, fd = 0;
                double md2 = 0, fd2 = 0;

                Matrix44d &MovTr=Mov->M;
                Matrix44d &FixTr=Fix->M;

                for (std::size_t i = 0; i < FixP.size(); ++i) {

                    mp = MovTr * MovP[i];
                    fp = FixTr * FixP[i];

                    md  +=        Distance(fp, M2F * mp);
                    md2 += SquaredDistance(fp, M2F * mp);

                    fd  +=        Distance(mp, F2M * fp);
                    fd2 += SquaredDistance(mp, F2M * fp);
                }

                int nn = static_cast<int>(MovP.size());

                std::fprintf(stdout, "Arc %3i -> %3i : %i pt\n", Fix->id, Mov->id, nn);
                std::fprintf(stdout, "SquaredSum Distance %7.3f %7.3f Avg %7.3f %7.3f\n", fd2, md2, fd2/nn, md2/nn);
                std::fprintf(stdout, "       Sum Distance %7.3f %7.3f Avg %7.3f %7.3f\n", fd , md ,  fd/nn, md/nn);
                return true;
            }
        };

        class Node {
        public:

            int id;  // id della mesh a cui corrisponde il nodo
            int sid; // Subgraph id;
            Matrix44d M; // La matrice che mette la mesh nella sua posizione di base;
            std::list<VirtAlign*> Adj;

            /***
             * True if the node is inside the active set
             */
            bool Active;

            /***
             * False if it's dormant
             */
            bool Queued;
            bool Discarded;

            Node() : id{-1}, Active{false}, Discarded{false}, Queued{false} {}

            // Allinea un nodo con tutti i suoi vicini
            double AlignWithActiveAdj(bool Rigid) {

                std::printf("--- AlignWithActiveAdj --- \nMoving node %i with respect to nodes:", id);

                for (auto li = std::begin(Adj); li != std::end(Adj); ++li) {
                    if ((*li)->Adj(this)->Active) {
                        std::printf(" %i,", (*li)->Adj(this)->id);
                    }
                }

                std::printf("\n");

                //printf("Base Matrix of Node %i\n",id);print(M);

                // Step 1; Costruiamo le due liste di punti da allineare
                std::vector<Point3d> FixP, MovP, FixN, MovN;
                Box3d FixBox, MovBox;
                FixBox.SetNull(); MovBox.SetNull();

                for (auto li = std::begin(Adj); li != std::end(Adj); ++li) {

                    // scorro tutti i nodi adiacenti attivi
                    if ((*li)->Adj(this)->Active) {
                        //printf("Base Matrix of Node %i adj to %i\n",id,(*li)->Adj(this)->id);
                        //print((*li)->Adj(this)->M);
                        std::vector<Point3d> &AP=(*li)->AdjP(this);   // Punti sul nodo adiacente corrente;
                        std::vector<Point3d> &AN=(*li)->AdjN(this);   // Normali sul nodo adiacente corrente;

                        //printf("Transf that bring points of %i onto %i\n",id,(*li)->Adj(this)->id);
                        //print((*li)->A2N(this));
                        //printf("Transf that bring points of %i onto %i\n",(*li)->Adj(this)->id,id);
                        //print((*li)->N2A(this));
                        vcg::Point3d pf, nf;
                        vcg::Point3d pm;

                        for (std::size_t i = 0; i < AP.size(); ++i) {

                            pf = (*li)->Adj(this)->M*AP[i]; // i punti fissi sono quelli sulla sup degli adiacenti messi nella loro pos corrente
                            FixP.push_back(pf);
                            FixBox.Add(pf);
                            nf=(*li)->Adj(this)->M*Point3d(AP[i]+AN[i])-pf;
                            nf.Normalize();
                            FixN.push_back(nf);

                            pm = (*li)->A2N(this)*pf;
                            MovP.push_back(pm); // i punti che si muovono sono quelli sul adj trasformati in modo tale da cascare sul nodo corr.
                            MovBox.Add(pm);
                        }
                    }
                }

                vcg::Matrix44d out;
                //if(Rigid) ret=ComputeRigidMatchMatrix(out,FixP,MovP);
                //else ret=ComputeMatchMatrix2(out,FixP,FixN,MovP);
                if (Rigid) {
                    ComputeRigidMatchMatrix(FixP,MovP,out);
                }
                else {
                    vcg::PointMatchingScale::computeRotoTranslationScalingMatchMatrix(out, FixP, MovP);
                }

                Matrix44d outIn=vcg::Inverse(out);
                //double maxdiff = MatrixNorm(out);
                double maxdiff = MatrixBoxNorm(out,FixBox);

                //	printf("Computed Transformation:\n"); print(out);printf("--\n");
                //	printf("Collected %i points . Err = %f\n",FixP.size(),maxdiff);

                // La matrice out calcolata e' quella che applicata ai punti MovP li porta su FixP, quindi i punti della mesh corrente
                // La nuova posizione di base della mesh diventa quindi
                // M * out
                // infatti se considero un punto della mesh originale applicarci la nuova matricie significa fare
                // p * M * out

                // M=M*out; //--Orig
                M = out * M;

                // come ultimo step occorre applicare la matrice trovata a tutti gli allineamenti in gioco.
                // scorro tutti i nodi adiacenti attivi
                for (auto li = std::begin(Adj); li != std::end(Adj); ++li) {
                    //print((*li)->N2A(this));
                    //print((*li)->A2N(this));printf("--\n");

                    (*li)->N2A(this)=(*li)->N2A(this)*outIn;
                    (*li)->A2N(this)=(*li)->A2N(this)*out  ;
                    //print((*li)->N2A(this));
                    //print((*li)->A2N(this));printf("--\n");
                }

                return maxdiff;
            }

            double MatrixNorm(vcg::Matrix44d &NewM) const {

                double maxDiff = 0;

                vcg::Matrix44d diff;
                diff.SetIdentity();
                diff = diff-NewM;

                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        maxDiff += (diff[i][j] * diff[i][j]);
                    }
                }

                return maxDiff;
            }

            double MatrixBoxNorm(vcg::Matrix44d &NewM, vcg::Box3d &bb) const {

                double maxDiff = 0;
                vcg::Point3d pt;

                pt = Point3d(bb.min[0], bb.min[1], bb.min[2]); maxDiff = std::max(maxDiff, Distance(pt, NewM * pt));
                pt = Point3d(bb.max[0], bb.min[1], bb.min[2]); maxDiff = std::max(maxDiff, Distance(pt, NewM * pt));
                pt = Point3d(bb.min[0], bb.max[1], bb.min[2]); maxDiff = std::max(maxDiff, Distance(pt, NewM * pt));
                pt = Point3d(bb.max[0], bb.max[1], bb.min[2]); maxDiff = std::max(maxDiff, Distance(pt, NewM * pt));
                pt = Point3d(bb.min[0], bb.min[1], bb.max[2]); maxDiff = std::max(maxDiff, Distance(pt, NewM * pt));
                pt = Point3d(bb.max[0], bb.min[1], bb.max[2]); maxDiff = std::max(maxDiff, Distance(pt, NewM * pt));
                pt = Point3d(bb.min[0], bb.max[1], bb.max[2]); maxDiff = std::max(maxDiff, Distance(pt, NewM * pt));
                pt = Point3d(bb.max[0], bb.max[1], bb.max[2]); maxDiff = std::max(maxDiff, Distance(pt, NewM * pt));

                return maxDiff;
            }

            int PushBackActiveAdj(std::queue<Node *> &Q) {

                assert(Active);

                int count = 0;
                AlignGlobal::Node *pt;

                for (auto li = std::begin(Adj); li != std::end(Adj); ++li) {

                    pt = (*li)->Adj(this);

                    if (pt->Active && !pt->Queued) {
                        ++count;
                        pt->Queued=true;
                        Q.push(pt);
                    }
                }
                return count;
            }

            int DormantAdjNum() {

                int count = 0;

                for (auto li = std::begin(Adj); li != std::end(Adj); ++li) {
                    if (!(*li)->Adj(this)->Active) ++count;
                }

                return count;
            }

            int ActiveAdjNum() {

                int count = 0;

                for (auto li = std::begin(Adj); li != std::end(Adj); ++li) {
                    if ((*li)->Adj(this)->Active) ++count;
                }

                return count;
            }
        };

        class SubGraphInfo {
        public:
            int sid;
            int size;
            Node *root;
        };

        std::list<Node> N;
        std::list<VirtAlign*> A;

        /**
         * Descrittori delle componenti connesse, riempito dalla ComputeConnectedComponents
         */
        std::list<SubGraphInfo> CC;

        static inline void LOG( FILE *fp, const char * f, ... ) {

            if (fp == 0) return;

            va_list marker;
            va_start(marker, f);
            std::vfprintf(fp, f, marker);
            va_end(marker);
            std::fflush(fp);
        }

        inline int DormantNum() const {

            int count = 0;

            for (auto li = std::begin(N); li != std::end(N); ++li) {
                if (!(*li).Active) ++count;
            }

            return count;
        }

        inline int ActiveNum() const {

            int count = 0;

            for (auto li = std::begin(N); li != std::end(N); ++li) {
                if ((*li).Active) ++count;
            }

            return count;
        }

        bool CheckGraph() {

            std::vector<bool> Visited(N.size(), false);
            std::stack<AlignGlobal::Node*> st;

            st.push(&(*N.begin()));
            while (!st.empty()) {

                AlignGlobal::Node *cur=st.top();
                st.pop();
                // std::printf("Visiting node %i\n",cur->id);

                for (auto li = std::begin(cur->Adj); li != std::end(cur->Adj); ++li) {
                    if (!Visited[(*li)->Adj(cur)->id]) {
                        Visited[(*li)->Adj(cur)->id] = true;
                        st.push((*li)->Adj(cur));
                    }
                }
            }

            size_t cnt = std::count(std::begin(Visited), std::end(Visited), true);
            std::printf("Nodes that can be reached from root %zu on %zu \n", cnt, N.size());

            return (cnt == N.size());
        }

        void Dump(FILE *fp) const {
            std::fprintf(fp, "Alignment Graph of %lu nodes and %lu arcs\n", N.size(), A.size());
            //  list<VirtAlign *>::iterator li;
            //	for(li=A.begin();li!=A.end();++li)
            //		printf("Arc : %3i ->%3i\n",(*li)->Fix->id,(*li)->Mov->id);
        }

        int ComputeConnectedComponents() {

            std::printf("Building Connected Components on a graph with %lu nodes and %lu arcs\n", N.size(), A.size());

            CC.clear();

            std::stack<AlignGlobal::Node *> ToReach; // nodi ancora da visitare
            std::stack<AlignGlobal::Node *> st;      // nodi che si stanno visitando

            for (auto li = std::begin(N); li != std::end(N); ++li) {
                (*li).sid = -1;
                ToReach.push(&*li);
            }

            int cnt = 0;

            while (!ToReach.empty()) {

                SubGraphInfo sg;
                st.push(&*ToReach.top());
                ToReach.pop();

                assert(st.top()->sid==-1);

                sg.root=st.top();
                sg.sid=cnt;
                sg.size=0;
                st.top()->sid=cnt;

                while (!st.empty()) {

                    AlignGlobal::Node *cur=st.top();
                    st.pop();
                    ++sg.size;

                    assert(cur->sid==cnt);
					// std::printf("Visiting node %2i %2i\n",cur->id,cur->sid);

                    for (auto li = std::begin(cur->Adj); li != std::end(cur->Adj); ++li) {

                        if ((*li)->Adj(cur)->sid == -1) {
                            (*li)->Adj(cur)->sid=cnt;
                            st.push((*li)->Adj(cur));
                        }
                        else {
                            assert((*li)->Adj(cur)->sid == cnt);
                        }
                    }

                }

                cnt++;
                CC.push_back(sg);

                while (!ToReach.empty() && ToReach.top()->sid != -1) ToReach.pop();
            }

            return cnt;
        }

        void Clear() {

            for (auto li = std::begin(A); li != std::end(A); ++li) {
                delete (*li);
            }

            N.clear();
            A.clear();
        }

        void MakeAllDormant() {
            for (auto li = std::begin(N); li != std::end(N); ++li) {
                (*li).Active=false;
            }
        }

        bool GlobalAlign(const std::map<int, std::string> &Names, const double epsilon, int maxiter, bool Rigid, FILE *elfp, vcg::CallBackPos* cb) {

            double change;
            int step = 0, localmaxiter;

            if (cb != NULL) cb(0, "Global Alignment...");
            AlignGlobal::LOG(elfp,"----------------\n----------------\nGlobalAlignment (target eps %7.3f)\n", epsilon);

            std::queue<AlignGlobal::Node *> Q;
            MakeAllDormant();
            AlignGlobal::Node *curr = ChooseDormantWithMostDormantLink();
            curr->Active = true;

            int cursid = curr->sid;
            AlignGlobal::LOG(elfp, "Root node %i '%s' with %i dormant link\n", curr->id, Names.find(curr->id)->second.c_str(), curr->DormantAdjNum());

            while (DormantNum() > 0) {

                AlignGlobal::LOG(elfp,"---------\nGlobalAlignment loop DormantNum = %i\n", DormantNum());

                curr = ChooseDormantWithMostActiveLink();
                if (!curr) {
                    // la componente connessa e' finita e si passa alla successiva cercando un dormant con tutti dormant.
                    AlignGlobal::LOG(elfp,"\nCompleted Connected Component %i\n", cursid);
                    AlignGlobal::LOG(elfp,"\nDormant Num: %i\n", DormantNum());

                    curr = ChooseDormantWithMostDormantLink();
                    if (curr == nullptr) {
                        AlignGlobal::LOG(elfp,"\nFailed ChooseDormantWithMostDormantLink, chosen id:%i\n" ,0);
                        break; // non ci sono piu' componenti connesse composte da piu' di una singola mesh.
                    }
                    else {
                        AlignGlobal::LOG(elfp,"\nCompleted ChooseDormantWithMostDormantLink, chosen id:%i\n" ,curr->id);
                    }

                    curr->Active = true;
                    cursid = curr->sid;
                    curr = ChooseDormantWithMostActiveLink ();
                    if (curr == nullptr) {
                        AlignGlobal::LOG(elfp, "\nFailed    ChooseDormantWithMostActiveLink, chosen id:%i\n", 0);
                    }
                    else {
                        AlignGlobal::LOG(elfp, "\nCompleted ChooseDormantWithMostActiveLink, chosen id:%i\n", curr->id);
                    }
                }

                AlignGlobal::LOG(elfp,"\nAdded node %i '%s' with %i/%i Active link\n",curr->id,Names.find(curr->id)->second.c_str(),curr->ActiveAdjNum(),curr->Adj.size());
                curr->Active=true;
                curr->Queued=true;

                // Si suppone, ad occhio, che per risistemare un insieme di n mesh servano al piu' 10n passi;
                localmaxiter = ActiveNum() * 10;
                Q.push(curr);
                step = 0;

                // Ciclo interno di allineamento
                while (!Q.empty()) {

                    curr = Q.front();
                    Q.pop();

                    curr->Queued=false;
                    change = curr->AlignWithActiveAdj(Rigid);
                    step++;

                    AlignGlobal::LOG(elfp, "     Step %5i Queue size %5i Moved %4i  err %10.4f\n", step, Q.size(), curr->id, change);

                    if (change > epsilon) {

                        curr->PushBackActiveAdj(Q);
                        AlignGlobal::LOG(elfp,"         Large Change pushing back active nodes adj to %i to Q (new size %i)\n",curr->id,Q.size());

                        if (change > epsilon * 1000) {
                            std::printf("Large Change Warning\n\n");
                        }
                    }
                    if (step > localmaxiter) return false;
                    if (step > maxiter) return false;
                }
            }

            if (!curr) {

                AlignGlobal::LOG(elfp,"Alignment failed for %i meshes:\n",DormantNum());
                for (auto li = std::begin(N); li != std::end(N); ++li){
                    if (!(*li).Active) {
                        //(*li).M.SetIdentity();
                        (*li).Discarded=true;
                        AlignGlobal::LOG(elfp, "%5i\n", (*li).id);
                    }
                }
            }

            AlignGlobal::LOG(elfp,"Completed Alignment in %i steps with error %f\n",step,epsilon);
            return true;
        }

        AlignGlobal::Node* ChooseDormantWithMostDormantLink() {

            int MaxAdjNum = 0;
            AlignGlobal::Node *BestNode = nullptr;

            for (auto li = std::begin(N); li != std::end(N); ++li) {
                if (!(*li).Active) {
                    int AdjNum = (*li).DormantAdjNum();
                    if (AdjNum > MaxAdjNum) {
                        MaxAdjNum = AdjNum;
                        BestNode = &(*li);
                    }
                }
            }

            if (!BestNode){
                std::printf("Warning! Unable to find a Node with at least a dormant link!!\n");
                return nullptr;
            }

            assert(BestNode);
            assert(!BestNode->Queued);
            assert(!BestNode->Active);

            return BestNode;
        }

        AlignGlobal::Node* ChooseDormantWithMostActiveLink() {

            int MaxAdjNum = 0;
            AlignGlobal::Node* BestNode = nullptr;

            for (auto li = std::begin(N); li != std::end(N); ++li) {
                if (!(*li).Active) {
                    int AdjNum = (*li).ActiveAdjNum();
                    if (AdjNum > MaxAdjNum) {
                        MaxAdjNum = AdjNum;
                        BestNode = &(*li);
                    }
                }
            }

            if (!BestNode){
                // Abbiamo finito di sistemare questa componente connessa.
                std::printf("Warning! Unable to find a Node with at least an active link!!\n");
                return nullptr;
            }

            assert(BestNode);
            assert(!BestNode->Queued);
            assert(!BestNode->Active);

            return BestNode;
        }

        void BuildGraph(std::vector<AlignPair::Result *> &Res, std::vector<Matrix44d> &Tr, std::vector<int> &Id) {

            Clear();
            // si suppone che la matrice Tr[i] sia relativa ad un nodo con id Id[i];
            int mn = static_cast<int>(Tr.size());

            //	printf("building graph\n");
            AlignGlobal::Node rgn;
            rgn.Active = false;
            rgn.Queued = false;
            rgn.Discarded = false;

            std::map<int, AlignGlobal::Node *> Id2N;
            std::map<int, int> Id2I;

            for (int i = 0; i < mn; ++i) {
                rgn.id = Id[i];
                rgn.M = Tr[i];
                N.push_back(rgn);
                Id2N[rgn.id] = &(N.back());
                Id2I[rgn.id] = i;
            }

            std::printf("building %zu graph arcs\n",Res.size());
            AlignGlobal::VirtAlign *tv;

            // Ciclo principale in cui si costruiscono i vari archi
            // Si assume che i result siano fatti nel sistema di riferimento della matrici fix.

            for (auto rii = std::begin(Res); rii != std::end(Res); ++rii) {

                AlignPair::Result *ri = *rii;
                tv = new VirtAlign();
                tv->Fix = Id2N[(*ri).FixName];
                tv->Mov = Id2N[(*ri).MovName];
                tv->Fix->Adj.push_back(tv);
                tv->Mov->Adj.push_back(tv);
                tv->FixP = (*ri).Pfix;
                tv->MovP = (*ri).Pmov;
                tv->FixN = (*ri).Nfix;
                tv->MovN = (*ri).Nmov;

                /*

                    Siano:
                    Pf e Pm   i punti sulle mesh fix e mov nei sist di rif originali
                    Pft e Pmt i punti sulle mesh fix e mov dopo le trasf correnti;
                    Mf e Mm  le trasf che portano le  mesh fix e mov nelle posizioni correnti;
                    If e Im  le trasf inverse di cui sopra
                    Vale:
                    Pft = Mf*Pf  e Pmt = Mm*Pm
                     Pf = If*Pft e Pm = Im*Pmt

                                 Res *   Pm     =         Pf;
                                 Res * Im * Pmt =      If * Pft
                        Mf * Res * Im * Pmt = Mf * If * Pft
                    (Mf * Res * Im) * Pmt = Pft

                */
                Matrix44d Mm = Tr[Id2I[(*ri).MovName]];
                Matrix44d Mf = Tr[Id2I[(*ri).FixName]];
                Matrix44d Im = Inverse(Mm);
                Matrix44d If = Inverse(Mf);

                Matrix44d NewTr = Mf * (*ri).Tr * Im; // --- orig

                tv->M2F = NewTr;
                tv->F2M = Inverse(NewTr);

                assert(tv->Check());
                A.push_back(tv);
            }

            ComputeConnectedComponents();
        }

        bool GetMatrixVector(std::vector<Matrix44d> &Tr, std::vector<int> &Id) {

            std::map<int, AlignGlobal::Node *> Id2N;

            Tr.clear();

            for (auto li = std::begin(N); li != std::end(N); ++li) {
                Id2N[(*li).id] = &*li;
            }

            Tr.resize(Id.size());

            for (std::size_t i = 0; i < Id.size(); ++i) {

                if (Id2N[Id[i]] == 0) return false;
                Tr[i] = Id2N[Id[i]]->M;
            }

            return false;
        }

    };
}

#endif //MESHLAB_ALIGNGLOBAL_H
