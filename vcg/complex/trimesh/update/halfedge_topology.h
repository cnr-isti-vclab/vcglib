#ifndef VCG_HEDGE_TOPOLOGY
#define VCG_HEDGE_TOPOLOGY

#include <vcg/connectors/halfedge_pos.h>
#include <vcg/complex/trimesh/allocate.h>

#include <vector>
#include <algorithm>

using namespace std;
using namespace vcg::hedge;
using namespace vcg::tri;

namespace vcg
{
    namespace tri
    {
        /*!
          * \brief Class containing functions to modify the topology of a halfedge based mesh
          *
          */
        template <class MeshType> class HalfEdgeTopology
        {
            public:

                typedef typename MeshType::VertexPointer VertexPointer;
                typedef typename MeshType::EdgePointer EdgePointer;
                typedef typename MeshType::HEdgePointer HEdgePointer;
                typedef typename MeshType::FacePointer FacePointer;

                typedef typename MeshType::VertexIterator VertexIterator;
                typedef typename MeshType::EdgeIterator EdgeIterator;
                typedef typename MeshType::HEdgeIterator HEdgeIterator;
                typedef typename MeshType::FaceIterator FaceIterator;

                /*!
                  * Collpases an edge shared by two quads, generating only quads.
                  * Made by a series of a vertex rotation and a diagonal collapse.
                  *
                  * \param m Mesh
                  * \param ep Edge to be collapsed
                  * \param vp Vertex that will be rotated
                  *
                  * \return Pointer to the new vertex
                  */
                static VertexPointer edge_collapse_quad(MeshType &m, EdgePointer ep, VertexPointer vp)
                {
                    assert(vp);
                    assert(ep);
                    assert(MeshType::EdgeType::HasEHAdjacency());
                    assert(MeshType::HEdgeType::HasHEAdjacency());
                    assert(MeshType::HEdgeType::HasHVAdjacency());
                    assert(ep->EHp()->HVp() == vp || ep->EHp()->HOp()->HVp() == vp);
                    assert(ep->EHp()->HFp()->VN() == 4);
                    assert(ep->EHp()->HOp()->HFp()->VN() == 4);

                    VertexPointer vp_opp = ep->EHp()->HOp()->HVp();

                    VertexPointer vp_rot = vertex_rotate( m, vp );

                    assert(vp_rot == vp);

                    FacePointer fp;

                    //retrieve right face
                    Pos<MeshType> p(vp->VHp(),true);

                    while(p.HE()->HNp->HOp()->HVp() != vp_opp)
                    {
                        p.FlipE();
                        p.FlipF();
                    }

                    fp = p.F();


                    return diagonal_collapse( m, fp, vp );

                }

                /*!
                  * Collpases a diagonal in a quad.
                  *
                  *
                  * \param m Mesh
                  * \param fp Face where diagonal resides
                  * \param vp One of the two vertices of the diagonal
                  *
                  * \return Pointer to the new vertex
                  */
                static VertexPointer diagonal_collapse(MeshType &m, FacePointer fp, VertexPointer vp)
                {

                    assert(MeshType::VertexType::HasVHAdjacency());
                    assert(MeshType::EdgeType::HasEHAdjacency());
                    assert(MeshType::FaceType::HasFHAdjacency());
                    assert(MeshType::HEdgeType::HasHVAdjacency());
                    assert(MeshType::HEdgeType::HasHEAdjacency());
                    assert(MeshType::HEdgeType::HasHFAdjacency());
                    assert(MeshType::HEdgeType::HasHOppAdjacency());
                    assert(MeshType::HEdgeType::HasHPrevAdjacency());

                    assert(fp);
                    assert(fp->FHp());
                    assert(fp->VN() == 4);

                    if( !can_remove_face(fp) )
                        return NULL;

                    HEdgePointer hp;

                    vector<VertexPointer> vps = getVertices(fp);
                    VertexPointer opp_vert = NULL;

                    for(unsigned int i = 0; i< vps.size(); i++)
                        if(vps[i] == vp)
                            opp_vert = vps[(i+2)%vps.size()];

                    assert(opp_vert);

                    if( fp->FHp()->HVp() == vp || fp->FHp()->HVp() == opp_vert)
                        hp = fp->FHp();
                    else
                        hp = fp->FHp()->HNp();

                    vector<HEdgePointer> hps = getHEdges(fp,hp);

                    int edge_cnt = 0;

                    if(hps[0]->HOp()->HFp() || hps[1]->HOp()->HFp())
                        edge_cnt++;
                    if(hps[2]->HOp()->HFp() || hps[3]->HOp()->HFp())
                        edge_cnt++;

                    VertexIterator vi;

                    if(edge_cnt > 0)
                    {
                        typename Allocator<MeshType>::template PointerUpdater<VertexPointer> puv;

                        if(m.vert.empty())
                            puv.oldBase = 0;
                        else
                        {
                            puv.oldBase = &*(m.vert.begin());
                            puv.oldEnd = &m.vert.back()+1;
                        }

                        vi = Allocator<MeshType>::AddVertices(m,1);

                        EdgeIterator ei = Allocator<MeshType>::AddEdges(m,edge_cnt);

                        puv.newBase = &*(m.vert.begin());
                        puv.newEnd = &m.vert.back()+1;

                        if( puv.NeedUpdate() )
                        {
                            puv.Update(vp);
                            puv.Update(opp_vert);
                            for(typename vector<VertexPointer>::iterator vpi = vps.begin(); vpi != vps.end(); ++vpi)
                                puv.Update(*vpi);
                        }

                        typename Allocator<MeshType>::template PointerUpdater<HEdgePointer> puh;

                        if(m.hedge.empty())
                            puh.oldBase = 0;
                        else
                        {
                            puh.oldBase = &*(m.hedge.begin());
                            puh.oldEnd = &m.hedge.back()+1;
                        }

                        HEdgeIterator hi = Allocator<MeshType>::AddHEdges(m,2*edge_cnt);

                        puh.newBase = &*(m.hedge.begin());
                        puh.newEnd = &m.hedge.back()+1;

                        if( puh.NeedUpdate() )
                            for(typename vector<HEdgePointer>::iterator hpi = hps.begin(); hpi != hps.end(); ++hpi)
                                    puh.Update(*hpi);



                        HEdgeIterator hi1 = hi;
                        HEdgeIterator hi2 = hi;
                        ++hi2;

                        (*vi).VHp() = &(*hi1);

                        change_vertex( hps[0]->HVp(), &(*vi));
                        change_vertex( hps[2]->HVp(), &(*vi));

                        for( int count = 0; count < 2; count++ )
                        {

                            int i = 2*count;

                            FacePointer fp1 = hps[i+1]->HOp()->HFp();
                            FacePointer fp2 = hps[i]->HOp()->HFp();

                            if( fp1 || fp2 )
                            {

                                // HOp
                                (*hi1).HOp() = &(*hi2);
                                (*hi2).HOp() = &(*hi1);

                                // EH
                                (*ei).EHp() = &(*hi1);

                                // HE
                                (*hi1).HEp() = &(*ei);
                                (*hi2).HEp() = &(*ei);

                                // HV
                                (*hi1).HVp() = &(*vi);
                                (*hi2).HVp() = hps[i+1]->HVp();

                                // FH
                                if( fp1 && ( fp1->FHp() == hps[i+1]->HOp() ) )
                                    fp1->FHp() = &(*hi1);

                                if( fp2 && ( fp2->FHp() == hps[i]->HOp() ) )
                                    fp2->FHp() = &(*hi2);

                                //HF
                                (*hi1).HFp() = fp1;
                                (*hi2).HFp() = fp2;

                                //HNp
                                (*hi1).HNp() = hps[i+1]->HOp()->HNp();
                                (*hi2).HNp() = hps[i]->HOp()->HNp();

                                (*hi1).HNp()->HPp() = &(*hi1);
                                (*hi2).HNp()->HPp() = &(*hi2);

                                //HPp
                                (*hi1).HPp() = hps[i+1]->HOp()->HPp();
                                (*hi2).HPp() = hps[i]->HOp()->HPp();

                                (*hi1).HPp()->HNp() = &(*hi1);
                                (*hi2).HPp()->HNp() = &(*hi2);

                                //VH
                                VertexPointer tmp = hps[i+1]->HVp();
                                if( tmp->VHp() == hps[i+1] ||  tmp->VHp() == hps[i]->HOp() )
                                    tmp->VHp() = &(*hi2);

                                ++ei;

                                ++hi1;
                                ++hi1;

                                ++hi2;
                                ++hi2;

                            }
                            else
                            {
                                hps[i+1]->HOp()->HPp()->HNp() = hps[i]->HOp()->HNp();
                                hps[i]->HOp()->HNp()->HPp() = hps[i+1]->HOp()->HPp();

                                hps[i+1]->HVp()->VHp() = NULL;
                            }


                        }

                    }


                    Allocator<MeshType>::DeleteFace(m, *(fp) );
                    Allocator<MeshType>::DeleteVertex(m, *(vp) );
                    Allocator<MeshType>::DeleteVertex(m, *(opp_vert) );

                    for(typename vector<HEdgePointer>::iterator hpi = hps.begin(); hpi != hps.end(); ++hpi)
                    {
                        if(! (*hpi)->HEp()->IsD() )
                            Allocator<MeshType>::DeleteEdge(m, *((*hpi)->HEp()) );

                        if(! (*hpi)->IsD())
                        {
                            Allocator<MeshType>::DeleteHEdge(m, *(*hpi) );
                            Allocator<MeshType>::DeleteHEdge(m, *((*hpi)->HOp()) );
                        }
                    }

                    if(edge_cnt > 0)
                        return &(*vi);

                    for(typename vector<VertexPointer>::iterator vpi = vps.begin(); vpi != vps.end(); ++vpi)
                        if(!(*vpi)->IsD())
                            Allocator<MeshType>::DeleteVertex(m, *(*vpi) );

                    return NULL;
                }

                /*!
                  * Removes a doublet merging the two quads in one
                  *
                  * \param m Mesh
                  * \param vp Vertex shared by the two consecutive edges of the doublet
                  *
                  * \return Pointer to the new face
                  */
                static FacePointer doublet_remove(MeshType &m, VertexPointer vp)
                {
                    assert(vp);

                    HEdgePointer hp = vp->VHp();

                    assert(hp);

                    FacePointer fp1 = hp->HFp();
                    FacePointer fp2 = hp->HOp()->HFp();

                    // check if face is a doublet

                    assert( fp1 );
                    assert( fp1 == hp->HPp()->HFp() );

                    assert( fp2 );
                    assert( fp2 == hp->HOp()->HNp()->HFp() );

                    assert( hp->HOp()->HNp()->HOp() == hp->HPp() );

                    assert( fp1->VN() == 4);
                    assert( fp2->VN() == 4);

                    // end of check

                    vector<VertexPointer> vert_face1 = getVertices(fp1, hp);
                    vector<VertexPointer> vert_face2 = getVertices(fp2, hp->HOp()->HNp());

                    remove_face_unsafe(m, fp1 );
                    remove_face_unsafe(m, fp2 );

                    Allocator<MeshType>::DeleteVertex(m, *vp);


                    vector<VertexPointer> new_face_vert;

                    new_face_vert.push_back( vert_face1[1] );
                    new_face_vert.push_back( vert_face1[2] );
                    new_face_vert.push_back( vert_face1[3] );
                    new_face_vert.push_back( vert_face2[2] );

                    return add_face_unsafe(m, new_face_vert );

                }

                /*!
                  * Removes a singlet replacing it with an edge
                  *
                  * \param m Mesh
                  * \param vp Vertex shared by the two consecutive edges inside the singlet
                  *
                  * \return Pointer to the new edge
                  */
                static EdgePointer singlet_remove(MeshType &m, VertexPointer vp)
                {
                    assert( vp );

                    HEdgePointer hp = vp->VHp();
                    assert( hp );

                    FacePointer fp1 = hp->HFp();
                    FacePointer fp2 = hp->HOp()->HFp();

                    assert( fp1 && fp2 && fp1 == fp2 ); // the faces pointed by the halfedges must be the same

                    HEdgePointer hp1 = hp->HNp()->HOp();
                    HEdgePointer hp2 = hp->HNp()->HNp()->HOp();

                    // pointers to the near faces
                    FacePointer fp3 = hp1->HFp();
                    FacePointer fp4 = hp2->HFp();

                    Allocator<MeshType>::DeleteFace( m, *fp1 );

                    Allocator<MeshType>::DeleteEdge( m, *(hp->HEp()) );

                    Allocator<MeshType>::DeleteHEdge( m, *hp );
                    Allocator<MeshType>::DeleteHEdge( m, *(hp->HOp()) );

                    Allocator<MeshType>::DeleteVertex(m, *(vp) );

                    Allocator<MeshType>::DeleteEdge(m, *(hp1->HEp()) );
                    Allocator<MeshType>::DeleteEdge(m, *(hp2->HEp()) );

                    Allocator<MeshType>::DeleteHEdge( m, *(hp1->HOp()) );
                    Allocator<MeshType>::DeleteHEdge( m, *(hp2->HOp()) );

                    if(!fp3 && !fp4) // there are no faces, nothing has to be created
                    {
                        Allocator<MeshType>::DeleteHEdge( m, *hp1 );
                        Allocator<MeshType>::DeleteHEdge( m, *hp2 );

                        return NULL;
                    }

                    EdgeIterator ei = Allocator<MeshType>::AddEdges(m,1);

                    (*ei).EHp() = hp1;

                    hp1->HEp() = &(*ei);
                    hp2->HEp() = &(*ei);

                    hp1->HOp() = hp2;
                    hp2->HOp() = hp1;

                    return &(*ei);

                }

                /*!
                  * Rotates a non-border edge shared by two quads
                  *
                  * \param m Mesh
                  * \param ep Edge to be rotated
                  * \param cw flag denoting a clockwise or counter-clockwise rotation
                  *
                  * \return Pointer to the rotated edge
                  */
                static EdgePointer edge_rotate(MeshType &m, EdgePointer ep, bool cw)
                {
                    assert( MeshType::EdgeType::HasEHAdjacency() );
                    assert( MeshType::HEdgeType::HasHFAdjacency() );
                    assert( MeshType::HEdgeType::HasHOppAdjacency() );
                    assert( MeshType::FaceType::HasFHAdjacency() );

                    assert( ep->EHp()->HFp() );

                    assert( ep->EHp()->HFp()->VN() == 4 );

                    assert( ep->EHp()->HOp()->HFp() );
                    assert( ep->EHp()->HOp()->HFp()->VN() == 4 );

                    FacePointer fp1 = ep->EHp()->HFp();
                    FacePointer fp2 = ep->EHp()->HOp()->HFp();

                    vector<VertexPointer> old_face1 = getVertices( fp1, ep->EHp() );
                    vector<VertexPointer> old_face2 = getVertices( fp2, ep->EHp()->HOp() );

                    remove_face_unsafe(m, fp1);
                    remove_face_unsafe(m, fp2);

                    vector<VertexPointer> new_face1;
                    vector<VertexPointer> new_face2;

                    if(cw)
                    {
                        new_face1.push_back( old_face1[3] );
                        new_face1.push_back( old_face2[3] );
                        new_face1.push_back( old_face1[1] );
                        new_face1.push_back( old_face1[2] );

                        new_face2.push_back( old_face2[3] );
                        new_face2.push_back( old_face1[3] );
                        new_face2.push_back( old_face2[1] );
                        new_face2.push_back( old_face2[2] );

                    }
                    else
                    {
                        new_face1.push_back( old_face1[2] );
                        new_face1.push_back( old_face1[3] );
                        new_face1.push_back( old_face1[0] );
                        new_face1.push_back( old_face2[2] );

                        new_face2.push_back( old_face2[2] );
                        new_face2.push_back( old_face2[3] );
                        new_face2.push_back( old_face2[0] );
                        new_face2.push_back( old_face1[2] );
                    }

                    fp1 = add_face_unsafe(m, new_face1);
                    fp2 = add_face_unsafe(m, new_face2);


                    // retrieve inserted edge
                    Pos<MeshType> p1( fp1->FHp(), true);

                    while( p1.HE()->HOp()->HFp() != fp2 )
                    {
                        p1.FlipV();
                        p1.FlipE();
                    }

                    return p1.E();

                }

                /*!
                  * Rotates a non-border vertex shared by only quads
                  *
                  * \param m Mesh
                  * \param vp Vertex to be rotated
                  *
                  * \return Pointer to the rotated vertex
                  */
                static VertexPointer vertex_rotate(MeshType &m, VertexPointer vp)
                {

                    assert(MeshType::VertexType::HasVHAdjacency());

                    vector<FacePointer> old_faces;

                    typedef vector<VertexPointer> vert_vect;
                    vector< vert_vect > old_face_verts;

                    assert( vp->VHp() );

                    Pos<MeshType> p(vp->VHp(), true);

                    HEdgePointer hep = p.HE();

                    do
                    {
                        assert( p.F() );
                        assert( p.F()->VN() == 4);

                        old_faces.push_back(p.F());

                        old_face_verts.push_back( getVertices( p.F(), p.HE() ) );

                        p.FlipE();
                        p.FlipF();

                    }while(p.HE() != hep);

                    assert( old_faces.size() == old_face_verts.size() );

                    int size = old_faces.size();

                    vector<vert_vect> new_face_verts;

                    for(int i = 0; i < size; i++)
                    {
                        new_face_verts.push_back(vector<VertexPointer>());

                        new_face_verts.back().push_back( old_face_verts[i][0] );
                        new_face_verts.back().push_back( old_face_verts[i][2] );
                        new_face_verts.back().push_back( old_face_verts[i][3] );
                        new_face_verts.back().push_back( old_face_verts[(i+1)%size][2] );

                    }

                    for(typename vector<FacePointer>::iterator fi = old_faces.begin(); fi != old_faces.end(); ++fi)
                        remove_face_unsafe( m, *fi );

                    for(typename vector<vert_vect>::iterator vi = new_face_verts.begin(); vi != new_face_verts.end(); ++vi)
                        add_face_unsafe(m, *vi);

                    return vp;

                }

                /*!
                  * Collapses a generic edge
                  *
                  * \param m Mesh
                  * \param ep Edge to be collapsed
                  * \param vp Vertex to be deleted
                  *
                  * \return Pointer to the other vertex belonging to the collapsed edge
                  */
                static VertexPointer edge_collapse(MeshType &m, EdgePointer ep, VertexPointer vp)
                {

                    assert(MeshType::EdgeType::HasEHAdjacency());
                    assert(MeshType::VertexType::HasVHAdjacency());
                    assert(MeshType::HEdgeType::HasHOppAdjacency());
                    assert(MeshType::HEdgeType::HasHVAdjacency());
                    assert(MeshType::HEdgeType::HasHPrevAdjacency());

                    if( ep->EHp()->HFp() )
                        assert(ep->EHp()->HFp()->VN() > 3);

                    if( ep->EHp()->HOp()->HFp())
                        assert(ep->EHp()->HOp()->HFp()->VN() > 3);

                    assert(ep->EHp()->HFp() || ep->EHp()->HOp()->HFp());
                    assert(ep->EHp()->HVp() == vp || ep->EHp()->HOp()->HVp() == vp);


                    HEdgePointer he = ep->EHp();
                    HEdgePointer hopp = he->HOp();

                    VertexPointer vp1;

                    if( he->HVp() == vp )
                        vp1 = hopp->HVp();
                    else
                        vp1 = he->HVp();

                    change_vertex( vp, vp1);

                    //HP
                    he->HNp()->HPp() = he->HPp();
                    hopp->HNp()->HPp() = hopp->HPp();

                    //HN
                    he->HPp()->HNp() = he->HNp();
                    hopp->HPp()->HNp() = hopp->HNp();

                    //FH
                    if( he->HFp() )
                        if( he->HFp()->FHp() == he )
                            he->HFp()->FHp() = he->HNp();

                    if( hopp->HFp() )
                        if( hopp->HFp()->FHp() == hopp )
                            hopp->HFp()->FHp() = hopp->HNp();

                    // VH
                    if( vp1->VHp() == hopp )
                        vp1->VHp() = hopp->HNp();

                    Allocator<MeshType>::DeleteEdge(m,*ep);
                    Allocator<MeshType>::DeleteHEdge(m,*he);
                    Allocator<MeshType>::DeleteHEdge(m,*hopp);
                    Allocator<MeshType>::DeleteVertex(m,*vp);

                    return vp1;

                }

                /*!
                  * Adds a face in a mesh, checking if the operation is possible.
                  *
                  * \param m Mesh
                  * \param vps Vector of vertices (in ccw order) that will belong to the new face
                  *
                  * \return Pointer to the new face if it has been inserted, NULL otherwise
                  */
                static FacePointer add_face(MeshType &m, vector<VertexPointer> &vps)
                {

                    assert(MeshType::VertexType::HasVHAdjacency());
                    assert(MeshType::EdgeType::HasEHAdjacency());
                    assert(MeshType::HEdgeType::HasHVAdjacency());
                    assert(MeshType::HEdgeType::HasHEAdjacency());
                    assert(MeshType::HEdgeType::HasHFAdjacency());
                    assert(MeshType::HEdgeType::HasHOppAdjacency());
                    assert(MeshType::HEdgeType::HasHPrevAdjacency());

                    unsigned int size = vps.size();

                    assert(size >= 3); //there must be at least 3 vertices

                    for(unsigned int i = 0; i< size; i++)
                    {
                        // all vertices must be different
                        assert( count(vps.begin(), vps.end(), vps[i]) == 1 );
                    }

                    vector<HEdgePointer> hps;

                    while(hps.size() < size)
                        if( !can_add_hedge(vps, hps) )
                            return NULL;

                    vector<bool> non_manifold_vertices(size, false);
                    return add_face_unsafe( m,vps, hps, non_manifold_vertices );

                }

                /*!
                  * Removes a face in a mesh, checking if the operation is possible
                  *
                  * \param m Mesh
                  * \param fp face to be removed
                  *
                  * \retval true if face has been removed
                  * \retval false otherwise
                  */
                static bool remove_face(MeshType &m, FacePointer fp)
                {

                    assert(MeshType::VertexType::HasVHAdjacency());
                    assert(MeshType::EdgeType::HasEHAdjacency());
                    assert(MeshType::FaceType::HasFHAdjacency());
                    assert(MeshType::HEdgeType::HasHVAdjacency());
                    assert(MeshType::HEdgeType::HasHEAdjacency());
                    assert(MeshType::HEdgeType::HasHFAdjacency());
                    assert(MeshType::HEdgeType::HasHOppAdjacency());
                    assert(MeshType::HEdgeType::HasHPrevAdjacency());

                    if( can_remove_face(fp) )
                    {
                        remove_face_unsafe(m, fp);
                        return true;
                    }

                    return false;
                }

            protected:

                /*!
                  * Adds a face in a mesh without any check
                  *
                  * \param m Mesh
                  * \param vps Vector of vertices (in ccw order) that will belong to the new face
                  *
                  * \return Pointer to the new face
                  */
                static FacePointer add_face_unsafe(MeshType &m, vector<VertexPointer> &vps)
                {
                    unsigned int size = vps.size();

                    vector<HEdgePointer> hps;
                    vector<bool> non_manifold_vertices;

                    while(hps.size() < size)
                    {
                        if( can_add_hedge(vps, hps) )
                            non_manifold_vertices.push_back( false );
                        else
                            non_manifold_vertices.push_back( hps.back() == NULL );
                    }

                    return add_face_unsafe(m,vps,hps, non_manifold_vertices);
                }

                /*!
                  * Adds a face in a mesh without any check
                  *
                  * \param m Mesh
                  * \param vps Vector of vertices (in ccw order) that will belong to the new face
                  * \param non_manifold_vertices Vector of booleans denoting on the i-th position if the i-th vertex is non-manifold
                  *
                  * \return Pointer to the new face
                  */
                static FacePointer add_face_unsafe(MeshType &m, vector<VertexPointer> &vps, vector<HEdgePointer> &hps, vector<bool> &non_manifold_vertices)
                {

                    assert(MeshType::VertexType::HasVHAdjacency());
                    assert(MeshType::EdgeType::HasEHAdjacency());
                    assert(MeshType::HEdgeType::HasHVAdjacency());
                    assert(MeshType::HEdgeType::HasHEAdjacency());
                    assert(MeshType::HEdgeType::HasHFAdjacency());
                    assert(MeshType::HEdgeType::HasHOppAdjacency());
                    assert(MeshType::HEdgeType::HasHPrevAdjacency());

                    unsigned int size = vps.size();

                    assert(size >= 3); //there must be at least 3 vertices

//                    for(unsigned int i = 0; i< size; i++)
//                    {
//                        // all vertices must be different
//                        assert( count(vps.begin(), vps.end(), vps[i]) == 1 );
//                    }

                    HEdgeIterator hi;

                    assert(hps.size() == size);

                    HEdgePointer nullPointer = NULL;
                    int edge_n = count(hps.begin(), hps.end(), nullPointer);

                    FaceIterator fi = Allocator<MeshType>::AddFaces(m,1);
                    m.face.back().Alloc( size );

                    if(edge_n > 0)
                    {
                        (*fi).SetD();

                        EdgeIterator ei = Allocator<MeshType>::AddEdges(m,edge_n);

                        for(EdgeIterator ei1 = ei; ei1 != m.edge.end(); ++ei1)
                            (*ei1).SetD();

                        typename Allocator<MeshType>::template PointerUpdater<HEdgePointer> pu;

                        if(m.hedge.empty())
                            pu.oldBase = 0;
                        else
                        {
                            pu.oldBase = &*(m.hedge.begin());
                            pu.oldEnd = &m.hedge.back()+1;
                        }

                        hi = Allocator<MeshType>::AddHEdges(m,2*edge_n);

                        pu.newBase = &*(m.hedge.begin());
                        pu.newEnd = &m.hedge.back()+1;


                        //undelete face
                        (*fi).ClearD();

                        //undelete edges
                        for(EdgeIterator ei1 = ei; ei1 != m.edge.end(); ++ei1)
                            (*ei1).ClearD();


                        HEdgeIterator hi1 = hi;
                        HEdgeIterator hi2 = hi;

                        ++hi2;


                        for(EdgeIterator ei1 = ei; ei1 != m.edge.end(); ++ei1, ++hi1, ++hi2)
                        {
                            // EH
                            (*ei1).EHp() = &(*hi1);

                            // HE
                            (*hi1).HEp() = &(*ei1);
                            (*hi2).HEp() = &(*ei1);

                            //HO
                            (*hi1).HOp() = &(*hi2);
                            (*hi2).HOp() = &(*hi1);

                            // HF
                            (*hi1).HFp() = &(*fi);

                            ++hi1;
                            ++hi2;
                        }

                        // update hedge pointers (if needed)
                        if( pu.NeedUpdate() )
                            for(typename vector<HEdgePointer>::iterator hpsi = hps.begin(); hpsi != hps.end(); ++hpsi)
                            {
                                if((*hpsi))
                                    pu.Update(*hpsi);
                            }



                    }

                    vector<HEdgePointer> hps1;

                    for(unsigned int i = 0; i < size; i++)
                    {
                        if(hps[i] == NULL)
                        {
                            hps1.push_back(&(*hi));
                            ++hi;
                            ++hi;
                        }
                        else
                            hps1.push_back(hps[i]);
                    }


                    assert( hps1.size() == size );


                    for(unsigned int i = 0; i < size; i++)
                    {

                        int next = (i+1)%size;

                        // hedge already exisitng
                        if(hps[i])
                        {
                            hps1[i]->HFp() = &(*fi);

                            // next hedge was disconnected
                            if(!hps[next])
                            {

                                hps1[next]->HOp()->HNp() = hps1[i]->HNp();

                                hps1[i]->HNp()->HPp() = hps1[next]->HOp();

                                hps1[i]->HNp() = hps1[next];

                                hps1[next]->HPp() = hps1[i];
                            }
                        }

                        // hedge wasn't existing, vertex was disconnected
                        else
                        {
                            //HV
                            hps1[i]->HVp() = vps[i];
                            hps1[i]->HOp()->HVp() = vps[next];


                            hps1[i]->HNp() = hps1[next];

                            // next hedge was existing (vertex was disconnected)
                            if(hps[next])
                            {
                                hps1[i]->HOp()->HPp() = hps1[next]->HPp();
                                hps1[next]->HPp()->HNp() = hps1[i]->HOp();
                            }

                            //vertex was detached
                            else
                            {
                                // after face insertion vertex will become non-manifold
                                if(non_manifold_vertices[next])
                                {
                                    Pos<MeshType> p(vps[next]->VHp(), true);

                                    while(p.F())
                                    {

                                        p.FlipE();
                                        p.FlipF();

                                        if(p.HE() == vps[next]->VHp())
                                            assert(0); //can't add a connection, there is no space
                                    }


                                    p.HE()->HPp()->HNp() = hps1[i]->HOp();
                                    hps1[i]->HOp()->HPp() = p.HE()->HPp();

                                    p.HE()->HPp() = hps1[next]->HOp();
                                    hps1[next]->HOp()->HNp() = p.HE();

                                }
                                else
                                {
                                    hps1[i]->HOp()->HPp() = hps1[next]->HOp();
                                    hps1[next]->HOp()->HNp() = hps1[i]->HOp();
                                }

                            }


                            hps1[next]->HPp() = hps1[i];

                            //VH
                            if( !vps[i]->VHp())
                                vps[i]->VHp() = hps1[i];
                        }
                    }

                    //FH
                    (*fi).FHp() = hps1.front();

                    return &(*fi);

                }

                /*!
                  * Removes a face in a mesh, without any check
                  *
                  * \param m Mesh
                  * \param fp Face to be removed
                  *
                  */
                static void remove_face_unsafe (MeshType &m, FacePointer fp)
                {

                    vector<HEdgePointer> hps = getHEdges(fp);

                    int size = hps.size();

                    for( int i = 0; i< size; i++ )
                    {
                        if( hps[i]->HOp()->HFp() )
                        {
                            hps[i]->HFp() = NULL;

                            if( !hps[(i+size-1)%size]->HOp()->HFp() )
                            {
                                // HP
                                hps[i]->HPp() = hps[(i+size-1)%size]->HOp()->HPp();
                                hps[(i+size-1)%size]->HOp()->HPp()->HNp() = hps[i];
                            }

                            if( !hps[(i+1)%size]->HOp()->HFp() )
                            {
                                // HN
                                hps[i]->HNp() = hps[(i+1)%size]->HOp()->HNp();
                                hps[(i+1)%size]->HOp()->HNp()->HPp() = hps[i];
                            }
                        }
                        else
                        {
                            Allocator<MeshType>::DeleteHEdge( m, *hps[i] );
                            Allocator<MeshType>::DeleteHEdge( m, *(hps[i]->HOp()) );
                            Allocator<MeshType>::DeleteEdge( m, *(hps[i]->HEp()) );


                            if( !hps[(i+size-1)%size]->HOp()->HFp() )
                            {
                                hps[i]->HOp()->HNp()->HPp() = hps[(i+size-1)%size]->HOp()->HPp();
                                hps[(i+size-1)%size]->HOp()->HPp()->HNp() = hps[i]->HOp()->HNp();
                            }

                        }

                    }

                    for( int i = 0; i< size; i++ )
                    {
                        if( hps[i]->HVp()->VHp()->IsD() )
                        {
                            if( !hps[i]->IsD() )
                                hps[i]->HVp()->VHp() = hps[i];

                            else if( !hps[(i+size-1)%size]->IsD() )
                                hps[i]->HVp()->VHp() = hps[(i+size-1)%size]->HOp();

                            else //search for a hedge (hedge can be found only if the vertex is non-manifold)
                            {
                                bool manifold = true;

                                Pos<MeshType> p(hps[i]->HVp()->VHp(), true);

                                p.HE()->SetV();

                                p.FlipE();
                                p.FlipF();

                                while( !p.HE()->IsV() )
                                {
                                    if( !p.HE()->IsD() )
                                    {
                                        manifold = false;
                                        hps[i]->HVp()->VHp() = p.HE();
                                        break;
                                    }

                                    p.FlipE();
                                    p.FlipF();
                                }

                                p.HE()->ClearV();

                                if(manifold)
                                    hps[i]->HVp()->VHp() = NULL;

                            }
                        }

                    }
                    Allocator<MeshType>::DeleteFace(m,*fp);

                }

                /*!
                  * Checks if the next hedge can be inserted into hps.
                  * If true, inserts the hedge into hps. If false, inserts NULL.
                  *
                  * \param vps Vector of vertices (in ccw order) that will belong to the new face
                  * \param hps Vector of hedges already checked
                  *
                  * \retval true if hedge can be inserted
                  * \retval false otherwise
                  */
                static bool can_add_hedge( vector<VertexPointer> &vps, vector<HEdgePointer> &hps )
                {

                    unsigned int i = hps.size();

                    assert( i < vps.size() );

                    HEdgePointer he = vps[i]->VHp();

                    if(!he) //vertex is detached
                    {
                        hps.push_back(NULL);
                        return true;
                    }
                    else
                    {
                        bool disconnected = false;

                        bool hasEdge = false;

                        unsigned int size = vps.size();

                        Pos<MeshType> p(he, false);

                        he->SetV();

                        while(p.V() != vps[(i+1)%size])
                        {
                            if(!hasEdge)
                                hasEdge= ( find( vps.begin(), vps.end(), p.V()) != (vps.end() ) );

                            p.FlipV();

                            p.FlipE();
                            p.FlipF();

                            p.FlipV();

                            if(p.HE()->IsV())
                            {
                                disconnected = true;
                                break;
                            }

                        }

                        he->ClearV();

                        if(disconnected) // edge does not exist
                        {
                            hps.push_back(NULL);

                            // if hasEdge is false after inserting the face there will be a non-manifold vertex
                            return hasEdge;
                        }

                        else //edge already existing
                        {
                            // try to insert consecutve hedges if they will belong to the new face
                            while( (p.V() ==  vps[(i+1)%size])  && (i < size) )
                            {
                                hps.push_back( p.HE() );

                                if(p.HE()->HFp() != NULL)
                                    return false;

                                i++;
                                p.FlipE();
                                p.FlipV();
                            }
                            return true;
                        }
                    }
                }

                /*!
                  * Checks if a face can be removed
                  *
                  * \param fp Face to check
                  *
                  * \retval true if the face can be removed
                  * \retval false otherwise
                  */
                static bool can_remove_face(FacePointer fp)
                {

                    Pos<MeshType> p(fp->FHp(), true);

                    do
                    {
                        vector<FacePointer> incident_faces = get_incident_faces( p.V() );

                        unsigned int size = incident_faces.size();

                        if(size > 2)
                        {
                            for(unsigned int i = 0; i < size; i++)
                            {
                                if(incident_faces[i] == NULL)
                                    if(incident_faces[(i+1)%size] != fp && incident_faces[((i+size)-1)%size] != fp )
                                        return false;
                            }
                        }

                        p.FlipV();
                        p.FlipE();

                    }while( p.HE() != fp->FHp() );

                    return true;
                }

                /*!
                  * Gets all vertices incident to a face
                  *
                  * \param fp Face
                  * \param starting_he A hedge in the face from which to start
                  *
                  * \return Vector containing the incident vertices
                  */
                static vector<VertexPointer> getVertices(FacePointer fp, HEdgePointer starting_he = NULL)
                {
                    if(starting_he)
                        assert( starting_he->HFp() == fp );
                    else
                        starting_he = fp->FHp();

                    Pos<MeshType> p( starting_he, true );

                    vector<VertexPointer> ret;

                    if( fp )
                    {
                        do
                        {
                            ret.push_back( p.V() );

                            p.FlipV();
                            p.FlipE();

                        }while(p.HE() != starting_he);
                    }
                    return ret;

                }

                /*!
                  * Gets all edges incident to a face
                  *
                  * \param fp Face
                  * \param starting_he A hedge in the face from which to start
                  *
                  * \return Vector containing the incident edges
                  */
                static vector<HEdgePointer> getHEdges(FacePointer fp, HEdgePointer starting_he = NULL)
                {
                    if(starting_he)
                        assert( starting_he->HFp() == fp );
                    else
                        starting_he = fp->FHp();

                    Pos<MeshType> p( starting_he, true );

                    vector<HEdgePointer> ret;

                    if( fp )
                    {
                        do
                        {
                            ret.push_back( p.HE() );

                            p.FlipV();
                            p.FlipE();

                        }while(p.HE() != starting_he);
                    }
                    return ret;

                }

                /*!
                  * Gets all faces incident to a vertex
                  *
                  * \param fp Vertex
                  * \param starting_he A hedge from which to start
                  *
                  * \return Vector containing the incident faces
                  */
                static vector<FacePointer> get_incident_faces(VertexPointer vp, HEdgePointer starting_he = NULL)
                {
                    assert(vp);

                    if(starting_he)
                        assert( starting_he->HVp() == vp );
                    else
                        starting_he = vp->VHp();

                    Pos<MeshType> p( starting_he, true );

                    vector<FacePointer> ret;

                    do
                    {
                        ret.push_back( p.F() );

                        p.FlipE();
                        p.FlipF();

                    }while(p.HE() != starting_he);

                    return ret;

                }

                /*!
                  * Connects to a new vertex all hedges incident to a vertex
                  *
                  * \param old_vp the old vertex to be disconnected
                  * \param new_vp the new vertex to be connected
                  *
                  */
                static void change_vertex(VertexPointer old_vp, VertexPointer new_vp)
                {
                    assert(old_vp);
                    assert(new_vp);
                    assert(old_vp != new_vp);

                    Pos<MeshType> p(old_vp->VHp(),true);

                    p.HE()->SetV();

                    do
                    {
                        p.HE()->HVp() = new_vp;

                        p.FlipE();
                        p.FlipF();

                    }while( !p.HE()->IsV() );

                    p.HE()->ClearV();

                    if( !new_vp->VHp() )
                        new_vp->VHp() = old_vp->VHp();

                }

            };

    }
}

#endif // VCG_HEDGE_TOPOLOGY

