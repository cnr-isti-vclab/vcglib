#ifndef VCG_PIVOT_H
#define VCG_PIVOT_H

#include <vector>
#include <list>
#include <wrap/callback.h>

#include "vcg/space/index/grid_static_ptr.h"
#include "vcg/complex/trimesh/closest.h"

namespace vcg {
namespace tri {



    
template <class MESH>
class Pivot {
  public:
    typedef GridStaticPtr<typename MESH::VertexType, typename MESH::ScalarType > StaticGrid;
    typedef typename MESH::VertexType CVertex;
    typedef typename MESH::FaceType CFace;
    typedef typename MESH::ScalarType ScalarType;
    typedef typename CVertex::CoordType Point3x;

    ScalarType radius;  //default 1 (not meaningful
    ScalarType mindist; //minimum distance between points in the mesh (% of radius)   
    ScalarType crease;  // -0.5 
    bool normals;       //default false
    Box3<ScalarType> box;
    
    MESH &mesh;    
    StaticGrid grid;
    
    template <class Coord> class Edge { 
      public:       
      int v0, v1, v2;   //v0, v1 represent the edge, v2 the other vertex in the face
                        //this edge belongs to
      int face;        //corresponding face
      Coord center;  //center of the sphere touching the face
      int count;   //test delay touch edges. 
      
      bool active; //keep tracks of wether it is in front or in deads
      float angle;
      int candidate;
      Coord newcenter;
      
      //the loops in the front are mantained as a double linked list
      typename std::list<Edge>::iterator next;            
      typename std::list<Edge>::iterator previous;
    
      Edge() {}
      Edge(int _v0, int _v1, int _v2, int _face, Point3f &_center): 
               v0(_v0), v1(_v1), v2(_v2), 
               face(_face), center(_center), count(-1), active(true) {
                 assert(v0 != v1 && v1 != v2 && v0 != v2);
               }
    };    
    typedef Edge<Point3x> Edgex;
    
    /* front edges of the mesh: 
       to expand the front we get the first edge
       if an edge cannot create a new triangle it is marked dead and moved
       to the end of the list
       the new edges are inserted to the back (before dead_begin) */
       
    std::list<Edgex> front;   
    std::list<Edgex> deads;
    std::vector<int> nb; //number of fronts a vertex is into,
                         //this is used for the Visited and Border flags
                         //but adding topology may not be needed anymode
    int last_seed;

      
Pivot(MESH &_mesh, ScalarType _radius, ScalarType _mindist = 0.1, ScalarType _crease = -0.5): 
       mesh(_mesh), radius(_radius), mindist(_mindist), crease(_crease), normals(true), last_seed(0) {
    
      //Compute bounding box. (this may be passed as a parameter?
      for(int i = 0; i < mesh.vert.size(); i++)
        box.Add(mesh.vert[i].P());

      //estimate radius if not provided     
      if(radius <= 0.0f)
        radius = sqrt((box.Diag()*box.Diag())/mesh.vn);
      
      /* we need to enlarge the grid to allow queries from little outside of the box 
         Someone is a bit lazy... */        
      box.Offset(4*radius);
      grid.Set(mesh.vert.begin(), mesh.vert.end(), box);
      nb.clear();
      nb.resize(mesh.vert.size(), 0);

      if(mesh.face.size()) {
        //init border from mesh
        Point3x center;
        CVertex *start = &*mesh.vert.begin();
        for(int i = 0; i < mesh.face.size(); i++) {
          CFace &face = mesh.face[i];          
          for(int k = 0; k < 3; k++) {
            if(!face.V(k)->IsB()) face.V(k)->SetV();
            cluster(face.V(k) - start);
            if(face.IsB(k)) {
              //compute center:
              findSphere(face.P(k), face.P((k+1)%3), face.P((k+2)%3), center);
              newEdge(Edgex(face.V((k)%3) -start, face.V((k+1)%3) - start, face.V((k+2)%3) - start,
                            i, center));
            }
          }
        }
        for(typename std::list<Edgex>::iterator s = front.begin(); s != front.end(); s++) {
          (*s).previous = front.end();
          (*s).next = front.end();
            printf("%d %d\n", (*s).v0, (*s).v1);
        }
        //now create loops:
        for(typename std::list<Edgex>::iterator s = front.begin(); s != front.end(); s++) {
          for(typename std::list<Edgex>::iterator j = front.begin(); j != front.end(); j++) {
            if(s == j) continue;
            if((*s).v1 != (*j).v0) continue;
            if((*j).previous != front.end()) continue;
            (*s).next = j;
            (*j).previous = s;

          }
        }
/*        for(typename std::list<Edgex>::iterator s = front.begin(); s != front.end(); s++) {
          assert((*s).next != front.end());
          assert((*s).previous != front.end());
        }   */
        for(int i = 0; i < mesh.face.size(); i++) {
          CFace &face = mesh.face[i];
          for(int k = 0; k < 3; k++) 
            face.V(k) = (CVertex *)(face.V(k) - start);
        }
        
      } else {
        for(int i = 0; i < mesh.vert.size(); i++) {
           mesh.vert[i].ClearFlags();
        }
        
      }      
      srand(time(NULL));
    }
    /* return false if you want to stop.\n */
void buildMesh(CallBackPos *call = NULL, int interval = 512) {
     bool done = false;
     float estimated_faces = mesh.vn*2;     
     while(!done) {
       //estimating progress
       float vdeleted = mesh.vert.size() - mesh.vn;
       float vused = mesh.face.size()/2.0f;       
       float progress = 100*(vdeleted + vused)/mesh.vert.size();
       if(progress > 99) progress = 99;
       if(call) call((int)progress, "Pivoting");
       for(int i = 0; i < interval; i++) {
         if(addFace() == -1) {       
           done = true;
           break;
         }
       }       
     }
     if(call) call(0, "Reindexing");
     for(int i = 0; i < mesh.face.size(); i++) {
       CFace &face = mesh.face[i];
       for(int k = 0; k < 3; k++)
         face.V(k) = &(mesh.vert[(int)face.V(k)]);
     }     

   }
    
    /* select a vertex at random, a small group of nearby vertices and looks
       for a sphere that touches 3 and contains none.
       Use the center of the box to get a sphere inside (or outside) the model 
       You may be unlucky... */
       
bool seed(bool outside = true, int start = -1) {   
           
      //pick a random point (well...)
      if(start == -1) start = 0;//rand()%mesh.vert.size();
      
      //get a sphere of neighbours
      std::vector<int> targets;
      std::vector<ScalarType> dists;
      int n = getInSphere(mesh.vert[start].P(), 2*radius, targets, dists);
      if(n < 3) { 
         mesh.vert[start].SetD();
        //bad luck. we should call seed again (assuming random pick) up to
        //some maximum tries. im lazy.
        return false;
      }
      int v0, v1, v2;
      bool found = false;
      //find a triplet that does not contains any other point
      Point3x center;
      for(int i = 0; i < n; i++) {
        v0 = targets[i];
        CVertex &vv0 = mesh.vert[v0];
        if(vv0.IsD() || vv0.IsB() || vv0.IsV()) continue;
        Point3x &p0 = mesh.vert[v0].P();
        Point3x out = (p0 - box.Center());
        if(!outside) out = -out;
    
        for(int k = i+1; k < n; k++) {
          v1 = targets[k];            
          CVertex &vv1 = mesh.vert[v1];
          if(vv1.IsD() || vv1.IsB() || vv1.IsV()) continue;
          Point3x &p1 = mesh.vert[v1].P();            
          if((p1 - p0).Norm() < mindist*radius) continue;
    
          for(int j = k+1; j < n; j++) {
            v2 = targets[j];
            CVertex &vv2 = mesh.vert[v2];
            if(vv2.IsD() || vv2.IsB() || vv2.IsV()) continue;
            Point3x &p2 = mesh.vert[v2].P();
            if((p2 - p0).Norm() < mindist*radius) continue;
            if((p2 - p1).Norm() < mindist*radius) continue;
            Point3x normal = (p1 - p0)^(p2 - p0);
            if(!normals) {
              //check normal pointing inside
              if(normal * out < 0) continue;
            } else {
              if(normal * vv0.N() < 0) continue;
              if(normal * vv1.N() < 0) continue;
              if(normal * vv2.N() < 0) continue;
            }
            if(!findSphere(p0, p1, p2, center)) continue;
            
            bool failed = false;
            //check no other point inside
            for(int t = 0; t < n; t++) {
              Point3x &p = mesh.vert[targets[t]].P();
              if((center - p).Norm() <= radius) {
                failed = true;
                break;
              }
            }
            //check on the other side there are not a surface
            Point3x recenter;
            if(!findSphere(p0, p2, p1, recenter)) continue;           
            for(int t = 0; t < n; t++) {
              CVertex &v = mesh.vert[targets[t]];
              Point3x &p = v.P();
              if((center - p).Norm() <= radius && (v.IsV() || v.IsB())) {
                failed = true;
                break;
              }
            }
            
            if(failed) continue;  
            found = true;
            i = k = j = n;
          }
        }
      }
      
      if(!found)  //see bad luck above
        return false;
      
      assert(!front.size());
      //TODO: should i check the edgex too?
      addFace(v0, v1, v2);
      
      //create the border of the first face  
      typename std::list<Edgex>::iterator e = front.end();
      typename std::list<Edgex>::iterator last;
      for(int i = 0; i < 3; i++) {
        int v0 = (int)(mesh.face.back().V0(i));
        int v1 = (int)(mesh.face.back().V1(i));    
        int v2 = (int)(mesh.face.back().V2(i));    
        nb[v0] = 1;
        assert(!mesh.vert[v0].IsB());
        mesh.vert[v0].SetB();
        Edgex edge(v0, v1, v2, 0, center);
        edge.previous = e;
        e = front.insert(front.begin(), edge);
        if(i == 0) last = e;
        (*edge.previous).next = e;
        
        cluster(v0);
      } 
    
      //connect last and first
      (*e).next = last;
      (*last).previous = e;
      return true;
    }
      
         
    /* expand the front adding 1 face. Return false on failure (id when
       all edges are dead  returns:
       1: added a face
       0: added nothing
       -1: finished         */
       
int addFace() {

      //We try to seed again
      if(!mesh.face.size()) {
        for(int i = 0; i < 100; i++) 
          if(seed()) return 1;
        return -1;
      }
      
      if(!front.size()) {
        //maybe there are unconnected parts of the mesh:
        //find a non D, V, B point and try to seed if failed D it.
        while(last_seed < mesh.vert.size()) {
          ++last_seed;
          CVertex &v = mesh.vert[last_seed-1];
          if(v.IsD() || v.IsV() || v.IsB()) continue;
          
          printf("seeding new: %i\n", last_seed-1);
          if(seed(true, last_seed-1)) return 1;
          
          v.SetD();
          --mesh.vn;
          return 0;
        }
        printf("done\n");
        return -1;
      }
      if(last_seed > 1) printf("frontsixe: %d\n", front.size());
      typename std::list<Edgex>::iterator ei = front.begin();
      Edgex &e = *ei;
      Edgex &previous = *e.previous;           
      Edgex &next = *e.next;        
 
      int v0 = e.v0, v1 = e.v1;
      
      assert(nb[v0] < 10 && nb[v1] < 10);
      int v2;
      Point3x center;
      bool success = pivot(e);
      v2 = e.candidate;
      center = e.newcenter;

      //if no pivoting or we are trying to connect to the inside of the mesh.
      if(!success || mesh.vert[v2].IsV()) { 
        killEdge(ei);
        return 0;
      } 
    
      //does v2 belongs to a front? (and which?)
      typename std::list<Edgex>::iterator touch = touches(ei);
    
      assert(v2 != v0 && v2 != v1);  
    
      int fn = mesh.face.size();
      if(touch != front.end()) {       

        //check for orientation and manifoldness    
        if(!checkEdge(v0, v2) || !checkEdge(v2, v1)) {                      
          killEdge(ei);
          return 0;
        }
        
        if(v2 == previous.v0) {   
               
          /*touching previous edge  (we reuse previous)        
                                    next
             ------->v2 -----> v1------>
                      \       /
                       \     /
               previous \   / e
                         \ /
                          v0           */
          
          detach(v0);
        
          typename std::list<Edgex>::iterator up = newEdge(Edgex(v2, v1, v0, fn, center));
          (*up).previous = previous.previous;
          (*up).next = e.next;
          (*previous.previous).next = up;
          next.previous = up;
          erase(e.previous);
          erase(ei);
          trovamiunnome(up);
          
           
        } else if(v2 == next.v1) {
        
               
        /*touching next edge  (we reuse next)        
          previous
             ------->v0 -----> v2------>
                      \       /
                       \     /
                        \   / next
                         \ /
                          v1           */      
    
          detach(v1);

          typename std::list<Edgex>::iterator up = newEdge(Edgex(v0, v2, v1, fn, center));
          (*up).previous = e.previous;
          (*up).next = (*e.next).next;
          previous.next = up;
          (*next.next).previous = up;
          erase(e.next);
          erase(ei);
          trovamiunnome(up);

              
        } else {     
    
        //touching some loop: split (or merge it is local does not matter.
        //like this 
        /*                 
                    left        right
                  <--------v2-<------
                          /|\
                         /   \
                     up /     \ down
                       /       \
                      /         V
                 ----v0 - - - > v1---------
                          e                         */           
          typename std::list<Edgex>::iterator left = touch;
          typename std::list<Edgex>::iterator right = (*touch).previous;      
          typename std::list<Edgex>::iterator up = ei;
          
          //this would be a really bad join
          if(v1 == (*right).v0 || v0 == (*left).v1) {
            killEdge(ei);
            return 0;
          }
          
          nb[v2]++;    
                      
          typename std::list<Edgex>::iterator down = newEdge(Edgex(v2, v1, v0, fn, center));      
    
          (*right).next = down;
          (*down).previous = right;
    
          (*down).next = e.next;
          next.previous = down;      
    
          (*left).previous = up;
          (*up).next = left;
          
          (*up).v2 = v1;      
          (*up).v1 = v2;
          (*up).face = fn;
          (*up).center = center;
          moveBack(ei);
        }                         
    
              
      
      } else {
    
        /*  adding a new vertex
                 
                           v2
                          /|\
                         /   \
                     up /     \ down
                       /       \
                      /         V
                 ----v0 - - - > v1--------- */
        assert(!mesh.vert[v2].IsB()); //fatal error! a new point is already a border?
        
        //clustering points aroundf v2                 
        cluster(v2);    
            
        nb[v2]++;                 
        mesh.vert[v2].SetB();
        typename std::list<Edgex>::iterator down = newEdge(Edgex(v2, v1, v0, fn, center));
        (*down).previous = ei;
        (*down).next = e.next;
        next.previous = down;
        
        e.v2 = v1;    
        e.v1 = v2;
        e.face = fn;
        e.center = center;
        e.next = down; 
        moveBack(ei);
      }
      addFace(v0, v2, v1);
      return 1;
    }        
    

    
    /* return new vertex and the center of the new sphere pivoting from edge
       if the vertex belongs to another edge, touch points to it. */
    bool pivot(Edgex &edge) {
        Point3x v0 = mesh.vert[edge.v0].P();
        Point3x v1 = mesh.vert[edge.v1].P();  
        Point3x v2 = mesh.vert[edge.v2].P();  
        /* TODO why using the face normals everything goes wrong? should be
           exactly the same................................................
           Check if the e.face is correct.
           Point3x &normal = mesh.face[edge.face].N();
        */
    
        Point3x normal = ((v1 - v0)^(v2 - v0)).Normalize();        
        Point3x middle = (v0 + v1)/2;    
        Point3x start_pivot = edge.center - middle;          
        Point3x axis = (v1 - v0);
        
        ScalarType axis_len = axis.SquaredNorm();
        if(axis_len > 4*radius*radius) return false;
        axis.Normalize();
        
        // r is the radius of the thorus of all possible spheres passing throug v0 and v1
        ScalarType r = sqrt(radius*radius - axis_len/4);
        
        std::vector<int> targets;
        std::vector<ScalarType> dists;    
        getInSphere(middle, r + radius, targets, dists);
    
        if(targets.size() == 0) return false; //this really would be strange but one never knows.
    
        edge.candidate = -1;
        ScalarType minangle = 0;
        Point3x center;  //to be computed for each sample
        for(int i = 0; i < targets.size(); i++) {      
          int id = targets[i];
          
          if(id == edge.v0 || id == edge.v1 || id == edge.v2) continue;
    
          if(mesh.vert[id].IsD()) 
            continue;                    

          Point3x p = mesh.vert[id].P();
    
                                
          /* Find the sphere through v0, p, v1 (store center on end_pivot */
          if(!findSphere(v0, p, v1, center)) {
            continue;      
          }
          
          /* Angle between old center and new center */
          ScalarType alpha = angle(start_pivot, center - middle, axis);
    
          /* adding a small bias to already chosen vertices.
             doesn't solve numerical problems, but helps. */
//          if(mesh.vert[id].IsB()) alpha -= 0.001;
          
          /* Sometimes alpha might be little less then M_PI while it should be 0,
             by numerical errors: happens for example pivoting 
             on the diagonal of a square. */
          
          if(alpha > 2*M_PI - 0.8) {               
            // Angle between old center and new *point* 
            //TODO is this really overshooting? shouldbe enough to alpha -= 2*M_PI
            Point3x proj = p - axis * (axis * p - axis * middle);
            ScalarType beta = angle(start_pivot, proj - middle, axis);
          
            if(alpha > beta) alpha -= 2*M_PI; 
          }
         //scale alpha by distance:
          if(edge.candidate == -1 || alpha < edge.angle) {
            edge.candidate = id; 
            edge.angle = alpha;
            edge.newcenter = center;
          }
        }
        
        if(edge.candidate == -1) return false;
        Point3x n = ((mesh.vert[edge.candidate].P() - v0)^(v1 - v0)).Normalize();
        //found no point suitable.
        if(normal * mesh.vert[edge.candidate].N() < 0 ||
           n * normal < crease ||
           nb[edge.candidate] >= 2) {
          return false;
        }
           
        assert(edge.candidate != edge.v0 && edge.candidate != edge.v1);
        return true;
    }         
    
  private:
     //front management:
     //Add a new edge to the back of the queue
     typename std::list<Edgex>::iterator newEdge(Edgex e) {  
       e.active = true;                
       return front.insert(front.end(), e);
     }     
     //move an Edge among the dead ones
     void killEdge(typename std::list<Edgex>::iterator e) {
       (*e).active = false;
       deads.splice(deads.end(), front, e);
     }

     void erase(typename std::list<Edgex>::iterator e) {
       if((*e).active) front.erase(e);
       else deads.erase(e);
     }
     //move an Edge to the back of the queue
     void moveBack(typename std::list<Edgex>::iterator e) {
       front.splice(front.end(), front, e);          
     }
     
     void moveFront(typename std::list<Edgex>::iterator e) {
       front.splice(front.begin(), front, e);
     }
    bool checkEdge(int v0, int v1) {
      int tot = 0;
      //HACK to speed up things until i can use a seach structure
      int i = mesh.face.size() - 4*(front.size());
      if(front.size() < 100) i = mesh.face.size() - 100;
//      i = 0;
      if(i < 0) i = 0;
      for(; i < mesh.face.size(); i++) { 
        CFace &f = mesh.face[i];
        for(int k = 0; k < 3; k++) {
          if(v1== (int)f.V(k) && v0 == (int)f.V((k+1)%3)) ++tot;
          else if(v0 == (int)f.V(k) && v1 == (int)f.V((k+1)%3)) { //orientation non constistent
             return false;
          }
        }
        if(tot >= 2) { //non manifold
          return false;
        }
      }
      return true;
    }        
    
               
    void cluster(int v) {
      /* clean up too close points */
        std::vector<int> targets;
        std::vector<ScalarType> dists;    
        getInSphere(mesh.vert[v].P(), mindist*radius, targets, dists);
        
        for(int i = 0; i < targets.size(); i++) {
          int id = targets[i];
          if(id == v) continue;
    
          CVertex &v = mesh.vert[id];
          if(v.IsD() || v.IsV() || v.IsB()) continue;
          v.SetD();
          --mesh.vn;
        }
            
    }
    
    bool trovamiunnome(typename std::list<Edgex>::iterator e) {
      return glue((*e).previous, e) || glue(e, (*e).next);
    }
    
    //glue toghether a and b (where a.next = b
    bool glue(typename std::list<Edgex>::iterator a, typename std::list<Edgex>::iterator b) {
      if((*a).v0 != (*b).v1) return false; 
      
      typename std::list<Edgex>::iterator previous = (*a).previous;
      typename std::list<Edgex>::iterator next = (*b).next;
      (*previous).next = next;
      (*next).previous = previous;
      detach((*a).v1);
      detach((*a).v0); 
      front.erase(a);
      front.erase(b);  
      return true;
    }
    
    void detach(int v) {
      assert(nb[v] > 0);
      if(--nb[v] == 0) {
        mesh.vert[v].SetV();
        mesh.vert[v].ClearB();      
      }
    }

          
    /* compute angle from p to q, using axis for orientation */
    ScalarType angle(Point3x p, Point3x q, Point3x &axis) {
      p.Normalize();
      q.Normalize();
      Point3x vec = p^q;
      ScalarType angle = acos(p*q);
      if(vec*axis < 0) angle = -angle;
      if(angle < 0) angle += 2*M_PI;
      return angle;
    }          
    /* add a new face. compute normals. */
    void addFace(int a, int b, int c) {
      CFace face;
      face.V(0) = (CVertex *)a;
      face.V(1) = (CVertex *)b;
      face.V(2) = (CVertex *)c;
      Point3x &p0 = mesh.vert[a].P();
      Point3x &p1 = mesh.vert[b].P();
      Point3x &p2 = mesh.vert[c].P();
      face.N() = ((p1 - p0)^(p2 - p0)).Normalize();
      
      mesh.face.push_back(face);
      mesh.fn++;
    }         
              
                             
    /* intersects segment [v0, v1] with the sphere of radius radius. */
    bool intersect(int v0, int v1, Point3x &center) {
      Point3x m =  mesh.vert[v1].P() -  mesh.vert[v0].P();
      ScalarType t = m*(center - mesh.vert[v0].P());
      if(t < 0) return false;
      if(t > m*m) return false;
      return true;
    }         


    ScalarType distance(int v0, int v1, Point3x &center) {
      Point3x m =  mesh.vert[v1].P() -  mesh.vert[v0].P();
      ScalarType t = m*(center - mesh.vert[v0].P())/(m*m);
      Point3x p = mesh.vert[v0].P() + m*t;
      return (p - center).Norm();
    }             


    /* return all point in a given ball, notice as i want the index
       of the vertices not the pointers... this may change in future */
    unsigned int getInSphere(Point3x &p, ScalarType distance, 
                             std::vector<int> &results,
                             std::vector<ScalarType> &dists) {
      std::vector<CVertex *> ptr;
      std::vector<Point3x> points;
      int n = trimesh::GetInSphereVertex(mesh, grid, p, distance, ptr, dists, points);
      for(int i = 0; i < ptr.size(); i++) 
        results.push_back(ptr[i] - &(mesh.vert[0]));
      return n;
    }                                                


    
    /* returns the sphere touching p0, p1, p2 of radius r such that
       the normal of the face points toward the center of the sphere */
    bool findSphere(Point3x &p0, Point3x &p1, Point3x &p2, Point3x &center) {
      Point3x q1 = p1 - p0;
      Point3x q2 = p2 - p0;  
    
      Point3x up = q1^q2;
      ScalarType uplen = up.Norm();
    
      //the three points are aligned
      if(uplen < 0.001*q1.Norm()*q2.Norm()) return false;
      up /= uplen;
      
    
      ScalarType a11 = q1*q1;
      ScalarType a12 = q1*q2;
      ScalarType a22 = q2*q2;
    
      ScalarType m = 4*(a11*a22 - a12*a12);
      ScalarType l1 = 2*(a11*a22 - a22*a12)/m;
      ScalarType l2 = 2*(a11*a22 - a12*a11)/m;
    
      center = q1*l1 + q2*l2;
      ScalarType circle_r = center.Norm();
      if(circle_r > radius) return false; //need too big a sphere
    
      ScalarType height = sqrt(radius*radius - circle_r*circle_r);
      center += p0 + up*height;
    
      return true;
    }         

    typename std::list<Edgex>::iterator touches(typename std::list<Edgex>::iterator e) {
      //TODO what happens when it touches more than one front?
      //might still work.
      int v = (*e).candidate;
      typename std::list<Edgex>::iterator touch = front.end();
      if(mesh.vert[v].IsB()) {
        //test nearby Edges: it is faster
        typename std::list<Edgex>::iterator p = e;
        p = (*e).previous;
        if(v == (*p).v0) return p;
        e = (*e).next;
        if(v == (*e).v0) return e;
    
        p = (*p).previous;
        if(v == (*p).v0) return p;
        e = (*e).next;
        if(v == (*e).v0) return e;
    
        //test all. sigh.
    
        for(typename std::list<Edgex>::iterator k = front.begin(); k != front.end(); k++) {
          if(v == (*k).v0) { 
            touch = k;
            break;
          }
        }
        for(typename std::list<Edgex>::iterator k = deads.begin(); k != deads.end(); k++) {
          if(v == (*k).v0) { 
            touch = k;
            break;
          }
        }
        assert(touch != front.end());
      }
    
      return touch;
    }                              
    
    
 public:
        
 };
        

}//namespace
}//namespace
/*  CODE FOR PIVOTING IN  A TOUCH SITUATION not used now.

    //if touch we want to check the ball could really pivot around that point
    if(touch != front.end() && touch != (*Edge.next).next && touch != Hinge.previous) {
      Point3x &hinge = mesh.vert[min].P();      
      Point3x target = (*touch).center - hinge;
      float d = (target * start_pivot)/(target.Norm()*start_pivot.Norm());
      if(d < -0.8) {
        return false;
      }
      
      

      if(d < 0.5) { //they are far enough so test .
        Point3x naxis = (target ^ start_pivot).Normalize();
        Point3x d1 = naxis^start_pivot;
        Point3x d2 = target^naxis;          
        
        for(int i = 0; i < targets.size(); i++) {
          int id = targets[i];
          if(id == Hinge.v0 || id == Hinge.v1 || id == Hinge.v2 || id == min) continue;
          if(mesh.vert[id].IsD()) {
            continue;
          }

          Point3x intruder = mesh.vert[targets[i]].P() - hinge;
          float h = intruder*naxis;
          if(fabs(h) > radius) continue;
          intruder -= naxis*h;
          assert(fabs(intruder *naxis) < 0.01);
          float off = radius - intruder.Norm(); //(distance from the center ring of the thorus
          if(h*h + off*off > radius*radius) continue; //outside of thorus
          if(d1*intruder < 0 || d2*intruder < 0) continue; //ouside of sector
          cout << "could not pivot while touching;\n";
          return false;
        }
        
      }          
    }*/


#endif
