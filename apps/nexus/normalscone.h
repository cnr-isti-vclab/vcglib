#ifndef NXS_NORMALS_CONE_H
#define NXS_NORMALS_CONE_H

#include <vector>

#include <vcg/space/point3.h>
#include <vcg/space/sphere3.h>

namespace nxs {

  
  //anchored normal cone.
  class ANCone3f {
  public:
    ANCone3f();

    bool Frontface(const vcg::Point3f &viewPoint);
    bool Backface(const vcg::Point3f &viewPoint);

    //threshold [0, 1]: percentage of normals to left out
    void AddNormals(std::vector<vcg::Point3f> &normals, 
		    float threshold = 1.0f);
    void AddNormals(std::vector<vcg::Point3f> &normals, 
		    std::vector<float> &areas, 
		    float threshold = 1.0f);
    void AddAnchors(std::vector<vcg::Point3f> &anchors);
  protected:
    vcg::Point3f scaledNormal;
    vcg::Point3f frontAnchor;
    vcg::Point3f backAnchor;
    
    friend class NCone3s;
  };



  //this is not anchored... need a bounding sphere to report something
  class NCone3s {
  public:
    void Import(const ANCone3f &c);
    
    bool Backface(const vcg::Sphere3f &sphere, 
		  const vcg::Point3f &view) const;
    bool Frontface(const vcg::Sphere3f &sphere, 
		   const vcg::Point3f &view) const;    

    //encode normal and sin(alpha) in n[3]
    short n[4];
  };

}

#endif
