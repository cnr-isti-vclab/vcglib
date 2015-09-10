#include <iostream>
#include <nanoply.hpp>

using namespace std;
using namespace nanoply;
struct Point3f
{
  float x,y,z;
};

int main()
{
  std::vector<Point3f> coordVec;
  std::vector<Point3f> normalVec;

  nanoply::Info info;
  nanoply::GetInfo("sphere13k.ply",info);
//  Open("pippo.ply",std::make_tuple(DataDescriptor<std::tuple<float,float,float> >(NNP_PXYZ,&*coordVec.begin()),
//                                   DataDescriptor<std::tuple<float,float,float> >(NNP_NXYZ,&*normalVec.begin())),
//        info );

  return 0;
}

