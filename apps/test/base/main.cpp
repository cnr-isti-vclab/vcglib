#include <stdio.h>
#include<vcg/math/base.h>
#include<vcg/space/point3.h>
#include<vcg/space/point4.h>
#include<vcg/space/color4.h>

using namespace vcg;

int main(int argc, char *argv[])
{
  printf("Hello Library!\n");
  Point3f pp0(0,1,2);
  Point3f pp1(2,1,0);
  Point3f pp2=pp1+pp0;
  Point3i ppi=Point3i::Construct(pp1+pp0);
  
  Point4i size(0,0,1,1);
  
  Color4b cb(Color4b::LightBlue);
  Color4f cf(Color4f::LightBlue);
  Color4b cbi; cbi.Import(cf);
  printf("ci %i %i %i %i\n",cbi.V(0),cbi.V(1),cbi.V(2),cbi.V(3));
  return -1;
}