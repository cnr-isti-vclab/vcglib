#include <iostream>
#include "crude.h"
using namespace nxs;
using namespace vcg;
using namespace std;

int main(int argc, char *argv[]) {
  if(argc != 2) {
    cerr << "Uso: " << argv[0] << " <side>\n";
    return -1;
  }

  unsigned int side = atoi(argv[1]);
  if(side > 1000 || side == 0) {
    cerr << "Invalid side: " << argv[1] << endl;
    return -1;
  }

  Crude crude;
  if(!crude.Create("square")) {
    cerr << "Could not create square" << endl;
    return -1;
  }
  
  int half = side/2;
  crude.Resize(side * side, (side-1) * (side-1) * 2);
  for(unsigned int x = 0; x < side; x++) 
    for(unsigned int y = 0; y < side; y++) {
      //      Point3f p(x*x*x/((float)side), 
      //		y*y*y/((float)side), x*y/((float)side));
      Point3f p(x, y, sqrt((float)x*x + y*y));
      crude.SetVertex(x + side * y, p);
      crude.GetBox().Add(p);
    }

  for(unsigned int x = 0; x < side-1; x++) 
    for(unsigned int y = 0; y < side-1; y++) {
      unsigned int pos = x + side*y;
      Crude::Face face(pos, pos + 1, pos + side);
      crude.SetFace(0 + 2*x + (side-1)*y*2, face);

      face = Crude::Face(pos + 1, pos + 1 + side, pos +side);
      crude.SetFace(1 + 2*x + (side-1)*y*2,  face);
    }
  
  crude.Close();
  return 0;
}
