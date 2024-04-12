#include "convexhull.h"
#include <vector>

int main()
{
  std::vector<PointnD> tests = \
    {{3,{3,0,0}}, {3, {0,3,0}}, {3, {0,0,3}}, {3, {3,3,3}}};
  
  cout<<tests[0]<<endl;
  // nConvexHull C(tests);
  // std::cout<<"vertices (x, y, z, intensity)\n"; C.Print("vertice");
  // std::cout<<"faces\n";  C.Print("face");

  return 0;
}