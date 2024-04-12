#include<iostream>
#include<cmath>
#include<vector>
#include "utility2.h"
#include "convexhull.h"

size_t ConvexHull::Key2Edge(const Point3D& a, const Point3D& b) const
{
  point_hash ph;
  return ph(a) ^ ph(b);
}

// size_t nConvexHull::Key2Edge(const PointnD& a, const PointnD& b) const
// {
//   npoint_hash ph;
//   return ph(a) ^ ph(b);
// }

double detc(const vector<vector<double>> &matG)
{
    const std::size_t nrows{matG.size()}, ncols{matG[0].size()};
    double detG = 0.0;
    if (nrows != ncols)
    {
        std::cout << "Error when using det: matrix is not square.\n";
        return detG;
    }

    const std::size_t nSize{nrows};
    std::size_t i{0}, j{0}, k{0};

    // ******************** Step 1: row permutation (swap diagonal zeros) ********************
    vector<vector<double>> matLU;
    std::vector<std::size_t> permuteLU; // Permute vector
    bool changeSign{false};

    for (i = 0; i < nSize; ++i)
    {
        permuteLU.push_back(i);
    }

    for (j = 0; j < nSize; ++j)
    {
        double maxv{0.0};
        for (i = j; i < nSize; ++i)
        {
            const double currentv{std::abs(matG[permuteLU[i]][j])};
            if (currentv > maxv)
            {
                maxv = currentv;
                if (permuteLU[i] != permuteLU[j]) // swap rows
                {
                    changeSign = !changeSign;
                    const std::size_t tmp{permuteLU[j]};
                    permuteLU[j] = permuteLU[i];
                    permuteLU[i] = tmp;
                }
            }
        }
    }

    for (i = 0; i < nSize; ++i)
    {
        matLU.push_back(matG[permuteLU[i]]);
    }
    
    // ******************** Step 2: LU decomposition (save both L & U in matLU) ********************
    if (matLU[0][0] == 0.0)
    {
        return detG; // Singular matrix, det(G) = 0
    }

    for(int i=1; i<nSize; i++){
        double curv = matLU[i-1][i-1];
        for(int j=i; j<nSize; j++){
            double f = matLU[j][i-1] / curv;
            for(int k=i-1; k<nSize; k++){
                matLU[j][k] -= f * matLU[i-1][k];
            }
        }
        if(matLU[i][i] == 0) {
            int flag = 0;
            for(int j=i+1; j<nSize; j++){
                if(matLU[j][i] != 0.0) {
                    swap(matLU[j], matLU[i]);
                    changeSign ^= 1;
                    flag = 1;
                    break;
                }
            }
            if(flag == 0) {
                return detG;
            }
        }
        // for(int j=0; j<nSize; j++){
        //     for(int k=0; k<nSize; k++){
        //         cout<<matLU[j][k]<<" ";
        //     }
        //     cout<<endl;
        // }
        if(matLU[i][i] == 0.0) return detG;
    }

    // for (i = 1; i < nSize; ++i)
    // {
    //     matLU[i][0] /= matLU[0][0];
    // }

    // for (i = 1; i < nSize; ++i)
    // {
    //     for (j = i; j < nSize; ++j)
    //     {
    //         for (k = 0; k < i; ++k)
    //         {
    //             matLU[i][j] -= matLU[i][k] * matLU[k][j]; // Calculate U matrix
    //         }
    //     }
    //     if (matLU[i][i] == 0.0)
    //     {
    //         // for(int p=0; p<nSize; p++){
    //         //     for(int q=0; q<nSize; q++){
    //         //         cout<<matLU[p][q]<<" ";
    //         //     }
    //         //     cout<<endl;
    //         // }
    //         return detG; // Singular matrix, det(G) = 0
    //     }
    //     for (k = i + 1; k < nSize; ++k)
    //     {
    //         for (j = 0; j < i; ++j)
    //         {
    //             matLU[k][i] -= matLU[k][j] * matLU[j][i]; // Calculate L matrix
    //         }
    //         matLU[k][i] /= matLU[i][i];
    //     }
    // }

    detG = 1.0;
    if (changeSign)
    {
        detG = -1.0; // Change the sign of the determinant
    }
    for (i = 0; i < nSize; ++i)
    {
        detG *= matLU[i][i]; // det(G) = det(L) * det(U). For triangular matrices, det(L) = prod(diag(L)) = 1, det(U) = prod(diag(U)), so det(G) = prod(diag(U))
    }

    return detG;
}

vector<vector<double>> dyn_inv(vector<vector<double>> &mat, vector<vector<double>> &invMat, vector<double> u, int col = 0) {
    int nSize = mat.size();
    vector<vector<double>> matG(nSize, vector<double>(nSize));

    vector<double> h1(nSize);
    for(int i=0; i<nSize; i++){
        h1[i] = 0.0;
        for(int j=0; j<nSize; j++){
            h1[i] += invMat[i][j] * (u[j] - mat[j][col]);
        }
    }

    for(int i=0; i<nSize; i++){
        h1[i] /= (1.0 + h1[col]);
    }

    vector<vector<double>> iMatG = invMat;
    for(int i=0; i<nSize; i++){
        for(int j=0; j<nSize; j++) {
            iMatG[i][j] -= h1[i] * invMat[col][j];
        }
    }
    return iMatG;
}

double dyn_det(vector<vector<double>> &mat, vector<vector<double>> &invMat, double det, vector<double> u, int col = 0) {
    int nSize = mat.size();
    double val = 0.0;
    for(int i=0; i<nSize; i++){
        val += invMat[col][i] * (u[i] - mat[i][col]);
    }
    return (1.0 + val) * det;
}

double ConvexHull::VolumeSign(const Face& f, const Point3D& p, bool debug) const
{
  double vol;
  vector<vector<double>> mat = {{f.vertices[0].x, f.vertices[1].x, f.vertices[2].x, p.x},
                                {f.vertices[0].y, f.vertices[1].y, f.vertices[2].y, p.y},
                                {f.vertices[0].z, f.vertices[1].z, f.vertices[2].z, p.z},
                                {1, 1, 1, 1}};
  if(debug){
    for(int p=0; p<4; p++){
    for(int q=0; q<4; q++){
      cout<<mat[p][q]<<" ";
    }
    cout<<'\n';
  }
  }
  vol = detc(mat);
  if(vol == 0 || vol == -0) return 0;
  return vol ;

}

// double nConvexHull::VolumeSign(const nFace& f, const PointnD& p, bool debug) const
// {
//   double vol;
//   vector<vector<double>> mat = {};
//   for(auto v: f.vertices){
//     vector<i_float_t> row = v.X;
//     row.push_back(1.0);
//     mat.push_back(row);
//   }
//   vector<i_float_t> row = p.X;
//   row.push_back(1.0);
//   mat.push_back(row);

//   vol = detc(mat);
//   return vol;
// }

void ConvexHull::AddOneFace(const Point3D& a, const Point3D& b, 
    const Point3D& c, const Point3D& inner_pt)
{
  // Make sure face is CCW with face normal pointing outward
  this->faces.emplace_back(a, b, c);
  auto& new_face = this->faces.back();
  double D = this->VolumeSign(this->faces.back(), inner_pt, 0);
  vector<vector<double>> mat;
  if(D < 0) {
    new_face.Reverse();
    mat = {{c.x, b.x, a.x, inner_pt.x},
          {c.y, b.y, a.y, inner_pt.y},
          {c.z, b.z, a.z, inner_pt.z},
          {1, 1, 1, 1}};
  } 
  else{
    mat = {{a.x, b.x, c.x, inner_pt.x},
          {a.y, b.y, c.y, inner_pt.y},
          {a.z, b.z, c.z, inner_pt.z},
          {1, 1, 1, 1}};
  }
  new_face.originalMatrix = mat;
  new_face.inverseMatrix = pinv(mat);
  new_face.determinant = detc(mat);

  // Create edges and link them to face pointer
  auto create_edge = [&](const Point3D& p1, const Point3D& p2)
  {
    size_t key = this->Key2Edge(p1, p2);
    if(!this->map_edges.count(key)) 
    { 
      this->edges.emplace_back(p1, p2);
      this->map_edges.insert({key, &this->edges.back()});
    }
    this->map_edges[key]->LinkAdjFace(&new_face);
  };
  create_edge(a, b);
  create_edge(a, c);
  create_edge(b, c);
}

// void nConvexHull::AddOneFace(const vector<PointnD> pts, const PointnD& inner_pt)
// {
//   // Make sure face is CCW with face normal pointing outward
//   this->faces.emplace_back(pts);
//   auto& new_face = this->faces.back();
//   int sign = this->VolumeSign(this->faces.back(), inner_pt, 1);
//   if(sign <= 0) new_face.Reverse();

//   // Create edges and link them to face pointer
//   auto create_edge = [&](const PointnD& p1, const PointnD& p2)
//   {
//     size_t key = this->Key2Edge(p1, p2);
//     if(!this->map_edges.count(key)) 
//     { 
//       this->edges.emplace_back(p1, p2);
//       this->map_edges.insert({key, &this->edges.back()});
//     }
//     this->map_edges[key]->LinkAdjFace(&new_face);
//   };
//   for(int i=0; i<pts.size(); i++){
//     for(int j=i+1; j<pts.size(); j++){
//       create_edge(pts[i], pts[j]);
//     }
//   }
//   // create_edge(a, b);
//   // create_edge(a, c);
//   // create_edge(b, c);
// }

bool ConvexHull::Colinear(const Point3D& a, const Point3D& b, const Point3D& c) const
{
  return ((c.z - a.z) * (b.y - a.y) - 
          (b.z - a.z) * (c.y - a.y)) == 0 &&\
         ((b.z - a.z) * (c.x - a.x) - 
          (b.x - a.x) * (c.z - a.z)) == 0 &&\
         ((b.x - a.x) * (c.y - a.y) - 
          (b.y - a.y) * (c.x - a.x)) == 0;
}

bool ConvexHull::BuildFirstHull(std::vector<Point3D>& pointcloud)
{
  const int n = pointcloud.size();
  if(n <= 3)
  {
    std::cout<<"Tetrahedron: points.size() < 4\n";
    return false;    
  }

  int i = 2;
  while(this->Colinear(pointcloud[i], pointcloud[i-1], pointcloud[i-2]))
  {
    if(i++ == n - 1)
    {
      std::cout<<"Tetrahedron: All points are colinear!\n";
      return false;
    }
  }

  Face face(pointcloud[i], pointcloud[i-1], pointcloud[i-2]);

  int j = i;
  while(!this->VolumeSign(face, pointcloud[j], 0))
  {
    if(j++ == n-1)
    {
      std::cout<<"Tetrahedron: All pointcloud are coplanar!\n";
      return false;    
    }
  }

  auto& p1 = pointcloud[i];    auto& p2 = pointcloud[i-1];
  auto& p3 = pointcloud[i-2];  auto& p4 = pointcloud[j];
  p1.processed = p2.processed = p3.processed = p4.processed = true;
  this->AddOneFace(p1, p2, p3, p4);
  this->AddOneFace(p1, p2, p4, p3);
  this->AddOneFace(p1, p3, p4, p2);
  this->AddOneFace(p2, p3, p4, p1);
  return true;

}
 
Point3D ConvexHull::FindInnerPoint(const Face* f, const Edge& e)
{
  for(int i = 0; i < 3; i++)
  {
    if(f->vertices[i] == e.endpoints[0]) continue;
    if(f->vertices[i] == e.endpoints[1]) continue;
    return f->vertices[i];
  } 
}

void ConvexHull::IncreHull(const Point3D& pt)
{
  // Find the illuminated faces (which will be removed later)
  bool vis = false;
  vector<double> col = {pt.x, pt.y, pt.z, 1.0};
  for(auto& face : this->faces)
  {
    // vector<vector<double>> newInv = dyn_inv(face.originalMatrix, face.inverseMatrix, col, 3);
    double newDet = dyn_det(face.originalMatrix, face.inverseMatrix, face.determinant, col, 3);
    // if(VolumeSign(face, pt, 1) < 0) 
    // {
    //   //std::cout<<"face illuminated by pt "<<pt<<" is \n"<<face<<"\n";
    //   face.visible = vis = true;
    // }
    // cout<<"Debug:\n";
    // showMatrix(face.originalMatrix);
    // cout<<"-------------------"<<endl;
    // cout<<VolumeSign(face, pt, 1)<<endl;
    // cout<<newDet<<endl;
    if(newDet < 0) 
    {
      //std::cout<<"face illuminated by pt "<<pt<<" is \n"<<face<<"\n";
      face.visible = vis = true;
    }
  }
  if(!vis) return;

  // Find the edges to make new tagent surface or to be removed
  for(auto it = this->edges.begin(); it != this->edges.end(); it++)
  {
    auto& edge = *it;
    auto& face1 = edge.adjface1;
    auto& face2 = edge.adjface2;

    // Newly added edge
    if(face1 == NULL || face2 == NULL)
    {
      continue;
    }

    // This edge is to be removed because two adjacent faces will be removed 
    else if(face1->visible && face2->visible) 
    {
      edge.remove = true;
    }

    // Edge on the boundary of visibility, which will be used to extend a tagent
    // cone surface.
    else if(face1->visible|| face2->visible) 
    {
      if(face1->visible) std::swap(face1, face2);
      auto inner_pt = this->FindInnerPoint(face2, edge);
      edge.Erase(face2);
      this->AddOneFace(edge.endpoints[0], edge.endpoints[1], pt, inner_pt);
    }
  }
}

void ConvexHull::ConstructHull(std::vector<Point3D>& pointcloud)
{
  if(!this->BuildFirstHull(pointcloud)) return;
  for(const auto& pt : pointcloud)
  {
    if(pt.processed) continue;
    this->IncreHull(pt);
    this->CleanUp();
  }
  this->ExtractExteriorPoints();
}

void ConvexHull::CleanUp()
{
  auto it_edge = this->edges.begin();
  while(it_edge != this->edges.end())
  {
    if(it_edge->remove)
    {
      auto pt1 = it_edge->endpoints[0];
      auto pt2 = it_edge->endpoints[1];
      auto key_to_evict = this->Key2Edge(pt1, pt2);
      this->map_edges.erase(key_to_evict);
      this->edges.erase(it_edge++);
    }
    else it_edge++;
  };
  auto it_face = this->faces.begin();
  while(it_face != this->faces.end())
  {
    if(it_face->visible) this->faces.erase(it_face++);
    else it_face++;
  }
}

void ConvexHull::ExtractExteriorPoints()
{
  std::unordered_set<Point3D, point_hash> exterior_set;
  for(const auto& f : this->faces)
  {
    for(int i =0; i < 3; i++)
      exterior_set.insert(f.vertices[i]);
  }
  this->exterior_points = \
      std::vector<Point3D>(exterior_set.begin(), exterior_set.end());
}

void ConvexHull::Print(const std::string mode = "none")
{
  if(mode == "vertice")
  {
    for(const auto& pt : this->exterior_points)  std::cout<<pt<<"\n";
  }    
  else if(mode == "edge")
  {
    for(const auto& e : this->edges)  std::cout<<(e)<<"\n";
  }
  else if( mode == "face")
  {
    for(const auto& f : this->faces)  std::cout<<f<<"\n";
  }
  else
  {
    std::cout<<"Print Usage: Print {'vertice', 'edge', 'face'}\n";
  }
}

int main()
{
  std::vector<Point3D> tests = \
    {{1,0,0}, {0,1,0}, {0,0,1}, {3,3,3}, {6, 6, 6}};
  ConvexHull C(tests);
  
  std::cout<<"vertices (x, y, z, intensity)\n"; C.Print("vertice");
  std::cout<<"faces\n";  C.Print("face");

  
//     std::vector<Point3D> tests1 = \
//     {{3,0,0}, {0,3,0}, {0,0,3}, {3,3,3}};
//     Face f1({3,0,0}, {0,3,0}, {0,0,3});
//     Point3D p1(3, 3, 3);
//   ConvexHull C1(tests1);
//     cout<<C1.VolumeSign(f1, p1, 1)<<endl;
  return 0;
}