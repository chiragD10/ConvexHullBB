#ifndef CONVEXHULL_H
#define CONVEXHULL_H

#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include "utility.h"


// Defined in CCW
struct Face
{
  Face(const Point3D& p1, const Point3D& p2, const Point3D& p3): visible(false)
      { 
        vertices[0] = p1; vertices[1] = p2; vertices[2] = p3;
        inverseMatrix.resize(3);
        for(int i=0; i<3; i++){
          inverseMatrix[i].resize(3);
        }
        originalMatrix.resize(3);
        for(int i=0; i<3; i++){
          originalMatrix[i].resize(3);
        }
      };

  void Reverse()
  {
    std::swap(vertices[0], vertices[2]); 
    // std::swap()
  };

  friend std::ostream& operator<<(std::ostream& os, const Face& f)
  {
    os << "[face pt1 = "   << f.vertices[0].ToString()
       << " | face pt2 = " << f.vertices[1].ToString()
       << " | face pt3 = " << f.vertices[2].ToString()<<"] ";
    return os;
  }

  bool visible;
  Point3D vertices[3];
  double determinant;
  vector<vector<double>> originalMatrix, inverseMatrix;
};

// struct nFace
// {
//     nFace(int d, const vector<PointnD> &F): visible(false), dim(d)
//     {
//         vertices.resize(d);
//         for(int i=0; i<d; i++){
//           vertices[i] = F[i];
//         }
//     }
//     void Reverse(){std::swap(vertices[0], vertices[1]); };

//     friend std::ostream& operator<<(std::ostream& os, const nFace& f)
//     {
//         for(int i=0; i<f.vertices.size(); i++){
//           os << "[face pt"<<i+1<<" = "<< f.vertices[i].ToString()<<" | ";
//         }
//         os << "] ";
//         return os;
//     }
    
//     int dim;
//     bool visible;
//     vector<PointnD> vertices;
// };

struct Edge
{
  Edge(const Point3D& p1, const Point3D& p2): 
      adjface1(nullptr), adjface2(nullptr), remove(false) 
      { endpoints[0] = p1; endpoints[1] = p2; };
  
  void LinkAdjFace(Face* face) 
  {
    if( adjface1 != NULL && adjface2 != NULL )
    {
      std::cout<<"warning: property violated!\n";
    }
    (adjface1 == NULL ? adjface1 : adjface2)= face;
  };

  void Erase(Face* face) 
  {
    if(adjface1 != face && adjface2 != face) return;
    (adjface1 == face ? adjface1 : adjface2) = nullptr;
  };

  friend std::ostream& operator<<(std::ostream& os, const Edge& e)
  {
    os << "[edge pt1 = "   << e.endpoints[0].ToString()
       << " | edge pt2 = " << e.endpoints[1].ToString() << " ]";
    return os;
  }

  bool remove; 
  Face *adjface1, *adjface2;
  Point3D endpoints[2];
};

// struct nEdge
// {
//   nEdge(const PointnD& p1, const PointnD& p2): adjface1(nullptr), adjface2(nullptr), remove(false)
//   {
//     endpoints[0] = p1; endpoints[1] = p2;
//   }

//   void LinkAdjFace(nFace* face) 
//   {
//     if( adjface1 != NULL && adjface2 != NULL )
//     {
//       std::cout<<"warning: property violated!\n";
//     }
//     (adjface1 == NULL ? adjface1 : adjface2)= face;
//   };

//   void Erase(nFace* face) 
//   {
//     if(adjface1 != face && adjface2 != face) return;
//     (adjface1 == face ? adjface1 : adjface2) = nullptr;
//   };

//   friend std::ostream& operator<<(std::ostream& os, const nEdge& e)
//   {
//     os << "[edge pt1 = "   << e.endpoints[0].ToString()
//        << " | edge pt2 = " << e.endpoints[1].ToString() << " ]";
//     return os;
//   }

//   bool remove; 
//   nFace *adjface1, *adjface2;
//   PointnD endpoints[2];
// };

class ConvexHull
{
  public:
    template<typename T> ConvexHull(const std::vector<T>& points);
    // All major works are conducted upon calling of constructor

    ~ConvexHull() = default;

    template<typename T> bool Contains(T p) const;
    // In out test for a query point (surface point is considered outside)

    const std::list<Face>& GetFaces() const {return this->faces;};
    
    const std::vector<Point3D> GetVertices() const \
        {return this->exterior_points;};
    // Return exterior vertices than defines the convell hull

    void Print(const std::string mode);
    // mode {face, edge, vertice}

    int Size() const {return this->exterior_points.size();};

    double VolumeSign(const Face& f, const Point3D& p, bool debug = 0) const;
    // A point is considered outside of a CCW face if the volume of tetrahedron
    // formed by the face and point is negative. Note that origin is set at p.
    std::list<Face> faces = {};
  private:

    bool Colinear(const Point3D& a, const Point3D& b, const Point3D& c) const;

    bool CoPlanar(Face& f, Point3D& p);


    size_t Key2Edge(const Point3D& a, const Point3D& b) const;
    // Hash key for edge. hash(a, b) = hash(b, a)

    void AddOneFace(const Point3D& a, const Point3D& b, 
        const Point3D& c, const Point3D& inner_pt);
    // Inner point is used to make the orientation of face consistent in counter-
    // clockwise direction

    bool BuildFirstHull(std::vector<Point3D>& pointcloud);
    // Build a tetrahedron as first convex hull

    void IncreHull(const Point3D& p);

    void ConstructHull(std::vector<Point3D>& pointcloud);

    void CleanUp();

    void ExtractExteriorPoints();

    Point3D FindInnerPoint(const Face* f, const Edge& e);
    // for face(a,b,c) and edge(a,c), return b

    std::vector<Point3D> pointcloud = {};
    std::vector<Point3D> exterior_points = {};
    std::list<Edge> edges = {};
    std::unordered_map<size_t, Edge*> map_edges;
};

// class nConvexHull
// {
//   public:
//     template<typename T> nConvexHull(const std::vector<T>& points);
//     // All major works are conducted upon calling of constructor

//     ~nConvexHull() = default;

//     template<typename T> bool Contains(T p) const;
//     // In out test for a query point (surface point is considered outside)

//     const std::list<nFace>& GetFaces() const {return this->faces;};
    
//     const std::vector<PointnD> GetVertices() const \
//         {return this->exterior_points;};
//     // Return exterior vertices than defines the convell hull

//     void Print(const std::string mode);
//     // mode {face, edge, vertice}

//     int Size() const {return this->exterior_points.size();};

//     double VolumeSign(const nFace& f, const PointnD& p, bool debug = 0) const;
//     // A point is considered outside of a CCW face if the volume of tetrahedron
//     // formed by the face and point is negative. Note that origin is set at p.
//   private:

//     bool Colinear(const PointnD& a, const PointnD& b, const PointnD& c) const;

//     bool CoPlanar(nFace& f, PointnD& p);


//     size_t Key2Edge(const PointnD& a, const PointnD& b) const;
//     // Hash key for edge. hash(a, b) = hash(b, a)

//     void AddOneFace(const vector<PointnD> pts, const PointnD& inner_pt);
//     // Inner point is used to make the orientation of face consistent in counter-
//     // clockwise direction

//     bool BuildFirstHull(std::vector<PointnD>& pointcloud);
//     // Build a tetrahedron as first convex hull

//     void IncreHull(const PointnD& p);

//     void ConstructHull(std::vector<PointnD>& pointcloud);

//     void CleanUp();

//     void ExtractExteriorPoints();

//     Point3D FindInnerPoint(const nFace* f, const nEdge& e);
//     // for face(a,b,c) and edge(a,c), return b

//     std::vector<PointnD> pointcloud = {};
//     std::vector<PointnD> exterior_points = {};
//     std::list<nFace> faces = {};
//     std::list<nEdge> edges = {};
//     std::unordered_map<size_t, nEdge*> map_edges;
// };

template<typename T> ConvexHull::ConvexHull(const std::vector<T>& points)
{
  const int n = points.size();
  this->pointcloud.resize(n);
  for(int i = 0; i < n; i++)
  { 
    this->pointcloud[i].x = points[i].x;
    this->pointcloud[i].y = points[i].y;
    this->pointcloud[i].z = points[i].z;
  }
  this->ConstructHull(this->pointcloud);
}

// template<typename T> nConvexHull::nConvexHull(const std::vector<T>& points)
// {
//   const int n = points.size();
//   this->pointcloud.resize(n);
//   for(int i = 0; i < n; i++)
//   { 
//     this->pointcloud[i] = points[i];
//   }
//   this->ConstructHull(this->pointcloud);
// }

template<typename T> bool ConvexHull::Contains(T p) const
{
  Point3D pt(p.x, p.y, p.z);
  for(auto& face : this->faces)
  {
    if(VolumeSign(face, pt) <= 0) return false;
  }
  return true;
}

// template<typename T> bool nConvexHull::Contains(T p) const
// {
//   PointnD pt(p.dim, p.X);
//   for(auto& face : this->faces)
//   {
//     if(VolumeSign(face, pt) <= 0) return false;
//   }
//   return true;
// }
#endif