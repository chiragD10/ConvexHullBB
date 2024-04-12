#include<iostream>
#include <vector>
#include<string>
using namespace std;

struct Point3D {
  float x;
  float y;
  float z;
  float intensity;
  bool processed;
  Point3D() = default;
  Point3D(float _x, float _y, float _z):\
      x(_x), y(_y), z(_z), intensity(0), processed(false) {}
  Point3D(float _x, float _y, float _z, float _i):\
      x(_x), y(_y), z(_z), intensity(_i), processed(false) {}

  bool operator ==(const Point3D& pt) const
  {
      return x == pt.x && y == pt.y && z == pt.z;
  }

  float operator *(const Point3D& pt) const
  {
      return x * pt.x + y * pt.y + z * pt.z;
  }

  Point3D operator *(const float factor) const
  {
      return Point3D(x*factor, y*factor, z*factor);
  }

  Point3D operator /(const float factor) const
  {
      return Point3D(x/factor, y/factor, z/factor);
  }

  Point3D operator +(const Point3D& pt) const
  {
      return Point3D(x+pt.x, y+pt.y, z+pt.z);
  }
  Point3D operator -(const Point3D& pt) const
  {
      return Point3D(x-pt.x, y-pt.y, z-pt.z);
  }

  std::string ToString() const
  {
    return to_string(x).substr(0,4) + ", " + 
           to_string(y).substr(0,4) + ", " + 
           to_string(z).substr(0,4) + ", " + 
           to_string(intensity).substr(0,4);
  };
  
  friend ostream& operator<<(ostream& os, const Point3D& p)
  {
    os << "["<< p.ToString() << "] ";
    return os;
  }
  
};

typedef vector<Point3D> PointStack;

struct point_hash
{
  size_t operator() (const Point3D& p) const
  {
      string sx, sy, sz;
      sx = to_string(p.x);
      sy = to_string(p.y);
      sz = to_string(p.z);
      return hash<string>{}(sx+sy+sz);
  }
};

struct PointnD {
    int dim;
    vector<double> X;
    float intensity;
    bool processed;
    PointnD() = default;
    PointnD(int d, vector<double> _X): dim(d), X(_X), intensity(0), processed(false) {}
    PointnD(int d, vector<double> _X, float _i):\
        dim(d), X(_X), intensity(_i), processed(false) {}
    
    bool operator ==(const PointnD& pt) const
    {
        for(int i=0; i<dim; i++){
            if(X[i] != pt.X[i]){
                return false;
            }
        }
        return true;
    }

    float operator *(const PointnD& pt) const
    {
        double res = 0;
        for(int i=0; i<dim; i++){
            res += (pt.X[i] * X[i]);
        }
        return res;
    }

    PointnD operator *(const float factor) const
    {
        vector<double> Y(dim);
        for(int i=0; i<dim; i++){
            Y[i] = factor * X[i];
        }
        return PointnD(dim, Y);
    }

    PointnD operator /(const float factor) const
    {
        vector<double> Y(dim);
        for(int i=0; i<dim; i++){
            Y[i] = X[i] / factor;
        }
        return PointnD(dim, Y);
    }

    PointnD operator +(const PointnD& pt) const
    {
        vector<double> Y(dim);
        for(int i=0; i<dim; i++){
            Y[i] = X[i] + pt.X[i];
        }
        return PointnD(dim, Y);
    }
    PointnD operator -(const PointnD& pt) const
    {
        vector<double> Y(dim);
        for(int i=0; i<dim; i++){
            Y[i] = X[i] - pt.X[i];
        }
        return PointnD(dim, Y);
    }

    std::string ToString() const
    {
        string str = "";
        for(int i=0; i<dim; i++){
            str += to_string(X[i]).substr(0,4) + ", ";
        }
        str += to_string(intensity).substr(0,4);
        return str;
    };
  
    friend ostream& operator<<(ostream& os, const PointnD& p)
    {
        os << "["<< p.ToString() << "] ";
        return os;
    }
};

typedef vector<PointnD> nPointStack;

struct npoint_hash
{
  size_t operator() (const PointnD& p) const
  {
    string str = "";
    for(int i=0; i<p.dim; i++){
        str += to_string(p.X[i]);
    }
    return hash<string>{}(str);
  }
};
