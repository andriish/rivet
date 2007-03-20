#include <stdexcept>
#include <string>
#include <ostream>
#include <iostream>


class Vector4 {
public:

  //  virtual ~Vector4() = 0;

public:
  double operator[](size_t index) {
    if (index < 0 || index > 3) {
      throw std::runtime_error("Tried to access an invalid 4-vector index.");
    } else {
      return _elements[index];
    }
  }

  double invariant() const {
    double metricdiag[4] = {1, -1, -1, -1};
    double invariant = 0;
    for (size_t i = 0; i < 4; ++i) {
      invariant += metricdiag[i] * _elements[i]*_elements[i];
    }
    return invariant;
  }

protected:
  double _elements[4];
};


class Vector4X : public Vector4 {
public:
  Vector4X(double t, double x, double y, double z) {
    this->t(t);
    this->x(x);
    this->y(y);
    this->z(z);
  }

  ~Vector4X() {}

public:
  double t() const { return _elements[0]; }
  double x() const { return _elements[1]; }
  double y() const { return _elements[2]; }
  double z() const { return _elements[3]; }
  Vector4X& t(double t) { _elements[0] = t; return *this; }
  Vector4X& x(double x) { _elements[1] = x; return *this; }
  Vector4X& y(double y) { _elements[2] = y; return *this; }
  Vector4X& z(double z) { _elements[3] = z; return *this; }
  double interval() const { return this->invariant(); }
};


class Vector4P : public Vector4 {
public:
  Vector4P(double E, double px, double py, double pz) {
    this->E(E);
    this->px(px);
    this->py(py);
    this->pz(pz);
  }

  ~Vector4P() {}

public:
  double E() const { return _elements[0]; }
  double px() const { return _elements[1]; }
  double py() const { return _elements[2]; }
  double pz() const { return _elements[3]; }
  Vector4P& E(double E) { _elements[0] = E; return *this; }
  Vector4P& px(double px) { _elements[1] = px; return *this; }
  Vector4P& py(double py) { _elements[2] = py; return *this; }
  Vector4P& pz(double pz) { _elements[3] = pz; return *this; }
  double mass() const { return this->invariant(); }
};


std::ostream& operator<<(std::ostream& out, Vector4& v4) {
  out << "("  << v4[0] 
      << ", " << v4[1]
      << ", " << v4[2]
      << ", " << v4[3] << ")";
  return out;
}


int main() {
  Vector4X a(1,0,0,0);
  std::cout << a << ": interval = " << a.interval() << std::endl;
  a.z(1);
  std::cout << a << ": interval = " << a.interval() << std::endl;
  a.y(2).z(3);
  std::cout << a << ": interval = " << a.interval() << std::endl;
  return EXIT_SUCCESS;
}
