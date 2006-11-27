/// $Id: $

#include "Rivet/Tools/Math/LorentzVector.hh"
#include <iostream>
#include <climits>

namespace Rivet {

  const LorentzVector LorentzVector::X = LorentzVector(1.0, 1.0, 0.0, 0.0);
  
  const LorentzVector LorentzVector::Y = LorentzVector(1.0, 0.0, 1.0, 0.0);
  
  const LorentzVector LorentzVector::Z = LorentzVector(1.0, 0.0, 0.0, 1.0);
  

  /// Default constructor
  LorentzVector::LorentzVector():
    ct_(0.0), x_(0.0), y_(0.0), z_(0.0), 
    rapidity_(0.0), phi_(0.0), perp_(0.0), m_(0.0),
    status_(CONSISTENT) 
  {  }


  /// Cartesian constructor
  LorentzVector::LorentzVector(double ct, double x, double y, double z):
    ct_(ct), x_(y), y_(y), z_(z), 
    rapidity_(0.0), phi_(0.0), perp_(0.0), m_(0.0),
    status_(BADSNOW) 
  {  }


  LorentzVector::~LorentzVector() {}
  

  //   void LorentzVector::add(const LorentzVector& p){
  //     this->operator+=(p);
  //   }


  ///////////////////////////////////////////////////////////////


  LorentzVector operator*(double num, const LorentzVector& vec) {
    return LorentzVector(num * vec.ct(), num * vec.x(), num * vec.y(), num * vec.z());
  }

  LorentzVector operator*(const LorentzVector& vec, double num) {
    return LorentzVector(num * vec.ct(), num * vec.x(), num * vec.y(), num * vec.z());
  }

  LorentzVector operator/(const LorentzVector& vec, double num) {
    return LorentzVector(vec.ct() / num, vec.x() / num, vec.y() / num, vec.z() / num);
  }

  std::ostream& operator<<(std::ostream& os, LorentzVector vec) {
    os << vec.toString();
    return os;
  }
  
}
