// -*- C++ -*-

#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Cmp.hh"
#include "Rivet/RivetCLHEP.hh"

using namespace Rivet;
using namespace std;


int ParisiTensor::compare(const Projection& p) const {
  const ParisiTensor& other = dynamic_cast<const ParisiTensor&>(p);
  int sphcmp = pcmp(_sphproj, other._sphproj);
  return sphcmp;
}


// RivetInfo ParisiTensor::getInfo() const {
//   return Projection::getInfo() + _sphproj.getInfo();
// }


void ParisiTensor::project(const Event & e) {
  // Apply sphericity projection to event
  Sphericity sph = e.applyProjection(_sphproj);

  // Set parameters
  _lambda[0] = sph.lambda1();
  _lambda[1] = sph.lambda2();
  _lambda[2] = sph.lambda3();
  _C = 3 * ( lambda1()*lambda2() + lambda1()*lambda3() + lambda2()*lambda3() );
  _D = 27 * lambda1() * lambda2() * lambda3();

}
