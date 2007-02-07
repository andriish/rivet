// -*- C++ -*-

#include "Rivet/Projections/Thrust.hh"
#include "Rivet/RivetCLHEP.hh"

using namespace Rivet;
using namespace CLHEP;


Thrust::~Thrust() {}


int Thrust::compare(const Projection & p) const {
  return 0;
}


double Thrust::calcThrust(const vector<Vector3>& momenta, const Vector3& n) {
  // Make sure the direction vector is normalised
  const Vector3 direction = n / n.modulus();

  // Run over momenta of particles in final state
  double sumOfMomComponentModuli = 0;
  double sumOfMomModuli = 0;
  for (vector<Vector3>::const_iterator p = momenta.begin(); p != momenta.end(); ++p) {
    sumOfMomComponentModuli += fabs( p->dot(direction) );
    sumOfMomModuli += p->modulus();
  }
  double thrust = sumOfMomComponentModuli / sumOfMomModuli;
  return thrust;
}


const vector<Vector3> Thrust::calcThrustAxes(const FinalState& fs) {
  // Make a vector of the three-momenta in the final state
  vector<Vector3> threeMomenta;
  for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
    Vector3 p3 = p->momentum.vec();
    threeMomenta.push_back(p3);
  }

  // Handle special case of thrust = 1 if there are only 2 particles

  // Handle special case where thrust = max{x1,x2,x3} if there are 3 particles

  // For more general case, build all 3-momentum subsets (well, half of them due to
  // symmetry of n and N-n subsets) and identify the thrust axis with the subset which 
  // maximises the thrust scalar. Identify the thrust major and minor axes by...?
  //Vector3 n = sum of 3-momenta
  double maxThrust = 0;
  Vector3 maxThrustAxis;
  for (...) {
    // Use std::bitfield
    Vector3 n = ...;
    thrust = calcThrust(threeMomenta, n);
    if (thrust > maxThrust) {
      maxThrust = thrust;
      maxThrustVector = n;
    }
  }

  vector<Vector3> thrustAxes = ...;
  return thrustAxes;
}


void Thrust::orderEigenvalues(vector<double>& eigenvalues) {
  // Check that there are 3 eigenvalues
  assert(eigenvalues.size() == 3);

  // Put the eigenvalues in the correct order
  rsort(order.begin(), order.end()); 
}


void Thrust::project(const Event & e) {
  // Reset parameters
  for (int i = 0; i < 3; ++i) lambdas[i] = 0;

  const FinalState& fs = e.addProjection(fsproj);

  thrustAxes_ = calcThrustAxes();

}
