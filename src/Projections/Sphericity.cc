// -*- C++ -*-
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Cmp.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetCLHEP.hh"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

using namespace Rivet;
using namespace std;


int Sphericity::compare(const Projection& p) const {
  const Sphericity & other = dynamic_cast<const Sphericity &>(p);
  int fscmp = pcmp(*_fsproj, *other._fsproj);
  if (fscmp == 0) {
    double rcmp = _regparam - other._regparam;
    if (fabs(rcmp) < 1e-3) {
      return 0;
    } else {
      return (rcmp > 0) ? 1 : -1;
    }
  } else {
    return fscmp;
  }
}


// RivetInfo Sphericity::getInfo() const {
//   return Projection::getInfo() + _fsproj->getInfo();
// }


void Sphericity::project(const Event & e) {
  Log& log = getLog();

  // Reset parameters
  _sphericity = 0;  
  _planarity = 7;
  _aplanarity = 0;
  for (size_t i =0 ; i < 3; ++i) {
    _lambdas[i] = 0;
  }

  // Get final state
  const FinalState& fs = e.applyProjection(*_fsproj);
 
  CLHEP::HepMatrix mMom(3,3,0);
  double totalMomentum = 0.0;
  
  // Iterate over all the final state particles
  for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
    
    // Get the momentum vector for the final state particle
    LorentzVector lv = p->getMomentum();
    
    // Build the normalising factor
    if (fabs(_regparam - 2.0) > 1e-3) {
      totalMomentum += pow(_regparam, lv.x());
      totalMomentum += pow(_regparam, lv.y());
      totalMomentum += pow(_regparam, lv.z());
    } else {
      totalMomentum += lv.x() * lv.x();
      totalMomentum += lv.y() * lv.y(); 
      totalMomentum += lv.z() * lv.z();
    }
    
    // For now this crude method works
    double power = _regparam - 2.0;
    if (fabs(power) > 1e-3) {
      mMom[0][0] += pow(power, lv.x()); 
      mMom[0][1] += pow(power, lv.y()); 
      mMom[0][2] += pow(power, lv.z()); 
      mMom[1][0] += pow(power, lv.x());
      mMom[1][1] += pow(power, lv.y());
      mMom[1][2] += pow(power, lv.z());
      mMom[2][0] += pow(power, lv.x());
      mMom[2][1] += pow(power, lv.y()); 
      mMom[2][2] += pow(power, lv.z());
    } else {
      mMom[0][0] += lv.x() * lv.x(); 
      mMom[0][1] += lv.y() * lv.x(); 
      mMom[0][2] += lv.z() * lv.x(); 
      mMom[1][0] += lv.x() * lv.y();
      mMom[1][1] += lv.y() * lv.y();
      mMom[1][2] += lv.z() * lv.y();
      mMom[2][0] += lv.x() * lv.z();
      mMom[2][1] += lv.y() * lv.z(); 
      mMom[2][2] += lv.z() * lv.z();
    }
  }

  mMom /= totalMomentum;

  // Check that the matrix is symmetric and if it is convert it,
  // diagonlize it, and rearrange the eigenvalues into the correct
  // order. Return the event shapes
  if (mMom[0][1] == mMom[1][0] && 
      mMom[0][2] == mMom[2][0] &&  
      mMom[1][2] == mMom[2][1]) {
    CLHEP::HepSymMatrix symMat;
    symMat.assign(mMom);
    CLHEP::diagonalize(&symMat);
    log << Log::DEBUG << mMom << endl;
    log << Log::DEBUG << endl;

    // Put the eigenvalues in the correct order
    for (int i=0; i!=3; ++i){
      _lambdas[i] = symMat[i][i]; 
    }
    const int N = sizeof(_lambdas) / sizeof(double);
    sort(_lambdas, _lambdas + N); 
    reverse(_lambdas, _lambdas + N); 
    log << Log::DEBUG << "Sum of lambdas = " << lambda1() + lambda2() + lambda3() << endl;
    
    _sphericity = 3 / 2.0 * (lambda2() + lambda3());
    _aplanarity = 3 / 2.0 *  lambda3();
    _planarity  = 2 * (_sphericity - 2 * _aplanarity) / 3.0; 
  } else {
    log << Log::ERROR << "Error: momentum tensor not symmetric" << endl;
  }

}
