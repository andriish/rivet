// -*- C++ -*-
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Cmp.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/Utils.hh"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

using namespace Rivet;
using namespace std;


int Sphericity::compare(const Projection& p) const {
  const Sphericity& other = dynamic_cast<const Sphericity&>(p);
  int fscmp = pcmp(_fsproj, other._fsproj);
  if (fscmp != 0) return fscmp;
  if (fuzzyEquals(_regparam, other._regparam)) return 0;
  return (_regparam > other._regparam) ? 1 : -1;
}


void Sphericity::project(const Event& e) {
  Log& log = getLog();
  log << Log::DEBUG << "Calculating sphericity with r = " << _regparam << endl;

  // Reset parameters
  _sphericity = -1;  
  _planarity = -1;
  _aplanarity = -1;
  for (size_t i =0 ; i < 3; ++i) _lambdas[i] = -1;

  // Get final state.
  const FinalState& fs = e.applyProjection(_fsproj);
  const ParticleVector prts = fs.particles();

  // Return (with "safe nonsense" sphericity params) if there are no final state particles.
  if (prts.empty()) {
    log << Log::DEBUG << "No particles in final state..." << endl; 
    return;
  }

  // Iterate over all the final state particles.
  CLHEP::HepMatrix mMom(3,3,0);
  double totalMomentum = 0.0;
  for (ParticleVector::const_iterator p = prts.begin(); p != prts.end(); ++p) {
    
    // Get the momentum vector for the final state particle.
    LorentzVector lv = p->getMomentum();
    
    // Build the (regulated) normalising factor.
    totalMomentum += pow(lv.vect().mag(), _regparam);
    
    // Build (regulated) quadratic momentum components.
    const double regfactor = pow(lv.vect().mag(), _regparam-2);
    if (!fuzzyEquals(regfactor, 1.0)) {
      log << Log::DEBUG << "Regfactor (r=" << _regparam << ") = " << regfactor << endl;
    }
    mMom[0][0] += regfactor * lv.x() * lv.x(); 
    mMom[1][1] += regfactor * lv.y() * lv.y();
    mMom[2][2] += regfactor * lv.z() * lv.z();
    mMom[0][1] += regfactor * lv.y() * lv.x(); 
    mMom[1][0] += regfactor * lv.x() * lv.y();
    mMom[0][2] += regfactor * lv.z() * lv.x(); 
    mMom[2][0] += regfactor * lv.x() * lv.z();
    mMom[1][2] += regfactor * lv.z() * lv.y();
    mMom[2][1] += regfactor * lv.y() * lv.z(); 
  }

  // Normalise to total (regulated) momentum.
  mMom /= totalMomentum;
  log << Log::DEBUG << "Momentum tensor = " << endl << mMom << endl;

  // Check that the matrix is symmetric.
  const bool isSymm = 
    fuzzyEquals(mMom[0][1], mMom[1][0]) && 
    fuzzyEquals(mMom[0][2], mMom[2][0]) && 
    fuzzyEquals(mMom[1][2], mMom[2][1]);
  if (!isSymm) {
    log << Log::ERROR << "Error: momentum tensor not symmetric (r=" << _regparam << ")" << endl;
    log << Log::ERROR << "[0,1] vs. [1,0]: " << mMom[0][1] << ", " << mMom[1][0] << endl;
    log << Log::ERROR << "[0,2] vs. [2,0]: " << mMom[0][2] << ", " << mMom[2][0] << endl;
    log << Log::ERROR << "[1,2] vs. [2,1]: " << mMom[1][2] << ", " << mMom[2][1] << endl;
  }
  // If not symmetric, something's wrong (we made sure the error msg appeared first).
  assert(isSymm); 

  // Convert to a SymMatrix and diagonalize.
  CLHEP::HepSymMatrix symMat;
  symMat.assign(mMom);
  CLHEP::diagonalize(&symMat);
  log << Log::DEBUG << "Diag momentum tensor = " << endl << symMat << endl;
  
  // Put the eigenvalues in the correct order.
  for (int i=0; i!=3; ++i){
    _lambdas[i] = symMat[i][i]; 
  }
  const int N = sizeof(_lambdas) / sizeof(double);
  sort(_lambdas, _lambdas + N); 
  reverse(_lambdas, _lambdas + N); 
  log << Log::DEBUG << "Sum of lambdas = " << lambda1() + lambda2() + lambda3() << endl;
  
  // Construct standard eigenvalue combinations.
  _sphericity = 3 / 2.0 * (lambda2() + lambda3());
  _aplanarity = 3 / 2.0 *  lambda3();
  _planarity  = 2 * (_sphericity - 2 * _aplanarity) / 3.0; 
  
}
