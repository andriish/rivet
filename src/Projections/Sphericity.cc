// -*- C++ -*-
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Cmp.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/Utils.hh"


namespace Rivet {

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

    // Get final state.
    const FinalState& fs = e.applyProjection(_fsproj);
    const ParticleVector prts = fs.particles();

    // Return (with "safe nonsense" sphericity params) if there are no final state particles.
    if (prts.empty()) {
      log << Log::DEBUG << "No particles in final state..." << endl; 
      return;
    }

    // Iterate over all the final state particles.
    Matrix3 mMom;
    double totalMomentum = 0.0;
    for (ParticleVector::const_iterator p = prts.begin(); p != prts.end(); ++p) {

      // Get the momentum vector for the final state particle.
      const FourMomentum lv = p->getMomentum();
      const Vector3 p3 = lv.vector3();

      // Build the (regulated) normalising factor.
      totalMomentum += pow(p3.mod(), _regparam);

      // Build (regulated) quadratic momentum components.
      const double regfactor = pow(p3.mod(), _regparam-2);
      if (!fuzzyEquals(regfactor, 1.0)) {
        log << Log::DEBUG << "Regfactor (r=" << _regparam << ") = " << regfactor << endl;
      }
      // mMom[0][0] += regfactor * lv.x() * lv.x(); 
      // mMom[1][1] += regfactor * lv.y() * lv.y();
      // mMom[2][2] += regfactor * lv.z() * lv.z();
      // mMom[0][1] += regfactor * lv.y() * lv.x(); 
      // mMom[1][0] += regfactor * lv.x() * lv.y();
      // mMom[0][2] += regfactor * lv.z() * lv.x(); 
      // mMom[2][0] += regfactor * lv.x() * lv.z();
      // mMom[1][2] += regfactor * lv.z() * lv.y();
      // mMom[2][1] += regfactor * lv.y() * lv.z(); 

      Matrix3 mMomPart;
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          mMomPart.set(i,j, p3[i]*p3[j]);
        }
      } 
      mMom += regfactor * mMomPart;
    }

    // Normalise to total (regulated) momentum.
    mMom /= totalMomentum;
    log << Log::DEBUG << "Momentum tensor = " << endl << mMom << endl;

    // Check that the matrix is symmetric.
    const bool isSymm = mMom.isSymm();
    if (!isSymm) {
      log << Log::ERROR << "Error: momentum tensor not symmetric (r=" << _regparam << ")" << endl;
      log << Log::ERROR << "[0,1] vs. [1,0]: " << mMom.get(0,1) << ", " << mMom.get(1,0) << endl;
      log << Log::ERROR << "[0,2] vs. [2,0]: " << mMom.get(0,2) << ", " << mMom.get(2,0) << endl;
      log << Log::ERROR << "[1,2] vs. [2,1]: " << mMom.get(1,2) << ", " << mMom.get(2,1) << endl;
    }
    // If not symmetric, something's wrong (we made sure the error msg appeared first).
    assert(isSymm); 

    // Diagonalize momentum matrix.
    const EigenSystem<3> eigen3 = diagonalize(mMom);
    log << Log::DEBUG << "Diag momentum tensor = " << endl << eigen3.getDiagMatrix() << endl;
    
    // Reset and set eigenvalue/vector parameters.
    _lambdas.clear();
    _sphAxes.clear();
    const EigenSystem<3>::EigenPairs epairs = eigen3.getEigenPairs();
    assert(epairs.size() == 3);
    for (size_t i = 0; i < 3; ++i) {
      _lambdas.push_back(epairs[i].first);
      _sphAxes.push_back(Vector3(epairs[i].second));
    }

    // Debug output.
    log << Log::DEBUG << "Sum of lambdas = " << lambda1() + lambda2() + lambda3() << endl;
    log << Log::DEBUG << "Lambdas = (" 
        << lambda1() << ", " << lambda2() << ", " << lambda3() << ")" << endl;
    log << Log::DEBUG << "Vectors = " 
        << sphericityAxis() << ", "
        << sphericityMajorAxis() << ", " 
        << sphericityMinorAxis() << ")" << endl;
  }

}
