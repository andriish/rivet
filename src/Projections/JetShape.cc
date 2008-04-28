// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  int JetShape::compare(const Projection& p) const {
    PCmp fscmp = mkNamedPCmp(p, "FS");
    if (fscmp == PCmp::EQUIVALENT) return PCmp::EQUIVALENT;
    const JetShape& other = dynamic_cast<const JetShape&>(p);
    return cmp(_jetaxes.size(), other._jetaxes.size());
    /// @todo Add proper comparison of vecs
    // for (size_t i = 0; i < _jetaxes.size(); ++i) {
    //   if (_jetaxes[i] != other._jetaxes[i]) return 1;
    // }
  }


  void JetShape::project(const Event& e) {
    // Clear for each event and resize with zero vectors
    _diffjetshapes.clear();
    _intjetshapes.clear();
    for (size_t i = 0; i < _jetaxes.size(); ++i) {
      vector<double> tmp(_nbins, 0.0);
      _diffjetshapes.push_back(tmp);
      _intjetshapes.push_back(tmp);
    }
    _PsiSlot.clear();
    _PsiSlot.resize(_jetaxes.size(), 0.0);

    // Determine jet shapes
    double y1, y2, eta1, eta2, phi1, phi2, drad;
    double dradmin = TWOPI;
    int dradminind = 0;
    if (_jetaxes.size() > 0) {
      const VetoedFinalState& vfs = applyProjection<VetoedFinalState>(e, "FS");
      for (ParticleVector::const_iterator p = vfs.particles().begin(); p != vfs.particles().end(); ++p) {
        for (size_t j = 0; j < _jetaxes.size(); ++j) {
          y1 = _jetaxes[j].rapidity();
          y2 = p->getMomentum().rapidity();
          eta1 = _jetaxes[j].vector3().pseudorapidity();
          eta2 = p->getMomentum().vector3().pseudorapidity();
          phi1 = _jetaxes[j].vector3().azimuthalAngle();
          phi2 = p->getMomentum().vector3().azimuthalAngle();
          
          if (_distscheme == SNOWMASS) {
            drad = delta_rad(eta1, phi1, eta2, phi2);
          } else { // _distscheme = ENERGY
            drad = delta_rad(y1, phi1, y2, phi2);
          }
          
          if (j == 0 || drad < dradmin) {
            dradminind = j;
            dradmin = drad;
          }
        }

        for (size_t i = 0; i < _nbins; ++i) {
          if (dradmin < _rmin+(i+1)*_interval) {
            _intjetshapes[dradminind][i] += p->getMomentum().vector3().polarRadius();
            if (dradmin > _rmin+i*_interval) {
              _diffjetshapes[dradminind][i] += p->getMomentum().vector3().polarRadius()/_interval;
            }
          }
        }
        
        if (dradmin < _r1minPsi) {
          _PsiSlot[dradminind] += p->getMomentum().vector3().polarRadius();
        }
      }
     
      
      // Normalize to total pT
      for (size_t j = 0; j < _jetaxes.size(); j++) {
        if (_intjetshapes[j][_nbins-1] > 0.) {
          _PsiSlot[j] /= _intjetshapes[j][_nbins-1];
          for (size_t i = 0; i < _nbins; ++i) {
            _diffjetshapes[j][i] /= _intjetshapes[j][_nbins-1];
            _intjetshapes[j][i] /= _intjetshapes[j][_nbins-1];
          }
        }
      }
      
    }
  }
  
  
}
