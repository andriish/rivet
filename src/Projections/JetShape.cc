// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/JetShape.hh"
#include "Rivet/Cmp.hh"


namespace Rivet {

  int JetShape::compare(const Projection& p) const {
    const JetShape& other = dynamic_cast<const JetShape&>(p);
    if (_jetaxes.size() != other._jetaxes.size()) return 1;
    /// @todo Add comparison of vecs
    // for (size_t i = 0; i < _jetaxes.size(); ++i) {
    //   if (_jetaxes[i] != other._jetaxes[i]) return 1;
    // }
    return 0;
  }


  void JetShape::project(const Event& e) {
    Log& log = getLog();
    
    // Project into vetoed final state if not done yet (supposed to be done by jet algorithm)
    const VetoedFinalState& vfs = e.applyProjection(_vfsproj);
    
    // Clear for each event and resize with zero vectors
    /// @todo Eliminate pointers to simplify this sort of thing:
    for (size_t i = 0; i < _diffjetshapes.size(); ++i) {
      _diffjetshapes[i].clear();
    }
    _diffjetshapes.clear();

    for (size_t i=0; i<_jetaxes.size(); ++i) {
      vector<double> diffjetshape(_nbins, 0.);
      _diffjetshapes.push_back(diffjetshape);
    }


    /// @todo Simplify
    for (size_t i = 0; i < _intjetshapes.size(); ++i) {
      _intjetshapes[i].clear();
    }
    _intjetshapes.clear();

    for (size_t i = 0; i < _jetaxes.size(); ++i) {
      vector<double> intjetshape(_nbins, 0.);
      _intjetshapes.push_back(intjetshape);
    }
    
    _PsiSlot.clear();
    _PsiSlot.resize(_jetaxes.size(), 0.);

    
    
    // Determine jet shapes
    double y1, y2, eta1, eta2, phi1, phi2, drad;
    /// @todo Fix!
    double dradmin = TWOPI; //< Dummy assignment, to avoid compile warning
    int dradminind = 0; //< Dummy asignment, to avoid compile warning
    if (_jetaxes.size() > 0) {
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

        /// @todo Clear up all these commented out bits. Version control means we can always get them back.
        for (size_t i = 0; i < _nbins; ++i) {
          if (dradmin<_rmin+(i+1)*_interval) {
            _intjetshapes[dradminind][i] += p->getMomentum().vector3().polarRadius();
            if (dradmin>_rmin+i*_interval) {
              _diffjetshapes[dradminind][i] += p->getMomentum().vector3().polarRadius()/_interval;
            }
          }
        }
        
        if (dradmin < _r1minPsi/_rmax) {
          _PsiSlot[dradminind] += p->getMomentum().vector3().polarRadius();
        }
        
      }
      
      
      // Normalize to total pT
      for (size_t j = 0; j < _jetaxes.size(); j++) {
        /// @todo Clear up all these commented out bits. Version control means we can always get them back.
        if (_intjetshapes[j][_nbins-1] > 0.) {
          _PsiSlot[j] /= _intjetshapes[j][_nbins-1];
          for (size_t i = 0; i < _nbins; ++i) {
            _diffjetshapes[j][i] /= _intjetshapes[j][_nbins-1];
            _intjetshapes[j][i] /= _intjetshapes[j][_nbins-1];
          }
        }
      }
      
    }
    
    log << Log::DEBUG << "Done" << endl;
  }
  
  
}
