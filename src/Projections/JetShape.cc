// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  /// Constructor.
  JetShape::JetShape(const VetoedFinalState& vfsp, 
                     const vector<FourMomentum>& jetaxes, 
                     double rmin, double rmax, double interval, 
                     double r1minPsi, DeltaRScheme distscheme)
    : _jetaxes(jetaxes), 
      _rmin(rmin), _rmax(rmax), _interval(interval), 
      _r1minPsi(r1minPsi), _distscheme(distscheme)
  {
    setName("JetShape");
    _nbins = int(round((rmax-rmin)/interval));
    addProjection(vfsp, "FS");
  }


  int JetShape::compare(const Projection& p) const {
    PCmp fscmp = mkNamedPCmp(p, "FS");
    if (fscmp == PCmp::EQUIVALENT) return PCmp::EQUIVALENT;
    const JetShape& other = dynamic_cast<const JetShape&>(p);
    PCmp sizecmp = cmp(_jetaxes.size(), other._jetaxes.size());
    if (sizecmp == PCmp::EQUIVALENT) return PCmp::EQUIVALENT;
    // Check that each axis vector has an equivalent in the other set
    foreach (const FourMomentum& j1, _jetaxes) {
      bool match = false;
      foreach (const FourMomentum& j2, other._jetaxes) {
        if (isZero(angle(j1, j2)) && fuzzyEquals(j1.mod2(), j1.mod2())) {
          match = true;
          break;
        }
        if (!match) return PCmp::UNDEFINED;
      }
    }
    return PCmp::EQUIVALENT;
  }


  void JetShape::clear() {
    // Reset vectors for each event
    _diffjetshapes.clear();
    _intjetshapes.clear();
    for (size_t i = 0; i < _jetaxes.size(); ++i) {
      const vector<double> tmp(_nbins, 0.0);
      _diffjetshapes.push_back(tmp);
      _intjetshapes.push_back(tmp);
    }
    _PsiSlot.clear();
    _PsiSlot.resize(_jetaxes.size(), 0.0);
  }


  void JetShape::project(const Event& e) {
    // Reset for new event
    clear();

    if (!_jetaxes.empty()) {
      const VetoedFinalState& vfs = applyProjection<VetoedFinalState>(e, "FS");
      foreach (const Particle& p, vfs.particles()) {
        double drad_min = TWOPI;
        size_t i_drad_min = 0;
        
        // Identify "best match" jet axis for this particle
        for (size_t j = 0; j < _jetaxes.size(); ++j) {
          const double drad = deltaR(_jetaxes[j], p.momentum(), _distscheme);
          if (drad < drad_min) {
            i_drad_min = j;
            drad_min = drad;
          }
        }

        // Fill diff & int jet shape histos for closest jet axis
        /// @todo Actually use histograms here, rather than doing the binning by hand
        /// @todo Calculate int jet shape from diff jet shape histo (YODA)
        for (size_t i = 0; i < _nbins; ++i) {
          if (drad_min < _rmin + (i+1)*_interval) {
            _intjetshapes[i_drad_min][i] += p.momentum().pT();
            if (drad_min > _rmin + i*_interval) {
              _diffjetshapes[i_drad_min][i] += p.momentum().pT()/_interval;
            }
          }
        }

        // Sum pT of closest match jet axes for dr < _r1minPsi
        /// @todo Calculate int [0.0, 0.3] jet shape from diff jet shape histo (YODA)
        if (drad_min < _r1minPsi) {
          _PsiSlot[i_drad_min] += p.momentum().pT();
        }

      }
     
      
      // Normalize to total pT
      for (size_t j = 0; j < _jetaxes.size(); j++) {
        const double psimax = _intjetshapes[j][_nbins-1];
        if (psimax > 0.0) {
          _PsiSlot[j] /= psimax;
          for (size_t i = 0; i < _nbins; ++i) {
            _diffjetshapes[j][i] /= psimax;
            _intjetshapes[j][i] /= psimax;
          }
        }
      }

      
    }
  }
  
  
}
