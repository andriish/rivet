// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  /// Constructor.
  JetShape::JetShape(const JetAlg& jetalg,
                     double rmin, double rmax, double interval,
                     DeltaRScheme distscheme)
    : _rmin(rmin), _rmax(rmax), _interval(interval),
      _distscheme(distscheme)
  {
    setName("JetShape");
    _nbins = int(round((rmax-rmin)/interval));
    addProjection(jetalg, "Jets");
  }


  int JetShape::compare(const Projection& p) const {
    PCmp jcmp = mkNamedPCmp(p, "Jets");
    return jcmp;
    /// @todo Also compare bin edges.
  }


  void JetShape::clear() {
    _diffjetshapes = vector<double>(_nbins, 0.0);
  }


  void JetShape::calc(const Jets& jets) {
    clear();

    foreach (const Jet& j, jets) {
      FourMomentum pj = j.momentum();
      foreach (const Particle& p, j.particles()) {
        const double dR = deltaR(pj, p.momentum());
        size_t dRindex = int(floor((dR-_rmin)/(_rmax - _rmin)));
        assert(dRindex < _nbins);
        _diffjetshapes[dRindex] += p.momentum().pT();
      }
    }

    // Normalize to total pT
    double integral = 0.0;
    for (size_t i = 0; i < _nbins; ++i) {
      integral += _diffjetshapes[i];
    }
    for (size_t i = 0; i < _nbins; ++i) {
      _diffjetshapes[i] /= integral;
    }

  }


  void JetShape::project(const Event& e) {
    const Jets jets = applyProjection<JetAlg>(e, "Jets").jets();
    calc(jets);
  }


}
