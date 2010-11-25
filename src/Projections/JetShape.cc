// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  // Constructor.
  JetShape::JetShape(const JetAlg& jetalg,
                     double rmin, double rmax, size_t nbins,
                     double ptmin, double ptmax,
                     double absrapmin, double absrapmax,
                     RapScheme rapscheme)
    : _rapscheme(rapscheme)
  {
    setName("JetShape");
    _binedges = linspace(rmin, rmax, nbins);
    _ptcuts = make_pair(ptmin, ptmax);
    _rapcuts = make_pair(absrapmin, absrapmax);
    addProjection(jetalg, "Jets");
  }


  // Constructor.
  JetShape::JetShape(const JetAlg& jetalg,
                     vector<double> binedges,
                     double ptmin, double ptmax,
                     double absrapmin, double absrapmax,
                     RapScheme rapscheme)
    : _binedges(binedges), _rapscheme(rapscheme)
  {
    setName("JetShape");
    _ptcuts = make_pair(ptmin, ptmax);
    _rapcuts = make_pair(absrapmin, absrapmax);
    addProjection(jetalg, "Jets");
  }


  int JetShape::compare(const Projection& p) const {
    const int jcmp = mkNamedPCmp(p, "Jets");
    if (jcmp != EQUIVALENT) return jcmp;
    const JetShape& other = pcast<JetShape>(p);
    const int ptcmp = cmp(ptMin(), other.ptMin()) || cmp(ptMax(), other.ptMax());
    if (ptcmp != EQUIVALENT) return ptcmp;
    int bincmp = cmp(numBins(), other.numBins());
    if (bincmp != EQUIVALENT) return bincmp;
    for (size_t i = 0; i < _binedges.size(); ++i) {
      bincmp = cmp(_binedges[i], other._binedges[i]);
      if (bincmp != EQUIVALENT) return bincmp;
    }
    return EQUIVALENT;
  }


  void JetShape::clear() {
    _diffjetshapes = vector<double>(numBins(), 0.0);
  }


  void JetShape::calc(const Jets& jets) {
    clear();

    foreach (const Jet& j, jets) {
      FourMomentum pj = j.momentum();
      if (!inRange(pj.pT(), _ptcuts)) continue;
      /// @todo Introduce a better (i.e. more safe and general) eta/y selection mechanism: MomentumFilter
      if (_rapscheme == PSEUDORAPIDITY && !inRange(fabs(pj.eta()), _rapcuts)) continue;
      if (_rapscheme == RAPIDITY && !inRange(fabs(pj.rapidity()), _rapcuts)) continue;
      foreach (const Particle& p, j.particles()) {
        const double dR = deltaR(pj, p.momentum(), _rapscheme);
        if (!inRange(dR, _binedges.front(), _binedges.back())) continue; //< Out of histo range
        size_t dRindex = -1;
        for (size_t i = 1; i < _binedges.size(); ++i) {
          if (dR < _binedges[i]) {
            dRindex = i-1;
            break;
          }
        }
        assert(inRange(dRindex, 0, numBins()));
        _diffjetshapes[dRindex] += p.momentum().pT();
      }
    }

    // Normalize to total pT
    double integral = 0.0;
    for (size_t i = 0; i < numBins(); ++i) {
      integral += _diffjetshapes[i];
    }
    for (size_t i = 0; i < numBins(); ++i) {
      _diffjetshapes[i] /= integral;
    }

  }


  void JetShape::project(const Event& e) {
    const Jets jets = applyProjection<JetAlg>(e, "Jets").jets(_ptcuts.first, _ptcuts.second,
                                                              -_rapcuts.second, _rapcuts.second, _rapscheme);
    calc(jets);
  }


}
