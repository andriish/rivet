// -*- C++ -*-
#include "Rivet/Projections/ParticleFinder.hh"

namespace Rivet {


  int ParticleFinder::compare(const Projection& p) const {
    const ParticleFinder& other = dynamic_cast<const ParticleFinder&>(p);

    //MSG_TRACE("PF::compare: " << 1 << " " << this << " " << &p);
    vector<pair<double, double> > eta1(_etaRanges);
    vector<pair<double, double> > eta2(other._etaRanges);
    sort(eta1.begin(), eta1.end());
    sort(eta2.begin(), eta2.end());

    //MSG_TRACE("PF::compare: " << 2 << " " << this << " " << &p);
    if (eta1 < eta2) return ORDERED;
    else if (eta2 < eta1) return ANTIORDERED;

    //MSG_TRACE("PF::compare: " << 3 << " " << this << " " << &p);
    return cmp(_ptmin, other._ptmin);
  }


  // void ParticleFinder::project(const Event& e) {
  //   _theParticles.clear();
  // }


}
