// -*- C++ -*-
//     Jet algorithm from the Field & Stuart minimum bias study.
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/TrackJet.hh"
#include "Rivet/Cmp.hh"
#include <algorithm>

namespace Rivet {


  int TrackJet::compare(const Projection& p) const {
    const TrackJet& other = dynamic_cast<const TrackJet&>(p);
    return pcmp(_fsproj, other._fsproj);
  }


  void TrackJet::project(const Event & e) {
    Log& log = getLog();
    _jets.clear();

    // Project into final state
    // NB to be true to the original, the final state projection should have 
    // specific cuts on eta, pT and require stable charged particles.
    log << Log::DEBUG << "About to apply the final state projection with eta and pt cuts" << endl;
    const FinalState& fs = e.applyProjection(_fsproj);

    // Put each particle into the collection of tracks and sort (decreasing in pT)
    vector<FourMomentum> tracksvector; // Need to use a vector because you can't use std::sort with lists
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      FourMomentum lv(p->getMomentum());
      tracksvector.push_back(lv);
    }
    // Now sort particles in pT
    sort(tracksvector.begin(), tracksvector.end(), compareVecsByPt);
    // Make it into a list
    Tracks tracks(tracksvector.begin(), tracksvector.end());

    // Now find jets using Field & Stuart criteria
    log << Log::DEBUG << "About to assign tracks into jets" << endl;
    while (!tracks.empty()) {

      Tracks::iterator t = tracks.begin();
      // Get eta and phi for this track
      const double eta = t->vector3().pseudorapidity();
      const double phi = t->vector3().azimuthalAngle();

      // Make a new jet and put this seed track into it
      Jet thisjet;
      thisjet.addParticle(*t);
      tracks.erase(t);

      // Compare with all unassociated tracks with a smaller pT measure
      Tracks::iterator t2 = tracks.begin();
      while (t2 != tracks.end()) {
        log << Log::DEBUG << "Building jet from tracks" << endl;

        // Compute Deta and Dphi, mapping Dphi into [0,pi]
        double Deta = eta - t2->vector3().pseudorapidity();
        double Dphi = phi - t2->vector3().azimuthalAngle();
        if (Dphi > PI) Dphi = 2*PI - Dphi;

        // Add to this jet if eta-phi distance < 0.7 
        if (sqrt(Deta*Deta + Dphi*Dphi) <= 0.7) {
          // Move this particle into the current jet (no extra sorting needed)
          thisjet.addParticle(*t2);
          t2 = tracks.erase(t2);
        } else {
          ++t2;
        }
      }

      _jets.push_back(thisjet);
    }

    // Sort the jets by pT.
    std::sort(_jets.begin(), _jets.end(), compareJetsByPt);

    if (log.isActive(Log::DEBUG)) {
      log << Log::DEBUG << "Number of jets = " << _jets.size() << endl;
      size_t njet = 1;
      for (Jets::const_iterator j = _jets.begin(); j != _jets.end(); ++j)
        log << Log::DEBUG << "Number of tracks in jet #" << njet++ << " = " << j->getNumParticles() << endl;
    }
  }

}
