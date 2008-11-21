// -*- C++ -*-
// Jet algorithm from the CDF Field & Stuart minimum bias study.

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/TrackJet.hh"
#include "Rivet/Cmp.hh"
#include <algorithm>

namespace Rivet {


  // Helper functions for sorting
  // For sorting Jet objects by pT.
  inline bool compareJetsByPt(const Rivet::Jet& first, const Rivet::Jet& second) {
    return first.ptSum() > second.ptSum();
  }
  

  // For sorting Lorentz four-vectors by pT.
  inline bool compareVecsByPt(const FourMomentum& first, const FourMomentum& second) {
    return pT2(first) > pT2(second);
  }


  /////////////////////////////////


  int TrackJet::compare(const Projection& p) const {
    return mkNamedPCmp(p, "FS");
  }


  void TrackJet::project(const Event& e) {
    _jets.clear();

    // Project into final state
    // NB to be true to the original, the final state projection should have 
    // specific cuts on eta, pT and require stable charged particles.
    getLog() << Log::DEBUG << "About to apply the final state projection with eta and pt cuts" << endl;
    const FinalState& fs = applyProjection<FinalState>(e, "FS");

    // Put each particle into the collection of tracks and sort (decreasing in pT)
    vector<FourMomentum> tracksvector; // Need to use a vector because you can't use std::sort with lists
    foreach (const Particle& p, fs.particles()) {
      const FourMomentum lv = p.momentum();
      tracksvector.push_back(lv);
    }
    // Now sort particles in pT
    sort(tracksvector.begin(), tracksvector.end(), compareVecsByPt);
    // Make it into a list
    Tracks tracks(tracksvector.begin(), tracksvector.end());

    // Now find jets using Field & Stuart criteria
    if (! tracks.empty()) {
      getLog() << Log::DEBUG << "About to assign tracks into jets" << endl;
    }
    while (! tracks.empty()) {

      // Make a new jet based on the highest pT independent track remaining.
      Jet thisjet;
      Tracks::iterator t = tracks.begin();
      thisjet.addParticle(*t);
      tracks.erase(t);

      // Compare with all unassociated tracks with a smaller pT measure
      Tracks::iterator t2 = tracks.begin();
      while (t2 != tracks.end()) {
        getLog() << Log::TRACE << "Building jet from tracks" << endl;

        // Get eta and phi for this jet
        const double jeteta = thisjet.ptWeightedEta();
        const double jetphi = thisjet.ptWeightedPhi();

        // Compute D(eta) and D(phi), mapping Dphi into [0,pi]
        const double Deta = fabs(jeteta - t2->pseudorapidity());
        double tmpDphi = fabs(jetphi - t2->azimuthalAngle());
        const double Dphi = (tmpDphi > PI) ? fabs( 2*PI - tmpDphi ) : tmpDphi;
        assert(Dphi >= 0 && Dphi <= PI);

        // Add to this jet if eta-phi distance < Rmax cutoff
        if (sqrt(Deta*Deta + Dphi*Dphi) < _Rmax) {
          // Move this particle into the current jet (no extra sorting needed)
          thisjet.addParticle(*t2);
          tracks.erase(t2);
          t2 = tracks.begin(); //< This was important!
        } else {
          ++t2;
        }
      }

      _jets.push_back(thisjet);
    }

    // Sort the jets by pT.
    std::sort(_jets.begin(), _jets.end(), compareJetsByPt);

    if (getLog().isActive(Log::DEBUG)) {
      getLog() << Log::DEBUG << "Number of jets = " << _jets.size() << endl;
      size_t njet = 1;
      foreach (const Jet& j, _jets) {
        getLog() << Log::DEBUG << "Number of tracks in jet #" << njet++ << " = " << j.size() << endl;
      }
    }
  }

}
