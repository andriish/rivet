// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/SVertex.hh"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenEvent.h"
#include "HepMC/SimpleVector.h"
#include "Rivet/Cmp.hh"

namespace Rivet {


  int SVertex::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "PV");
    if (fscmp != PCmp::EQUIVALENT) return fscmp;
    const SVertex& other = pcast<SVertex>(p);
    return \
      cmp(_detEta, other._detEta) ||
      cmp(_IPres, other._IPres) ||
      cmp(_DLS, other._DLS) ||
      cmp(_DLSres, _DLSres);
  }



  void SVertex::project(const Event& e) {
    const PVertex& pvtx = applyProjection<PVertex>(e, "PV");
    const Vector3 pvpos = pvtx.position();
    const ChargedFinalState& chfs = applyProjection<ChargedFinalState>(e, "FS");
    
    // Produce vector of vertices, each containing a vector of all charged 
    // final state particles belonging to this vertex
    typedef map<GenVertex*,ParticleVector> VtxPartsMap;
    VtxPartsMap vtxparts;
    foreach (const Particle& p, chfs.particles()) {
      // Consider only charged particles in tracker geometrical acceptance
      /// @todo Use acceptance from the FinalState instead
      if (fabs(p.momentum().pseudorapidity()) > _detEta) continue;
      HepMC::GenVertex* pvtx = p.genParticle().production_vertex();
      vtxparts[pvtx].push_back(p);
    }
  
    // Check if jets are tagged, by means of selected vertices fulfilling track criteria
    _taggedjets.clear();
    for (VtxPartsMap::const_iterator vp = vtxparts.begin(); vp != vtxparts.end(); ++vp) {
      FourMomentum vtxVisMom;
      if (! _applyVtxTrackCuts(vp->second, pvpos, vtxVisMom) ) break;
      for (vector<FourMomentum>::const_iterator jaxis = _jetaxes.begin(); jaxis != _jetaxes.end(); ++jaxis) {
        // Delta{R} requirement between jet and visible vector sum of vertex tracks
        if (deltaR(*jaxis, vtxVisMom) < _deltaR) {
          const double dls = get2dDecayLength(vp->first->position(), pvpos, *jaxis) / _DLSres;
          if (dls > _DLS) _taggedjets.push_back(*jaxis);
        }
      }
    }
    
  }
  
  
  
  /// Analysis dependent cuts on vertex tracks in SVertex projection 
  /// Since the analysis specific cuts are very complex, they are not 
  /// implemented in the projection and are instead passed via a function (object).
  /// SVertex member function implementation below 
  /// in: reference to instance of SVertex projection, ParticleVector of
  ///     vertex to be analyzed, primary (Gen)Vertex
  /// out: FourMomentum = visible Momentum of vertex (selected tracks), 
  /// return bool: cuts passed? 1 : 0 
  /// @todo Move this into the projection concrete class.
  bool SVertex::_applyVtxTrackCuts(const ParticleVector& vtxparts, 
                                   const Vector3& pvtxpos, 
                                   FourMomentum vtxVisMom) 
  {
    // Check vertex final state charged particles, if fulfilling track criteria
    size_t pass1trk1pTdcaSig25(0), pass1trk05pTdcaSig25(0), 
      pass2trk15pTdcaSig3(0), pass2trk1pTdcaSig3(0);
    
    foreach (const Particle& vp, vtxparts) {
      const double IPsig = get2dClosestApproach(vp.genParticle(), pvtxpos) / _IPres;
      
      // Update "visible momentum" vector (returned by reference).
      if (vp.momentum().pT() > 0.5) {
        vtxVisMom += vp.momentum();
      }
      // 1st pass
      if (vtxparts.size() >= 3 && IPsig > 2.5) {
        if (vp.momentum().pT() > 1.0) pass1trk1pTdcaSig25++;
        else if (vp.momentum().pT() > 0.5) pass1trk05pTdcaSig25++;
      }
      // 2nd pass
      if (vtxparts.size() >= 2 && IPsig > 3.) {
        if (vp.momentum().pT() > 1.5) pass2trk15pTdcaSig3++;
        else if (vp.momentum().pT() > 1.0) pass2trk1pTdcaSig3++;
      } 
    }

    // Combine info from passes to make yes/no decision about whether this is significant:
    if (pass1trk1pTdcaSig25 >= 1 && pass1trk1pTdcaSig25 + pass1trk05pTdcaSig25>=3) return true;
    if (pass2trk15pTdcaSig3 >= 1 && pass2trk15pTdcaSig3 + pass2trk1pTdcaSig3>=2) return true;
    return false;
  }

  
 
}
