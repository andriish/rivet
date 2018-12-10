// -*- C++ -*-
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @todo Add a FIRST/LAST/ANY enum to specify the mode for uniquifying replica chains (default = LAST)


  void UnstableFinalState::project(const Event& e) {
    _theParticles.clear();

    /// @todo Replace PID veto list with PID:: functions?
    vector<PdgId> vetoIds;
    vetoIds += 22; // status 2 photons don't count!
    vetoIds += 110; vetoIds += 990; vetoIds += 9990; // Reggeons
    //vetoIds += 9902210; // something weird from PYTHIA6

    for (ConstGenParticlePtr p : Rivet::particles(e.genEvent())) {
      const int st = p->status();
      bool passed =
        (st == 1 || (st == 2 && !contains(vetoIds, abs(p->pdg_id())))) &&
        !PID::isParton(p->pdg_id()) && ///< Always veto partons
        p->status() !=4 && // Filter beam particles
        _cuts->accept(p->momentum());

      // Avoid double counting by re-marking as unpassed if ID == (any) parent ID
      ConstGenVertexPtr pv = p->production_vertex();
      // Avoid double counting by re-marking as unpassed if ID == any child ID
      ConstGenVertexPtr dv = p->end_vertex();
      if (passed && dv) {
        for (ConstGenParticlePtr pp : ::Rivet::particles(dv, Relatives::CHILDREN)) {
          if (p->pdg_id() == pp->pdg_id() && pp->status() == 2) {
            passed = false;
            break;
          }
        }
      }

      // Add to output particles collection
      if (passed) _theParticles.push_back(Particle(p));

      // Log parents and children
      if (getLog().isActive(Log::TRACE)) {
        MSG_TRACE("ID = " << p->pdg_id()
                  << ", status = " << st
                  << ", pT = " << p->momentum().perp()
                  << ", eta = " << p->momentum().eta()
                  << ": result = " << std::boolalpha << passed);
        if (pv) {
          for (ConstGenParticlePtr pp : ::Rivet::particles(pv, Relatives::PARENTS))
            MSG_TRACE("  parent ID = " << pp->pdg_id());
        }
        if (dv) {
          for (ConstGenParticlePtr pp : Rivet::particles(dv, Relatives::CHILDREN))
            MSG_TRACE("  child ID  = " << pp->pdg_id());
        }
      }
    }
    MSG_DEBUG("Number of unstable final-state particles = " << _theParticles.size());
  }


}
