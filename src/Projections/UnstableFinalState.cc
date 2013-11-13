// -*- C++ -*-
#include "Rivet/Projections/UnstableFinalState.hh"

#define IS_PARTON_PDGID(id) ( abs(id) <= 100 && abs(id) != 22 && (abs(id) < 11 || abs(id) > 18) )

namespace Rivet {


  void UnstableFinalState::project(const Event& e) {
    _theParticles.clear();

    vector<PdgId> vetoIds;
    vetoIds += 22; // status 2 photons don't count!
    vetoIds += 110; vetoIds += 990; vetoIds += 9990; // Reggeons
    //vetoIds += 9902210; // something weird from PYTHIA6
    foreach (GenParticle* p, Rivet::particles(e.genEvent())) {
      const int st = p->status();
      bool passed =
        (st == 1 || (st == 2 && find(vetoIds.begin(), vetoIds.end(), abs(p->pdg_id())) == vetoIds.end())) &&
        !IS_PARTON_PDGID(p->pdg_id()) && //< Always veto partons?
        !isZero(p->momentum().perp()) && 
	_cuts->accept(p->momentum());

      // Avoid double counting by re-marking as unpassed if particle ID == parent ID
      const GenVertex* pv = p->production_vertex();
      const GenVertex* dv = p->end_vertex();
      if (passed && pv) {
        for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin() ;
             pp != pv->particles_in_const_end() ; ++pp) {
          if ( p->pdg_id() == (*pp)->pdg_id() && (*pp)->status() == 2 ) {
            passed = false;
            break;
          }
        }
      }

      // Add to output particles collection
      if (passed) {
        _theParticles.push_back(Particle(*p));
      }

      // Log parents and children
      if (getLog().isActive(Log::TRACE)) {
        MSG_TRACE("ID = " << p->pdg_id()
                  << ", status = " << st
                  << ", pT = " << p->momentum().perp()
                  << ", eta = " << p->momentum().eta()
                  << ": result = " << std::boolalpha << passed);
        if (pv) {
          for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin() ;
               pp != pv->particles_in_const_end() ; ++pp) {
            MSG_TRACE("  parent ID = " << (*pp)->pdg_id());
          }
        }
        if (dv) {
          for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin() ;
               pp != dv->particles_out_const_end() ; ++pp) {
            MSG_TRACE("  child ID  = " << (*pp)->pdg_id());
          }
        }
      }
    }
    MSG_DEBUG("Number of unstable final-state particles = " << _theParticles.size());
  }


}
