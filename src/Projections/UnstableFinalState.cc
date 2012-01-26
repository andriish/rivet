// -*- C++ -*-
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Cmp.hh"

#define IS_PARTON_PDGID(id) ( abs(id) <= 100 && abs(id) != 22 && (abs(id) < 11 || abs(id) > 18) )

namespace Rivet {

  void UnstableFinalState::project(const Event& e) {
    _theParticles.clear();

    // \todo We should really implement the FinalState algorithm here instead
    double etamin, etamax;
    if ( _etaRanges.empty() ) {
      etamin = -MAXRAPIDITY;
      etamax = MAXRAPIDITY;
    }
    else {
      // with our current constructor choice, we can only ever have one entry
      assert( _etaRanges.size() == 1 );
      etamin = _etaRanges[0].first;
      etamax = _etaRanges[0].second;
    }

    foreach (GenParticle* p, Rivet::particles(e.genEvent())) {
      const int st = p->status();
      bool passed = \
        ( st == 1 || (st == 2 && abs(p->pdg_id()) != 22) ) &&
        !isZero(p->momentum().perp()) && p->momentum().perp() >= _ptmin &&
        p->momentum().eta() > etamin && p->momentum().eta() < etamax &&
        !IS_PARTON_PDGID(p->pdg_id());

      // Avoid double counting by re-marking as unpassed if particle ID == parent ID
      const GenVertex* pv = p->production_vertex();
      const GenVertex* dv = p->end_vertex();
      if (passed && pv) {
        for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin() ;
             pp != pv->particles_in_const_end() ; ++pp) {
          if ( p->pdg_id() == (*pp)->pdg_id() ) {
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
