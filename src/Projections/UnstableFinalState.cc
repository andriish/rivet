// -*- C++ -*-
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Cmp.hh"

#define IS_PARTON_PDGID(id) ( abs(id) <= 100 && abs(id) != 22 && (abs(id) < 11 || abs(id) > 18) )

namespace Rivet {

  int UnstableFinalState::compare(const Projection& p) const {
    const UnstableFinalState& other = dynamic_cast<const UnstableFinalState&>(p);
    return \
      cmp(_etamin, other._etamin) || 
      cmp(_etamax, other._etamax) || 
      cmp(_ptmin, other._ptmin);
  }


  void UnstableFinalState::project(const Event& e) {
    Log& log = getLog();
    _theParticles.clear();

    for (GenEvent::particle_const_iterator p = e.genEvent().particles_begin();
         p != e.genEvent().particles_end(); ++p) {
      // Only include particles which are final state (status = 1) and which
      // pass the eta and phi cuts. The eta cut is pre-tested by checking if the
      // x and y components of the momentum are non-zero since the vectors might
      // throw an exception otherwise.
      const int st = (*p)->status();
      const double pT = (*p)->momentum().perp();
      const double eta = (*p)->momentum().eta();
      const GenVertex* pv = (*p)->production_vertex();
      const GenVertex* dv = (*p)->end_vertex();
      bool passed = ( st == 1 || (st == 2 && abs((*p)->pdg_id()) != 22) ) &&
                    !isZero(pT) &&
                    pT >= _ptmin &&
                    eta > _etamin &&
                    eta < _etamax &&
                    !IS_PARTON_PDGID((*p)->pdg_id());
//      if (passed) {
//        if (pv!=NULL) {
//          for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin() ;
//              pp != pv->particles_in_const_end() ; ++pp) {
//                passed = passed && IS_PARTON_PDGID((*pp)->pdg_id());
//          }
//        }
//      }

      if (log.isActive(Log::TRACE)) {
        log << Log::TRACE << std::boolalpha 
            << "ID = " << (*p)->pdg_id() << ", status = " << st << ", pT = " << pT 
            << ", eta = " << eta << ": result = " << passed << endl;
        if (pv!=NULL) {
          for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin() ;
              pp != pv->particles_in_const_end() ; ++pp) {
            log << Log::TRACE << std::boolalpha
                << "     parent ID = " << (*pp)->pdg_id() << endl;
          }
        }
        if (dv!=NULL) {
          for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin() ;
              pp != dv->particles_out_const_end() ; ++pp) {
            log << Log::TRACE << std::boolalpha
                << "     child ID  = " << (*pp)->pdg_id() << endl;
          }
        }
      }
      if (passed) _theParticles.push_back(Particle(**p));
    }
    log << Log::DEBUG << "Number of final-state particles = " 
        << _theParticles.size() << endl;
  }

}
