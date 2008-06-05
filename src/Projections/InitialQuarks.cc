// -*- C++ -*-
#include "Rivet/Projections/InitialQuarks.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Cmp.hh"

#define IS_PARTON_PDGID(id) ( abs(id) <= 100 && abs(id) != 22 && (abs(id) < 11 || abs(id) > 18) )

namespace Rivet {

  int InitialQuarks::compare(const Projection& p) const {
    const InitialQuarks& other = dynamic_cast<const InitialQuarks&>(p);
    return true;
  }


  void InitialQuarks::project(const Event& e) {
    Log& log = getLog();
    _theParticles.clear();

    for (GenEvent::particle_const_iterator p = e.genEvent().particles_begin();
         p != e.genEvent().particles_end(); ++p) {
      const GenVertex* pv = (*p)->production_vertex();
      const GenVertex* dv = (*p)->end_vertex();
      bool passed = abs((*p)->pdg_id()) >= 1 && abs((*p)->pdg_id()) <= 5;
      if (passed) {
        if (pv!=NULL) {
          for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin() ;
              pp != pv->particles_in_const_end() ; ++pp) {
            // Only accept if parent is electron or Z0
            if ( abs((*pp)->pdg_id()) != 11 && abs((*pp)->pdg_id()) != 23 )
              passed = false;
          }
        }
        else {
          passed = false;
        }
      }

      if (log.isActive(Log::TRACE)) {
        const int st = (*p)->status();
        const double pT = (*p)->momentum().perp();
        const double eta = (*p)->momentum().eta();
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
    log << Log::DEBUG << "Number of initial quarks = " 
        << _theParticles.size() << endl;
    if (not _theParticles.empty())
      for (size_t i=0 ; i<_theParticles.size() ; i++)
        log << Log::DEBUG << "Initial quark[" << i << "] = " 
            << _theParticles[i].getPdgId() << std::endl;
  }

}
