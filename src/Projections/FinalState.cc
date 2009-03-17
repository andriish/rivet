// -*- C++ -*-
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {

  int FinalState::compare(const Projection& p) const {
    const FinalState& other = dynamic_cast<const FinalState&>(p);

    std::vector<std::pair<double, double> > eta1(_etaRanges);
    std::vector<std::pair<double, double> > eta2(other._etaRanges);
    std::sort(eta1.begin(), eta1.end());
    std::sort(eta2.begin(), eta2.end());

    if (eta1 < eta2) return PCmp::ORDERED;
    else if (eta2 < eta1) return PCmp::UNORDERED;

    return cmp(_ptmin, other._ptmin);
  }


  void FinalState::project(const Event& e) {
    Log& log = getLog();
    _theParticles.clear();

    for (GenEvent::particle_const_iterator p = e.genEvent().particles_begin();
         p != e.genEvent().particles_end(); ++p) {
      // Only include particles which are final state (status = 1) and which
      // pass the eta and phi cuts. The eta cut is pre-tested by checking if the
      // x and y components of the momentum are non-zero since the vectors might
      // throw an exception otherwise.
      const bool passed = accept(**p);  //GIULIO
      if (log.isActive(Log::TRACE)) {
        log << Log::TRACE << std::boolalpha 
            << "ID = " << (*p)->pdg_id() << ", status = " << (*p)->status() << ", pT = " << (*p)->momentum().perp() 
            << ", eta = " << (*p)->momentum().eta() << ": result = " << passed << endl;
      }
      if (passed) _theParticles.push_back(Particle(**p));
    }
    log << Log::DEBUG << "Number of final-state particles = " 
        << _theParticles.size() << endl;
  }

  /// Decide if a particle is to be accepted or not.
  bool FinalState::accept(const GenParticle& p) const {
    const int st = p.status();
    if (st!=1) return false;
    
    if (_ptmin>0.0) {
      const double pT = p.momentum().perp();
      if (pT<_ptmin) return false;
    }
    
    if (_etaRanges.size()>0) {
      bool eta_pass=false;
      const double eta = p.momentum().eta();
      for (size_t i=0; i<_etaRanges.size(); ++i) {
        if (eta>_etaRanges[i].first && eta<_etaRanges[i].second) {
          eta_pass=true;
          break;
        }
      }
      if (!eta_pass) return false;
    }
    
    return true;
  }

}
