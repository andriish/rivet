// -*- C++ -*-
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  // @deprecated, keep for backwards compatibility for now.
  FinalState::FinalState(double mineta, double maxeta, double minpt)
  {
    setName("FinalState");
    const bool openpt = isZero(minpt);
    const bool openeta = (mineta <= -MAXDOUBLE && maxeta >= MAXDOUBLE);
    MSG_TRACE("Check for open FS conditions:" << std::boolalpha
              << " eta=" << openeta
              << ", pt=" << openpt);
    if ( openpt && openeta ) {
      _cuts = Cuts::open();
    }
    else {
      addProjection(FinalState(), "OpenFS");
      if ( openeta )
	_cuts = Cuts::pT >= minpt;
      else if ( openpt )
	_cuts = EtaIn(mineta,maxeta);
      else
	_cuts = EtaIn(mineta,maxeta) & (Cuts::pT >= minpt);
    }
  }


  FinalState::FinalState(Cut c)
    : ParticleFinder(c)
  {
    setName("FinalState");
    const bool open = ( c == Cuts::open() );
    MSG_TRACE("Check for open FS conditions: " << std::boolalpha << open);
    if ( !open ) {
      addProjection(FinalState(), "OpenFS");
    }
  }




  /// @todo HOW DO WE COMPARE CUTS OBJECTS?
  int FinalState::compare(const Projection& p) const {
    const FinalState& other = dynamic_cast<const FinalState&>(p);

    //MSG_TRACE("FS::compare: " << 1 << " " << this << " " << &p);
    //    std::vector<std::pair<double, double> > eta1(_etaRanges);
    //std::vector<std::pair<double, double> > eta2(other._etaRanges);
    //std::sort(eta1.begin(), eta1.end());
    //std::sort(eta2.begin(), eta2.end());

    //MSG_TRACE("FS::compare: " << 2 << " " << this << " " << &p);
    //if (eta1 < eta2) return ORDERED;
    //else if (eta2 < eta1) return UNORDERED;

    //MSG_TRACE("FS::compare: " << 3 << " " << this << " " << &p);
    //return cmp(_ptmin, other._ptmin);
    return _cuts == other._cuts ? EQUIVALENT : UNDEFINED;
  }



  void FinalState::project(const Event& e) {
    _theParticles.clear();

    // Handle "open FS" special case
    if ( _cuts == Cuts::open() ) {
      //MSG_TRACE("Open FS processing: should only see this once per event ("
      //           << e.genEvent().event_number() << ")");
      foreach (const GenParticle* p, Rivet::particles(e.genEvent())) {
        if (p->status() == 1) {
          //MSG_TRACE("FS GV = " << p->production_vertex());
          _theParticles.push_back(Particle(*p));
        }
      }
      return;
    }

    // If this is not itself the "open" FS, base the calculations on the open FS' results
    /// @todo In general, we'd like to calculate a restrictive FS based on the most restricted superset FS.
    const Particles allstable = applyProjection<FinalState>(e, "OpenFS").particles();
    foreach (const Particle& p, allstable) {
      const bool passed = accept(p);
      MSG_TRACE("Choosing: ID = " << p.pdgId()
                << ", pT = " << p.pT()
                << ", eta = " << p.eta()
                << ": result = " << std::boolalpha << passed);
      if (passed) _theParticles.push_back(p);
    }
    //MSG_DEBUG("Number of final-state particles = " << _theParticles.size());
  }


  /// Decide if a particle is to be accepted or not.
  bool FinalState::accept(const Particle& p) const {
    // Not having s.c. == 1 should never happen!
    assert(p.genParticle() == NULL || p.genParticle()->status() == 1);

    return _cuts->accept(p);
  }


}
