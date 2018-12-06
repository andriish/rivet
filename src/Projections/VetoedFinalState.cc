// -*- C++ -*-
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  int VetoedFinalState::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != EQUIVALENT) return fscmp;
    /// @todo We can do better than this...
    if (_vetofsnames.size() != 0) return UNDEFINED;
    const VetoedFinalState& other = dynamic_cast<const VetoedFinalState&>(p);
    return \
      cmp(_vetoCodes, other._vetoCodes) ||
      cmp(_compositeVetoes, other._compositeVetoes) ||
      cmp(_parentVetoes, other._parentVetoes);
  }


  void VetoedFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size());

    for (const Particle& p : fs.particles()) {
      if (getLog().isActive(Log::TRACE)) {
        vector<long> codes;
        for (VetoDetails::const_iterator code = _vetoCodes.begin(); code != _vetoCodes.end(); ++code) {
          codes.push_back(code->first);
        }
        const string codestr = "{ " + join(codes) + " }";
        MSG_TRACE(p.pid() << " vs. veto codes = " << codestr << " (" << codes.size() << ")");
      }
      VetoDetails::iterator iter = _vetoCodes.find(p.pid());
      if (iter == _vetoCodes.end()) {
        MSG_TRACE("Storing with PDG code = " << p.pid() << ", pT = " << p.pT());
        _theParticles.push_back(p);
      } else {
        // This particle code is listed as a possible veto... check pT.
        // Make sure that the pT range is sensible:
        BinaryCut ptrange = iter->second;
        assert(ptrange.first <= ptrange.second);
        stringstream rangess;
        if (ptrange.first < numeric_limits<double>::max()) rangess << ptrange.second;
        rangess << " - ";
        if (ptrange.second < numeric_limits<double>::max()) rangess << ptrange.second;
        MSG_TRACE("ID = " << p.pid() << ", pT range = " << rangess.str());
        stringstream debugline;
        debugline << "with PDG code = " << p.pid() << " pT = " << p.pT();
        if (p.pT() < ptrange.first || p.pT() > ptrange.second) {
          MSG_TRACE("Storing " << debugline.str());
          _theParticles.push_back(p);
        } else {
          MSG_TRACE("Vetoing " << debugline.str());
        }
      }
    }

    set<Particles::iterator> toErase;
    for (set<int>::iterator nIt = _nCompositeDecays.begin();
         nIt != _nCompositeDecays.end() && !_theParticles.empty(); ++nIt) {
      map<set<Particles::iterator>, FourMomentum> oldMasses;
      map<set<Particles::iterator>, FourMomentum> newMasses;
      set<Particles::iterator> start;
      start.insert(_theParticles.begin());
      oldMasses.insert(pair<set<Particles::iterator>, FourMomentum>
                       (start, _theParticles.begin()->momentum()));

      for (int nParts = 1; nParts != *nIt; ++nParts) {
        for (map<set<Particles::iterator>, FourMomentum>::iterator mIt = oldMasses.begin();
             mIt != oldMasses.end(); ++mIt) {
          Particles::iterator pStart = *(mIt->first.rbegin());
          for (Particles::iterator pIt = pStart + 1; pIt != _theParticles.end(); ++pIt) {
            FourMomentum cMom = mIt->second + pIt->momentum();
            set<Particles::iterator> pList(mIt->first);
            pList.insert(pIt);
            newMasses[pList] = cMom;
          }
        }
        oldMasses = newMasses;
        newMasses.clear();
      }
      for (map<set<Particles::iterator>, FourMomentum>::iterator mIt = oldMasses.begin();
           mIt != oldMasses.end(); ++mIt) {
        double mass2 = mIt->second.mass2();
        if (mass2 >= 0.0) {
          double mass = sqrt(mass2);
          for (CompositeVeto::iterator cIt = _compositeVetoes.lower_bound(*nIt);
               cIt != _compositeVetoes.upper_bound(*nIt); ++cIt) {
            BinaryCut massRange = cIt->second;
            if (mass < massRange.second && mass > massRange.first) {
              for (set<Particles::iterator>::iterator lIt = mIt->first.begin();
                   lIt != mIt->first.end(); ++lIt) {
                toErase.insert(*lIt);
              }
            }
          }
        }
      }
    }

    for (set<Particles::iterator>::reverse_iterator p = toErase.rbegin(); p != toErase.rend(); ++p) {
      _theParticles.erase(*p);
    }

    // Remove particles whose parents match entries in the parent veto PDG ID codes list
    for (PdgId vetoid : _parentVetoes)
      ifilter_discard(_theParticles, [&](const Particle& p) {
          // Loop over parents and test their IDs
          for (ConstGenParticlePtr parent : Rivet::particles(p.genParticle(), Relatives::ANCESTORS))
            if (parent->pdg_id() == vetoid) return true;
          return false; });

    // Now veto on FS particles from other projections
    for (const string& ifs : _vetofsnames) {
      const Particles& vfsp = applyProjection<FinalState>(e, ifs).particles();
      ifilter_discard(_theParticles, [&](const Particle& p) {
          for (const Particle& pcheck : vfsp) {
            MSG_TRACE("Comparing with veto particle: " << p << " vs. " << pcheck);
            if (p.genParticle() == pcheck.genParticle()) return true;
          }
          return false; });
    }

  }


}
