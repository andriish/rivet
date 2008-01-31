// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Cmp.hh"
#include "Rivet/Tools/Utils.hh"
#include <algorithm>

namespace Rivet {

  int VetoedFinalState::compare(const Projection& p) const {
    const VetoedFinalState& other = dynamic_cast<const VetoedFinalState&>(p);
    const int fscmp = pcmp(_fsproj, other._fsproj);
    if (fscmp != 0) return fscmp;
    if (_vetoCodes == other._vetoCodes && 
        _compositeVetoes == other._compositeVetoes && 
        _parentVetoes == other._parentVetoes) return 0;
    int ret = 1;
    if (_vetoCodes < other._vetoCodes) ret = -1;
    if (_compositeVetoes < other._compositeVetoes) ret *= -1;
    if (_parentVetoes < other._parentVetoes) ret *=-1;
    return ret;
  }


  void VetoedFinalState::project(const Event& e) {
    Log log = getLog();
    const FinalState& fs = e.applyProjection(_fsproj);
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size());
    const ParticleVector& fsps = fs.particles();
    for (ParticleVector::const_iterator p = fsps.begin(); p != fsps.end(); ++p) {
      if (log.isActive(Log::DEBUG)) {
        vector<long> codes; 
        for (VetoDetails::const_iterator code = _vetoCodes.begin(); code != _vetoCodes.end(); ++code) {
          codes.push_back(code->first);
        }
        const string codestr = "{ " + join(codes) + " }";
        log << Log::DEBUG << p->getPdgId() << " vs. veto codes = " 
            << codestr << " (" << codes.size() << ")" << endl;
      }
      const long pdgid = p->getPdgId();
      VetoDetails::iterator iter = _vetoCodes.find(pdgid);
      if ( (iter == _vetoCodes.end())) {
        log << Log::DEBUG << "Storing with PDG code " << pdgid << " pt " 
            << p->getMomentum().vector3().polarRadius() << endl;
        _theParticles.push_back(*p);
      } else {
        // This particle code is listed as a possible veto... check pT.
        BinaryCut ptrange = iter->second;
        // Make sure that the pT range is sensible.
        assert(ptrange.getHigherThan() <= ptrange.getLowerThan());
        double pt = p->getMomentum().vector3().polarRadius();
        stringstream rangess;
        if (ptrange.getHigherThan() < numeric_limits<double>::max()) rangess << ptrange.getHigherThan();
        rangess << " - ";
        if (ptrange.getLowerThan() < numeric_limits<double>::max()) rangess << ptrange.getLowerThan();
        log << Log::DEBUG << "ID = " << pdgid << ", pT range = " << rangess.str();
        stringstream debugline;
        debugline << "with PDG code = " << pdgid << " pT = " << p->getMomentum().vector3().polarRadius();
        if (pt < ptrange.getHigherThan() || pt > ptrange.getLowerThan()) { 
          log << Log::DEBUG << "Storing " << debugline << endl;
          _theParticles.push_back(*p);
        } else {
          log << Log::DEBUG << "Vetoing " << debugline << endl;
        }
      }
    }
  
    
    for(set<int>::iterator nIt = _nCompositeDecays.begin();
        nIt != _nCompositeDecays.end()&&_theParticles.size()!=0; ++nIt){
      map<set<ParticleVector::iterator>, FourMomentum> oldMasses;
      map<set<ParticleVector::iterator>, FourMomentum> newMasses;
      set<ParticleVector::iterator> start;
      start.insert(_theParticles.begin());
      oldMasses.insert(pair<set<ParticleVector::iterator>, FourMomentum>
       (start, _theParticles.begin()->getMomentum()));

      for(int nParts=1; nParts != *nIt; ++nParts){

        for(map<set<ParticleVector::iterator>, FourMomentum>::iterator mIt = oldMasses.begin();
            mIt != oldMasses.end(); ++mIt){
          ParticleVector::iterator pStart = *(mIt->first.rbegin());
          for(ParticleVector::iterator pIt = pStart + 1; 
              pIt != _theParticles.end(); ++pIt){
            //Particle p = (*pIt);
//            FourMomentum momIncr = pIt->getMomentum();
            //double tmp = momIncr.p().mod2();
  
            FourMomentum cMom = mIt->second + pIt->getMomentum();//momIncr;
            set<ParticleVector::iterator> pList(mIt->first);
            pList.insert(pIt);
            newMasses[pList] = cMom;
 
          }
        }
        oldMasses = newMasses;
        newMasses.clear();
      }

      for(map<set<ParticleVector::iterator>, FourMomentum>::iterator mIt = oldMasses.begin();
          mIt != oldMasses.end(); ++mIt){
        double mass2 = mIt->second.mass2();
        if(mass2 >= 0.0){
          double mass = sqrt(mass2);
          for(CompositeVeto::iterator cIt = _compositeVetoes.lower_bound(*nIt);
              cIt != _compositeVetoes.upper_bound(*nIt); ++cIt){
            BinaryCut massRange = cIt->second;

            if(mass < massRange.getLowerThan() && mass > massRange.getHigherThan()){
              for(set<ParticleVector::iterator>::iterator lIt = mIt->first.begin();
                  lIt != mIt->first.end(); ++lIt){
                _theParticles.erase(*lIt);
              }
            }
          }
        }
      }
    }
    
    for(vector<long>::iterator vIt = _parentVetoes.begin();
        vIt != _parentVetoes.end(); ++vIt){
      for(ParticleVector::iterator p = _theParticles.begin();
          p != _theParticles.end(); ++p){
        GenVertex *startVtx=((*p).getHepMCParticle()).production_vertex();
        bool veto = false;
        GenParticle HepMCP = (*p).getHepMCParticle();
        if(startVtx!=0){
          for(GenVertex::particle_iterator pIt = startVtx->particles_begin(HepMC::ancestors);
              pIt != startVtx->particles_end(HepMC::ancestors) && !veto; ++pIt){

            if(*vIt == (*pIt)->pdg_id()){
              veto = true;
              p = _theParticles.erase(p);
              --p;
            }
          }
        }
      }
    }
    return;
  } 
}
