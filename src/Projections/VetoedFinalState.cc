// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Cmp.hh"
#include "Rivet/Tools/Utils.hh"
#include <algorithm>

namespace Rivet {


  int VetoedFinalState::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != PCmp::EQUIVALENT) return fscmp;
    const VetoedFinalState& other = dynamic_cast<const VetoedFinalState&>(p);
    int vfssize = cmp(_vetofsnames.size(), other._vetofsnames.size());
    if (vfssize != PCmp::EQUIVALENT) return vfssize;
    //if the size is the same retrieve the FS projections, store them in two ordered set,
    //compare the sets element by element
    set<const Projection*> vfs;
    set<const Projection*> other_vfs;
    for (FSNames::const_iterator ivfs = _vetofsnames.begin(); ivfs != _vetofsnames.end(); ++ivfs){
      vfs.insert(&(getProjection(*ivfs)));
      other_vfs.insert(&(other.getProjection(*ivfs)));
    }
    int isetcmp = cmp(vfs, other_vfs);
    if (isetcmp != PCmp::EQUIVALENT) return isetcmp;
    return \
      cmp(_vetoCodes, other._vetoCodes) ||
      cmp(_compositeVetoes, other._compositeVetoes) ||
      cmp(_parentVetoes, other._parentVetoes);
  }


  void VetoedFinalState::project(const Event& e) {
    Log log = getLog();
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
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
        log << Log::DEBUG << p->pdgId() << " vs. veto codes = " 
            << codestr << " (" << codes.size() << ")" << endl;
      }
      const long pdgid = p->pdgId();
      VetoDetails::iterator iter = _vetoCodes.find(pdgid);
      if ( (iter == _vetoCodes.end())) {
        log << Log::DEBUG << "Storing with PDG code " << pdgid << " pt " 
            << p->momentum().vector3().polarRadius() << endl;
        _theParticles.push_back(*p);
      } else {
        // This particle code is listed as a possible veto... check pT.
        BinaryCut ptrange = iter->second;
        // Make sure that the pT range is sensible.
        assert(ptrange.getHigherThan() <= ptrange.getLowerThan());
        double pt = p->momentum().vector3().polarRadius();
        stringstream rangess;
        if (ptrange.getHigherThan() < numeric_limits<double>::max()) rangess << ptrange.getHigherThan();
        rangess << " - ";
        if (ptrange.getLowerThan() < numeric_limits<double>::max()) rangess << ptrange.getLowerThan();
        log << Log::DEBUG << "ID = " << pdgid << ", pT range = " << rangess.str();
        stringstream debugline;
        debugline << "with PDG code = " << pdgid << " pT = " << p->momentum().vector3().polarRadius();
        if (pt < ptrange.getHigherThan() || pt > ptrange.getLowerThan()) { 
          log << Log::DEBUG << "Storing " << debugline.str() << endl;
          _theParticles.push_back(*p);
        } else {
          log << Log::DEBUG << "Vetoing " << debugline.str() << endl;
        }
      }
    }
  
    set<ParticleVector::iterator> toErase;
    for (set<int>::iterator nIt = _nCompositeDecays.begin();
         nIt != _nCompositeDecays.end() && !_theParticles.empty(); ++nIt) {
      map<set<ParticleVector::iterator>, FourMomentum> oldMasses;
      map<set<ParticleVector::iterator>, FourMomentum> newMasses;
      set<ParticleVector::iterator> start;
      start.insert(_theParticles.begin());
      oldMasses.insert(pair<set<ParticleVector::iterator>, FourMomentum>
                       (start, _theParticles.begin()->momentum()));

      for (int nParts = 1; nParts != *nIt; ++nParts) {
        for (map<set<ParticleVector::iterator>, FourMomentum>::iterator mIt = oldMasses.begin();
             mIt != oldMasses.end(); ++mIt) {
          ParticleVector::iterator pStart = *(mIt->first.rbegin());
          for (ParticleVector::iterator pIt = pStart + 1; pIt != _theParticles.end(); ++pIt) { 
            FourMomentum cMom = mIt->second + pIt->momentum();
            set<ParticleVector::iterator> pList(mIt->first);
            pList.insert(pIt);   
            newMasses[pList] = cMom;
          }
	}
        oldMasses = newMasses;
        newMasses.clear();
      }
      for (map<set<ParticleVector::iterator>, FourMomentum>::iterator mIt = oldMasses.begin();
           mIt != oldMasses.end(); ++mIt) {
        double mass2 = mIt->second.mass2();
        if (mass2 >= 0.0) {
          double mass = sqrt(mass2);
	  for (CompositeVeto::iterator cIt = _compositeVetoes.lower_bound(*nIt);
	       cIt != _compositeVetoes.upper_bound(*nIt); ++cIt) {
	    BinaryCut massRange = cIt->second;
            if (mass < massRange.getLowerThan() && mass > massRange.getHigherThan()) {
              for (set<ParticleVector::iterator>::iterator lIt = mIt->first.begin();
                   lIt != mIt->first.end(); ++lIt) {
		toErase.insert(*lIt);                
              }
	     }
	  }
        }
      }
    }
    
    for(set<ParticleVector::iterator>::reverse_iterator p = toErase.rbegin();
	p != toErase.rend(); ++p){
      _theParticles.erase(*p);
    }
    

    for (ParentVetos::const_iterator vIt = _parentVetoes.begin(); vIt != _parentVetoes.end(); ++vIt) {
      for (ParticleVector::iterator p = _theParticles.begin(); p != _theParticles.end(); ++p) {
        GenVertex *startVtx=((*p).getHepMCParticle()).production_vertex();
        bool veto = false;
        GenParticle HepMCP = (*p).getHepMCParticle();
        if (startVtx!=0) {
          for (GenVertex::particle_iterator pIt = startVtx->particles_begin(HepMC::ancestors);
               pIt != startVtx->particles_end(HepMC::ancestors) && !veto; ++pIt) {
            
            if (*vIt == (*pIt)->pdg_id()) {
              veto = true;
              p = _theParticles.erase(p);
              --p;
            }
          }
        }
      }
    }

    // Now veto on the FS
    for (FSNames::const_iterator ivfs = _vetofsnames.begin(); ivfs != _vetofsnames.end(); ++ivfs){
      const FinalState& vfs = applyProjection<FinalState>(e, *ivfs);
      const ParticleVector& vfsp = vfs.particles();
      for (ParticleVector::iterator icheck = _theParticles.begin(); icheck != _theParticles.end(); ++icheck){
        if (!icheck->hasHepMCParticle()) continue;
        bool found = false;
        for (ParticleVector::const_iterator ipart = vfsp.begin(); ipart != vfsp.end(); ++ipart){
          if (!ipart->hasHepMCParticle()) continue;
          log << Log::DEBUG << "comparing barcode " << icheck->getHepMCParticle().barcode() << " with veto particle " << ipart->getHepMCParticle().barcode() << endl; 
          if (ipart->getHepMCParticle().barcode() == icheck->getHepMCParticle().barcode()){
            found = true;
            break;
          }
        }
        if (found) {
          _theParticles.erase(icheck);
          --icheck;
        }	
      }	
    }
  }
  

}
