// -*- C++ -*-
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {


  InvMassFinalState::InvMassFinalState(const FinalState& fsp,
                                       const std::pair<long, long>& idpair, // pair of decay products
                                       double minmass, // min inv mass
                                       double maxmass) // max inv mass
    : _minmass(minmass), _maxmass(maxmass)
  {
    setName("InvMassFinalState");
    addProjection(fsp, "FS");
    _decayids.push_back(idpair);
  }


  InvMassFinalState::InvMassFinalState(const FinalState& fsp,
                                       const std::vector<std::pair<long, long> >& idpairs,  // vector of pairs of decay products
                                       double minmass, // min inv mass
                                       double maxmass) // max inv mass
    : _decayids(idpairs), _minmass(minmass), _maxmass(maxmass)
  { 
    setName("InvMassFinalState");
    addProjection(fsp, "FS");
  }
  

  int InvMassFinalState::compare(const Projection& p) const {
    // First compare the final states we are running on
    int fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != EQUIVALENT) return fscmp;

    // Then compare the two as final states
    const InvMassFinalState & other = dynamic_cast <const InvMassFinalState&>(p);
    fscmp = FinalState::compare(other);
    if (fscmp != EQUIVALENT) return fscmp;

    // Then compare the mass limits
    int massllimcmp = cmp(_minmass, other._minmass);
    if (massllimcmp != EQUIVALENT) return massllimcmp;
    int masshlimcmp = cmp(_maxmass, other._maxmass);
    if (masshlimcmp != EQUIVALENT) return masshlimcmp;

    // Compare the decay species
    int decaycmp = cmp(_decayids, other._decayids);
    if (decaycmp != EQUIVALENT) return decaycmp;

    // Finally compare them as final states 
    return FinalState::compare(other);
  } 
  


  void InvMassFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _theParticles.clear();

    // Containers for the particles of type specified in the pair
    vector<const Particle*> type1;
    vector<const Particle*> type2;
    // Get all the particles of the type specified in the pair from the particle list
    foreach (const Particle& ipart, fs.particles()) {
      // Loop around possible particle pairs
      typedef pair<long,long> longpair;
      foreach (const longpair& ipair, _decayids) {
        if (ipart.pdgId() == ipair.first) {
          if (accept(ipart.genParticle())) {
            type1 += &ipart;
          }
        } else if (ipart.pdgId() == ipair.second) {
          if (accept(ipart.genParticle())) {
            type2 += &ipart;
          }
        }
      }
    }

    // Temporary container of selected particles iterators
    // Useful to compare iterators and avoid double occurrences of the same 
    // particle in case it matches with more than another particle
    vector<const Particle*> tmp;

    // Now calculate the inv mass
    foreach (const Particle* i1, type1) {
      foreach (const Particle* i2, type2) {
        FourMomentum v4 = i1->momentum() + i2->momentum();
        if (v4.mass() > _minmass && v4.mass() < _maxmass) {
          // Avoid duplicates
          if (find(tmp.begin(), tmp.end(), i1) == tmp.end()) {
            tmp.push_back(i1);
            _theParticles.push_back(*i1);
          }
          if (find(tmp.begin(), tmp.end(), i2) == tmp.end()) {
            tmp.push_back(i2);
            _theParticles.push_back(*i2);
          }
          getLog() << Log::DEBUG << "Selecting particles with IDs " 
                   << i1->pdgId() << " & " << i2->pdgId()
                   << " and mass = " << v4.mass()/GeV << " GeV" << endl;
        }
      }
    }
    
    getLog() << Log::DEBUG << "Selected " << _theParticles.size() << " particles." << endl;
    if (getLog().isActive(Log::TRACE)) {
      foreach (const Particle& p, _theParticles) {
        getLog() << Log::TRACE << "ID: " << p.pdgId()
                 << ", barcode: " << p.genParticle().barcode() << endl;
      }
    }
  }
 
 
}
