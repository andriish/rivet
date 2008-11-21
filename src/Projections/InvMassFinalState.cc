// -*- C++ -*-
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {


  int InvMassFinalState::compare(const Projection& p) const {
    // First compare the final states we are running on
    int fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != PCmp::EQUIVALENT) return fscmp;

    // Then compare the two as final states
    const InvMassFinalState & other = dynamic_cast <const InvMassFinalState&>(p);
    fscmp = FinalState::compare(other);
    if (fscmp != PCmp::EQUIVALENT) return fscmp;

    // Then compare the mass limits
    int massllimcmp = cmp(_minmass, other._minmass);
    if (massllimcmp != PCmp::EQUIVALENT) return massllimcmp;
    int masshlimcmp = cmp(_maxmass, other._maxmass);
    if (masshlimcmp != PCmp::EQUIVALENT) return masshlimcmp;

    // Compare the decay species
    int decaycmp = cmp(_decayids, other._decayids);
    if (decaycmp != PCmp::EQUIVALENT) return decaycmp;

    // Finally compare them as final states 
    return FinalState::compare(other);
  } 
  


  void InvMassFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _theParticles.clear();

    // Containers for the particles of type specified in the pair
    vector <ParticleVector::const_iterator> type1;
    vector <ParticleVector::const_iterator> type2;
    // Get all the particles of the type specified in the pair from the particle list
    const ParticleVector& allparticles = fs.particles();
    for (ParticleVector::const_iterator ipart = allparticles.begin(); ipart != allparticles.end(); ++ipart) {
      // Loop around possible particle pairs
      for (vector<pair<long,long> >::const_iterator ipair = _decayids.begin(); ipair != _decayids.end() ; ++ipair) {
        if (ipart->pdgId() == ipair->first) {
          if (accept(ipart->genParticle())) {
            type1.push_back(ipart);
          }
        } else if (ipart->pdgId() == ipair->second) {
          if (accept(ipart->genParticle())) {
            type2.push_back(ipart);
          }
        }
      }
    }

    // Temporary container of selected particles iterators
    // Useful to compare iterators and avoid double occurrences of the same 
    // particle in case it matches with more than another particle
    vector <ParticleVector::const_iterator> tmp;

    // Now calculate the inv mass
    vector <ParticleVector::const_iterator>::const_iterator i1;
    vector <ParticleVector::const_iterator>::const_iterator i2;
    vector <ParticleVector::const_iterator>::const_iterator begin1;
    vector <ParticleVector::const_iterator>::const_iterator end1;
    begin1 = type1.begin();
    end1 = type1.end();

    for (i1 = begin1; i1 < end1; ++i1) {
      vector <ParticleVector::const_iterator>::const_iterator begin2;
      vector <ParticleVector::const_iterator>::const_iterator end2;
      begin2 = type2.begin();
      end2 = type2.end();
      
      for (i2 = begin2; i2 < end2; ++i2) {
        FourMomentum v4 = (*i1)->momentum() + (*i2)->momentum();
        if (v4.mass() > _minmass && v4.mass() < _maxmass) {
          // Avoid duplicates
          if (find(tmp.begin(), tmp.end(), *i1) == tmp.end()) {
            tmp.push_back(*i1);
            _theParticles.push_back(**i1);
          }
          if (find(tmp.begin(), tmp.end(), *i2) == tmp.end()) {
            tmp.push_back(*i2);
            _theParticles.push_back(**i2);
          }
          getLog() << Log::DEBUG << "Selecting particles with IDs " 
                   << (*i1)->pdgId() << " & " << (*i2)->pdgId()
                   << " and mass = " << v4.mass()/GeV << " GeV" << endl;
        }
      }
    }
    
    getLog() << Log::DEBUG << "Selected " << _theParticles.size() << " particles (" << endl;
    if (getLog().isActive(Log::TRACE)) {  
      foreach (const Particle& p, _theParticles) {
        getLog() << Log::TRACE << "ID: " << p.pdgId() 
                 << ", barcode: " << p.genParticle().barcode() << endl;
      }
    }
  }
 
 
}
