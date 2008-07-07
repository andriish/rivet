// -*- C++ -*-
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {

  int InvMassFinalState::compare(const Projection & p) const {
    //first compare the final states we are running on
    int fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != PCmp::EQUIVALENT)
       return fscmp;
    const InvMassFinalState & other = dynamic_cast < const InvMassFinalState & >(p);
    //then compare the two as final states
     fscmp = FinalState::compare(other);
    if (fscmp != PCmp::EQUIVALENT)
       return fscmp;
    // then compare the mass limits
    int massllimcmp = cmp(_minmass, other._minmass);
    if (massllimcmp != PCmp::EQUIVALENT)
       return massllimcmp;
    int masshlimcmp = cmp(_maxmass, other._maxmass);
    if (masshlimcmp != PCmp::EQUIVALENT)
       return masshlimcmp;
    //compare the decay species
    int decaycmp = cmp(_decayids, other._decayids);
    if (decaycmp != PCmp::EQUIVALENT)
       return decaycmp;
    //finally compare them as final states 
     return FinalState::compare(other);
  } 
  
  void InvMassFinalState::project(const Event & e) {
    Log & log = getLog();
    const FinalState & fs = applyProjection < FinalState > (e, "FS");
     _theParticles.clear();

    //containers for the particles of type specified in the pair (if the id is different)
     vector < ParticleVector::const_iterator > type1;
     vector < ParticleVector::const_iterator > type2;
    //containers for the particles of type specified in the pair (if the id is same) 
     vector < ParticleVector::const_iterator > type3;
    // get all the particles of the type specified in the pair from the particle list
    const ParticleVector & allparticles = fs.particles();
     ParticleVector::const_iterator ipart;
    for (ipart = allparticles.begin(); ipart != allparticles.end(); ++ipart) {
      if (_decayids.first != _decayids.second) {
        if (ipart->getPdgId() == _decayids.first) {
          if (accept(ipart->getHepMCParticle()))
            type1.push_back(ipart);
        } else if (ipart->getPdgId() == _decayids.second) {
          if (accept(ipart->getHepMCParticle()))
            type2.push_back(ipart);
        }
      } else {
        if (ipart->getPdgId() == _decayids.first) {
          if (accept(ipart->getHepMCParticle()))
            type3.push_back(ipart);
        }
      }
    }

    /// temporaty container of selected particles iterators
    /// useful to compare iterators and avoid double occurrences of the same particle
    /// in case it matches with more than another particle
    vector < ParticleVector::const_iterator > tmp;

    //now calculate the inv mass
    bool equalDecayProducts = (_decayids.first == _decayids.second);
    vector < ParticleVector::const_iterator >::const_iterator i1;
    vector < ParticleVector::const_iterator >::const_iterator i2;
    vector < ParticleVector::const_iterator >::const_iterator begin1;
    vector < ParticleVector::const_iterator >::const_iterator end1;
    if (equalDecayProducts) {
      begin1 = type3.begin();
      end1 = type3.end() - 1;
    } else {
      begin1 = type1.begin();
      end1 = type1.end();
    }
    for (i1 = begin1; i1 < end1; ++i1) {
      vector < ParticleVector::const_iterator >::const_iterator begin2;
      vector < ParticleVector::const_iterator >::const_iterator end2;
      if (equalDecayProducts) {
        begin2 = i1 + 1;
        end2 = type3.end();
      } else {
        begin2 = type2.begin();
        end2 = type2.end();
      }
      for (i2 = begin2; i2 < end2; ++i2) {
        FourMomentum vector = (*i1)->getMomentum() + (*i2)->getMomentum();
        if (vector.mass() > _minmass && vector.mass() < _maxmass) {
          //avoid duplicates
          if (find(tmp.begin(), tmp.end(), *i1) == tmp.end()) {
            tmp.push_back(*i1);
            _theParticles.push_back(**i1);
          }
          if (find(tmp.begin(), tmp.end(), *i2) == tmp.end()) {
            tmp.push_back(*i2);
            _theParticles.push_back(**i2);
          }
          log << Log::DEBUG << "InvMassFinalState is selecting particles with id " << (*i1)->getPdgId() << " and " << (*i2)->getPdgId()
            << " with mass " << vector.mass() << endl;
        }
      }
    }
    if (log.isActive(Log::DEBUG)) {
      stringstream msg;
      for (ParticleVector::const_iterator ipart = _theParticles.begin(); ipart != _theParticles.end(); ++ipart) {
        msg << "ID " << ipart->getPdgId() << " barcode " << ipart->getHepMCParticle().barcode() << ", ";
      }
      log << Log::DEBUG << "The following " << _theParticles.size() << " particles have been selected by InvMassFinalState: " << msg.str() << endl;
    }
  }

}
