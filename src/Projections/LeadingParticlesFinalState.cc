#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Particle.hh"

namespace Rivet {

  int LeadingParticlesFinalState::compare(const Projection & p) const {
    //first compare the final states we are running on
    int fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != PCmp::EQUIVALENT)
       return fscmp;
    //then compare the two as final states
    const LeadingParticlesFinalState & other = dynamic_cast < const LeadingParticlesFinalState & >(p);
     fscmp = FinalState::compare(other);
    if (fscmp != PCmp::EQUIVALENT)
       return fscmp;
    //finally compare the ids
    if (_ids < other._ids)
      return PCmp::ORDERED;
    else if (other._ids < _ids)
      return PCmp::UNORDERED;
    return PCmp::EQUIVALENT;
  } 
  
  void LeadingParticlesFinalState::project(const Event & e) {

    Log & log = getLog();

    _theParticles.clear();
    const FinalState & fs = applyProjection < FinalState > (e, "FS");

    /// temporary container
     map < long, ParticleVector::const_iterator > tmp;

    const ParticleVector & particles = fs.particles();
     log << Log::DEBUG << "Original final state particles size " << particles.size() << endl;

     ParticleVector::const_iterator ifs;
    for (ifs = particles.begin(); ifs != particles.end(); ++ifs) {
      if (inList(*ifs) && FinalState::accept(ifs->getHepMCParticle())) {
        //look for an existing particle in tmp container
        map < long, ParticleVector::const_iterator >::const_iterator itmp = tmp.find(ifs->getPdgId());
        if (itmp != tmp.end()) {  // if a particle with this type has been already selected 
          //if the new pT is higher than the previous one substitute
          if (ifs->getMomentum().pT() > itmp->second->getMomentum().pT())
            tmp[ifs->getPdgId()] = ifs;
          // else insert in the container
        } else
           tmp[ifs->getPdgId()] = ifs;
      }
    }

    // loop on the tmp container and fill _theParticles
    map < long, ParticleVector::const_iterator >::const_iterator i;
    for (i = tmp.begin(); i != tmp.end(); ++i) {
      log << Log::DEBUG << "LeadingParticlesFinalState is accepting particle ID " << i->second->getPdgId()
        << " with momentum " << i->second->getMomentum() << endl;
      _theParticles.push_back(*(i->second));
    }
  }

  bool LeadingParticlesFinalState::inList(const Particle & particle) const {
    std::set < long >::const_iterator ilist = _ids.find(particle.getPdgId());
    if (ilist != _ids.end())
       return true;
    return false;
  }
}
