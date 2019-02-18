// -*- C++ -*-

#include "Rivet/Tools/RivetHepMC.hh"

namespace Rivet{
  
  std::vector<ConstGenParticlePtr> particles(ConstGenEventPtr ge){
    return ge->particles();
  }

  std::vector<ConstGenParticlePtr> particles(const GenEvent *ge){
    assert(ge);
    return ge->particles();
  }
  
  std::vector<ConstGenVertexPtr> vertices(ConstGenEventPtr ge){
    return ge->vertices();
  }

  std::vector<ConstGenVertexPtr> vertices(const GenEvent *ge){
    assert(ge);
    return ge->vertices();
  }
  
  std::vector<ConstGenParticlePtr> particles(ConstGenVertexPtr gv, const Relatives &relo){
    return relo(gv);
  }

  std::vector<ConstGenParticlePtr> particles(ConstGenParticlePtr gp, const Relatives &relo){
    return relo(gp);
  }
  
  int uniqueId(ConstGenParticlePtr gp){
    return gp->id();
  }
  
  std::vector<ConstGenParticlePtr> beams(const GenEvent *ge){
    return ge->beams();
  }
  
}
