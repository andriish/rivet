// -*- C++ -*-

#include "Rivet/Tools/RivetHepMC.hh"

namespace Rivet{
  
  std::vector<ConstGenParticlePtr> particles(const GenEvent *ge){
    assert(ge);
    return ge->particles();
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
  
}
