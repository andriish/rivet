// -*- C++ -*-
#ifndef RIVET_RivetHepMC_HH
#define RIVET_RivetHepMC_HH

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"

#include <vector>

namespace Rivet {


  using HepMC::GenEvent;
  using HepMC::GenParticle;
  using HepMC::GenVertex;


  inline vector<GenParticle*> particles(const GenEvent& ge) {
    vector<GenParticle*> rtn;
    for (GenEvent::particle_const_iterator pi = ge.particles_begin(); pi != ge.particles_end(); ++pi) {
      rtn.push_back(*pi);
    }
    return rtn;
  }
  inline vector<GenParticle*> particles(const GenEvent* ge) {
    assert(ge);
    return particles(*ge);
  }


  inline vector<GenVertex*> vertices(const GenEvent& ge) {
    vector<GenVertex*> rtn;
    for (GenEvent::vertex_const_iterator vi = ge.vertices_begin(); vi != ge.vertices_end(); ++vi) {
      rtn.push_back(*vi);
    }
    return rtn;
  }
  inline vector<GenVertex*> vertices(const GenEvent* ge) {
    assert(ge);
    return vertices(*ge);
  }


  inline const vector<GenParticle*> particles_in(const GenVertex* gv) {
    vector<GenParticle*> rtn;
    for (GenVertex::particles_in_const_iterator pi = gv->particles_in_const_begin(); pi != gv->particles_in_const_end(); ++pi) {
      rtn.push_back(*pi);
    }
    return rtn;
  }


  inline const vector<GenParticle*> particles_out(const GenVertex* gv) {
    vector<GenParticle*> rtn;
    for (GenVertex::particles_out_const_iterator pi = gv->particles_out_const_begin(); pi != gv->particles_out_const_end(); ++pi) {
      rtn.push_back(*pi);
    }
    return rtn;
  }


  inline vector<GenParticle*> particles(GenVertex* gv, HepMC::IteratorRange range=HepMC::relatives) {
    vector<GenParticle*> rtn;
    for (GenVertex::particle_iterator pi = gv->particles_begin(range); pi != gv->particles_end(range); ++pi) {
      rtn.push_back(*pi);
    }
    return rtn;
  }


}

#endif
