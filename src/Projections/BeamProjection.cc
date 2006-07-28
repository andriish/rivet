// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BeamProjection class.
//

#include "BeamProjection.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BeamProjection.tcc"
#endif


using namespace Rivet;

BeamProjection::~BeamProjection() {}

void BeamProjection::project(const Event & e) {
  vector<GenParticle*> inc =
    e.genEvent().signal_process_vertex()->listParents();
  if ( inc.size() != 2 )
    throw runtime_error("Wrong number of beams.");
  // *** ATTENTION *** Maybe we should have our own exception classes.

  theBeams.first = Particle(*inc[0]);
  theBeams.second = Particle(*inc[1]);

}

int BeamProjection::compare(const Projection &) const {
  return 0;
}


