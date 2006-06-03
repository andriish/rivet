// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Event class.
//

#include "Event.h"

using namespace Rivet;

Event::Event(const GenEvent & geneve)
  : theGenEvent(&geneve), theWeight(1.0) {
  if ( geneve.weights().size() ) theWeight = geneve.weights()[0];
}

Event::~Event() {}


