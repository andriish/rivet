// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Projection class.
//

#include "Rivet/Projections/Projection.h"

using namespace Rivet;

Projection::~Projection() {}

RivetInfo Projection::getInfo() const {
  return info;
}



