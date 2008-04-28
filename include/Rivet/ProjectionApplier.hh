// -*- C++ -*-
#ifndef RIVET_ProjectionApplier_HH
#define RIVET_ProjectionApplier_HH

#include "Rivet/Rivet.hh"

namespace Rivet {

  /// Empty interface used for storing Projection and Analysis pointers in the
  /// same container (used by the ProjectionHandler)
  class ProjectionApplier {
  public:
    virtual ~ProjectionApplier() { }
  };

}

#endif
