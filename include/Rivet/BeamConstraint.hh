// -*- C++ -*-
#ifndef RIVET_BeamConstraint_HH
#define RIVET_BeamConstraint_HH

#include "Rivet/Rivet.hh"
#include "Rivet/ParticleName.hh"
#include <iostream>


namespace Rivet {

  /// Find whether ParticleName @a p is compatible with the
  /// template ParticleName @a allowed. Effectively this is
  /// asking whether @a p is a subset of @a allowed.
  inline bool compatible(ParticleName p, ParticleName allowed) {
    //assert(p != ANY);
    return (allowed == ANY || p == allowed);
  }

  /// Find whether BeamPair @a pair is compatible with the template
  /// BeamPair @a allowedpair. This assesses whether either of the 
  /// two possible pairings of @a pair's constituents is compatible.
  inline bool compatible(BeamPair pair, BeamPair allowedpair) {
    bool oneToOne = compatible(pair.first, allowedpair.first);
    bool twoToTwo = compatible(pair.second, allowedpair.second);
    bool oneToTwo = compatible(pair.first, allowedpair.second);
    bool twoToOne = compatible(pair.second, allowedpair.first);
    return (oneToOne && twoToTwo) || (oneToTwo && twoToOne);
  }

  /// Find whether a BeamPair @a pair is compatible with at least one template
  /// beam pair in a set @a allowedpairs.
  inline bool subset(BeamPair pair, set<BeamPair> allowedpairs) {
    for (set<BeamPair>::const_iterator bp = allowedpairs.begin(); bp != allowedpairs.end(); ++bp) {
      if (compatible(pair, *bp)) return true;
    }
    return false;
  }


}

#endif
