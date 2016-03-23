// -*- C++ -*-
#ifndef RIVET_SmearingFunctions_HH
#define RIVET_SmearingFunctions_HH

#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"

namespace Rivet {


  /// Return a uniformly sampled random number between 0 and 1
  /// @todo Move to (math?)utils
  double rand01() { return rand() / (double)RAND_MAX; }


  /// @name Standard particle efficiency and smearing functions
  //@{

  double PARTICLE_FN0(const Particle& p) { return 0; }
  double PARTICLE_FN1(const Particle& p) { return 1; }

  double P4_FN0(const FourMomentum& p) { return 0; }
  double P4_FN1(const FourMomentum& p) { return 1; }

  Particle PARTICLE_SMEAR_IDENTITY(const Particle& p) { return p; }

  //@}


  /// @name Standard jet efficiency and smearing functions
  //@{

  /// Return a constant 0 given a Jet as argument
  double JET_EFF_ZERO(const Jet& p) { return 0; }
  /// Return a constant 1 given a Jet as argument
  double JET_EFF_ONE(const Jet& p) { return 1; }

  /// Return 1 if the given Jet contains a b, otherwise 0
  double JET_BTAG_PERFECT(const Jet& j) { return j.bTagged() ? 1 : 0; }
  /// Return the ATLAS Run 1c jet flavour tagging efficiency for the given Jet
  double JET_BTAG_ATLAS_RUN1C(const Jet& j) {
    if (j.bTagged()) return 0.80*tanh(0.003*j.pT()/GeV)*(30/(1+0.086*j.pT()/GeV));
    if (j.cTagged()) return 0.20*tanh(0.02*j.pT()/GeV)*(1/(1+0.0034*j.pT()/GeV));
    return 0.002 + 7.3e-6*j.pT()/GeV;
  }

  /// Return 1 if the given Jet contains a c, otherwise 0
  double JET_CTAG_PERFECT(const Jet& j) { return j.cTagged() ? 1 : 0; }

  /// Take a jet and return an unmodified copy
  /// @todo Modify constituent particle vectors for consistency
  /// @todo Set a null PseudoJet if the Jet is smeared?
  Jet JET_SMEAR_IDENTITY(const Jet& j) { return j; }

  //@}


}

#endif
