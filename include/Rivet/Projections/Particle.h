// -*- C++ -*-
#ifndef RIVET_Particle_H
#define RIVET_Particle_H
//
// This is the declaration of the Particle class.
//

#include "Rivet/Rivet.h"
#include "Rivet/CLHEPWrap/LorentzVector.h"
#include "CLHEP/HepMC/GenParticle.h"
#include "Particle.fh"

namespace Rivet {

/** Alias for the HepMC namespace. */
namespace CLHEPMC = HepMC;

/** Forward typedefs from CLHEPMC. */
typedef CLHEPMC::GenParticle GenParticle;

/**
 * This class is used to represent particles projected from a
 * HepMC::GenEvent. It is currently just a struct to access the
 * members directly, but it should probably be promoted to a proper
 * class.
 */
struct Particle {

  /**
   * Simple constructor for copying a GenParticle.
   */
  inline Particle(const GenParticle & p);

  /**
   * Default constructor.
   */
  inline Particle();

  /**
   * Copy-constructor.
   */
  inline Particle(const Particle & p);

  /**
   * A pointer to the original GenParticle from which this Particle is
   * projected.
   */
  const GenParticle * original;

  /**
   * The PDG id number for this Particle.
   */
  long id;

  /**
   * The momentum of this projection of the Particle.
   */
  LorentzVector momentum;

  /**
   * The mass of this Particle.
   */
  double mass;

};

/** Typedef a vector of Particle objects. */
typedef std::vector<Particle> PVector;

/** Typedef a pair of Particle objects. */
typedef std::pair<Particle, Particle> PPair;

}

#include "Rivet/Projections/Particle.icc"

#endif /* RIVET_Particle_H */
