// -*- C++ -*-
#ifndef RIVET_BeamConstraint_HH
#define RIVET_BeamConstraint_HH

#include "Rivet/Tools/Beams.hh"
#include "Rivet/Tools/ParticleName.hh"

namespace Rivet {


  /// @defgroup beam_constraints Standalone beam-compatibility functions
  /// @{


  /// @brief Find whether ParticleName @a p is compatible with the template ParticleName @a allowed
  ///
  /// Effectively this is asking whether @a p is a subset of @a allowed.
  inline bool compatibleBeamID(PdgId p, PdgId allowed) {
    return (allowed == PID::ANY || p == allowed);
  }

  /// @brief Find whether PdgIdPair @a pair is compatible with the template PdgIdPair @a allowedpair
  ///
  /// This assesses whether either of the two possible pairings of @a pair's
  /// constituents is compatible.
  inline bool compatibleBeamIDs(const PdgIdPair& pair, const PdgIdPair& allowedpair) {
    bool oneToOne = compatibleBeamID(pair.first, allowedpair.first);
    bool twoToTwo = compatibleBeamID(pair.second, allowedpair.second);
    bool oneToTwo = compatibleBeamID(pair.first, allowedpair.second);
    bool twoToOne = compatibleBeamID(pair.second, allowedpair.first);
    return (oneToOne && twoToTwo) || (oneToTwo && twoToOne);
  }

  /// Check particle compatibility of two Particle pairs
  inline bool compatibleBeamIDs(const ParticlePair& ppair, const PdgIdPair& allowedpair) {
    return compatibleBeamIDs(PID::make_pdgid_pair(ppair.first.pid(), ppair.second.pid()), allowedpair);
  }

  /// Check particle compatibility of two Particle pairs
  inline bool compatibleBeamIDs(const ParticlePair& a, const ParticlePair& b) {
    return compatibleBeamIDs(PID::make_pdgid_pair(a.first.pid(), a.second.pid()),
                             PID::make_pdgid_pair(b.first.pid(), b.second.pid()) );
  }


  /// Check the energy compatibility of two pairs of particle energies
  inline bool compatibleBeamEnergies(const pair<double,double>& energies,
                                     const pair<double,double>& allowedenergies, double reltol=1e-3) {
    return
      fuzzyEquals(energies.first, allowedenergies.first, reltol) &&
      fuzzyEquals(energies.second, allowedenergies.second, reltol);
  }

  /// Check the energy compatibility of a pair of particles
  inline bool compatibleBeamEnergies(const ParticlePair& ppair,
                                     const pair<double,double>& allowedenergies, double reltol=1e-3) {
    return compatibleBeamEnergies(make_pair(ppair.first.E(), ppair.second.E()), allowedenergies, reltol);
  }

  /// Check the energy compatibility of two Particle pairs
  inline bool compatibleBeamEnergies(const ParticlePair& ppair,
                                     const ParticlePair& allowedppair, double reltol=1e-3) {
    return compatibleBeamEnergies(make_pair(ppair.first.E(), ppair.second.E()),
                                  make_pair(allowedppair.first.E(), allowedppair.second.E()), reltol);
  }


  /// Check the sqrt(s) compatibility of two massless opposing beams
  inline bool compatibleBeamEnergy(const pair<double,double>& energies, double allowedsqrts, double reltol=1e-3) {
    return fuzzyEquals(sqrtS(energies.first, energies.second), allowedsqrts, reltol);
  }

  /// Check the sqrt(s) compatibility of two massless opposing beams
  inline bool compatibleBeamEnergy(const ParticlePair& ppair,
                                   const pair<double,double>& allowedenergies, double reltol=1e-3) {
    return fuzzyEquals(sqrtS(ppair), sqrtS(allowedenergies), reltol);
  }

  /// Check the sqrt(s) compatibility of a Particle pair
  inline bool compatibleBeamEnergy(const ParticlePair& ppair,
                                   const double allowedsqrts, double reltol=1e-3) {
    return fuzzyEquals(sqrtS(ppair), allowedsqrts, reltol);
  }


  /// Check the particle ID and beam-energy compatibility of two Particle pairs
  inline bool compatibleBeams(const ParticlePair& ppair, const PdgIdPair& allowedpids, double allowedsqrts) {
    return compatibleBeamIDs(ppair, allowedpids) && compatibleBeamEnergy(ppair, allowedsqrts);
  }

  /// Check the particle ID and beam-energy compatibility of two Particle pairs
  inline bool compatibleBeams(const ParticlePair& ppair, const PdgIdPair& allowedpids, const pair<double,double>& allowedenergies) {
    return compatibleBeamIDs(ppair, allowedpids) && compatibleBeamEnergies(ppair, allowedenergies);
  }

  /// Check the particle ID and beam-energy compatibility of two Particle pairs
  inline bool compatibleBeams(const ParticlePair& a, const ParticlePair& b) {
    return compatibleBeamIDs(a, b) && compatibleBeamEnergies(a, b);
  }

  /// @}



  // /// Find whether PdgIdPair @a pair is compatible with at least one template beam pair in a set @a allowedpairs
  // inline bool anyCompatibleBeamIDs(const PdgIdPair& pair, const set<PdgIdPair>& allowedpairs) {
  //   for (const PdgIdPair& bp : allowedpairs)
  //     if (compatibleBeamIDs(pair, bp)) return true;
  //   return false;
  // }

  // /// Return the intersection of two sets of {PdgIdPair}s.
  // inline set<PdgIdPair> intersection(const set<PdgIdPair>& a, const set<PdgIdPair>& b) {
  //   set<PdgIdPair> ret;
  //   for (set<PdgIdPair>::const_iterator bp = a.begin(); bp != a.end(); ++bp) {
  //     if (compatibleBeamIDs(*bp, b)) ret.insert(*bp);
  //   }
  //   return ret;
  // }



}

#endif
