// -*- C++ -*-
#ifndef RIVET_PartonicTops_HH
#define RIVET_PartonicTops_HH

#include "Rivet/Projections/ParticleFinder.hh"

namespace Rivet {


  /// @brief Convenience finder of partonic top quarks
  ///
  /// @warning Requires there to be tops in the event record. A fiducial pseudo-top
  /// analysis approach is strongly recommended instead of this.
  class PartonicTops : public ParticleFinder {
  public:


    /// @brief Enum for categorising top quark decay modes
    ///
    /// More specifically, the decay mode of the W from the top. We presume top decay to a W and b quark.
    ///
    /// @note E_MU mode does not include intermediate taus, while LEPTONIC does.
    ///   Similarly the QUARKS mode does not include hadronic taus, while HADRONIC does.
    enum DecayMode { ELECTRON, MUON, E_MU, TAU, TAU_LEPTONIC, TAU_HADRONIC, LEPTONIC, QUARKS, HADRONIC, ALL };


    /// @name Constructors
    //@{

    /// Constructor optionally taking cuts object
    PartonicTops(const Cut& c=Cuts::OPEN)
      : ParticleFinder(c), _decaymode(ALL)
    {  }

    /// Constructor taking decay mode (and an optional cuts object)
    PartonicTops(DecayMode decaymode, const Cut& c=Cuts::OPEN)
      : ParticleFinder(c), _decaymode(decaymode)
    {  }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(PartonicTops);

    //@}


    /// Access to the found partonic tops
    const Particles& tops() const { return _theParticles; }


    /// Clear the projection
    void clear() {
      _theParticles.clear();
    }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& event) {
      // Find partonic tops
      _theParticles = filter_select(event.allParticles(_cuts), isLastWith(isTop));
      // Filtering by decay mode
ELECTRON, MUON, E_MU, TAU, TAU_LEPTONIC, TAU_HADRONIC, LEPTONIC, QUARKS, HADRONIC, ALL
    }

    /// Compare projections.
    int compare(const Projection& p) const {
      const PartonicTops& other = dynamic_cast<const PartonicTops&>(p);
      return cmp(_cut, other._cut);
    }


  private:

    DecayMode _decaymode;


  };


}


#endif
