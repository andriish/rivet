// -*- C++ -*-
#ifndef RIVET_PartonicTopFinder_HH
#define RIVET_PartonicTopFinder_HH

#include "Rivet/Projections/ParticleFinder.hh"

namespace Rivet {


  /// @brief Convenience finder of partonic top quarks
  ///
  /// @warning Requires there to be tops in the event record. A fiducial pseudo-top
  /// analysis approach is strongly recommended instead of this.
  class PartonicTopFinder : public ParticleFinder {
  public:

    /// @name Constructors
    //@{

    /// Constructor taking cuts object
    /// @todo Parametrise on accepted decay types: ELECTRONIC, MUONIC, LEPTONIC, HADRONIC, ALL
    PartonicTopFinder(const Cut& c=Cuts::OPEN)
      : _cut(c)
    {  }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(PartonicTopFinder);

    //@}


    /// Access to the found partonic tops
    const Particles& tops() const { return _tops; }

    // /// Find the particles not descended from the found tops
    // const VetoedFinalState& remainingFinalState() const;


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& event) {
      // Find partonic tops
      _tops.clear();
      for (const Particle& p : event.allParticles(_cut)) {
        if (!isTop(p)) continue;
        if (any(p.children(), isTop)) continue;
        _tops += p;
      }
    }

    /// Compare projections.
    int compare(const Projection& p) const {
      const PartonicTopFinder& other = dynamic_cast<const PartonicTopFinder&>(p);
      return cmp(_cut, other._cut);
    }


  public:

    /// Clear the projection
    void clear() {
      _tops.clear();
    }


  private:

    /// The cut
    Cut _cut;

    /// The discovered partonic tops
    Particles _tops;

  };


}


#endif
