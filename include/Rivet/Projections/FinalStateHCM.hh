// -*- C++ -*-
#ifndef RIVET_FinalStateHCM_H
#define RIVET_FinalStateHCM_H
//
// This is the declaration of the FinalStateHCM class.

#include "Rivet/Tools/Event/Particle.hh"
#include "Rivet/Tools/Event/Event.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"


namespace Rivet {

/**
 * Project all final state particles (except the scattered lepton)
 * boosted to the hadronic center of mass system.
 */
class FinalStateHCM: public Projection {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor. Must specify the PDG id of the incoming
   * and scattered lepton and of the incoming hadron. May also specify
   * the minimum and maximum pseudorapidity (in the lab-system).
   */
  inline FinalStateHCM(long inid, long outid, long hadid,
		       double mineta = -MaxRapidity,
		       double maxeta = MaxRapidity);

  /**
   * The copy constructor.
   */
  inline FinalStateHCM(const FinalStateHCM &);

  /**
   * The destructor.
   */
  virtual ~FinalStateHCM();
  //@}

protected:

  /**
   * Take the information available in the Event and make the
   * calculations necessary to obtain the projection. Note that this
   * function must never be called except inside the
   * Event::addProjection(Projection *) function. If the information
   * from other projections are necessary, their project(const Event
   * &) should not be called, rather the corresponding objects should
   * be added to the Event using the Even::addProjection(Projection *)
   * function.
   */
  void project(const Event & e);

  /**
   * This function is used to define a unique ordering between
   * different Projection objects of the same class. If this is
   * considered to be equivalent to the Projector object, \a p, in the
   * argument the function should return 0. If this object should be
   * ordered before \a p a negative value should be returned,
   * otherwise a positive value should be returned. This function must
   * never be called explicitly, but should only be called from the
   * operator<(const Projection &). When implementing the function in
   * concrete sub-classes, it is then guarranteed that the Projection
   * object \a p in the argument is of the same class as the sub-class
   * and can be safely dynamically casted to that class.
   *
   * When implementing this function in a sub-class, the immediate
   * base class version of the function should be called first. If the
   * base class function returns a non-zero value, that value should
   * be returned immediately. Only if zero is returned should this
   * function check the member variables of the sub-class to determine
   * whether this should be ordered before or after \a p, or if it is
   * equivalent with \a p.
   */
  int compare(const Projection & p) const;

public:

  /**
   * Access the projected final-state particles.
   */
  inline const ParticleVector & particles() const;

  /**
   * Return the RivetInfo object of this Projection. Derived classes
   * should re-implement this function to return the combined
   * RivetInfo object of this and of any other Projection upon which
   * this depends.
   */
  virtual RivetInfo getInfo() const;

private:

  /**
   * The projector for the DIS kinematics.
   */
  DISLepton lepton;

  /**
   * The projector for the DIS kinematics.
   */
  DISKinematics kinematics;

  /**
   * The projector for the full final state.
   */
  FinalState fsproj;

  /**
   * The final-state particles.
   */
  ParticleVector theParticles;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FinalStateHCM & operator=(const FinalStateHCM &);

};

}

#include "Rivet/Projections/FinalStateHCM.icc"

#endif /* RIVET_FinalStateHCM_H */
