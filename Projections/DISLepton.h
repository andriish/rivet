// -*- C++ -*-
#ifndef RIVET_DISLepton_H
#define RIVET_DISLepton_H
//
// This is the declaration of the DISLepton class.
//

#include "Rivet/Projections/BeamProjection.h"
#include "Rivet/Projections/Particle.h"
#include "Rivet/Projections/Event.h"

namespace Rivet {

/**
 * This class projects out the incoming and outgoing leptons in a DIS
 * event. The incoming incoming lepton is assumed to be along the
 * positive z-axis.
 */
class DISLepton: public Projection {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor. Must specify the incoming and outgoing
   * PDG codes of the leptons to project.  If \a inid is the
   * anti-particle of \a outid, either a scattered lepton or
   * anti-lepton is searched for.
   */
  inline DISLepton(long inid, long outid);

  /**
   * The copy constructor.
   */
  inline DISLepton(const DISLepton &);

  /**
   * The destructor.
   */
  virtual ~DISLepton();
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
  virtual void project(const Event & e);

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
  virtual int compare(const Projection & p) const;

public:

  /**
   * The incoming lepton.
   */
  inline const Particle & in() const;

  /**
   * The outgoing lepton.
   */
  inline const Particle & out() const;

  /**
   * Return the RivetInfo object of this Projection. Derived classes
   * should re-implement this function to return the combined
   * RivetInfo object of this and of any other Projection upon which
   * this depends.
   */
  virtual RivetInfo getInfo() const;

private:

  /**
   * The BeamProjector object defining the incoming beam particles.
   */
  BeamProjection beams;

  /**
   * The PDG id of the incoming lepton.
   */
  long idin;

  /**
   * The PDG id of the outcoming lepton.
   */
  long idout;

  /**
   * The incoming lepton.
   */
  Particle incoming;

  /**
   * The incoming lepton.
   */
  Particle outgoing;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DISLepton & operator=(const DISLepton &);

};

}

#include "DISLepton.icc"

#endif /* RIVET_DISLepton_H */
