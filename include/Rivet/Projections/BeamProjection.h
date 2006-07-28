// -*- C++ -*-
#ifndef RIVET_BeamProjection_H
#define RIVET_BeamProjection_H
//
// This is the declaration of the BeamProjection class.
//

#include "Rivet/Projections/Projection.h"
#include "Rivet/Projections/Event.h"
#include "Rivet/Projections/Particle.h"

namespace Rivet {

/**
 * This class is used to project out the beams in a HepMC::GenEvent.
 */
class BeamProjection: public Projection {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline BeamProjection();

  /**
   * The copy constructor.
   */
  inline BeamProjection(const BeamProjection &);

  /**
   * The destructor.
   */
  virtual ~BeamProjection();
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
   * The pair of beam particles in the current collision in GenEvent 
   */
  inline const PPair & operator()() const;

private:

  /**
   * The beam particles in the current collision in GenEvent 
   */
  PPair theBeams;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BeamProjection & operator=(const BeamProjection &);

};

}

#include "Rivet/Projections/BeamProjection.icc"

#endif /* RIVET_BeamProjection_H */
