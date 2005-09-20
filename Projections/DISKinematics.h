// -*- C++ -*-
#ifndef RIVET_DISKinematics_H
#define RIVET_DISKinematics_H
//
// This is the declaration of the DISKinematics class.
//

#include "Rivet/Projections/Projection.h"
#include "Rivet/Projections/Particle.h"
#include "Rivet/Projections/Event.h"
#include "Rivet/Projections/DISLepton.h"

namespace Rivet {

/**
 * This class projects out the DIS kinematic variables and relevant
 * boosts for an event.
 */
class DISKinematics: public Projection {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor. Must specify, the incoming and outgoing
   * (\a inid and \a outid respectively) PDG codes of the scattered
   * lepton as well as the PDG code of the incoming hadron (\a hadid).
   */
  inline DISKinematics(long inid, long outid, long hadid);

  /**
   * The copy constructor.
   */
  inline DISKinematics(const DISKinematics &);

  /**
   * The destructor.
   */
  virtual ~DISKinematics();
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
  virtual int cmp(const Projection & p) const;

public:

  /**
   * The \f$Q^2\f$.
   */
  inline double Q2() const;

  /**
   * The \f$W^2\f$.
   */
  inline double W2() const;

  /**
   * The Bjorken \f$x\f$.
   */
  inline double x() const;

  /**
   * The LorentzRotation needed to boost a particle to the hadronic CM
   * frame.
   */
  inline const LorentzRotation & boostHCM() const;

  /**
   * The LorentzRotation needed to boost a particle to the hadronic Breit
   * frame.
   */
  inline const LorentzRotation & boostBreit() const;

private:

  /**
   * The projector for the scattered lepton.
   */
  DISLepton lepton;

  /**
   * The PDG id of the incoming hadron.
   */
  long idhad;

  /**
   * The \f$Q^2\f$.
   */
  double theQ2();

  /**
   * The \f$W^2\f$.
   */
  double theW2;

  /**
   * The Bjorken \f$x\f$.
   */
  double theX;

  /**
   * The LorentzRotation needed to boost a particle to the hadronic CM
   * frame.
   */
  LorentzRotation hcm;

  /**
   * The LorentzRotation needed to boost a particle to the hadronic Breit
   * frame.
   */
  LorentzRotation breit;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DISKinematics & operator=(const DISKinematics &);

};

}

#include "DISKinematics.icc"

#endif /* RIVET_DISKinematics_H */
