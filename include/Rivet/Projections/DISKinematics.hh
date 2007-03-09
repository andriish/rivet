// -*- C++ -*-
#ifndef RIVET_DISKinematics_H
#define RIVET_DISKinematics_H
//
// This is the declaration of the DISKinematics class.
//

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Tools/Event/Particle.hh"
#include "Rivet/Tools/Event/Event.hh"
#include "Rivet/Projections/DISLepton.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/RivetCLHEP.hh"

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
   * The default constructor. Must specify Beam and DISLepton
   * projection objects which are guaranteed to live throughout the
   * run. Also the PDG code of the incoming hadron (\a hadid) must be
   * specified.
   */
  inline DISKinematics(Beam & beamp, DISLepton & leptonp, long hadid)
    : beams(&beamp), lepton(&leptonp), idhad(hadid), theQ2(-1.0), theW2(-1.0),
      theX(-1.0) {
    info.declareParameter("BeamB", hadid);
  }

  /**
   * The copy constructor.
   */
  inline DISKinematics(const DISKinematics & x)
    : Projection(x), beams(x.beams), lepton(x.lepton), idhad(x.idhad),
      theQ2(x.theQ2), theW2(x.theW2), theX(x.theX), hcm(x.hcm), breit(x.breit) 
  { }
  
  /**
   * The destructor.
   */
  virtual ~DISKinematics() { };
  //@}

protected:

  /**
   * Take the information available in the Event and make the
   * calculations necessary to obtain the projection. Note that this
   * function must never be called except inside the
   * Event::applyProjection(Projection *) function. If the information
   * from other projections are necessary, their project(const Event
   * &) should not be called, rather the corresponding objects should
   * be added to the Event using the Event::applyProjection(Projection *)
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
   * The \f$Q^2\f$.
   */
  inline double Q2() const { return theQ2; }

  /**
   * The \f$W^2\f$.
   */
  inline double W2() const { return theW2; }

  /**
   * The Bjorken \f$x\f$.
   */
  inline double x() const { return theX; }

  /**
   * The LorentzRotation needed to boost a particle to the hadronic CM
   * frame.
   */
  inline const LorentzRotation & boostHCM() const {
    return hcm; 
  }

  /**
   * The LorentzRotation needed to boost a particle to the hadronic Breit
   * frame.
   */
  inline const LorentzRotation & boostBreit() const {
    return breit;
  }

  /**
   * Return the RivetInfo object of this Projection. Derived classes
   * should re-implement this function to return the combined
   * RivetInfo object of this and of any other Projection upon which
   * this depends.
   */
  virtual RivetInfo getInfo() const;

private:

  /**
   * The Beam projector object defining the incoming beam particles.
   */
  Beam * beams;

  /**
   * The projector for the scattered lepton.
   */
  DISLepton * lepton;

  /**
   * The PDG id of the incoming hadron.
   */
  long idhad;

  /**
   * The \f$Q^2\f$.
   */
  double theQ2;

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

#endif /* RIVET_DISKinematics_H */
