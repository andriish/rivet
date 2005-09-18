// -*- C++ -*-
#ifndef RIVET_Event_H
#define RIVET_Event_H
//
// This is the declaration of the Event class.
//

#include "Rivet/Config/Rivet.h"
#include "Event.fh"
#include "Rivet/Projections/Projection.h"
#include "CLHEP/HepMC/GenEvent.h"
namespace Rivet {

/**
 * Event is a concrete class representing an generated event in
 * Rivet. It is constructed given a HepMC::GenEvent, a pointer to
 * which is kept by the Event object throughout its lifetime. The user
 * must therefore make sure that the corresponding HepMC::GenEvent
 * will persist at least as long as the Event object.
 *
 * In addition to the HepMC::GenEvent object the Event also keeps
 * track of all Projections object which have been applied to the
 * Event so far.
 */
class Event {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Event(const HepMC::GenEvent & geneve);

  /**
   * The copy constructor.
   */
  inline Event(const Event &);

  /**
   * The destructor.
   */
  virtual ~Event();
  //@}

public:


  /**
   * Return the generated event obtained from an external event
   * generator.
   */
  inline const HepMC::GenEvent & genEvent() const;

  /**
   * Add a projection \a p to this Event. If an equivalent Projection
   * has been applied before, the Projection::project(const Event &)
   * of \a p is not called and a pointer to the previous equivalent
   * projection is returned. If no previous Projection was found, the
   * Projection::project(const Event &) of \a p is called and a
   * pointer to p is returned.
   */
  template <typename PROJ>
  inline const PROJ * addProjection(Proj & p) const;

private:

  /**
   * A pointer to the generated event obtained from an external event
   * generator.
   */
  const HepMC::GenEvent * theGenEvent;

  /**
   * The set of Projection objects applied so far.
   */
  set<const Projection *> theProjections;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Event & operator=(const Event &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Event. */
template <>
struct BaseClassTrait<Rivet::Event,1> {
  /** Typedef of the first base class of Event. */
  typedef  NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Event class and the shared object where it is defined. */
template <>
struct ClassTraits<Rivet::Event>
  : public ClassTraitsBase<Rivet::Event> {
  /** Return a platform-independent class name */
  static string className() { return "Rivet::Event"; }
  /** Return the name of the shared library be loaded to get
   *  access to the Event class and every other class it uses
   *  (except the base class). */
  static string library() { return "Event.so"; }
};

/** @endcond */

}

#include "Event.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Event.tcc"
#endif

#endif /* RIVET_Event_H */
