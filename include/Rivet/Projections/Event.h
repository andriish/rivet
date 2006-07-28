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

/** Alias for the HepMC namespace. */
namespace CLHEPMC = HepMC;

/** Forward typedef from CLHEPMC. */
typedef CLHEPMC::GenEvent GenEvent;

/** Forward typedef from CLHEPMC. */
typedef CLHEPMC::GenVertex GenVertex;

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
  Event(const HepMC::GenEvent & geneve);

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
  inline const GenEvent & genEvent() const;

  /**
   * Add a projection \a p to this Event. If an equivalent Projection
   * has been applied before, the Projection::project(const Event &)
   * of \a p is not called and a pointer to the previous equivalent
   * projection is returned. If no previous Projection was found, the
   * Projection::project(const Event &) of \a p is called and a
   * pointer to p is returned.
   */
  template <typename PROJ>
  inline const PROJ * addProjection(PROJ & p) const;

  /**
   * Same as addProjection() except that a reference to the projected
   * object is returned.
   */
  template <typename PROJ>
  inline const PROJ & operator()(PROJ & p) const;

  /**
   * The weight associated with the event.
   */
  inline double weight() const;

private:

  /**
   * A pointer to the generated event obtained from an external event
   * generator.
   */
  const GenEvent * theGenEvent;

  /**
   * The set of Projection objects applied so far.
   */
  mutable std::set<const Projection *> theProjections;

  /**
   * The weight associated with the event.
   */
  double theWeight;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Event & operator=(const Event &);

};

}

#include "Rivet/Projections/Event.icc"

#endif /* RIVET_Event_H */
