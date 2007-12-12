// -*- C++ -*-
#ifndef RIVET_Event_HH
#define RIVET_Event_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Projection.hh"
#include "HepMC/GenEvent.h"
#include "Event.fhh"


namespace Rivet {

  /// Forward typedef from HepMC.
  typedef HepMC::GenEvent GenEvent;
  
  /// Forward typedef from HepMC.
  typedef HepMC::GenVertex GenVertex;


  /// Event is a concrete class representing an generated event in
  /// Rivet. It is constructed given a HepMC::GenEvent, a pointer to
  /// which is kept by the Event object throughout its lifetime. The user
  /// must therefore make sure that the corresponding HepMC::GenEvent
  /// will persist at least as long as the Event object.
  ///
  /// In addition to the HepMC::GenEvent object the Event also keeps
  /// track of all Projections object which have been applied to the
  /// Event so far.
  class Event {

  public:

    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor.
    Event(const GenEvent& geneve)
      : theGenEvent(&geneve), theWeight(1.0) {
      if ( !geneve.weights().empty() ) theWeight = geneve.weights()[0];
    }

    /// The copy constructor.
    Event(const Event& x)  
      : theGenEvent(x.theGenEvent), theWeight(x.theWeight) 
    { }

    /// The destructor.
    ~Event() { }
    //@}

  public:

    /// Return the generated event obtained from an external event generator.
    const GenEvent& genEvent() const { return *theGenEvent; }

    /// Add a projection \a p to this Event. If an equivalent Projection
    /// has been applied before, the Projection::project(const Event &)
    /// of \a p is not called and a reference to the previous equivalent
    /// projection is returned. If no previous Projection was found, the
    /// Projection::project(const Event &) of \a p is called and a
    /// reference to p is returned.
    template <typename PROJ>
    const PROJ& applyProjection(PROJ& p) const {
      ConstProjectionPtr cpp(&p);
      std::set<ConstProjectionPtr>::const_iterator old = theProjections.find(cpp);
      if (old != theProjections.end()) {
        //return *( dynamic_cast<const PROJ*>(old->get()) );
        return *( dynamic_cast<const PROJ*>(*old) );
      }
      // Add the projection via the Projection base class (only 
      // possible because Event is a friend of Projection)
      ProjectionPtr pp(&p);
      pp->project(*this);
      theProjections.insert(pp);
      return p;
    }

    /// The weight associated with the event.
    double weight() const { return theWeight; }

  private:

    /// A pointer to the generated event obtained from an external generator.
    const GenEvent* theGenEvent;

    /// The set of Projection objects applied so far.
    mutable std::set<ConstProjectionPtr> theProjections;

    /// The weight associated with the event.
    double theWeight;

  private:

    /// The assignment operator is private and must never be called.
    /// In fact, it should not even be implemented.
    Event& operator=(const Event&);

  };
  
}

#endif
