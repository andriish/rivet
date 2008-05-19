// -*- C++ -*-
#ifndef RIVET_Event_HH
#define RIVET_Event_HH

#include "Rivet/Rivet.hh"
#include "Event.fhh"
#include "Rivet/Projection.hh"


namespace Rivet {

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
      : _genEvent(&geneve), _weight(1.0) {
      if ( !geneve.weights().empty() ) _weight = geneve.weights()[0];
    }

    /// The copy constructor.
    Event(const Event& e)  
      : _genEvent(e._genEvent), _weight(e._weight) 
    { }

    /// The destructor.
    ~Event() { }
    //@}

  public:

    /// Return the generated event obtained from an external event generator.
    const GenEvent& genEvent() const { return *_genEvent; }

    /// Add a projection \a p to this Event. If an equivalent Projection
    /// has been applied before, the Projection::project(const Event &)
    /// of \a p is not called and a reference to the previous equivalent
    /// projection is returned. If no previous Projection was found, the
    /// Projection::project(const Event &) of \a p is called and a
    /// reference to p is returned.
    template <typename PROJ>
    const PROJ& applyProjection(PROJ& p) const {
      const Projection* cpp(&p);
      std::set<const Projection*>::const_iterator old = _projections.find(cpp);
      if (old != _projections.end()) {
        const Projection& pRef = **old;
        return pcast<PROJ>(pRef);
      }
      // Add the projection via the Projection base class (only 
      // possible because Event is a friend of Projection)
      Projection* pp = const_cast<Projection*>(cpp);
      pp->project(*this);
      _projections.insert(pp);
      return p;
    }

    template <typename PROJ>
    const PROJ& applyProjection(PROJ* pp) const {
      if (!pp) throw Error("Event::applyProjection(PROJ*): Projection pointer is null.");
      return applyProjection(*pp);
    }

    /// The weight associated with the event.
    double weight() const { return _weight; }

  private:

    /// A pointer to the generated event obtained from an external generator.
    const GenEvent* _genEvent;

    /// The set of Projection objects applied so far.
    mutable std::set<ConstProjectionPtr> _projections;

    /// The weight associated with the event.
    double _weight;

  private:

    /// The assignment operator is private and must never be called.
    /// In fact, it should not even be implemented.
    Event& operator=(const Event&);

  };
  
}

#endif
