// -*- C++ -*-
#ifndef RIVET_PVertex_H
#define RIVET_PVertex_H

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Event.hh"
#include "Rivet/Particle.hh"

namespace Rivet {
  
  
  /// @todo This is wrong: This class is used to project out the beams in a HepMC::GenEvent.
  class PVertex: public Projection {

  public:

    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor.
    inline PVertex() { 
    }

    /// The copy constructor.
    inline PVertex(const PVertex& x)
      : Projection(x), thePVertex(x.thePVertex) { }

    /// The destructor.
    virtual ~PVertex() { }
    //@}

    public:
      /// Return the name of the projection
      inline string getName() const {
        return "PVertex";
      }

  protected:

    /*
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
     * The pair of beam particles in the current collision in GenEvent 
     */
    inline const GenVertex & operator()() const {
     return *thePVertex;
    }

  private:

    /**
     * The Primary Vertex in the current collision in GenEvent 
     */
    GenVertex * thePVertex;

  private:

    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    PVertex & operator=(const PVertex &);

  };

}


#endif /* RIVET_PVertex_H */
