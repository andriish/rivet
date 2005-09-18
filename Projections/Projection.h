// -*- C++ -*-
#ifndef RIVET_Projection_H
#define RIVET_Projection_H
//
// This is the declaration of the Projection class.
//

// #include "Rivet/Config/Rivet.h"
#include "Projection.fh"
#include "Rivet/Analysis/Event.fh"

namespace Rivet {

/**
 * Projection is the base class of all Projections to be used by
 * Rivet. A Projection object can be assigned to an Event object and
 * will then define a processed part of the information available in
 * the Event, which then can be used by other Projection objects
 * and/or Analysis objects.
 *
 * The main virtual functions to be overridden by concrete sub-classes
 * are project(const Event &) and cmp(const Projection &).
 *
 */
class Projection {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Projection();

  /**
   * The copy constructor.
   */
  inline Projection(const Projection &);

  /**
   * The destructor.
   */
  virtual ~Projection();
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
  int cmp(const Projection & p) const;

public:

  
  /**
   * Determine whether this object should be ordered before the object
   * \a p given as argument. If \a p is of a different class than
   * this, the before() function of the corresponding type_info
   * objects is used. Otherwise, if the objects are of the same class,
   * the virtual cmp(const Projection &) will be returned.
   */
  inline bool before(const Projection & p) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Projection & operator=(const Projection &);

};

}

namespace std {
template <>
struct less<const Rivet::Projection *> :
    public binary_function<const Rivet::Projection *,
                           const Rivet::Projection *, bool> 
{
  /**
   * This is the function called when comparing two pointers to
   * Rivet::Projection.
   */
  bool operator()(const Rivet::Projection * x,
		  const Rivet::Projection * y) const {
    return x->before(*y);
  }
};

}

#endif /* RIVET_Projection_H */
