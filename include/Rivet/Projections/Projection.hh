// -*- C++ -*-
#ifndef RIVET_Projection_HH
#define RIVET_Projection_HH

#include "Rivet/Rivet.hh"
#include "Projection.fhh"
#include "Rivet/Constraints.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/Event.fhh"
#include "Rivet/Tools/Logging.fhh"
#include "Rivet/Projections/Cmp.fhh"


namespace Rivet {

  /**
   * Projection is the base class of all Projections to be used by
   * Rivet. A Projection object can be assigned to an Event object and
   * will then define a processed part of the information available in
   * the Event, which then can be used by other Projection objects
   * and/or Analysis objects.
   *
   * The main virtual functions to be overridden by concrete sub-classes
   * are project(const Event &) and compare(const Projection &).
   */
  class Projection {
    
  public:
    
    /// Event is a friend.
    friend class Event;
    
    /// The Cmp specialization for Projection is a friend.
    friend class Cmp<Projection>;
    
  public:
    
    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor.
    inline Projection() { };
    
    /// The destructor.
    inline virtual ~Projection() { };
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
    virtual void project(const Event & e) = 0;
    
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
    virtual int compare(const Projection & p) const = 0;
    
  public:
    
    /**
     * Determine whether this object should be ordered before the object
     * \a p given as argument. If \a p is of a different class than
     * this, the before() function of the corresponding type_info
     * objects is used. Otherwise, if the objects are of the same class,
     * the virtual compare(const Projection &) will be returned.
     */
    inline bool before(const Projection& p) const {
      const std::type_info& thisid = typeid(*this);
      const std::type_info& otherid = typeid(p);
      if (thisid == otherid) {
        return compare(p) < 0;
      } else {
        return thisid.before(otherid);
      }
    }
    
    /// Return the Cuts objects for this projection. Derived
    /// classes should re-implement this function to return the combined
    /// RivetInfo object of this object and of any Projection objects
    /// upon which this depends.
    inline virtual const Cuts getCuts() const {
      return _cuts;
    }

    /// Return the BeamConstraints for this analysis. Derived
    /// classes should re-implement this function to return the combined
    /// allowed beam pairs for this object and for any other Projections
    /// upon which it depends.
    inline virtual const set<BeamPair> getBeamPairs() const {
      return _beamPairs;
    }

    /// Get the name of the projection
    inline virtual string getName() const {
      return "";
    }
    
  protected:

    /// Get a Log object based on the getName() property of the calling projection object.
    Log& getLog();
    
    /// Parameter constraints
    Cuts _cuts;

    /// Beam-type constraint
    set<BeamPair> _beamPairs;


  private:
    
    /// The assignment operator is private and must never be called.
    /// In fact, it should not even be implemented.
    Projection & operator=(const Projection &);
    
  };
  
}


namespace std {
  
  /// This is the function called when comparing two pointers to Rivet::Projection.
  template <>
  struct less<const Rivet::Projection*> 
    : public binary_function<const Rivet::Projection*, const Rivet::Projection*, bool> {
    bool operator()(const Rivet::Projection* x, const Rivet::Projection* y) const {
      return x->before(*y);
    }
  };

}


#endif
