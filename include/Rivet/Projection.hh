// -*- C++ -*-
#ifndef RIVET_Projection_HH
#define RIVET_Projection_HH

#include "Rivet/Rivet.hh"
#include "Projection.fhh"
#include "Rivet/Constraints.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/Event.fhh"
#include "Rivet/Tools/Logging.hh"


namespace Rivet {

  /// Projection is the base class of all Projections to be used by
  /// Rivet. A Projection object can be assigned to an Event object and
  /// will then define a processed part of the information available in
  /// the Event, which then can be used by other Projection objects
  /// and/or Analysis objects.
  ///
  /// The main virtual functions to be overridden by concrete sub-classes
  /// are project(const Event &) and compare(const Projection &).
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
    Projection() {
      addBeamPair(ANY, ANY);
      getLog() << Log::DEBUG << "Creating " << this->getName() << " at " << this << endl;
    }
    

    /// The destructor.
    virtual ~Projection() {
      getLog() << Log::DEBUG << "Destroying " << this->getName() << " at " << this << endl;
    }
    //@}
    
  protected:
    
    /// Take the information available in the Event and make the
    /// calculations necessary to obtain the projection. Note that this
    /// function must never be called except inside the
    /// Event::applyProjection(Projection *) function. If the information
    /// from other projections are necessary, their project(const Event&)
    /// should not be called, rather the corresponding objects should
    /// be added to the Event using the Event::applyProjection(Projection *)
    /// function.
    virtual void project(const Event& e) = 0;
    

    /// This function is used to define a unique ordering between
    /// different Projection objects of the same class. If this is
    /// considered to be equivalent to the Projector object, \a p, in the
    /// argument the function should return 0. If this object should be
    /// ordered before \a p a negative value should be returned,
    /// otherwise a positive value should be returned. This function must
    /// never be called explicitly, but should only be called from the
    /// operator<(const Projection &). When implementing the function in
    /// concrete sub-classes, it is then guarranteed that the Projection
    /// object \a p in the argument is of the same class as the sub-class
    /// and can be safely dynamically casted to that class.
    ///
    /// When implementing this function in a sub-class, the immediate
    /// base class version of the function should be called first. If the
    /// base class function returns a non-zero value, that value should
    /// be returned immediately. Only if zero is returned should this
    /// function check the member variables of the sub-class to determine
    /// whether this should be ordered before or after \a p, or if it is
    /// equivalent with \a p.
    virtual int compare(const Projection& p) const = 0;
    
  public:
    

    /// Determine whether this object should be ordered before the object
    /// \a p given as argument. If \a p is of a different class than
    /// this, the before() function of the corresponding type_info
    /// objects is used. Otherwise, if the objects are of the same class,
    /// the virtual compare(const Projection &) will be returned.
    bool before(const Projection& p) const {
      const std::type_info& thisid = typeid(*this);
      const std::type_info& otherid = typeid(p);
      if (thisid == otherid) {
        return compare(p) < 0;
      } else {
        return thisid.before(otherid);
      }
    }
    
    /// Return the Cuts objects for this projection. Derived
    /// classes should ensure that all contained projections are
    /// registered in the @a _projections set for the cut chaining 
    /// to work.
    virtual const Cuts getCuts() const {
      Cuts totalCuts = _cuts;
      for (set<ProjectionPtr>::const_iterator p = _projections.begin(); p != _projections.end(); ++p) {
        totalCuts.addCuts((*p)->getCuts());
      }
      return totalCuts;
    }

    /// Return the BeamConstraints for this projection, not including
    /// recursion. Derived classes should ensure that all contained projections
    /// are registered in the @a _projections set for the beam constraint
    /// chaining to work.
    virtual const set<BeamPair> getBeamPairs() const {
      set<BeamPair> ret = _beamPairs;
      for (set<ProjectionPtr>::const_iterator ip = _projections.begin(); ip != _projections.end(); ++ip) {
        ProjectionPtr p = *ip;
        getLog() << Log::DEBUG << "Proj addr = " << p << endl;
        if (p) ret = intersection(ret, p->getBeamPairs());
      }
      return ret;
    }

    /// Get the name of the projection.
    virtual string getName() const {
      return "BaseProjection";
    }

    /// Get the contained projections, including recursion.
    set<ProjectionPtr> getProjections() const {
      set<ProjectionPtr> allProjections = _projections;
      for (set<ProjectionPtr>::const_iterator p = _projections.begin(); p != _projections.end(); ++p) {
        allProjections.insert((*p)->getProjections().begin(), (*p)->getProjections().end());
      }
      return allProjections;
    }
    
  protected:

    /// Add a projection dependency to the projection list.
    Projection& addProjection(Projection& proj) {
      getLog() << Log::DEBUG << " Inserting projection at: " << &proj << endl;
      getLog() << Log::DEBUG << " Inserter/insertee: " << this->getName() << " inserts " << proj.getName() << endl;
      ProjectionPtr pp(& proj);
      _projections.insert(pp);
      return *this;
    }

    /// Add a colliding beam pair.
    Projection& addBeamPair(const ParticleName& beam1, const ParticleName& beam2) {
      _beamPairs.insert(BeamPair(beam1, beam2));
      return *this;
    }

    /// Add a cut.
    Projection& addCut(const string& quantity, const Comparison& comparison, const double value) {
      getLog() << Log::DEBUG << getName() << "::addCut(): " << quantity << " " << comparison << " " << value << endl;
      _cuts.addCut(quantity, comparison, value);
      return *this;
    }

    /// Get a Log object based on the getName() property of the calling projection object.
    Log& getLog() const {
      string logname = "Rivet.Projection." + getName();
      return Log::getLog(logname);
    }
    
    /// Parameter constraints
    Cuts _cuts;

    /// Beam-type constraint
    set<BeamPair> _beamPairs;
 
    /// Collection of pointers to projections, for automatically combining constraints.
    set<ProjectionPtr> _projections;

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
