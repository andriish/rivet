// -*- C++ -*-
#ifndef RIVET_Projection_HH
#define RIVET_Projection_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Projection.fhh"
#include "Rivet/ProjectionApplier.hh"
#include "Rivet/ProjectionHandler.hh"
#include "Rivet/Constraints.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/Event.fhh"
#include "Rivet/Tools/Logging.hh"
//#include "Rivet/Tools/TypeTraits.hh"
#include "Rivet/Cmp.fhh"


namespace std {
  
  /// This is the function called when comparing two (const) pointers to Rivet::Projection.
  template <>
  struct less<const Rivet::Projection*>
    : public binary_function<const Rivet::Projection*, const Rivet::Projection*, bool> {
    bool operator()(const Rivet::Projection* x, const Rivet::Projection* y) const;
  };

}


namespace Rivet {


  /// Convenience method for casting to a const Projection reference.
  template <typename PROJ>
  inline const PROJ& pcast(const Projection& p) {
    return dynamic_cast<const PROJ&>(p);
  }
  
  
  /// Convenience method for casting to a const Projection pointer.
  template <typename PROJ>
  inline const PROJ* pcast(const Projection* p) {
    return dynamic_cast<const PROJ*>(p);
  }


  /// Projection is the base class of all Projections to be used by
  /// Rivet. A Projection object can be assigned to an Event object and
  /// will then define a processed part of the information available in
  /// the Event, which then can be used by other Projection objects
  /// and/or Analysis objects.
  ///
  /// The main virtual functions to be overridden by concrete sub-classes
  /// are project(const Event &) and compare(const Projection &).
  class Projection : public ProjectionApplier {
    
  public:
    
    /// Event is a friend.
    friend class Event;
    
    /// The Cmp specialization for Projection is a friend.
    friend class Cmp<Projection>;
    
  public:
    
    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor.
    Projection();
    
    /// The destructor.
    virtual ~Projection();
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
    ///
    /// By default, this function returns the result of a comparison between two
    /// requested child projections named "FS", since almost all projections
    /// should have such a child.
    virtual int compare(const Projection& p) const;
    
  public:

    /// Determine whether this object should be ordered before the object
    /// \a p given as argument. If \a p is of a different class than
    /// this, the before() function of the corresponding type_info
    /// objects is used. Otherwise, if the objects are of the same class,
    /// the virtual compare(const Projection &) will be returned.
    bool before(const Projection& p) const;
    
    /// Return the Cuts objects for this projection. Derived
    /// classes should ensure that all contained projections are
    /// registered in the @a _projections set for the cut chaining 
    /// to work.
    virtual const Cuts getCuts() const;

    /// Return the BeamConstraints for this projection, not including
    /// recursion. Derived classes should ensure that all contained projections
    /// are registered in the @a _projections set for the beam constraint
    /// chaining to work.
    virtual const set<BeamPair> getBeamPairs() const;

    /// Get the name of the projection.
    virtual string getName() const {
      return _name;
    }

    /// Get the contained projections, including recursion.
    set<ConstProjectionPtr> getProjections() const {
      return getProjHandler().getChildProjections(*this, ProjectionHandler::DEEP);
    }

    /// Get the named projection, specifying return type via a template argument.
    template <typename PROJ>
    const PROJ& getProjection(const string& name) const {
      const Projection& p = getProjHandler().getProjection(*this, name);
      return pcast<PROJ>(p);
    }

    /// Get the named projection (non-templated, so returns as a reference to a
    /// Projection base class).
    const Projection& getProjection(const string& name) const {
      return getProjHandler().getProjection(*this, name);
    }
    

    /// Apply the supplied projection on @a event.
    template <typename PROJ>
    const PROJ& applyProjection(const Event& evt, const PROJ& proj) const {
      return pcast<PROJ>(_applyProjection(evt, proj));
    }


    /// Apply the supplied projection on @a event.
    template <typename PROJ>
    const PROJ& applyProjection(const Event& evt, const Projection& proj) const {
      return pcast<PROJ>(_applyProjection(evt, proj));
    }


    /// Apply the named projection on @a event.
    template <typename PROJ>
    const PROJ& applyProjection(const Event& evt, const string& name) const {
      return pcast<PROJ>(_applyProjection(evt, name));
    }


   
    /// Shortcut to make a named Cmp<Projection> comparison with the @c *this
    /// object automatically passed as one of the parent projections.
    Cmp<Projection> mkNamedPCmp(const Projection& otherparent, const string& pname) const;

    /// Shortcut to make a named Cmp<Projection> comparison with the @c *this
    /// object automatically passed as one of the parent projections.
    Cmp<Projection> mkPCmp(const Projection&, const string& pname) const;

    
  protected:

    /// Get a reference to the ProjectionHandler for this thread.
    ProjectionHandler& getProjHandler() const {
      assert(_projhandler);
      return *_projhandler;
    }

    // /// Register a contained projection.
    // template <typename PROJ>
    // PROJ addProjection(PROJ proj, const string& name) {
    //   // Use traits to determine ref/ptr type of argument:
    //   return _addProjHelper(proj, name, TypeTraits<PROJ>::ArgType());
    // }

    /// Register a contained projection (via reference).
    template <typename PROJ>
    const PROJ& addProjection(const PROJ& proj, const string& name) {
      // getLog() << Log::TRACE << this->getName() << " inserts " 
      //          << proj.getName() << " at: " << &proj << endl;
      return dynamic_cast<const PROJ&>(getProjHandler().registerProjection(*this, proj, name));
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

    /// Used by derived classes to set their name.
    void setName(const string& name) {
      _name = name;
    }

  private:

    /// Non-templated version of string-based applyProjection, to work around
    /// header dependency issue.
    const Projection& _applyProjection(const Event& evt, const string& name) const;
    
    /// Non-templated version of proj-based applyProjection, to work around
    /// header dependency issue.
    const Projection& _applyProjection(const Event& evt, const Projection& proj) const;
    

  private:

    /// Name variable is used by the base class messages to identify
    /// which derived class is being handled.
    string _name;

    /// Parameter constraints.
    Cuts _cuts;
    
    /// Beam-type constraint.
    set<BeamPair> _beamPairs;
    
    /// Pointer to projection handler.
    ProjectionHandler* _projhandler;
  };


}


/// Define "less" operator for Projection* containers in terms of the Projection::before virtual method.
inline bool std::less<const Rivet::Projection *>::operator()(const Rivet::Projection* x, 
                                                             const Rivet::Projection* y) const {
  return x->before(*y);
}


// Definition of the comparison objects and functions
#include "Rivet/Cmp.hh"


#endif
