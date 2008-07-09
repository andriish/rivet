// -*- C++ -*-
#ifndef RIVET_ProjectionApplier_HH
#define RIVET_ProjectionApplier_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Event.fhh"
#include "Rivet/Projection.fhh"
#include "Rivet/ProjectionApplier.fhh"
#include "Rivet/ProjectionHandler.hh"
#include "Rivet/Tools/Logging.hh"


namespace Rivet {



  /// Empty interface used for storing Projection and Analysis pointers in the
  /// same container (used by the ProjectionHandler)
  class ProjectionApplier {
  public:
    ProjectionApplier();

    // Ensure that inheritance is possible.
    virtual ~ProjectionApplier() { }


  public:

    /// @name Metadata functions
    //@{
    /// Get the name of this Projection or Analysis class
    virtual string getName() const = 0;
    //@}

    /// @name Projection "getting" functions
    //@{
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
    //@}


    /// @name Projection applying functions 
    //@{
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
    //@}
   

  protected:

    Log& getLog() const {
      return Log::getLog("Rivet.ProjectionHandler");
    }

    /// Get a reference to the ProjectionHandler for this thread.
    ProjectionHandler& getProjHandler() const {
      assert(_projhandler);
      return *_projhandler;
    }


  protected:


    /// @name Projection registration functions 
    //@{

    /// Register a contained projection. The type of the argument is used to
    /// instantiate a new projection internally: this new object is applied to 
    /// events rather than the argument object. Hence you are advised to only use 
    /// locally-scoped Projection objects in your Projection and Analysis
    /// constructors, and to avoid polymorphism (e.g. handling @c ConcreteProjection
    /// via a pointer or reference to type @c Projection) since this will screw
    /// up the internal type management.
    template <typename PROJ>
    const PROJ& addProjection(const PROJ& proj, const string& name) {
      getLog() << Log::TRACE << "Cloning projection " << proj.getName() << endl;
      const Projection* newpproj = proj.clone();
      getLog() << Log::TRACE << "Cloned projection " << proj.getName() << " at " << newpproj << endl;
      const Projection* reg = getProjHandler().registerClonedProjection(*this, &proj, newpproj, name);
      assert(reg);
      if (reg != newpproj) delete newpproj;
      return dynamic_cast<const PROJ&>(*reg);
    }

    //@}
    
    
  private:
    
    /// Non-templated version of string-based applyProjection, to work around
    /// header dependency issue.
    const Projection& _applyProjection(const Event& evt, const string& name) const;
    
    /// Non-templated version of proj-based applyProjection, to work around
    /// header dependency issue.
    const Projection& _applyProjection(const Event& evt, const Projection& proj) const;
    
    
  private:
    
    /// Pointer to projection handler.
    ProjectionHandler* _projhandler;
    
  };

}

#endif
