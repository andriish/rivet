// -*- C++ -*-
#ifndef RIVET_ProjectionApplier_HH
#define RIVET_ProjectionApplier_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Event.fhh"
#include "Rivet/Projection.fhh"
#include "Rivet/ProjectionApplier.fhh"
#include "Rivet/ProjectionHandler.hh"


namespace Rivet {



  /// Empty interface used for storing Projection and Analysis pointers in the
  /// same container (used by the ProjectionHandler)
  class ProjectionApplier {
  public:
    ProjectionApplier();

    // Ensure that inheritance is possible.
    virtual ~ProjectionApplier() { }


  public:
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
  
  


    // /// Get the contained projections, including recursion.
    // set<ConstProjectionPtr> getProjections() const {
    //   return getProjHandler().getChildProjections(*this, ProjectionHandler::DEEP);
    // }

    // /// Get the named projection, specifying return type via a template argument.
    // template <typename PROJ>
    // const PROJ& getProjection(const string& name) const {
    //   const Projection& p = getProjHandler().getProjection(*this, name);
    //   return pcast<PROJ>(p);
    // }

    // /// Get the named projection (non-templated).
    // const Projection& getProjection(const string& name) const {
    //   return getProjHandler().getProjection(*this, name);
    // }

    // /// Get a reference to the ProjectionHandler for this thread.
    // ProjectionHandler& getProjHandler() const {
    //   assert(_projhandler);
    //   return *_projhandler;
    // }

    // /// Register a contained projection (via reference).
    // template <typename PROJ>
    // const PROJ& addProjection(const PROJ& proj, const string& name) {
    //   getLog() << Log::TRACE << this->getName() << " inserts " 
    //            << proj.getName() << " at: " << &proj << endl;
    //   const Projection& p = getProjHandler().registerProjection(*this, proj, name);
    //   return pcast<PROJ>(p);
    // }

    // /// @todo Discriminate with templated pointer type?
    // // /// Register a contained projection (via pointer).
    // // template <typename PROJ>
    // // const PROJ* addProjection(const PROJ* pproj, const string& name) {
    // //   getLog() << Log::TRACE << this->getName() << " inserts " 
    // //            << pproj->getName() << " at: " << pproj << endl;
    // //   const Projection* p = getProjHandler().registerProjection(*this, pproj, name);
    // //   return dynamic_cast<const PROJ*>(p);
    // // }

    // /// Register a contained projection (via pointer, untemplated).
    // const Projection* addProjection(const Projection* pproj, const string& name) {
    //   getLog() << Log::TRACE << this->getName() << " inserts " 
    //            << pproj->getName() << " at: " << pproj << endl;
    //   return getProjHandler().registerProjection(*this, pproj, name);
    // }

    // /// Apply the named projection on @a event.
    // template <typename PROJ>
    // const PROJ& applyProjection(const Event& evt, const string& name) {
    //   return pcast<PROJ>(evt.applyProjection(getProjection(name)));
    // }

    // /// Apply the supplied projection on @a event.
    // template <typename PROJ>
    // const PROJ& applyProjection(const Event& evt, const PROJ& proj) {
    //   return pcast<PROJ>(evt.applyProjection(proj));
    // }

    // /// Apply the supplied projection on @a event.
    // template <typename PROJ>
    // const PROJ& applyProjection(const Event& evt, const Projection& proj) {
    //   return pcast<PROJ>(evt.applyProjection(proj));
    // }


}

#endif
