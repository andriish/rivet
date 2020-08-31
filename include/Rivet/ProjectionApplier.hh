// -*- C++ -*-
#ifndef RIVET_ProjectionApplier_HH
#define RIVET_ProjectionApplier_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Projection.fhh"
#include "Rivet/ProjectionHandler.hh"
#include "Rivet/Tools/Logging.hh"
#include<queue>

using namespace std; // @@@AK
namespace Rivet {


  // Forward declarations
  class Event;


  /// @brief Common base class for Projection and Analysis, used for internal polymorphism
  ///
  /// Empty interface used for storing Projection and Analysis pointers in the
  /// same container (used by the ProjectionHandler)
  class ProjectionApplier {
  public:

    // The proj handler needs access to reset the _allowProjReg flag before calling a.init()
    // friend class ProjectionHandler;

    /// Constructor
    ProjectionApplier();

    // Virtual destructor: ensure that inheritance is possible.
    virtual ~ProjectionApplier();


    /// @name Metadata functions
    //@{
    /// Get the name of this Projection or Analysis class
    virtual std::string name() const = 0;
    //@}

    /// @name Projection "getting" functions
    //@{

    /// Get the contained projections, including recursion.
    std::set<ConstProjectionPtr> getProjections() const {
      return getProjHandler().getChildProjections(*this, ProjectionHandler::DEEP);
    }

    /// Does this applier have a projection registered under the name @a name?
    bool hasProjection(const std::string& name) const {
      return getProjHandler().hasProjection(*this, name);
    }

    /// Get the named projection, specifying return type via a template argument.
    /// @todo Add SFINAE to require that PROJ inherit from Projection
    template <typename PROJ>
    const PROJ& getProjection(const std::string& name) const {
      const Projection& p = getProjHandler().getProjection(*this, name);
      return pcast<PROJ>(p);
    }
    /// Get the named projection, specifying return type via a template argument (user-facing alias).
    /// @todo Add SFINAE to require that PROJ inherit from Projection
    template <typename PROJ>
    const PROJ& get(const std::string& name) const { return getProjection<PROJ>(name); }

    /// Get the named projection (non-templated, so returns as a reference to a
    /// Projection base class).
    const Projection& getProjection(const std::string& name) const {
      return getProjHandler().getProjection(*this, name);
    }

    //@}


    /// @name Projection applying functions
    //@{

    /// Apply the supplied projection on event @a evt.
    ///
    /// @deprecated Prefer the simpler apply<> form
    template <typename PROJ=Projection>
    typename std::enable_if_t<std::is_base_of<Projection, PROJ>::value, const PROJ&>
    applyProjection(const Event& evt, const Projection& proj) const {
      return pcast<PROJ>(_applyProjection(evt, proj));
    }

    /// Apply the supplied projection on event @a evt (user-facing alias).
    template <typename PROJ=Projection>
    typename std::enable_if_t<std::is_base_of<Projection, PROJ>::value, const PROJ&>
    apply(const Event& evt, const Projection& proj) const { return applyProjection<PROJ>(evt, proj); }


    /// Apply the supplied projection on event @a evt.
    ///
    /// @deprecated Prefer the simpler apply<> form
    template <typename PROJ=Projection>
    typename std::enable_if_t<std::is_base_of<Projection, PROJ>::value, const PROJ&>
    applyProjection(const Event& evt, const PROJ& proj) const {
      return pcast<PROJ>(_applyProjection(evt, proj));
    }

    /// Apply the supplied projection on event @a evt (user-facing alias).
    template <typename PROJ=Projection>
    typename std::enable_if_t<std::is_base_of<Projection, PROJ>::value, const PROJ&>
    apply(const Event& evt, const PROJ& proj) const { return applyProjection<PROJ>(evt, proj); }


    /// Apply the named projection on event @a evt.
    ///
    /// @deprecated Prefer the simpler apply<> form
    template <typename PROJ=Projection>
    typename std::enable_if_t<std::is_base_of<Projection, PROJ>::value, const PROJ&>
    applyProjection(const Event& evt, const std::string& name) const {
      return pcast<PROJ>(_applyProjection(evt, name));
    }

    /// Apply the supplied projection on event @a evt (user-facing alias).
    template <typename PROJ=Projection>
    typename std::enable_if_t<std::is_base_of<Projection, PROJ>::value, const PROJ&>
    apply(const Event& evt, const std::string& name) const { return applyProjection<PROJ>(evt, name); }

    /// Apply the supplied projection on event @a evt (convenience arg-reordering alias).
    template <typename PROJ=Projection>
    typename std::enable_if_t<std::is_base_of<Projection, PROJ>::value, const PROJ&>
    apply(const std::string& name, const Event& evt) const { return applyProjection<PROJ>(evt, name); }

    //@}


    /// Mark this object as owned by a proj-handler
    void markAsOwned() const { _owned = true; }


  protected:

    Log& getLog() const {
      return Log::getLog("Rivet.ProjectionHandler");
    }


    /// Get a reference to the ProjectionHandler for this thread.
    ProjectionHandler& getProjHandler() const {
      return *_projhandler;
    }


    /// @name Projection registration functions
    //@{

    /// @brief Register a contained projection
    ///
    /// The type of the argument is used to instantiate a new projection
    /// internally: this new object is applied to events rather than the
    /// argument object. Hence you are advised to only use locally-scoped
    /// Projection objects in your Projection and Analysis constructors, and to
    /// avoid polymorphism (e.g. handling @c ConcreteProjection via a pointer or
    /// reference to type @c Projection) since this will screw up the internal
    /// type management.
    ///
    /// @todo Add SFINAE to require that PROJ inherit from Projection
    template <typename PROJ>
    const PROJ& declareProjection(const PROJ& proj, const std::string& name) {
      const Projection& reg = _declareProjection(proj, name);
      const PROJ& rtn = dynamic_cast<const PROJ&>(reg);
      return rtn;
    }

    /// @brief Register a contained projection (user-facing version)
    /// @todo Add SFINAE to require that PROJ inherit from Projection
    template <typename PROJ>
    const PROJ& declare(const PROJ& proj, const std::string& name) {
      //{@@@AK
      if(_projhandler){
        return declareProjection(proj, name);
      }
      std::unique_ptr<Projection> projClone = proj.clone();
      _declQueue.push(make_pair(projClone.get(), name));
      return proj;
    }
    /// @brief Register a contained projection (user-facing, arg-reordered version)
    /// @todo Add SFINAE to require that PROJ inherit from Projection
    template <typename PROJ>
    const PROJ& declare(const std::string& name, const PROJ& proj) {
      if(_projhandler){
        return declareProjection(proj, name);
      }
      std::unique_ptr<Projection> projClone = proj.clone();
      _declQueue.push(make_pair(projClone.get(), name));
      return proj;
    }
//@@@AK}

    /// Untemplated function to do the work...
    const Projection& _declareProjection(const Projection& proj, const std::string& name);

    //@}


    /// Non-templated version of string-based applyProjection, to work around
    /// header dependency issue.
    const Projection& _applyProjection(const Event& evt, const std::string& name) const;

    /// Non-templated version of proj-based applyProjection, to work around
    /// header dependency issue.
    const Projection& _applyProjection(const Event& evt, const Projection& proj) const;

///{ @@@AK
    /// @todo AB: Add Doxygen comment, follow surrounding coding style
    /// @todo AB: Changes to this header force a rebuild of everything: put the implementation in the .cc file
    void setProjectionHandler(ProjectionHandler& projectionHandler);
// @@@AK}

    /// Flag to forbid projection registration in analyses until the init phase
    bool _allowProjReg;


  private:

    /// Mark object as owned by the _projhandler
    mutable bool _owned;

    /// Pointer to projection handler.
    /// @todo AB: I don't think this can work: what's the null value on construction?
    //ProjectionHandler& _projhandler;
    ProjectionHandler* _projhandler;
    //shared_ptr<ProjectionHandler&> _projhandler;

    /// @todo AB: You can't store references... how does this work????
    /// @todo AB: What's the string for? - declare receives reference to a Projection and name, so we need to store both as long as we delays the declareProjection
    std::queue<pair<Projection*, string>> _declQueue; // @@@AK

  };


}

#endif
