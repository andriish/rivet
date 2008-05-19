// -*- C++ -*-
#ifndef RIVET_ProjectionHandler_HH
#define RIVET_ProjectionHandler_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.fhh"
#include "Rivet/ProjectionApplier.fhh"
#include "Rivet/Projection.fhh"


namespace Rivet {

  class ProjectionApplier;

  /// @brief The projection handler is a central repository for projections to be used
  /// in a Rivet analysis run.
  ///
  /// Without centralised projections, it can be hard to know which of an
  /// equivalent set of projections will be run on a particular event. In turn,
  /// this may mean that certain projections in the chain can go out of scope
  /// unexpectedly. There were originally also the issues that projections may
  /// need to be held as member pointers to an abstract base class, since
  /// post-construction setup is needed; that projections contained pointers to
  /// their own dependency chain, which could go out of scope; and that
  /// projection members could be modified after being applied to an event
  /// which, due to the caching model, would have unpredictable consequences.
  ///
  /// By centralising all the projections, these issues are eliminated, as well
  /// as allowing analysis classes to contain fewer data members (since
  /// projections are now better accessed by name than by storing a data member
  /// reference or pointer).
  /// 
  /// The core of the ProjectionHandler design is that it is a singleton class,
  /// essentially a wrapper around a map of @c Projection*, indexed by a hash of
  /// the registering object and its local name for the registered projection.
  ///
  // To avoid the need for Projections to implement a copy-constructor,
  // copy-assignment operator etc. (cf. the Big Three rule), any projection
  // which must be set up after construction
  //
  /// TODO:
  /// @todo Thread safety... do we need to be able to have multiple proj handlers, one per thread?
  /// @todo Replace Projection* containers in Projection and Analysis classes with calls to this.
  /// @todo Replace registerProjection calls in Projection and Analysis classes with calls to this.
  class ProjectionHandler {
  private:
    
    /// @name Construction. */
    //@{
    /// The standard constructor.
    ProjectionHandler() { }
    //@}

    /// Singleton instance
    /// @todo Threading?
    static ProjectionHandler* _instance;


  public:
    /// Singleton creation function
    static ProjectionHandler* create();


  public:
    /// @name Projection registration. */
    //@{
    /// Attach and retrieve a projection as a reference.
    const Projection& registerProjection(const ProjectionApplier& parent, 
                                         const Projection& proj, const string& name);

    /// Attach and retrieve a projection as a pointer.
    const Projection* registerProjection(const ProjectionApplier& parent, 
                                         const Projection* proj, const string& name) {
      if (!proj) return 0;
      return &registerProjection(parent, *proj, name);
    }
    //@}


    /// Enum to specify depth of projection search.
    enum ProjDepth { SHALLOW, DEEP };


    /// @name Projection retrieval. */
    //@{
    /// Retrieve a named projection for the given parent. Returning as a
    /// reference is partly to discourage ProjectionApplier classes from storing
    /// pointer members to the registered projections, since that can lead to
    /// problems and there is no need to do so.
    const Projection& getProjection(const ProjectionApplier& parent,
                                    const string& name) const;
    
    /// Get child projections for the given parent. By default this will just
    /// return the projections directly contained by the @a parent, but the @a
    /// depth argument can be changed to do a deep retrieval, which will recurse
    /// through the whole projection chain. In this case, there is no protection
    /// against getting stuck in a circular projection dependency loop.
    set<const Projection*> getChildProjections(const ProjectionApplier& parent,
                                               ProjDepth depth=SHALLOW) const;
    //@}

    /// Projection clearing method: deletes all known projections and empties
    /// the reference collections.
    void clear();


  private:

    /// Get a logger.
    Log& getLog() const;

    /// Typedef for Projection pointer, to allow conversion to a smart pointer in this context.
    typedef const Projection* ProjHandle;

    /// Typedef for a vector of Projection pointers.
    typedef vector<ProjHandle> ProjHandles;

    /// @brief Typedef for the structure used to contain named projections for a
    /// particular containing Analysis or Projection.
    /// @todo Use a shared_pointer class?
    typedef map<const string, ProjHandle> NamedProjs;

    /// Structure used to map a containing Analysis or Projection to its set of
    /// contained projections.
    typedef map<const ProjectionApplier*, NamedProjs> NamedProjsMap;

    /// Core data member, associating a given containing class (via a
    /// ProjectionApplier pointer) to its contained projections.
    NamedProjsMap _namedprojs;

    /// Cache of {@link Projection}s for reverse lookup, to speed up registering
    /// new projections as @c _namedprojs gets large.
    ProjHandles _projs;

  private:

    /// Private destructor means no inheritance from this class.
    ~ProjectionHandler();

    /// The assignment operator is hidden.
    ProjectionHandler& operator=(const ProjectionHandler&);

    /// The copy constructor is hidden.
    ProjectionHandler(const ProjectionHandler&);
  };


}

#endif
