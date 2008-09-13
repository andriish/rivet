// -*- C++ -*-
#ifndef RIVET_HistoHandler_HH
#define RIVET_HistoHandler_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.fhh"
#include "Rivet/Analysis.fhh"


namespace Rivet {


  /// @brief The projection handler is a central repository for histograms (and
  /// other analysis stats objects) to be used in a Rivet analysis run. This
  /// eliminates the need for analysis classes to contain large numbers of
  /// histogram pointer members, and allows histograms to be accessed via more
  /// user-friendly names than C++ variable names allow.
  ///
  /// The core of the HistoHandler design is that it is a singleton class,
  /// essentially a wrapper around a map of @c AnalysisObject*, indexed by a
  /// hash of the registering object and its local name for the registered
  /// projection.
  ///
  class HistoHandler {
  private:
    
    /// @name Construction. */
    //@{
    /// The standard constructor.
    HistoHandler() { }
    //@}

    /// Singleton instance
    /// @todo Threading?
    static HistoHandler* _instance;

    /// Private destructor means no inheritance from this class.
    ~HistoHandler();

    /// The assignment operator is hidden.
    HistoHandler& operator=(const HistoHandler&);

    /// The copy constructor is hidden.
    HistoHandler(const HistoHandler&);

  public:
    /// Singleton creation function
    static HistoHandler* create();


    ////////////////////////////////////////////////////////


  public:
    /// @name Projection registration. */
    //@{
    /// Attach and retrieve a projection as a reference.
    const AnalysisObject& registerHisto(const Analysis& parent, 
                                        const AnalysisObject& histo, const string& name);


    /// @name Projection retrieval. */
    //@{
    /// Retrieve a named projection for the given parent. Returning as a
    /// reference is partly to discourage ProjectionApplier classes from storing
    /// pointer members to the registered projections, since that can lead to
    /// problems and there is no need to do so.
    const AnalysisObject& getHisto(const Analysis& parent,
                                   const string& name) const;

    /// Projection clearing method: deletes all known projections and empties
    /// the reference collections.
    void clear();


  private:

    /// Get a logger.
    Log& getLog() const;


  private:

    /// Typedef for Projection pointer, to allow conversion to a smart pointer in this context.
    typedef const AnalysisObject* HistoHandle;

    /// Typedef for a vector of histo pointers.
    typedef vector<HistoHandle> HistoHandles;

    /// @brief Typedef for the structure used to contain named projections for a
    /// particular containing Analysis or Projection.
    /// @todo Use a shared_pointer class?
    typedef map<const string, HistoHandle> NamedHistos;

    /// Structure used to map a containing Analysis or Projection to its set of
    /// contained projections.
    typedef map<const Analysis*, NamedHistos> NamedHistosMap;

    /// Core data member, associating a given containing class (via a
    /// ProjectionApplier pointer) to its contained projections.
    NamedHistosMap _namedhistos;

    /// Cache of histos for reverse lookup, to speed up registering
    /// new projections as @c _namedhistos gets large.
    HistoHandles _histos;

  };


}

#endif
