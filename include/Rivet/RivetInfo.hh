// -*- C++ -*-
#ifndef RIVET_RivetInfo_HH
#define RIVET_RivetInfo_HH

#include "Rivet/Rivet.hh"
#include "Rivet/BeamParticle.hh"
#include <iostream>

namespace Rivet {

  /// RivetInfo contains information which can be passed from the
  /// different Projection and Analysis objects in Rivet to the
  /// outside world. The main purpose of the RivetInfo objects is
  /// to allow applications which use Rivet to determine whether or
  /// not a given analysis "makes sense" on the provided events. An
  /// Analysis' RivetInfo object should be the combination of those
  /// from its projections, plus the additional constraints (experiment
  /// conditions and cuts) specific to that analysis.
  ///
  /// To make the output constraints predictable, the valid quantities
  /// and their comparison operations are defined via enums. If 
  /// additional enums are required for an analysis that you are 
  /// implementing, please put a request in to the Rivet developers 
  /// to extend the Quantity enum.
  class RivetInfo {
    
  public:
    
    enum Comparison {
      LESS, GREATER, EQUAL
    };
    
    inline static const string comparisonToStr(const Comparison& c) {
      switch(c) {
      case LESS:
        return "<";
      case GREATER:
        return ">";
      case EQUAL:
        return "=";
      }
    }

    enum Quantity {
      SQRTS, PT, Q2
    };

    inline static const string quantityToStr(const Quantity& q) {
      switch(q) {
      case SQRTS:
        return "sqrt(s)";
      case PT:
        return "pT";
      case Q2:
        return "Q^2";
      }
    }

  public:

    /// Typedef for a param name and comparison-type pair.
    typedef pair<Quantity, Comparison> ParamKey;

    /// Typedef for a collection of params.
    typedef map<ParamKey, double> Params;

    /// Typedef for a beam particle and momentum.
    typedef pair<BeamParticle, double> Beam;

    /// Typedef for a pair of beam particles.
    typedef pair<Beam, Beam> BeamPair;

    /// Typedef for a set of pairs of beams.
    typedef set<BeamPair> BeamsSet;


  public:

    /** @name Standard constructors, destructors and assignment. */
    //@{
    /// The default constructor.
    inline RivetInfo() {
      addValidBeamPair(ANY, ANY);
    }

    /// The copy constructor.
    //inline RivetInfo(const RivetInfo &);

    /// The destructor.
    //virtual ~RivetInfo();

    /// The assignment operator.
    //RivetInfo & operator=(const RivetInfo &);
    //@}

  public:

    /// Add a parameter, defined by a quantity to be constrained, the comparison type
    /// of the constraint and a value.
    RivetInfo& addParam(Quantity quantity, Comparison comparison, double value);

    /// Get the collection of parameters.
    inline const Params& getParams() const { return _params; }

    /// Add a valid pair of beams, with optional momentum constraints in GeV. 
    /// Negative momenta will be counted as declaring that there is no momentum constraint.
    RivetInfo& addValidBeamPair(BeamParticle b1, BeamParticle b2, double mom1 = -1, double mom2 = -1);

    /// Get the set of allowed beam pairs
    inline const BeamsSet& getValidBeamPairs() const { return _validBeams; }

    /// Append all the parameters from \a inf.
    inline RivetInfo& operator+=(const RivetInfo& inf) {
      _append(inf);
      return *this;
    }

    /// Return a new RivetInfo object with all information in this and
    /// \a inf added together.
    inline RivetInfo operator+(const RivetInfo& inf) const {
      RivetInfo ret(*this);
      ret += inf;
      return ret;
    }

    /// Print the parameters to the given \a stream.
    ostream& print(ostream& stream) const;

  protected:

    /// Check consistency and purge multiple parameter definitions. 
    /// @throw runtime_error if there were conflicting parameters.
    void _check() const;
    
    /// Combine params and beam type constraints and check consistency.
    /// If inconsistent, throw a runtime_error.
    inline void _append(const RivetInfo& inf) {
      _combineParams(inf.getParams());
      _combineValidBeams(inf.getValidBeamPairs());
      _check();
    }

    /// Combine parameters into a new minimal set.
    void _combineParams(const Params& otherps);

    /// Combine beam types into the intersection of the two input sets.
    void _combineValidBeams(const BeamsSet& otherbs);

  private:

    /// Map of parameter names, types and values.
    Params _params;

    /// Set of allowed beam pairs.
    BeamsSet _validBeams;

  };


  /// Allow RivetInfo to be passed to an iostream
  inline ostream& operator<<(ostream& os, const RivetInfo& i) {
    return i.print(os);
  }

}

#endif
