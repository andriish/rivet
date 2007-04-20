// -*- C++ -*-
#ifndef RIVET_ParamConstraint_HH
#define RIVET_ParamConstraint_HH

#include "Rivet/Rivet.hh"
#include "Rivet/ParamConstraint.fhh"
#include <iostream>

namespace Rivet {

  /// ParamConstraint contains information which can be passed from the
  /// different Projection and Analysis objects in Rivet to the
  /// outside world. The main purpose of the ParamConstraint objects is
  /// to allow applications which use Rivet to determine whether or
  /// not a given analysis "makes sense" on the provided events. An
  /// Analysis' ParamConstraint object should be the combination of those
  /// from its projections, plus the additional constraints (experiment
  /// conditions and cuts) specific to that analysis.
  ///
  /// To make the output constraints predictable, the valid quantities
  /// and their comparison operations are defined via enums. If 
  /// additional enums are required for an analysis that you are 
  /// implementing, please put a request in to the Rivet developers 
  /// to extend the Quantity enum.
  class ParamConstraint {
    
  public:

    /** @name Standard constructors, destructors and assignment. */
    //@{
    /// The default constructor.
    //inline ParamConstraint() {}

    /// The copy constructor.
    //inline ParamConstraint(const ParamConstraint &);

    /// The destructor.
    //virtual ~ParamConstraint();

    /// The assignment operator.
    //ParamConstraint & operator=(const ParamConstraint &);
    //@}

  public:

    /// Add a parameter, defined by a quantity to be constrained, the comparison type
    /// of the constraint and a value.
    inline ParamConstraint& set(Quantity quantity, Comparison comparison, double value) {
      _quantity = quantity;
      _comp = comparison;
      _value = value;
      return *this;
    }

//    /// Append all the parameters from \a inf.
//     inline ParamConstraint& operator+=(const ParamConstraint& inf) {
//       _append(inf);
//       return *this;
//     }

//     /// Return a new ParamConstraint object with all information in this and
//     /// \a inf added together.
//     inline ParamConstraint operator+(const ParamConstraint& inf) const {
//       ParamConstraint ret(*this);
//       ret += inf;
//       return ret;
//     }

    /// Print the parameters to the given \a stream.
    ostream& print(ostream& stream) const;

  protected:

    /// Check consistency and purge multiple parameter definitions. 
    /// @throw runtime_error if there were conflicting parameters.
    //    void _check() const;
    
    /// Combine params and beam type constraints and check consistency.
    /// If inconsistent, throw a runtime_error.
//     inline void _append(const ParamConstraint& inf) {
//       _combineParams(inf.getParams());
//       _check();
//     }

    /// Combine parameters into a new minimal set.
    //    void _combineParams(const Params& otherps);

  private:

    Quantity _quantity;
    Comparison _comp;
    double _value;

  };


  /// Allow ParamConstraint to be passed to an iostream
  inline ostream& operator<<(ostream& os, const ParamConstraint& p) {
    return p.print(os);
  }

}

#endif
