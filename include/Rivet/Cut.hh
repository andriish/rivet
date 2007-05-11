// -*- C++ -*-
#ifndef RIVET_Cut_HH
#define RIVET_Cut_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Cut.fhh"
#include <iostream>


namespace Rivet {

  /// A Cuts object contains information which can be passed from the
  /// different Projection and Analysis objects in Rivet to the
  /// outside world. The main purpose of the Cut objects is
  /// to allow applications which use Rivet to determine whether or
  /// not a given analysis "makes sense" on the provided events. An
  /// Analysis' Cuts object should be the combination of Cuts
  /// from its projections, plus the additional constraints (experiment
  /// conditions and cuts) specific to that analysis.
  class Cuts {
    
  public:

    /// @name Standard constructors, destructors and assignment.
    //@{
    /// The default constructor.
    inline Cuts() { }
    //@}

  public:

    /// Define this cut by a quantity to be constrained and the comparison type of the constraint.
    Cuts& addCut(const string& quantity, const Comparison& comparison, const double value);

    /// Set the value of this cut, bypassing the combination mechanism.
    // inline Cuts& setCut(const string& quantity, const BinaryCut& bincut) {
    //   _cuts[quantity] = bincut;
    //   return *this;
    // }

    /// Make sure that the cuts are internally consistent. Returns no
    /// value, but an exception will be thrown if there is an inconsistency.
    bool checkConsistency() const;

    /// Print the parameters to the given \a stream.
    ostream& print(ostream& stream) const;


  public:
    /// A minimal wrapper class to define >=, <= cut pairs.
    class BinaryCut {
    public:
      /// Default constructor.
      BinaryCut() {
        _raw = pair<double, double>(numeric_limits<double>::max(), numeric_limits<double>::min());
      }

      /// Valued constructor.
      BinaryCut(double lowerthan, double higherthan) {
        _raw = pair<double, double>(lowerthan, higherthan);
      }

      /// @name Accessing the low and high cut values
      //@{
      inline double& lowerthan() { return _raw.first; }
      inline double& higherthan() { return _raw.second; }
      inline const double& lowerthan() const { return _raw.first; }
      inline const double& higherthan() const { return _raw.second; }
      //@}

    private:
      pair<double, double> _raw;
    };

    /// Typedef for a named collection of binary cut objects.
    typedef map<string, BinaryCut> NamedBinaryCuts;


  public:

    /// @name Non-const iterators over the cuts
    //@{
    typedef NamedBinaryCuts::iterator iterator;
    inline iterator begin() { return _cuts.begin(); }
    inline iterator end() { return _cuts.end(); }
    inline iterator find(const string& quantity) { return _cuts.find(quantity); }
    //@}

    /// @name Non-const iterators over the cuts
    //@{
    typedef NamedBinaryCuts::const_iterator const_iterator;
    inline const const_iterator begin() const { return _cuts.begin(); }
    inline const const_iterator end() const { return _cuts.end(); }
    inline const const_iterator find(const string& quantity) const { return _cuts.find(quantity); }
    //@}


  private:
    NamedBinaryCuts _cuts;

  };


  /// Allow Cut to be passed to an iostream
  inline ostream& operator<<(ostream& os, const Cuts& cuts) {
    return cuts.print(os);
  }

}

#endif
