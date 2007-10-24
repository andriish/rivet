// -*- C++ -*-
#ifndef RIVET_Cuts_HH
#define RIVET_Cuts_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Cuts.fhh"
#include <iostream>


namespace Rivet {


  /// A minimal wrapper class to define >=, <= cut pairs.
  class BinaryCut {
  public:
    /// Default constructor.
    BinaryCut() {
      _raw = pair<double, double>(-numeric_limits<double>::max(), numeric_limits<double>::max());
    }
    
    /// Valued constructor.
    BinaryCut(double higherthan, double lowerthan) {
      _raw = pair<double, double>(higherthan, lowerthan);
    }
    
    /// @name Accessing the low and high cut values
    //@{
    BinaryCut& setLowerThan(double val) { _raw.second = val; return *this; }
    BinaryCut& setHigherThan(double val) { _raw.first = val; return *this; }
    double& getLowerThan() { return _raw.second; }
    double& getHigherThan() { return _raw.first; }
    const double& getLowerThan() const { return _raw.second; }
    const double& getHigherThan() const { return _raw.first; }
    //@}
    
  private:
    pair<double, double> _raw;
  };

  inline bool operator==(BinaryCut a, BinaryCut b) {
    return a.getLowerThan() == b.getLowerThan() && a.getHigherThan() == b.getHigherThan();
  }  

  inline bool operator<(BinaryCut a, BinaryCut b) {
    const bool lowersEqual = a.getLowerThan() == b.getLowerThan();
    if (!lowersEqual) return a.getLowerThan() < b.getLowerThan();
    return a.getHigherThan() < b.getHigherThan();
  }  

  inline bool operator>(BinaryCut a, BinaryCut b) {
    const bool lowersEqual = a.getLowerThan() == b.getLowerThan();
    if (!lowersEqual) return a.getLowerThan() > b.getLowerThan();
    return a.getHigherThan() > b.getHigherThan();
  }  




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
    Cuts() { }
    //@}

  public:

    /// Define this cut by a quantity to be constrained and the comparison type of the constraint.
    Cuts& addCut(const string& quantity, const Comparison& comparison, const double value);


    /// Combine with another set of Cuts, using addCut iternally.
    Cuts& addCuts(const Cuts& other) {
      for (const_iterator cut = other.begin(); cut != other.end(); ++cut) {
        addCut(cut->first, MORE_EQ, cut->second.getHigherThan());
        addCut(cut->first, LESS_EQ, cut->second.getLowerThan());
      }
      return *this;
    }

    /// Set the value of this cut, bypassing the combination mechanism.
    // Cuts& setCut(const string& quantity, const BinaryCut& bincut) {
    //   _cuts[quantity] = bincut;
    //   return *this;
    // }

    /// Make sure that the cuts are internally consistent. Returns no
    /// value, but an exception will be thrown if there is an inconsistency.
    bool checkConsistency() const;

    /// Print the parameters to the given \a stream.
    ostream& print(ostream& stream) const;


  public:
    /// Typedef for a named collection of binary cut objects.
    typedef map<string, BinaryCut> NamedBinaryCuts;


  public:

    /// @name Non-const iterators over the cuts
    //@{
    typedef NamedBinaryCuts::iterator iterator;
    iterator begin() { return _cuts.begin(); }
    iterator end() { return _cuts.end(); }
    iterator find(const string& quantity) { return _cuts.find(quantity); }
    //@}

    /// @name Non-const iterators over the cuts
    //@{
    typedef NamedBinaryCuts::const_iterator const_iterator;
    const const_iterator begin() const { return _cuts.begin(); }
    const const_iterator end() const { return _cuts.end(); }
    const const_iterator find(const string& quantity) const { return _cuts.find(quantity); }
    //@}


  private:
    NamedBinaryCuts _cuts;

  };


  /// Allow Cuts to be passed to an ostream.
  inline ostream& operator<<(ostream& os, const Cuts& cuts) {
    return cuts.print(os);
  }

}

#endif
