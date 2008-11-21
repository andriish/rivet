// -*- C++ -*-
#ifndef RIVET_InvMassFinalState_HH
#define RIVET_InvMassFinalState_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// Identify particles which can be paired to make an invariant mass within a given mass window
  class InvMassFinalState : public FinalState {

  public:
    
    // Constructor for a single inv-mass pair
    InvMassFinalState(FinalState& fsp,
                      std::pair<long, long> idpair, // pair of decay products
                      double minmass, // min inv mass
                      double maxmass) // max inv mass
      : _minmass(minmass), _maxmass(maxmass)
    {
      setName("InvMassFinalState");
      addProjection(fsp, "FS");
      _decayids.push_back(idpair);
    }


    InvMassFinalState(FinalState& fsp,
                      std::vector<std::pair<long, long> > idpairs,  // vector of pairs of decay products
                      double minmass, // min inv mass
                      double maxmass) // max inv mass
      : _decayids(idpairs), _minmass(minmass), _maxmass(maxmass)
    { 
      setName("InvMassFinalState");
      addProjection(fsp, "FS");
    }


    /// Clone on the heap.
    virtual const Projection* clone() const {
    	return new InvMassFinalState(*this);
    }
		

  protected:
    
    /// Apply the projection on the supplied event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;


  private:
    
    /// ids of the decay products
    std::vector<std::pair<long, long> > _decayids;
   
    /// min inv mass
    double _minmass;

    /// max inv mass
    double _maxmass;
    
  };

  
}


#endif
