// -*- C++ -*-
#ifndef RIVET_InvMassFinalState_HH
#define RIVET_InvMassFinalState_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {

  /// Project only charged final state particles.
  class InvMassFinalState : public FinalState {

  public:
    
    /// Constructor: the supplied FinalState projection is assumed to live through the run.
    // constructor for a single pair
  InvMassFinalState(FinalState& fsp,
                      std::pair<long, long> ids,  //vector of pairs of decay products
                      double minmass,     //min inv mass
                      double maxmass,     //max inv mass
                      double mineta = -MaxRapidity,
                      double maxeta = MaxRapidity,
                      double minpt = 0.0*GeV)
    { 
      vector<pair<long,long> > vids;
      vids.push_back(ids);
      InvMassFinalState(fsp,vids,minmass,maxmass,mineta,maxeta,minpt);   
    }

    InvMassFinalState(FinalState& fsp,
                      std::vector<std::pair<long, long> > ids,  //vector of pairs of decay products
                      double minmass,     //min inv mass
                      double maxmass,     //max inv mass
                      double mineta = -MaxRapidity,
                      double maxeta = MaxRapidity,
                      double minpt = 0.0*GeV)
      : FinalState(mineta, maxeta, minpt), _decayids(ids), _minmass(minmass), _maxmass(maxmass)
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
