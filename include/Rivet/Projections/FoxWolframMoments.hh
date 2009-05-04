// -*- C++ -*-
#ifndef RIVET_FoxWolframMoments_HH
#define RIVET_FoxWolframMoments_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/FastJets.hh"

#include <gsl/gsl_sf_legendre.h>

#define MAXMOMENT 7

namespace Rivet {


  /// Project out the total visible energy vector, allowing missing 
  /// \f$ E_T \f$ etc. to be calculated.
  class FoxWolframMoments : public Projection {
    
  public:
    
    /// Constructor.
    FoxWolframMoments(const FinalState& fsp)
    { 
        setName("FoxWolframMoments");
        addProjection(fsp, "FS");
        //addProjection(TotalVisibleMomentum(fsp), "SumET");

        VetoedFinalState vfs(fsp);
        vfs
        .addVetoPairId(NU_E)
        .addVetoPairId(NU_MU)
        .addVetoPairId(NU_TAU)
        .addVetoDetail(MUON, 1.0*GeV, MAXDOUBLE);
        addProjection(vfs, "VFS");
        
        addProjection(FastJets(vfs, FastJets::CDFJETCLU, 0.4), "JetsC4");

        // initialize moments vector
        for ( int i = 0; i < MAXMOMENT ; ++i) {
            _fwmoments.push_back(0.0);
        }
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new FoxWolframMoments(*this);
    }
    
  public:

    /// The projected Fox-Wolfram Moment of order l
      const double getFoxWolframMoment(unsigned int l) const { 
        if ( l < MAXMOMENT )
            return _fwmoments[l]; 
        else return -666.0;
      }
      
  protected:
    
    /// Apply the projection to the event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;
        
  private:
      vector<double> _fwmoments;

  };
  
}


#endif
