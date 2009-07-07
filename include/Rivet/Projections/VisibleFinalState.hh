// -*- C++ -*-
#ifndef RIVET_VisibleFinalState_HH
#define RIVET_VisibleFinalState_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// Final state modifier which excludes any particles which are not experimentally visible
  class VisibleFinalState : public FinalState {

  public:
    
    /// @name Constructors
    //@{
    /// Default constructor.
    VisibleFinalState() {
      setName("VisibleFinalState");
      VetoedFinalState vfs;
      vfs.vetoNeutrinos();
      addProjection(vfs, "VFS");
    }

    /// Constructor with min and max pseudorapidity \f$ \eta \f$ and min \f$ p_T
    /// \f$ (in GeV).
    VisibleFinalState(double mineta = -MaxRapidity,
                      double maxeta =  MaxRapidity,
                      double minpt  =  0.0*GeV) {
      setName("VisibleFinalState");
      VetoedFinalState vfs(FinalState(mineta, maxeta, minpt));
      vfs.vetoNeutrinos();
      addProjection(vfs, "VFS");
    }

    /// Constructor with specific FinalState.
    VisibleFinalState(const FinalState& fsp)
    {
      setName("VisibleFinalState");
      VetoedFinalState vfs(fsp);
      vfs.vetoNeutrinos();
      addProjection(vfs, "VFS");
    }


    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new VisibleFinalState(*this);
    }
    //@}
    

  protected:
    
    /// Apply the projection on the supplied event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;

  };

  
}


#endif
