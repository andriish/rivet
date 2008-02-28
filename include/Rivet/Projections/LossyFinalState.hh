// -*- C++ -*-
#ifndef RIVET_LossyFinalState_HH
#define RIVET_LossyFinalState_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {

  /// Randomly lose a fraction of the particles from the supplied final state projection.
  /// @todo This needs an overhaul to make the base projections work properly. Slicing + inheritance again.
  class LossyFinalState : public FinalState {

  public:
    
    /// @name Constructors
    //@{

    /// Constructor from FinalState. The supplied FinalState projection is assumed 
    /// to live through the run.
    /// @todo Does this work? Slicing problem with using "FinalState(fsp)"?
    LossyFinalState(FinalState& fsp, double lossfraction)
      //: FinalState(fsp), _lossFraction(lossfraction)
      : _fsproj(&fsp), _lossFraction(lossfraction)
    { 
      assert(_lossFraction >= 0);
    }
    
    /// Stand-alone constructor. Initialises the base FinalState projection.
    LossyFinalState(double lossfraction,
                    double mineta = -MaxRapidity,
                    double maxeta = MaxRapidity,
                    double minpt = 0.0)
      //: FinalState(mineta, maxeta, minpt), _lossFraction(lossfraction)
      : _fsproj(0), _lossFraction(lossfraction)
    { 
      assert(_lossFraction >= 0);
      _fsproj = new FinalState(mineta, maxeta, minpt);
    }

    /// Return the name of the projection.
    string getName() const {
      return "LossyFinalState";
    }

  protected:
    
    /// Apply the projection on the supplied event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;

  private:

    /// Inner functor used to implement the random lossiness.
    struct RandomFilter {
      RandomFilter(double lossFraction) : _lossFraction(lossFraction) {}
      bool operator()(const Particle& p) {
        /// @todo Use a better RNG
        return (rand()/static_cast<double>(RAND_MAX) < _lossFraction);
      }
      double _lossFraction;
    };

    /// @todo For now, we'll hold a constituent FinalState pointer...
    FinalState* _fsproj;

    /// Fraction of particles to lose.
    const double _lossFraction;
    
  };

  
}


#endif
