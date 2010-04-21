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


  /// @brief Randomly lose a fraction of the particles from the supplied final state projection.
  class LossyFinalState : public FinalState {
  public:

    /// @name Constructors
    //@{

    /// Constructor from FinalState.
    LossyFinalState(const FinalState& fsp, double lossfraction)
      : _lossFraction(lossfraction)
    {
      setName("LossyFinalState");
      addProjection(fsp, "FS");
      assert(_lossFraction >= 0);
    }

    /// Stand-alone constructor. Initialises the base FinalState projection.
    LossyFinalState(double lossfraction,
                    double mineta = -MAXRAPIDITY,
                    double maxeta = MAXRAPIDITY,
                    double minpt = 0.0)
      : _lossFraction(lossfraction)
    {
      setName("LossyFinalState");
      addProjection(FinalState(mineta, maxeta, minpt), "FS");
      assert(_lossFraction >= 0);
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new LossyFinalState(*this);
    }

    //@}


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

    /// Fraction of particles to lose.
    const double _lossFraction;

  };


}

#endif
