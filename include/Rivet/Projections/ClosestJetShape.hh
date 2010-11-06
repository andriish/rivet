// -*- C++ -*-
#ifndef RIVET_ClosestJetShape_HH
#define RIVET_ClosestJetShape_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Tools/Utils.hh"

namespace Rivet {


  /**
     @brief Calculate the jet shape.

     Calculate the differential and integral jet shapes in \f$P_{\perp}\f$ for a given
     set of jet axes each event.

     The rapidity scheme (\f$ \eta \f$ or \f$ y \f$) has to be specified when
     invoking the constructor.

     The differential jet shape around a given jet axis at distance interval
     \f$ r \pm \delta{r}/2 \f$ is defined as
     \f[
     \rho(r) =
       \frac{1}{\delta r} \frac{1}{N_\mathrm{jets}}
       \sum_\mathrm{jets} \frac{P_\perp(r - \delta r/2, r+\delta r/2)}{p_\perp(0, R)}
     \f]
     with \f$ 0 \le r \le R \f$ and \f$ P_\perp(r_1, r_2) = \sum_{\in [r_1, r_2)} p_\perp \f$.

     The integral jet shape around a given jet axes until distance \f$ r \f$ is defined as
     \f[
     \Psi(r) =
       \frac{1}{N_\mathrm{jets}}
       \sum_\mathrm{jets} \frac{P_\perp(0, r)}{p_\perp(0, R)}
     \f]
     with \f$ 0 \le r \le R \f$ and \f$ P_\perp(r_1, r_2) = \sum_{\in [r_1, r_2)} p_\perp \f$.

     The constructor expects also the equidistant binning in radius \f$ r \f$ to produce the
     jet shape of all bins in a vector and this separately for each jet to allow
     post-selection.

     In this implementation, the jet axes are passed for each event.

     @deprecated The closest-axis jet shape algorithm is incorrect and should not be used.
  */
  class ClosestJetShape : public Projection {

    /// @todo Review: remove external jet axes, binning, etc.

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor.
    ClosestJetShape(const FinalState& fs, const vector<FourMomentum>& jetaxes,
                    double rmin=0.0, double rmax=0.7, double interval=0.1,
                    double r1minPsi=0.3, DeltaRScheme distscheme=RAPIDITY);

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new ClosestJetShape(*this);
    }

    //@}


    /// Reset projection between events
    void clear();


  public:


    /// Number of equidistant radius bins.
    double numBins() const {
      return _nbins;
    }

    /// \f$ r_\text{min} \f$ value.
    double rMin() const {
      return _rmin;
    }

    /// \f$ r_\text{max} \f$ value.
    double rMax() const {
      return _rmax;
    }

    /// Radius interval size.
    double interval() const {
      return _interval;
    }

    /// Return value of differential jet shape profile histo bin.
    /// @todo Remove this external indexing thing
    double diffJetShape(size_t pTbin, size_t rbin) const {
      return _diffjetshapes[pTbin][rbin];
    }

    /// Return value of integrated jet shape profile histo bin.
    /// @todo Remove this external indexing thing
    double intJetShape(size_t pTbin, size_t rbin) const {
      return _intjetshapes[pTbin][rbin];
    }

    /// Return value of \f$ \Psi \f$ (integrated jet shape) at given radius for a \f$ p_T \f$ bin.
    /// @todo Remove this external indexing thing
    double psi(size_t pTbin) const {
      return _PsiSlot[pTbin];
    }


  protected:

    /// Apply the projection to the event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;


  private:

    /// The jet axes of the jet algorithm projection
    const vector<FourMomentum>& _jetaxes;


    /// @name The projected jet shapes
    //@{

    /// Jet shape histo
    vector<vector<double> > _diffjetshapes;
    vector<vector<double> > _intjetshapes;
    vector<double> _PsiSlot;

    //@}


    /// @name Jet shape parameters
    //@{

    /// Min radius (typically r=0)
    double _rmin;

    /// Max radius
    double _rmax;

    /// Length of radius interval
    double _interval;

    /// One minus Psi radius
    double _r1minPsi;

    /// Rapidity scheme
    DeltaRScheme _distscheme;

    /// Number of bins
    size_t _nbins;

    //@}
  };


}

#endif
