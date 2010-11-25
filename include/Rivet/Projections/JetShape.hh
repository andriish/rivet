// -*- C++ -*-
#ifndef RIVET_JetShape_HH
#define RIVET_JetShape_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Tools/Utils.hh"

namespace Rivet {


  /**
     @brief Calculate the jet shape.

     Calculate the differential and integral jet shapes in \f$P_{\perp}\f$ for a given
     set of jets. This particular jet shape projection calculates jet shapes relative
     to jet centroids, using only the particles associated to each jet, for the hardest
     \f$ n \f$ jets.

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
  */
  class JetShape : public Projection {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor from histo range and number of bins.
    JetShape(const JetAlg& jetalg,
             double rmin, double rmax, size_t nbins,
             double ptmin=0, double ptmax=MAXDOUBLE,
             double absrapmin=-MAXDOUBLE, double absrapmax=-MAXDOUBLE,
             RapScheme rapscheme=RAPIDITY);

    /// Constructor from vector of bin edges.
    JetShape(const JetAlg& jetalg, vector<double> binedges,
             double ptmin=0, double ptmax=MAXDOUBLE,
             double absrapmin=-MAXDOUBLE, double absrapmax=-MAXDOUBLE,
             RapScheme rapscheme=RAPIDITY);

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new JetShape(*this);
    }

    //@}


    /// Reset projection between events.
    void clear();


    /// Do the calculation directly on a supplied collection of Jet objects.
    void calc(const Jets& jets);


  public:


    /// Number of equidistant radius bins.
    size_t numBins() const {
      return _binedges.size() - 1;
    }

    /// \f$ r_\text{min} \f$ value.
    double rMin() const {
      return _binedges.front();
    }

    /// \f$ r_\text{max} \f$ value.
    double rMax() const {
      return _binedges.back();
    }

    /// \f$ p_\perp^\text{min} \f$ value.
    double ptMin() const {
      return _ptcuts.first;
    }

    /// \f$ p_\perp^\text{max} \f$ value.
    double ptMax() const {
      return _ptcuts.second;
    }

    /// Central \f$ r \f$ value for bin @a rbin.
    /// @todo This external indexing thing is a bit nasty...
    double rBinMin(size_t rbin) const {
      assert(inRange(rbin, 0, numBins()));
      return _binedges[rbin];
    }

    /// Central \f$ r \f$ value for bin @a rbin.
    /// @todo This external indexing thing is a bit nasty...
    double rBinMax(size_t rbin) const {
      assert(inRange(rbin, 0, numBins()));
      return _binedges[rbin+1];
    }

    /// Central \f$ r \f$ value for bin @a rbin.
    /// @todo This external indexing thing is a bit nasty...
    double rBinMid(size_t rbin) const {
      assert(inRange(rbin, 0, numBins()));
      return (_binedges[rbin] + _binedges[rbin+1])/2.0;
    }

    /// Return value of differential jet shape profile histo bin.
    /// @todo This external indexing thing is a bit nasty...
    double diffJetShape(size_t rbin) const {
      assert(inRange(rbin, 0, numBins()));
      return _diffjetshapes[rbin];
    }

    /// Return value of integrated jet shape profile histo bin.
    /// @todo This external indexing thing is a bit nasty...
    double intJetShape(size_t rbin) const {
      assert(inRange(rbin, 0, numBins()));
      double rtn  = 0;
      for (size_t i = 0; i <= rbin; ++i) {
        rtn += _diffjetshapes[i];
      }
      return rtn;
    }

    /// @todo Provide int and diff jet shapes with some sort of area normalisation?

    // /// Return value of \f$ \Psi \f$ (integrated jet shape) at given radius for a \f$ p_T \f$ bin.
    // /// @todo Remove this external indexing thing
    // double psi(size_t pTbin) const {
    //   return _PsiSlot[pTbin];
    // }


  protected:

    /// Apply the projection to the event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;


  private:

    /// @name Jet shape parameters
    //@{

    /// Vector of radius bin edges
    vector<double> _binedges;

    /// Lower and upper cuts on contributing jet \f$ p_\perp \f$.
    pair<double, double> _ptcuts;

    /// Lower and upper cuts on contributing jet (pseudo)rapidity.
    pair<double, double> _rapcuts;

    /// Rapidity scheme
    RapScheme _rapscheme;

    //@}


    /// @name The projected jet shapes
    //@{

    /// Jet shape histo
    vector<double> _diffjetshapes;

    //@}

  };


}

#endif
