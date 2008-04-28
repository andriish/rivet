// -*- C++ -*-
#ifndef RIVET_JetShape_HH
#define RIVET_JetShape_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Tools/Utils.hh"

namespace Rivet {


  /// Enums to distinguish between different recombination schemes
  enum Scheme { ENERGY, SNOWMASS };


  /**
     @brief Calculate the jet shape.

     Calculate the differential and integral jet shapes in \f$P_{\perp}\f$ for a given 
     set of jet axes each event.
     
     The recombination scheme (ENERGY or SNOWMASS) has to be specified when invoking the 
     constructor.

     The differential jet shape around a given jet axes at distance interval 
     \f$ r/pm\delta r/2 \f$ is defined as
     \f[
     \rho(r) = 
     \frac{1}{\delta r} \frac{1}{N_\mathrm{jets}}
     \sum_\mathrm{jets} \frac{P_{\perp}(r - \delta r/2, r+\delta r/2)}{P_{\perp}(0,R)}
     \f]
     with \f$ 0 \le r \le R \f$

     The integral jet shape around a given jet axes until distance r is defined as
     \f[
     \Psi(r) = 
       \frac{1}{N_\mathrm{jets}} \sum_\mathrm{jets} 
       \frac{P_{\perp}(0, r)}{P_{\perp}(0,R)} 
     \f]
     with \f$ 0 \le r \le R \f$
     
     The constructor expects also the equidistant binning in radius r to produce the
     jet shape of all bins in a vector and this separately for each jet to allow
     post selection.

     Internally, this projection uses the VetoedFinalState projection to determine the
     jet shapes around the jet axes.

     The jet axes are passed for each event.
  */
  class JetShape: public Projection {
    
  public:

    /// Constructor.
    JetShape(const VetoedFinalState& vfsp, const vector<FourMomentum>& jetaxes, 
             double rmin=0.0, double rmax=0.7, double interval=0.1, 
             double r1minPsi=0.3, Scheme distscheme=ENERGY)
      : _jetaxes(jetaxes), 
        _rmin(rmin), _rmax(rmax), 
        _interval(interval), _r1minPsi(r1minPsi), _distscheme(distscheme)
    { 
      setName("JetShape");
      _nbins = int(round((rmax-rmin)/interval));
      addProjection(vfsp, "FS");
    }

    
  public:
    
    /// Return number of equidistant radius bins.
    double getNbins() const {
      return _nbins;
    }
    
    /// Return \f$ r_\text{min} \f$ value.
    double getRmin() const {
      return _rmin;
    }
    
    /// Return \f$ r_\text{max} \f$ value.
    double getRmax() const {
      return _rmax;
    }
    
    /// Return Rad interval size.
    double getInterval() const {
      return _interval;
    }


    /// Return value of differential jet shape profile histo bin.
    double getDiffJetShape(size_t pTbin, size_t rbin) const {
      return _diffjetshapes[pTbin][rbin];
    }
    
    /// Return value of integrated jet shape profile histo bin.
    double getIntJetShape(size_t pTbin, size_t rbin) const {
      return _intjetshapes[pTbin][rbin];
    }
    
    /// Return value of \f$ \Psi \f$ (integrated jet shape) at given radius for a \f$ p_T \f$ bin.
    double getPsi(size_t pTbin) const {
      return _PsiSlot[pTbin];
    }
    

  protected:
    
    /// Apply the projection to the event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;
 
       
  private:

    ///The jet axes of the jet algorithm projection
    //vector<FourMomentum> _jetaxes;
    const vector<FourMomentum>& _jetaxes;

    /// @name The projected jet shapes
    /// @{
    ///Output jets vector containg jet shape vectors
    vector<vector<double> >  _diffjetshapes;
    vector<vector<double> >  _intjetshapes;
    vector<double> _PsiSlot;
    /// @}

    ///Jet shape parameters
    ///min radius (typycally r=0)
    double _rmin;
    ///max radius
    double _rmax;
    ///length of radius interval
    double _interval;
    ///One minus Psi radius
    double _r1minPsi;
    ///ENERGY or SNOWMASS recombination scheme
    Scheme _distscheme;
    size_t _nbins;
  };

  
}

#endif
