// -*- C++ -*-
#ifndef RIVET_JetShape_HH
#define RIVET_JetShape_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/RivetCLHEP.hh"
#include "Rivet/Event.hh"

/// @todo: we should have a global math library to define 
/// such things as Pi, deltaPhi, deltaR, etc.
/// @todo Remove use of D0-specific code.
/// The inline_maths.h code is not D0 specific, but globally useful
#include "Rivet/Tools/D0RunIIcone/inline_maths.h"


namespace Rivet {

  /// Enums to distinguish between different recombination schemes
  enum schemelist {ENERGY, SNOWMASS};


  /**
     @brief Calculate the jetshape.

     Calculate the differential and integral jet shapes in \f$P_{\perp}\f$ for a given 
     set of jet axes each event.
     
     The recombination scheme (ENERGY or SNOWMASS) has to be specified when invoking the 
     constructor.

     The differential jet shape around a given jet axes at distance interval 
     \f$ r/pm\delta r/2 \f$ is defined as
     \f$ \rho(r)=\frac{1}{\delta r}\frac{1}{N_{jets}}
     \sum_jets\frac{P_{\perp}(r-\delta r/2, r+\delta r/2)}{P_{\perp}(0,R)} \f$
     with \f$ 0\lesseq r \lesseq R \f$

     The integral jet shape around a given jet axes until distance r is defined as
     \f$ \Psi(r)=\frac{1}{N_{jets}}\sum_jets\frac{P_{\perp}(0, r)}{P_{\perp}(0,R)} \f$
     with \f$ 0\lesseq r \lesseq R \f$
     
     The constructor expects also the equidistant binning in radius r to produce the
     jet shape of all bins in a vector and this separately for each jet to allow
     post selection.
     Internally, this projection uses the Vetoed Final State projection to determine the
     jet shapes around the jet axes.
     The jet axes are passed for each event.
  */
  class JetShape: public Projection {
    
  public:


    /// Constructor. The provided FinalState projection must live throughout the run.
    inline JetShape(VetoedFinalState& vfsp, vector<LorentzVector>& jetaxes, 
		    vector<vector<double> >&  diffjetshapes,
		    vector<vector<double> >&  intjetshapes,
		    vector<double>& oneminPsiShape,
		    double rmin=0.0, double rmax=0.7, double interval=0.1, 
		    double r1minPsi=0.3, schemelist distscheme=ENERGY)
      : _vfsproj(vfsp), _jetaxes(jetaxes), _diffjetshapes(diffjetshapes),
	_intjetshapes(intjetshapes), _oneminPsiShape(oneminPsiShape), _rmin(rmin), _rmax(rmax), 
	_interval(interval), _r1minPsi(r1minPsi), _distscheme(distscheme)
    { 
      _nbins = int(round((rmax-rmin)/interval));
      addProjection(_vfsproj);
    }
    
  public:
    /// Return the name of the projection
    inline string getName() const {
      return "JetShape";
    }
    
    /// Return number of equidistant radius bins
    inline double getNbins() const {
      return _nbins;
    }
    
    /// Return rmin value
    inline double getRmin() const {
      return _rmin;
    }
    
    /// Return rmax value
    inline double getRmax() const {
      return _rmax;
    }
    
    /// Return Rad interval size
    inline double getInterval() const {
      return _interval;
    }
    

  protected:
    
    /// Apply the projection to the event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;
 
       
  private:
        
    /// The VetoedFinalState projection used by this projection
    VetoedFinalState _vfsproj;

    ///The jet axes of the jet algorithm projection
    vector<LorentzVector>& _jetaxes;

    /// @name The projected jet shapes
    /// @{
    ///Output jets vector containg jet shape vectors
    vector<vector<double> >&  _diffjetshapes;
    vector<vector<double> >&  _intjetshapes;
    vector<double>& _oneminPsiShape;
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
    schemelist _distscheme;
    int _nbins;

    
  };
  
}

#endif
