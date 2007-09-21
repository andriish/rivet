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
#include "Rivet/Tools/D0RunIIcone/inline_maths.h"


namespace Rivet {

  /// @todo Enum consts should be ALL-CAPS
  enum schemelist {energy, snowmass};

  /// Project out the differential and integral \f$ p_T \f$ shapes of jets 
  /// by means of visible final state particles around jet axes.
  class JetShape: public Projection {
    
  public:


    /// Constructor. The provided FinalState projection must live throughout the run.
    inline JetShape(VetoedFinalState& vfsp, double rmin=0.0, double rmax=0.7, 
                    double interval=0.1, double r1minPsi=0.3, schemelist distscheme=energy)
      : _vfsproj(&vfsp), _rmin(rmin), _rmax(rmax), _interval(interval), _r1minPsi(r1minPsi), 
        _distscheme(distscheme)
    { 
      _nbins = int(round((rmax-rmin)/interval));
      addProjection(*_vfsproj);
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
    
    /// Set the list of pre-selected jet axes
    inline void setJetAxes(vector<LorentzVector>& jetaxes) 
    {
      _jetaxes = &jetaxes;
    }
    
    /// @name The projected jet shapes
    /// @{
    //inline const list<vector> getDiffJetShapes() const { return _diffjetshapes; }
    //inline const list<vector> getIntJetShapes() const { return _intjetshapes; }
    inline void setDiffJetShapes(vector<vector<double> >& diffjetshapes) { 
      _diffjetshapes = &diffjetshapes; 
    }

    inline void setIntJetShapes(vector<vector<double> >& intjetshapes) { 
      _intjetshapes = &intjetshapes; 
    }

    inline void setOneMinPsiShape(vector<double>& oneminPsishape) { 
      _oneminPsishape = &oneminPsishape; 
    }
    /// @}

  protected:
    
    /// Apply the projection to the event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;
 
       
  private:
        
    /// The VetoedFinalState projection used by this projection
    /// @todo Are the pointers really needed?
    VetoedFinalState* _vfsproj;

    ///Jet shape parameters
    double _rmin;
    double _rmax;
    double _interval;
    double _r1minPsi;
    schemelist _distscheme;
    int _nbins;

    ///The jet axes of the jet algorithm projection
    /// @todo Are the pointers really needed?
    vector<LorentzVector>* _jetaxes;
    ///Output lists containg jet shapes
    /// @todo Are the pointers really needed?
    vector<vector<double> >*  _diffjetshapes;
    /// @todo Are the pointers really needed?
    vector<vector<double> >*  _intjetshapes;
    /// @todo Are the pointers really needed?
    vector<double>* _oneminPsishape;
    
  };
  
}

#endif
