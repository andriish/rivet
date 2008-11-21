// -*- C++ -*-
#ifndef RIVET_D0ILConeJets_HH
#define RIVET_D0ILConeJets_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/D0RunIIcone/HepEntity.h"
#include "Rivet/Tools/D0RunIIcone/energycluster/ILConeAlgorithm.hpp"

namespace Rivet {

  
  /// @brief Find jets according to the D0 "improved legacy" cone algorithm.
  ///
  /// The jet algorithm is described in:
  ///
  /// D0 Note Number: 003750, Date: 4/20/00,
  /// Title: Run II Jet Physics,
  /// Proceedings of the jet physics group of the Run II QCD Workshop
  /// Author(s): Gerald C. Blazey, Jay R. Dittmann, Stephen D. Ellis, 
  /// V. Daniel Elvira, K. Frame, S. Grinstein, Robert Hirosky, R.Peigaia, 
  /// H. Schellman, R. Snihur,V. Sorin, Dieter Zeppenfeld.
  ///
  /// The actual implementation differs in details: mid points are only 
  /// considered between 4-vectors above threshold and the midpoint is 
  /// determined \f$ p_T \f$-weighted.
  ///
  /// The cone radius has to be specified according to the analysis
  /// D0 JCCA: _cone_radius=0.7, D0 JCCB: _cone_radius=0.5    
  class D0ILConeJets : public JetAlg {

  public:

    typedef HepEntity entity_type;
    typedef list<HepEntity> collection_type;

    
    /// @name Constructors and destructors.
    //@{

    /// @brief Default constructor. 
    ///
    /// The algorithm parameters are supposed to be set as used by D0 in Run II.
    /// this constructor will initialise the correct parameter values.
    D0ILConeJets(const FinalState& fsp, float cone_radius=0.7, float etMin=0.0)
      : _cone_radius(cone_radius), _min_jet_Et(etMin)
    { 
      setName("D0ILConeJets");
      addProjection(fsp, "FS");
      _init_extra_params();
    }

    // Internal helper
    /// @ignore
    void _init_extra_params();

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new D0ILConeJets(*this);
    }
    //@}


  protected:   

    /// Perform the projection on the Event.
    void project(const Event& e);
    
    /// This function defines a unique ordering between different 
    /// Projection objects of the same class. Should only be called from 
    /// operator<(const Projection&).
    int compare(const Projection& p) const;  


  public:

    /// Do the calculation locally (no caching).
    void calc(const ParticleVector& ps);


  public:

    /// @todo Kill off the HepEntity crap.
    
    /// Number of jets.
    size_t size() const { return _jets.size(); }
    
    /// Get a reference to the jets collection.
    /// @deprecated Analyses should not use HepEntities: prefer the JetAlg interface
    //list<HepEntity>& getJets() { return _jets; }

    /// Get a reference to the jets collection (const version).
    /// @deprecated Analyses should not use HepEntities: prefer the JetAlg interface
    //const list<HepEntity>& getJets() const { return _jets; }

    /// Common interface to FinalState and JetAlg
    /// @deprecated Analyses should not use HepEntities: prefer the JetAlg interface
    const list<HepEntity>& entities() const { return _jets; }

    /// Get jets in the standard way offered by the JetAlg interface.
    Jets jets(double ptmin=0.0*GeV) const;

    /// Get a reference to the lorentzvecjets collection.
    const list<FourMomentum>& getLorentzJets() const {
      return _lorentzvecjets; 
    }

    /// Clear the jets list.
    void reset();


  private:

    /// The collection of jets.
    list<HepEntity> _jets;

    /// Collection of jets converted to list of Lorentz vectors
    list<FourMomentum> _lorentzvecjets;

    /// List of the event particles
    list<HepEntity> _particlelist;
    /// List of the event particles (as pointers)
    /// @todo Why both? Can we eliminate the pointers?
    list<const HepEntity*> _particlepointerlist;

    /// @name Cone algorithm parameters
    //@{
    /// Cone radius
    float _cone_radius;
    /// Jet \f$ E_T \f$ threshold
    float _min_jet_Et;
    /// Split/merge fraction
    float _split_ratio;
    //@}

    /// @name More algorithm parameters
    //@{
    /// The original author Laurent Duflot might be able to explain these.
    float _far_ratio;
    float _et_min_ratio;
    bool _kill_duplicate;
    float _duplicate_dR; 
    float _duplicate_dPT; 
    float _search_factor; 
    float _pT_min_leading_protojet; 
    float _pT_min_second_protojet;
    int _merge_max; 
    float _pT_min_nomerge;
    //@}

    /// The jet algorithm function itself
    ILConeAlgorithm<HepEntity> _algo;
  };

}

#endif
