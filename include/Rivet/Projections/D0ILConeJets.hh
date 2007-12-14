// -*- C++ -*-
#ifndef RIVET_D0ILConeJets_HH
#define RIVET_D0ILConeJets_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Projection.hh"
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
  class D0ILConeJets : public Projection {

  public:
    
    /// @name Constructors and destructors.
    //@{

    /// @brief Default constructor. 
    ///
    /// The algorithm parameters are supposed to be set as used by D0 in RunII -
    /// this constructor will initialise the correct parameter values.
    /// Must specify a FinalState projection which is assumed to live throughout the run.    
    D0ILConeJets(FinalState& fsp)
      : _fsproj(fsp), _cone_radius(0.7), _min_jet_Et(0.0), 
        _split_ratio(0.5), _far_ratio(0.5), 
        _et_min_ratio(0.5), _kill_duplicate(true), 
        _duplicate_dR(0.005), _duplicate_dPT(0.01), 
        _search_factor(1.0), 
        _pT_min_leading_protojet(0.0), 
        _pT_min_second_protojet(0.0), 
        _merge_max(1000), _pT_min_nomerge(0.0), 
        _algo(_cone_radius, _min_jet_Et, _split_ratio,
              _far_ratio, _et_min_ratio, _kill_duplicate, _duplicate_dR, 
              _duplicate_dPT, _search_factor, _pT_min_leading_protojet, 
              _pT_min_second_protojet, _merge_max, _pT_min_nomerge)
    { 
      addProjection(fsp);
    }

        
    /// @brief Alternative constructor.
    ///
    /// Added so that same projection can be ran but with different parameters.
    /// Must specify a FinalState projection which is assumed to live throughout the run. 
    D0ILConeJets(FinalState& fsp, float r, float etMin, float split,
			float farRatio, float etMinRatio, bool killDuplicate,
			float duplicateDR, float duplicateDPT, float searchFactor,
			float pTMinLeadingProtojet, float pTMinSecondProtojet,
			int mergeMax, float pTMinNomerge)
      : _fsproj(fsp), _cone_radius(r), _min_jet_Et(etMin), 
        _split_ratio(split), _far_ratio(farRatio), 
        _et_min_ratio(etMinRatio), _kill_duplicate(killDuplicate),
        _duplicate_dR(duplicateDR), _duplicate_dPT(duplicateDPT), 
        _search_factor(searchFactor),
        _pT_min_leading_protojet(pTMinLeadingProtojet), 
        _pT_min_second_protojet(pTMinSecondProtojet), 
        _merge_max(mergeMax), _pT_min_nomerge(pTMinNomerge),
        _algo(_cone_radius, _min_jet_Et, _split_ratio,
              _far_ratio, _et_min_ratio, _kill_duplicate, _duplicate_dR, 
              _duplicate_dPT, _search_factor, _pT_min_leading_protojet, 
              _pT_min_second_protojet, _merge_max, _pT_min_nomerge)
    {  
      addProjection(fsp);
    }
    

    /// Destructor.
    virtual ~D0ILConeJets() {};

    //@}

  public:
    /// Return the name of the projection
    string getName() const {
      return "D0ILConeJets";
    }

  protected:   

    /// Perform the projection on the Event.
    void project(const Event& e);

    /// This function defines a unique ordering between different 
    /// Projection objects of the same class. Should only be called from 
    /// operator<(const Projection&).
    int compare(const Projection& p) const;  

  public:
    
    /// Get the number of jets.
    int getNJets() const { return _jets.size(); }
  
    /// Get a reference to the jets collection.
    list<HepEntity>& getJets() { return _jets; }
    /// Get a reference to the jets collection (const version).
    const list<HepEntity>& getJets() const { return _jets; }

    /// Get a reference to the lorentzvecjets collection.
    const list<FourMomentum>& getLorentzJets() const {
      return _lorentzvecjets; 
    }


    /// Clear the jets list.
    D0ILConeJets& clearJets() { 
      _jets.clear(); 
      return *this; 
    }



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

    /// The FinalState projection used by this projection.
    FinalState _fsproj;


    /// @name Cone algorithm parameters
    //@{
    /// Cone radius
    const float _cone_radius;
    /// Jet \f$ E_T \f$ threshold
    const float _min_jet_Et;
    /// Split/merge fraction
    const float _split_ratio;
    //@}


    /// @name More algorithm parameters
    //@{
    /// The original author Laurent Duflot might be able to explain these.
    const float _far_ratio;
    const float _et_min_ratio;
    const bool _kill_duplicate;
    const float _duplicate_dR; 
    const float _duplicate_dPT; 
    const float _search_factor; 
    const float _pT_min_leading_protojet; 
    const float _pT_min_second_protojet;
    const int _merge_max; 
    const float _pT_min_nomerge;
    //@}

    /// The jet algorithm function itself
    ILConeAlgorithm<HepEntity> _algo;

  };

}

#endif
