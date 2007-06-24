// -*- C++ -*-
#ifndef RIVET_D0ILConeJets_HH
#define RIVET_D0ILConeJets_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

#include "Rivet/Tools/D0RunIIcone/HepEntity.h"
#include "Rivet/Tools/D0RunIIcone/energycluster/ILConeAlgorithm.hpp"


namespace Rivet {
  
  /// Find jets according to the D0 "improved legacy" cone algorithm.
  class D0ILConeJets : public Projection {

    ///@name description
    //@{ 
    /// The jet algorithm is described in:
    /// D0 Note Number: 003750, Date: 4/20/00,
    /// Title: Run II Jet Physics,
    /// Proceedings of the jet physics group of the Run II QCD Workshop
    /// Author(s): Gerald C. Blazey, Jay R. Dittmann, Stephen D. Ellis, V. Daniel Elvira, K. Frame, S. Grinstein, Robert Hirosky, R.Peigaia, H. Schellman, R. Snihur,V. Sorin, Dieter Zeppenfeld
    ///The actual implementation differs in details like: mid points are only considered between
    /// 4-vectors above threshold and the midpoint is determined pT weighted.
    
  public:
    
    /// @name Standard constructors and destructors.
    //@{
    /// Default constructor. Must specify a FinalState projection which is
    ///  assumed to live throughout the run.    
    inline D0ILConeJets(FinalState& fsp)
      /// @todo Why is _fsproj a pointer?
      : _fsproj(&fsp), _cone_radius(0.7), _min_jet_Et(0.0), _split_ratio(0.5),
      _far_ratio(0.5), _et_min_ratio(0.5), _kill_duplicate(true), _duplicate_dR(0.005), 
      _duplicate_dPT(0.01), _search_factor(1.0), _pT_min_leading_protojet(0.0), 
      _pT_min_second_protojet(0.0), _merge_max(1000), _pT_min_nomerge(0.0)
	, _algo(_cone_radius, _min_jet_Et, _split_ratio,
		_far_ratio, _et_min_ratio, _kill_duplicate, _duplicate_dR, 
		_duplicate_dPT, _search_factor, _pT_min_leading_protojet, 
		_pT_min_second_protojet, _merge_max, _pT_min_nomerge)

    { 
      addProjection(fsp);
    }


        
    /// Argument constructor.
    /// Added so that same projection can be ran but with different parameters.
    /// Must specify a FinalState projection which is
    /// assumed to live throughout the run. 
    inline D0ILConeJets(FinalState& fsp, float R, float Etmin, float split,
			float Far_Ratio, float Et_Min_Ratio, bool Kill_Duplicate,
			float Duplicate_DR, float Duplicate_DPT, float Search_Factor,
			float PT_Min_Leading_Protojet, float PT_Min_Second_Protojet,
			int Merge_Max, float PT_Min_Nomerge)
      : _fsproj(&fsp), _cone_radius(R), _min_jet_Et(Etmin), _split_ratio(split),
	_far_ratio(Far_Ratio), _et_min_ratio(Et_Min_Ratio), _kill_duplicate(Kill_Duplicate),
	_duplicate_dR(Duplicate_DR), _duplicate_dPT(Duplicate_DPT), _search_factor(Search_Factor),
	_pT_min_leading_protojet(PT_Min_Leading_Protojet), 
	_pT_min_second_protojet(PT_Min_Second_Protojet), _merge_max(Merge_Max),
	_pT_min_nomerge(PT_Min_Nomerge)
	, _algo(_cone_radius, _min_jet_Et, _split_ratio,
		_far_ratio, _et_min_ratio, _kill_duplicate, _duplicate_dR, 
		_duplicate_dPT, _search_factor, _pT_min_leading_protojet, 
		_pT_min_second_protojet, _merge_max, _pT_min_nomerge)
 
    {  }
    

    
    /// Destructor.
    virtual ~D0ILConeJets() { 

    };
    //@}

  public:
    /// Return the name of the projection
    inline string getName() const {
      return "D0ILConeJets";
    }

  protected:   

    /// Perform the projection on the Event: only to be called by 
    /// Event::applyProjection(Projection &).
    void project(const Event& e);

    /// This function defines a unique ordering between different 
    /// Projection objects of the same class. Should only be called from 
    /// operator<(const Projection&).
    int compare(const Projection& p) const;  

  public:
    
    inline int getNJets() const { return _jets.size(); }
  
    //Here is a problem: 
    //If I stick with the upper version alone I get into trouble in D0IIConeJets.cc
    //If I stick with the lower version alone I get into trouble in ILConeAlgorithm.hpp
    //I tried to dig into both possibilities but get still more trouble
    //Notice that this only happens since the public list jets has been made private: _jets
    inline list<HepEntity>& getJets() { return _jets; }
    inline const list<HepEntity>& getJets() const { return _jets; }

    inline void clearJets() { _jets.clear(); return; }

  private:
    list<HepEntity> _jets;


//   private:
//     /// The assignment operator is private and must never be called.
//     /// In fact, it shouldn't even be implemented.
//     D0ILConeJets& operator=(const D0ILConeJets&);
  
    list<HepEntity> _particlelist;
    list<const HepEntity*> _particlepointerlist;

    /// The FinalState projection used by this projection.
    /// @todo Why a pointer?
    /// That's how I learned it from the KtJet projection
    FinalState* _fsproj;



    // Initialize D0RunII cone algorithm (what does this mean?)
    // This means that when the constructer gets instantiated the cone radius has to be specified
    // D0 JCCA: _cone_radius=0.7, D0 JCCB: _cone_radius=0.5
    const float _cone_radius;
    const float _min_jet_Et;
    const float _split_ratio;

    /// @todo Documentation! 
    /// The original author Laurent Duflot might be able to explain those.
    /// The parameters are supposed to be set as used by D0 in RunII,
    /// correct values are preset in the default constructor
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

    ///The jet algorithm function itself
     ILConeAlgorithm<HepEntity> _algo;

  };

}

#endif
