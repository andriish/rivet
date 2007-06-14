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
    
  public:
    
    /// @name Standard constructors and destructors.
    //@{
    /// Default constructor. Must specify a FinalState projection which is
    ///  assumed to live throughout the run.    
    inline D0ILConeJets(FinalState& fsp)
      /// @todo Why is _fsproj a pointer?
      : _cone_radius(0.7), _min_jet_Et(0.0), _split_ratio(0.5), _fsproj(&fsp)
    { 
      addProjection(fsp);
      
      // The parameters are supposed to be set as used by D0 in RunII
      /// @todo Set these in the argument list and make the members const
      _far_ratio = 0.5;
      _et_min_ratio = 0.5;
      _kill_duplicate = true;
      _duplicate_dR = 0.005; 
      _duplicate_dPT = 0.01; 
      _search_factor = 1.0; 
      _pT_min_leading_protojet = 0.0; 
      _pT_min_second_protojet = 0.0;
      _merge_max = 10000;
      _pT_min_nomerge = 0.0;
      
      /// @todo Why not use the stack?
      _algo = new ILConeAlgorithm<HepEntity>(_cone_radius, _min_jet_Et, _split_ratio,
                                             _far_ratio, _et_min_ratio, _kill_duplicate, _duplicate_dR, 
                                             _duplicate_dPT, _search_factor, _pT_min_leading_protojet, 
                                             _pT_min_second_protojet, _merge_max, _pT_min_nomerge);
      /// @todo Why not use the stack?
      jets = new list<HepEntity>;
    }

    
    /// Argument constructor.
    /// Added so that same projection can be ran but with different parameters.
    /// Must specify a FinalState projection which is
    /// assumed to live throughout the run. 
    inline D0ILConeJets(FinalState& fsp, float R, float Etmin, float split) 
      : _cone_radius(R), _min_jet_Et(Etmin), _split_ratio(split), _fsproj(&fsp) 
    {
      /// @todo Why not use the stack?
      _algo = new ILConeAlgorithm<HepEntity>(_cone_radius, _min_jet_Et, _split_ratio,
                                             _far_ratio, _et_min_ratio, _kill_duplicate, _duplicate_dR, 
                                             _duplicate_dPT, _search_factor, _pT_min_leading_protojet, 
                                             _pT_min_second_protojet, _merge_max, _pT_min_nomerge);
      /// @todo Most of the variables in the previous line are undefined!

      /// @todo Why not use the stack?
      jets = new list<HepEntity>;
    }

    
    /// Destructor.
    virtual ~D0ILConeJets() { 
      /// @todo Memory leaks are very likely...
      delete _particlelist; 
      delete _algo;
      delete jets;
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
    
    inline int getNJets() const { return jets->size(); }
  

  public:
    /// @todo Why public? And why a pointer to a list? Name should get 
    /// an underscore prefix and there should be relevant reference-based 
    /// accessor methods.
    list<HepEntity>* jets;


  private:
    /// The assignment operator is private and must never be called.
    /// In fact, it shouldn't even be implemented.
    D0ILConeJets& operator=(const D0ILConeJets&);
  
    /// @todo Why a pointer to the list?
    list<const HepEntity*>* _particlelist;
    
    // Initialize D0RunII cone algorithm (what does this mean?)
    float _cone_radius;
    float _min_jet_Et;
    float _split_ratio;

    /// @todo Documentation! 
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

    /// The FinalState projection used by this projection.
    /// @todo Why a pointer?
    FinalState* _fsproj;

    /// @todo Why a pointer?
    ILConeAlgorithm<HepEntity>* _algo;

  };

}

#endif
