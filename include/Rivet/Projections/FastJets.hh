// -*- C++ -*-
#ifndef RIVET_FastJets_HH
#define RIVET_FastJets_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "fastjet/ClusterSequence.hh"

namespace Rivet {
  
  /// Project out kT jets found using fastJet package.
  class FastJets : public Projection {
    
  public:
    
    /// @name Standard constructors and destructors.
    //@{
    /// Default constructor. Must specify a FinalState projection which is
    //  assumed to live throughout the run.
    inline FastJets(FinalState& fsp)
      : _cseq(0), _type(fastjet::kt_algorithm), _recom(fastjet::E_scheme), 
        _rparameter(1.0), _fsproj(&fsp) 
    { 
      addProjection(fsp);
      _jdef = fastjet::JetDefinition(_type,_rparameter,_recom); 
    }

    // @todo implement configurable constructor for fastJet
    /// Argument constructor. Allows the to be run with different parameters.
    /// Must specify a FinalState projection which is assumed to live throughout the run. 
    //    inline FastJets(FinalState& fsp, int type, int angle, int recom, double rparameter)
    //  : _pktev(0), _type(type), _angle(angle), _recom(recom),
    //    _rparameter(rparameter), _fsproj(&fsp)
    //{ 
    //  addProjection(fsp);
    //}
    
    /// Destructor.
    inline virtual ~FastJets() { 
      if (_cseq) delete _cseq; 
    }
    //@}

  public:
    /// Return the name of the projection
    inline string getName() const {
      return "FastJets";
    }
    
  protected:   

    /// Perform the projection on the Event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;  

  public:
    
    // @todo for more efficiency, pass the pt cut.
    inline int getNJets() const {
      return _cseq->inclusive_jets().size();
    }
    inline vector<fastjet::PseudoJet> getJetsPt() const {
      return sorted_by_pt(_cseq->inclusive_jets());
    }

    /// return the cluster sequence.
    inline fastjet::ClusterSequence getCSeq() const {
      return *_cseq;
    }

    /// return the jet definitions
    inline fastjet::JetDefinition getJetDef() const {
      return _jdef;
    }

    /// Get the subjet splitting variables for the given jet.
    vector<double> getYSubJet(const fastjet::PseudoJet& jet) const; 

  private:
    

    /// Internal KtEvent, rebuilt every time an event is projected, but not otherwise.
    fastjet::ClusterSequence* _cseq;

    fastjet::JetFinder _type;
    fastjet::RecombinationScheme _recom;
    double _rparameter;  

    /// The FinalState projection used by this projection.
    FinalState* _fsproj;

    /// Jet definition
    fastjet::JetDefinition  _jdef;

    fastjet::Strategy _strat;

    /// Map of vectors of y scales. This is mutable so we can use caching/lazy evaluation.
    mutable map<int, vector<double> > _yscales;
  
  };

}

#endif
