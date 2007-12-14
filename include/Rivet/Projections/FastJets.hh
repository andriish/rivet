// -*- C++ -*-
#ifndef RIVET_FastJets_HH
#define RIVET_FastJets_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"
#include "fastjet/CDFJetCluPlugin.hh"
#include "fastjet/CDFMidPointPlugin.hh"


namespace Rivet {
  
  /// Typedef for a collection of PseudoJets.
  typedef vector<fastjet::PseudoJet> PseudoJets;


  /// Project out jets found using fastJet package.
  class FastJets : public JetAlg {

  public:

    /// @name Standard constructors and destructors.
    //@{
    /// Default constructor. Must specify a FinalState projection which is
    //  assumed to live throughout the run.
    FastJets(FinalState& fsp)
      : _cseq(), _type(fastjet::kt_algorithm), _recom(fastjet::E_scheme), 
        _rparameter(1.0), _fsproj(fsp) 
    { 
      addProjection(fsp);
      _jdef = fastjet::JetDefinition(_type,_rparameter,_recom); 
    }

    // @todo implement configurable constructor for fastJet
    /// Argument constructor. Allows the to be run with different parameters.
    /// Must specify a FinalState projection which is assumed to live throughout the run. 
    inline FastJets(fastjet::JetFinder type,
                    fastjet::RecombinationScheme recom, double rparameter, 
		    FinalState& fsp)
      : _type(type), _recom(recom),
	_rparameter(rparameter), _fsproj(fsp)
    {
      addProjection(fsp);
      _jdef = fastjet::JetDefinition(_type, _rparameter, _recom);
    }
    


    enum ExtJetType {SISCone, CDFJetClu, CDFMidPoint};
    /// Argument constructor for external plugin algorithms
    inline FastJets(ExtJetType type, double rparameter, FinalState& fsp)
      : _extype(type), _rparameter(rparameter), _fsproj(fsp)
    {
      addProjection(fsp);
      if (_extype == SISCone) {
	double _overlapthreshold = 0.5;
	_plugin = new fastjet::SISConePlugin(_rparameter,_overlapthreshold);
	_jdef = fastjet::JetDefinition(_plugin);
      }
      else if (_extype == CDFJetClu) { //use external plugin jet algorithm
	_plugin = new fastjet::CDFJetCluPlugin(_rparameter);
	_jdef = fastjet::JetDefinition(_plugin);
      }
      else if (_extype == CDFMidPoint) {
	_plugin = new fastjet::CDFMidPointPlugin(_rparameter);
	_jdef = fastjet::JetDefinition(_plugin);
      }
    }
    


    /// Destructor.
    ~FastJets() { }
    //@}


  public:
    /// Return the name of the projection
    string getName() const {
      return "FastJets";
    }
    
    /// @todo For more efficiency, pass the \f$ p_T \f$ cut. (What?)
    size_t getNumJets() const {
      return _cseq.inclusive_jets().size();
    }

    /// Get the pseudo jets (unordered).
    Jets getJets() const {
      Jets rtn;
      const PseudoJets pjets = _cseq.inclusive_jets();
      for (PseudoJets::const_iterator pj = pjets.begin(); pj != pjets.end(); ++pj) {
        Jet j;
        const PseudoJets parts = getClusterSeq().constituents(*pj);
        for (PseudoJets::const_iterator p = parts.begin(); p != parts.end(); ++p) {
          const FourMomentum particle(p->E(), p->px(), p->py(), p->pz());
          j.addParticle(particle);
        }
        rtn.push_back(j);
      }
      return rtn;
    }


    /// Get the pseudo jets (unordered).
    PseudoJets getPseudoJets() const {
      return _cseq.inclusive_jets();
    }

    /// Get the jets, ordered by \f$ p_T \f$.
    PseudoJets getPseudoJetsPt() const {
      return sorted_by_pt(_cseq.inclusive_jets());
    }

    /// Return the cluster sequence (FastJet-specific).
    const fastjet::ClusterSequence& getClusterSeq() const {
      return _cseq;
    }

    /// Return the jet definition (FastJet-specific).
    const fastjet::JetDefinition& getJetDef() const {
      return _jdef;
    }

    /// Get the subjet splitting variables for the given jet.
    vector<double> getYSubJet(const fastjet::PseudoJet& jet) const;


  protected:   

    /// Perform the projection on the Event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;  



  private:
    
    ///FastJet external jetalgo parameters
    ExtJetType _extype;
    fastjet::JetDefinition::Plugin  * _plugin; 

    /// FastJet parameters
    fastjet::ClusterSequence _cseq;
    fastjet::JetFinder _type;
    fastjet::RecombinationScheme _recom;
    double _rparameter;  


    /// The FinalState projection used by this projection.
    FinalState _fsproj;

    /// Jet definition
    fastjet::JetDefinition  _jdef;
    fastjet::Strategy _strat;

    /// Map of vectors of y scales. This is mutable so we can use caching/lazy evaluation.
    mutable map<int, vector<double> > _yscales;
  
  };

}

#endif
