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

    /// @name Constructors and destructors.
    //@{
    /// Default constructor uses Kt algorithm and E-scheme recombination. 
    /// The FinalState projection must live throughout the run. 
    FastJets(FinalState* fsp) : _plugin(0)
    { 
      defaultConstructor(fsp);
    }

    FastJets(FinalState& fsp) : _plugin(0)
    { 
      defaultConstructor(&fsp);
    }
    
  private:
    void defaultConstructor(FinalState* fsp){
      _fsproj = fsp;
      addProjection(*_fsproj);
      const double RPARAM = 1.0;
      _jdef = fastjet::JetDefinition(fastjet::kt_algorithm, RPARAM, fastjet::E_scheme);
      return;
    }
    
  public:
    
    /// Wrapper enum for selected Fastjet jet algorithms.
    enum JetAlg { KT, CAM, SISCONE, CDFJETCLU, CDFMIDPOINT };


    /// "Wrapped" argument constructor using Rivet enums for most common
    /// jet alg choices (including some plugins). For the built-in algs,
    /// E-scheme recombination is used. For full control of
    /// FastJet built-in jet algs, use the native arg constructor.
    /// The FinalState projection must live throughout the run. 
    FastJets(FinalState* fsp, JetAlg alg, double rparameter)
    : _plugin(0)
    {
      wrappedConstructor(fsp, alg, rparameter);
    }

    FastJets(FinalState& fsp, JetAlg alg, double rparameter)
    //: _fsproj(fsp), _plugin(0)
    //FastJets(FinalState* fsp, JetAlg alg, double rparameter)
    : _plugin(0)
    {
      wrappedConstructor(&fsp, alg, rparameter);
    }
    
  private:
    void wrappedConstructor(FinalState* fsp, JetAlg alg, double rparameter){
      _fsproj = fsp;
      addProjection(*_fsproj);
      if (alg == KT) {
        _jdef = fastjet::JetDefinition(fastjet::kt_algorithm, rparameter, fastjet::E_scheme);
      } else if (alg == CAM) {
        _jdef = fastjet::JetDefinition(fastjet::cambridge_algorithm, rparameter, fastjet::E_scheme);
      } else if (alg == SISCONE) {
        const double overlapthreshold = 0.5;
        _plugin = new fastjet::SISConePlugin(rparameter, overlapthreshold);
        _jdef = fastjet::JetDefinition(_plugin);
      } else if (alg == CDFJETCLU) {
        _plugin = new fastjet::CDFJetCluPlugin(rparameter);
        _jdef = fastjet::JetDefinition(_plugin);
      } else if (alg == CDFMIDPOINT) {
        _plugin = new fastjet::CDFMidPointPlugin(rparameter);
        _jdef = fastjet::JetDefinition(_plugin);
      }
      return;
    }
    
  public:
    
    /// Native argument constructor, using FastJet alg/scheme enums.
    /// The FinalState projection must live throughout the run. 
    FastJets(FinalState* fsp, fastjet::JetAlgorithm type,
             fastjet::RecombinationScheme recom, double rparameter)
    : _plugin(0)
    {
      nativeConstructor(fsp, type, recom, rparameter);
    }

    FastJets(FinalState& fsp, fastjet::JetAlgorithm type,
             fastjet::RecombinationScheme recom, double rparameter)
    : _plugin(0)
    {
      nativeConstructor(&fsp, type, recom, rparameter);
    }
    
  private:
    
    void nativeConstructor(FinalState* fsp, fastjet::JetAlgorithm type,
                           fastjet::RecombinationScheme recom, double rparameter){
      _fsproj = fsp;
      addProjection(*_fsproj);
      _jdef = fastjet::JetDefinition(type, rparameter, recom);
      return;
    }

  public:
    
    
    /// Destructor.
    ~FastJets() { 
      // Delete plugin
      if (_plugin) delete _plugin;
    }
    //@}


  public:
    /// Return the name of the projection
    string getName() const {
      return "FastJets";
    }
    
    /// Return the number of jets above the pt cut 
    size_t getNumJets(double ptmin = 0.) const {
      return _cseq.inclusive_jets(ptmin).size();
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
    PseudoJets getPseudoJets(double ptmin = 0.) const {
      return _cseq.inclusive_jets(ptmin);
    }

    /// Get the jets, ordered by \f$ p_T \f$.
    PseudoJets getPseudoJetsPt(double ptmin = 0.) const {
      return sorted_by_pt(_cseq.inclusive_jets(ptmin));
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

    /// The FinalState projection used by this projection.
    FinalState* _fsproj;
    
    /// Jet definition
    fastjet::ClusterSequence _cseq;
    fastjet::JetDefinition _jdef;

    /// FastJet external plugin
    fastjet::JetDefinition::Plugin* _plugin; 

    /// Map of vectors of y scales. This is mutable so we can use caching/lazy evaluation.
    mutable map<int, vector<double> > _yscales;
  
  };

}

#endif
