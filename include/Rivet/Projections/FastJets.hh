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
#ifdef HAVE_JADE
#include "fastjet/JadePlugin.hh"
#endif


namespace Rivet {
  
  /// Typedef for a collection of PseudoJets.
  typedef vector<fastjet::PseudoJet> PseudoJets;


  /// Project out jets found using fastJet package.
  class FastJets : public JetAlg {

  public:
    /// Wrapper enum for selected Fastjet jet algorithms.
    #ifdef HAVE_JADE
    enum JetAlg { KT, CAM, SISCONE, CDFJETCLU, CDFMIDPOINT, JADE, DURHAM };
    #else
    enum JetAlg { KT, CAM, SISCONE, CDFJETCLU, CDFMIDPOINT };
    #endif


    /// @name Constructors and destructors.
    //@{
    /// Default constructor uses Kt algorithm and E-scheme recombination. 
    /// The FinalState projection must live throughout the run. 
    FastJets(const FinalState& fsp) 
      : _plugin(0)
    { 
      setName("FastJets");
      addProjection(fsp, "FS");
      const double RPARAM = 1.0;
      _jdef = fastjet::JetDefinition(fastjet::kt_algorithm, RPARAM, fastjet::E_scheme);
    }


    /// "Wrapped" argument constructor using Rivet enums for most common
    /// jet alg choices (including some plugins). For the built-in algs,
    /// E-scheme recombination is used. For full control of
    /// FastJet built-in jet algs, use the native arg constructor.
    /// The FinalState projection must live throughout the run. 
    FastJets(const FinalState& fsp, JetAlg alg, double rparameter)
      : _plugin(0)
    {
      setName("FastJets");
      addProjection(fsp, "FS");
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
      #ifdef HAVE_JADE
      } else if (alg == JADE) {
        _plugin = new fastjet::JadePlugin("jade");
        _jdef = fastjet::JetDefinition(_plugin);
      } else if (alg == DURHAM) {
        _plugin = new fastjet::JadePlugin("durham");
        _jdef = fastjet::JetDefinition(_plugin);
      #endif
      }
    }
    
    
    /// Native argument constructor, using FastJet alg/scheme enums.
    FastJets(const FinalState& fsp, fastjet::JetAlgorithm type,
             fastjet::RecombinationScheme recom, double rparameter)
      : _plugin(0)
    {
      setName("FastJets");
      addProjection(fsp, "FS");
      _jdef = fastjet::JetDefinition(type, rparameter, recom);
    }

  
    /// Destructor.
    ~FastJets() { 
      // Delete plugin
      if (_plugin) delete _plugin;
    }
    //@}


  public:
    /// Return the number of jets above the pt cut 
    size_t getNumJets(double ptmin = 0.) const {
      return _cseq.inclusive_jets(ptmin).size();
    }

    /// Get the pseudo jets (unordered).
    Jets getJets() const;

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
