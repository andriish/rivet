// -*- C++ -*-
#ifndef RIVET_D0_2008_S6879055_HH
#define RIVET_D0_2008_S6879055_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/IsolationTools.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {


  /// @brief Measurement of the ratio sigma(Z/gamma* + n jets)/sigma(Z/gamma*)
  class D0_2008_S6879055 : public Analysis {

    typedef IsolationProjection<D0ILConeJets, FinalState> D0JetFromParticleIso;
    typedef MultiplicityInConeEstimator< D0ILConeJets::entity_type, FinalState::collection_type > D0JetIsoEstimator;

  public:

    /// Default constructor.
     D0_2008_S6879055();

    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S6879055();
    }

    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "6879055";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Measurement of the ratio sigma(Z/gamma* + n jets)/sigma(Z/gamma*)";
    }
    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "D0";
    }
    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2008";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("hep-ex/0608052");
      return ret;
    }

    //@}
    /// @name Analysis methods
    //@{ 
    void init();
    void analyze(const Event & event);
    void finalize();
    //@}

  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _crossSectionRatio;
    AIDA::IHistogram1D * _pTjet1;
    AIDA::IHistogram1D * _pTjet2;
    AIDA::IHistogram1D * _pTjet3;
    //@}

  };

}

#endif
