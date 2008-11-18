// -*- C++ -*-
#ifndef RIVET_D0_2008_S7863608_HH
#define RIVET_D0_2008_S7863608_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/IsolationTools.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {


  /// @brief Measurement differntial Z/gamma* + jet +X cross sections
  /// @author Gavin Hesketh
  class D0_2008_S7863608:public Analysis {

    typedef IsolationProjection<D0ILConeJets, FinalState> D0JetFromParticleIso;
    typedef MultiplicityInConeEstimator< D0ILConeJets::entity_type, FinalState::collection_type > D0JetIsoEstimator;

  public:

    /// Default constructor.
     D0_2008_S7863608();

    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7863608();
    }
    /// @name Publication metadata
    //@{
    /// Get a description of the analysis. 
    string getSpiresId() const {
      return "7863608";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "Measurement differential Z/gamma* + jet + X cross sections";
    }
    
    /// Experiment which performed and published this analysis. 
    string getExpt() const {
      return "D0";
    }
    /// When published (preprint year according to SPIRES). 
    string getYear() const {
      return "2008";
    }
    /// Publication references.
    vector<string> getReferences() const {
      vector<string> ret;
      ret.push_back("hep-ex/08081296");
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

    /// Hide the assignment operator
    D0_2008_S7863608 & operator=(const D0_2008_S7863608 & x);

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _h_jet_pT_cross_section;
    AIDA::IHistogram1D * _h_jet_y_cross_section;
    AIDA::IHistogram1D * _h_Z_pT_cross_section;
    AIDA::IHistogram1D * _h_Z_y_cross_section;
    AIDA::IHistogram1D * _h_total_cross_section;

    //@}

    double _events;

  };

}

#endif
