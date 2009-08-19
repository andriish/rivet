// -*- C++ -*-
#ifndef RIVET_MC_JetAnalysis_HH
#define RIVET_MC_JetAnalysis_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  class MC_JetAnalysis : public Analysis {

  public:

    /// Default constructor.
    MC_JetAnalysis(const std::string& name, const double& sqrts,
                   const size_t& njet, const std::string& jetpro_name);


    /// @name Analysis methods
    //@{ 
    virtual void init();
    virtual void analyze(const Event& event);
    virtual void finalize();
    //@}
  
  protected:
    
    /// The energy scale and number of jets for which histograms are to be
    /// initialised
    double m_sqrts;
    size_t m_njet;
    
    /// The name of the jet projection to be used for this analysis
    /// (this projection has to be registered by the derived analysis!)
    const std::string m_jetpro_name;

    /// @name Histograms
    //@{
    std::vector<AIDA::IHistogram1D *> _h_log10_d;
    std::vector<AIDA::IDataPointSet *> _h_log10_R;
    std::vector<AIDA::IHistogram1D *> _h_pT_jet;
    std::vector<AIDA::IHistogram1D *> _h_eta_jet;
    std::map<std::pair<size_t, size_t>, AIDA::IHistogram1D*> _h_deta_jets;
    std::map<std::pair<size_t, size_t>, AIDA::IHistogram1D*> _h_dR_jets;
    AIDA::IHistogram1D * _h_jet_multi_exclusive;
    AIDA::IHistogram1D * _h_jet_multi_inclusive;
    AIDA::IDataPointSet * _h_jet_multi_ratio;
    //@}

  };

}

#endif
