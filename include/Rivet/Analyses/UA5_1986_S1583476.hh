// -*- C++ -*-
#ifndef RIVET_UA5_1986_S1583476_HH
#define RIVET_UA5_1986_S1583476_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  class UA5_1986_S1583476 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    string name() const {
        return "UA5_1986_S1583476";
    }
    /// Constructor
    UA5_1986_S1583476();

    /// Factory method
    static Analysis* create() {
      return new UA5_1986_S1583476();
    }
    //@}


  public:

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    // Data members like post-cuts event weight counters go here

  private:

    /// @name Histograms
    //@{
    // Histos of Figure 1 (HepData Table 1)
    AIDA::IHistogram1D *_hist_eta_nsd_200;
    AIDA::IHistogram1D *_hist_eta_inelastic_200;
    AIDA::IHistogram1D *_hist_eta_nsd_900;
    AIDA::IHistogram1D *_hist_eta_inelastic_900;

    // Histos of Figure 3a (HepData Table 2)
    AIDA::IHistogram1D *_hist_eta_nsd_n_2_200;
    AIDA::IHistogram1D *_hist_eta_nsd_n_12_200;
    AIDA::IHistogram1D *_hist_eta_nsd_n_22_200;
    AIDA::IHistogram1D *_hist_eta_nsd_n_32_200;
    AIDA::IHistogram1D *_hist_eta_nsd_n_42_200;
    AIDA::IHistogram1D *_hist_eta_nsd_n_52_200;

    // Histos of Figure 3b (HepData Table 3)
    AIDA::IHistogram1D *_hist_eta_nsd_n_2_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_12_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_22_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_32_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_42_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_52_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_62_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_72_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_82_900;
    //@}

  };


}

#endif
