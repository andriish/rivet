// -*- C++ -*-
#ifndef RIVET_CDF_2005_S6080774_HH
#define RIVET_CDF_2005_S6080774_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  class CDF_2005_S6080774 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_2005_S6080774();

    /// Factory method
    static Analysis* create() {
      return new CDF_2005_S6080774();
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

    /// @name Histograms
    //@{

    std::vector<AIDA::IHistogram1D*> _h_m_PP;
    std::vector<AIDA::IHistogram1D*> _h_pT_PP;
    std::vector<AIDA::IHistogram1D*> _h_dphi_PP;
    //@}

  };


}

#endif

