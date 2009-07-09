// -*- C++ -*-
#ifndef RIVET_CDF_2008_S7540469_HH
#define RIVET_CDF_2008_S7540469_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Measurement differential Z/gamma* + jet + X cross sections
  /// @author Frank Siegert
  class CDF_2008_S7540469 : public Analysis {

  public:

    /// Default constructor.
    CDF_2008_S7540469();


    /// Factory method 
    static Analysis* create() {
      return new CDF_2008_S7540469();
    }


    /// @name Analysis methods
    //@{ 
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _h_jet_multiplicity;
    AIDA::IHistogram1D * _h_jet_pT_cross_section_incl_1jet;
    AIDA::IHistogram1D * _h_jet_pT_cross_section_incl_2jet;
    //@}

  };

}

#endif
