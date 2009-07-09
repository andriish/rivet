// -*- C++ -*-
#ifndef RIVET_D0_2004_S5992206_HH
#define RIVET_D0_2004_S5992206_HH

#include "Rivet/Analysis.hh"

namespace Rivet {  


  /* @brief D0 Run II jet analysis
   * @author Lars Sonnenschein
   * 
   * Measurement of angular correlations in di-jet events.
   * 
   * 
   * @par Run conditions
   * 
   * @arg \f$ \sqrt{s} = \f$ 1960 GeV
   * @arg Run with generic QCD events.
   * @arg Several \f$ p_\perp^\text{min} \f$ cutoffs are probably required to fill the histograms:
   *   @arg \f$ p_\perp^\text{min} = \f$ 50, 75, 100, 150 GeV for the four pT ranges respecively
   * 
   */ 
  class D0_2004_S5992206 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor.
    D0_2004_S5992206();

    /// Factory method
    static Analysis* create() { 
      return new D0_2004_S5992206(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histJetAzimuth_pTmax75_100;
    AIDA::IHistogram1D* _histJetAzimuth_pTmax100_130;
    AIDA::IHistogram1D* _histJetAzimuth_pTmax130_180;
    AIDA::IHistogram1D* _histJetAzimuth_pTmax180_;
    //@}

  };

}

#endif
