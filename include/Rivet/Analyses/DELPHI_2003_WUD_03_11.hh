// -*- C++ -*-
#ifndef RIVET_DELPHI_2003_WUD_03_11_HH
#define RIVET_DELPHI_2003_WUD_03_11_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /**
   * @brief DELPHI 4-jet angular distributions
   * @author Hendrik Hoeth
   *
   * This is Hendrik Hoeth's Diploma thesis, measuring the 4-jet angular
   * distributions at LEP-1.
   *
   *
   * @par Run conditions
   *
   * @arg LEP1 beam energy: \f$ \sqrt{s} = \$f 91.2 GeV
   * @arg Run with generic QCD events.
   * @arg No \f$ p_\perp^\text{min} \f$ cutoff is required
   */
  class DELPHI_2003_WUD_03_11 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    DELPHI_2003_WUD_03_11();

    /// Factory method.
    static Analysis* create() { 
      return new DELPHI_2003_WUD_03_11(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event& event);
    virtual void finalize();

    double calc_BZ(std::vector<fastjet::PseudoJet> jets);
    double calc_KSW(std::vector<fastjet::PseudoJet> jets);
    double calc_NR(std::vector<fastjet::PseudoJet> jets);
    double calc_ALPHA34(std::vector<fastjet::PseudoJet> jets);
    //@}


  private:

    int _numdurjets;
    int _numjadejets;

    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_histDurhamBZ;
    AIDA::IHistogram1D *_histDurhamKSW;
    AIDA::IHistogram1D *_histDurhamNR;
    AIDA::IHistogram1D *_histDurhamALPHA34;
    AIDA::IHistogram1D *_histJadeBZ;
    AIDA::IHistogram1D *_histJadeKSW;
    AIDA::IHistogram1D *_histJadeNR;
    AIDA::IHistogram1D *_histJadeALPHA34;
    //@}

  };

}

#endif
