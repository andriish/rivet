// -*- C++ -*-
#ifndef RIVET_PDG_Hadron_Multiplicities_Ratios_HH
#define RIVET_PDG_Hadron_Multiplicities_Ratios_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief Implementation of PDG hadron multiplicities
  /// @author Hendrik Hoeth
  class PDG_Hadron_Multiplicities_Ratios : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor.
    PDG_Hadron_Multiplicities_Ratios() 
    {
      setBeams(ELECTRON, POSITRON); 
      addProjection(Beam(), "Beams");
      addProjection(ChargedFinalState(), "FS");
      addProjection(UnstableFinalState(), "UFS");
      _weightedTotalNumPiPlus10       = 0.;
      _weightedTotalNumPiPlus32       = 0.;
      _weightedTotalNumPiPlus91       = 0.;
      _weightedTotalNumPiPlus165      = 0.;
    }


    /// Factory method.
    static Analysis* create() { 
      return new PDG_Hadron_Multiplicities_Ratios(); 
    }
    //@}


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getName() const {
      return "PDG_Hadron_Multiplicities_Ratios";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "Ratios (w.r.t. piplus/piminus) of hadron multiplicities in hadronic e+e- events, taken from Particle Data Book";
    }
    /// Experiment which performed and published this analysis.
    //string getExpt() const {
    //  return "PDG";
    //}
    ///// When published (preprint year according to SPIRES).
    //string getYear() const {
    //  return "2006";
    //}
    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event& event);
    virtual void finalize();
    //@}


  private:

    /// Hide the assignment operator
    PDG_Hadron_Multiplicities_Ratios& operator=(const PDG_Hadron_Multiplicities_Ratios&);

    double _weightedTotalNumPiPlus10;
    double _weightedTotalNumPiPlus32;
    double _weightedTotalNumPiPlus91;
    double _weightedTotalNumPiPlus165;

    AIDA::IHistogram1D *_hist10MeanMultiPi0;
    AIDA::IHistogram1D *_hist10MeanMultiKPlus;
    AIDA::IHistogram1D *_hist10MeanMultiK0;
    AIDA::IHistogram1D *_hist10MeanMultiEta;
    AIDA::IHistogram1D *_hist10MeanMultiEtaPrime;
    AIDA::IHistogram1D *_hist10MeanMultiDPlus;
    AIDA::IHistogram1D *_hist10MeanMultiD0;
    AIDA::IHistogram1D *_hist10MeanMultiDPlus_s;
    AIDA::IHistogram1D *_hist10MeanMultiF0_980;
    AIDA::IHistogram1D *_hist10MeanMultiRho770_0;
    AIDA::IHistogram1D *_hist10MeanMultiOmega782;
    AIDA::IHistogram1D *_hist10MeanMultiKStar892Plus;
    AIDA::IHistogram1D *_hist10MeanMultiKStar892_0;
    AIDA::IHistogram1D *_hist10MeanMultiPhi1020;
    AIDA::IHistogram1D *_hist10MeanMultiDStar2010Plus;
    AIDA::IHistogram1D *_hist10MeanMultiDStar2007_0;
    AIDA::IHistogram1D *_hist10MeanMultiDStar_s2112Plus;
    AIDA::IHistogram1D *_hist10MeanMultiJPsi1S;
    AIDA::IHistogram1D *_hist10MeanMultiF2_1270;
    AIDA::IHistogram1D *_hist10MeanMultiP;
    AIDA::IHistogram1D *_hist10MeanMultiLambda;
    AIDA::IHistogram1D *_hist10MeanMultiSigma0;
    AIDA::IHistogram1D *_hist10MeanMultiXiMinus;
    AIDA::IHistogram1D *_hist10MeanMultiDelta1232PlusPlus;
    AIDA::IHistogram1D *_hist10MeanMultiSigma1385Minus;
    AIDA::IHistogram1D *_hist10MeanMultiSigma1385Plus;
    AIDA::IHistogram1D *_hist10MeanMultiSigma1385PlusMinus;
    AIDA::IHistogram1D *_hist10MeanMultiXi1530_0;
    AIDA::IHistogram1D *_hist10MeanMultiOmegaMinus;
    AIDA::IHistogram1D *_hist10MeanMultiLambda_c_Plus;
    AIDA::IHistogram1D *_hist10MeanMultiSigma_c_PlusPlus_0;
    AIDA::IHistogram1D *_hist10MeanMultiLambda1520;

    AIDA::IHistogram1D *_hist32MeanMultiPi0;
    AIDA::IHistogram1D *_hist32MeanMultiKPlus;
    AIDA::IHistogram1D *_hist32MeanMultiK0;
    AIDA::IHistogram1D *_hist32MeanMultiEta;
    AIDA::IHistogram1D *_hist32MeanMultiEtaPrime;
    AIDA::IHistogram1D *_hist32MeanMultiDPlus;
    AIDA::IHistogram1D *_hist32MeanMultiD0;
    AIDA::IHistogram1D *_hist32MeanMultiDPlus_s;
    AIDA::IHistogram1D *_hist32MeanMultiF0_980;
    AIDA::IHistogram1D *_hist32MeanMultiRho770_0;
    AIDA::IHistogram1D *_hist32MeanMultiKStar892Plus;
    AIDA::IHistogram1D *_hist32MeanMultiKStar892_0;
    AIDA::IHistogram1D *_hist32MeanMultiPhi1020;
    AIDA::IHistogram1D *_hist32MeanMultiDStar2010Plus;
    AIDA::IHistogram1D *_hist32MeanMultiDStar2007_0;
    AIDA::IHistogram1D *_hist32MeanMultiF2_1270;
    AIDA::IHistogram1D *_hist32MeanMultiK2Star1430Plus;
    AIDA::IHistogram1D *_hist32MeanMultiK2Star1430_0;
    AIDA::IHistogram1D *_hist32MeanMultiP;
    AIDA::IHistogram1D *_hist32MeanMultiLambda;
    AIDA::IHistogram1D *_hist32MeanMultiXiMinus;
    AIDA::IHistogram1D *_hist32MeanMultiSigma1385Minus;
    AIDA::IHistogram1D *_hist32MeanMultiSigma1385Plus;
    AIDA::IHistogram1D *_hist32MeanMultiSigma1385PlusMinus;
    AIDA::IHistogram1D *_hist32MeanMultiOmegaMinus;
    AIDA::IHistogram1D *_hist32MeanMultiLambda_c_Plus;

    AIDA::IHistogram1D *_hist91MeanMultiPi0;
    AIDA::IHistogram1D *_hist91MeanMultiKPlus;
    AIDA::IHistogram1D *_hist91MeanMultiK0;
    AIDA::IHistogram1D *_hist91MeanMultiEta;
    AIDA::IHistogram1D *_hist91MeanMultiEtaPrime;
    AIDA::IHistogram1D *_hist91MeanMultiDPlus;
    AIDA::IHistogram1D *_hist91MeanMultiD0;
    AIDA::IHistogram1D *_hist91MeanMultiDPlus_s;
    AIDA::IHistogram1D *_hist91MeanMultiBPlus_B0_d;
    AIDA::IHistogram1D *_hist91MeanMultiBPlus_u;
    AIDA::IHistogram1D *_hist91MeanMultiB0_s;
    AIDA::IHistogram1D *_hist91MeanMultiF0_980;
    AIDA::IHistogram1D *_hist91MeanMultiA0_980Plus;
    AIDA::IHistogram1D *_hist91MeanMultiRho770_0;
    AIDA::IHistogram1D *_hist91MeanMultiRho770Plus;
    AIDA::IHistogram1D *_hist91MeanMultiOmega782;
    AIDA::IHistogram1D *_hist91MeanMultiKStar892Plus;
    AIDA::IHistogram1D *_hist91MeanMultiKStar892_0;
    AIDA::IHistogram1D *_hist91MeanMultiPhi1020;
    AIDA::IHistogram1D *_hist91MeanMultiDStar2010Plus;
    AIDA::IHistogram1D *_hist91MeanMultiDStar_s2112Plus;
    AIDA::IHistogram1D *_hist91MeanMultiBStar;
    AIDA::IHistogram1D *_hist91MeanMultiJPsi1S;
    AIDA::IHistogram1D *_hist91MeanMultiPsi2S;
    AIDA::IHistogram1D *_hist91MeanMultiUpsilon1S;
    AIDA::IHistogram1D *_hist91MeanMultiF1_1285;
    AIDA::IHistogram1D *_hist91MeanMultiF1_1420;
    AIDA::IHistogram1D *_hist91MeanMultiChi_c1_3510;
    AIDA::IHistogram1D *_hist91MeanMultiF2_1270;
    AIDA::IHistogram1D *_hist91MeanMultiF2Prime1525;
    AIDA::IHistogram1D *_hist91MeanMultiK2Star1430_0;
    AIDA::IHistogram1D *_hist91MeanMultiBStarStar;
    AIDA::IHistogram1D *_hist91MeanMultiDs1Plus;
    AIDA::IHistogram1D *_hist91MeanMultiDs2Plus;
    AIDA::IHistogram1D *_hist91MeanMultiP;
    AIDA::IHistogram1D *_hist91MeanMultiLambda;
    AIDA::IHistogram1D *_hist91MeanMultiSigma0;
    AIDA::IHistogram1D *_hist91MeanMultiSigmaMinus;
    AIDA::IHistogram1D *_hist91MeanMultiSigmaPlus;
    AIDA::IHistogram1D *_hist91MeanMultiSigmaPlusMinus;
    AIDA::IHistogram1D *_hist91MeanMultiXiMinus;
    AIDA::IHistogram1D *_hist91MeanMultiDelta1232PlusPlus;
    AIDA::IHistogram1D *_hist91MeanMultiSigma1385Minus;
    AIDA::IHistogram1D *_hist91MeanMultiSigma1385Plus;
    AIDA::IHistogram1D *_hist91MeanMultiSigma1385PlusMinus;
    AIDA::IHistogram1D *_hist91MeanMultiXi1530_0;
    AIDA::IHistogram1D *_hist91MeanMultiOmegaMinus;
    AIDA::IHistogram1D *_hist91MeanMultiLambda_c_Plus;
    AIDA::IHistogram1D *_hist91MeanMultiLambda_b_0;
    AIDA::IHistogram1D *_hist91MeanMultiLambda1520;

    AIDA::IHistogram1D *_hist165MeanMultiKPlus;
    AIDA::IHistogram1D *_hist165MeanMultiK0;
    AIDA::IHistogram1D *_hist165MeanMultiP;
    AIDA::IHistogram1D *_hist165MeanMultiLambda;

    //@}

  };

}

#endif
