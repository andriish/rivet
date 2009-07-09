// -*- C++ -*-
#ifndef RIVET_DELPHI_1996_S3430090_HH
#define RIVET_DELPHI_1996_S3430090_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /**
   * @brief DELPHI event shapes and identified particle spectra
   * @author Andy Buckley
   * @author Hendrik Hoeth
   *
   * This is the paper which was used for the original PROFESSOR MC tuning
   * study. It studies a wide range of e+ e- event shape variables, differential
   * jet rates in the Durham and JADE schemes, and incorporates identified
   * particle spectra, from other LEP analyses.
   *
   *
   * @par Run conditions
   *
   * @arg LEP1 beam energy: \f$ \sqrt{s} = \$f 91.2 GeV
   * @arg Run with generic QCD events.
   * @arg No \f$ p_\perp^\text{min} \f$ cutoff is required
   */
  class DELPHI_1996_S3430090 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    DELPHI_1996_S3430090(); 

    /// Factory method.
    static Analysis* create() { 
      return new DELPHI_1996_S3430090(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event& event);
    virtual void finalize();
    //@}


  private:

    /// Little helper functions for the axis labels
    string dsigbyd(const string&);
    string Ndsigbyd(const string&);
    string unitdsigbyd(const string&);
    string texmath(const string&);

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the 
    /// inclusive single particle distributions' normalisations.
    double _weightedTotalPartNum;

    double _passedCutWeightSum;

    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_histPtTIn;
    AIDA::IHistogram1D *_histPtTOut;
    AIDA::IHistogram1D *_histPtSIn;
    AIDA::IHistogram1D *_histPtSOut;

    AIDA::IHistogram1D *_histRapidityT;
    AIDA::IHistogram1D *_histRapidityS;

    AIDA::IHistogram1D *_histScaledMom, *_histLogScaledMom;

    AIDA::IProfile1D   *_histPtTOutVsXp, *_histPtVsXp;

    AIDA::IHistogram1D *_hist1MinusT; 
    AIDA::IHistogram1D *_histTMajor; 
    AIDA::IHistogram1D *_histTMinor; 
    AIDA::IHistogram1D *_histOblateness; 

    AIDA::IHistogram1D *_histSphericity;
    AIDA::IHistogram1D *_histAplanarity;
    AIDA::IHistogram1D *_histPlanarity;

    AIDA::IHistogram1D *_histCParam;
    AIDA::IHistogram1D *_histDParam;

    AIDA::IHistogram1D *_histHemiMassD;
    AIDA::IHistogram1D *_histHemiMassH;
    AIDA::IHistogram1D *_histHemiMassL;
               
    AIDA::IHistogram1D *_histHemiBroadW;
    AIDA::IHistogram1D *_histHemiBroadN;
    AIDA::IHistogram1D *_histHemiBroadT;
    AIDA::IHistogram1D *_histHemiBroadD;

    AIDA::IHistogram1D *_histDiffRate2Durham;
    AIDA::IHistogram1D *_histDiffRate2Jade; 
    AIDA::IHistogram1D *_histDiffRate3Durham;
    AIDA::IHistogram1D *_histDiffRate3Jade;
    AIDA::IHistogram1D *_histDiffRate4Durham;
    AIDA::IHistogram1D *_histDiffRate4Jade;

    AIDA::IHistogram1D *_histEEC, *_histAEEC;

    AIDA::IHistogram1D *_histMultiCharged;

    AIDA::IHistogram1D *_histMultiPiPlus;
    AIDA::IHistogram1D *_histMultiPi0;
    AIDA::IHistogram1D *_histMultiKPlus;
    AIDA::IHistogram1D *_histMultiK0;
    AIDA::IHistogram1D *_histMultiEta;
    AIDA::IHistogram1D *_histMultiEtaPrime;
    AIDA::IHistogram1D *_histMultiDPlus;
    AIDA::IHistogram1D *_histMultiD0;
    AIDA::IHistogram1D *_histMultiBPlus0;

    AIDA::IHistogram1D *_histMultiF0;

    AIDA::IHistogram1D *_histMultiRho;
    AIDA::IHistogram1D *_histMultiKStar892Plus;
    AIDA::IHistogram1D *_histMultiKStar892_0;
    AIDA::IHistogram1D *_histMultiPhi;
    AIDA::IHistogram1D *_histMultiDStar2010Plus;

    AIDA::IHistogram1D *_histMultiF2;
    AIDA::IHistogram1D *_histMultiK2Star1430_0;

    AIDA::IHistogram1D *_histMultiP;
    AIDA::IHistogram1D *_histMultiLambda0;
    AIDA::IHistogram1D *_histMultiXiMinus;
    AIDA::IHistogram1D *_histMultiOmegaMinus;
    AIDA::IHistogram1D *_histMultiDeltaPlusPlus;
    AIDA::IHistogram1D *_histMultiSigma1385Plus;
    AIDA::IHistogram1D *_histMultiXi1530_0;
    AIDA::IHistogram1D *_histMultiLambdaB0;
    //@}

  };

}

#endif
