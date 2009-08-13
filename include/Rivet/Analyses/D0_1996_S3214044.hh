// -*- C++ -*-
#ifndef RIVET_D0_1996_S3214044_HH
#define RIVET_D0_1996_S3214044_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Jet.hh"

namespace Rivet {


  class D0_1996_S3214044 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    D0_1996_S3214044();

    /// Factory method
    static Analysis* create() {
      return new D0_1996_S3214044();
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

    void threeJetAnalysis(const Jets& jets, const double& weight);
    void fourJetAnalysis(const Jets& jets, const double& weight);


  private:

    /// @name Histograms
    //@{

    AIDA::IHistogram1D *_h_3j_x3;
    AIDA::IHistogram1D *_h_3j_x5;
    AIDA::IHistogram1D *_h_3j_costheta3;
    AIDA::IHistogram1D *_h_3j_psi;
    AIDA::IHistogram1D *_h_3j_mu34;
    AIDA::IHistogram1D *_h_3j_mu35;
    AIDA::IHistogram1D *_h_3j_mu45;
    
    AIDA::IHistogram1D *_h_4j_x3;
    AIDA::IHistogram1D *_h_4j_x4;
    AIDA::IHistogram1D *_h_4j_x5;
    AIDA::IHistogram1D *_h_4j_x6;
    AIDA::IHistogram1D *_h_4j_costheta3;
    AIDA::IHistogram1D *_h_4j_costheta4;
    AIDA::IHistogram1D *_h_4j_costheta5;
    AIDA::IHistogram1D *_h_4j_costheta6;
    AIDA::IHistogram1D *_h_4j_cosomega34;
    AIDA::IHistogram1D *_h_4j_cosomega35;
    AIDA::IHistogram1D *_h_4j_cosomega36;
    AIDA::IHistogram1D *_h_4j_cosomega45;
    AIDA::IHistogram1D *_h_4j_cosomega46;
    AIDA::IHistogram1D *_h_4j_cosomega56;
    AIDA::IHistogram1D *_h_4j_mu34;
    AIDA::IHistogram1D *_h_4j_mu35;
    AIDA::IHistogram1D *_h_4j_mu36;
    AIDA::IHistogram1D *_h_4j_mu45;
    AIDA::IHistogram1D *_h_4j_mu46;
    AIDA::IHistogram1D *_h_4j_mu56;
    AIDA::IHistogram1D *_h_4j_theta_BZ;
    AIDA::IHistogram1D *_h_4j_costheta_NR;
    //@}

  };


}

#endif

