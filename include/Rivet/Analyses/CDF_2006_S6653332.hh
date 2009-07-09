// -*- C++ -*-
#ifndef RIVET_CDF_2006_S6653332_HH
#define RIVET_CDF_2006_S6653332_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /* @brief CDF Run II analysis: jet \f$ p_T \f$ and \f$ \eta \f$ 
   *   distributions in Z + (b) jet production
   * @author Lars Sonnenschein
   *
   * This CDF analysis provides \f$ p_T \f$ and \f$ \eta \f$ distributions of
   * jets in Z + (b) jet production, before and after tagging.
   */
  class CDF_2006_S6653332 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_2006_S6653332();

    /// Factory method.
    static Analysis* create() { 
      return new CDF_2006_S6653332(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}


  private:


    double _Rjet;
    double _JetPtCut;
    double _JetEtaCut;

    double _sumWeightsWithZ;
    double _sumWeightsWithZJet;


    //@{
    /// Histograms
    AIDA::IHistogram1D* _sigmaBJet;
    AIDA::IHistogram1D* _ratioBJetToZ;
    AIDA::IHistogram1D* _ratioBJetToJet;
    //@}

  };


}

#endif
