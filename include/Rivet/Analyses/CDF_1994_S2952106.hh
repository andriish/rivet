// -*- C++ -*-
#ifndef RIVET_CDF_1994_S2952106_HH
#define RIVET_CDF_1994_S2952106_HH

#include "Rivet/Analysis.hh"

namespace Rivet {  

  /* @brief CDF Run I color coherence analysis
   * @author Lars Sonnenschein
   */
  class CDF_1994_S2952106 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_1994_S2952106();

    /// Factory method
    static Analysis* create() { 
      return new CDF_1994_S2952106(); 
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    /// Counter for the number of events analysed
    double _eventsTried;

    /// Counter for the number of  3jet events passed
    double _events3jPassed;


    /// @name Analysis cuts
    //@{
    ///Cut on primary vertex z-position (z(PV) < 60 cm)
    const double _pvzmax;

    /// Min \f$ p_T \f$ of the leading and 3rd leading jets.
    //@{
    const double _leadJetPt;
    const double _3rdJetPt;
    //@}

    /// Max pseudorapidity range of 2nd and 3rd leading jets.
    const double _etamax;

    /// Delta phi (azimuthal angle) requirement (transverse back to back'ness).
    const double _phimin;

    /// MET over sqrt(scalar \f$ E_T \f$) cut requirement.
    const double _metsetmax;
    //@}


  private:

    /// @name Histogram collections
    //@{
    // AIDA::IHistogram2D* _histHvsDphi;
    // AIDA::IHistogram2D* _histRvsAlpha;
    AIDA::IHistogram1D* _histJet1Et;
    AIDA::IHistogram1D* _histJet2Et;
    AIDA::IHistogram1D* _histR23;
    AIDA::IHistogram1D* _histJet3eta;
    AIDA::IHistogram1D* _histAlpha;
    // AIDA::IHistogram1D* _histAlphaMCvsDat;
    AIDA::IHistogram1D* _histAlpaIdeal;
    AIDA::IHistogram1D* _histAlpaCDF;
    AIDA::IHistogram1D* _histR23Ideal;
    AIDA::IHistogram1D* _histR23CDF;
    AIDA::IHistogram1D* _histJet3etaIdeal;
    AIDA::IHistogram1D* _histJet3etaCDF;
    //@}

  };


}

#endif
