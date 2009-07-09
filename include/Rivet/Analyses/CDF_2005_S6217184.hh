// -*- C++ -*-
#ifndef RIVET_CDF_2005_S6217184_HH
#define RIVET_CDF_2005_S6217184_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /* CDF Run II jet shape analysis
   * @author Lars Sonnenschein
   * @author Andy Buckley
   */	
  class CDF_2005_S6217184 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    // Constructor
    CDF_2005_S6217184();

    /// Factory method
    static Analysis* create() { 
      return new CDF_2005_S6217184(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}


  private:

    /// @name Analysis data
    //@{

    /// Vector of jet axes
    vector<FourMomentum> _jetaxes;

    /// \f$p_\perp\f$ bins to be distinguished during analysis
    vector<double> _pTbins;
    //@}


    /// @name Histograms
    //@{
    AIDA::IProfile1D* _profhistRho_pT[18];
    AIDA::IProfile1D* _profhistPsi_pT[18];
    AIDA::IProfile1D* _profhistPsi;
    //@}

  };


}

#endif
