// -*- C++ -*-
#ifndef RIVET_CDF_2008_S7782535_HH
#define RIVET_CDF_2008_S7782535_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// Implementation of CDF RunII b-jet shape paper
  class CDF_2008_S7782535 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_2008_S7782535();

    /// Factory method
    static Analysis* create() { 
      return new CDF_2008_S7782535(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event & event);
    void finalize();
    //@}


  private:

    /// @name Analysis data
    //@{
    vector<FourMomentum> _jetaxes;
    double _Rjet;
    vector<double> _pTbins;
    int _NpTbins;
    //@}


    /// @name Histograms
    //@{
    AIDA::IProfile1D* _h_Psi_pT[4];
    AIDA::IDataPointSet* _h_OneMinusPsi_vs_pT;
    //@}

  };

}

#endif
