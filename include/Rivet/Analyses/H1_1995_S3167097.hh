// -*- C++ -*-
#ifndef RIVET_H1_1995_S3167097_HH
#define RIVET_H1_1995_S3167097_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Measures energy flow in DIS? To be checked!
  /// @todo Check this analysis!
  /// @author Leif Lonnblad
  class H1_1995_S3167097 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor.
    H1_1995_S3167097();

    /// Factory method
    static Analysis* create() { 
      return new H1_1995_S3167097(); 
    }

    //@}    


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

    
  private:

    /// Some integer constants used.
    static const size_t _nb = 24, _nbin = 9;
    
    /// Some double constants used.
    static const double _xmin, _xmax;

    /// Histograms for the \f$ E_T \f$ flows
    vector<AIDA::IHistogram1D*> _hEtFlow, _hEtFlowStat;

    /// Histograms for averages in different kinematical bins.
    AIDA::IHistogram1D *_hAvEt, *_hAvX, *_hAvQ2, *_hN;

    /// Helper vector;
    vector<double> _nev;

  };

}

#endif
