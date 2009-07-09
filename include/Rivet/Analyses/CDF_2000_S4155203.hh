// -*- C++ -*-
#ifndef RIVET_CDF_2000_S4155203_HH
#define RIVET_CDF_2000_S4155203_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /*
   * @brief CDF Run I Z pT in Drell-Yan events
   * @author Hendrik Hoeth
   */ 
  class CDF_2000_S4155203 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.5 \f$ GeV.
    CDF_2000_S4155203();

    /// Factory method
    static Analysis* create() {
      return new CDF_2000_S4155203();
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}


  private:

    AIDA::IHistogram1D *_hist_zpt;

  };


}

#endif
