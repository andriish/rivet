// -*- C++ -*-
#ifndef RIVET_UA5_1989_S1926373_HH
#define RIVET_UA5_1989_S1926373_HH

#include "Rivet/Analysis.hh"

namespace Rivet {

    
  class UA5_1989_S1926373 : public Analysis {

  public:

    /// @name Constructor etc.
    //@{

    /// constructor
    UA5_1989_S1926373();

    /// Factory method
    static Analysis* create() { 
      return new UA5_1989_S1926373(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:
    
    //@{
    /// Histograms
    AIDA::IHistogram1D* _hist_nch200;
    AIDA::IHistogram1D* _hist_nch900;
    AIDA::IHistogram1D* _hist_nch200eta0point5;
    AIDA::IHistogram1D* _hist_nch200eta1point5;
    AIDA::IHistogram1D* _hist_nch200eta3;
    AIDA::IHistogram1D* _hist_nch200eta5;
    AIDA::IHistogram1D* _hist_nch900eta0point5;
    AIDA::IHistogram1D* _hist_nch900eta1point5;
    AIDA::IHistogram1D* _hist_nch900eta3;
    AIDA::IHistogram1D* _hist_nch900eta5;
    AIDA::IHistogram1D* _hist_mean_nch_200;
    AIDA::IHistogram1D* _hist_mean_nch_900;
    //@}

  };

}

#endif
