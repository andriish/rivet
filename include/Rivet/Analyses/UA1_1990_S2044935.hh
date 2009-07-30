// -*- C++ -*-
#ifndef RIVET_UA1_1990_S2044935_HH
#define RIVET_UA1_1990_S2044935_HH

#include "Rivet/Analysis.hh"

namespace Rivet {
  
  class UA1_1990_S2044935 : public Analysis {

  public:

    /// @name Constructor etc.
    //@{

    /// constructor
    UA1_1990_S2044935();

    /// Factory method
    static Analysis* create() { 
      return new UA1_1990_S2044935(); 
    }

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:
    
    /// @name Histogram collections
    //@{
    AIDA::IHistogram1D* _hist_sigma200;
    AIDA::IHistogram1D* _hist_sigma500;
    AIDA::IHistogram1D* _hist_sigma900;
    AIDA::IHistogram1D* _hist_Esigma200;
    AIDA::IHistogram1D* _hist_Esigma500;
    AIDA::IHistogram1D* _hist_Esigma900;
    AIDA::IHistogram1D* _hist_Esigmapoint8;
    AIDA::IHistogram1D* _hist_Esigma4;
    AIDA::IHistogram1D* _hist_Esigma8;
        AIDA::IProfile1D* _hist_Pt63;
    AIDA::IProfile1D* _hist_Pt200;
    AIDA::IProfile1D* _hist_Pt900;
    AIDA::IProfile1D* _hist_Etavg200;
    AIDA::IProfile1D* _hist_Etavg500;
    AIDA::IProfile1D* _hist_Etavg900;
    AIDA::IHistogram1D* _hist_Et200;
    AIDA::IHistogram1D* _hist_Et500;
    AIDA::IHistogram1D* _hist_Et900;
    //@}

  };

}

#endif

