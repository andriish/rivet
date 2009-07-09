// -*- C++ -*-
#ifndef RIVET_CDF_1990_S2089246_HH
#define RIVET_CDF_1990_S2089246_HH

#include "Rivet/Analysis.hh"

namespace Rivet {  

  /* @brief CDF pseudorapidity analysis
   * @author Andy Buckley
   */
  class CDF_1990_S2089246 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_1990_S2089246();

    /// Factory method
    static Analysis* create() { 
      return new CDF_1990_S2089246(); 
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}


  private:

    /// @name Histogram collections
    //@{
    AIDA::IHistogram1D* _hist_eta630;
    AIDA::IHistogram1D* _hist_eta1800;
    //@}

  };


}

#endif
