// -*- C++ -*-
#ifndef RIVET_D0_2006_S6438750_HH
#define RIVET_D0_2006_S6438750_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Inclusive isolated photon cross-section, differential in \f$ p_\perp(gamma) \f$.
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2006_S6438750 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Default constructor.
     D0_2006_S6438750();

    /// Factory method 
    static Analysis* create() {
      return new D0_2006_S6438750();
    }

    //@}


    /// @name Analysis methods
    //@{ 
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _h_pTgamma;
    //@}

  };

}

#endif
