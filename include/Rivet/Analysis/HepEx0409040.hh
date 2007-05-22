// -*- C++ -*-
#ifndef RIVET_HepEx0409040_H
#define RIVET_HepEx0409040_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/D0RunIIconeJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/CalMET.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {

  /// Analysis based on the D0 Run II jet analysis described in hep-ex/0409040.
  /// @author Lars Sonnenschein
  class HepEx0409040 : public Analysis {

  public:

    /// Default constructor.
    inline HepEx0409040(bool xrad=true)
      : fs(-3.0, 3.0), conejets(fs), p_vertex(), p_calmet(fs)
    { 
      setBeams(PROTON, ANTIPROTON);
      addProjection(fs);
      addProjection(conejets);
      addProjection(p_vertex);
      addProjection(p_calmet);
      /// @todo Remove - the analysis must be absolutely repeatable.
      if (xrad) 
        xscale = PI/128.; //histo x-axis range in [rad]
      else
        xscale = 1.; //histo x-axis range in [rad*128/Pi]
    }


    /// The name of this analysis is "HepEx0409040"
    inline string getName() const {
      return "HepEx0409040";
    }

  public:

    void init();
    
    void analyze(const Event & event);
    
    void finalize();

    /// Return the RivetInfo object of this analysis object.
    //RivetInfo getInfo() const;

  private:

    /// The FinalState projector used by this analysis.
    FinalState fs;

    /// The D0RunIIconeJets projector used by this analysis.
    D0RunIIconeJets conejets;

    /// The Primary Vertex projector
    PVertex p_vertex;

    /// The Calorimeter Missing Et projector
    CalMET p_calmet;

    /// Hide the assignment operator
    HepEx0409040 & operator=(const HepEx0409040& x);

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* histJetAzimuth_pTmax75_100;
    AIDA::IHistogram1D* histJetAzimuth_pTmax100_130;
    AIDA::IHistogram1D* histJetAzimuth_pTmax130_180;
    AIDA::IHistogram1D* histJetAzimuth_pTmax180_;
    //@}
    double xscale;

  };

}

#endif
