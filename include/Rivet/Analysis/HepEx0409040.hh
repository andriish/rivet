// -*- C++ -*-
#ifndef RIVET_HepEx0409040_H
#define RIVET_HepEx0409040_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/D0RunIIconeJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/CalMET.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {

  /**
   *  This class is under construction. It will be a reproduction of the HZTool routine
   *  for the ZEUS dijet photoproduction paper which was used in the ZEUS Jets PDF fit.
   *  And hopefully it will teach me how to write a rivet analysis.
   *  @author Jon Butterworth
   */
  class HepEx0409040 : public Analysis {

  public:

    /// Default constructor.
    inline HepEx0409040()
      : fs(-3., 3.), conejets(fs), p_vertex(), p_calmet(fs)
    { }

    /// Copy constructor.
    inline HepEx0409040(const HepEx0409040& x) 
      : Analysis(x), conejets(x.conejets), p_vertex(x.p_vertex), p_calmet(fs)
    { }

    /// Destructor
    inline ~HepEx0409040()
    { }

    /// The name of this analysis is "HepEx0409040"
    inline std::string name() const {
      return "HepEx0409040";
    }

  public:

    void init();
    
    void analyze(const Event & event);
    
    void finalize();

    /// Return the RivetInfo object of this analysis object.
    RivetInfo getInfo() const;

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

    //@{
    /// Histograms
    AIDA::IHistogram1D* histJetAzimuthpTmax75_100;
    AIDA::IHistogram1D* histJetAzimuthpTmax100_130;
    AIDA::IHistogram1D* histJetAzimuthpTmax130_180;
    AIDA::IHistogram1D* histJetAzimuthpTmax180_;
    //@}

  };

}

#endif
