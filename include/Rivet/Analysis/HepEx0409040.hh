// -*- C++ -*-
#ifndef RIVET_HepEx0409040_H
#define RIVET_HepEx0409040_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/D0RunIIconeJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/CalMET.hh"
#include "Rivet/RivetAIDA.fhh"

//#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {

  

  /// Analysis based on the D0 Run II jet analysis described in hep-ex/0409040.
  /// @author Lars Sonnenschein
  class HepEx0409040 : public Analysis {

  public:

    /// Default constructor.
    inline HepEx0409040()
      //: fs(-3.0, 3.0), conejets(fs), vertex(), calmet(fs)
      //: fs(-3.0, 3.0), vfs(fs, vetopids), conejets(fs, vfs), vertex(), calmet(vfs) //ls
      : fs(-3.0, 3.0), vertex(), calmet(fs) //ls
    { 
      //veto pids: 12=nu_e, 14=nu_mu, 16=nu_tau, 13=mu
      /*
      vetopids.insert(-12);
      vetopids.insert(12);
      vetopids.insert(-14);
      vetopids.insert(14);
      vetopids.insert(-16);
      vetopids.insert(16);
      vetopids.insert(-13);
      vetopids.insert(13);
      */

      setBeams(PROTON, ANTIPROTON);
      addProjection(fs);
      addProjection(vertex);

      vfs = new VetoedFinalState(fs, vetopids);
      addProjection(*vfs); //ls

      conejets = new D0RunIIconeJets(*vfs);
      addProjection(*conejets);

      //calmet = new CalMET(fs);
      addProjection(calmet);
 

   }


    /// The name of this analysis is "HepEx0409040"
    inline string getName() const {
      return "HepEx0409040";
    }

  public:

    void init();
    
    void analyze(const Event & event);
    
    void finalize();

  private:

    /// The final state projector used by this analysis.
    FinalState fs;

    ///The veto against final state particles
    set<long> vetopids; //12=nu_e, 14=nu_mu, 16=nu_tau, 13=mu

    ///The vetoed final state projector needed by the jet algorithm
    //VetoedFinalState vfs(FinalState fs, set<long> vetopids); //ls
    VetoedFinalState* vfs; //ls



    /// The D0RunIIconeJets projector used by this analysis.
    D0RunIIconeJets* conejets;

    /// The Primary Vertex projector
    PVertex vertex;

    /// The Calorimeter Missing Et projector
    CalMET calmet;

    /// Hide the assignment operator
    HepEx0409040 & operator=(const HepEx0409040& x);

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* histJetAzimuth_pTmax75_100;
    AIDA::IHistogram1D* histJetAzimuth_pTmax100_130;
    AIDA::IHistogram1D* histJetAzimuth_pTmax130_180;
    AIDA::IHistogram1D* histJetAzimuth_pTmax180_;
    //@}

  };

}

#endif
