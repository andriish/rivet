// -*- C++ -*-
#ifndef RIVET_HepEx0409040_H
#define RIVET_HepEx0409040_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/D0RunIIconeJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {  

  /// Analysis based on the D0 Run II jet analysis described in hep-ex/0409040.
  /// @author Lars Sonnenschein
  class HepEx0409040 : public Analysis {

  public:

    /// Default constructor.
    inline HepEx0409040()
      : fs(-3.0, 3.0), vfs(fs), vertex()
    { 
      setBeams(PROTON, ANTIPROTON);
      addProjection(fs);
      addProjection(vertex);
      addProjection(vfs);

      vfs.addVetoId(12);
      vfs.addVetoId(14);
      vfs.addVetoId(16);
      vfs.addVetoId(-12);
      vfs.addVetoId(-14);
      vfs.addVetoId(-16);
      // Veto muons with pT above 1.0 GeV
      vfs.addVetoDetail(13, 1.0, numeric_limits<double>::max());

      // Put all particles into the jet finder.
      /// @todo Spot the memory leak!
      conejets = new D0RunIIconeJets(fs);
      addProjection(*conejets);

      // Don't put neutrinos or low pT muons into the cal missing ET.
      /// @todo Spot the memory leak!
      calmet = new TotalVisibleMomentum(vfs);
      addProjection(*calmet);
   }


    /// Return the name of this analysis.
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
    //map<long,double*> vetopids; //12=nu_e, 14=nu_mu, 16=nu_tau, 13=mu

    ///The vetoed final state projector needed by the jet algorithm
    VetoedFinalState vfs; 

    /// The D0RunIIconeJets projector used by this analysis.
    D0RunIIconeJets* conejets;

    /// The Primary Vertex projector
    PVertex vertex;

    /// The Calorimeter Missing Et projector
    TotalVisibleMomentum* calmet;

    /// Hide the assignment operator
    HepEx0409040& operator=(const HepEx0409040& x);

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
