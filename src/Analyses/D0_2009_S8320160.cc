// -*- C++ -*-
#include "Rivet/Analyses/D0_2009_S8320160.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2009_S8320160::D0_2009_S8320160()
    : Analysis("D0_2009_S8320160")
  {
    setBeams(PROTON, ANTIPROTON);
    
    ChargedFinalState cfs;
    FastJets conefinder(cfs, FastJets::D0ILCONE, 0.7);
    addProjection(conefinder, "ConeFinder");
  } 


  // Book histograms
  void D0_2009_S8320160::init() {
    _h_chi_dijet.addHistogram(250., 300., bookHistogram1D(1, 1, 1));
    _h_chi_dijet.addHistogram(300., 400., bookHistogram1D(2, 1, 1));
    _h_chi_dijet.addHistogram(400., 500., bookHistogram1D(3, 1, 1));
    _h_chi_dijet.addHistogram(500., 600., bookHistogram1D(4, 1, 1));
    _h_chi_dijet.addHistogram(600., 700., bookHistogram1D(5, 1, 1));
    _h_chi_dijet.addHistogram(700., 800., bookHistogram1D(6, 1, 1));
    _h_chi_dijet.addHistogram(800., 900., bookHistogram1D(7, 1, 1));
    _h_chi_dijet.addHistogram(900., 1000., bookHistogram1D(8, 1, 1));
    _h_chi_dijet.addHistogram(1000., 1100., bookHistogram1D(9, 1, 1));
    _h_chi_dijet.addHistogram(1100., 1960, bookHistogram1D(10, 1, 1));
  }



  // Do the analysis 
  void D0_2009_S8320160::analyze(const Event & e) {
    double weight = e.weight();

    const JetAlg& jetpro = applyProjection<JetAlg>(e, "ConeFinder");
    const Jets& jets = jetpro.jetsByPt();
    
    if(jets.size()<2) {
      vetoEvent;
    }
    
    FourMomentum j0(jets[0].momentum());
    FourMomentum j1(jets[1].momentum());
    double mjj = FourMomentum(j0+j1).mass();
    double chi = exp(fabs(j0.rapidity()-j1.rapidity()));
    _h_chi_dijet.fill(mjj, chi, weight);
  }



  // Finalize
  void D0_2009_S8320160::finalize() {
    foreach (AIDA::IHistogram1D* hist, _h_chi_dijet.getHistograms()) {
      scale(hist, 1.0/sumOfWeights());
    }
  }

}
