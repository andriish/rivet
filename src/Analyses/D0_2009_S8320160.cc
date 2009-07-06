// -*- C++ -*-
#include "Rivet/Analyses/D0_2009_S8320160.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2009_S8320160::D0_2009_S8320160()
  {
    setBeams(PROTON, ANTIPROTON);
    
    ChargedFinalState cfs;
    FastJets conefinder(cfs, FastJets::D0ILCONE, 0.7);
    addProjection(conefinder, "ConeFinder");
  } 


  // Book histograms
  void D0_2009_S8320160::init() {
    _h_chi_dijet.addHistogram(0.25, 0.3, bookHistogram1D(
        1, 1, 1, "$0.25 \\leq M_{jj}/\\text{TeV} \\leq 0.3$",
        "\\chi_{\\text{dijet}}=\\exp(|y_1 - y_2|)",
        "1/\\sigma_{\\text{dijet}} d\\sigma/d\\chi_{\\text{dijet}}"));
    _h_chi_dijet.addHistogram(0.3, 0.4, bookHistogram1D(
        2, 1, 1, "$0.3 \\leq M_{jj}/\\text{TeV} \\leq 0.4$",
        "\\chi_{\\text{dijet}}=\\exp(|y_1 - y_2|)",
        "1/\\sigma_{\\text{dijet}} d\\sigma/d\\chi_{\\text{dijet}}"));
    _h_chi_dijet.addHistogram(0.4, 0.5, bookHistogram1D(
        3, 1, 1, "$0.4 \\leq M_{jj}/\\text{TeV} \\leq 0.5$",
        "\\chi_{\\text{dijet}}=\\exp(|y_1 - y_2|)",
        "1/\\sigma_{\\text{dijet}} d\\sigma/d\\chi_{\\text{dijet}}"));
    _h_chi_dijet.addHistogram(0.5, 0.6, bookHistogram1D(
        4, 1, 1, "$0.5 \\leq M_{jj}/\\text{TeV} \\leq 0.6$",
        "\\chi_{\\text{dijet}}=\\exp(|y_1 - y_2|)",
        "1/\\sigma_{\\text{dijet}} d\\sigma/d\\chi_{\\text{dijet}}"));
    _h_chi_dijet.addHistogram(0.6, 0.7, bookHistogram1D(
        5, 1, 1, "$0.6 \\leq M_{jj}/\\text{TeV} \\leq 0.7$",
        "\\chi_{\\text{dijet}}=\\exp(|y_1 - y_2|)",
        "1/\\sigma_{\\text{dijet}} d\\sigma/d\\chi_{\\text{dijet}}"));
    _h_chi_dijet.addHistogram(0.7, 0.8, bookHistogram1D(
        6, 1, 1, "$0.7 \\leq M_{jj}/\\text{TeV} \\leq 0.8$",
        "\\chi_{\\text{dijet}}=\\exp(|y_1 - y_2|)",
        "1/\\sigma_{\\text{dijet}} d\\sigma/d\\chi_{\\text{dijet}}"));
    _h_chi_dijet.addHistogram(0.8, 0.9, bookHistogram1D(
        7, 1, 1, "$0.8 \\leq M_{jj}/\\text{TeV} \\leq 0.9$",
        "\\chi_{\\text{dijet}}=\\exp(|y_1 - y_2|)",
        "1/\\sigma_{\\text{dijet}} d\\sigma/d\\chi_{\\text{dijet}}"));
    _h_chi_dijet.addHistogram(0.9, 1.0, bookHistogram1D(
        8, 1, 1, "$0.9 \\leq M_{jj}/\\text{TeV} \\leq 1.0$",
        "\\chi_{\\text{dijet}}=\\exp(|y_1 - y_2|)",
        "1/\\sigma_{\\text{dijet}} d\\sigma/d\\chi_{\\text{dijet}}"));
    _h_chi_dijet.addHistogram(1.0, 1.1, bookHistogram1D(
        9, 1, 1, "$1.0 \\leq M_{jj}/\\text{TeV} \\leq 1.1$",
        "\\chi_{\\text{dijet}}=\\exp(|y_1 - y_2|)",
        "1/\\sigma_{\\text{dijet}} d\\sigma/d\\chi_{\\text{dijet}}"));
    _h_chi_dijet.addHistogram(1.1, 1.960, bookHistogram1D(
        10, 1, 1, "$1.1 \\leq M_{jj}/\\text{TeV}$",
        "\\chi_{\\text{dijet}}=\\exp(|y_1 - y_2|)",
        "1/\\sigma_{\\text{dijet}} d\\sigma/d\\chi_{\\text{dijet}}"));
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
      normalize(hist);
    }
  }

}
