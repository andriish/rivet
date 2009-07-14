// -*- C++ -*-
#include "Rivet/Analyses/CDF_2008_S8093652.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  CDF_2008_S8093652::CDF_2008_S8093652()
    : Analysis("CDF_2008_S8093652")
  {
    setBeams(PROTON, ANTIPROTON);
    setNeedsCrossSection(true);
    
    ChargedFinalState cfs;
    FastJets conefinder(cfs, FastJets::CDFMIDPOINT, 0.7);
    addProjection(conefinder, "ConeFinder");
  } 


  // Book histograms
  void CDF_2008_S8093652::init() {
    _h_m_dijet = bookHistogram1D(1, 1, 1);
  }



  // Do the analysis 
  void CDF_2008_S8093652::analyze(const Event & e) {
    double weight = e.weight();

    const JetAlg& jetpro = applyProjection<JetAlg>(e, "ConeFinder");
    const Jets& jets = jetpro.jetsByPt();
    
    if(jets.size()<2) {
      vetoEvent;
    }
    
    FourMomentum j0(jets[0].momentum());
    FourMomentum j1(jets[1].momentum());
    
    if (fabs(j1.rapidity())>1.0 || fabs(j0.rapidity())>1.0) {
      vetoEvent;
    }
    
    double mjj = FourMomentum(j0+j1).mass();
    
    _h_m_dijet->fill(mjj, weight);
  }



  // Finalize
  void CDF_2008_S8093652::finalize() {
    /// Scale by L_eff = sig_MC * L_exp / num_MC
    const double lumi_mc = sumOfWeights() / crossSection();
    const double scalefactor =  1 / lumi_mc;
    scale(_h_m_dijet, scalefactor);
  }

}
