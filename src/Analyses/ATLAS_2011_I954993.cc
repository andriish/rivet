// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief MC validation analysis for WZ events
  class ATLAS_2011_I954993 : public Analysis {
  public:

    /// Default constructor
    ATLAS_2011_I954993()
      : Analysis("ATLAS_2011_I954993")
    {
      setNeedsCrossSection(true);
    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {

      //// ZFinder: etaMin,etaMax,pTmin,pid,m2_min,m2_max,dRmax,clusterPhotons,excludePhotonsFromRFS
      ZFinder zfinder_e( -2.5, 2.5, 15.0*GeV, ELECTRON, 81.1876*GeV, 101.1876*GeV, 0.1, true, true);
      addProjection(zfinder_e, "ZFinder_e");
      ZFinder zfinder_mu(-2.5, 2.5, 15.0*GeV, MUON, 81.1876*GeV, 101.1876*GeV, 0.1, true, true);
      addProjection(zfinder_mu, "ZFinder_mu");

      //// WFinder: etaRanges,pTmin,pid,m2_min,m2_max,missingET,dRmax
      VetoedFinalState weinput;
      weinput.addVetoOnThisFinalState(zfinder_e);
      WFinder wfinder_e(weinput, -2.5, 2.5, 15.0*GeV, ELECTRON, 0.0*GeV, 1000.0*GeV, 25.0*GeV, 0.1);
      addProjection(wfinder_e, "WFinder_e");

      VetoedFinalState wminput;
      wminput.addVetoOnThisFinalState(zfinder_mu);
      WFinder wfinder_mu(wminput,-2.5, 2.5, 15.0*GeV, MUON, 0.0*GeV, 1000.0*GeV, 25.0*GeV, 0.1);
      addProjection(wfinder_mu, "WFinder_mu");

      //// Histograms
      _h_fiducial = bookHistogram1D(1,1,1);

    }
    
    /// Do the analysis
    void analyze(const Event & e) {

      const double weight = e.weight();

      const ZFinder& zfinder_e = applyProjection<ZFinder>(e, "ZFinder_e");
      const ZFinder& zfinder_mu = applyProjection<ZFinder>(e, "ZFinder_mu");
      const WFinder& wfinder_e = applyProjection<WFinder>(e, "WFinder_e");
      const WFinder& wfinder_mu = applyProjection<WFinder>(e, "WFinder_mu");
      
    
      // Looking for a Z
      if (zfinder_e.bosons().size()!= 1 && zfinder_mu.bosons().size() != 1) {
        MSG_DEBUG("No Z boson found, vetoing event");
        vetoEvent;
      }

      // Looking for a W
      if (wfinder_e.bosons().size()!= 1 && wfinder_mu.bosons().size() != 1) {
        MSG_DEBUG("No W boson found, vetoing event");
        vetoEvent;
      }

      // If we find a W...
      FourMomentum wmom_e(0.0,0.0,0.0,0.0), We(0.0,0.0,0.0,0.0), Wenu(0.0,0.0,0.0,0.0);
      FourMomentum wmom_mu(0.0,0.0,0.0,0.0), Wmu(0.0,0.0,0.0,0.0), Wmunu(0.0,0.0,0.0,0.0);
      if(wfinder_e.bosons().size()== 1){
	wmom_e = wfinder_e.bosons().front().momentum(); 
	We = wfinder_e.constituentLeptons()[0].momentum();
	Wenu = wfinder_e.constituentNeutrinos()[0].momentum();    
      }
      if(wfinder_mu.bosons().size()== 1){
	wmom_mu = wfinder_mu.bosons().front().momentum(); 
	Wmu = wfinder_mu.constituentLeptons()[0].momentum();
	Wmunu = wfinder_mu.constituentNeutrinos()[0].momentum();
      }

      // Applying remaining fiducial phase space requirements
      double mT = 0;
      if(wfinder_e.bosons().size() == 1){
	mT = sqrt(2*We.pT()*Wenu.Et()*(1.0-cos(We.phi()-Wenu.phi())));
	if (Wenu.pT()/GeV < 25.0 || We.pT()/GeV < 20.0 || mT/GeV < 20.0) {
	  MSG_DEBUG(" Wnu.pT()/GeV:" << Wenu.pT()/GeV<<" Wl.pT()/GeV:" << We.pT()/GeV<<" mT/GeV:" << mT/GeV);
	  vetoEvent;
	}
      }
      else if(wfinder_mu.bosons().size() == 1){
	mT = sqrt(2*Wmu.pT()*Wmunu.Et()*(1.0-cos(Wmu.phi()-Wmunu.phi())));
	if (Wmunu.pT()/GeV < 25.0 || Wmu.pT()/GeV < 20.0 || mT/GeV < 20.0) {
	  MSG_DEBUG(" Wnu.pT()/GeV:" << Wmunu.pT()/GeV<<" Wl.pT()/GeV:" << Wmu.pT()/GeV<<" mT/GeV:" << mT/GeV);
	  vetoEvent;
	}
      }
      else{
	MSG_DEBUG("No W boson found, can't make a transverse mass, vetoing event");
	vetoEvent;
      }

      _h_fiducial->fill(7000.0, weight);

    }


    /// Finalize
    void finalize() {
      
      scale(_h_fiducial, crossSection()/femtobarn/sumOfWeights());

    }

    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _h_fiducial;

    //@}

  };


  //// The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_I954993);

}
