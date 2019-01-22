// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
namespace Rivet {

  /// Primary particles + resonances projection.
  class PrimaryResonances : public Rivet::PrimaryParticles {
    public:
      PrimaryResonances(const Cut& c=Cuts::open())
	: Rivet::PrimaryParticles({},c) {}
      
      virtual int compare(const Projection& p) const {
	return UNDEFINED;
      }
      
      virtual std::unique_ptr<Rivet::Projection> clone() const {
	return std::unique_ptr<Projection>(new PrimaryResonances(*this));
      }
    protected:
      bool isPrimaryPID(const HepMC::GenParticle* p) const
      {
	int pdg = PID::abspid(p->pdg_id());
	if (pdg > 1000000000) return true;
	
	switch (pdg) {
	case Rivet::PID::MUON:
	case Rivet::PID::ELECTRON:
	case Rivet::PID::GAMMA:
	case Rivet::PID::PIPLUS:
	case Rivet::PID::K0S:
	case Rivet::PID::K0L:
	case Rivet::PID::KPLUS:
	case 313: // K*0 (892)
	case 333: // phi
	case Rivet::PID::PROTON:
	case Rivet::PID::NEUTRON:
	case Rivet::PID::LAMBDA:
	case Rivet::PID::SIGMAMINUS:
	case Rivet::PID::SIGMAPLUS:
	case Rivet::PID::XIMINUS:
	case Rivet::PID::XI0:
	case Rivet::PID::OMEGAMINUS:
	case Rivet::PID::NU_E:
	case Rivet::PID::NU_MU:
	case Rivet::PID::NU_TAU:
	  return true;
	}
	return false;
      }
    };

  /// @brief Strangeness enhancement and <pT> in pp 7 TeV, no spectra.
  class ALICE_2018_I1684320 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2018_I1684320);

    int profileIndex(vector<double> cBins, double c) {
      int index = 100;
      if (c > 0 && c <= cBins[0]) return cBins.size() - 1;
      for (size_t i = 0; i < cBins.size() - 1; ++i) {
        if (c > cBins[i] && c <= cBins[i + 1]) {
	  index = i;
	  break;
	} 
      }
      // Catch low fluctuation.
      return max(0, int(cBins.size() - index - 2));
    }

    /// Book histograms and initialise projections before the run
    void init() {
      // Centrality projection.
      declareCentrality(ALICE::V0MMultiplicity(), 
        "ALICE_2015_PPCentrality","V0M","V0M");
      // Central primary particles
      declare(PrimaryResonances(Cuts::abseta < 0.5),"PP");
      declare(PrimaryResonances(Cuts::absrap < 0.5),"PPy");
      centralityBins = {1.,5.,10.,15.,20., 30., 40., 50., 70., 100.};
      centralityBinsPhi = {1.,5.,10.,20., 30., 40., 50., 70., 100.};
      centralityBinsOmega = {5.,15.,30.,50.,100.};
      // Book histograms
      
      pipT = bookProfile1D(1,1,1);
      kpmpT = bookProfile1D(2,1,1);
      ppT = bookProfile1D(3,1,1);
      kzeropT = bookProfile1D(4,1,1);
      lambdapT = bookProfile1D(5,1,1);
      kstarpT = bookProfile1D(6,1,1);
      phipT = bookProfile1D(7,1,1);
      xipT = bookProfile1D(8,1,1);
      omegapT = bookProfile1D(9,1,1);

      piYield = bookProfile1D(1,2,1);
      kpmYield = bookProfile1D(2,2,1);
      kstarYield = bookProfile1D(3,2,1);
      phiYield = bookProfile1D(4,2,1);
      piRebinnedK = shared_ptr<YODA::Profile1D>(kstarYield->newclone());
      piRebinnedK->setTitle("piRebinnedK");
      piRebinnedK->setPath("/" + name() + "/piRebinnedK");
      addAnalysisObject(piRebinnedK);

      piRebinnedP = shared_ptr<YODA::Profile1D>(phiYield->newclone());
      piRebinnedP->setTitle("piRebinnedP");
      piRebinnedP->setPath("/" + name() + "/piRebinnedP");
      addAnalysisObject(piRebinnedP);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (apply<PrimaryResonances>(event,"PP").particles().size() < 1) vetoEvent;
      const PrimaryResonances& prim = apply<PrimaryResonances>(event,"PPy");
      const double weight = event.weight();
      const CentralityProjection& cent = apply<CentralityProjection>(event,"V0M");
      double c  = cent();
      // Fill the pt histograms and count yields.
      int npi = 0, nkpm = 0, nkstar = 0, nphi = 0;
      // Profile centrality index. 
      int index = profileIndex(centralityBins,c);
      int indexPhi = profileIndex(centralityBinsPhi, c);
      int indexOmega = profileIndex(centralityBinsOmega, c);
      for (auto p : prim.particles()) {
        const double pT = p.pT();
	const int pid = abs(p.pid());
	if (pid == 211) {
	  ++npi;
	  pipT->fillBin(index, pT, weight);
	}
	else if (pid == 321) {
	  ++nkpm;
	  kpmpT->fillBin(index, pT, weight);
	}
	else if (pid == 2212) {
	  ppT->fillBin(index, pT, weight);
	}
	else if (pid == 310) {
	  kzeropT->fillBin(index, pT, weight);
	}
	else if (pid == 3122) {
	  lambdapT->fillBin(index, pT, weight);
	}
	else if (pid == 313) {
	  ++nkstar;
	  kstarpT->fillBin(indexPhi, pT, weight);
	}
	else if (pid == 333) {
	  ++nphi;
	  phipT->fillBin(indexPhi, pT, weight);
	}
	else if (pid == 3312) {
	  xipT->fillBin(index, pT, weight);
	}
	else if (pid == 3334) {
	  omegapT->fillBin(indexOmega, pT, weight);
	}
      }
      // Fill the profiles of yields.
      piYield->fillBin(index, double(npi)/2., weight);
      kpmYield->fillBin(index, double(nkpm)/2., weight);
      kstarYield->fillBin(indexPhi, double(nkstar)/2., weight);
      phiYield->fillBin(indexPhi, 2.*double(nphi), weight);
      piRebinnedK->fillBin(indexPhi,double(npi)/2.,weight);
      piRebinnedP->fillBin(indexPhi,double(npi)/2.,weight);
    
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Make the ratios
      kpmpi = bookScatter2D(1, 1, 2, true);
      kstarpi = bookScatter2D(2, 1, 2, true);
      phipi = bookScatter2D(3, 1, 2, true);

      divide(kpmYield, piYield, kpmpi);
      divide(kstarYield, piRebinnedK, kstarpi);
      divide(phiYield, piRebinnedP, phipi);
    }

    //@}


    /// @name Histograms
    //@{
    // Histograms ordered in centrality classes
    vector<double> centralityBins;
    vector<double> centralityBinsPhi;
    vector<double> centralityBinsOmega;

    // Average pT
    Profile1DPtr pipT;
    Profile1DPtr kpmpT;
    Profile1DPtr ppT;
    Profile1DPtr kzeropT;
    Profile1DPtr lambdapT;
    Profile1DPtr kstarpT;
    Profile1DPtr phipT;
    Profile1DPtr xipT;
    Profile1DPtr omegapT;

    // Total yields
    Profile1DPtr piYield;
    Profile1DPtr kpmYield;
    Profile1DPtr kstarYield;
    Profile1DPtr phiYield;
    Profile1DPtr piRebinnedK;
    Profile1DPtr piRebinnedP;

    // Ratios
    Scatter2DPtr kpmpi;
    Scatter2DPtr kstarpi;
    Scatter2DPtr phipi;
    //@}
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2018_I1684320);


}
