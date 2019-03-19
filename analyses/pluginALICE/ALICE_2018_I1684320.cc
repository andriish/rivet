// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"
namespace Rivet {

  /// Used resonances.
  class PrimaryPK : public Rivet::PrimaryParticles {
    public:
      PrimaryPK(const Cut& c=Cuts::open())
	: Rivet::PrimaryParticles({},c) {}
      
      virtual int compare(const Projection& p) const {
	return UNDEFINED;
      }
      
      virtual std::unique_ptr<Rivet::Projection> clone() const {
	return std::unique_ptr<Projection>(new PrimaryPK(*this));
      }
    protected:
      bool isPrimaryPID(const HepMC::GenParticle* p) const
      {
	int pdg = PID::abspid(p->pdg_id());
	if (pdg > 1000000000) return true;
	
	switch (pdg) {
	case 313: // K*0 (892)
	case 333: // phi
	  return true;
	}
	return false;
      }
    };
  /// @brief Strangeness spectra
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
      return max(0, int(cBins.size() - index - 2));
    }

    /// Book histograms and initialise projections before the run
    void init() {
      // Centrality projection.
      declareCentrality(ALICE::V0MMultiplicity(), 
        "ALICE_2015_PPCentrality","V0M","V0M");
      // Central primary particles
      declare(ChargedFinalState(Cuts::abseta < 1.0),"PP");
      declare(ALICE::PrimaryParticles(Cuts::absrap < 0.5),"PPy");
      declare(PrimaryPK(Cuts::absrap < 0.5),"PRy");
      centralityBins = {1.,5.,10.,15.,20., 30., 40., 50., 70., 100.};
      centralityBinsPhi = {1.,5.,10.,20., 30., 40., 50., 70., 100.};
      
      // Book histograms
      for (int i = 0; i < 10; ++i) {
	chargedpT[centralityBins[i]] = bookHisto1D(i+1,1,1);
	pipT[centralityBins[i]] = bookHisto1D(i+11,1,1);
	KpmpT[centralityBins[i]] = bookHisto1D(i+21,1,1);
	ppT[centralityBins[i]] = bookHisto1D(i+31,1,1);
	sow[centralityBins[i]] = bookCounter("sow_" + toString(i));
      }
      for (int i = 0; i < 9; ++i) {
	KstarpT[centralityBinsPhi[i]] = bookHisto1D(i+41,1,1);
	phipT[centralityBinsPhi[i]] = bookHisto1D(i+49,1,1);
	sowPhi[centralityBinsPhi[i]] = bookCounter("sowPhi_" + toString(i));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (apply<ChargedFinalState>(event,"PP").particles().size() < 1) vetoEvent;
      const ALICE::PrimaryParticles& prim = apply<ALICE::PrimaryParticles>(event,"PPy");
      const PrimaryPK& primR = apply<PrimaryPK>(event,"PRy");
      const double weight = event.weight();
      const CentralityProjection& cent = apply<CentralityProjection>(event,"V0M");
      double c  = cent();
      // Find the correct histograms
      auto cItr = chargedpT.upper_bound(c);
      if (cItr == chargedpT.end()) return;
      auto piItr = pipT.upper_bound(c);
      if (piItr == pipT.end()) return;
      auto kItr = KpmpT.upper_bound(c);
      if (kItr == KpmpT.end()) return;
      auto pItr = ppT.upper_bound(c);
      if (pItr == ppT.end()) return;
      auto ksItr = KstarpT.upper_bound(c);
      if (ksItr == KstarpT.end()) return;
      auto phiItr = phipT.upper_bound(c);
      if (phiItr == phipT.end()) return;
      // Fill the sow.
      auto sowItr = sow.upper_bound(c);
      if (sowItr == sow.end()) return;
      auto sowPhiItr = sowPhi.upper_bound(c);
      if (sowPhiItr == sowPhi.end()) return;
      sowItr->second->fill(weight);
      sowPhiItr->second->fill(weight);
      // Fill the pt histograms.
      for (auto p : prim.particles()) {
        const double pT = p.pT();
	const int pid = abs(p.pid());
	if (p.charge() != 0) cItr->second->fill(pT, weight);
	if (pid == 211) piItr->second->fill(pT, weight);
        else if (pid == 321) kItr->second->fill(pT, weight);
	else if (pid == 2212) pItr->second->fill(pT, weight);
      }
      for (auto p : primR.particles()) {
        const double pT = p.pT();
	const int pid = abs(p.pid());
        if (pid == 313) ksItr->second->fill(pT, weight);
	else if (pid == 333) phiItr->second->fill(pT, weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Normalize the spectra
      for (int i = 0; i < 10; ++i) {
        chargedpT[centralityBins[i]]->scaleW(1./sow[centralityBins[i]]->sumW());
        pipT[centralityBins[i]]->scaleW(1./sow[centralityBins[i]]->sumW());
        KpmpT[centralityBins[i]]->scaleW(1./sow[centralityBins[i]]->sumW());
        ppT[centralityBins[i]]->scaleW(1./sow[centralityBins[i]]->sumW());
      }
      for (int i = 0; i < 9; ++i) {
        KstarpT[centralityBinsPhi[i]]->scaleW(1./sowPhi[centralityBinsPhi[i]]->sumW());
        phipT[centralityBinsPhi[i]]->scaleW(1./sowPhi[centralityBinsPhi[i]]->sumW());
      }

    }

    //@}


    /// @name Histograms
    //@{
    // Histograms ordered in centrality classes
    vector<double> centralityBins;
    vector<double> centralityBinsPhi;

    // pT spectra
    map<double, Histo1DPtr> chargedpT;
    map<double, Histo1DPtr> pipT;
    map<double, Histo1DPtr> KpmpT;
    map<double, Histo1DPtr> ppT;
    map<double, Histo1DPtr> KstarpT;
    map<double, Histo1DPtr> phipT;
    map<double, CounterPtr> sow;
    map<double, CounterPtr> sowPhi;

    //@}
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2018_I1684320);


}
