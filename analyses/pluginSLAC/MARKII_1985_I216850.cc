// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief 3 jet events at 29 GeV
  class MARKII_1985_I216850 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1985_I216850);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      const FinalState fs;
      declare(fs, "FS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      Sphericity sphere(fs);
      declare(sphere, "Sphericity");
      FastJets jadeJets = FastJets(fs, FastJets::JADE, 0.7, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      declare(jadeJets, "JadeJets");
      // histos
      book(_h,1,1,1);
      book(_njet,"TMP/njet");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // at least 5 charged particles
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      if(cfs.particles().size()<5) vetoEvent;
      // sphericity cuts
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      if(sphericity.lambda3()<0.06 || sphericity.lambda2()-sphericity.lambda3()<0.05) vetoEvent;
      // cluster particles using jet
      const FastJets& jadejet = apply<FastJets>(event, "JadeJets");
      Particles fs = apply<FinalState>(event, "FS").particles();
      PseudoJets pjets = jadejet.clusterSeq()->exclusive_jets(3);
      Jets jets = jadejet.mkJets(pjets,fs);
      Vector3 n[3];
      // loop over the jets
      for(unsigned int ix=0;ix<jets.size();++ix) {
	if(jets[ix].constituents().size()<3) vetoEvent;
	double energy=0.;
	for(const Particle & p : jets[ix].constituents()) energy+=p.E();
	if(energy<2.*GeV) vetoEvent;
	n[ix] = jets[ix].mom().p3().unit();
      }
      vector<double> dot; dot.reserve(3);
      for(unsigned int ix=0;ix<jets.size();++ix) {
	for(unsigned int iy=ix+1;iy<jets.size();++iy) {
	  dot.push_back(n[ix].dot(n[iy]));
	  double angle = acos(dot.back())*180./M_PI;
	  if(angle<100 || angle>140) vetoEvent;
	}
      }
      double Ejet[3];
      Ejet[0] = (sqrtS()*(dot[1] - dot[0]*dot[2]))/((-1 + dot[0])*(1 + dot[0] - dot[1] - dot[2]));
      Ejet[1] = (sqrtS()*(-(dot[0]*dot[1]) + dot[2]))/((-1 + dot[0])*(1 + dot[0] - dot[1] - dot[2]));
      Ejet[2] = (sqrtS()*(1 + dot[0]))/(1 + dot[0] - dot[1] - dot[2]);
      // now loop over jets
      for(unsigned int ix=0;ix<3;++ix) {
	_njet->fill();
	for(const Particle & p : jets[ix].constituents()) {
	  if(p.isCharged()) {
	    _h->fill(p.momentum().p3().mod()/Ejet[ix]);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h,1./ *_njet);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    CounterPtr _njet;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1985_I216850);

}
