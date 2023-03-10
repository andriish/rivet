// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {

  /// ATLAS same-sign WW production at 8 TeV (PRL)
  class ATLAS_2014_I1298023 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2014_I1298023);


    /// @name Analysis methods
    //@{

    void init() {

      const FinalState fs;

      // bare leptons
      ChargedLeptons bare_leptons(fs);

      // dressed leptons
      Cut cuts = (Cuts::abseta < 2.5) & (Cuts::pT > 25*GeV);
      DressedLeptons leptons(fs, bare_leptons, 0.1, cuts);
      declare(leptons, "leptons");

      // MET
      declare(MissingMomentum(fs), "MissingET");

      // jets
      VetoedFinalState vfs(fs);
      vfs.addVetoPairId(PID::MUON);
      vfs.vetoNeutrinos();
      declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "jets");

      // book histogram
      book(_hist ,1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const vector<DressedLepton>& leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();
      if ( leptons.size() < 2 )  vetoEvent;

      double minDR_ll = DBL_MAX, mll = -1.0;
      for (unsigned int i = 0; i < leptons.size(); ++i) {
        DressedLepton lep1 = leptons[i];
        for (unsigned int j = i + 1; j < leptons.size(); ++j) {
          DressedLepton lep2 = leptons[j];
          double dr = deltaR(lep1, lep2);
          if ( dr < minDR_ll )  minDR_ll = dr;
          double m = (lep1.momentum() + lep2.momentum()).mass();
          if ( mll < 0. || m < mll )  mll = m;
        }
      }
      if ( minDR_ll <= 0.3 || mll <= 20*GeV )  vetoEvent;

      if ( leptons[0].charge() * leptons[1].charge() < 0.0 )  vetoEvent;

      const MissingMomentum& met = apply<MissingMomentum>(event, "MissingET");
      if ( met.vectorEt().mod() <= 40*GeV )  vetoEvent;

      const Jets& all_jets = apply<FastJets>(event, "jets").jetsByPt( (Cuts::abseta < 4.5) && (Cuts::pT > 30*GeV) );
      Jets jets;
      double minDR_overall = DBL_MAX;
      for (const Jet& jet : all_jets) {
        double minDR_jet = DBL_MAX, minDR_electrons = DBL_MAX;
        for( DressedLepton lep : leptons ) {
          double dr = deltaR(jet, lep);
          if ( dr < minDR_jet )  minDR_jet = dr;
          if ( lep.abspid() == 11 && dr < minDR_electrons )  minDR_electrons = dr;
        }
        if ( minDR_electrons < 0.05 )  continue; // veto jet if it overlaps with electron
        if ( minDR_jet < minDR_overall )  minDR_overall = minDR_jet;
        jets += jet;
      }
      if ( jets.size() < 2 || minDR_overall <= 0.3 )  vetoEvent;

      FourMomentum dijet = jets[0].momentum() + jets[1].momentum();
      if ( dijet.mass() <= 500*GeV )  vetoEvent;

      // inclusive region
      _hist->fill(1);

      // VBS region
      if ( deltaRap(jets[0], jets[1]) > 2.4 )  _hist->fill(2);
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      const double scalefactor( crossSection() / sumOfWeights() );
      scale(_hist, scalefactor);

    }
    //@}

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _hist;
    //@}

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2014_I1298023);

}
