#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  class ALICE_2011_S8945144 : public Analysis {
  public:

    ALICE_2011_S8945144()
      : Analysis("ALICE_2011_S8945144")
    {}


  public:

    void init() {
      const ChargedFinalState cfs(-15, 15);
      addProjection(cfs, "CFS");

      _histPtPions        = bookHistogram1D("d01-x01-y01");
      _histPtAntiPions    = bookHistogram1D("d01-x01-y02");
      _histPtKaons        = bookHistogram1D("d02-x01-y01");
      _histPtAntiKaons    = bookHistogram1D("d02-x01-y02");
      _histPtProtons      = bookHistogram1D("d03-x01-y01");
      _histPtAntiProtons  = bookHistogram1D("d03-x01-y02");
      _histAveragePt      = bookProfile1D("d04-x01-y01");
    }


    void analyze(const Event& event) {
      const double weight = event.weight();
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");
      foreach (const Particle& p, cfs.particles()) {
        if(fabs(p.momentum().rapidity())<0.5) {
          switch (p.pdgId()) {
            case 211:
              _histPtPions->fill(p.momentum().pT()/GeV, weight);
              _histAveragePt->fill(p.mass()/GeV, p.momentum().pT()/GeV, weight);
              break;
            case -211:
              _histPtAntiPions->fill(p.momentum().pT()/GeV, weight);
              _histAveragePt->fill(p.mass()/GeV, p.momentum().pT()/GeV, weight);
              break;
            case 2212:
              _histPtProtons->fill(p.momentum().pT()/GeV, weight);
              _histAveragePt->fill(p.mass()/GeV, p.momentum().pT()/GeV, weight);
              break;
            case -2212:
              _histPtAntiProtons->fill(p.momentum().pT()/GeV, weight);
              _histAveragePt->fill(p.mass()/GeV, p.momentum().pT()/GeV, weight);
              break;
            case 321:
              _histPtKaons->fill(p.momentum().pT()/GeV, weight);
              _histAveragePt->fill(p.mass()/GeV, p.momentum().pT()/GeV, weight);
              break;
            case -321:
              _histPtAntiKaons->fill(p.momentum().pT()/GeV, weight);
              _histAveragePt->fill(p.mass()/GeV, p.momentum().pT()/GeV, weight);
              break;
          }
        }
      }
    }


    void finalize() {
      scale(_histPtPions,       1./sumOfWeights());
      scale(_histPtProtons,     1./sumOfWeights());
      scale(_histPtKaons,       1./sumOfWeights());
      scale(_histPtAntiPions,   1./sumOfWeights());
      scale(_histPtAntiProtons, 1./sumOfWeights());
      scale(_histPtAntiKaons,   1./sumOfWeights());
    }


  private:

    AIDA::IHistogram1D *_histPtPions;
    AIDA::IHistogram1D *_histPtProtons;
    AIDA::IHistogram1D *_histPtKaons;
    AIDA::IHistogram1D *_histPtAntiPions;
    AIDA::IHistogram1D *_histPtAntiProtons;
    AIDA::IHistogram1D *_histPtAntiKaons;
    AIDA::IProfile1D   *_histAveragePt;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2011_S8945144);

}
