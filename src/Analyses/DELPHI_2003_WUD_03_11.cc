// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/DELPHI_2003_WUD_03_11.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  // Constructor
  DELPHI_2003_WUD_03_11::DELPHI_2003_WUD_03_11() 
  {
    const ChargedFinalState cfs;
    addProjection(cfs, "FS");
    addProjection(FastJets(cfs, FastJets::JADE, 0.7), "JadeJets");
    addProjection(FastJets(cfs, FastJets::DURHAM, 0.7), "DurhamJets");
    _numdurjets = 0;
    _numjadejets = 0;
  }


  double DELPHI_2003_WUD_03_11::calc_BZ(std::vector<fastjet::PseudoJet> jets) {
    assert(jets.size()==4);
    Vector3 p12 = cross( momentum3(jets[0]), momentum3(jets[1]));
    Vector3 p34 = cross( momentum3(jets[2]), momentum3(jets[3]));
    return dot(p12,p34) / (p12.mod()*p34.mod());
  }

  double DELPHI_2003_WUD_03_11::calc_KSW(std::vector<fastjet::PseudoJet> jets) {
    assert(jets.size()==4);
    Vector3 p13 = cross( momentum3(jets[0]), momentum3(jets[2]));
    Vector3 p24 = cross( momentum3(jets[1]), momentum3(jets[3]));
    Vector3 p14 = cross( momentum3(jets[0]), momentum3(jets[3]));
    Vector3 p23 = cross( momentum3(jets[1]), momentum3(jets[2]));
    return cos (0.5*( acos (dot(p14,p23) / (p14.mod()*p23.mod())) +
                      acos (dot(p13,p24) / (p13.mod()*p24.mod())) ));
  }

  double DELPHI_2003_WUD_03_11::calc_NR(std::vector<fastjet::PseudoJet> jets) {
    assert(jets.size()==4);
    Vector3 p12 = momentum3(jets[0]) - momentum3(jets[1]);
    Vector3 p34 = momentum3(jets[2]) - momentum3(jets[3]);
    return dot(p12,p34) / (p12.mod()*p34.mod());
  }

  double DELPHI_2003_WUD_03_11::calc_ALPHA34(std::vector<fastjet::PseudoJet> jets) {
    assert(jets.size()==4);
    Vector3 p3 = momentum3(jets[2]);
    Vector3 p4 = momentum3(jets[3]);
    return dot(p3,p4) / (p3.mod()*p4.mod());
  }

  void DELPHI_2003_WUD_03_11::analyze(const Event& e) {
    // First, veto on leptonic events by requiring at least 4 charged FS particles
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const size_t numParticles = fs.particles().size();

    // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
    if (numParticles < 2) {
      getLog() << Log::DEBUG << "Failed multiplicity cut" << endl;
      vetoEvent(e);
    }
    getLog() << Log::DEBUG << "Passed multiplicity cut" << endl;

    // Get event weight for histo filling
    const double weight = e.weight();

    // Jets
    getLog() << Log::DEBUG << "Using FastJet JADE patch to make diff jet rate plots:" << endl;
    const FastJets& durjet = applyProjection<FastJets>(e, "DurhamJets");
    std::vector<fastjet::PseudoJet> jets_durham;
    if (durjet.clusterSeq()) {
      jets_durham = fastjet::sorted_by_E(durjet.clusterSeq()->exclusive_jets(0.008));
      if (jets_durham.size()==4) {
        _histDurhamBZ->fill(fabs(calc_BZ(jets_durham)), weight);
        _histDurhamKSW->fill(calc_KSW(jets_durham), weight);
        _histDurhamNR->fill(fabs(calc_NR(jets_durham)), weight);
        _histDurhamALPHA34->fill(calc_ALPHA34(jets_durham), weight);
      }
      if (durjet.clusterSeq()->exclusive_dmerge(3) > 0.008 && durjet.clusterSeq()->exclusive_dmerge(4) < 0.008)
        _numdurjets++;
    }

    const FastJets& jadejet = applyProjection<FastJets>(e, "JadeJets");
    std::vector<fastjet::PseudoJet> jets_jade;
    if (jadejet.clusterSeq()) {
      jets_jade = fastjet::sorted_by_E(jadejet.clusterSeq()->exclusive_jets(0.015));
      if (jets_jade.size()==4) {
        _histJadeBZ->fill(fabs(calc_BZ(jets_jade)), weight);
        _histJadeKSW->fill(calc_KSW(jets_jade), weight);
        _histJadeNR->fill(fabs(calc_NR(jets_jade)), weight);
        _histJadeALPHA34->fill(calc_ALPHA34(jets_jade), weight);
      }
      if (jadejet.clusterSeq()->exclusive_dmerge(3) > 0.015 && jadejet.clusterSeq()->exclusive_dmerge(4) < 0.015)
        _numjadejets++;
    }

  }



  void DELPHI_2003_WUD_03_11::init() {
    _histDurhamBZ      = bookHistogram1D(1, 1, 1, "Bengtsson-Zerwas $|\\cos(\\chi_\\text{BZ})|$, Durham $y_\\text{cut}=0.008$",
                                                  "$|\\cos(\\chi_\\text{BZ})|$", "$1/\\sigma \\, \\text{d}{\\sigma}/\\text{d}|\\cos(\\chi_\\text{BZ})|$");
    _histDurhamKSW     = bookHistogram1D(2, 1, 1, "K\\\"orner-Schierholz-Willrodt $\\cos(\\phi_\\text{KSW})$, Durham $y_\\text{cut}=0.008$",
                                                  "$\\cos(\\phi_\\text{KSW})$", "$1/\\sigma \\, \\text{d}{\\sigma}/\\text{d}\\,\\cos(\\phi_\\text{KSW})$");
    _histDurhamNR      = bookHistogram1D(3, 1, 1, "Nachtmann-Reiter (mod.) $|\\cos(\\theta^*_\\text{NR})|$, Durham $y_\\text{cut}=0.008$",
                                                  "$|\\cos(\\theta^*_\\text{NR})|$", "$1/\\sigma \\, \\text{d}{\\sigma}/\\text{d}|\\cos(\\theta^*_\\text{NR})|$");
    _histDurhamALPHA34 = bookHistogram1D(4, 1, 1, "$\\cos(\\alpha_{34})$, Durham $y_\\text{cut}=0.008$",
                                                  "$\\cos(\\alpha_{34})$", "$1/\\sigma \\, \\text{d}{\\sigma}/\\text{d}\\,\\cos(\\alpha_{34})$");
    _histJadeBZ        = bookHistogram1D(1, 2, 1, "Bengtsson-Zerwas $|\\cos(\\chi_\\text{BZ})|$, Jade $y_\\text{cut}=0.015$",
                                                  "$|\\cos(\\chi_\\text{BZ})|$", "$1/\\sigma \\, \\text{d}{\\sigma}/\\text{d}|\\cos(\\chi_\\text{BZ})|$");
    _histJadeKSW       = bookHistogram1D(2, 2, 1, "K\\\"orner-Schierholz-Willrodt $\\cos(\\phi_\\text{KSW})$, Jade $y_\\text{cut}=0.015$",
                                                  "$\\cos(\\phi_\\text{KSW})$", "$1/\\sigma \\, \\text{d}{\\sigma}/\\text{d}\\,\\cos(\\phi_\\text{KSW})$");
    _histJadeNR        = bookHistogram1D(3, 2, 1, "Nachtmann-Reiter (mod.) $|\\cos(\\theta^*_\\text{NR})|$, Jade $y_\\text{cut}=0.015$",
                                                  "$|\\cos(\\theta^*_\\text{NR})|$", "$1/\\sigma \\, \\text{d}{\\sigma}/\\text{d}|\\cos(\\theta^*_\\text{NR})|$");
    _histJadeALPHA34   = bookHistogram1D(4, 2, 1, "$\\cos(\\alpha_{34})$, Jade $y_\\text{cut}=0.015$",
                                                  "$\\cos(\\alpha_{34})$", "$1/\\sigma \\, \\text{d}{\\sigma}/\\text{d}\\,\\cos(\\alpha_{34})$");
  }



  // Finalize
  void DELPHI_2003_WUD_03_11::finalize() { 
    // Normalize inclusive single particle distributions to the average number 
    // of charged particles per event.

    getLog() << Log::WARN << "numdurjets   = " << _numdurjets << endl;
    getLog() << Log::WARN << "numjadejets  = " << _numjadejets << endl;
    getLog() << Log::WARN << "sumofweights = " << sumOfWeights() << endl;
    normalize(_histDurhamBZ      , 0.0785);
    normalize(_histDurhamKSW     , 0.0785);
    normalize(_histDurhamNR      , 0.0785);
    normalize(_histDurhamALPHA34 , 0.0785);
    normalize(_histJadeBZ        , 0.0277);
    normalize(_histJadeKSW       , 0.0277);
    normalize(_histJadeNR        , 0.0277);
    normalize(_histJadeALPHA34   , 0.0277);
  }

}
