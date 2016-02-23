#include "Rivet/Analysis.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "GeneratorInterface/RivetInterface/interface/PseudoTop.hh"
//#include "TopMonteCarlo/RivetTop/interface/PseudoTop.hh"

namespace Rivet {

class CMS_2015_I1370682 : public Analysis {
public:
  CMS_2015_I1370682() : Analysis("CMS_2015_I1370682"),
    _applyCorrection(false), _doShapeOnly(false) {
  }

  void init() {
    addProjection(PseudoTop(0.1, 20, 2.4, 0.5, 30, 2.4), "ttbar");

    // Lepton + Jet channel
    _hSL_topPt         = bookHisto1D("d15-x01-y01"); // 1/sigma dsigma/dpt(top)
    _hSL_topPtTtbarSys = bookHisto1D("d16-x01-y01"); // 1/sigma dsigma/dpt*(top)
    _hSL_topY          = bookHisto1D("d17-x01-y01"); // 1/sigma dsigma/dy(top)
    _hSL_ttbarDelPhi   = bookHisto1D("d18-x01-y01"); // 1/sigma dsigma/ddeltaphi(t,tbar)
    _hSL_topPtLead     = bookHisto1D("d19-x01-y01"); // 1/sigma dsigma/dpt(t1)
    _hSL_topPtSubLead  = bookHisto1D("d20-x01-y01"); // 1/sigma dsigma/dpt(t2)
    _hSL_ttbarPt       = bookHisto1D("d21-x01-y01"); // 1/sigma dsigma/dpt(ttbar)
    _hSL_ttbarY        = bookHisto1D("d22-x01-y01"); // 1/sigma dsigma/dy(ttbar)
    _hSL_ttbarMass     = bookHisto1D("d23-x01-y01"); // 1/sigma dsigma/dm(ttbar)

    // Dilepton channel
    _hDL_topPt         = bookHisto1D("d24-x01-y01"); // 1/sigma dsigma/dpt(top)
    _hDL_topPtTtbarSys = bookHisto1D("d25-x01-y01"); // 1/sigma dsigma/dpt*(top)
    _hDL_topY          = bookHisto1D("d26-x01-y01"); // 1/sigma dsigma/dy(top)
    _hDL_ttbarDelPhi   = bookHisto1D("d27-x01-y01"); // 1/sigma dsigma/ddeltaphi(t,tbar)
    _hDL_topPtLead     = bookHisto1D("d28-x01-y01"); // 1/sigma dsigma/dpt(t1)
    _hDL_topPtSubLead  = bookHisto1D("d29-x01-y01"); // 1/sigma dsigma/dpt(t2)
    _hDL_ttbarPt       = bookHisto1D("d30-x01-y01"); // 1/sigma dsigma/dpt(ttbar)
    _hDL_ttbarY        = bookHisto1D("d31-x01-y01"); // 1/sigma dsigma/dy(ttbar)
    _hDL_ttbarMass     = bookHisto1D("d32-x01-y01"); // 1/sigma dsigma/dm(ttbar)

  };

  void analyze(const Event& event) {
    const double weight = event.weight();

    // Get the ttbar candidate
    const PseudoTop& ttbar = applyProjection<PseudoTop>(event, "ttbar");
    if ( ttbar.mode() == PseudoTop::CH_NONE ) vetoEvent;

    const FourMomentum& t1P4 = ttbar.t1().momentum();
    const FourMomentum& t2P4 = ttbar.t2().momentum();
    const double pt1 = std::max(t1P4.pT(), t2P4.pT());
    const double pt2 = std::min(t1P4.pT(), t2P4.pT());
    const double dPhi = deltaPhi(t1P4, t2P4);
    const FourMomentum ttP4 = t1P4+t2P4;
    const FourMomentum t1P4AtCM = LorentzTransform(-ttP4.boostVector()).transform(t1P4);

    if ( ttbar.mode() == PseudoTop::CH_SEMILEPTON ) {
      const Particle lCand1 = ttbar.wDecays1()[0]; // w1 dau0 is the lepton in the PseudoTop
      if ( lCand1.pt() < 33 or std::abs(lCand1.eta()) > 2.1 ) vetoEvent;

      _hSL_topPt->fill(t1P4.pT(), weight);
      _hSL_topPt->fill(t2P4.pT(), weight);
      _hSL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
      _hSL_topY->fill(t1P4.rapidity(), weight);
      _hSL_topY->fill(t2P4.rapidity(), weight);
      _hSL_ttbarDelPhi->fill(dPhi, weight);
      _hSL_topPtLead->fill(pt1, weight);
      _hSL_topPtSubLead->fill(pt2, weight);
      _hSL_ttbarPt->fill(ttP4.pT(), weight);
      _hSL_ttbarY->fill(ttP4.rapidity(), weight);
      _hSL_ttbarMass->fill(ttP4.mass(), weight);
    }
    else if ( ttbar.mode() == PseudoTop::CH_FULLLEPTON ) {
      const Particle lCand1 = ttbar.wDecays1()[0]; // dau0 are the lepton in the PseudoTop
      const Particle lCand2 = ttbar.wDecays2()[0]; // dau0 are the lepton in the PseudoTop
      if ( lCand1.pt() < 20 or std::abs(lCand1.eta()) > 2.4 or
           lCand2.pt() < 20 or std::abs(lCand2.eta()) > 2.4 ) vetoEvent;

      _hDL_topPt->fill(t1P4.pT(), weight);
      _hDL_topPt->fill(t2P4.pT(), weight);
      _hDL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
      _hDL_topY->fill(t1P4.rapidity(), weight);
      _hDL_topY->fill(t2P4.rapidity(), weight);
      _hDL_ttbarDelPhi->fill(dPhi, weight);
      _hDL_topPtLead->fill(pt1, weight);
      _hDL_topPtSubLead->fill(pt2, weight);
      _hDL_ttbarPt->fill(ttP4.pT(), weight);
      _hDL_ttbarY->fill(ttP4.rapidity(), weight);
      _hDL_ttbarMass->fill(ttP4.mass(), weight);
    }

  };

  void finalize() {
    if ( _applyCorrection ) {
      // Correction functions for TOP-12-028 paper, (parton bin height)/(pseudotop bin height)
      const double ch15[] = { 5.473609, 4.941048, 4.173346, 3.391191, 2.785644, 2.371346, 2.194161, 2.197167, };
      const double ch16[] = { 5.470905, 4.948201, 4.081982, 3.225532, 2.617519, 2.239217, 2.127878, 2.185918, };
      const double ch17[] = { 10.003667, 4.546519, 3.828115, 3.601018, 3.522194, 3.524694, 3.600951, 3.808553, 4.531891, 9.995370, };
      const double ch18[] = { 4.406683, 4.054041, 3.885393, 4.213646, };
      const double ch19[] = { 6.182537, 5.257703, 4.422280, 3.568402, 2.889408, 2.415878, 2.189974, 2.173210, };
      const double ch20[] = { 5.199874, 4.693318, 3.902882, 3.143785, 2.607877, 2.280189, 2.204124, 2.260829, };
      const double ch21[] = { 6.053523, 3.777506, 3.562251, 3.601356, 3.569347, 3.410472, };
      const double ch22[] = { 11.932351, 4.803773, 3.782709, 3.390775, 3.226806, 3.218982, 3.382678, 3.773653, 4.788191, 11.905338, };
      const double ch23[] = { 7.145255, 5.637595, 4.049882, 3.025917, 2.326430, 1.773824, 1.235329, };

      const double ch24[] = { 2.268193, 2.372063, 2.323975, 2.034655, 1.736793, };
      const double ch25[] = { 2.231852, 2.383086, 2.341894, 2.031318, 1.729672, 1.486993, };
      const double ch26[] = { 3.993526, 2.308249, 2.075136, 2.038297, 2.036302, 2.078270, 2.295817, 4.017713, };
      const double ch27[] = { 2.205978, 2.175010, 2.215376, 2.473144, };
      const double ch28[] = { 2.321077, 2.371895, 2.338871, 2.057821, 1.755382, };
      const double ch29[] = { 2.222707, 2.372591, 2.301688, 1.991162, 1.695343, };
      const double ch30[] = { 2.599677, 2.026855, 2.138620, 2.229553, };
      const double ch31[] = { 5.791779, 2.636219, 2.103642, 1.967198, 1.962168, 2.096514, 2.641189, 5.780828, };
      const double ch32[] = { 2.006685, 2.545525, 2.477745, 2.335747, 2.194226, 2.076500, };

      applyCorrection(_hSL_topPt, ch15);
      applyCorrection(_hSL_topPtTtbarSys, ch16);
      applyCorrection(_hSL_topY, ch17);
      applyCorrection(_hSL_ttbarDelPhi, ch18);
      applyCorrection(_hSL_topPtLead, ch19);
      applyCorrection(_hSL_topPtSubLead, ch20);
      applyCorrection(_hSL_ttbarPt, ch21);
      applyCorrection(_hSL_ttbarY, ch22);
      applyCorrection(_hSL_ttbarMass, ch23);

      applyCorrection(_hDL_topPt, ch24);
      applyCorrection(_hDL_topPtTtbarSys, ch25);
      applyCorrection(_hDL_topY, ch26);
      applyCorrection(_hDL_ttbarDelPhi, ch27);
      applyCorrection(_hDL_topPtLead, ch28);
      applyCorrection(_hDL_topPtSubLead, ch29);
      applyCorrection(_hDL_ttbarPt, ch30);
      applyCorrection(_hDL_ttbarY, ch31);
      applyCorrection(_hDL_ttbarMass, ch32);
    }

    if ( _doShapeOnly ) {
      normalize(_hSL_topPt        );
      normalize(_hSL_topPtTtbarSys);
      normalize(_hSL_topY         );
      normalize(_hSL_ttbarDelPhi  );
      normalize(_hSL_topPtLead    );
      normalize(_hSL_topPtSubLead );
      normalize(_hSL_ttbarPt      );
      normalize(_hSL_ttbarY       );
      normalize(_hSL_ttbarMass    );

      normalize(_hDL_topPt        );
      normalize(_hDL_topPtTtbarSys);
      normalize(_hDL_topY         );
      normalize(_hDL_ttbarDelPhi  );
      normalize(_hDL_topPtLead    );
      normalize(_hDL_topPtSubLead );
      normalize(_hDL_ttbarPt      );
      normalize(_hDL_ttbarY       );
      normalize(_hDL_ttbarMass    );
    }
    else {
      const double s = 1./sumOfWeights();
      scale(_hSL_topPt        , s);
      scale(_hSL_topPtTtbarSys, s);
      scale(_hSL_topY         , s);
      scale(_hSL_ttbarDelPhi  , s);
      scale(_hSL_topPtLead    , s);
      scale(_hSL_topPtSubLead , s);
      scale(_hSL_ttbarPt      , s);
      scale(_hSL_ttbarY       , s);
      scale(_hSL_ttbarMass    , s);

      scale(_hDL_topPt        , s);
      scale(_hDL_topPtTtbarSys, s);
      scale(_hDL_topY         , s);
      scale(_hDL_ttbarDelPhi  , s);
      scale(_hDL_topPtLead    , s);
      scale(_hDL_topPtSubLead , s);
      scale(_hDL_ttbarPt      , s);
      scale(_hDL_ttbarY       , s);
      scale(_hDL_ttbarMass    , s);
    }

  };

  void applyCorrection(Histo1DPtr h, const double* cf) {
    std::vector<YODA::HistoBin1D>& bins = h->bins();
    for ( int i=0, n=bins.size(); i<n; ++i ) {
      const double s = cf[i];
      YODA::HistoBin1D& bin = bins[i];
      bin.scaleW(s);
    }
  };

private:
  const bool _applyCorrection, _doShapeOnly;

  Histo1DPtr _hSL_topPt        ;
  Histo1DPtr _hSL_topPtTtbarSys;
  Histo1DPtr _hSL_topY         ;
  Histo1DPtr _hSL_ttbarDelPhi  ;
  Histo1DPtr _hSL_topPtLead    ;
  Histo1DPtr _hSL_topPtSubLead ;
  Histo1DPtr _hSL_ttbarPt      ;
  Histo1DPtr _hSL_ttbarY       ;
  Histo1DPtr _hSL_ttbarMass    ;

  Histo1DPtr _hDL_topPt        ;
  Histo1DPtr _hDL_topPtTtbarSys;
  Histo1DPtr _hDL_topY         ;
  Histo1DPtr _hDL_ttbarDelPhi  ;
  Histo1DPtr _hDL_topPtLead    ;
  Histo1DPtr _hDL_topPtSubLead ;
  Histo1DPtr _hDL_ttbarPt      ;
  Histo1DPtr _hDL_ttbarY       ;
  Histo1DPtr _hDL_ttbarMass    ;

};

// This global object acts as a hook for the plugin system
AnalysisBuilder<CMS_2015_I1370682> plugin_CMS_2015_I1370682;

}
