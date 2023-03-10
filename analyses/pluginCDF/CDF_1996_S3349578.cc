// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/SmearedJets.hh"

namespace Rivet {


  /// @brief CDF properties of high-mass multi-jet events
  class CDF_1996_S3349578 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_1996_S3349578);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections here
      const FinalState fs(Cuts::abseta < 4.2);
      FastJets fj(fs, FastJets::CDFJETCLU, 0.7);
      declare(fj, "Jets");

      // Smear energy and mass with the 10% uncertainty quoted in the paper
      SmearedJets sj_E(fj, [](const Jet& jet){ return P4_SMEAR_MASS_GAUSS(P4_SMEAR_E_GAUSS(jet, 0.1*jet.E()), 0.1*jet.mass()); });
      declare(sj_E, "SmearedJets");

      /// Book histograms here, e.g.:
      book(_h_3_mNJ ,1, 1, 1);
      book(_h_3_X3 ,2, 1, 1);
      book(_h_3_X4 ,3, 1, 1);
      book(_h_3_costheta3 ,8, 1, 1);
      book(_h_3_psi3 ,9, 1, 1);
      book(_h_3_f3 ,14, 1, 1);
      book(_h_3_f4 ,14, 1, 2);
      book(_h_3_f5 ,14, 1, 3);

      book(_h_4_mNJ ,1, 1, 2);
      book(_h_4_X3 ,4, 1, 1);
      book(_h_4_X4 ,5, 1, 1);
      book(_h_4_costheta3 ,10, 1, 1);
      book(_h_4_psi3 ,11, 1, 1);
      book(_h_4_f3 ,15, 1, 1);
      book(_h_4_f4 ,15, 1, 2);
      book(_h_4_f5 ,15, 1, 3);
      book(_h_4_XA ,17, 1, 1);
      book(_h_4_psiAB ,19, 1, 1);
      book(_h_4_fA ,21, 1, 1);
      book(_h_4_fB ,21, 1, 2);

      book(_h_5_mNJ ,1, 1, 3);
      book(_h_5_X3 ,6, 1, 1);
      book(_h_5_X4 ,7, 1, 1);
      book(_h_5_costheta3 ,12, 1, 1);
      book(_h_5_psi3 ,13, 1, 1);
      book(_h_5_f3 ,16, 1, 1);
      book(_h_5_f4 ,16, 1, 2);
      book(_h_5_f5 ,16, 1, 3);
      book(_h_5_XA ,18, 1, 1);
      book(_h_5_XC ,18, 1, 2);
      book(_h_5_psiAB ,20, 1, 1);
      book(_h_5_psiCD ,20, 1, 2);
      book(_h_5_fA ,22, 1, 1);
      book(_h_5_fB ,23, 1, 1);
      book(_h_5_fC ,24, 1, 1);
      book(_h_5_fD ,25, 1, 1);
    }


    void analyze(const Event& event) {
      Jets jets;
      FourMomentum jetsystem(0.0, 0.0, 0.0, 0.0);
      for (const Jet& jet : apply<JetAlg>(event, "SmearedJets").jets(Cuts::Et > 20.0*GeV, cmpMomByEt)) {
        bool separated = true;
        for (const Jet& ref : jets) {
          if (deltaR(jet, ref) < 0.9) {
            separated = false;
            break;
          }
        }
        if (!separated) continue;
        jets.push_back(jet);
        jetsystem += jet.momentum();
        if (jets.size() >= 5) break;
      }

      if (jets.size() > 4) {
        _fiveJetAnalysis(jets);
        jets.resize(4);
      }
      if (jets.size() > 3) {
        _fourJetAnalysis(jets);
        jets.resize(3);
      }
      if (jets.size() > 2) {
        _threeJetAnalysis(jets);
      }
    }


    void _threeJetAnalysis(const Jets& jets) {
      MSG_DEBUG("3 jet analysis");

      double sumEt = 0.0;
      FourMomentum jetsystem(0.0, 0.0, 0.0, 0.0);
      for (const Jet& jet : jets) {
        sumEt += jet.Et();
        jetsystem += jet.momentum();
      }
      if (sumEt < 420.0*GeV) return;

      const double m3J = _safeMass(jetsystem);
      if (m3J < 600*GeV) return;

      const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(jetsystem.betaVec());
      vector<FourMomentum> jets3;
      for (Jet jet : jets) {
        jets3.push_back(cms_boost.transform(jet.momentum()));
      }
      std::sort(jets3.begin(), jets3.end(), FourMomentum::byEDescending());
      FourMomentum p3(jets3[0]), p4(jets3[1]), p5(jets3[2]);

      FourMomentum pAV = cms_boost.transform(_avg_beam_in_lab(m3J, jetsystem.rapidity()));
      double costheta3 = pAV.p3().unit().dot(p3.p3().unit());
      if (fabs(costheta3) > 0.6) return;

      double X3 = 2.0*p3.E()/m3J;
      if (X3 > 0.9) return;

      const double X4 = 2.0*p4.E()/m3J;
      const double psi3 = _psi(p3, pAV, p4, p5);
      const double f3 = _safeMass(p3)/m3J;
      const double f4 = _safeMass(p4)/m3J;
      const double f5 = _safeMass(p5)/m3J;

      _h_3_mNJ->fill(m3J);
      _h_3_X3->fill(X3);
      _h_3_X4->fill(X4);
      _h_3_costheta3->fill(costheta3);
      _h_3_psi3->fill(psi3);
      _h_3_f3->fill(f3);
      _h_3_f4->fill(f4);
      _h_3_f5->fill(f5);

    }


    void _fourJetAnalysis(const Jets& jets) {
      MSG_DEBUG("4 jet analysis");

      double sumEt=0.0;
      FourMomentum jetsystem(0.0, 0.0, 0.0, 0.0);
      for (const Jet& jet : jets) {
        sumEt+=jet.Et();
        jetsystem+=jet.momentum();
      }
      if (sumEt < 420.0*GeV) return;

      const double m4J = _safeMass(jetsystem);
      if (m4J < 650*GeV) return;

      const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(jetsystem.betaVec());
      vector<FourMomentum> jets4;
      for (Jet jet : jets) {
        jets4.push_back(cms_boost.transform(jet.momentum()));
      }
      std::sort(jets4.begin(), jets4.end(), FourMomentum::byEDescending());

      FourMomentum pA, pB;
      vector<FourMomentum> jets3(_reduce(jets4, pA, pB));
      std::sort(jets3.begin(), jets3.end(), FourMomentum::byEDescending());
      FourMomentum p3(jets3[0]);
      FourMomentum p4(jets3[1]);
      FourMomentum p5(jets3[2]);

      FourMomentum pAV = cms_boost.transform(_avg_beam_in_lab(m4J, jetsystem.rapidity()));
      double costheta3=pAV.p3().unit().dot(p3.p3().unit());
      if (fabs(costheta3)>0.8) {
        return;
      }

      const double X3 = 2.0*p3.E()/m4J;
      if (X3>0.9) {
        return;
      }

      // fill histograms
      const double X4 = 2.0*p4.E()/m4J;
      const double psi3 = _psi(p3, pAV, p4, p5);
      const double f3 = _safeMass(p3)/m4J;
      const double f4 = _safeMass(p4)/m4J;
      const double f5 = _safeMass(p5)/m4J;
      const double fA = _safeMass(pA)/m4J;
      const double fB = _safeMass(pB)/m4J;
      const double XA = pA.E()/(pA.E()+pB.E());
      const double psiAB = _psi(pA, pB, pA+pB, pAV);

      _h_4_mNJ->fill(m4J);
      _h_4_X3->fill(X3);
      _h_4_X4->fill(X4);
      _h_4_costheta3->fill(costheta3);
      _h_4_psi3->fill(psi3);
      _h_4_f3->fill(f3);
      _h_4_f4->fill(f4);
      _h_4_f5->fill(f5);
      _h_4_XA->fill(XA);
      _h_4_psiAB->fill(psiAB);
      _h_4_fA->fill(fA);
      _h_4_fB->fill(fB);
    }


    void _fiveJetAnalysis(const Jets& jets) {
      MSG_DEBUG("5 jet analysis");

      double sumEt=0.0;
      FourMomentum jetsystem(0.0, 0.0, 0.0, 0.0);
      for (const Jet& jet : jets) {
        sumEt+=jet.Et();
        jetsystem+=jet.momentum();
      }
      if (sumEt < 420.0*GeV) return;

      const double m5J = _safeMass(jetsystem);
      if (m5J < 750*GeV) return;

      const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(jetsystem.betaVec());
      vector<FourMomentum> jets5;
      for (Jet jet : jets) {
        jets5.push_back(cms_boost.transform(jet.momentum()));
      }
      std::sort(jets5.begin(), jets5.end(), FourMomentum::byEDescending());

      FourMomentum pC, pD;
      vector<FourMomentum> jets4(_reduce(jets5, pC, pD));
      std::sort(jets4.begin(), jets4.end(), FourMomentum::byEDescending());

      FourMomentum pA, pB;
      vector<FourMomentum> jets3(_reduce(jets4, pA, pB));
      std::sort(jets3.begin(), jets3.end(), FourMomentum::byEDescending());
      FourMomentum p3(jets3[0]);
      FourMomentum p4(jets3[1]);
      FourMomentum p5(jets3[2]);

      // fill histograms
      FourMomentum pAV = cms_boost.transform(_avg_beam_in_lab(m5J, jetsystem.rapidity()));
      const double costheta3 = pAV.p3().unit().dot(p3.p3().unit());
      const double X3 = 2.0*p3.E()/m5J;
      const double X4 = 2.0*p4.E()/m5J;
      const double psi3 = _psi(p3, pAV, p4, p5);
      const double f3 = _safeMass(p3)/m5J;
      const double f4 = _safeMass(p4)/m5J;
      const double f5 = _safeMass(p5)/m5J;
      const double fA = _safeMass(pA)/m5J;
      const double fB = _safeMass(pB)/m5J;
      const double XA = pA.E()/(pA.E()+pB.E());
      const double psiAB = _psi(pA, pB, pA+pB, pAV);
      const double fC = _safeMass(pC)/m5J;
      const double fD = _safeMass(pD)/m5J;
      const double XC = pC.E()/(pC.E()+pD.E());
      const double psiCD = _psi(pC, pD, pC+pD, pAV);

      _h_5_mNJ->fill(m5J);
      _h_5_X3->fill(X3);
      _h_5_X4->fill(X4);
      _h_5_costheta3->fill(costheta3);
      _h_5_psi3->fill(psi3);
      _h_5_f3->fill(f3);
      _h_5_f4->fill(f4);
      _h_5_f5->fill(f5);
      _h_5_XA->fill(XA);
      _h_5_psiAB->fill(psiAB);
      _h_5_fA->fill(fA);
      _h_5_fB->fill(fB);
      _h_5_XC->fill(XC);
      _h_5_psiCD->fill(psiCD);
      _h_5_fC->fill(fC);
      _h_5_fD->fill(fD);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      /// Normalise, scale and otherwise manipulate histograms here
      normalize(_h_3_mNJ, 1.0);
      normalize(_h_3_X3, 1.0);
      normalize(_h_3_X4, 1.0);
      normalize(_h_3_costheta3, 1.0);
      normalize(_h_3_psi3, 1.0);
      normalize(_h_3_f3, 1.0);
      normalize(_h_3_f4, 1.0);
      normalize(_h_3_f5, 1.0);

      normalize(_h_4_mNJ, 1.0);
      normalize(_h_4_X3, 1.0);
      normalize(_h_4_X4, 1.0);
      normalize(_h_4_costheta3, 1.0);
      normalize(_h_4_psi3, 1.0);
      normalize(_h_4_f3, 1.0);
      normalize(_h_4_f4, 1.0);
      normalize(_h_4_f5, 1.0);
      normalize(_h_4_XA, 1.0);
      normalize(_h_4_psiAB, 1.0);
      normalize(_h_4_fA, 1.0);
      normalize(_h_4_fB, 1.0);

      normalize(_h_5_mNJ, 1.0);
      normalize(_h_5_X3, 1.0);
      normalize(_h_5_X4, 1.0);
      normalize(_h_5_costheta3, 1.0);
      normalize(_h_5_psi3, 1.0);
      normalize(_h_5_f3, 1.0);
      normalize(_h_5_f4, 1.0);
      normalize(_h_5_f5, 1.0);
      normalize(_h_5_XA, 1.0);
      normalize(_h_5_XC, 1.0);
      normalize(_h_5_psiAB, 1.0);
      normalize(_h_5_psiCD, 1.0);
      normalize(_h_5_fA, 1.0);
      normalize(_h_5_fB, 1.0);
      normalize(_h_5_fC, 1.0);
      normalize(_h_5_fD, 1.0);

    }

    //@}


  private:

    vector<FourMomentum> _reduce(const vector<FourMomentum>& jets, FourMomentum& combined1, FourMomentum& combined2) {
      double minMass2 = 1e9;
      size_t idx1(jets.size()), idx2(jets.size());
      for (size_t i=0; i<jets.size(); ++i) {
        for (size_t j=i+1; j<jets.size(); ++j) {
          double mass2 = FourMomentum(jets[i]+jets[j]).mass2();
          if (mass2<minMass2) {
            idx1=i;
            idx2=j;
          }
        }
      }
      vector<FourMomentum> newjets;
      for (size_t i=0; i<jets.size(); ++i) {
        if (i!=idx1 && i!=idx2) newjets.push_back(jets[i]);
      }
      newjets.push_back(jets[idx1]+jets[idx2]);
      combined1 = jets[idx1];
      combined2 = jets[idx2];
      return newjets;
    }


    FourMomentum _avg_beam_in_lab(const double& m, const double& y) {
      const double mt = m/2.0;
      FourMomentum beam1(mt, 0, 0, mt);
      FourMomentum beam2(mt, 0, 0, -mt);
      if (fabs(y)>1e-3) {
        FourMomentum boostvec(cosh(y), 0.0, 0.0, sinh(y));
        const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(boostvec.betaVec()).inverse();
        beam1 = cms_boost.transform(beam1);
        beam2 = cms_boost.transform(beam2);
      }
      return (beam1.E() > beam2.E()) ? beam1-beam2 : beam2-beam1;
    }


    double _psi(const FourMomentum& p1, const FourMomentum& p2,
                const FourMomentum& p3, const FourMomentum& p4) {
      Vector3 p1xp2 = p1.p3().cross(p2.p3());
      Vector3 p3xp4 = p3.p3().cross(p4.p3());
      return mapAngle0ToPi(acos(p1xp2.unit().dot(p3xp4.unit())));
    }


    double _safeMass(const FourMomentum& p) {
      double mass2=p.mass2();
      if (mass2>0.0) return sqrt(mass2);
      else if (mass2<-1.0e-5) {
        MSG_WARNING("m2 = " << m2 << ". Assuming m2=0.");
        return 0.0;
      }
      else return 0.0;
    }


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_3_mNJ;
    Histo1DPtr _h_3_X3;
    Histo1DPtr _h_3_X4;
    Histo1DPtr _h_3_costheta3;
    Histo1DPtr _h_3_psi3;
    Histo1DPtr _h_3_f3;
    Histo1DPtr _h_3_f4;
    Histo1DPtr _h_3_f5;

    Histo1DPtr _h_4_mNJ;
    Histo1DPtr _h_4_X3;
    Histo1DPtr _h_4_X4;
    Histo1DPtr _h_4_costheta3;
    Histo1DPtr _h_4_psi3;
    Histo1DPtr _h_4_f3;
    Histo1DPtr _h_4_f4;
    Histo1DPtr _h_4_f5;
    Histo1DPtr _h_4_XA;
    Histo1DPtr _h_4_psiAB;
    Histo1DPtr _h_4_fA;
    Histo1DPtr _h_4_fB;

    Histo1DPtr _h_5_mNJ;
    Histo1DPtr _h_5_X3;
    Histo1DPtr _h_5_X4;
    Histo1DPtr _h_5_costheta3;
    Histo1DPtr _h_5_psi3;
    Histo1DPtr _h_5_f3;
    Histo1DPtr _h_5_f4;
    Histo1DPtr _h_5_f5;
    Histo1DPtr _h_5_XA;
    Histo1DPtr _h_5_XC;
    Histo1DPtr _h_5_psiAB;
    Histo1DPtr _h_5_psiCD;
    Histo1DPtr _h_5_fA;
    Histo1DPtr _h_5_fB;
    Histo1DPtr _h_5_fC;
    Histo1DPtr _h_5_fD;
    //@}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_1996_S3349578, CDF_1996_I418504);

}
