// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class CDF_1996_S3349578 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_1996_S3349578()
      : Analysis("CDF_1996_S3349578") 
    {
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);
    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections here
      const FinalState fs(-4.2, 4.2);
      addProjection(FastJets(fs, FastJets::CDFJETCLU, 0.7), "Jets");

      /// Book histograms here, e.g.:
      _h_3_mNJ = bookHistogram1D(1, 1, 1);
      _h_3_X3 = bookHistogram1D(2, 1, 1);
      _h_3_X4 = bookHistogram1D(3, 1, 1);
      _h_3_costheta3 = bookHistogram1D(8, 1, 1);
      _h_3_psi3 = bookHistogram1D(9, 1, 1);
      _h_3_f3 = bookHistogram1D(14, 1, 1);
      _h_3_f4 = bookHistogram1D(14, 1, 2);
      _h_3_f5 = bookHistogram1D(14, 1, 3);
      
      _h_4_mNJ = bookHistogram1D(1, 1, 2);
      _h_4_X3 = bookHistogram1D(4, 1, 1);
      _h_4_X4 = bookHistogram1D(5, 1, 1);
      _h_4_costheta3 = bookHistogram1D(10, 1, 1);
      _h_4_psi3 = bookHistogram1D(11, 1, 1);
      _h_4_f3 = bookHistogram1D(15, 1, 1);
      _h_4_f4 = bookHistogram1D(15, 1, 2);
      _h_4_f5 = bookHistogram1D(15, 1, 3);
      _h_4_XA = bookHistogram1D(17, 1, 1);
      _h_4_psiAB = bookHistogram1D(19, 1, 1);
      _h_4_fA = bookHistogram1D(21, 1, 1);
      _h_4_fB = bookHistogram1D(21, 1, 2);
      
      _h_5_mNJ = bookHistogram1D(1, 1, 3);
      _h_5_X3 = bookHistogram1D(6, 1, 1);
      _h_5_X4 = bookHistogram1D(7, 1, 1);
      _h_5_costheta3 = bookHistogram1D(12, 1, 1);
      _h_5_psi3 = bookHistogram1D(13, 1, 1);
      _h_5_f3 = bookHistogram1D(16, 1, 1);
      _h_5_f4 = bookHistogram1D(16, 1, 2);
      _h_5_f5 = bookHistogram1D(16, 1, 3);
      _h_5_XA = bookHistogram1D(18, 1, 1);
      _h_5_XC = bookHistogram1D(18, 1, 2);
      _h_5_psiAB = bookHistogram1D(20, 1, 1);
      _h_5_psiCD = bookHistogram1D(20, 1, 2);
      _h_5_fA = bookHistogram1D(22, 1, 1);
      _h_5_fB = bookHistogram1D(23, 1, 1);
      _h_5_fC = bookHistogram1D(24, 1, 1);
      _h_5_fD = bookHistogram1D(25, 1, 1);
      
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      /// Do the event by event analysis here
      Jets jets;
      double sumEt = 0.0;
      FourMomentum jetsystem(0.0, 0.0, 0.0, 0.0);
      foreach (const Jet& jet, applyProjection<FastJets>(event, "Jets").jetsByEt()) {
        double Et = jet.momentum().Et();
        if (Et > 20.0*GeV) {
          jets.push_back(jet);
          sumEt += Et;
          jetsystem += jet.momentum();
        }
      }
      /// @todo include gaussian jet energy resolution smearing?
      
      if (jets.size() < 3) {
        vetoEvent;
      }
      
      if (sumEt < 420.0*GeV) {
        vetoEvent;
      }
      
      if (jets.size() > 2) _threeJetAnalysis(jets, weight);
      if (jets.size() > 3) _fourJetAnalysis(jets, weight);
      if (jets.size() > 4) _fiveJetAnalysis(jets, weight);
    }
    
    void _threeJetAnalysis(const Jets& jets, const double& weight) {
      getLog() << Log::DEBUG << "3 jet analysis" << std::endl;
      FourMomentum jjj(jets[0].momentum()+jets[1].momentum()+jets[2].momentum());
      const double m3J = jjj.mass();
      if (m3J<600*GeV) {
        return;
      }
    
      LorentzTransform cms_boost(-jjj.boostVector());
      vector<FourMomentum> jets_boosted;
      foreach (Jet jet, jets) {
        jets_boosted.push_back(cms_boost.transform(jet.momentum()));
      }
      std::sort(jets_boosted.begin(), jets_boosted.end(), FourMomentum::byEDescending());
      FourMomentum p3(jets_boosted[0]);
      FourMomentum p4(jets_boosted[1]);
      FourMomentum p5(jets_boosted[2]);
      
      double costheta3 = cos(p3.theta());
      if (fabs(costheta3)>0.6) {
        return;
      }
      
      double X3 = 2.0*p3.E()/m3J;
      if (X3>0.9) {
        return;
      }
      
      
      // fill histograms
      const double X4 = 2.0*p4.E()/m3J;
      Vector3 beam1(0.0, 0.0, 1.0);
      Vector3 p1xp3 = beam1.cross(p3.vector3());
      Vector3 p4xp5 = p4.vector3().cross(p5.vector3());
      const double cospsi3 = p1xp3.dot(p4xp5)/p1xp3.mod()/p4xp5.mod();
      const double f3 = p3.mass()/m3J;
      const double f4 = p4.mass()/m3J;
      const double f5 = p5.mass()/m3J;
      
      _h_3_mNJ->fill(m3J, weight);
      _h_3_X3->fill(X3, weight);
      _h_3_X4->fill(X4, weight);
      _h_3_costheta3->fill(costheta3, weight);
      _h_3_psi3->fill(mapAngle0ToPi(acos(cospsi3)), weight);
      _h_3_f3->fill(f3, weight);
      _h_3_f4->fill(f4, weight);
      _h_3_f5->fill(f5, weight);
      
    }

    void _fourJetAnalysis(const Jets& jets, const double& weight) {
      getLog() << Log::DEBUG << "4 jet analysis" << std::endl;
      FourMomentum jjjj(0.0, 0.0, 0.0, 0.0);
      vector<FourMomentum> jetmoms;
      for (size_t i=0; i<4; ++i) {
        jetmoms.push_back(jets[i].momentum());
        jjjj += jets[i].momentum();
      }
      const double m4J = jjjj.mass();
      if (m4J < 650*GeV) return;
      
      FourMomentum pA, pB;
      vector<FourMomentum> jetmoms3(_reduce(jetmoms, pA, pB));
      LorentzTransform cms_boost(-jjjj.boostVector());
      vector<FourMomentum> jetmoms3_boosted;
      foreach (FourMomentum mom, jetmoms3) {
        jetmoms3_boosted.push_back(cms_boost.transform(mom));
      }
      pA = cms_boost.transform(pA);
      pB = cms_boost.transform(pB);
      
      sort(jetmoms3_boosted.begin(), jetmoms3_boosted.end(), FourMomentum::byEDescending());
      if (pB.E()>pA.E()) std::swap(pA, pB);
      FourMomentum p3(jetmoms3_boosted[0]);
      FourMomentum p4(jetmoms3_boosted[1]);
      FourMomentum p5(jetmoms3_boosted[2]);
      
      const double costheta3 = cos(p3.theta());
      if (fabs(costheta3)>0.8) {
        return;
      }
      
      const double X3 = 2.0*p3.E()/m4J;
      if (X3>0.9) {
        return;
      }
      
      // fill histograms
      const double X4 = 2.0*p4.E()/m4J;
      Vector3 beam1(0.0, 0.0, 1.0);
      Vector3 p1xp3 = beam1.cross(p3.vector3());
      Vector3 p4xp5 = p4.vector3().cross(p5.vector3());
      const double cospsi3 = p1xp3.dot(p4xp5)/p1xp3.mod()/p4xp5.mod();
      const double f3 = p3.mass()/m4J;
      const double f4 = p4.mass()/m4J;
      const double f5 = p5.mass()/m4J;
      const double fA = pA.mass()/m4J;
      const double fB = pB.mass()/m4J;
      const double XA = pA.E()/(pA.E()+pB.E());
      FourMomentum pAB = pA+pB;
      Vector3 pABxp1 = pAB.vector3().cross(beam1);
      Vector3 pAxpB = pA.vector3().cross(pB.vector3());
      const double cospsiAB = pAxpB.dot(pABxp1)/pAxpB.mod()/pABxp1.mod();
      
      _h_4_mNJ->fill(m4J, weight);
      _h_4_X3->fill(X3, weight);
      _h_4_X4->fill(X4, weight);
      _h_4_costheta3->fill(costheta3, weight);
      _h_4_psi3->fill(mapAngle0ToPi(acos(cospsi3)), weight);
      _h_4_f3->fill(f3, weight);
      _h_4_f4->fill(f4, weight);
      _h_4_f5->fill(f5, weight);
      _h_4_XA->fill(XA, weight);
      _h_4_psiAB->fill(mapAngle0ToPi(acos(cospsiAB)), weight);
      _h_4_fA->fill(fA, weight);
      _h_4_fB->fill(fB, weight);
    }
      
      
    void _fiveJetAnalysis(const Jets& jets, const double& weight) {
      getLog() << Log::DEBUG << "5 jet analysis" << std::endl;
      FourMomentum jjjjj(0.0, 0.0, 0.0, 0.0);
      vector<FourMomentum> jetmoms;
      for (size_t i=0; i<5; ++i) {
        jetmoms.push_back(jets[i].momentum());
        jjjjj += jets[i].momentum();
      }
      const double m5J = jjjjj.mass();
      if (m5J < 750*GeV) return;
      
      FourMomentum pA, pB, pC, pD;
      vector<FourMomentum> jetmoms4(_reduce(jetmoms, pC, pD));
      vector<FourMomentum> jetmoms3(_reduce(jetmoms4, pA, pB));
      
      LorentzTransform cms_boost(-jjjjj.boostVector());
      vector<FourMomentum> jetmoms3_boosted;
      foreach (FourMomentum mom, jetmoms3) {
        jetmoms3_boosted.push_back(cms_boost.transform(mom));
      }
      pA = cms_boost.transform(pA);
      pB = cms_boost.transform(pB);
      pC = cms_boost.transform(pC);
      pD = cms_boost.transform(pD);
      
      sort(jetmoms3_boosted.begin(), jetmoms3_boosted.end(), FourMomentum::byEDescending());
      if (pB.E()>pA.E()) std::swap(pA, pB);
      if (pD.E()>pC.E()) std::swap(pD, pC);
      FourMomentum p3(jetmoms3_boosted[0]);
      FourMomentum p4(jetmoms3_boosted[1]);
      FourMomentum p5(jetmoms3_boosted[2]);
      
      // fill histograms
      const double costheta3 = cos(p3.theta());
      const double X3 = 2.0*p3.E()/m5J;
      const double X4 = 2.0*p4.E()/m5J;
      Vector3 beam1(0.0, 0.0, 1.0);
      Vector3 p1xp3 = beam1.cross(p3.vector3());
      Vector3 p4xp5 = p4.vector3().cross(p5.vector3());
      const double cospsi3 = p1xp3.dot(p4xp5)/p1xp3.mod()/p4xp5.mod();
      const double f3 = p3.mass()/m5J;
      const double f4 = p4.mass()/m5J;
      const double f5 = p5.mass()/m5J;
      const double fA = pA.mass()/m5J;
      const double fB = pB.mass()/m5J;
      const double XA = pA.E()/(pA.E()+pB.E());
      FourMomentum pAB = pA+pB;
      Vector3 pABxp1 = pAB.vector3().cross(beam1);
      Vector3 pAxpB = pA.vector3().cross(pB.vector3());
      const double cospsiAB = pAxpB.dot(pABxp1)/pAxpB.mod()/pABxp1.mod();
      const double fC = pC.mass()/m5J;
      const double fD = pD.mass()/m5J;
      const double XC = pC.E()/(pC.E()+pD.E());
      FourMomentum pCD = pC+pD;
      Vector3 pCDxp1 = pCD.vector3().cross(beam1);
      Vector3 pCxpD = pC.vector3().cross(pD.vector3());
      const double cospsiCD = pCxpD.dot(pCDxp1)/pCxpD.mod()/pCDxp1.mod();
      
      _h_5_mNJ->fill(m5J, weight);
      _h_5_X3->fill(X3, weight);
      _h_5_X4->fill(X4, weight);
      _h_5_costheta3->fill(costheta3, weight);
      _h_5_psi3->fill(mapAngle0ToPi(acos(cospsi3)), weight);
      _h_5_f3->fill(f3, weight);
      _h_5_f4->fill(f4, weight);
      _h_5_f5->fill(f5, weight);
      _h_5_XA->fill(XA, weight);
      _h_5_psiAB->fill(mapAngle0ToPi(acos(cospsiAB)), weight);
      _h_5_fA->fill(fA, weight);
      _h_5_fB->fill(fB, weight);
      _h_5_XC->fill(XC, weight);
      _h_5_psiCD->fill(mapAngle0ToPi(acos(cospsiCD)), weight);
      _h_5_fC->fill(fC, weight);
      _h_5_fD->fill(fD, weight);
    }
      
      
    /// Normalise histograms etc., after the run
    void finalize() {
      
      /// Normalise, scale and otherwise manipulate histograms here
      scale(_h_3_mNJ, crossSection()/sumOfWeights());
      scale(_h_3_X3, crossSection()/sumOfWeights());
      scale(_h_3_X4, crossSection()/sumOfWeights());
      scale(_h_3_costheta3, crossSection()/sumOfWeights());
      scale(_h_3_psi3, crossSection()/sumOfWeights());
      scale(_h_3_f3, crossSection()/sumOfWeights());
      scale(_h_3_f4, crossSection()/sumOfWeights());
      scale(_h_3_f5, crossSection()/sumOfWeights());
      
      scale(_h_4_mNJ, crossSection()/sumOfWeights());
      scale(_h_4_X3, crossSection()/sumOfWeights());
      scale(_h_4_X4, crossSection()/sumOfWeights());
      scale(_h_4_costheta3, crossSection()/sumOfWeights());
      scale(_h_4_psi3, crossSection()/sumOfWeights());
      scale(_h_4_f3, crossSection()/sumOfWeights());
      scale(_h_4_f4, crossSection()/sumOfWeights());
      scale(_h_4_f5, crossSection()/sumOfWeights());
      scale(_h_4_XA, crossSection()/sumOfWeights());
      scale(_h_4_psiAB, crossSection()/sumOfWeights());
      scale(_h_4_fA, crossSection()/sumOfWeights());
      scale(_h_4_fB, crossSection()/sumOfWeights());
      
      scale(_h_5_mNJ, crossSection()/sumOfWeights());
      scale(_h_5_X3, crossSection()/sumOfWeights());
      scale(_h_5_X4, crossSection()/sumOfWeights());
      scale(_h_5_costheta3, crossSection()/sumOfWeights());
      scale(_h_5_psi3, crossSection()/sumOfWeights());
      scale(_h_5_f3, crossSection()/sumOfWeights());
      scale(_h_5_f4, crossSection()/sumOfWeights());
      scale(_h_5_f5, crossSection()/sumOfWeights());
      scale(_h_5_XA, crossSection()/sumOfWeights());
      scale(_h_5_XC, crossSection()/sumOfWeights());
      scale(_h_5_psiAB, crossSection()/sumOfWeights());
      scale(_h_5_psiCD, crossSection()/sumOfWeights());
      scale(_h_5_fA, crossSection()/sumOfWeights());
      scale(_h_5_fB, crossSection()/sumOfWeights());
      scale(_h_5_fC, crossSection()/sumOfWeights());
      scale(_h_5_fD, crossSection()/sumOfWeights());
      
    }

    //@}


  private:
    vector<FourMomentum> _reduce(const vector<FourMomentum>& jets,
                                 FourMomentum& combined1,
                                 FourMomentum& combined2) {
      double minMass2 = 1e9;
      FourMomentum combined;
      size_t idx1(jets.size()), idx2(jets.size());
      for (size_t i=0; i<jets.size(); ++i) {
        for (size_t j=i+1; j<jets.size(); ++j) {
          double mass2 = FourMomentum(jets[i]+jets[j]).mass2();
          if (mass2<minMass2) {
            combined = jets[i]+jets[j];
            idx1=i;
            idx2=j;
          }
        }
      }
      vector<FourMomentum> newjets;
      for (size_t i=0; i<jets.size(); ++i) {
        if (i!=idx1 && i!=idx2) newjets.push_back(jets[i]);
      }
      newjets.push_back(combined);
      combined1 = jets[idx1];
      combined2 = jets[idx2];
      return newjets;
    }
    

  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_h_3_mNJ;
    AIDA::IHistogram1D *_h_3_X3;
    AIDA::IHistogram1D *_h_3_X4;
    AIDA::IHistogram1D *_h_3_costheta3;
    AIDA::IHistogram1D *_h_3_psi3;
    AIDA::IHistogram1D *_h_3_f3;
    AIDA::IHistogram1D *_h_3_f4;
    AIDA::IHistogram1D *_h_3_f5;
    
    AIDA::IHistogram1D *_h_4_mNJ;
    AIDA::IHistogram1D *_h_4_X3;
    AIDA::IHistogram1D *_h_4_X4;
    AIDA::IHistogram1D *_h_4_costheta3;
    AIDA::IHistogram1D *_h_4_psi3;
    AIDA::IHistogram1D *_h_4_f3;
    AIDA::IHistogram1D *_h_4_f4;
    AIDA::IHistogram1D *_h_4_f5;
    AIDA::IHistogram1D *_h_4_XA;
    AIDA::IHistogram1D *_h_4_psiAB;
    AIDA::IHistogram1D *_h_4_fA;
    AIDA::IHistogram1D *_h_4_fB;
    
    AIDA::IHistogram1D *_h_5_mNJ;
    AIDA::IHistogram1D *_h_5_X3;
    AIDA::IHistogram1D *_h_5_X4;
    AIDA::IHistogram1D *_h_5_costheta3;
    AIDA::IHistogram1D *_h_5_psi3;
    AIDA::IHistogram1D *_h_5_f3;
    AIDA::IHistogram1D *_h_5_f4;
    AIDA::IHistogram1D *_h_5_f5;
    AIDA::IHistogram1D *_h_5_XA;
    AIDA::IHistogram1D *_h_5_XC;
    AIDA::IHistogram1D *_h_5_psiAB;
    AIDA::IHistogram1D *_h_5_psiCD;
    AIDA::IHistogram1D *_h_5_fA;
    AIDA::IHistogram1D *_h_5_fB;
    AIDA::IHistogram1D *_h_5_fC;
    AIDA::IHistogram1D *_h_5_fD;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_1996_S3349578> plugin_CDF_1996_S3349578;


}
