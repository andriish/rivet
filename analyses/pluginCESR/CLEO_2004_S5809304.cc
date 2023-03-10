// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief CLEO charmed mesons and baryons from fragmentation
  ///
  /// @author Peter Richardson
  class CLEO_2004_S5809304 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2004_S5809304);

    void init() {
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      // continuum cross sections
      book(_sigmaDPlus      ,1,1,1);
      book(_sigmaD0A        ,1,1,2);
      book(_sigmaD0B        ,1,1,3);
      book(_sigmaDStarPlusA ,1,1,4);
      book(_sigmaDStarPlusB ,1,1,5);
      book(_sigmaDStar0A    ,1,1,6);
      book(_sigmaDStar0B    ,1,1,7);

       // histograms for continuum data
      book(_histXpDplus      ,2, 1, 1);
      book(_histXpD0A        ,3, 1, 1);
      book(_histXpD0B        ,4, 1, 1);
      book(_histXpDStarPlusA ,5, 1, 1);
      book(_histXpDStarPlusB ,6, 1, 1);
      book(_histXpDStar0A    ,7, 1, 1);
      book(_histXpDStar0B    ,8, 1, 1);
      book(_histXpTotal      ,9, 1, 1);

    }


    void analyze(const Event& e) {
      // Loop through unstable FS particles and look for charmed mesons/baryons
      const UnstableParticles& ufs = apply<UnstableParticles>(e, "UFS");

      const Beam beamproj = apply<Beam>(e, "Beams");
      const ParticlePair& beams = beamproj.beams();
      const FourMomentum mom_tot = beams.first.momentum() + beams.second.momentum();
      LorentzTransform cms_boost;
      if (mom_tot.p3().mod() > 1*MeV)
        cms_boost = LorentzTransform::mkFrameTransformFromBeta(mom_tot.betaVec());
      const double s = sqr(beamproj.sqrtS());

      // Particle masses from PDGlive (accessed online 16. Nov. 2009).
      for (const Particle& p : ufs.particles()) {

        double xp = 0.0;
        double mH2 = 0.0;
        // 3-momentum in CMS frame
        const double mom = cms_boost.transform(p.momentum()).vector3().mod();

        const int pdgid = p.abspid();
        MSG_DEBUG("pdgID = " << pdgid << "  mom = " << mom);
        switch (pdgid) {

        case 421:
          MSG_DEBUG("D0 found");
          mH2 = 3.47763; // 1.86484^2
          xp = mom/sqrt(s/4.0 - mH2);
          _sigmaD0A->fill(10.6);
          _sigmaD0B->fill(10.6);
          _histXpD0A->fill(xp);
          _histXpD0B->fill(xp);
          _histXpTotal->fill(xp);
          break;

        case 411:
          MSG_DEBUG("D+ found");
          mH2 = 3.49547; // 1.86962^2
          xp = mom/sqrt(s/4.0 - mH2);
          _sigmaDPlus->fill(10.6);
          _histXpDplus->fill(xp);
          _histXpTotal->fill(xp);
          break;

        case 413:
          MSG_DEBUG("D*+ found");
          mH2 = 4.04119; // 2.01027^2
          xp = mom/sqrt(s/4.0 - mH2);
          _sigmaDStarPlusA->fill(10.6);
          _sigmaDStarPlusB->fill(10.6);
          _histXpDStarPlusA->fill(xp);
          _histXpDStarPlusB->fill(xp);
          _histXpTotal->fill(xp);
          break;

        case 423:
          MSG_DEBUG("D*0 found");
          mH2 = 4.02793; // 2.00697**2
          xp = mom/sqrt(s/4.0 - mH2);
          _sigmaDStar0A->fill(10.6);
          _sigmaDStar0B->fill(10.6);
          _histXpDStar0A->fill(xp);
          _histXpDStar0B->fill(xp);
          _histXpTotal->fill(xp);
          break;
        }

      }
    }


    void finalize() {

      scale(_sigmaDPlus     , crossSection()/picobarn/sumOfWeights());
      scale(_sigmaD0A       , crossSection()/picobarn/sumOfWeights());
      scale(_sigmaD0B       , crossSection()/picobarn/sumOfWeights());
      scale(_sigmaDStarPlusA, crossSection()/picobarn/sumOfWeights());
      scale(_sigmaDStarPlusB, crossSection()/picobarn/sumOfWeights());
      scale(_sigmaDStar0A   , crossSection()/picobarn/sumOfWeights());
      scale(_sigmaDStar0B   , crossSection()/picobarn/sumOfWeights());

      scale(_histXpDplus     , crossSection()/picobarn/sumOfWeights());
      scale(_histXpD0A       , crossSection()/picobarn/sumOfWeights());
      scale(_histXpD0B       , crossSection()/picobarn/sumOfWeights());
      scale(_histXpDStarPlusA, crossSection()/picobarn/sumOfWeights());
      scale(_histXpDStarPlusB, crossSection()/picobarn/sumOfWeights());
      scale(_histXpDStar0A   , crossSection()/picobarn/sumOfWeights());
      scale(_histXpDStar0B   , crossSection()/picobarn/sumOfWeights());
      scale(_histXpTotal     , crossSection()/picobarn/sumOfWeights()/4.);
    }


  private:

    // Histograms for the continuum cross sections
    Histo1DPtr _sigmaDPlus     ;
    Histo1DPtr _sigmaD0A       ;
    Histo1DPtr _sigmaD0B       ;
    Histo1DPtr _sigmaDStarPlusA;
    Histo1DPtr _sigmaDStarPlusB;
    Histo1DPtr _sigmaDStar0A   ;
    Histo1DPtr _sigmaDStar0B   ;

    // Histograms for continuum data
    Histo1DPtr _histXpDplus     ;
    Histo1DPtr _histXpD0A       ;
    Histo1DPtr _histXpD0B       ;
    Histo1DPtr _histXpDStarPlusA;
    Histo1DPtr _histXpDStarPlusB;
    Histo1DPtr _histXpDStar0A   ;
    Histo1DPtr _histXpDStar0B   ;
    Histo1DPtr _histXpTotal     ;

  };


  RIVET_DECLARE_ALIASED_PLUGIN(CLEO_2004_S5809304, CLEO_2004_I645209);

}
