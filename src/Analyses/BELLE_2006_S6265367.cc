// -*- C++ -*-
#include <iostream>
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief BELLE charmed mesons and baryons from fragmentation
  /// @author Eike von Seggern
  class BELLE_2006_S6265367 : public Analysis {
  public:

    BELLE_2006_S6265367()
      : Analysis("BELLE_2006_S6265367")
    { }


    void analyze(const Event& e) {
      const double weight = e.weight();

      // Loop through unstable FS particles and look for charmed mesons/baryons
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");

      const Beam beamproj = applyProjection<Beam>(e, "Beams");
      const ParticlePair& beams = beamproj.beams();
      FourMomentum mom_tot = beams.first.momentum() + beams.second.momentum();
      LorentzTransform cms_boost(-mom_tot.boostVector());
      const double s = sqr(beamproj.sqrtS());

      const bool onresonance = fuzzyEquals(beamproj.sqrtS(), 10.58, 5E-3);

      // Particle masses from PDGlive (accessed online 16. Nov. 2009).
      foreach (const Particle& p, ufs.particles()) {

        double xp = 0.0;
        double mH2 = 0.0;
        // 3-momentum in CMS frame
        const double mom = cms_boost.transform(p.momentum()).vector3().mod();

        const int PdgId = abs(p.pdgId());
        MSG_DEBUG("pdgID = " << PdgId << "  mom = " << mom);
        switch (PdgId) {

          case 421:
            MSG_DEBUG("D0 found");
            mH2 = 3.47763; // 1.86484^2
            xp = mom/sqrt(s/4.0 - mH2);
	    // off-resonance cross section
	    if(!onresonance) _sigmaD0->fill(10.6,weight);
	    if(checkDecay(p.genParticle())) {
	      if (onresonance)
		_histXpD0_R->fill(xp, s*weight);
	      else
		_histXpD0_C->fill(xp, s*weight);
	    }
	    if (onresonance)
	      _histXpD0_R_N->fill(xp, weight);
	    else
	      _histXpD0_C_N->fill(xp, weight);
            break;

          case 411:
            MSG_DEBUG("D+ found");
            mH2 = 3.49547; // 1.86962^2
            xp = mom/sqrt(s/4.0 - mH2);
	    if(!onresonance) _sigmaDPlus->fill(10.6,weight);
	    if(checkDecay(p.genParticle())) {
	      if (onresonance)
		_histXpDplus_R->fill(xp, s*weight);
	      else 
		_histXpDplus_C->fill(xp, s*weight);
	    }
	    if (onresonance)
	      _histXpDplus_R_N->fill(xp, weight);
	    else 
	      _histXpDplus_C_N->fill(xp, weight);
            break;
          case 431:
            MSG_DEBUG("D+_s found");
            mH2 = 3.87495; // 1.96849^2
            xp = mom/sqrt(s/4.0 - mH2);
	    if(!onresonance) _sigmaDs->fill(10.6,weight);
	    if(checkDecay(p.genParticle())) {
	      if (onresonance)
		_histXpDplus_s_R->fill(xp, s*weight);
	      else
		_histXpDplus_s_C->fill(xp, s*weight);
	    }
	    if (onresonance)
	      _histXpDplus_s_R_N->fill(xp, weight);
	    else
	      _histXpDplus_s_C_N->fill(xp, weight);
            break;

          case 4122:
            MSG_DEBUG("Lambda_c found");
            mH2 = 5.22780; // 2.28646^2
            xp = mom/sqrt(s/4.0 - mH2);
	    if(!onresonance) _sigmaLambdac   ->fill(10.6,weight);
	    if(checkDecay(p.genParticle())) {
	      if (onresonance)
		_histXpLambda_c_R->fill(xp, s*weight);
	      else
		_histXpLambda_c_C->fill(xp, s*weight);
	    }
	    if (onresonance)
	      _histXpLambda_c_R_N->fill(xp, weight);
	    else
	      _histXpLambda_c_C_N->fill(xp, weight);
            break;


          case 413: {
            MSG_DEBUG("D*+ found");
            mH2 = 4.04119; // 2.01027^2
            xp = mom/sqrt(s/4.0 - mH2);
	    if(!onresonance) {
	      _sigmaDStarPlusA->fill(10.6,weight);
	      _sigmaDStarPlusB->fill(10.6,weight);
	      _sigmaDStarPlusC->fill(10.6,weight);
	    }
	    const GenParticle * Dmeson = &p.genParticle();
            const GenVertex* dv = p.genParticle().end_vertex();
            bool D0decay(false), Pi0decay(false), Piplusdecay(false), Dplusdecay(false);

            for (GenVertex::particles_out_const_iterator
		   pp = dv->particles_out_const_begin();
                 pp != dv->particles_out_const_end(); ++pp) {
              if (abs((*pp)->pdg_id()) == 421) {
		Dmeson = *pp;
                D0decay = true;
              } else if (abs((*pp)->pdg_id()) == 411) {
		Dmeson = *pp;
                Dplusdecay = true;
              } else if (abs((*pp)->pdg_id()) == 111) {
                Pi0decay = true;
              } else if (abs((*pp)->pdg_id()) == 211) {
                Piplusdecay = true;
              }
            }

            if (D0decay && Piplusdecay && checkDecay(*Dmeson)) {
              if (onresonance)
                _histXpDstarplus2D0_R->fill(xp, s*weight);
              else
                _histXpDstarplus2D0_C->fill(xp, s*weight);
            } 
	    else if (Dplusdecay && Pi0decay && checkDecay(*Dmeson)) {
              if (onresonance)
                _histXpDstarplus2Dplus_R->fill(xp, s*weight);
              else
                _histXpDstarplus2Dplus_C->fill(xp, s*weight);
            } 
	    else {
              MSG_WARNING("Unexpected D* decay!");
            }
	    if (onresonance) {
	      _histXpDstarplus2D0_R_N->fill(xp, weight);
	      _histXpDstarplus2Dplus_R_N->fill(xp, weight);
	    }
	    else {
	      _histXpDstarplus2D0_C_N->fill(xp, weight);
	      _histXpDstarplus2Dplus_C_N->fill(xp, weight);
	    }
            break;
            }

          case 423:
            MSG_DEBUG("D*0 found");
            mH2 = 4.02793; // 2.00697**2
            xp = mom/sqrt(s/4.0 - mH2);
	    if(!onresonance) _sigmaDStar0    ->fill(10.6,weight);
            MSG_DEBUG("xp = " << xp);
	    const GenParticle * Dmeson = &p.genParticle();
	    const GenVertex* dv = p.genParticle().end_vertex();
	    bool D0decay(false), Pi0decay(false);
            for (GenVertex::particles_out_const_iterator
		   pp = dv->particles_out_const_begin();
                 pp != dv->particles_out_const_end(); ++pp) {
              if (abs((*pp)->pdg_id()) == 421) {
		Dmeson = *pp;
                D0decay = true;
              }
	      else if (abs((*pp)->pdg_id()) == 111) {
		Pi0decay = true;
	      }
	    }
            if (D0decay && Pi0decay && checkDecay(*Dmeson)) {
	      if (onresonance)
		_histXpDstar0_R->fill(xp, s*weight);
	      else {
		_histXpDstar0_C->fill(xp, s*weight);
	      }
	    }
	    if (onresonance)
	      _histXpDstar0_R_N->fill(xp, weight);
	    else {
	      _histXpDstar0_C_N->fill(xp, weight);
	    }
            break;
        }
      }
    } // analyze


    void finalize() {
      scale(_histXpDstarplus2D0_R, crossSection()/nanobarn/sumOfWeights());
      scale(_histXpD0_R, crossSection()/nanobarn/sumOfWeights());
      scale(_histXpDplus_R, crossSection()/nanobarn/sumOfWeights());
      scale(_histXpDplus_s_R, crossSection()/nanobarn/sumOfWeights());
      scale(_histXpLambda_c_R, crossSection()/nanobarn/sumOfWeights());
      scale(_histXpDstarplus2Dplus_R, crossSection()/nanobarn/sumOfWeights());
      scale(_histXpDstar0_R, crossSection()/nanobarn/sumOfWeights());

      scale(_histXpDstarplus2D0_C, crossSection()/nanobarn/sumOfWeights());
      scale(_histXpD0_C, crossSection()/nanobarn/sumOfWeights());
      scale(_histXpDplus_C, crossSection()/nanobarn/sumOfWeights());
      scale(_histXpDplus_s_C, crossSection()/nanobarn/sumOfWeights());
      scale(_histXpLambda_c_C, crossSection()/nanobarn/sumOfWeights());
      scale(_histXpDstarplus2Dplus_C, crossSection()/nanobarn/sumOfWeights());
      scale(_histXpDstar0_C, crossSection()/nanobarn/sumOfWeights());

      normalize(_histXpDstarplus2D0_R_N);
      normalize(_histXpD0_R_N);
      normalize(_histXpDplus_R_N);
      normalize(_histXpDplus_s_R_N);
      normalize(_histXpLambda_c_R_N);
      normalize(_histXpDstarplus2Dplus_R_N);
      normalize(_histXpDstar0_R_N);

      normalize(_histXpDstarplus2D0_C_N);
      normalize(_histXpD0_C_N);
      normalize(_histXpDplus_C_N);
      normalize(_histXpDplus_s_C_N);
      normalize(_histXpLambda_c_C_N);
      normalize(_histXpDstarplus2Dplus_C_N);
      normalize(_histXpDstar0_C_N);

      scale(_sigmaD0, crossSection()/picobarn/sumOfWeights());
      scale(_sigmaDPlus, crossSection()/picobarn/sumOfWeights());
      scale(_sigmaDs, crossSection()/picobarn/sumOfWeights());
      scale(_sigmaLambdac, crossSection()/picobarn/sumOfWeights());
      scale(_sigmaDStar0, crossSection()/picobarn/sumOfWeights());
      scale(_sigmaDStarPlusA, crossSection()/picobarn/sumOfWeights());
      scale(_sigmaDStarPlusB, crossSection()/picobarn/sumOfWeights());
      scale(_sigmaDStarPlusC, crossSection()/picobarn/sumOfWeights());
    } // finalize


    void init() {
      addProjection(Beam(), "Beams");
      addProjection(UnstableFinalState(-1.3170, 1.9008), "UFS");

      // continuum cross sections
      _sigmaD0         = bookHistogram1D(1,1,1);
      _sigmaDPlus      = bookHistogram1D(1,1,2);
      _sigmaDs         = bookHistogram1D(1,1,3);
      _sigmaLambdac    = bookHistogram1D(1,1,4);
      _sigmaDStar0     = bookHistogram1D(1,1,5);
      _sigmaDStarPlusA = bookHistogram1D(1,1,6);
      _sigmaDStarPlusB = bookHistogram1D(1,1,7);
      _sigmaDStarPlusC = bookHistogram1D(1,1,8);

      // histograms for continuum data (sqrt(s) = 10.52 GeV)
      _histXpDstarplus2D0_C = bookHistogram1D(2, 1, 1);
      _histXpD0_C = bookHistogram1D(3, 1, 1);
      _histXpDplus_C = bookHistogram1D(4, 1, 1);
      _histXpDplus_s_C = bookHistogram1D(5, 1, 1);
      _histXpLambda_c_C = bookHistogram1D(6, 1, 1);
      _histXpDstarplus2Dplus_C = bookHistogram1D(7, 1, 1);
      _histXpDstar0_C = bookHistogram1D(8, 1, 1);

      // histograms for on-resonance data (sqrt(s) = 10.58 GeV)
      _histXpDstarplus2D0_R = bookHistogram1D(9, 1, 1);
      _histXpD0_R = bookHistogram1D(10, 1, 1);
      _histXpDplus_R = bookHistogram1D(11, 1, 1);
      _histXpDplus_s_R = bookHistogram1D(12, 1, 1);
      _histXpLambda_c_R = bookHistogram1D(13, 1, 1);
      _histXpDstarplus2Dplus_R = bookHistogram1D(14, 1, 1);
      _histXpDstar0_R = bookHistogram1D(15, 1, 1);

      // histograms for continuum data (sqrt(s) = 10.52 GeV)
      _histXpDstarplus2D0_C_N = bookHistogram1D(2, 1, 2);
      _histXpD0_C_N = bookHistogram1D(3, 1, 2);
      _histXpDplus_C_N = bookHistogram1D(4, 1, 2);
      _histXpDplus_s_C_N = bookHistogram1D(5, 1, 2);
      _histXpLambda_c_C_N = bookHistogram1D(6, 1, 2);
      _histXpDstarplus2Dplus_C_N = bookHistogram1D(7, 1, 2);
      _histXpDstar0_C_N = bookHistogram1D(8, 1, 2);

      // histograms for on-resonance data (sqrt(s) = 10.58 GeV)
      _histXpDstarplus2D0_R_N = bookHistogram1D(9, 1, 2);
      _histXpD0_R_N = bookHistogram1D(10, 1, 2);
      _histXpDplus_R_N = bookHistogram1D(11, 1, 2);
      _histXpDplus_s_R_N = bookHistogram1D(12, 1, 2);
      _histXpLambda_c_R_N = bookHistogram1D(13, 1, 2);
      _histXpDstarplus2Dplus_R_N = bookHistogram1D(14, 1, 2);
      _histXpDstar0_R_N = bookHistogram1D(15, 1, 2);

    } // init

  private:

    //@{
    /// Histograms

    // Histograms for the continuum cross sections
    AIDA::IHistogram1D* _sigmaD0;
    AIDA::IHistogram1D* _sigmaDPlus;
    AIDA::IHistogram1D* _sigmaDs;
    AIDA::IHistogram1D* _sigmaLambdac;
    AIDA::IHistogram1D* _sigmaDStar0;
    AIDA::IHistogram1D* _sigmaDStarPlusA;
    AIDA::IHistogram1D* _sigmaDStarPlusB;
    AIDA::IHistogram1D* _sigmaDStarPlusC;

    // histograms for continuum data (sqrt(s) = 10.52 GeV)
    AIDA::IHistogram1D* _histXpDstarplus2D0_C;
    AIDA::IHistogram1D* _histXpD0_C;
    AIDA::IHistogram1D* _histXpDplus_C;
    AIDA::IHistogram1D* _histXpDplus_s_C;
    AIDA::IHistogram1D* _histXpLambda_c_C;
    AIDA::IHistogram1D* _histXpDstarplus2Dplus_C;
    AIDA::IHistogram1D* _histXpDstar0_C;
    AIDA::IHistogram1D* _histXpDstarplus2D0_C_N;
    AIDA::IHistogram1D* _histXpD0_C_N;
    AIDA::IHistogram1D* _histXpDplus_C_N;
    AIDA::IHistogram1D* _histXpDplus_s_C_N;
    AIDA::IHistogram1D* _histXpLambda_c_C_N;
    AIDA::IHistogram1D* _histXpDstarplus2Dplus_C_N;
    AIDA::IHistogram1D* _histXpDstar0_C_N;

    // histograms for on-resonance data (sqrt(s) = 10.58 GeV)
    AIDA::IHistogram1D* _histXpDstarplus2D0_R;
    AIDA::IHistogram1D* _histXpD0_R;
    AIDA::IHistogram1D* _histXpDplus_R;
    AIDA::IHistogram1D* _histXpDplus_s_R;
    AIDA::IHistogram1D* _histXpLambda_c_R;
    AIDA::IHistogram1D* _histXpDstarplus2Dplus_R;
    AIDA::IHistogram1D* _histXpDstar0_R;
    AIDA::IHistogram1D* _histXpDstarplus2D0_R_N;
    AIDA::IHistogram1D* _histXpD0_R_N;
    AIDA::IHistogram1D* _histXpDplus_R_N;
    AIDA::IHistogram1D* _histXpDplus_s_R_N;
    AIDA::IHistogram1D* _histXpLambda_c_R_N;
    AIDA::IHistogram1D* _histXpDstarplus2Dplus_R_N;
    AIDA::IHistogram1D* _histXpDstar0_R_N;
    //@}

    bool checkDecay(const GenParticle & p) {
      unsigned int nstable=0,npip=0,npim=0;
      unsigned int np=0,nap=0,nKp=0,nKm=0,nPhi=0;
      findDecayProducts(p,nstable,npip,npim,
			np,nap,nKp,nKm,nPhi);
      int id = p.pdg_id();
      //D0 
      if(id==421) {
	if(nstable==2&&nKm==1&&npip==1) return true;
      }
      //Dbar0
      else if(id==-421) {
	if(nstable==2&&nKp==1&&npim==1) return true;
      }
      // D+
      else if(id==411) {
	if(nstable==3&&nKm==1&&npip==2) return true;
      }
      // D-
      else if(id==-411) {
	if(nstable==3&&nKp==1&&npim==2) return true;
      }
      // D_s+
      else if(id==431) {
	if(nstable==1&&nPhi==1&&npip==1) return true;
      }
      // D_s-
      else if(id==-431) {
	if(nstable==1&&nPhi==1&&npim==1) return true;
      }
      // Lambda_c
      else if(id==4122) {
	if(nstable==3&&np==1&&npip==1&&nKm==1) return true;
      }
      // Lambda_c bar
      else if(id==-4122) {
	if(nstable==3&&nap==1&&npim==1&&nKp==1) return true;
      }
      return false;
    }

    void findDecayProducts(const GenParticle & p,
			   unsigned int & nstable, unsigned int & npip,
			   unsigned int & npim   , unsigned int & np,
			   unsigned int & nap    , unsigned int & nKp,
			   unsigned int & nKm    , unsigned int & nPhi) {
      const GenVertex* dv = p.end_vertex();
      for (GenVertex::particles_out_const_iterator
	     pp = dv->particles_out_const_begin();
	   pp != dv->particles_out_const_end(); ++pp) {
	int id = (*pp)->pdg_id();
	if(id==333) 
	  ++nPhi;
	else if(id==111||id==221)
	  ++nstable;
	else if((*pp)->end_vertex())
	  findDecayProducts(**pp,nstable,npip,npim,np,nap,nKp,nKm,nPhi);
	else {
	  if(id!=22) ++nstable;
	  if     (id ==   211) ++npip;
	  else if(id ==  -211) ++npim;
	  else if(id ==  2212) ++np;
	  else if(id == -2212) ++nap;
	  else if(id ==   321) ++nKp;
	  else if(id ==  -321) ++nKm;
	}
      }
    }
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BELLE_2006_S6265367);

}
