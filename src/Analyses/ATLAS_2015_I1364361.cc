/*
 * Rivet routine for the Higgs differential cross sections combination
 * between the ATLAS measurements in the yy and 4l channels for
 * Higgs transverse momentum, rapidity, jet multiplicity and leading jet pT
 *
 * Author: Michaela Queitsch-Maitland (ATLAS Collaboration)
 * Contact: Michaela Queitsch-Maitland <michaela.queitsch-maitland@cern.ch>,
 *          Dag Gillberg <dag.gillberg@cern.ch>,
 *          Florian Bernlochner <florian.bernlochner@cern.ch>,
 *          Sarah Heim <sarah.heim@cern.ch>
 */

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/Vector4.hh"
#include "Rivet/Particle.hh"
#include "HepMC/GenEvent.h"

namespace Rivet {


  class ATLAS_2015_I1364361 : public Analysis {
  public:

    /// Constructor
    ATLAS_2015_I1364361()
      : Analysis("ATLAS_2015_I1364361")
    {    }



    /// Book histograms and initialise projections before the run
    void init() {

      // All final state particles
      const FinalState FS(-MAXDOUBLE, MAXDOUBLE, 0.0*GeV);
      addProjection(FS,"FS");


      // Histograms with data bins
      // Inclusive
      _h_pTH_incl = bookHisto1D(1,1,1);// "pTH_incl", vector<double>{0, 20, 30, 40, 50, 60, 80, 100, 200}, "", "", "" );
      _h_yH_incl = bookHisto1D(2,1,1);// "yH_incl", vector<double>{0.0, 0.3, 0.6, 0.9, 1.2, 1.6, 2.4}, "", "", "" );
      _h_Njets_incl = bookHisto1D(3,1,1);// "Njets_incl", 4, 0, 4, "", "", "" );
      _h_pTj1_incl = bookHisto1D(4,1,1);// "pTj1_incl", vector<double>{0, 30, 50, 70, 100, 140}, "", "", "" );
      //// |y_H| < 2.0
      //_h_pTH_yH2 = bookHisto1D( "pTH_yH2", vector<double>{0, 20, 30, 40, 50, 60, 80, 100, 200}, "", "", "" );
      //_h_yH_yH2 = bookHisto1D( "yH_yH2", vector<double>{0.0, 0.3, 0.6, 0.9, 1.2, 1.6, 2.4}, "", "", "" );
      //_h_Njets_yH2 = bookHisto1D( "Njets_yH2", 4, 0, 4, "", "", "" );
      //_h_pTj1_yH2 = bookHisto1D( "pTj1_yH2", vector<double>{0, 30, 50, 70, 100, 140}, "", "", "" );
      //// |y_H| < 1.2
      //_h_pTH_yH12 = bookHisto1D( "pTH_yH12", vector<double>{0, 20, 30, 40, 50, 60, 80, 100, 200}, "", "", "" );
      //_h_yH_yH12 = bookHisto1D( "yH_yH12", vector<double>{0.0, 0.3, 0.6, 0.9, 1.2, 1.6, 2.4}, "", "", "" );
      //_h_Njets_yH12 = bookHisto1D( "Njets_yH12", 4, 0, 4, "", "", "" );
      //_h_pTj1_yH12 = bookHisto1D( "pTj1_yH12", vector<double>{0, 30, 50, 70, 100, 140}, "", "", "" );

      //// Histograms with fine bins
      //_h_pTH = bookHisto1D( "pTH", 100, 0, 500 , "", "", "" );
      //_h_yH = bookHisto1D( "yH", 80, -4, 4, "", "", "" );
      //_h_pTj1 = bookHisto1D( "pTj1", 100, 0, 500, "", "", "" );
      //_h_pTj2 = bookHisto1D( "pTj2", 100, 0, 500, "", "", "" );
      //_h_pTj3 = bookHisto1D( "pTj3", 100, 0, 500, "", "", "" );
      //// H+dijet kinematics
      //_h_m_jj     = bookHisto1D( "m_jj", 80, 0, 1000, "", "", "" );
      //_h_dphi_jj = bookHisto1D( "dph_jj", 60, 0, PI, "", "", "" );
      //_h_dy_jj    = bookHisto1D( "dy_jj", 100, 0, 10, "", "", "" );
      //_h_dphi_H_jj    = bookHisto1D( "dphi_H_jj", 60, 0, PI, "", "", "" );

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get event weight
      const double weight = event.weight();
      _weight = weight;

      // Get the final state particles ordered by pT
      const ParticleVector& FS = applyProjection<FinalState>(event, "FS").particlesByPt();

      // Find the Higgs
      bool stable_higgs = false;
      const Particle* higgs=0;
      foreach ( const Particle& p, FS ) {
	if ( p.pid()==25 ) {
	  //	  printf("Found the Higgs!!");
	  stable_higgs = true;
	  higgs = &p;
	  break;
	}
      }

      // If no stable Higgs found in event record, can't do anything (abort)
      if ( !stable_higgs ) {
	MSG_WARNING("FATAL: No stable Higgs found in event record.\n");
	vetoEvent;
      }

      ParticleVector leptons;
      ParticleVector photons; // for dressing
      ParticleVector jet_ptcls;

      // Loop over final state particles and fill jet particles vector
      foreach ( const Particle& ptcl, FS ) {
	// Do not include the Higgs in jet finding!
	if ( ptcl.pid()==25 ) continue;
	// Neutrinos not from hadronisation
	if ( ptcl.isNeutrino() && !fromHadronDecay(ptcl) ) continue;
	// Electrons and muons not from hadronisation
	if ( ( ptcl.abspid() == 11 || ptcl.abspid() == 13 ) && !fromHadronDecay(ptcl) ) {
	  leptons.push_back(ptcl);
	  continue;
	}
	// Photons not from hadronisation
	if ( ptcl.abspid() == 22 && !fromHadronDecay(ptcl) ) {
	  photons.push_back(ptcl);
	  continue;
	}
	// Add particle to jet inputs
	jet_ptcls.push_back(ptcl);
      }

      // Match FS photons to leptons within cone R=0.1
      // If they are not 'dressing' photons, add to jet particle vector
      foreach ( const Particle& ph, photons ) {
	bool fsr_photon = false;
	foreach ( const Particle& lep, leptons ) {
	  if ( deltaR(ph.momentum(),lep.momentum()) < 0.1 ){
	    fsr_photon=true;
	    continue;
	  }
	}
	if ( !fsr_photon ) jet_ptcls.push_back(ph);
      }

      // Let's build the jets!
      FastJets jet_pro(FastJets::ANTIKT, 0.4);
      jet_pro.calc(jet_ptcls);
      Jets jets = jet_pro.jetsByPt(Cuts::pT>30*GeV && Cuts::absrap<4.4);

      _pTH = higgs->momentum().pT();
      _yH = higgs->momentum().rapidity();
      _Njets = jets.size() > 3 ? 3 : jets.size();
      _pTj1 = jets.size() > 0 ? jets[0].momentum().pT() : 0.;
      _pTj2 = jets.size() > 1 ? jets[1].momentum().pT() : 0.;
      _pTj3 = jets.size() > 2 ? jets[2].momentum().pT() : 0.;

      _h_pTH_incl->fill(_pTH,weight);
      //_h_pTH->fill(_pTH,weight);
      _h_yH_incl->fill( fabs(_yH),weight);
      //_h_yH->fill(_yH,weight);
      _h_Njets_incl->fill(_Njets,weight);
      _h_pTj1_incl->fill(_pTj1,weight);

      //_h_pTj1->fill(_pTj1,weight);
      //_h_pTj2->fill(_pTj2,weight);
      //_h_pTj3->fill(_pTj3,weight);

      // H+dijet kinematics
      //if ( jets.size()>1 ) {
	//FourMomentum jet1 = jets[0].momentum();
	//FourMomentum jet2 = jets[1].momentum();
	//_mjj = ( jet1+jet2 ).mass();
	//_dyjj = fabs( jet1.rapidity()-jet2.rapidity() );
	//_dphijj = fabs( deltaPhi(jet1,jet2) );
	//_dphiHjj = fabs( deltaPhi(higgs->momentum(),jet1+jet2) );

	//_h_m_jj->fill(_mjj,weight);
	//_h_dy_jj->fill(_dyjj,weight);
	//_h_dphi_jj->fill(_dphijj,weight);
	//_h_dphi_H_jj->fill(_dphiHjj,weight);
      //}

      //if ( fabs(_yH) < 1.2 ) {
	//_h_pTH_yH12->fill(_pTH,weight);
	//_h_yH_yH12->fill( fabs(_yH),weight);
	//_h_Njets_yH12->fill(_Njets,weight);
	//_h_pTj1_yH12->fill(_pTj1,weight);
      //}

      //if ( fabs(_yH) < 2.0 ) {
	//_h_pTH_yH2->fill(_pTH,weight);
	//_h_yH_yH2->fill( fabs(_yH),weight);
	//_h_Njets_yH2->fill(_Njets,weight);
	//_h_pTj1_yH2->fill(_pTj1,weight);
      //}

    }


    /// Normalise histograms etc., after the run
    void finalize() {


      double xs = crossSectionPerEvent();
      scale( _h_pTH_incl, xs );
      scale( _h_yH_incl, xs );
      scale( _h_Njets_incl, xs );
      scale( _h_pTj1_incl, xs );
      //scale( _h_pTH_yH2, xs );
      //scale( _h_yH_yH2, xs );
      //scale( _h_Njets_yH2, xs );
      //scale( _h_pTj1_yH2, xs );
      //scale( _h_pTH_yH12, xs );
      //scale( _h_yH_yH12, xs );
      //scale( _h_Njets_yH12, xs );
      //scale( _h_pTj1_yH12, xs );

      //scale( _h_pTH, xs );
      //scale( _h_yH, xs );
      //scale( _h_pTj1, xs );
      //scale( _h_pTj2, xs );
      //scale( _h_pTj3, xs );
      //scale( _h_m_jj, xs );
      //scale( _h_dphi_jj, xs );
      //scale( _h_dy_jj, xs );
      //scale( _h_dphi_H_jj, xs );
    }


    bool fromHadronDecay(const Particle& p ) {
      return p.fromHadron();
      // const GenVertex* prodVtx = p.genParticle()->production_vertex();
      // if (prodVtx == NULL) return false;
      // foreach (const GenParticle* ancestor, particles(prodVtx, HepMC::ancestors)) {
      //   const PdgId pid = ancestor->pdg_id();
      //   if (ancestor->status() == 2 && PID::isHadron(pid)) return true;
      //   if (ancestor->status() == 2 && (abs(pid) == PID::TAU && fromHadronDecay(ancestor))) return true;
      // }
      // return false;
    }


  private:

    double _weight;
    double _pTH;
    double _yH;
    double _Njets;
    double _pTj1;
    double _pTj2;
    double _pTj3;
    double _mjj;
    double _dphijj;
    double _dphiHjj;
    double _dyjj;


    Histo1DPtr _h_pTH_incl;
    Histo1DPtr _h_yH_incl;
    Histo1DPtr _h_Njets_incl;
    Histo1DPtr _h_pTj1_incl;

    //Histo1DPtr _h_pTH_yH2;
    //Histo1DPtr _h_yH_yH2;
    //Histo1DPtr _h_Njets_yH2;
    //Histo1DPtr _h_pTj1_yH2;

    //Histo1DPtr _h_pTH_yH12;
    //Histo1DPtr _h_yH_yH12;
    //Histo1DPtr _h_Njets_yH12;
    //Histo1DPtr _h_pTj1_yH12;

    //Histo1DPtr _h_pTH;
    //Histo1DPtr _h_yH;
    //Histo1DPtr _h_pTj1;
    //Histo1DPtr _h_pTj2;
    //Histo1DPtr _h_pTj3;
    //Histo1DPtr _h_m_jj;
    //Histo1DPtr _h_dphi_jj;
    //Histo1DPtr _h_dy_jj;
    //Histo1DPtr _h_dphi_H_jj;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1364361);

}
