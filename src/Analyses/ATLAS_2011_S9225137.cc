// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/RivetMT2.hh"

namespace Rivet {


  class ATLAS_2011_S9225137 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ATLAS_2011_S9225137()
      : Analysis("ATLAS_2011_S9225137")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // veto region electrons
      std::vector<std::pair<double, double> > eta_v_e;
      eta_v_e.push_back(make_pair(-1.52,-1.37));
      eta_v_e.push_back(make_pair( 1.37, 1.52));
      IdentifiedFinalState veto_elecs(eta_v_e, 10.0*GeV);
      veto_elecs.acceptIdPair(ELECTRON);
      addProjection(veto_elecs, "veto_elecs");

      // projection to find the electrons
      std::vector<std::pair<double, double> > eta_e;
      eta_e.push_back(make_pair(-2.47,2.47));
      IdentifiedFinalState elecs(eta_e, 20.0*GeV);
      elecs.acceptIdPair(ELECTRON);
      addProjection(elecs, "elecs");

      // projection to find the muons
      std::vector<std::pair<double, double> > eta_m;
      eta_m.push_back(make_pair(-2.4,2.4));
      IdentifiedFinalState muons(eta_m, 10.0*GeV);
      muons.acceptIdPair(MUON);
      addProjection(muons, "muons");

      // for pTmiss
      addProjection(VisibleFinalState(-4.9,4.9),"vfs");

      VetoedFinalState vfs;
      vfs.addVetoPairId(MUON);

      /// Jet finder
      addProjection(FastJets(vfs, FastJets::ANTIKT, 0.4),
		    "AntiKtJets04");

      // all tracks (to do deltaR with leptons)
      addProjection(ChargedFinalState(-3.0,3.0),"cfs");

      /// Book histograms
      _etmissHTA = bookHistogram1D("etmissHTA", 64, 0., 16.);
      _etmissHTB = bookHistogram1D("etmissHTB", 64, 0., 16.);

      _njet55A = bookHistogram1D("njet55A", 14, 0.5, 14.5);
      _njet55B = bookHistogram1D("njet55B", 14, 0.5, 14.5);
      _njet80A = bookHistogram1D("njet80A", 14, 0.5, 14.5);
      _njet80B = bookHistogram1D("njet80B", 14, 0.5, 14.5);

      _count_7j55 = bookHistogram1D("count_7j55", 1, 0., 1.);
      _count_8j55 = bookHistogram1D("count_8j55", 1, 0., 1.);
      _count_6j80 = bookHistogram1D("count_6j80", 1, 0., 1.);
      _count_7j80 = bookHistogram1D("count_7j80", 1, 0., 1.);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // apply electron veto region
      ParticleVector veto_e
	= applyProjection<IdentifiedFinalState>(event, "veto_elecs").particles();
      if ( ! veto_e.empty() ) {
	MSG_DEBUG("electrons in veto region");
	vetoEvent;
      }

      // get the jet candidates
      Jets cand_jets;
      foreach (const Jet& jet,
	       applyProjection<FastJets>(event, "AntiKtJets04").jetsByPt(20.0*GeV) ) {
        if ( fabs( jet.momentum().eta() ) < 4.9 ) {
          cand_jets.push_back(jet);
        }
      }

      // candidate muons
      ParticleVector cand_mu;
      ParticleVector chg_tracks = 
	applyProjection<ChargedFinalState>(event, "cfs").particles();
      foreach ( const Particle & mu,
		applyProjection<IdentifiedFinalState>(event, "muons").particlesByPt() ) {
	double pTinCone = -mu.momentum().pT();
	foreach ( const Particle & track, chg_tracks ) {
	  if ( deltaR(mu.momentum(),track.momentum()) <= 0.2 )
	    pTinCone += track.momentum().pT();
	}
	if ( pTinCone < 1.8*GeV )
	  cand_mu.push_back(mu);
      }

      // candidate electrons

      ParticleVector cand_e  = 
	applyProjection<IdentifiedFinalState>(event, "elecs").particlesByPt();

      // resolve jet/lepton ambiguity 
      Jets cand_jets_2;
      foreach ( const Jet& jet, cand_jets ) {
	// candidates above eta=2.8 are jets
	if ( fabs( jet.momentum().eta() ) >= 2.8 )
	  cand_jets_2.push_back( jet );
	// otherwise more the R=0.2 from an electrons
	else {
	  bool away_from_e = true;
	  foreach ( const Particle & e, cand_e ) {
	    if ( deltaR(e.momentum(),jet.momentum()) <= 0.2 ) {
	      away_from_e = false;
	      break;
	    }
	  }
	  if ( away_from_e )
	    cand_jets_2.push_back( jet );
	}
      }

      // only keep electrons more than R=0.4 from jets
      ParticleVector recon_e;
      foreach ( const Particle & e, cand_e ) {
	bool away = true;
	foreach ( const Jet& jet, cand_jets_2 ) {
	  if ( deltaR(e.momentum(),jet.momentum()) < 0.4 ) {
	    away = false;
	    break;
	  }
	}
	if ( away )
	  recon_e.push_back( e );
      }

      // only keep muons more than R=0.4 from jets
      ParticleVector recon_mu;
      foreach ( const Particle & mu, cand_mu ) {
	bool away = true;
	foreach ( const Jet& jet, cand_jets_2 ) {
	  if ( deltaR(mu.momentum(),jet.momentum()) < 0.4 ) {
	    away = false;
	    break;
	  }
	}
	if ( away )
	  recon_mu.push_back( mu );
      }

      // pTmiss
      ParticleVector vfs_particles = 
	applyProjection<VisibleFinalState>(event, "vfs").particles();
      FourMomentum pTmiss;
      foreach ( const Particle & p, vfs_particles ) {
	pTmiss -= p.momentum();
      }
      double eTmiss = pTmiss.pT();

      // final jet filter
      Jets recon_jets;
      foreach ( const Jet& jet, cand_jets_2 ) {
	if ( fabs( jet.momentum().eta() ) <= 2.8 )
	  recon_jets.push_back( jet );
      }

      // now only use recon_jets, recon_mu, recon_e

      // reject events with electrons and muons
      if ( ! ( recon_mu.empty() && recon_e.empty() ) ) {
	MSG_DEBUG("Charged leptons left after selection");
	vetoEvent;
      }

      // calculate H_T
      double HT=0;
      foreach ( const Jet& jet, recon_jets ) {
	if ( jet.momentum().pT() > 40 * GeV )
	  HT += jet.momentum().pT() ;
      }

      // number of jets and deltaR
      bool pass55DeltaR=true;
      unsigned int njet55=0;
      bool pass80DeltaR=true;
      unsigned int njet80=0;
      for (unsigned int ix=0;ix<recon_jets.size();++ix) {
	if(recon_jets[ix].momentum().pT()>80.*GeV) ++njet80;
	if(recon_jets[ix].momentum().pT()>55.*GeV) ++njet55;

	for (unsigned int iy=ix+1;iy<recon_jets.size();++iy) {
	  if(recon_jets[ix].momentum().pT()>55.*GeV &&
	     recon_jets[iy].momentum().pT()>55.*GeV &&
	     deltaR(recon_jets[ix],recon_jets[ix]) <0.6 ) 
	    pass55DeltaR = false;
	  if(recon_jets[ix].momentum().pT()>80.*GeV &&
	     recon_jets[iy].momentum().pT()>80.*GeV &&
	     deltaR(recon_jets[ix],recon_jets[ix]) <0.6 ) 
	    pass80DeltaR = false;
	}
      }

      // plots of etmiss/ht
      double etht = eTmiss/sqrt(HT);
      if(njet55==6) _etmissHTA->fill(etht,weight);
      if(njet80==5) _etmissHTB->fill(etht,weight);

      if(etht>1.5&&etht<2. ) {
	_njet55A->fill(njet55,weight);
	_njet80A->fill(njet80,weight);
      }
      if(etht>2. &&etht<3. ) {
	_njet55B->fill(njet55,weight);
	_njet80B->fill(njet80,weight);
      }

      // apply E_T/sqrt(H_T) cut
      if(etht<=3.5*GeV) {
	MSG_DEBUG("Fails ET/sqrt(HT) cut ");
	vetoEvent;
      }

      // check passes at least one delta5/ njet number cut
      if(!(pass55DeltaR && njet55 >= 7) &&
	 !(pass80DeltaR && njet80 >= 6) ) {
	MSG_DEBUG("Fails DeltaR cut or jet number cuts");
	vetoEvent;
      }

      // 7j55
      if(njet55>=7&&pass55DeltaR)
	_count_7j55->fill( 0.5, weight) ;
      // 8j55
      if(njet55>=8&&pass55DeltaR)
	_count_8j55->fill( 0.5, weight) ;
      // 6j80
      if(njet80>=6&&pass80DeltaR)
	_count_6j80->fill( 0.5, weight) ;
      // 7j80
      if(njet80>=7&&pass80DeltaR)
	_count_7j80->fill( 0.5, weight) ;

    }

    //@}

    void finalize() {}

  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _etmissHTA;
    AIDA::IHistogram1D* _etmissHTB;
    AIDA::IHistogram1D* _njet55A;
    AIDA::IHistogram1D* _njet55B;
    AIDA::IHistogram1D* _njet80A;
    AIDA::IHistogram1D* _njet80B;
    AIDA::IHistogram1D* _count_7j55;
    AIDA::IHistogram1D* _count_8j55;
    AIDA::IHistogram1D* _count_6j80;
    AIDA::IHistogram1D* _count_7j80;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_S9225137);

}
