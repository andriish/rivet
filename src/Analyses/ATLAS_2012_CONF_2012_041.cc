// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/RivetYODA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2012_CONF_2012_041 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor

    ATLAS_2012_CONF_2012_041()
      : Analysis("ATLAS_2012_CONF_2012_041")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialize projections before the run
    void init() {

      // projection to find the electrons
      std::vector<std::pair<double, double> > eta_e;
      eta_e.push_back(make_pair(-2.47,2.47));
      IdentifiedFinalState elecs(eta_e, 7.0*GeV);
      elecs.acceptIdPair(ELECTRON);
      addProjection(elecs, "elecs");

      // projection to find the muons
      std::vector<std::pair<double, double> > eta_m;
      eta_m.push_back(make_pair(-2.4,2.4));
      IdentifiedFinalState muons(eta_m, 6.0*GeV);
      muons.acceptIdPair(MUON);
      addProjection(muons, "muons");

      // Jet finder
      VetoedFinalState vfs;
      vfs.addVetoPairId(MUON);
      addProjection(FastJets(vfs, FastJets::ANTIKT, 0.4),
                   "AntiKtJets04");

      // all tracks (to do deltaR with leptons)
      addProjection(ChargedFinalState(-3.0,3.0,0.5*GeV),"cfs");

      // for pTmiss
      addProjection(VisibleFinalState(-4.9,4.9),"vfs");

      // Book histograms
      _count_3jet_channel = bookHisto1D("count_3jet_channel", 1, 0., 1.);
      _count_4jet_channel = bookHisto1D("count_4jet_channel", 1, 0., 1.);
      _count_soft_channel = bookHisto1D("count_soft_channel", 1, 0., 1.);

      _hist_m_eff_3jet        = bookHisto1D("hist_m_eff_3jet"       ,  6, 400., 1600.);
      _hist_m_eff_4jet        = bookHisto1D("hist_m_eff_4jet"       ,  6, 400., 1600.);
      _hist_eTmiss_m_eff_soft = bookHisto1D("hist_eTmiss_m_eff_soft",  6, 0.1 , 0.7  );

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // get the candiate jets
      Jets cand_jets;
      foreach ( const Jet& jet,
                applyProjection<FastJets>(event, "AntiKtJets04").jetsByPt(20.0*GeV) ) {
        if ( fabs( jet.momentum().eta() ) < 4.5 ) {
          cand_jets.push_back(jet);
        }
      }

      // get the candidate "medium" leptons without isolation
      ParticleVector cand_soft_e,cand_hard_e;
      foreach( const Particle & e,
               applyProjection<IdentifiedFinalState>(event, "elecs").particlesByPt()) {
        double pT  = e.momentum().perp();
        double eta = e.momentum().eta();
        // remove any leptons within 0.4 of any candidate jets
        bool e_near_jet = false;
        foreach ( const Jet& jet, cand_jets ) {
          double dR = deltaR(e.momentum(),jet.momentum());
          if ( dR < 0.4 && dR > 0.2 ) {
            e_near_jet = true;
            break;
          }
        }
        if ( e_near_jet ) continue;
        // soft selection
        if(pT>7.&&!(fabs(eta)>1.37&&fabs(eta)<1.52)) {
          cand_soft_e.push_back(e);
        }
        // hard selection
        if(pT>10.) cand_hard_e.push_back(e);
      }
      ParticleVector cand_soft_mu,cand_hard_mu;
      foreach( const Particle & mu,
               applyProjection<IdentifiedFinalState>(event, "muons").particlesByPt()) {
        double pT  = mu.momentum().perp();
        double eta = mu.momentum().eta();
        // remove any leptons within 0.4 of any candidate jets
        bool mu_near_jet = false;
        foreach ( const Jet& jet, cand_jets ) {
          if ( deltaR(mu.momentum(),jet.momentum()) < 0.4 ) {
            mu_near_jet = true;
            break;
          }
        }
        if ( mu_near_jet ) continue;
        // soft selection
        if(pT>6.&&!(fabs(eta)>1.37&&fabs(eta)<1.52)) {
          cand_soft_mu.push_back(mu);
        }
        // hard selection
        if(pT>10.) cand_hard_mu.push_back(mu);
      }
      // apply the isolation
      ParticleVector chg_tracks =
        applyProjection<ChargedFinalState>(event, "cfs").particles();
      // pTcone around muon track (hard)
      ParticleVector recon_hard_mu;
      foreach ( const Particle & mu, cand_hard_mu ) {
        double pTinCone = -mu.momentum().pT();
        if(-pTinCone<20.) continue;
        foreach ( const Particle & track, chg_tracks ) {
          if ( deltaR(mu.momentum(),track.momentum()) < 0.2 )
            pTinCone += track.momentum().pT();
        }
        if ( pTinCone < 1.8*GeV ) recon_hard_mu.push_back(mu);
      }
      // pTcone around muon track (soft)
      ParticleVector recon_soft_mu;
      foreach ( const Particle & mu, cand_soft_mu ) {
        double pTinCone = -mu.momentum().pT();
        if(-pTinCone>20.) continue;
        foreach ( const Particle & track, chg_tracks ) {
          if ( deltaR(mu.momentum(),track.momentum()) < 0.2 )
            pTinCone += track.momentum().pT();
        }
        if ( pTinCone < 1.8*GeV ) recon_soft_mu.push_back(mu);
      }
      // pTcone around electron track (hard)
      ParticleVector recon_hard_e;
      foreach ( const Particle & e, cand_hard_e ) {
        double pTinCone = -e.momentum().pT();
        if(-pTinCone<25.) continue;
        foreach ( const Particle & track, chg_tracks ) {
          if ( deltaR(e.momentum(),track.momentum()) < 0.2 )
            pTinCone += track.momentum().pT();
        }
        if ( pTinCone < 0.1 * e.momentum().pT() ) recon_hard_e.push_back(e);
      }
      // pTcone around electron track (soft)
      ParticleVector recon_soft_e;
      foreach ( const Particle & e, cand_soft_e ) {
        double pTinCone = -e.momentum().pT();
        if(-pTinCone>25.) continue;
        foreach ( const Particle & track, chg_tracks ) {
          if ( deltaR(e.momentum(),track.momentum()) < 0.2 )
            pTinCone += track.momentum().pT();
        }
        if ( pTinCone < 0.1 * e.momentum().pT() ) recon_soft_e.push_back(e);
      }

      // discard jets that overlap with electrons
      Jets recon_jets;
      foreach ( const Jet& jet, cand_jets ) {
        if(fabs(jet.momentum().eta())>2.5||
           jet.momentum().perp()<25.) continue;
        bool away_from_e = true;
        foreach ( const Particle & e, cand_hard_e ) {
          if ( deltaR(e.momentum(),jet.momentum()) < 0.2 ) {
            away_from_e = false;
            break;
          }
        }
        if ( away_from_e ) recon_jets.push_back( jet );
      }

      // pTmiss
      FourMomentum pTmiss;
      foreach ( const Particle & p,
                applyProjection<VisibleFinalState>(event, "vfs").particles() ) {
        pTmiss -= p.momentum();
      }
      double eTmiss = pTmiss.pT();

      // both selections require at least 2 jets
      if(recon_jets.size()<2) vetoEvent;

      // start of meff calculation
      double HT=0.;
      foreach( const Jet & jet, recon_jets) {
        HT += jet.momentum().perp();
      }
      double m_eff_inc  = HT+eTmiss;

      // hard selection exactly one candidate
      // and 1 recon and at least 3 jets
      if( cand_hard_e.size()  +  cand_hard_mu.size() == 1 &&
          recon_hard_e.size() + recon_hard_mu.size() == 1 &&
          recon_jets.size() >= 3 ) {
        // get the lepton
        Particle lepton = recon_hard_e.empty() ?
          recon_hard_mu[0] : recon_hard_e[0];
        // lepton variables
        double pT = lepton.momentum().perp();
        double mT  = 2.*(pT*eTmiss -
                         lepton.momentum().x()*pTmiss.x() -
                         lepton.momentum().y()*pTmiss.y());
        mT = sqrt(mT);
        HT += pT;
        m_eff_inc += pT;
        double m_eff = pT+eTmiss+recon_jets[0].momentum().perp()+
          recon_jets[1].momentum().perp()+recon_jets[2].momentum().perp();
        // three jet selection
        if(recon_jets[0].momentum().perp()>100. &&
           (recon_jets.size() == 3 ||
            recon_jets[3].momentum().perp() < 80. ) &&
           mT>100. && eTmiss>250. && eTmiss/m_eff>0.3) {
          if(m_eff_inc>1200.) _count_3jet_channel->fill(0.5,weight);
          _hist_m_eff_3jet->fill(min(1599.,m_eff_inc),weight);
        }
        // four jet selecton
        if(recon_jets.size() >= 4) {
          m_eff += recon_jets[3].momentum().perp();
          if(recon_jets[3].momentum().perp() > 80.  &&
             mT>100. && eTmiss>250. && eTmiss/m_eff>0.2) {
            if(m_eff_inc>800.) _count_4jet_channel->fill(0.5,weight);
            _hist_m_eff_4jet->fill(min(1599.,m_eff_inc),weight);
          }
        }
      }

      // soft selection exactly one candidate
      // and 1 recon and 4 jets
      if( cand_soft_e.size()  +  cand_soft_mu.size() == 1 &&
          recon_soft_e.size() + recon_soft_mu.size() == 1 &&
          recon_jets.size() >= 2 &&
          recon_jets[0].momentum().perp()>130.) {
        // get the lepton
        Particle lepton = recon_soft_e.empty() ?
          recon_soft_mu[0] : recon_soft_e[0];
        // lepton variables
        double pT = lepton.momentum().perp();
        double mT  = 2.*(pT*eTmiss -
                         lepton.momentum().x()*pTmiss.x() -
                         lepton.momentum().y()*pTmiss.y());
        mT = sqrt(mT);
        HT += pT;
        m_eff_inc += pT;
        double m_eff = pT+eTmiss+recon_jets[0].momentum().perp()+
          recon_jets[1].momentum().perp();
        if (mT>100. && eTmiss>250.) {
          if( eTmiss/m_eff>0.3 ) _count_soft_channel->fill(0.5,weight);
          _hist_eTmiss_m_eff_soft->fill( eTmiss/m_eff_inc,weight);
        }
      }
    }
    //@}


    void finalize() {

      double norm = 4.7* crossSection()/sumOfWeights()/femtobarn;
      scale(_count_3jet_channel    ,norm);
      scale(_count_4jet_channel    ,norm);
      scale(_count_soft_channel    ,norm);
      scale(_hist_m_eff_3jet       ,200.*norm);
      scale(_hist_m_eff_4jet       ,200.*norm);
      scale(_hist_eTmiss_m_eff_soft,0.1*norm);

    }

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _count_3jet_channel;
    Histo1DPtr _count_4jet_channel;
    Histo1DPtr _count_soft_channel;

    Histo1DPtr _hist_m_eff_3jet;
    Histo1DPtr _hist_m_eff_4jet;
    Histo1DPtr _hist_eTmiss_m_eff_soft;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_CONF_2012_041);

}
