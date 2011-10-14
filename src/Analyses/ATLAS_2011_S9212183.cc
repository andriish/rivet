// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

    /// @author Chris Wymant
    class ATLAS_2011_S9212183 : public Analysis {
    public:

        /// @name Constructors etc.
        //@{

        /// Constructor
        ATLAS_2011_S9212183()
            : Analysis("ATLAS_2011_S9212183")
        {
            /// Set whether your finalize method needs the generator cross section
            setNeedsCrossSection(false);
        }

    //@}


    public:

        /// @name Analysis methods
        //@{

        /// Book histograms and initialise projections before the run
        void init() {

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

            /// Jet finder
            addProjection(FastJets(FinalState(), FastJets::ANTIKT, 0.4), 
            "AntiKtJets04");

            // all tracks (to do deltaR with leptons)
            addProjection(ChargedFinalState(-3.0,3.0),"cfs");

            // used for pTmiss (N.B. the real 'vfs' extends beyond 4.5 to |eta| = 4.9)
            addProjection(VisibleFinalState(-4.5,4.5),"vfs");


            /// Book histograms
            _count_2j =   bookHistogram1D("count_2j", 1, 0., 1.);
            _count_3j =   bookHistogram1D("count_3j", 1, 0., 1.);
            _count_4j5 =  bookHistogram1D("count_4j5", 1, 0., 1.);
            _count_4j10 = bookHistogram1D("count_4j10", 1, 0., 1.);
            _count_HM =   bookHistogram1D("count_HM", 1, 0., 1.);

            _hist_meff_2j  = bookHistogram1D("m_eff_2j", 30, 0., 3000.);
            _hist_meff_3j = bookHistogram1D("m_eff_3j", 30, 0., 3000.);
            _hist_meff_4j = bookHistogram1D("m_eff_4j", 30, 0., 3000.);
            _hist_meff_HM = bookHistogram1D("m_eff_HM", 20, 0., 3000.);

            _hist_eTmiss  = bookHistogram1D("Et_miss", 20, 0., 1000.);
        }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
        const double weight = event.weight();

        Jets cand_jets;
        foreach (const Jet& jet, 
            applyProjection<FastJets>(event, "AntiKtJets04").jetsByPt(20.0*GeV) ) {
            if ( fabs( jet.momentum().eta() ) < 4.9 ) {
                cand_jets.push_back(jet);
            }
        } 

        ParticleVector cand_e  = applyProjection<IdentifiedFinalState>(event, "elecs").particlesByPt();
      
        // muon isolation not mentioned in hep-exp 1109.6572, unlike in 1102.5290,
        // but assumed to still be applicable
        ParticleVector cand_mu;
        ParticleVector chg_tracks = applyProjection<ChargedFinalState>(event, "cfs").particles();
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

        // resolve jet-lepton overlap for jets with |eta| < 2.8 
        Jets cand_jets_2;
        foreach ( const Jet& jet, cand_jets ) {
            if ( fabs( jet.momentum().eta() ) >= 2.8 )
                cand_jets_2.push_back( jet );
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


        ParticleVector recon_e, recon_mu;

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
        // based on all candidate electrons, muons and jets, plus everything else with |eta| < 4.5
        // i.e. everything in our projection "vfs" plus the jets with |eta| > 4.5
        ParticleVector vfs_particles = applyProjection<VisibleFinalState>(event, "vfs").particles();
        FourMomentum pTmiss;
        foreach ( const Particle & p, vfs_particles ) {
            pTmiss -= p.momentum();
        }
        foreach ( const Jet& jet, cand_jets_2 ) {
            if ( fabs( jet.momentum().eta() ) > 4.5 )
                pTmiss -= jet.momentum();
        }
        double eTmiss = pTmiss.pT();


        // final jet filter
        Jets recon_jets;
        foreach ( const Jet& jet, cand_jets_2 ) {
            if ( fabs( jet.momentum().eta() ) <= 2.8 )
                recon_jets.push_back( jet );
        }
        // NB it seems that jets with |eta| > 2.8 could have been thrown away at
        // the start; we don't do so, in order to follow both the structure of 
        // the paper and the similar Rivet analysis ATLAS_2011_S8983313

        // 'candidate' muons needed only 10 GeV, to cause a veto they need 20 GeV
        ParticleVector veto_mu;
        foreach ( const Particle & mu, cand_mu ) {
            if ( mu.momentum().pT() >= 20.0*GeV )
                veto_mu.push_back(mu);
        }

        if ( ! ( veto_mu.empty() && recon_e.empty() ) ) {
            MSG_DEBUG("Charged leptons left after selection");
            vetoEvent;
        }

        if ( eTmiss <= 130 * GeV ) {
            MSG_DEBUG("Not enough eTmiss: " << eTmiss << " < 130");
            vetoEvent;
        }


        if ( recon_jets.empty() || recon_jets[0].momentum().pT() <= 130.0 * GeV ) {
            MSG_DEBUG("No hard leading jet in " << recon_jets.size() << " jets");
            vetoEvent;
        }
	
        // ==================== observables ====================

        int Njets = 0;
        double min_dPhi = 999.999;
        double pTmiss_phi = pTmiss.phi();
        foreach ( const Jet& jet, recon_jets ) {
            if ( jet.momentum().pT() > 40 * GeV ) {
                if ( Njets < 3 )
                    min_dPhi = min( min_dPhi, 
        deltaPhi( pTmiss_phi, jet.momentum().phi() ) );
                ++Njets;
            }
        }

        int NjetsHighMass = 0;
        foreach ( const Jet& jet, recon_jets ) {
            if ( jet.momentum().pT() > 80.0 * GeV ) {
                ++NjetsHighMass;
            }
        }

        if ( Njets < 2 ) {
            MSG_DEBUG("Only " << Njets << " >40 GeV jets left");
            vetoEvent;
        }

        if ( min_dPhi <= 0.4 ) {
            MSG_DEBUG("dPhi too small");
            vetoEvent;
        }

        // m_eff

        double m_eff_2j = eTmiss 
            + recon_jets[0].momentum().pT() 
            + recon_jets[1].momentum().pT();

        double m_eff_3j = recon_jets.size() < 3 ? -999.0 : 
            m_eff_2j + recon_jets[2].momentum().pT();

        double m_eff_4j = recon_jets.size() < 4 ? -999.0 : 
            m_eff_3j + recon_jets[3].momentum().pT();

        double m_eff_HM = eTmiss;
        foreach ( const Jet& jet, recon_jets ) {
            if ( jet.momentum().pT() > 40.0 * GeV )
                m_eff_HM += jet.momentum().pT();
        }

        double et_meff_2j = eTmiss / m_eff_2j;
        double et_meff_3j = eTmiss / m_eff_3j;
        double et_meff_4j = eTmiss / m_eff_4j;
        double et_meff_HM = eTmiss / m_eff_HM;

  
        // ==================== FILL ====================

        MSG_DEBUG( "Trying to fill " 
            << Njets << ' '
            << m_eff_2j << ' '
            << et_meff_2j << ' '
            << m_eff_3j << ' '
            << et_meff_3j << ' '
            << m_eff_4j << ' '
            << et_meff_4j << ' '
            << m_eff_HM << ' '
            << et_meff_HM );

        _hist_eTmiss->fill(eTmiss, weight);


        // 2j region
        if ( et_meff_2j > 0.3 ) {

            _hist_meff_2j->fill(m_eff_2j, weight);

            if ( m_eff_2j > 1000 * GeV ) {
                MSG_DEBUG("Hits 2j");
                _count_2j->fill(0.5, weight);
            }
        }


        // 3j region
        if ( Njets >= 3 && et_meff_3j > 0.25 ) {
	
            _hist_meff_3j->fill(m_eff_3j, weight);
	
            if ( m_eff_3j > 1000 * GeV ) {
                MSG_DEBUG("Hits 3j");
                _count_3j->fill(0.5, weight);
            }
        }


        // 4j5 & 4j10 regions
        if ( Njets >= 4 && et_meff_4j > 0.25 ) {
	
            _hist_meff_4j->fill(m_eff_4j, weight);
	
            if ( m_eff_4j > 500 * GeV ) {
                MSG_DEBUG("Hits 4j5");
                _count_4j5->fill(0.5, weight);
            }
            if ( m_eff_4j > 1000 * GeV ) {
                MSG_DEBUG("Hits 4j10");
                _count_4j10->fill(0.5, weight);
            }
        }


        // High Mass region
        if ( NjetsHighMass >= 4 && et_meff_HM > 0.2 ) {
	
            _hist_meff_HM->fill(m_eff_HM, weight);
	
            if ( m_eff_HM > 1100 * GeV ) {
                MSG_DEBUG("Hits HM");
                _count_HM->fill(0.5, weight);
            }
        }


    }

    //@}
    
    void finalize() {

        // Two, three and four jet channels have bin width = 100 (GeV)
        // High mass channel has bin width = 150 (GeV)
        // Integrated luminosity = 1040 (pb)
	
        scale( _hist_meff_2j, 100. * 1040 * crossSection()/sumOfWeights() );
        scale( _hist_meff_3j, 100. * 1040 * crossSection()/sumOfWeights() );
        scale( _hist_meff_4j, 100. * 1040 * crossSection()/sumOfWeights() );
        scale( _hist_meff_HM, 150. * 1040 * crossSection()/sumOfWeights() );


    }

    private:

        /// @name Histograms
        //@{
        AIDA::IHistogram1D* _count_2j;
        AIDA::IHistogram1D* _count_3j;
        AIDA::IHistogram1D* _count_4j5;
        AIDA::IHistogram1D* _count_4j10;
        AIDA::IHistogram1D* _count_HM;
        AIDA::IHistogram1D* _hist_meff_2j;
        AIDA::IHistogram1D* _hist_meff_3j;
        AIDA::IHistogram1D* _hist_meff_4j;
        AIDA::IHistogram1D* _hist_meff_HM;
        AIDA::IHistogram1D* _hist_eTmiss;
        //@}


    };


    // This global object acts as a hook for the plugin system
    AnalysisBuilder<ATLAS_2011_S9212183> plugin_ATLAS_2011_S9212183;

}
