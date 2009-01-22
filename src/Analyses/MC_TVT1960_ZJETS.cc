// -*- C++ -*-
#include "Rivet/Analyses/MC_TVT1960_ZJETS.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  MC_TVT1960_ZJETS::MC_TVT1960_ZJETS()
  {
    setBeams(PROTON, ANTIPROTON);
    
    //full final state
    FinalState fs(-5.0, 5.0);
    addProjection(fs, "FS");

    // leading leptons (for Z candidates)
    IdentifiedFinalState lfs(-5.0, 5.0, 25.0*GeV);
    lfs.acceptIdPair(ELECTRON);
    addProjection(lfs, "Leptons");
  } 



  // Book histograms
  void MC_TVT1960_ZJETS::init() {
    _h_Z_mass = bookHistogram1D("Z_mass", "$m_{\\mathrm{Z}}$", 50, 66.0, 116.0);
    _h_jet1_pT = bookHistogram1D("jet1_pT", "$p_{\\perp}^{\\mathrm{1st jet}}$", 50, 0.0, 500.0);
    _h_jet2_pT = bookHistogram1D("jet2_pT", "$p_{\\perp}^{\\mathrm{2nd jet}}$", 30, 0.0, 300.0);
    _h_jet3_pT = bookHistogram1D("jet3_pT", "$p_{\\perp}^{\\mathrm{3rd jet}}$", 20, 0.0, 200.0);
    _h_jet4_pT = bookHistogram1D("jet4_pT", "$p_{\\perp}^{\\mathrm{4th jet}}$", 10, 0.0, 100.0);
    _h_deta_Z_jet1 = bookHistogram1D("deta_Z_jet2", "$|\\Delta{\\eta}(\\mathrm{Z, 1st jet})|$", 20, 0.0, 5.0);
    _h_dR_jet2_jet3 = bookHistogram1D("dR_jet2_jet3", "$|\\Delta{R}(\\mathrm{2nd jet, 3rd jet})|$", 20, 0.0, 5.0);
    for (size_t i=0; i<4; ++i) {
      stringstream name, title;
      name<<"log10_d_"<<i<<i+1;
      title<<"$\\log_{10}$($k_\\perp$ jet resolution $"<<i<<" \\to "<<i+1<<"$ [GeV])";
      _h_log10_d[i] = bookHistogram1D(name.str(), title.str(), 100, 0.2, 2.6);
    }
    for (size_t i=0; i<5; ++i) {
      stringstream name, title;
      name<<"log10_R_"<<i;
      title<<"$\\log_{10}$(Integrated $"<<i<<"$ jet rate in $k_\\perp$ [GeV])";
      _h_log10_R[i] = bookHistogram1D(name.str(), title.str(), 100, 0.2, 2.6);
    }
  }



  // Do the analysis 
  void MC_TVT1960_ZJETS::analyze(const Event & event) {
    double weight = event.weight();

    // Skip if the event is empty
    const FinalState& fs = applyProjection<FinalState>(event, "FS");
    if (fs.isEmpty()) {
      getLog() << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
               << " because no final state pair found " << endl;
      vetoEvent(event);
    }
    
    // Find the Z candidates
    const FinalState & lfs = applyProjection<FinalState>(event, "Leptons");
    std::vector<std::pair<Particle, Particle> > Z_candidates;
    ParticleVector all_leptons=lfs.particles();
    for (size_t i=0; i<all_leptons.size(); ++i) {
      for (size_t j=i+1; j<all_leptons.size(); ++j) {
        double mZ=FourMomentum(all_leptons[i].momentum()+all_leptons[j].momentum()).mass()/GeV;
        if (mZ>66.0 && mZ<116.0) Z_candidates.push_back(make_pair(all_leptons[i], all_leptons[j]));
      }
    }
    if (Z_candidates.size() != 1) {
      getLog() << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
               << " because no unique lepton pair found " << endl;
      vetoEvent(event);
    }

    // Now build the jets on a FS without the electrons from the Z and their QED radiation
    ParticleVector jetparts;
    foreach (const Particle& p, fs.particles()) {
      bool copy = true;
      if (p.pdgId() == PHOTON) {
        FourMomentum p_e0=Z_candidates[0].first.momentum();
        FourMomentum p_e1=Z_candidates[0].second.momentum();
        FourMomentum p_P=p.momentum();
        if (deltaR(p_e0.pseudorapidity(), p_e0.azimuthalAngle(),
                   p_P.pseudorapidity(), p_P.azimuthalAngle()) < 0.2) {
            copy = false;
            Z_candidates[0].first.momentum()+=p_P;
        }
        if (deltaR(p_e1.pseudorapidity(), p_e1.azimuthalAngle(),
                   p_P.pseudorapidity(), p_P.azimuthalAngle()) < 0.2) {
            copy = false;
            Z_candidates[0].second.momentum()+=p_P;
        }
      }
      else {
        if (p.genParticle().barcode()==Z_candidates[0].first.genParticle().barcode()) {
          copy = false;
        }
        if (p.genParticle().barcode()==Z_candidates[0].second.genParticle().barcode()) {
          copy = false;
        }
      }
      if (copy) jetparts.push_back(p);
    }

    FastJets jetpro(fs); // fs only as dummy here
    jetpro.calc(jetparts);

    // Take jets with pt > 20
    /// @todo Make this neater, using the JetAlg interface and the built-in sorting
    const Jets& jets = jetpro.jets();
    Jets jets_cut;
    foreach (const Jet& j, jets) {
      if (j.momentum().pT()/GeV > 20.0) {
        jets_cut.push_back(j);
      }
    }
    getLog() << Log::DEBUG << "Num jets above 20 GeV = " << jets_cut.size() << endl;

    // Sort by pT:
    sort(jets_cut.begin(), jets_cut.end(), cmpJetsByPt);

    // cut on Delta R between jet and electrons
    foreach (const Jet& j, jets_cut) {
      Particle el=Z_candidates[0].first;
      if (deltaR(el.momentum().pseudorapidity(), el.momentum().azimuthalAngle(),
                 j.momentum().pseudorapidity(), j.momentum().azimuthalAngle()) < 0.5) {
        vetoEvent(event);
      }
      el=Z_candidates[0].second;
      if (deltaR(el.momentum().pseudorapidity(), el.momentum().azimuthalAngle(),
                 j.momentum().pseudorapidity(), j.momentum().azimuthalAngle()) < 0.5) {
        vetoEvent(event);
      }
    }

    FourMomentum zmom(Z_candidates[0].first.momentum()+Z_candidates[0].second.momentum());
    _h_Z_mass->fill(zmom.mass());
    if (jets_cut.size()>0) {
      _h_jet1_pT->fill(jets_cut[0].momentum().pT(), weight);
      double deta=fabs(zmom.pseudorapidity()-jets_cut[0].momentum().pseudorapidity());
      _h_deta_Z_jet1->fill(deta, weight);
    }
    if (jets_cut.size()>1) {
      _h_jet2_pT->fill(jets_cut[1].momentum().pT(), weight);
    }
    if (jets_cut.size()>2) {
      _h_jet3_pT->fill(jets_cut[2].momentum().pT(), weight);
      double dR23=deltaR(jets_cut[1].momentum().pseudorapidity(), jets_cut[1].momentum().azimuthalAngle(),
                         jets_cut[2].momentum().pseudorapidity(), jets_cut[2].momentum().azimuthalAngle());
      _h_dR_jet2_jet3->fill(dR23, weight);
    }
    if (jets_cut.size()>3) {
      _h_jet4_pT->fill(jets_cut[3].momentum().pT(), weight);
    }

    // jet resolutions
    const fastjet::ClusterSequence* seq = jetpro.clusterSeq();
    if (seq!=NULL) {
      for (size_t i=0; i<4; ++i) {
        double d_min=seq->exclusive_dmerge_max(i);
        if (d_min>0.0) _h_log10_d[i]->fill(log10(sqrt(d_min)), weight);
      }
    }
  }


  // Finalize
  void MC_TVT1960_ZJETS::finalize() {
    // integrated jet rates
    /// @todo: can probably be done much smarter using DPS or yoda?
    const AIDA::IAxis& axis(_h_log10_d[0]->axis());
    const int N=axis.bins();
    // R_0 (d) = 1 - int_d^dmax d_{01}
    double fullint0=integral(_h_log10_d[0]);
    for (int ibin=0; ibin<N; ++ibin) {
      double integral=0.0;
      for (int jbin=ibin; jbin<N; ++jbin) {
        integral+=_h_log10_d[0]->binHeight(jbin)*_h_log10_d[0]->axis().binWidth(jbin);
      }
      double binmean=axis.binLowerEdge(ibin)+0.5*axis.binWidth(ibin);
      _h_log10_R[0]->fill(binmean, fullint0-integral);
    }
    // R_n (d) = int_d^dmax d_{n-1,n} - d_{n,n+1}
    for (int n=1; n<4; ++n) {
      for (int ibin=0; ibin<N; ++ibin) {
        double integral=0.0;
        for (int jbin=ibin; jbin<N; ++jbin) {
          integral+=_h_log10_d[n-1]->binHeight(jbin)*_h_log10_d[n-1]->axis().binWidth(jbin)-
            _h_log10_d[n]->binHeight(jbin)*_h_log10_d[n]->axis().binWidth(jbin);
        }
        double binmean=axis.binLowerEdge(ibin)+0.5*axis.binWidth(ibin);
        _h_log10_R[n]->fill(binmean, integral);
      }
    }
    // R_4 (d) = int_d^dmax d_{3,4}
    for (int ibin=0; ibin<N; ++ibin) {
      double integral=0.0;
      for (int jbin=ibin; jbin<N; ++jbin) {
        integral+=_h_log10_d[3]->binHeight(jbin)*_h_log10_d[3]->axis().binWidth(jbin);
      }
      double binmean=axis.binLowerEdge(ibin)+0.5*axis.binWidth(ibin);
      _h_log10_R[4]->fill(binmean, integral);
    }

    normalize(_h_Z_mass,1.0);
    normalize(_h_jet1_pT,1.0);
    normalize(_h_jet2_pT,1.0);
    normalize(_h_jet3_pT,1.0);
    normalize(_h_jet4_pT,1.0);
    normalize(_h_deta_Z_jet1,1.0);
    normalize(_h_dR_jet2_jet3,1.0);
    for (size_t i=0; i<4; ++i) {
      normalize(_h_log10_d[i],1.0);
    }
  }

}
