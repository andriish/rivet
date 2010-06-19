// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  MC_JetAnalysis::MC_JetAnalysis(const string& name,
                                 size_t njet,
                                 const string& jetpro_name,
                                 double jetptcut)
    : Analysis(name), m_njet(njet), m_jetpro_name(jetpro_name), m_jetptcut(jetptcut),
      _h_log10_d(njet, NULL), _h_log10_R(njet+1, NULL), _h_pT_jet(njet, NULL),
      _h_eta_jet(njet, NULL), _h_rap_jet(njet, NULL), _h_mass_jet(njet, NULL)
  {
    setNeedsCrossSection(true);
  }



  // Book histograms
  void MC_JetAnalysis::init() {

    for (size_t i=0; i<m_njet; ++i) {
      stringstream dname;
      dname<<"log10_d_"<<i<<i+1;

      _h_log10_d[i] = bookHistogram1D(dname.str(), 50, 0.2, log10(0.5*sqrtS()));

      stringstream Rname;
      Rname<<"log10_R_"<<i;
      _h_log10_R[i] = bookDataPointSet(Rname.str(), 50, 0.2, log10(0.5*sqrtS()));

      stringstream pTname;
      pTname<<"jet_pT_"<<i+1;
      double pTmax = 1.0/(double(i)+2.0)*sqrtS()/GeV/2.0;
      int nbins_pT = 100/(i+1);
      _h_pT_jet[i] = bookHistogram1D(pTname.str(), logBinEdges(nbins_pT, 10.0, pTmax));

      stringstream massname;
      massname<<"jet_mass_"<<i+1;
      double mmax = 100.0;
      int nbins_m = 100/(i+1);
      _h_mass_jet[i] = bookHistogram1D(massname.str(), logBinEdges(nbins_m, 1.0, mmax));

      stringstream etaname;
      etaname<<"jet_eta_"<<i+1;
      _h_eta_jet[i] = bookHistogram1D(etaname.str(), i>1 ? 25 : 50, -5.0, 5.0);

      stringstream rapname;
      rapname<<"jet_y_"<<i+1;
      _h_rap_jet[i] = bookHistogram1D(rapname.str(), i>1 ? 25 : 50, -5.0, 5.0);

      for (size_t j=i+1; j<m_njet; ++j) {
        std::pair<size_t, size_t> ij(std::make_pair(i, j));

        stringstream detaname;
        detaname<<"jets_deta_"<<i+1<<j+1;
        _h_deta_jets.insert(make_pair(ij, bookHistogram1D(detaname.str(), 25, -5.0, 5.0)));

        stringstream dRname;
        dRname<<"jets_dR_"<<i+1<<j+1;
        _h_dR_jets.insert(make_pair(ij, bookHistogram1D(dRname.str(), 25, 0.0, 5.0)));
      }
    }
    stringstream Rname;
    Rname<<"log10_R_"<<m_njet;
    _h_log10_R[m_njet] = bookDataPointSet(Rname.str(), 50, 0.2, log10(0.5*sqrtS()));

    _h_jet_multi_exclusive = bookHistogram1D("jet_multi_exclusive", m_njet+3, -0.5, m_njet+3-0.5);
    _h_jet_multi_inclusive = bookHistogram1D("jet_multi_inclusive", m_njet+3, -0.5, m_njet+3-0.5);
    _h_jet_multi_ratio = bookDataPointSet("jet_multi_ratio", m_njet+2, 0.5, m_njet+3-0.5);
  }



  // Do the analysis
  void MC_JetAnalysis::analyze(const Event & e) {
    const double weight = e.weight();

    const FastJets& jetpro = applyProjection<FastJets>(e, m_jetpro_name);

    // Jet resolutions and integrated jet rates
    const fastjet::ClusterSequence* seq = jetpro.clusterSeq();
    if (seq!=NULL) {
      double previous_dij = 10.0;
      for (size_t i=0; i<m_njet; ++i) {
        // Jet resolution i -> j
        double d_ij=log10(sqrt(seq->exclusive_dmerge_max(i)));

        // Fill differential jet resolution
        _h_log10_d[i]->fill(d_ij, weight);

        // Fill integrated jet resolution
        for (int ibin=0; ibin<_h_log10_R[i]->size(); ++ibin) {
          IDataPoint* dp=_h_log10_R[i]->point(ibin);
          double dcut=dp->coordinate(0)->value();
          if (d_ij<dcut && previous_dij>dcut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value()+weight);
          }
        }
        previous_dij = d_ij;
      }
      // One remaining integrated jet resolution
      for (int ibin=0; ibin<_h_log10_R[m_njet]->size(); ++ibin) {
        IDataPoint* dp=_h_log10_R[m_njet]->point(ibin);
        double dcut=dp->coordinate(0)->value();
        if (previous_dij>dcut) {
          dp->coordinate(1)->setValue(dp->coordinate(1)->value()+weight);
        }
      }
    }

    const Jets& jets = jetpro.jetsByPt(m_jetptcut);

    // The remaining direct jet observables
    for (size_t i=0; i<m_njet; ++i) {
      if (jets.size()<i+1) continue;
      _h_pT_jet[i]->fill(jets[i].momentum().pT()/GeV, weight);
      _h_mass_jet[i]->fill(jets[i].momentum().mass()/GeV, weight);
      _h_eta_jet[i]->fill(jets[i].momentum().eta(), weight);
      _h_rap_jet[i]->fill(jets[i].momentum().rapidity(), weight);
      // cout << "Jet mass [" << i+1 << "] = " << jets[i].momentum().mass()/GeV << " GeV" << endl;

      for (size_t j=i+1; j<m_njet; ++j) {
        if (jets.size()<j+1) continue;
        std::pair<size_t, size_t> ij(std::make_pair(i, j));
        double deta = jets[i].momentum().eta()-jets[j].momentum().eta();
        double dR = deltaR(jets[i].momentum(), jets[j].momentum());
        _h_deta_jets[ij]->fill(deta, weight);
        _h_dR_jets[ij]->fill(dR, weight);
      }
    }
    _h_jet_multi_exclusive->fill(jets.size(), weight);

    for (size_t i=0; i<m_njet+2; ++i) {
      if (jets.size()>=i) {
        _h_jet_multi_inclusive->fill(i, weight);
      }
    }
  }


  // Finalize
  void MC_JetAnalysis::finalize() {
    for (size_t i=0; i<m_njet; ++i) {
      scale(_h_log10_d[i], crossSection()/sumOfWeights());
      for (int ibin=0; ibin<_h_log10_R[i]->size(); ++ibin) {
        IDataPoint* dp=_h_log10_R[i]->point(ibin);
        dp->coordinate(1)->setValue(dp->coordinate(1)->value()*crossSection()/sumOfWeights());
      }

      scale(_h_pT_jet[i], crossSection()/sumOfWeights());
      scale(_h_mass_jet[i], crossSection()/sumOfWeights());
      scale(_h_eta_jet[i], crossSection()/sumOfWeights());
      scale(_h_rap_jet[i], crossSection()/sumOfWeights());

    }
    for (int ibin=0; ibin<_h_log10_R[m_njet]->size(); ++ibin) {
      IDataPoint* dp=_h_log10_R[m_njet]->point(ibin);
      dp->coordinate(1)->setValue(dp->coordinate(1)->value()*crossSection()/sumOfWeights());
    }

    // Scale the d{eta,R} histograms
    map<pair<size_t, size_t>, AIDA::IHistogram1D*>::iterator it;
    for (it=_h_deta_jets.begin(); it!=_h_deta_jets.begin(); ++it) {
      scale(it->second, crossSection()/sumOfWeights());
    }
    for (it=_h_dR_jets.begin(); it!=_h_dR_jets.begin(); ++it) {
      scale(it->second, crossSection()/sumOfWeights());
    }

    // Fill inclusive jet multi ratio
    int Nbins=_h_jet_multi_inclusive->axis().bins();
    std::vector<double> ratio(Nbins-1, 0.0);
    std::vector<double> err(Nbins-1, 0.0);
    for (int i=0; i<Nbins-1; ++i) {
      if (_h_jet_multi_inclusive->binHeight(i)>0.0 && _h_jet_multi_inclusive->binHeight(i+1)>0.0) {
        ratio[i]=_h_jet_multi_inclusive->binHeight(i+1)/_h_jet_multi_inclusive->binHeight(i);
        double relerr_i=_h_jet_multi_inclusive->binError(i)/_h_jet_multi_inclusive->binHeight(i);
        double relerr_j=_h_jet_multi_inclusive->binError(i+1)/_h_jet_multi_inclusive->binHeight(i+1);
        err[i]=ratio[i]*(relerr_i+relerr_j);
      }
    }
    _h_jet_multi_ratio->setCoordinate(1, ratio, err);

    scale(_h_jet_multi_exclusive, crossSection()/sumOfWeights());
    scale(_h_jet_multi_inclusive, crossSection()/sumOfWeights());
  }


}
