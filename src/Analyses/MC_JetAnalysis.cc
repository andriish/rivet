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
      _h_pT_jet(njet, NULL),
      _h_eta_jet(njet, NULL), _h_eta_jet_plus(njet), _h_eta_jet_minus(njet),
      _h_rap_jet(njet, NULL), _h_rap_jet_plus(njet), _h_rap_jet_minus(njet),
      _h_mass_jet(njet, NULL)
  {
    setNeedsCrossSection(true); // legitimate use, since a base class has no .info file!
  }



  // Book histograms
  void MC_JetAnalysis::init() {

    for (size_t i=0; i < m_njet; ++i) {
      stringstream pTname;
      pTname << "jet_pT_" << i+1;
      double pTmax = 1.0/(double(i)+2.0) * sqrtS()/GeV/2.0;
      int nbins_pT = 100/(i+1);
      _h_pT_jet[i] = bookHistogram1D(pTname.str(), logBinEdges(nbins_pT, 10.0, pTmax));

      stringstream massname;
      massname << "jet_mass_" << i+1;
      double mmax = 100.0;
      int nbins_m = 100/(i+1);
      _h_mass_jet[i] = bookHistogram1D(massname.str(), logBinEdges(nbins_m, 1.0, mmax));

      stringstream etaname;
      etaname << "jet_eta_" << i+1;
      _h_eta_jet[i] = bookHistogram1D(etaname.str(), i > 1 ? 25 : 50, -5.0, 5.0);
      _h_eta_jet_plus[i].reset(new LWH::Histogram1D(i > 1 ? 15 : 25, 0, 5));
      _h_eta_jet_minus[i].reset(new LWH::Histogram1D(i > 1 ? 15 : 25, 0, 5));

      stringstream rapname;
      rapname << "jet_y_" << i+1;
      _h_rap_jet[i] = bookHistogram1D(rapname.str(), i>1 ? 25 : 50, -5.0, 5.0);
      _h_rap_jet_plus[i].reset(new LWH::Histogram1D(i > 1 ? 15 : 25, 0, 5));
      _h_rap_jet_minus[i].reset(new LWH::Histogram1D(i > 1 ? 15 : 25, 0, 5));

      for (size_t j = i+1; j < min(size_t(3),m_njet); ++j) {
        std::pair<size_t, size_t> ij = std::make_pair(i, j);

        stringstream detaname;
        detaname << "jets_deta_" << i+1 << j+1;
        _h_deta_jets.insert(make_pair(ij, bookHistogram1D(detaname.str(), 25, -5.0, 5.0)));

        stringstream dphiname;
        dphiname << "jets_dphi_" << i+1 << j+1;
        _h_dphi_jets.insert(make_pair(ij, bookHistogram1D(dphiname.str(), 25, 0.0, M_PI)));

        stringstream dRname;
        dRname << "jets_dR_" << i+1 << j+1;
        _h_dR_jets.insert(make_pair(ij, bookHistogram1D(dRname.str(), 25, 0.0, 5.0)));
      }
    }

    _h_jet_multi_exclusive = bookHistogram1D("jet_multi_exclusive", m_njet+3, -0.5, m_njet+3-0.5);
    _h_jet_multi_inclusive = bookHistogram1D("jet_multi_inclusive", m_njet+3, -0.5, m_njet+3-0.5);
    _h_jet_multi_ratio = bookDataPointSet("jet_multi_ratio", m_njet+2, 0.5, m_njet+3-0.5);
    _h_jet_HT = bookHistogram1D("jet_HT", logBinEdges(50, m_jetptcut, sqrtS()/GeV/2.0));
  }



  // Do the analysis
  void MC_JetAnalysis::analyze(const Event & e) {
    const double weight = e.weight();

    const Jets& jets = applyProjection<FastJets>(e, m_jetpro_name).jetsByPt(m_jetptcut);

    for (size_t i = 0; i < m_njet; ++i) {
      if (jets.size() < i+1) continue;
      _h_pT_jet[i]->fill(jets[i].momentum().pT()/GeV, weight);
      // Check for numerical precision issues with jet masses
      double m2_i = jets[i].momentum().mass2();
      if (m2_i < 0) {
        if (m2_i < -1e-4) {
          MSG_WARNING("Jet mass2 is negative: " << m2_i << " GeV^2. "
                      << "Truncating to 0.0, assuming numerical precision is to blame.");
        }
        m2_i = 0.0;
      }

      // Jet mass
      _h_mass_jet[i]->fill(sqrt(m2_i)/GeV, weight);

      // Jet eta
      const double eta_i = jets[i].momentum().eta();
      _h_eta_jet[i]->fill(eta_i, weight);
      if (eta_i > 0.0) {
        _h_eta_jet_plus[i]->fill(eta_i, weight);
      } else {
        _h_eta_jet_minus[i]->fill(fabs(eta_i), weight);
      }

      // Jet rapidity
      const double rap_i = jets[i].momentum().rapidity();
      _h_rap_jet[i]->fill(rap_i, weight);
      if (rap_i > 0.0) {
        _h_rap_jet_plus[i]->fill(rap_i, weight);
      } else {
        _h_rap_jet_minus[i]->fill(fabs(rap_i), weight);
      }

      // Inter-jet properties
      for (size_t j = i+1; j < min(size_t(3),m_njet); ++j) {
        if (jets.size() < j+1) continue;
        std::pair<size_t, size_t> ij = std::make_pair(i, j);
        double deta = jets[i].momentum().eta()-jets[j].momentum().eta();
        double dphi = deltaPhi(jets[i].momentum(),jets[j].momentum());
        double dR = deltaR(jets[i].momentum(), jets[j].momentum());
        _h_deta_jets[ij]->fill(deta, weight);
        _h_dphi_jets[ij]->fill(dphi, weight);
        _h_dR_jets[ij]->fill(dR, weight);
      }
    }
    _h_jet_multi_exclusive->fill(jets.size(), weight);

    for (size_t i = 0; i < m_njet+2; ++i) {
      if (jets.size() >= i) {
        _h_jet_multi_inclusive->fill(i, weight);
      }
    }

    double HT=0.0;
    foreach (const Jet& jet, jets) {
      HT += jet.momentum().pT();
    }
    _h_jet_HT->fill(HT, weight);
  }


  // Finalize
  void MC_JetAnalysis::finalize() {
    for (size_t i = 0; i < m_njet; ++i) {
      scale(_h_pT_jet[i], crossSection()/sumOfWeights());
      scale(_h_mass_jet[i], crossSection()/sumOfWeights());
      scale(_h_eta_jet[i], crossSection()/sumOfWeights());
      scale(_h_rap_jet[i], crossSection()/sumOfWeights());

      // Create eta/rapidity ratio plots
      stringstream etaname;
      etaname << "jet_eta_pmratio_" << i+1;
      histogramFactory().divide(histoPath(etaname.str()), *_h_eta_jet_plus[i], *_h_eta_jet_minus[i]);
      stringstream rapname;
      rapname << "jet_y_pmratio_" << i+1;
      histogramFactory().divide(histoPath(rapname.str()), *_h_rap_jet_plus[i], *_h_rap_jet_minus[i]);
    }

    // Scale the d{eta,R} histograms
    map<pair<size_t, size_t>, AIDA::IHistogram1D*>::iterator it;
    for (it = _h_deta_jets.begin(); it != _h_deta_jets.end(); ++it) {
      scale(it->second, crossSection()/sumOfWeights());
    }
    for (it = _h_dphi_jets.begin(); it != _h_dphi_jets.end(); ++it) {
      scale(it->second, crossSection()/sumOfWeights());
    }
    for (it = _h_dR_jets.begin(); it != _h_dR_jets.end(); ++it) {
      scale(it->second, crossSection()/sumOfWeights());
    }

    // Fill inclusive jet multi ratio
    int Nbins = _h_jet_multi_inclusive->axis().bins();
    std::vector<double> ratio(Nbins-1, 0.0);
    std::vector<double> err(Nbins-1, 0.0);
    for (int i = 0; i < Nbins-1; ++i) {
      if (_h_jet_multi_inclusive->binHeight(i) > 0.0 && _h_jet_multi_inclusive->binHeight(i+1) > 0.0) {
        ratio[i] = _h_jet_multi_inclusive->binHeight(i+1)/_h_jet_multi_inclusive->binHeight(i);
        double relerr_i = _h_jet_multi_inclusive->binError(i)/_h_jet_multi_inclusive->binHeight(i);
        double relerr_j = _h_jet_multi_inclusive->binError(i+1)/_h_jet_multi_inclusive->binHeight(i+1);
        err[i] = ratio[i] * (relerr_i + relerr_j);
      }
    }
    _h_jet_multi_ratio->setCoordinate(1, ratio, err);

    scale(_h_jet_multi_exclusive, crossSection()/sumOfWeights());
    scale(_h_jet_multi_inclusive, crossSection()/sumOfWeights());
    scale(_h_jet_HT, crossSection()/sumOfWeights());
  }


}
