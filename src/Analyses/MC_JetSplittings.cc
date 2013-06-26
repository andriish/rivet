// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  MC_JetSplittings::MC_JetSplittings(const string& name,
                                     size_t njet,
                                     const string& jetpro_name)
    : Analysis(name), m_njet(njet), m_jetpro_name(jetpro_name),
      _h_log10_d(njet), _h_log10_R(njet+1)
  {
    setNeedsCrossSection(true); // legitimate use, since a base class has no .info file!
  }



  // Book histograms
  void MC_JetSplittings::init() {
    for (size_t i = 0; i < m_njet; ++i) {
      string dname = "log10_d_" + to_str(i) + to_str(i+1);
      _h_log10_d[i] = bookHisto1D(dname, 100, 0.2, log10(0.5*sqrtS()));
      string Rname = "log10_R_" + to_str(i);
      _h_log10_R[i] = bookScatter2D(Rname, 50, 0.2, log10(0.5*sqrtS()));
    }
    string Rname = "log10_R_" + to_str(m_njet);
    _h_log10_R[m_njet] = bookScatter2D(Rname, 50, 0.2, log10(0.5*sqrtS()));
  }



  // Do the analysis
  void MC_JetSplittings::analyze(const Event & e) {
    const double weight = e.weight();

    const FastJets& jetpro = applyProjection<FastJets>(e, m_jetpro_name);

    // Jet resolutions and integrated jet rates
    const fastjet::ClusterSequence* seq = jetpro.clusterSeq();
    if (seq != NULL) {
      double previous_dij = 10.0;
      for (size_t i = 0; i < min(m_njet,(size_t)seq->n_particles()); ++i) {
        // Jet resolution i -> j
        double d_ij = log10(sqrt(seq->exclusive_dmerge_max(i)));

        // Fill differential jet resolution
        _h_log10_d[i]->fill(d_ij, weight);

        // Fill integrated jet resolution
        for (size_t ibin = 0; ibin < _h_log10_R[i]->numPoints(); ++ibin) {
          Point2D & dp = _h_log10_R[i]->point(ibin);
          /// @todo Could inRange be used here (avoiding the temporary variable)?
          double dcut = dp.x();
          if (d_ij < dcut && previous_dij > dcut) {
            dp.setY(dp.y() + weight);
          }
        }
        previous_dij = d_ij;
      }
      // One remaining integrated jet resolution
      for (size_t ibin = 0; ibin<_h_log10_R[m_njet]->numPoints(); ++ibin) {
        Point2D & dp = _h_log10_R[m_njet]->point(ibin);
        double dcut = dp.x();
        if (previous_dij > dcut) {
          dp.setY(dp.y() + weight);
        }
      }
    }

  }


  // Finalize
  void MC_JetSplittings::finalize() {
    for (size_t i = 0; i < m_njet; ++i) {
      scale(_h_log10_d[i], crossSection()/sumOfWeights());
      for (size_t ibin = 0; ibin<_h_log10_R[i]->numPoints(); ++ibin) {
        Point2D & dp = _h_log10_R[i]->point(ibin);
        dp.setY(dp.y()*crossSection()/sumOfWeights());
      }
    }

    for (size_t ibin = 0; ibin < _h_log10_R[m_njet]->numPoints(); ++ibin) {
      Point2D & dp =_h_log10_R[m_njet]->point(ibin);
      dp.setY(dp.y()*crossSection()/sumOfWeights());
    }
  }


}
