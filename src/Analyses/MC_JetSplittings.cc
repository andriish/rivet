// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  MC_JetSplittings::MC_JetSplittings(const string& name,
                                     size_t njet,
                                     const string& jetpro_name)
    : Analysis(name), m_njet(njet), m_jetpro_name(jetpro_name),
      _h_log10_d(njet, NULL), _h_log10_R(njet+1, NULL)
  {
    setNeedsCrossSection(true); // legitimate use, since a base class has no .info file!
  }



  // Book histograms
  void MC_JetSplittings::init() {

    for (size_t i=0; i < m_njet; ++i) {
      stringstream dname;
      dname << "log10_d_" << i << i+1;

      _h_log10_d[i] = bookHistogram1D(dname.str(), 100, 0.2, log10(0.5*sqrtS()));

      stringstream Rname;
      Rname << "log10_R_" << i;
      _h_log10_R[i] = bookDataPointSet(Rname.str(), 100, 0.2, log10(0.5*sqrtS()));
    }
    stringstream Rname;
    Rname << "log10_R_" << m_njet;
    _h_log10_R[m_njet] = bookDataPointSet(Rname.str(), 100, 0.2, log10(0.5*sqrtS()));
  }



  // Do the analysis
  void MC_JetSplittings::analyze(const Event & e) {
    const double weight = e.weight();

    const FastJets& jetpro = applyProjection<FastJets>(e, m_jetpro_name);

    // Jet resolutions and integrated jet rates
    const fastjet::ClusterSequence* seq = jetpro.clusterSeq();
    if (seq != NULL) {
      double previous_dij = 10.0;
      for (size_t i = 0; i < m_njet; ++i) {
        // Jet resolution i -> j
        double d_ij = log10(sqrt(seq->exclusive_dmerge_max(i)));

        // Fill differential jet resolution
        _h_log10_d[i]->fill(d_ij, weight);

        // Fill integrated jet resolution
        for (int ibin = 0; ibin < _h_log10_R[i]->size(); ++ibin) {
          IDataPoint* dp = _h_log10_R[i]->point(ibin);
          double dcut = dp->coordinate(0)->value();
          if (d_ij < dcut && previous_dij > dcut) {
            dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
          }
        }
        previous_dij = d_ij;
      }
      // One remaining integrated jet resolution
      for (int ibin = 0; ibin<_h_log10_R[m_njet]->size(); ++ibin) {
        IDataPoint* dp = _h_log10_R[m_njet]->point(ibin);
        double dcut = dp->coordinate(0)->value();
        if (previous_dij > dcut) {
          dp->coordinate(1)->setValue(dp->coordinate(1)->value() + weight);
        }
      }
    }

  }


  // Finalize
  void MC_JetSplittings::finalize() {
    for (size_t i = 0; i < m_njet; ++i) {
      scale(_h_log10_d[i], crossSection()/sumOfWeights());
      for (int ibin = 0; ibin<_h_log10_R[i]->size(); ++ibin) {
        IDataPoint* dp = _h_log10_R[i]->point(ibin);
        dp->coordinate(1)->setValue(dp->coordinate(1)->value()*crossSection()/sumOfWeights());
      }
    }

    for (int ibin = 0; ibin < _h_log10_R[m_njet]->size(); ++ibin) {
      IDataPoint* dp =_h_log10_R[m_njet]->point(ibin);
      dp->coordinate(1)->setValue(dp->coordinate(1)->value()*crossSection()/sumOfWeights());
    }
  }


}
