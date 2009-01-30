// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/JADE_OPAL_2000_S4300807.hh"

namespace Rivet {

  void JADE_OPAL_2000_S4300807::analyze(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const double weight = e.weight();

    // are we running with a compatible CMS energy?
    const ParticlePair& beams = applyProjection<Beam>(e, "Beams").beams();
    double sbeams = beams.first.momentum().vector3().mod() + beams.second.momentum().vector3().mod();

    if (sbeams<_sqrts-0.5 || sbeams>_sqrts+0.5) {
      getLog() << Log::ERROR 
               << "CMS energy of events sqrt(s)=" << sbeams
               <<" doesn't match analysis energy sqrt(s)=" << _sqrts << endl;
      exit(1);
    }

    // Jets
    #ifdef HAVE_JADE
    getLog() << Log::DEBUG << "Using FastJet JADE patch to make diff jet rate plots:" << endl;

    const FastJets& jadejet = applyProjection<FastJets>(e, "JadeJets");
    if (jadejet.clusterSeq()) {
      double y_23=jadejet.clusterSeq()->exclusive_dmerge_max(2);
      double y_34=jadejet.clusterSeq()->exclusive_dmerge_max(3);
      double y_45=jadejet.clusterSeq()->exclusive_dmerge_max(4);
      double y_56=jadejet.clusterSeq()->exclusive_dmerge_max(5);

      const AIDA::IAxis& axis(_h_R_Jade[0]->axis());
      const int N=axis.bins();
      for (int ibin=axis.coordToIndex(y_23); ibin<N; ++ibin) {
        double binmean=axis.binLowerEdge(ibin)+0.5*axis.binWidth(ibin);
        if (y_23<binmean)
          _h_R_Jade[0]->fill(binmean, weight*axis.binWidth(ibin));
      }

      for (int ibin=axis.coordToIndex(y_34); ibin<axis.coordToIndex(y_23); ++ibin) {
        double binmean=axis.binLowerEdge(ibin)+0.5*axis.binWidth(ibin);
        if (y_34<binmean && y_23>binmean)
          _h_R_Jade[1]->fill(binmean, weight*axis.binWidth(ibin));
      }

      for (int ibin=axis.coordToIndex(y_45); ibin<axis.coordToIndex(y_34); ++ibin) {
        double binmean=axis.binLowerEdge(ibin)+0.5*axis.binWidth(ibin);
        if (y_45<binmean && y_34>binmean)
          _h_R_Jade[2]->fill(binmean, weight*axis.binWidth(ibin));
      }

      for (int ibin=axis.coordToIndex(y_56); ibin<axis.coordToIndex(y_45); ++ibin) {
        double binmean=axis.binLowerEdge(ibin)+0.5*axis.binWidth(ibin);
        if (y_56<binmean && y_45>binmean)
          _h_R_Jade[3]->fill(binmean, weight*axis.binWidth(ibin));
      }

      for (int ibin=0.0; ibin<axis.coordToIndex(y_56); ++ibin) {
        double binmean=axis.binLowerEdge(ibin)+0.5*axis.binWidth(ibin);
        if (y_56>binmean)
          _h_R_Jade[4]->fill(binmean, weight*axis.binWidth(ibin));
      }
    }

    const FastJets& durjet = applyProjection<FastJets>(e, "DurhamJets");
    if (durjet.clusterSeq()) {
      double y_23=durjet.clusterSeq()->exclusive_dmerge_max(2);
      double y_34=durjet.clusterSeq()->exclusive_dmerge_max(3);
      double y_45=durjet.clusterSeq()->exclusive_dmerge_max(4);
      double y_56=durjet.clusterSeq()->exclusive_dmerge_max(5);

      _h_y_Durham[0]->fill(y_23, weight);
      _h_y_Durham[1]->fill(y_34, weight);
      _h_y_Durham[2]->fill(y_45, weight);
      _h_y_Durham[3]->fill(y_56, weight);

      const AIDA::IAxis& axis(_h_R_Durham[0]->axis());
      const int N=axis.bins();
      for (int ibin=axis.coordToIndex(y_23); ibin<N; ++ibin) {
        double binmean=axis.binLowerEdge(ibin)+0.5*axis.binWidth(ibin);
        if (y_23<binmean)
          _h_R_Durham[0]->fill(binmean, weight*axis.binWidth(ibin));
      }

      for (int ibin=axis.coordToIndex(y_34); ibin<axis.coordToIndex(y_23); ++ibin) {
        double binmean=axis.binLowerEdge(ibin)+0.5*axis.binWidth(ibin);
        if (y_34<binmean && y_23>binmean)
          _h_R_Durham[1]->fill(binmean, weight*axis.binWidth(ibin));
      }

      for (int ibin=axis.coordToIndex(y_45); ibin<axis.coordToIndex(y_34); ++ibin) {
        double binmean=axis.binLowerEdge(ibin)+0.5*axis.binWidth(ibin);
        if (y_45<binmean && y_34>binmean)
          _h_R_Durham[2]->fill(binmean, weight*axis.binWidth(ibin));
      }

      for (int ibin=axis.coordToIndex(y_56); ibin<axis.coordToIndex(y_45); ++ibin) {
        double binmean=axis.binLowerEdge(ibin)+0.5*axis.binWidth(ibin);
        if (y_56<binmean && y_45>binmean)
          _h_R_Durham[3]->fill(binmean, weight*axis.binWidth(ibin));
      }

      for (int ibin=0.0; ibin<axis.coordToIndex(y_56); ++ibin) {
        double binmean=axis.binLowerEdge(ibin)+0.5*axis.binWidth(ibin);
        if (y_56>binmean)
          _h_R_Durham[4]->fill(binmean, weight*axis.binWidth(ibin));
      }
    }
    #endif
  }



  void JADE_OPAL_2000_S4300807::init() {
    #ifndef HAVE_JADE
    getLog() << Log::WARN << "No FastJet JADE patch, so not making any diff jet rate plots." << endl;
    #endif
    stringstream ss;
    ss<<_sqrts;
    _h_R_Jade[0]=bookHistogram1D(_nr_R_Jade, 1, 1, "Integrated 2-jet rate with Jade algorithm, $R_2^\\text{Jade}$ ("+ss.str()+" GeV)");
    _h_R_Jade[1]=bookHistogram1D(_nr_R_Jade, 1, 2, "Integrated 3-jet rate with Jade algorithm, $R_3^\\text{Jade}$ ("+ss.str()+" GeV)");
    _h_R_Jade[2]=bookHistogram1D(_nr_R_Jade, 1, 3, "Integrated 4-jet rate with Jade algorithm, $R_4^\\text{Jade}$ ("+ss.str()+" GeV)");
    _h_R_Jade[3]=bookHistogram1D(_nr_R_Jade, 1, 4, "Integrated 5-jet rate with Jade algorithm, $R_5^\\text{Jade}$ ("+ss.str()+" GeV)");
    _h_R_Jade[4]=bookHistogram1D(_nr_R_Jade, 1, 5, "Integrated $>$6-jet rate with Jade algorithm, $R_6^\\text{Jade}$ ("+ss.str()+" GeV)");

    _h_R_Durham[0]=bookHistogram1D(_nr_R_Durham, 1, 1, "Integrated 2-jet rate with Durham algorithm, $R_2^\\text{Durham}$ ("+ss.str()+" GeV)");
    _h_R_Durham[1]=bookHistogram1D(_nr_R_Durham, 1, 2, "Integrated 3-jet rate with Durham algorithm, $R_3^\\text{Durham}$ ("+ss.str()+" GeV)");
    _h_R_Durham[2]=bookHistogram1D(_nr_R_Durham, 1, 3, "Integrated 4-jet rate with Durham algorithm, $R_4^\\text{Durham}$ ("+ss.str()+" GeV)");
    _h_R_Durham[3]=bookHistogram1D(_nr_R_Durham, 1, 4, "Integrated 5-jet rate with Durham algorithm, $R_5^\\text{Durham}$ ("+ss.str()+" GeV)");
    _h_R_Durham[4]=bookHistogram1D(_nr_R_Durham, 1, 5, "Integrated $>$6-jet rate with Durham algorithm, $R_6^\\text{Durham}$ ("+ss.str()+" GeV)");

    _h_y_Durham[0]=bookHistogram1D(_nr_y_Durham, 1, 1, "Differential 2-jet rate with Durham algorithm, $y_{23}^\\text{Durham}$ ("+ss.str()+" GeV)");
    _h_y_Durham[1]=bookHistogram1D(_nr_y_Durham, 1, 2, "Differential 3-jet rate with Durham algorithm, $y_{34}^\\text{Durham}$ ("+ss.str()+" GeV)");
    _h_y_Durham[2]=bookHistogram1D(_nr_y_Durham, 1, 3, "Differential 4-jet rate with Durham algorithm, $y_{45}^\\text{Durham}$ ("+ss.str()+" GeV)");
    _h_y_Durham[3]=bookHistogram1D(_nr_y_Durham, 1, 4, "Differential 5-jet rate with Durham algorithm, $y_{56}^\\text{Durham}$ ("+ss.str()+" GeV)");
  }



  // Finalize
  void JADE_OPAL_2000_S4300807::finalize() {
    for (size_t n=0; n<4; ++n) {
      scale(_h_y_Durham[n], 1.0/sumOfWeights());
    }
    
    for (size_t n=0; n<5; ++n) {
      scale(_h_R_Durham[n], 100.0/sumOfWeights());
      scale(_h_R_Jade[n], 100.0/sumOfWeights());
    }
  }


}
