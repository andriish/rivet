// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "LWH/Profile1D.h"
#include "LWH/Histogram1D.h"

namespace Rivet {

  class ATLAS_2010_CONF_2010_081 : public Analysis {
  public:

    ATLAS_2010_CONF_2010_081() : Analysis("ATLAS_2010_CONF_2010_081") {
      setNeedsCrossSection(false);
    }


    void init() {
      ChargedFinalState cfs(-2.5, 2.5, 0.5*GeV);
      addProjection(cfs, "CFS");

      // We need at least one track with pT > 1GeV
      ChargedFinalState cfslead(-2.5, 2.5, 1.0*GeV);
      addProjection(cfslead, "CFSlead");

      int offset;
      if (sqrtS()/GeV >= 890 && sqrtS()/GeV <= 910)
        offset=0;
      else if (sqrtS()/GeV >= 6990 && sqrtS()/GeV <= 7010)
        offset=10;
      else
        return;

      _h_nch_toward    = bookProfile1D(offset+1, 1, 1);
      _h_nch_trans     = bookProfile1D(offset+1, 2, 1);
      _h_nch_away      = bookProfile1D(offset+1, 3, 1);
      _h_sumpt_toward  = bookProfile1D(offset+2, 1, 1);
      _h_sumpt_trans   = bookProfile1D(offset+2, 2, 1);
      _h_sumpt_away    = bookProfile1D(offset+2, 3, 1);
      _h_avgpt_toward  = bookProfile1D(offset+3, 1, 1);
      _h_avgpt_trans   = bookProfile1D(offset+3, 2, 1);
      _h_avgpt_away    = bookProfile1D(offset+3, 3, 1);
      _h_ptnch_toward  = bookProfile1D(offset+4, 1, 1);
      _h_ptnch_trans   = bookProfile1D(offset+4, 2, 1);
      _h_ptnch_away    = bookProfile1D(offset+4, 3, 1);
      for (size_t i=0; i<4; i++) {
        _h_n_vs_dphi[i]  = bookProfile1D(offset+5, 1, i+1);
        _h_pt_vs_dphi[i] = bookProfile1D(offset+6, 1, i+1);
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      // Require at least one track in the event with pT >= 1 GeV
      const ChargedFinalState& cfslead = applyProjection<ChargedFinalState>(event, "CFSlead");
      if (cfslead.size() < 1) {
        vetoEvent;
      }

      // These are the charged particles (tracks) with pT > 500 MeV
      const ChargedFinalState& charged = applyProjection<ChargedFinalState>(event, "CFS");

      // Identify leading track and its phi and pT
      ParticleVector particles = charged.particlesByPt();
      Particle p_lead = particles[0];
      const double philead = p_lead.momentum().phi();
      const double pTlead  = p_lead.momentum().perp();

      size_t numToward(0),     numTrans(0),     numAway(0);
      double ptSumToward(0.0), ptSumTrans(0.0), ptSumAway(0.0);

      // Loop over particles
      foreach (const Particle& p, particles) {
        const double dPhi = deltaPhi(philead, p.momentum().phi());
        const double pT = p.momentum().pT();
        if (dPhi < PI/3.0) {
          ptSumToward += pT;
          ++numToward;
        }
        else if (dPhi < 2*PI/3.0) {
          ptSumTrans += pT;
          ++numTrans;
        }
        else {
          ptSumAway += pT;
          ++numAway;
        }
      }

      // Temporary histos that bin N_ch in dPhi
      std::vector<LWH::Histogram1D> h_n_dphi;
      std::vector<LWH::Histogram1D> h_pt_dphi;
      for (size_t i=0; i<4; i++) {
        h_n_dphi.push_back(LWH::Histogram1D(binEdges(5, 1, i+1)));
        h_pt_dphi.push_back(LWH::Histogram1D(binEdges(6, 1, i+1)));
      }

      // Iterate over charged particles again - this is neccessary since the <pT> vs Nch
      // histos don't fill the Nch of the whole event as it was in the CDF analyses but
      // the number of charged particles in the transverse, toward and away region, respectively.
      // And these numbers are calculated in the for loop above
      double ptmin[4];
      if (sqrtS()/GeV >= 890 && sqrtS()/GeV <= 910) {
        ptmin[0] = 1.0;
        ptmin[1] = 1.5;
        ptmin[2] = 2.0;
        ptmin[3] = 2.5;
      }
      else {
        ptmin[0] = 1.0;
        ptmin[1] = 2.0;
        ptmin[2] = 3.0;
        ptmin[3] = 5.0;
      }

      for (ParticleVector::const_iterator p = particles.begin(); p != particles.end(); ++p) {
        const double dPhi = deltaPhi(philead, p->momentum().phi());
        const double pT = p->momentum().pT();
        if (dPhi < PI/3.0) {
          _h_ptnch_toward->fill(numToward, pT/GeV, weight);
        }
        else if (dPhi < 2*PI/3.0) {
          _h_ptnch_trans->fill(numTrans, pT/GeV, weight);
        }
        else {
          _h_ptnch_away->fill(numAway, pT/GeV, weight);
        }

        // Fill temp histos to bin N_ch in dPhi
        if (p != particles.begin()) { // We don't want to fill all those zeros from the leading track
          for (size_t i=0; i<4; i++) {
            if (pTlead/GeV >= 1.0) {
              h_n_dphi[i].fill( dPhi, 1);
              h_n_dphi[i].fill(-dPhi, 1);
              h_pt_dphi[i].fill( dPhi, pT/GeV);
              h_pt_dphi[i].fill(-dPhi, pT/GeV);
            }
          }
        }
      }

      const double avgpt_toward = numToward>0?ptSumToward/numToward:0;
      const double avgpt_trans  = numTrans>0?ptSumTrans/numTrans:0;
      const double avgpt_away   = numAway>0?ptSumAway/numAway:0;
      _h_nch_toward   ->fill(pTlead/GeV, numToward/(2*2.5 * 2*PI/3.0),    weight);
      _h_nch_trans    ->fill(pTlead/GeV, numTrans/(2*2.5 * 2*PI/3.0),     weight);
      _h_nch_away     ->fill(pTlead/GeV, numAway/(2*2.5 * 2*PI/3.0),      weight);
      _h_sumpt_toward ->fill(pTlead/GeV, ptSumToward/(2*2.5 * 2*PI/3.0),  weight);
      _h_sumpt_trans  ->fill(pTlead/GeV, ptSumTrans/(2*2.5 * 2*PI/3.0),   weight);
      _h_sumpt_away   ->fill(pTlead/GeV, ptSumAway/(2*2.5 * 2*PI/3.0),    weight);
      _h_avgpt_toward ->fill(pTlead/GeV, avgpt_toward, weight);
      if (numTrans>0) _h_avgpt_trans  ->fill(pTlead/GeV, avgpt_trans,  weight);
      if (numAway>0)  _h_avgpt_away   ->fill(pTlead/GeV, avgpt_away,   weight);
      _h_ptnch_toward ->fill(numToward,  avgpt_toward, weight);
      _h_ptnch_trans  ->fill(numTrans,   avgpt_trans,  weight);
      _h_ptnch_away   ->fill(numAway,    avgpt_away,   weight);

      for (size_t i=0; i<4; i++) {
        if (pTlead/GeV >= ptmin[i]) { // @todo: Do we need particles.size()>=2?
          const size_t nbins = _h_n_vs_dphi[i]->axis().bins();
          // deta*dphi, the range is not full 2*PI, thus the "+2". Extra factor of 2 because plot is symmetrised.
          const double scale = 2 * 2*2.5 * 2*PI/(nbins+2);
          for (size_t n=0; n<nbins; n++) {
            _h_n_vs_dphi[i]->fill(h_n_dphi[i].binMean(n), h_n_dphi[i].binHeight(n)/scale, weight);
            _h_pt_vs_dphi[i]->fill(h_pt_dphi[i].binMean(n), h_pt_dphi[i].binHeight(n)/scale, weight);
          }
        }
      }
    }

    void finalize() {
    }


  private:

    AIDA::IProfile1D* _h_nch_toward;
    AIDA::IProfile1D* _h_nch_trans;
    AIDA::IProfile1D* _h_nch_away;
    AIDA::IProfile1D* _h_sumpt_toward;
    AIDA::IProfile1D* _h_sumpt_trans;
    AIDA::IProfile1D* _h_sumpt_away;
    AIDA::IProfile1D* _h_avgpt_toward;
    AIDA::IProfile1D* _h_avgpt_trans;
    AIDA::IProfile1D* _h_avgpt_away;
    AIDA::IProfile1D* _h_ptnch_toward;
    AIDA::IProfile1D* _h_ptnch_trans;
    AIDA::IProfile1D* _h_ptnch_away;
    AIDA::IProfile1D* _h_n_vs_dphi[4];
    AIDA::IProfile1D* _h_pt_vs_dphi[4];
  };

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<ATLAS_2010_CONF_2010_081> plugin_ATLAS_2010_CONF_2010_081;

}

