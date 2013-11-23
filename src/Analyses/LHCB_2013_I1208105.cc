// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "LWH/Profile1D.h"
#include "LWH/Histogram1D.h"
#include "LWH/AIMeasurement.h"

namespace Rivet {


  class LHCB_2013_I1208105 : public Analysis {
  public:

    LHCB_2013_I1208105()
      : Analysis("LHCB_2013_I1208105")
    {   }


    void init() {
      // Projections
      addProjection(FinalState(1.9, 4.9), "forwardFS");
      addProjection(FinalState(-3.5,-1.5), "backwardFS");
      addProjection(ChargedFinalState(1.9, 4.9), "forwardCFS");
      addProjection(ChargedFinalState(-3.5,-1.5), "backwardCFS");

      // Histos
      /// @todo Why manual booking?
      _dps_chEF_minbias = bookDataPointSet("d01-x01-y01", 10, 1.9, 4.9);
      _dps_chEF_hard = bookDataPointSet("d02-x01-y01", 10, 1.9, 4.9);
      _dps_chEF_diff = bookDataPointSet("d03-x01-y01", 10, 1.9, 4.9);
      _dps_chEF_nondiff = bookDataPointSet("d04-x01-y01", 10, 1.9, 4.9);
      _dps_totEF_minbias = bookDataPointSet("d05-x01-y01", 10, 1.9, 4.9);
      _dps_totEF_hard = bookDataPointSet("d06-x01-y01", 10, 1.9, 4.9);
      _dps_totEF_diff = bookDataPointSet("d07-x01-y01", 10, 1.9, 4.9);
      _dps_totEF_nondiff = bookDataPointSet("d08-x01-y01", 10, 1.9, 4.9);
      //
      _tp_chEF_minbias.reset(new LWH::Profile1D(10, 1.9, 4.9));
      _tp_chEF_hard.reset(new LWH::Profile1D(10, 1.9, 4.9));
      _tp_chEF_diff.reset(new LWH::Profile1D(10, 1.9, 4.9));
      _tp_chEF_nondiff.reset(new LWH::Profile1D(10, 1.9, 4.9));
      _tp_totEF_minbias.reset(new LWH::Profile1D(10, 1.9, 4.9));
      _tp_totEF_hard.reset(new LWH::Profile1D(10, 1.9, 4.9));
      _tp_totEF_diff.reset(new LWH::Profile1D(10, 1.9, 4.9));
      _tp_totEF_nondiff.reset(new LWH::Profile1D(10, 1.9, 4.9));
      //
      _h_chN_minbias.reset(new LWH::Histogram1D(10, 1.9, 4.9));
      _h_chN_hard.reset(new LWH::Histogram1D(10, 1.9, 4.9));
      _h_chN_diff.reset(new LWH::Histogram1D(10, 1.9, 4.9));
      _h_chN_nondiff.reset(new LWH::Histogram1D(10, 1.9, 4.9));
      _h_totN_minbias.reset(new LWH::Histogram1D(10, 1.9, 4.9));
      _h_totN_hard.reset(new LWH::Histogram1D(10, 1.9, 4.9));
      _h_totN_diff.reset(new LWH::Histogram1D(10, 1.9, 4.9));
      _h_totN_nondiff.reset(new LWH::Histogram1D(10, 1.9, 4.9));

      // Counters
      _mbSumW = 0.0; _hdSumW = 0.0; _dfSumW = 0.0; _ndSumW = 0.0;
      _mbchSumW = 0.0; _hdchSumW = 0.0; _dfchSumW = 0.0; _ndchSumW = 0.0;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const FinalState& ffs = applyProjection<FinalState>(event, "forwardFS");
      const FinalState& bfs = applyProjection<FinalState>(event, "backwardFS");
      const ChargedFinalState& fcfs = applyProjection<ChargedFinalState>(event, "forwardCFS");
      const ChargedFinalState& bcfs = applyProjection<ChargedFinalState>(event, "backwardCFS");

      // Veto this event completely if there are no forward *charged* particles
      if (fcfs.empty()) vetoEvent;

      // Charged and neutral version
      {
        // Decide empirically if this is a "hard" or "diffractive" event
        bool ishardEvt = false;
        foreach (const Particle& p, ffs.particles()) {
          if (p.momentum().pT() > 3.0) { ishardEvt = true; break; }
        }
        // Decide empirically if this is a "diffractive" event
        /// @todo Can be "diffractive" *and* "hard"?
        bool isdiffEvt = (bfs.size() == 0);

        // Update event-type weight counters
        _mbSumW += weight;
        (isdiffEvt ? _dfSumW : _ndSumW) += weight;
        if (ishardEvt) _hdSumW += weight;

        // Plot energy flow
        foreach (const Particle& p, ffs.particles()) {
          const double eta = p.momentum().eta();
          const double energy = p.momentum().E();
          _tp_totEF_minbias->fill(eta, energy, weight);
          _h_totN_minbias->fill(eta, weight);
          if (ishardEvt) {
            _tp_totEF_hard->fill(eta, energy, weight);
            _h_totN_hard->fill(eta, weight);
          }
          if (isdiffEvt) {
            _tp_totEF_diff->fill(eta, energy, weight);
            _h_totN_diff->fill(eta, weight);
          } else {
            _tp_totEF_nondiff->fill(eta, energy, weight);
            _h_totN_nondiff->fill(eta, weight);
          }
        }
      }


      // Charged-only version
      {
        bool ishardEvt = false;
        foreach (const Particle& p, fcfs.particles()) {
          if (p.momentum().pT() > 3.0) { ishardEvt = true; break; }
        }
        // Decide empirically if this is a "diffractive" event
        /// @todo Can be "diffractive" *and* "hard"?
        bool isdiffEvt = (bcfs.size() == 0);

        // Update event-type weight counters
        _mbchSumW += weight;
        (isdiffEvt ? _dfchSumW : _ndchSumW) += weight;
        if (ishardEvt) _hdchSumW += weight;

        // Plot energy flow
        foreach (const Particle& p, fcfs.particles()) {
          const double eta = p.momentum().eta();
          const double energy = p.momentum().E();
          _tp_chEF_minbias->fill(eta, energy, weight);
          _h_chN_minbias->fill(eta, weight);
          if (ishardEvt) {
            _tp_chEF_hard->fill(eta, energy, weight);
            _h_chN_hard->fill(eta, weight);
          }
          if (isdiffEvt) {
            _tp_chEF_diff->fill(eta, energy, weight);
            _h_chN_diff->fill(eta, weight);
          } else {
            _tp_chEF_nondiff->fill(eta, energy, weight);
            _h_chN_nondiff->fill(eta, weight);
          }
        }
      }

    }


    void finalize() {
      double norm = -1.0;
      double err = -1.0;
      double val = -1.0;
      IMeasurement* mpt;
      for (int i = 0; i < _dps_totEF_minbias->size(); ++i) {
        mpt = _dps_totEF_minbias->point(i)->coordinate(1);
        norm = 1.0/0.3/_mbSumW;
        val = _tp_totEF_minbias->binHeight(i)*norm*_h_totN_minbias->binHeight(i);
        err = _tp_totEF_minbias->binHeight(i)*norm*_h_totN_minbias->binError(i) + _tp_totEF_minbias->binError(i)*norm*_h_totN_minbias->binHeight(i);
        mpt->setValue(val);
        mpt->setErrorPlus(err);
        mpt->setErrorMinus(err);

        mpt = _dps_totEF_hard->point(i)->coordinate(1);
        norm = 1.0/0.3/_hdSumW;
        val = _tp_totEF_hard->binHeight(i)*norm*_h_totN_hard->binHeight(i);
        err = _tp_totEF_hard->binHeight(i)*norm*_h_totN_hard->binError(i) + _tp_totEF_hard->binError(i)*norm*_h_totN_hard->binHeight(i);
        mpt->setValue(val);
        mpt->setErrorPlus(err);
        mpt->setErrorMinus(err);

        mpt = _dps_totEF_diff->point(i)->coordinate(1);
        norm = 1.0/0.3/_dfSumW;
        val = _tp_totEF_diff->binHeight(i)*norm*_h_totN_diff->binHeight(i);
        err = _tp_totEF_diff->binHeight(i)*norm*_h_totN_diff->binError(i) + _tp_totEF_diff->binError(i)*norm*_h_totN_diff->binHeight(i);
        mpt->setValue(val);
        mpt->setErrorPlus(err);
        mpt->setErrorMinus(err);

        mpt = _dps_totEF_nondiff->point(i)->coordinate(1);
        norm = 1.0/0.3/_ndSumW;
        val = _tp_totEF_nondiff->binHeight(i)*norm*_h_totN_nondiff->binHeight(i);
        err = _tp_totEF_nondiff->binHeight(i)*norm*_h_totN_nondiff->binError(i) + _tp_totEF_nondiff->binError(i)*norm*_h_totN_nondiff->binHeight(i);
        mpt->setValue(val);
        mpt->setErrorPlus(err);
        mpt->setErrorMinus(err);

        mpt = _dps_chEF_minbias->point(i)->coordinate(1);
        norm = 1.0/0.3/_mbchSumW;
        val = _tp_chEF_minbias->binHeight(i)*norm*_h_chN_minbias->binHeight(i);
        err = _tp_chEF_minbias->binHeight(i)*norm*_h_chN_minbias->binError(i) + _tp_chEF_minbias->binError(i)*norm*_h_chN_minbias->binHeight(i);
        mpt->setValue(val);
        mpt->setErrorPlus(err);
        mpt->setErrorMinus(err);

        mpt = _dps_chEF_hard->point(i)->coordinate(1);
        norm = 1.0/0.3/_hdchSumW;
        val = _tp_chEF_hard->binHeight(i)*norm*_h_chN_hard->binHeight(i);
        err = _tp_chEF_hard->binHeight(i)*norm*_h_chN_hard->binError(i) + _tp_chEF_hard->binError(i)*norm*_h_chN_hard->binHeight(i);
        mpt->setValue(val);
        mpt->setErrorPlus(err);
        mpt->setErrorMinus(err);

        mpt = _dps_chEF_diff->point(i)->coordinate(1);
        norm = 1.0/0.3/_dfchSumW;
        val = _tp_chEF_diff->binHeight(i)*norm*_h_chN_diff->binHeight(i);
        err = _tp_chEF_diff->binHeight(i)*norm*_h_chN_diff->binError(i) + _tp_chEF_diff->binError(i)*norm*_h_chN_diff->binHeight(i);
        mpt->setValue(val);
        mpt->setErrorPlus(err);
        mpt->setErrorMinus(err);

        mpt = _dps_chEF_nondiff->point(i)->coordinate(1);
        norm = 1.0/0.3/_ndchSumW;
        val = _tp_chEF_nondiff->binHeight(i)*norm*_h_chN_nondiff->binHeight(i);
        err = _tp_chEF_nondiff->binHeight(i)*norm*_h_chN_nondiff->binError(i) + _tp_chEF_nondiff->binError(i)*norm*_h_chN_nondiff->binHeight(i);
        mpt->setValue(val);
        mpt->setErrorPlus(err);
        mpt->setErrorMinus(err);
      }
    }


  private:

    // Histograms correspond to charged and total EF for each class of events:
    // minimum bias, hard scattering, diffractive enriched and non-diffractive enriched

    AIDA::IDataPointSet* _dps_chEF_minbias; // will be filled in finalize with 1/d_eta <N(eta)><E(eta)>
    shared_ptr<LWH::Profile1D> _tp_chEF_minbias;  // contains <E(eta)>
    shared_ptr<LWH::Histogram1D> _h_chN_minbias; // contains N(eta)
    double _mbSumW; // sum of weights (~ #events) in each event class

    AIDA::IDataPointSet* _dps_chEF_hard;
    shared_ptr<LWH::Profile1D> _tp_chEF_hard;
    shared_ptr<LWH::Histogram1D> _h_chN_hard;
    double _hdSumW;

    AIDA::IDataPointSet* _dps_chEF_diff;
    shared_ptr<LWH::Profile1D> _tp_chEF_diff;
    shared_ptr<LWH::Histogram1D> _h_chN_diff;
    double _dfSumW;

    AIDA::IDataPointSet* _dps_chEF_nondiff;
    shared_ptr<LWH::Profile1D> _tp_chEF_nondiff;
    shared_ptr<LWH::Histogram1D> _h_chN_nondiff;
    double _ndSumW;

    AIDA::IDataPointSet* _dps_totEF_minbias;
    shared_ptr<LWH::Profile1D> _tp_totEF_minbias;
    shared_ptr<LWH::Histogram1D> _h_totN_minbias;
    double _mbchSumW;

    AIDA::IDataPointSet* _dps_totEF_hard;
    shared_ptr<LWH::Profile1D> _tp_totEF_hard;
    shared_ptr<LWH::Histogram1D> _h_totN_hard;
    double _hdchSumW;

    AIDA::IDataPointSet* _dps_totEF_diff;
    shared_ptr<LWH::Profile1D> _tp_totEF_diff;
    shared_ptr<LWH::Histogram1D> _h_totN_diff;
    double _dfchSumW;

    AIDA::IDataPointSet* _dps_totEF_nondiff;
    shared_ptr<LWH::Profile1D> _tp_totEF_nondiff;
    shared_ptr<LWH::Histogram1D> _h_totN_nondiff;
    double _ndchSumW;
  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(LHCB_2013_I1208105);

}
