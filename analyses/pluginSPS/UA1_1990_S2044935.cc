// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// UA1 minbias track multiplicities, \f$ p_\perp \f$ and \f$ E_\perp \f$
  class UA1_1990_S2044935 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(UA1_1990_S2044935);


    /// @name Analysis methods
    //@{

    /// Book projections and histograms
    void init() {
      declare(ChargedFinalState((Cuts::etaIn(-5.5, 5.5))), "TriggerFS");
      declare(ChargedFinalState((Cuts::etaIn(-2.5, 2.5))), "TrackFS");
      const FinalState trkcalofs((Cuts::etaIn(-2.5, 2.5)));
      declare(MissingMomentum(trkcalofs), "MET25");
      const FinalState calofs((Cuts::etaIn(-6.0, 6.0)));
      declare(MissingMomentum(calofs), "MET60");

      if (isCompatibleWithSqrtS(63*GeV)) {
        book(_hist_Pt ,8,1,1);
      } else if (isCompatibleWithSqrtS(200*GeV)) {
        book(_hist_Nch ,1,1,1);
        book(_hist_Esigd3p ,2,1,1);
        book(_hist_Pt ,6,1,1);
        book(_hist_Et ,9,1,1);
        book(_hist_Etavg ,12,1,1);
      } else if (isCompatibleWithSqrtS(500*GeV)) {
        book(_hist_Nch ,1,1,2);
        book(_hist_Esigd3p ,2,1,2);
        book(_hist_Et ,10,1,1);
        book(_hist_Etavg ,12,1,2);
      } else if (isCompatibleWithSqrtS(900*GeV)) {
        book(_hist_Nch ,1,1,3);
        book(_hist_Esigd3p ,2,1,3);
        book(_hist_Pt ,7,1,1);
        book(_hist_Et ,11,1,1);
        book(_hist_Etavg ,12,1,3);
        book(_hist_Esigd3p08 ,3,1,1);
        book(_hist_Esigd3p40 ,4,1,1);
        book(_hist_Esigd3p80 ,5,1,1);
      }
      book(_sumwTrig, "TMP/sumwTrig");
      // book(_sumwTrig08, "TMP/sumwTrig08");
      // book(_sumwTrig40, "TMP/sumwTrig40");
      // book(_sumwTrig80, "TMP/sumwTrig80");

    }


    void analyze(const Event& event) {
      // Trigger
      const FinalState& trigfs = apply<FinalState>(event, "TriggerFS");
      unsigned int n_minus(0), n_plus(0);
      for (const Particle& p : trigfs.particles()) {
        const double eta = p.eta();
        if (inRange(eta, -5.5, -1.5)) n_minus++;
        else if (inRange(eta, 1.5, 5.5)) n_plus++;
      }
      MSG_DEBUG("Trigger -: " << n_minus << ", Trigger +: " << n_plus);
      if (n_plus == 0 || n_minus == 0) vetoEvent;
      _sumwTrig->fill();

      // Use good central detector tracks
      const FinalState& cfs = apply<FinalState>(event, "TrackFS");
      const double Et25 = apply<MissingMomentum>(event, "MET25").scalarEt();
      const double Et60 = apply<MissingMomentum>(event, "MET60").scalarEt();
      const unsigned int nch = cfs.size();

      // Event level histos
      if (!isCompatibleWithSqrtS(63*GeV)) {
        _hist_Nch->fill(nch);
        _hist_Et->fill(Et60/GeV);
        _hist_Etavg->fill(nch, Et25/GeV);
      }

      // Particle/track level histos
      const double deta = 2 * 5.0;
      const double dphi = TWOPI;
      const double dnch_deta = nch/deta;
      for (const Particle& p : cfs.particles()) {
        const double pt = p.pT();
        const double scaled_weight = 1.0/(deta*dphi*pt/GeV);
        if (!isCompatibleWithSqrtS(500*GeV)) {
          _hist_Pt->fill(nch, pt/GeV);
        }
        if (!isCompatibleWithSqrtS(63*GeV)) {
          _hist_Esigd3p->fill(pt/GeV, scaled_weight);
        }
        // Also fill for specific dn/deta ranges at 900 GeV
        if (isCompatibleWithSqrtS(900*GeV, 1E-3)) {
          if (inRange(dnch_deta, 0.8, 4.0)) {
            //_sumwTrig08 ->fill();
            _hist_Esigd3p08->fill(pt/GeV, scaled_weight);
          } else if (inRange(dnch_deta, 4.0, 8.0)) {
            //_sumwTrig40 ->fill();
            _hist_Esigd3p40->fill(pt/GeV, scaled_weight);
          } else {
            //MSG_WARNING(dnch_deta);
            if (dnch_deta > 8.0) {
              //_sumwTrig80 ->fill();
              _hist_Esigd3p80->fill(pt/GeV, scaled_weight);
            }
          }
        }
      }

    }


    void finalize() {
      if (_sumwTrig->val() <= 0) {
        MSG_WARNING("No events passed the trigger!");
        return;
      }
      const double xsec = crossSectionPerEvent();
      if (!isCompatibleWithSqrtS(63*GeV)) {
        scale(_hist_Nch, 2*xsec/millibarn); ///< Factor of 2 for Nch bin widths?
        scale(_hist_Esigd3p, xsec/millibarn);
        scale(_hist_Et, xsec/millibarn);
      }
      if (isCompatibleWithSqrtS(900*GeV)) {
        // NB. Ref data is normalised to a fixed value not reproducible from MC.
        const double scale08 =  (_hist_Esigd3p08->bin(0).area() > 0) ?
          0.933e5/_hist_Esigd3p08->bin(0).height() : 0;
        scale(_hist_Esigd3p08, scale08);
        const double scale40 = (_hist_Esigd3p40->bin(0).area() > 0) ?
          1.369e5/_hist_Esigd3p40->bin(0).height() : 0;
        scale(_hist_Esigd3p40, scale40);
        const double scale80 = (_hist_Esigd3p80->bin(0).area() > 0) ?
          1.657e5/_hist_Esigd3p80->bin(0).height() : 0;
        scale(_hist_Esigd3p80, scale80);
      }
    }

    //@}


  private:

    /// @name Weight counters
    //@{
    CounterPtr _sumwTrig; //, _sumwTrig08, _sumwTrig40, _sumwTrig80;
    //@}

    /// @name Histogram collections
    //@{
    Histo1DPtr _hist_Nch;
    Histo1DPtr _hist_Esigd3p;
    Histo1DPtr _hist_Esigd3p08;
    Histo1DPtr _hist_Esigd3p40;
    Histo1DPtr _hist_Esigd3p80;
    Profile1DPtr _hist_Pt;
    Profile1DPtr _hist_Etavg;
    Histo1DPtr _hist_Et;
    //@}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(UA1_1990_S2044935, UA1_1990_I280412);

}
