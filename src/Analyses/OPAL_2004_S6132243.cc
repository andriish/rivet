// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"

namespace Rivet {


  class OPAL_2004_S6132243 : public Analysis {
  public:

    /// Constructor
    OPAL_2004_S6132243() : Analysis("OPAL_2004_S6132243") {
      //
    }


    /// @name Analysis methods
    //@{

    void init() {
      // Projections
      addProjection(Beam(), "Beams");
      const ChargedFinalState cfs;
      addProjection(cfs, "FS");
      addProjection(FastJets(cfs, FastJets::DURHAM, 0.7), "DurhamJets");
      addProjection(Sphericity(cfs), "Sphericity");
      addProjection(ParisiTensor(cfs), "Parisi");
      const Thrust thrust(cfs);
      addProjection(thrust, "Thrust");
      addProjection(Hemispheres(thrust), "Hemispheres");

      // Counters
      _sumPassedWeights = 0;

      // Book histograms
      // Energies: 91, 133, 177 (161-183), 197 (189-209) => index 0..4
      for (int isqrts = 0; isqrts < 4; ++isqrts) {
        _hist1MinusT[isqrts]    = bookHistogram1D(1, 1, isqrts+1);
        _histHemiMassH[isqrts]  = bookHistogram1D(2, 1, isqrts+1);
        _histCParam[isqrts]     = bookHistogram1D(3, 1, isqrts+1);
        _histHemiBroadT[isqrts] = bookHistogram1D(4, 1, isqrts+1);
        _histHemiBroadW[isqrts] = bookHistogram1D(5, 1, isqrts+1);
        _histY23Durham[isqrts]  = bookHistogram1D(6, 1, isqrts+1);
        _histTMajor[isqrts]     = bookHistogram1D(7, 1, isqrts+1);
        _histTMinor[isqrts]     = bookHistogram1D(8, 1, isqrts+1);
        _histAplanarity[isqrts] = bookHistogram1D(9, 1, isqrts+1);
        _histSphericity[isqrts] = bookHistogram1D(10, 1, isqrts+1);
        _histOblateness[isqrts] = bookHistogram1D(11, 1, isqrts+1);
        _histHemiMassL[isqrts]  = bookHistogram1D(12, 1, isqrts+1);
        _histHemiBroadN[isqrts] = bookHistogram1D(13, 1, isqrts+1);
        _histDParam[isqrts]     = bookHistogram1D(14, 1, isqrts+1);
        //
        _hist1MinusTMom[isqrts]    = bookHistogram1D(15, 1, isqrts+1);
        _histHemiMassHMom[isqrts]  = bookHistogram1D(16, 1, isqrts+1);
        _histCParamMom[isqrts]     = bookHistogram1D(17, 1, isqrts+1);
        _histHemiBroadTMom[isqrts] = bookHistogram1D(18, 1, isqrts+1);
        _histHemiBroadWMom[isqrts] = bookHistogram1D(19, 1, isqrts+1);
        _histY23DurhamMom[isqrts]  = bookHistogram1D(20, 1, isqrts+1);
        _histTMajorMom[isqrts]     = bookHistogram1D(21, 1, isqrts+1);
        _histTMinorMom[isqrts]     = bookHistogram1D(22, 1, isqrts+1);
        _histSphericityMom[isqrts] = bookHistogram1D(23, 1, isqrts+1);
        _histOblatenessMom[isqrts] = bookHistogram1D(24, 1, isqrts+1);
        _histHemiMassLMom[isqrts]  = bookHistogram1D(25, 1, isqrts+1);
        _histHemiBroadNMom[isqrts] = bookHistogram1D(26, 1, isqrts+1);
      }
    }


    void analyze(const Event& event) {
      const FinalState& cfs = applyProjection<FinalState>(event, "FS");
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (cfs.size() < 2) vetoEvent;

      // Increment passed-cuts weight sum
      const double weight = event.weight();
      _sumPassedWeights += weight;

      // Get beams and average beam momentum
      const double sqrts = applyProjection<Beam>(event, "Beams").sqrtS();

      // Translate sqrt(s) into a histo index for the majority of histos
      size_t ih = 5;
      if (inRange(sqrts/GeV, 89.9, 91.5)) {
        ih = 0;
      } else if (fuzzyEquals(sqrts/GeV, 133)) {
        ih = 1;
      } else if (fuzzyEquals(sqrts/GeV, 177)) { // (161-183)
        ih = 2;
      } else if (fuzzyEquals(sqrts/GeV, 197)) { // (189-209)
        ih = 3;
      } else {
        throw Error("Invalid energy for OPAL_2004 analysis! Require 91, 133, 177, or 197 GeV");
      }

      // Thrusts
      const Thrust& thrust = applyProjection<Thrust>(event, "Thrust");
      _hist1MinusT[ih]->fill(1-thrust.thrust(), weight);
      _histTMajor[ih]->fill(thrust.thrustMajor(), weight);
      _histTMinor[ih]->fill(thrust.thrustMinor(), weight);
      _histOblateness[ih]->fill(thrust.oblateness(), weight);
      for (int n = 1; n <= 5; ++n) {
        _hist1MinusTMom[ih]->fill(n, pow(1-thrust.thrust(), n)*weight);
        _histTMajorMom[ih]->fill(n, pow(thrust.thrustMajor(), n)*weight);
        _histTMinorMom[ih]->fill(n, pow(thrust.thrustMinor(), n)*weight);
        _histOblatenessMom[ih]->fill(n, pow(thrust.oblateness(), n)*weight);
      }

      // Jets
      const FastJets& durjet = applyProjection<FastJets>(event, "DurhamJets");
      if (durjet.clusterSeq()) {
        /// @todo Need separate normalisation due to clusterseq / 3 jet requirement?
        const double y23 = durjet.clusterSeq()->exclusive_ymerge(3);
        _histY23Durham[ih]->fill(y23, weight);
        for (int n = 1; n <= 5; ++n) {
          _histY23DurhamMom[ih]->fill(n, pow(y23, n)*weight);
        }
      }

      // Sphericities
      const Sphericity& sphericity = applyProjection<Sphericity>(event, "Sphericity");
      const double sph = sphericity.sphericity();
      const double apl = sphericity.aplanarity();
      _histSphericity[ih]->fill(sph, weight);
      _histAplanarity[ih]->fill(apl, weight);
      for (int n = 1; n <= 5; ++n) {
        _histSphericityMom[ih]->fill(n, pow(sph, n)*weight);
      }

      // C & D params
      const ParisiTensor& parisi = applyProjection<ParisiTensor>(event, "Parisi");
      const double cparam = parisi.C();
      const double dparam = parisi.D();
      _histCParam[ih]->fill(cparam, weight);
      _histDParam[ih]->fill(dparam, weight);
      for (int n = 1; n <= 5; ++n) {
        _histCParamMom[ih]->fill(n, pow(cparam, n)*weight);
      }
   
      // Hemispheres
      const Hemispheres& hemi = applyProjection<Hemispheres>(event, "Hemispheres");
      const double hemi_mh = hemi.scaledM2high();
      const double hemi_ml = hemi.scaledM2low();
      const double hemi_bmax = hemi.Bmax();
      const double hemi_bmin = hemi.Bmin();
      const double hemi_bsum = hemi.Bsum();
      _histHemiMassH[ih]->fill(hemi_mh, weight);
      _histHemiMassL[ih]->fill(hemi_ml, weight);
      _histHemiBroadW[ih]->fill(hemi_bmax, weight);
      _histHemiBroadN[ih]->fill(hemi_bmin, weight);
      _histHemiBroadT[ih]->fill(hemi_bsum, weight);
      for (int n = 1; n <= 5; ++n) {
        _histHemiMassHMom[ih]->fill(n, pow(hemi_mh, n)*weight);
        _histHemiMassLMom[ih]->fill(n, pow(hemi_ml, n)*weight);
        _histHemiBroadWMom[ih]->fill(n, pow(hemi_bmax, n)*weight);
        _histHemiBroadNMom[ih]->fill(n, pow(hemi_bmin, n)*weight);
        _histHemiBroadTMom[ih]->fill(n, pow(hemi_bsum, n)*weight);
      }
    }


    void finalize() {
      /// @todo Normalisations / scalings, etc.
      for (int isqrts = 0; isqrts < 4; ++isqrts) {
        normalize(_hist1MinusT[isqrts]);
        normalize(_histTMajor[isqrts]);
        normalize(_histTMinor[isqrts]);
        normalize(_histOblateness[isqrts]);
        normalize(_histSphericity[isqrts]);
        normalize(_histAplanarity[isqrts]);
        normalize(_histHemiMassH[isqrts]);
        normalize(_histHemiMassL[isqrts]);
        normalize(_histHemiBroadW[isqrts]);
        normalize(_histHemiBroadN[isqrts]);
        normalize(_histHemiBroadT[isqrts]);
        normalize(_histCParam[isqrts]);
        normalize(_histDParam[isqrts]);
        normalize(_histY23Durham[isqrts]);
        //
        scale(_hist1MinusTMom[isqrts], 1.0/_sumPassedWeights);
        scale(_histTMajorMom[isqrts], 1.0/_sumPassedWeights);
        scale(_histTMinorMom[isqrts], 1.0/_sumPassedWeights);
        scale(_histOblatenessMom[isqrts], 1.0/_sumPassedWeights);
        scale(_histSphericityMom[isqrts], 1.0/_sumPassedWeights);
        scale(_histHemiMassHMom[isqrts], 1.0/_sumPassedWeights);
        scale(_histHemiMassLMom[isqrts], 1.0/_sumPassedWeights);
        scale(_histHemiBroadWMom[isqrts], 1.0/_sumPassedWeights);
        scale(_histHemiBroadNMom[isqrts], 1.0/_sumPassedWeights);
        scale(_histHemiBroadTMom[isqrts], 1.0/_sumPassedWeights);
        scale(_histCParamMom[isqrts], 1.0/_sumPassedWeights);
        scale(_histY23DurhamMom[isqrts], 1.0/_sumPassedWeights);
      }
    }

    //@}


  private:

    // Counter of event weights passing the cuts
    double _sumPassedWeights;

    // Event shape histos at 4 energies
    AIDA::IHistogram1D* _hist1MinusT[4];
    AIDA::IHistogram1D* _histHemiMassH[4];
    AIDA::IHistogram1D* _histCParam[4];
    AIDA::IHistogram1D* _histHemiBroadT[4];
    AIDA::IHistogram1D* _histHemiBroadW[4];
    AIDA::IHistogram1D* _histY23Durham[4];
    AIDA::IHistogram1D* _histTMajor[4];
    AIDA::IHistogram1D* _histTMinor[4];
    AIDA::IHistogram1D* _histAplanarity[4];
    AIDA::IHistogram1D* _histSphericity[4];
    AIDA::IHistogram1D* _histOblateness[4];
    AIDA::IHistogram1D* _histHemiMassL[4];
    AIDA::IHistogram1D* _histHemiBroadN[4];
    AIDA::IHistogram1D* _histDParam[4];

    // Event shape moment histos at 4 energies
    AIDA::IHistogram1D* _hist1MinusTMom[4];
    AIDA::IHistogram1D* _histHemiMassHMom[4];
    AIDA::IHistogram1D* _histCParamMom[4];
    AIDA::IHistogram1D* _histHemiBroadTMom[4];
    AIDA::IHistogram1D* _histHemiBroadWMom[4];
    AIDA::IHistogram1D* _histY23DurhamMom[4];
    AIDA::IHistogram1D* _histTMajorMom[4];
    AIDA::IHistogram1D* _histTMinorMom[4];
    AIDA::IHistogram1D* _histSphericityMom[4];
    AIDA::IHistogram1D* _histOblatenessMom[4];
    AIDA::IHistogram1D* _histHemiMassLMom[4];
    AIDA::IHistogram1D* _histHemiBroadNMom[4];

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<OPAL_2004_S6132243> plugin_OPAL_2004_S6132243;

}
