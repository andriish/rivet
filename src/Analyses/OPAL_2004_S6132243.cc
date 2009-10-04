// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"

namespace Rivet {


  class OPAL_2004_S6132243 : public Analysis { 
  public:

    /// Constructor
    OPAL_2004_S6132243() : Analysis("OPAL_2004_S6132243") { 

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

      // Book histograms
      // Energies: 91, 133, 177 (161-183), 197 (189-209) => index 0..4
      for (int isqrts = 0; isqrts < 4; ++isqrts) {
        _hist1MinusT[isqrts] = bookHistogram(1, 1, isqrts+1);
        _histHemiMassH[isqrts] = bookHistogram(2, 1, isqrts+1);
        _histCParam[isqrts] = bookHistogram(3, 1, isqrts+1);
        _histHemiBroadT[isqrts] = bookHistogram(4, 1, isqrts+1);
        _histHemiBroadW[isqrts] = bookHistogram(5, 1, isqrts+1);
        _histY23Durham[isqrts] = bookHistogram(6, 1, isqrts+1);
        _histTMajor[isqrts] = bookHistogram(7, 1, isqrts+1);
        _histTMinor[isqrts] = bookHistogram(8, 1, isqrts+1);
        _histAplanarity[isqrts] = bookHistogram(9, 1, isqrts+1);
        _histSphericity[isqrts] = bookHistogram(10, 1, isqrts+1);
        _histOblateness[isqrts] = bookHistogram(11, 1, isqrts+1);
        _histHemiMassL[isqrts] = bookHistogram(12, 1, isqrts+1);
        _histHemiBroadN[isqrts] = bookHistogram(13, 1, isqrts+1);
        _histDParam[isqrts] = bookHistogram(14, 1, isqrts+1);
        //
        _hist1MinusTMom[isqrts] = bookHistogram(15, 1, isqrts+1);
        _histHemiMassHMom[isqrts] = bookHistogram(16, 1, isqrts+1);
        _histCParamMom[isqrts] = bookHistogram(17, 1, isqrts+1);
        _histHemiBroadTMom[isqrts] = bookHistogram(18, 1, isqrts+1);
        _histHemiBroadWMom[isqrts] = bookHistogram(19, 1, isqrts+1);
        _histY23DurhamMom[isqrts] = bookHistogram(20, 1, isqrts+1);
        _histTMajorMom[isqrts] = bookHistogram(21, 1, isqrts+1);
        _histTMinorMom[isqrts] = bookHistogram(22, 1, isqrts+1);
        _histSphericityMom[isqrts] = bookHistogram(23, 1, isqrts+1);
        _histOblatenessMom[isqrts] = bookHistogram(24, 1, isqrts+1);
        _histHemiMassLMom[isqrts] = bookHistogram(25, 1, isqrts+1);
        _histHemiBroadNMom[isqrts] = bookHistogram(26, 1, isqrts+1);
      }
      //
      _hist1MinusTVar[isqrts] = bookHistogram(27, 1, 1);
      _histHemiMassHVar[isqrts] = bookHistogram(27, 1, 2);
      _histCParamVar[isqrts] = bookHistogram(27, 1, 3);
      _histHemiBroadTVar[isqrts] = bookHistogram(27, 1, 4);
      _histHemiBroadWVar[isqrts] = bookHistogram(27, 1, 5);
      _histY23DurhamVar[isqrts] = bookHistogram(27, 1, 6);
      _histTMajorVar[isqrts] = bookHistogram(27, 1, 7);
      _histTMinorVar[isqrts] = bookHistogram(27, 1, 8);
      _histSphericityVar[isqrts] = bookHistogram(27, 1, 9);
      _histOblatenessVar[isqrts] = bookHistogram(27, 1, 10);
      _histHemiMassLVar[isqrts] = bookHistogram(27, 1, 11);
      _histHemiBroadNVar[isqrts] = bookHistogram(27, 1, 12);
    }


    void analyze(const Event & event) { 
      const FinalState& cfs = applyProjection<FinalState>(e, "CFS");
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (cfs.size() < 2) vetoEvent;
      const double weight = e.weight();

      // Get beams and average beam momentum
      const ParticlePair& beams = applyProjection<Beam>(e, "Beams").beams();
      const double sqrts = beams.sqrtS();

      // Translate sqrt(s) into a histo index for the majority of histos
      size_t ih = 5;
      if (fuzzyEquals(sqrts/GeV, 91)) {
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
      const Thrust& thrust = applyProjection<Thrust>(e, "Thrust");
      _hist1MinusT[ih]->fill(1-thrust.thrust(), weight);
      _histTMajor[ih]->fill(thrust.thrustMajor(), weight); 
      _histTMinor[ih]->fill(thrust.thrustMinor(), weight); 
      _histOblateness[ih]->fill(thrust.oblateness(), weight);
      for (int n = 1; n <= 5; ++n) {
        _hist1MinusTMom[ih]->fill(n, pow(n, 1-thrust.thrust()), weight);
        _histTMajorMom[ih]->fill(n, pow(n, 1-thrust.thrustMajor()), weight);
        _histTMinorMom[ih]->fill(n, pow(n, 1-thrust.thrustMinor()), weight);
        _histOblatenessMom[ih]->fill(n, pow(n, 1-thrust.oblateness()), weight);
      }
      /// @todo Variance counters
      // var(1-T) d27-x01-y01
      // var(T_maj) d27-x01-y07
      // var(T_min) d27-x01-y08
      // var(O) d27-x01-y10

      // Jets
      const FastJets& durjet = applyProjection<FastJets>(e, "DurhamJets");
      if (durjet.clusterSeq()) {
        _histY23Durham[ih]->fill(durjet.clusterSeq()->exclusive_ymerge(3), weight); 
      }
      // Moments 1..5 of Y23 (Dur)
      // d20-x01-y01
      // d20-x01-y02
      // d20-x01-y03
      // d20-x01-y04
      // var(Y23) d27-x01-y06

      // Sphericities
      const Sphericity& sphericity = applyProjection<Sphericity>(e, "Sphericity");
      _histSphericity[ih]->fill(sphericity.sphericity(), weight);
      _histAplanarity[ih]->fill(sphericity.aplanarity(), weight);
      // Moments 1..5 of S
      // d23-x01-y01
      // d23-x01-y02
      // d23-x01-y03
      // d23-x01-y04
      // var(S) d27-x01-y09
      
      // C & D params
      const ParisiTensor& parisi = applyProjection<ParisiTensor>(e, "Parisi");
      _histCParam[ih]->fill(parisi.C(), weight);
      _histDParam[ih]->fill(parisi.D(), weight);
      // Moments 1..5 of C
      // d17-x01-y01
      // d17-x01-y02
      // d17-x01-y03
      // d17-x01-y04
      // var(C) d27-x01-y03
      
      // Hemispheres
      const Hemispheres& hemi = applyProjection<Hemispheres>(e, "Hemispheres");
      _histHemiMassH[ih]->fill(hemi.getScaledM2high(), weight);
      _histHemiMassL[ih]->fill(hemi.getScaledM2low(), weight);
      _histHemiBroadW[ih]->fill(hemi.getBmax(), weight);
      _histHemiBroadN[ih]->fill(hemi.getBmin(), weight);
      _histHemiBroadT[ih]->fill(hemi.getBsum(), weight);
      // Moments 1..5 of M_H
      // d16-x01-y01
      // d16-x01-y02
      // d16-x01-y03
      // d16-x01-y04
      // Moments 1..5 of B_T
      // d18-x01-y01
      // d18-x01-y02
      // d18-x01-y03
      // d18-x01-y04
      // Moments 1..5 of B_W
      // d19-x01-y01
      // d19-x01-y02
      // d19-x01-y03
      // d19-x01-y04
      // Moments 1..5 of M_L
      // d25-x01-y01
      // d25-x01-y02
      // d25-x01-y03
      // d25-x01-y04
      // Moments 1..5 of B_N
      // d26-x01-y01
      // d26-x01-y02
      // d26-x01-y03
      // d26-x01-y04
      // var(M_H) d27-x01-y02
      // var(B_T) d27-x01-y04
      // var(B_W) d27-x01-y05
      // var(M_L) d27-x01-y11
      // var(B_N) d27-x01-y12 
    }


    void finalize() { 
      /// @todo Normalisations / scalings, etc.    
    }

    //@}


  private:

    // Event shape histos at 4 energies
    AIDA::IHistogram1D _hist1MinusT[4];
    AIDA::IHistogram1D _histHemiMassH[4];
    AIDA::IHistogram1D _histCParam[4];
    AIDA::IHistogram1D _histHemiBroadT[4];
    AIDA::IHistogram1D _histHemiBroadW[4];
    AIDA::IHistogram1D _histY23Durham[4];
    AIDA::IHistogram1D _histTMajor[4];
    AIDA::IHistogram1D _histTMinor[4];
    AIDA::IHistogram1D _histAplanarity[4];
    AIDA::IHistogram1D _histSphericity[4];
    AIDA::IHistogram1D _histOblateness[4];
    AIDA::IHistogram1D _histHemiMassL[4];
    AIDA::IHistogram1D _histHemiBroadN[4];
    AIDA::IHistogram1D _histDParam[4];

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<OPAL_2004_S6132243> plugin_OPAL_2004_S6132243;

}
