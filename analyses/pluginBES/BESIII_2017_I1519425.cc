// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Differential decay rates for $D^+\to \{\bar{K}^0,\pi^0\} e^+\nu_e$ from BES
  class BESIII_2017_I1519425 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1519425);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_q2_K , 1, 1, 1);
      book(_h_q2_pi, 2, 1, 1);

    }

    // Calculate the Q2 using mother and daugher meson
    double q2(const Particle& B, int mesonID) {
      FourMomentum q = B.mom() - filter_select(B.children(), Cuts::pid==mesonID)[0];
      return q*q;
    }

    // Check for explicit decay into pdgids
    bool isSemileptonicDecay(const Particle& mother, vector<int> ids) {
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      // Check for the explicit decay
      return all(ids, [&](int i){return count(children, hasPID(i))==1;});
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over D+/- mesons
      for (const Particle& p :  apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::DPLUS)) {
        if (isSemileptonicDecay(p, {PID::PI0, PID::POSITRON, PID::NU_E}) ||
            isSemileptonicDecay(p, {PID::PI0, PID::ELECTRON, PID::NU_EBAR})) {
          _h_q2_pi->fill(q2(p, PID::PI0));
        }
        else if(isSemileptonicDecay(p, {-311, PID::POSITRON, PID::NU_E})) {
          _h_q2_K ->fill(q2(p, -311));
        }
        else if(isSemileptonicDecay(p, { 311, PID::ELECTRON, PID::NU_EBAR})) {
          _h_q2_K ->fill(q2(p, 311));
        }
        else if(isSemileptonicDecay(p, {PID::K0S, PID::POSITRON, PID::NU_E}) ||
                isSemileptonicDecay(p, {PID::K0S, PID::ELECTRON, PID::NU_EBAR})) {
          _h_q2_K ->fill(q2(p, PID::K0S));
        }
        else if(isSemileptonicDecay(p, {PID::K0L, PID::POSITRON, PID::NU_E}) ||
                isSemileptonicDecay(p, {PID::K0L, PID::ELECTRON, PID::NU_EBAR})) {
          _h_q2_K ->fill(q2(p, PID::K0L));
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize to unity
      normalize(_h_q2_K );
      normalize(_h_q2_pi);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_q2_K,_h_q2_pi;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BESIII_2017_I1519425);


}
