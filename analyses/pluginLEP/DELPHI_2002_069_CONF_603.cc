// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief DELPHI b-fragmentation measurement
  /// @author Hendrik Hoeth
  class DELPHI_2002_069_CONF_603 : public Analysis {
  public:

    /// Constructor
    DELPHI_2002_069_CONF_603()
      : Analysis("DELPHI_2002_069_CONF_603")
    {    }


    /// @name Helper functions
    /// @note The PID:: namespace functions would be preferable, but don't have exactly the same behaviour. Preserving the original form.
    //@{
    bool isParton(int id) { return abs(id) <= 100 && abs(id) != 22 && (abs(id) < 11 || abs(id) > 18); }
    // bool isBHadron(int id) { return ((abs(id)/100)%10 == 5) || (abs(id) >= 5000 && abs(id) <= 5999); }
    //@}


    /// @name Analysis methods
    //@{

    /// Book projections and histograms
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");

      _histXbprim     = bookHisto1D(1, 1, 1);
      _histXbweak     = bookHisto1D(2, 1, 1);
      _histMeanXbprim = bookProfile1D(4, 1, 1);
      _histMeanXbweak = bookProfile1D(5, 1, 1);
    }


    void analyze(const Event& e) {
      const FinalState& fs = apply<FinalState>(e, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed ncharged cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed ncharged cut");

      // Get event weight for histo filling
      const double weight = e.weight();

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);


      for (ConstGenParticlePtr p : particles(e.genEvent())) {
        ConstGenVertexPtr pv = p->production_vertex();
        ConstGenVertexPtr dv = p->end_vertex();
        if (PID::isBottomHadron(p->pdg_id())) {
          const double xp = p->momentum().e()/meanBeamMom;

          // If the B-hadron has a parton as parent, call it primary B-hadron:
          if (pv) {
            bool is_primary = false;
            for (ConstGenParticlePtr pp: pv->particles_in()){
              if (isParton(pp->pdg_id())) is_primary = true;
            }
            if (is_primary) {
              _histXbprim->fill(xp, weight);
              _histMeanXbprim->fill(_histMeanXbprim->bin(0).xMid(), xp, weight);
            }
          }

          // If the B-hadron has no B-hadron as a child, it decayed weakly:
          if (dv) {
            bool is_weak = true;
            for (ConstGenParticlePtr pp: dv->particles_out()){
              if (PID::isBottomHadron(pp->pdg_id())) {
                is_weak = false;
              }
            }
            if (is_weak) {
              _histXbweak->fill(xp, weight);
              _histMeanXbweak->fill(_histMeanXbweak->bin(0).xMid(), xp, weight);
            }
          }

        }
      }
    }


    // Finalize
    void finalize() {
      normalize(_histXbprim);
      normalize(_histXbweak);
    }


  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the
    /// inclusive single particle distributions' normalisations.

    Histo1DPtr _histXbprim;
    Histo1DPtr _histXbweak;

    Profile1DPtr _histMeanXbprim;
    Profile1DPtr _histMeanXbweak;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_2002_069_CONF_603);

}
