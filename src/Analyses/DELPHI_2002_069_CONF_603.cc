// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/DELPHI_2002_069_CONF_603.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"

/// @todo Use inline function rather than preprocessor, and use ParticleName enum values.
/// @todo Put these in a PID utility header? (isParton(id) would be a good candidate for inclusion in HepMC)
#define IS_PARTON_PDGID(id) ( abs(id) <= 100 && abs(id) != 22 && (abs(id) < 11 || abs(id) > 18) )
#define IS_BHADRON_PDGID(id) ( ((abs(id)/100)%10 == 5) || (abs(id) >= 5000 && abs(id) <= 5999) )

namespace Rivet {


  // Constructor
  DELPHI_2002_069_CONF_603::DELPHI_2002_069_CONF_603() 
  {
    setBeams(ELECTRON, POSITRON); 
    addProjection(Beam(), "Beams");
    addProjection(ChargedFinalState(), "FS");
    addProjection(InitialQuarks(), "IQF");
  }


  void DELPHI_2002_069_CONF_603::analyze(const Event& e) {
    // First, veto on leptonic events by requiring at least 4 charged FS particles
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const size_t numParticles = fs.particles().size();

    // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
    if (numParticles < 2) {
      getLog() << Log::DEBUG << "Failed ncharged cut" << endl;
      vetoEvent(e);
    }
    getLog() << Log::DEBUG << "Passed ncharged cut" << endl;

    // Get event weight for histo filling
    const double weight = e.weight();

    // Get beams and average beam momentum
    const ParticlePair& beams = applyProjection<Beam>(e, "Beams").beams();
    const double meanBeamMom = ( beams.first.momentum().vector3().mod() + 
                                 beams.second.momentum().vector3().mod() ) / 2.0;
    getLog() << Log::DEBUG << "Avg beam momentum = " << meanBeamMom << endl;


    for (GenEvent::particle_const_iterator p = e.genEvent().particles_begin();
         p != e.genEvent().particles_end(); ++p) {
      const GenVertex* pv = (*p)->production_vertex();
      const GenVertex* dv = (*p)->end_vertex();
      if (IS_BHADRON_PDGID((*p)->pdg_id())) {
        const double xp = (*p)->momentum().e()/meanBeamMom;

        // If the B-hadron has a parton as parent, call it primary B-hadron:
        if (pv!=NULL) {
          bool is_primary = false;
          for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin() ;
              pp != pv->particles_in_const_end() ; ++pp) {
            if (IS_PARTON_PDGID((*pp)->pdg_id()))
              is_primary = true;
          }
          if (is_primary) {
            _histXbprim->fill(xp, weight);
            _histMeanXbprim->fill(_histMeanXbprim->binMean(0), xp, weight);
          }
        }

        // If the B-hadron has no B-hadron as a child, it decayed weakly:
        if (dv!=NULL) {
          bool is_weak = true;
          for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin() ;
              pp != dv->particles_out_const_end() ; ++pp) {
            if (IS_BHADRON_PDGID((*pp)->pdg_id())) {
              is_weak = false;
            }
          }
          if (is_weak) {
            _histXbweak->fill(xp, weight);
            _histMeanXbweak->fill(_histMeanXbweak->binMean(0), xp, weight);
          }
        }

      }
    }
  }



  void DELPHI_2002_069_CONF_603::init() {
    _histXbprim     = bookHistogram1D(1, 1, 1, "$b$ quark fragmentation function $f(x_B^\\text{prim})$", "$x_B$", "$1/N \\, \\text{d}{N}/\\text{d}x_B$");
    _histXbweak     = bookHistogram1D(2, 1, 1, "$b$ quark fragmentation function $f(x_B^\\text{weak})$", "$x_B$", "$1/N \\, \\text{d}{N}/\\text{d}x_B$");
    _histMeanXbprim = bookProfile1D(4, 1, 1, "Mean of $b$ quark fragmentation function $f(x_B^\\text{prim})$", "", "$\\langle x_B \\rangle$");
    _histMeanXbweak = bookProfile1D(5, 1, 1, "Mean of $b$ quark fragmentation function $f(x_B^\\text{weak})$", "", "$\\langle x_B \\rangle$");
  }

  // Finalize
  void DELPHI_2002_069_CONF_603::finalize() {
    normalize(_histXbprim);
    normalize(_histXbweak);
  }

}
