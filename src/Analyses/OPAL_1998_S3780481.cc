// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/OPAL_1998_S3780481.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"

namespace Rivet {


  void OPAL_1998_S3780481::analyze(const Event& e) {
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
    _weightedTotalPartNum += numParticles * weight;

    // Get beams and average beam momentum
    const ParticlePair& beams = applyProjection<Beam>(e, "Beams").beams();
    const double meanBeamMom = ( beams.first.momentum().vector3().mod() + 
                                 beams.second.momentum().vector3().mod() ) / 2.0;
    getLog() << Log::DEBUG << "Avg beam momentum = " << meanBeamMom << endl;

    int flavour = 0;
    const InitialQuarks& iqf = applyProjection<InitialQuarks>(e, "IQF");

    // If we only have two quarks (qqbar), just take the flavour.
    // If we have more than two quarks, look for the highest energetic q-qbar pair.
    if (iqf.particles().size() == 2) {
      flavour = abs( iqf.particles().front().pdgId() );
    }
    else {
      std::map<int, double> quarkmap;
      foreach (const Particle& p, iqf.particles()) {
        if (quarkmap[p.pdgId()] < p.momentum().E()) {
          quarkmap[p.pdgId()] = p.momentum().E();
        }
      }
      double maxenergy = 0.;
      for (int i = 1; i <= 5; ++i) {
        if (quarkmap[i]+quarkmap[-i] > maxenergy) {
          flavour = i;
        }
      }
    }

    switch (flavour) {
      case 1:
      case 2:
      case 3:
        _SumOfudsWeights += weight;
        break;
      case 4:
        _SumOfcWeights += weight;
        break;
      case 5:
        _SumOfbWeights += weight;
        break;
    }

    foreach (const Particle& p, fs.particles()) {
      const double xp = p.momentum().vector3().mod()/meanBeamMom;
      const double logxp = -std::log(xp);
      _histXpall->fill(xp, weight);
      _histLogXpall->fill(logxp, weight);
      _histMultiChargedall->fill(_histMultiChargedall->binMean(0), weight);
      switch (flavour) {
        /// @todo Use PDG code enums
        case 1:
        case 2:
        case 3:
          _histXpuds->fill(xp, weight);
          _histLogXpuds->fill(logxp, weight);
          _histMultiChargeduds->fill(_histMultiChargeduds->binMean(0), weight);
          break;
        case 4:
          _histXpc->fill(xp, weight);
          _histLogXpc->fill(logxp, weight);
          _histMultiChargedc->fill(_histMultiChargedc->binMean(0), weight);
          break;
        case 5:
          _histXpb->fill(xp, weight);
          _histLogXpb->fill(logxp, weight);
          _histMultiChargedb->fill(_histMultiChargedb->binMean(0), weight);
          break;
      }
    }

  }


  void OPAL_1998_S3780481::init() {
    _histXpuds           = bookHistogram1D(1, 1, 1, "$uds$ events scaled momentum");
    _histXpc             = bookHistogram1D(2, 1, 1, "$c$ events scaled momentum");
    _histXpb             = bookHistogram1D(3, 1, 1, "$b$ events scaled momentum");
    _histXpall           = bookHistogram1D(4, 1, 1, "All events scaled momentum");
    _histLogXpuds        = bookHistogram1D(5, 1, 1, "$uds$ events $\\ln(1/x_p)$");
    _histLogXpc          = bookHistogram1D(6, 1, 1, "$c$ events $\\ln(1/x_p)$");
    _histLogXpb          = bookHistogram1D(7, 1, 1, "$b$ events $\\ln(1/x_p)$");
    _histLogXpall        = bookHistogram1D(8, 1, 1, "All events $\\ln(1/x_p)$");
    _histMultiChargeduds = bookHistogram1D(9, 1, 1, "$uds$ events mean charged multiplicity");
    _histMultiChargedc   = bookHistogram1D(9, 1, 2, "$c$ events mean charged multiplicity");
    _histMultiChargedb   = bookHistogram1D(9, 1, 3, "$b$ events mean charged multiplicity");
    _histMultiChargedall = bookHistogram1D(9, 1, 4, "All events mean charged multiplicity");
  }


  // Finalize
  void OPAL_1998_S3780481::finalize() {
    const double avgNumParts = _weightedTotalPartNum / sumOfWeights();
    normalize(_histXpuds    , avgNumParts);
    normalize(_histXpc      , avgNumParts);
    normalize(_histXpb      , avgNumParts);
    normalize(_histXpall    , avgNumParts);
    normalize(_histLogXpuds , avgNumParts);
    normalize(_histLogXpc   , avgNumParts);
    normalize(_histLogXpb   , avgNumParts);
    normalize(_histLogXpall , avgNumParts);

    scale(_histMultiChargeduds, 1.0/_SumOfudsWeights);
    scale(_histMultiChargedc  , 1.0/_SumOfcWeights);
    scale(_histMultiChargedb  , 1.0/_SumOfbWeights);
    scale(_histMultiChargedall, 1.0/sumOfWeights());
  }


}
