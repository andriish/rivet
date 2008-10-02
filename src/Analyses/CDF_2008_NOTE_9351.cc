// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_NOTE_9351.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {


  // Book histograms
  void CDF_2008_NOTE_9351::init() {
    _hist_tnchg      = bookProfile1D( 1, 1, 1, "Toward Region Charged Particle Density");
    _hist_pnchg      = bookProfile1D( 2, 1, 1, "Transverse Region Charged Particle Density");
    _hist_pmaxnchg   = bookProfile1D( 3, 1, 1, "TransMAX Region Charged Particle Density");
    _hist_pminnchg   = bookProfile1D( 4, 1, 1, "TransMIN Region Charged Particle Density");
    _hist_pdifnchg   = bookProfile1D( 5, 1, 1, "TransDIF Region Charged Particle Density");
    _hist_anchg      = bookProfile1D( 6, 1, 1, "Away Region Charged Particle Density");
    _hist_tcptsum    = bookProfile1D( 7, 1, 1, "Toward Region Charged pT Sum Density");
    _hist_pcptsum    = bookProfile1D( 8, 1, 1, "Transverse Region Charged pT Sum Density");
    _hist_pmaxcptsum = bookProfile1D( 9, 1, 1, "TransMAX Region Charged pT Sum Density");
    _hist_pmincptsum = bookProfile1D(10, 1, 1, "TransMIN Region Charged pT Sum Density");
    _hist_pdifcptsum = bookProfile1D(11, 1, 1, "TransDIF Region Charged pT Sum Density");
    _hist_acptsum    = bookProfile1D(12, 1, 1, "Away Region Charged pT Sum Density");
    _hist_tcptave    = bookProfile1D(13, 1, 1, "Toward Region Charged pT Average");
    _hist_pcptave    = bookProfile1D(14, 1, 1, "Transverse Region Charged pT Average");
    _hist_acptave    = bookProfile1D(15, 1, 1, "Away Region Charged pT Average");
    _hist_tcptmax    = bookProfile1D(16, 1, 1, "Toward Region Charged pT Maximum");
    _hist_pcptmax    = bookProfile1D(17, 1, 1, "Transverse Region Charged pT Maximum");
    _hist_acptmax    = bookProfile1D(18, 1, 1, "Away Region Charged pT Maximum");
    _hist_zptvsnchg  = bookProfile1D(19, 1, 1, "Average Lepton Pair pT versus Charged Multiplicity");
    _hist_cptavevsnchg = bookProfile1D(20, 1, 1, "Average Charged pT versus Charged Multiplicity");
    _hist_cptavevsnchgsmallzpt = bookProfile1D(21, 1, 1, "Average Charged pT versus Charged Multiplicity, pT(Z) less than 10 GeV");
  }


  // Do the analysis
  void CDF_2008_NOTE_9351::analyze(const Event& e) {
    Log log = getLog();

    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const size_t numParticles = fs.particles().size();

    // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
    if (numParticles < 1) {
      getLog() << Log::DEBUG << "Failed multiplicity cut" << endl;
      vetoEvent(e);
    }

    // Get the event weight
    const double weight = e.weight();

    // Get the leptons
    const ParticleVector& leptons = applyProjection<ChargedLeptons>(e, "CL").chargedLeptons();

    // We want exactly two leptons of the same flavour.
    getLog() << Log::DEBUG << "lepton multiplicity = " << leptons.size() << endl;
    if (leptons.size() != 2 || leptons[0].getPdgId() != -leptons[1].getPdgId() )
      vetoEvent(e);

    // Lepton pT > 20 GeV
    if (leptons[0].momentum().pT() <= 20 || leptons[1].momentum().pT() <= 20)
      vetoEvent(e);

    FourVector dilepton = leptons[0].momentum() + leptons[1].momentum();

    // Lepton pair should have an invariant mass between 70 and 110 and |eta|<6
    if (mass(dilepton) < 70 || mass(dilepton) > 110 || fabs(pseudorapidity(dilepton)) >= 6)
      vetoEvent(e);
    getLog() << Log::DEBUG << "dilepton mass = " << mass(dilepton) << endl;
    getLog() << Log::DEBUG << "dilepton pT   = " << pT(dilepton) << endl;


    // Calculate the observables
    size_t   numToward(0),     numTrans1(0),     numTrans2(0),     numAway(0);
    double ptSumToward(0.0), ptSumTrans1(0.0), ptSumTrans2(0.0), ptSumAway(0.0);
    double ptMaxToward(0.0), ptMaxTrans1(0.0), ptMaxTrans2(0.0), ptMaxAway(0.0);
    const double phiZ = azimuthalAngle(dilepton);
    const double pTZ  = pT(dilepton);
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      // don't use the leptons:
      if (abs(p->getPdgId()) < 20)
        continue;

      const double deltaPhi = delta_phi(p->momentum().azimuthalAngle(), phiZ);
      const double pT = p->momentum().pT();
      double rotatedphi = p->momentum().azimuthalAngle() - phiZ;
      while (rotatedphi < 0) rotatedphi += 2*PI;

      if (deltaPhi < PI/3.0) {
        ptSumToward += pT;
        ++numToward;
        if (pT > ptMaxToward)
          ptMaxToward = pT;
      } else if (deltaPhi < 2*PI/3.0) {
        if (rotatedphi <= PI) {
          ptSumTrans1 += pT;
          ++numTrans1;
          if (pT > ptMaxTrans1)
            ptMaxTrans1 = pT;
        }
        else {
          ptSumTrans2 += pT;
          ++numTrans2;
          if (pT > ptMaxTrans2)
            ptMaxTrans2 = pT;
        }
      } else {
        ptSumAway += pT;
        ++numAway;
        if (pT > ptMaxAway)
          ptMaxAway = pT;
      }
      // We need to subtract the two leptons from the number of particles to get the correct multiplicity
      _hist_cptavevsnchg->fill(numParticles-2, pT, weight);
      if (pTZ < 10)
        _hist_cptavevsnchgsmallzpt->fill(numParticles-2, pT, weight);
    }

    // Fill the histograms
    _hist_tnchg->fill(pTZ, numToward/(4*PI/3), weight);
    _hist_pnchg->fill(pTZ, (numTrans1+numTrans2)/(4*PI/3), weight);
    _hist_pmaxnchg->fill(pTZ, (numTrans1>numTrans2 ? numTrans1 : numTrans2)/(2*PI/3), weight);
    _hist_pminnchg->fill(pTZ, (numTrans1<numTrans2 ? numTrans1 : numTrans2)/(2*PI/3), weight);
    _hist_pdifnchg->fill(pTZ, abs(numTrans1-numTrans2)/(2*PI/3), weight);
    _hist_anchg->fill(pTZ, numAway/(4*PI/3), weight);

    _hist_tcptsum->fill(pTZ, ptSumToward/(4*PI/3), weight);
    _hist_pcptsum->fill(pTZ, (ptSumTrans1+ptSumTrans2)/(4*PI/3), weight);
    _hist_pmaxcptsum->fill(pTZ, (ptSumTrans1>ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/(2*PI/3), weight);
    _hist_pmincptsum->fill(pTZ, (ptSumTrans1<ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/(2*PI/3), weight);
    _hist_pdifcptsum->fill(pTZ, fabs(ptSumTrans1-ptSumTrans2)/(2*PI/3), weight);
    _hist_acptsum->fill(pTZ, ptSumAway/(4*PI/3), weight);

    if (numToward > 0) {
      _hist_tcptave->fill(pTZ, ptSumToward/numToward, weight);
      _hist_tcptmax->fill(pTZ, ptMaxToward, weight);
    }
    if ((numTrans1+numTrans2) > 0) {
      _hist_pcptave->fill(pTZ, (ptSumTrans1+ptSumTrans2)/(numTrans1+numTrans2), weight);
      _hist_pcptmax->fill(pTZ, (ptMaxTrans1 > ptMaxTrans2 ? ptMaxTrans1 : ptMaxTrans2), weight);
    }
    if (numAway > 0) {
      _hist_acptave->fill(pTZ, ptSumAway/numAway, weight);
      _hist_acptmax->fill(pTZ, ptMaxAway, weight);
    }

    // We need to subtract the two leptons from the number of particles to get the correct multiplicity
    _hist_zptvsnchg->fill(numParticles-2, pTZ, weight);
  }


  void CDF_2008_NOTE_9351::finalize() { 
  }


}
