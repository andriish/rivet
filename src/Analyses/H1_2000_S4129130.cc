// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/H1_2000_S4129130.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  // Constructor
  H1_2000_S4129130::H1_2000_S4129130() 
    : Analysis("H1_2000_S4129130")
  {
    setBeams(ELECTRON, PROTON);
    addProjection(DISLepton(), "Lepton");
    addProjection(DISKinematics(), "Kinematics");
    addProjection(FinalState(), "FS");
  }


  void H1_2000_S4129130::analyze(const Event& event) {
    // Get the projections
    const FinalState & fs = applyProjection<FinalState>(event, "FS");
    const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");
    const DISLepton    & dl = applyProjection<DISLepton>(event,"Lepton");

    // Get the DIS kinematics
    double q2  = dk.Q2();
    double x   = dk.x();
    double y   = dk.y();
    double w2  = dk.W2();

    // Momentum of the scattered lepton
    FourMomentum leptonMom = dl.out().momentum();
    // pT energy and angle
    double enel = leptonMom.E();
    double thel = 180.-leptonMom.angle(dl.in().momentum())/degree;

    // Extract the particles other than the lepton
    ParticleVector particles;
    particles.reserve(fs.particles().size());
    const GenParticle& dislepGP = dl.out().genParticle();
    for (ParticleVector::const_iterator p = fs.particles().begin();
         p != fs.particles().end(); ++p) {
      const GenParticle& loopGP = p->genParticle(); 
      if (&loopGP == &dislepGP) continue;
      particles.push_back(*p);
    }

    // Cut on the forward energy
    double efwd = 0.;
    foreach (const Particle& p, particles) {
      double th = 180.-p.momentum().angle(dl.in().momentum())/degree;
//      double th = beamAngle(p.momentum(),order);
      if(th > 4.4 && th < 15.0) efwd += p.momentum().E();
    }
    // There are four possible selections for events
    bool evcut[4];
    // Low  Q2 selection a
    /// @todo Units and inRange
    evcut[0] = enel/GeV > 12. && w2 >= 4400. && efwd/GeV > 0.5 && 
               thel > 157. && thel < 176.0;
    // Low  Q2 selection b
    /// @todo Units and inRange
    evcut[1] = enel/GeV > 12. && y > 0.3 && y < 0.5;
    // High Q2 selection a
    /// @todo Units and inRange
    evcut[2] = thel > 12. && thel < 150.0 && y > 0.05 && y < 0.6 && 
               w2 >= 4400. && efwd > 0.5;
    // High Q2 selection b
    /// @todo Units and inRange
    evcut[3] = thel > 12. && thel < 150.0 && y > 0.05 && y < 0.6 &&
               w2 > 27110. && w2 < 45182.;

    // Veto if fails all cuts
    if (! (evcut[0] || evcut[1] || evcut[2] || evcut[3]) ) {
      vetoEvent;
    }

    // Find the bins
    int bin[4] = {-1,-1,-1,-1};
    // For the low Q2 selection a)
    /// @todo Units
    if (q2 > 2.5 && q2 <= 5.) {
      if (x > 0.00005 && x <= 0.0001 ) bin[0] = 0;
      if (x > 0.0001  && x <= 0.0002 ) bin[0] = 1;
      if (x > 0.0002  && x <= 0.00035) bin[0] = 2;
      if (x > 0.00035 && x <= 0.0010 ) bin[0] = 3;
    }
    /// @todo Units
    else if(q2 > 5. && q2 <= 10.) {
      if (x > 0.0001  && x <= 0.0002 ) bin[0] = 4;
      if (x > 0.0002  && x <= 0.00035) bin[0] = 5;
      if (x > 0.00035 && x <= 0.0007 ) bin[0] = 6;
      if (x > 0.0007  && x <= 0.0020 ) bin[0] = 7;
    }
    /// @todo Units
    else if(q2 > 10. && q2 <= 20.) {
      if (x > 0.0002 && x <= 0.0005) bin[0] = 8;
      if (x > 0.0005 && x <= 0.0008) bin[0] = 9;
      if (x > 0.0008 && x <= 0.0015) bin[0] = 10;
      if (x > 0.0015 && x <= 0.040 ) bin[0] = 11;
    }
    /// @todo Units
    else if(q2 > 20. && q2 <= 50.) {
      if (x > 0.0005 && x <= 0.0014) bin[0] = 12;
      if (x > 0.0014 && x <= 0.0030) bin[0] = 13;
      if (x > 0.0030 && x <= 0.0100) bin[0] = 14;
    }
    /// @todo Units
    else if (q2 > 50. && q2 <= 100.) {
      if (x >0.0008 && x <= 0.0030) bin[0] = 15;
      if (x >0.0030 && x <= 0.0200) bin[0] = 16;
    }
    /// @todo Huh?
    evcut[0] &= bin[0] >= 0;
    // For the low Q2 selection b)
    if (q2 > 2.5 && q2 <= 5.  ) bin[1] = 0;
    if (q2 > 5.  && q2 <= 10. ) bin[1] = 1;
    if (q2 > 10. && q2 <= 20. ) bin[1] = 2;
    if (q2 > 20. && q2 <= 50. ) bin[1] = 3;
    if (q2 > 50. && q2 <= 100.) bin[1] = 4;
    evcut[1] &= bin[1] >= 0;
    // for the high Q2 selection a)
    /// @todo Units
    if (q2 > 100. && q2 <= 400.) {
      if (x > 0.00251 && x <= 0.00631) bin[2] = 0;
      if (x > 0.00631 && x <= 0.0158 ) bin[2] = 1;
      if (x > 0.0158  && x <= 0.0398 ) bin[2] = 2;
    }
    /// @todo Units
    else if (q2 > 400 && q2 <= 1100.) {
      if (x > 0.00631 && x <= 0.0158 ) bin[2] = 3;
      if (x > 0.0158  && x <= 0.0398 ) bin[2] = 4;
      if (x > 0.0398  && x <= 1.     ) bin[2] = 5;
    }
    /// @todo Units
    else if (q2 > 1100. && q2 <= 100000.) {
      if (x > 0. && x <= 1.) bin[2] = 6;
    }
    evcut[2] &= bin[2] >= 0;
    // for the high Q2 selection b)
    /// @todo Units
    if      (q2 > 100. && q2 <= 220.) bin[3] = 0;
    else if (q2 > 220. && q2 <= 400.) bin[3] = 1;
    else if (q2 > 400.              ) bin[3] = 2;
    evcut[3] &= bin[3] >= 0;

    // Veto if fails all cuts after bin selection
    if (! (evcut[0] || evcut[1] || evcut[2] || evcut[3]));

    // Increment the count for normalisation
    const double weight = event.weight();
    if (evcut[0]) _weightETLowQa [bin[0]] += weight;
    if (evcut[1]) _weightETLowQb [bin[1]] += weight;
    if (evcut[2]) _weightETHighQa[bin[2]] += weight;
    if (evcut[3]) _weightETHighQb[bin[3]] += weight;

    // Boost to hadronicCM
    const LorentzTransform hcmboost = dk.boostHCM();

    // Loop over the particles
    double etcent = 0;
    double etfrag = 0;
    foreach (const Particle& p, particles) {
      // Boost momentum to CMS
      const FourMomentum hcmMom = hcmboost.transform(p.momentum());
      double et = fabs(Et(hcmMom));
      double eta = -hcmMom.pseudorapidity();
      // Averages in central and forward region
      if (fabs(eta) < .5 ) etcent += et;
      if (eta > 2 && eta <= 3.) etfrag += et;
      // Histograms of Et flow
      if (evcut[0]) _histETLowQa [bin[0]]->fill(eta, et*weight);
      if (evcut[1]) _histETLowQb [bin[1]]->fill(eta, et*weight);
      if (evcut[2]) _histETHighQa[bin[2]]->fill(eta, et*weight);
      if (evcut[3]) _histETHighQb[bin[3]]->fill(eta, et*weight);
    }
    // Fill histograms for the average quantities
    if (evcut[1] || evcut[3]) {
      _histAverETCentral->fill(q2, etcent*weight,weight);
      _histAverETFrag   ->fill(q2, etfrag*weight,weight);
    }
  }


  void H1_2000_S4129130::init() {

    string t = "Transverse energy flow for ";
    IHistogram1D* h = 0;

    string xlabel = "$\\eta$";
    /// @todo What is "N"?
    string ylabel = "$1/N \\, \\mathrm{d}{E_\\perp}/\\mathrm{d}{\\eta}$ / GeV";

    const string xt = "\\langle x \\rangle";
    const string Q2t = "\\langle Q^2 \\rangle";

    // Histograms and weight vectors for low Q^2 a
    _histETLowQa.reserve(17);
    _weightETLowQa.reserve(17);
    for (size_t ix = 0; ix < 17; ++ix) {
      string title = t + "$" + xt;
      if      (ix ==  0) title += " = 0.08\\cdot 10^{-3}, " + Q2t + " =  3.2";
      else if (ix ==  1) title += " = 0.14\\cdot 10^{-3}, " + Q2t + " =  3.8";
      else if (ix ==  2) title += " = 0.26\\cdot 10^{-3}, " + Q2t + " =  3.9";
      else if (ix ==  3) title += " = 0.57\\cdot 10^{-3}, " + Q2t + " =  4.2";
      else if (ix ==  4) title += " = 0.16\\cdot 10^{-3}, " + Q2t + " =  6.3";
      else if (ix ==  5) title += " = 0.27\\cdot 10^{-3}, " + Q2t + " =  7.0";
      else if (ix ==  6) title += " = 0.50\\cdot 10^{-3}, " + Q2t + " =  7.0";
      else if (ix ==  7) title += " = 1.10\\cdot 10^{-3}, " + Q2t + " =  7.3";
      else if (ix ==  8) title += " = 0.36\\cdot 10^{-3}, " + Q2t + " = 13.1";
      else if (ix ==  9) title += " = 0.63\\cdot 10^{-3}, " + Q2t + " = 14.1";
      else if (ix == 10) title += " = 1.10\\cdot 10^{-3}, " + Q2t + " = 14.1";
      else if (ix == 11) title += " = 2.30\\cdot 10^{-3}, " + Q2t + " = 14.9";
      else if (ix == 12) title += " = 0.93\\cdot 10^{-3}, " + Q2t + " = 28.8";
      else if (ix == 13) title += " = 2.10\\cdot 10^{-3}, " + Q2t + " = 31.2";
      else if (ix == 14) title += " = 4.70\\cdot 10^{-3}, " + Q2t + " = 33.2";
      else if (ix == 15) title += " = 2.00\\cdot 10^{-3}, " + Q2t + " = 59.4";
      else if (ix == 16) title += " = 7.00\\cdot 10^{-3}, " + Q2t + " = 70.2";
      title += " \\text{ GeV}^2$";
      h = bookHistogram1D(ix+1, 1, 1, title, xlabel, ylabel);
      _histETLowQa.push_back(h);
      _weightETLowQa.push_back(0.);
    }

    // Histograms and weight vectors for high Q^2 a
    _histETHighQa.reserve(7);
    _weightETHighQa.reserve(7);
    for (size_t ix = 0; ix < 7; ++ix) {
      string title = t + "$" + xt;
      if      (ix == 0) title += " = 0.0043, " + Q2t +	" = 175";
      else if (ix == 1) title += " = 0.01, "   + Q2t +	" = 253";
      else if (ix == 2) title += " = 0.026, "  + Q2t +	" = 283";
      else if (ix == 3) title += " = 0.012, "  + Q2t +	" = 511";
      else if (ix == 4) title += " = 0.026, "  + Q2t +	" = 617";
      else if (ix == 5) title += " = 0.076, "  + Q2t +	" = 682";
      else if (ix == 6) title += " = 0.11, "   + Q2t + " = 2200";
      title += " \\text{ GeV}^2$";
      h = bookHistogram1D(ix+18, 1, 1, title, xlabel, ylabel);
      _histETHighQa.push_back(h);
      _weightETHighQa.push_back(0.);
    }

    // Histograms and weight vectors for low Q^2 b
    _histETLowQb.reserve(5);
    _weightETLowQb.reserve(5);
    for (size_t ix = 0; ix < 5; ++ix) {
      string title = t + "$" + Q2t;
      if      (ix == 0) title +=  " = 2.5-5";
      else if (ix == 1) title +=   " = 5-10";
      else if (ix == 2) title +=  " = 10-20";
      else if (ix == 3) title +=  " = 20-50";
      else if (ix == 4) title += " = 50-100";
      title += " \\text{ GeV}^2$";
      h = bookHistogram1D(ix+25, 1, 1, title, xlabel, ylabel);
      _histETLowQb.push_back(h);
      _weightETLowQb.push_back(0.);
    }

    // Histograms and weight vectors for high Q^2 b
    _histETHighQb.reserve(3);
    _weightETHighQb.reserve(3);
    for (size_t ix = 0; ix < 3; ++ix) {
      string title = t + "$" + Q2t;
      if      (ix == 0) title += " = 100-220";
      else if (ix == 1) title += " = 220-400";
      else if (ix == 1) title += " > 400";
      title += " \\text{ GeV}^2$";
      h = bookHistogram1D(30+ix, 1, 1, title, xlabel, ylabel);
      _histETHighQb.push_back(h);
      _weightETHighQb.push_back(0.0);
    }

    // Histograms for the averages
    xlabel = "$Q^2$ / $\\text{GeV}^2$";
    ylabel = "$\\langle E_\\perp \\rangle$ / GeV";
    _histAverETCentral = 
      bookProfile1D(33,  1, 1, 
                    "Average $E_\\perp$ in the central region", xlabel, ylabel);
    _histAverETFrag =
      bookProfile1D(34,  1, 1, 
                    "Average $E_\\perp$ in the forward region", xlabel, ylabel);
  }


  // Finalize
  void H1_2000_S4129130::finalize() { 
    // Normalization of the Et distributions
    for (size_t ix=0; ix<17; ++ix) {
      scale(_histETLowQa[ix], 1./_weightETLowQa[ix]);
    }
    for(size_t ix=0; ix<7; ++ix) {
      scale(_histETHighQa[ix], 1./_weightETHighQa[ix]);
    }
    for(size_t ix=0; ix<5; ++ix) {
      scale(_histETLowQb[ix], 1./_weightETLowQb[ix]);
    }
    for(size_t ix=0; ix<3; ++ix) {
      scale(_histETHighQb[ix], 1./_weightETHighQb[ix]);
    }
  }


}
