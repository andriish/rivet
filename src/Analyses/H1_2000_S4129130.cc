// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/H1_2000_S4129130.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"

namespace Rivet {

  void H1_2000_S4129130::analyze(const Event& event) {
    // get the projections
    const FinalState & fs = applyProjection<FinalState>(event, "FS");
    const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");
    const DISLepton    & dl = applyProjection<DISLepton>(event,"Lepton");
    // require outgoing lepton same type as incoming
    if(dl.in().getPdgId()!=dl.out().getPdgId()) vetoEvent(event);
    // get the DIS kinematics
    double q2  = dk.Q2();
    double x   = dk.x();
    double y   = dk.y();
    double w2  = dk.W2();
    // whether electron in + (false) or - (true) z
    bool order = dl.in().momentum().z()<0.;
    // momentum of the scattered lepton
    FourMomentum leptonMom = dl.out().momentum();
    // pt energy and angle
    double enel = leptonMom.E();
    double thel = beamAngle(leptonMom,order);
    // extract the particles other than the lepton
    ParticleVector particles;
    particles.reserve(fs.particles().size());
    const GenParticle& dislepGP = dl.out().getHepMCParticle();
    for (ParticleVector::const_iterator p = fs.particles().begin();
 	 p != fs.particles().end(); ++p) {
      const GenParticle& loopGP = p->getHepMCParticle(); 
      if (&loopGP == &dislepGP) continue;
      particles.push_back(*p);
    }
    // cut on the forward energy
    double efwd = 0.;
    for (ParticleVector::const_iterator p = particles.begin(); 
	 p != particles.end(); ++p) {
      double th = beamAngle(p->momentum(),order);
      if(th>4.4&&th<15.) efwd += p->momentum().E();
    }
    // there are four possible selections for events
    bool evcut[4];
    // low  Q2 selection a
    evcut[0] = enel>12. && w2>=4400. && efwd>0.5 && 
               thel>157. && thel<176.0;
    // low  Q2 selection b
    evcut[1] = enel>12. && y>0.3 && y<0.5;
    // high Q2 selection a
    evcut[2] = thel>12. && thel<150.0 && y>0.05 && y<0.6 && 
               w2>=4400. && efwd>0.5;
    // high Q2 selection b
    evcut[3] = thel>12. && thel<150.0 && y>0.05 && y<0.6 &&
               w2>27110. && w2<45182.;
    // veto if fails all cuts
    if(!evcut[0]&&!evcut[1]&&!evcut[2]&&!evcut[3])
      vetoEvent(event);
    // find the bins
    int bin[4]={-1,-1,-1,-1};
    // for the low Q2 selection a)
    if(q2>2.5 && q2<=5.) {
      if(x>0.00005 && x<=0.0001 ) bin[0] = 0;
      if(x>0.0001  && x<=0.0002 ) bin[0] = 1;
      if(x>0.0002  && x<=0.00035) bin[0] = 2;
      if(x>0.00035 && x<=0.0010 ) bin[0] = 3;
    }
    else if(q2>5. && q2<=10.) {
      if(x>0.0001  && x<=0.0002 ) bin[0] = 4;
      if(x>0.0002  && x<=0.00035) bin[0] = 5;
      if(x>0.00035 && x<=0.0007 ) bin[0] = 6;
      if(x>0.0007  && x<=0.0020 ) bin[0] = 7;
    }
    else if(q2>10. && q2<=20.) {
      if(x>0.0002 && x<=0.0005) bin[0] = 8;
      if(x>0.0005 && x<=0.0008) bin[0] = 9;
      if(x>0.0008 && x<=0.0015) bin[0] = 10;
      if(x>0.0015 && x<=0.040 ) bin[0] = 11;
    }
    else if(q2>20. && q2<=50.) {
      if(x>0.0005 && x<=0.0014) bin[0] = 12;
      if(x>0.0014 && x<=0.0030) bin[0] = 13;
      if(x>0.0030 && x<=0.0100) bin[0] = 14;
    }
    else if(q2>50. && q2<=100.) {
      if(x>0.0008 && x<=0.0030) bin[0] = 15;
      if(x>0.0030 && x<=0.0200) bin[0] = 16;
    }
    evcut[0] &= bin[0]>=0;
    // for the low Q2 selection b)
    if(q2>2.5 && q2<=5.  ) bin[1] = 0;
    if(q2>5.  && q2<=10. ) bin[1] = 1;
    if(q2>10. && q2<=20. ) bin[1] = 2;
    if(q2>20. && q2<=50. ) bin[1] = 3;
    if(q2>50. && q2<=100.) bin[1] = 4;
    evcut[1] &= bin[1]>=0;
    // for the high Q2 selection a)
    if(q2>100. && q2<=400.) {
      if(x>0.00251 && x<=0.00631) bin[2] = 0;
      if(x>0.00631 && x<=0.0158 ) bin[2] = 1;
      if(x>0.0158  && x<=0.0398 ) bin[2] = 2;
    }
    else if(q2>400 && q2<=1100.) {
      if(x>0.00631 && x<=0.0158 ) bin[2] = 3;
      if(x>0.0158  && x<=0.0398 ) bin[2] = 4;
      if(x>0.0398  && x<=1.     ) bin[2] = 5;
    }
    else if(q2>1100. && q2<=100000.) {
      if(x>0.      && x<=1.     ) bin[2] = 6;
    }
    evcut[2] &= bin[2]>=0;
    // for the high Q2 selection b)
    if     (q2>100. && q2<=220.) bin[3] = 0;
    else if(q2>220. && q2<=400.) bin[3] = 1;
    else if(q2>400.            ) bin[3] = 2;
    evcut[3] &= bin[3]>=0;
    // veto if fails all cuts after bin selection
    if(!evcut[0]&&!evcut[1]&&!evcut[2]&&!evcut[3]);
    // increment the count for normalisation
    const double weight = event.weight();
    if(evcut[0]) _weightETLowQa [bin[0]] += weight;
    if(evcut[1]) _weightETLowQb [bin[1]] += weight;
    if(evcut[2]) _weightETHighQa[bin[2]] += weight;
    if(evcut[3]) _weightETHighQb[bin[3]] += weight;
    // boost to hadronicCM
    const LorentzTransform hcmboost = dk.boostHCM();
    // loop over the particles
    double etcent=0;
    double etfrag=0;
    for(ParticleVector::const_iterator p = particles.begin();
 	p!= particles.end(); ++p ) {
      // boost momentum to CMS
      const FourMomentum hcmMom = hcmboost.transform(p->momentum());
      double et = fabs(Et(hcmMom));
      double eta = -hcmMom.pseudorapidity();
      // averages in central and forward region
      if(fabs(eta)< .5 )  etcent += et;
      if(eta>2&&eta<=3.) etfrag += et;
      // histograms of et flow
      if(evcut[0]) _histETLowQa [bin[0]]->fill(eta,et*weight);
      if(evcut[1]) _histETLowQb [bin[1]]->fill(eta,et*weight);
      if(evcut[2]) _histETHighQa[bin[2]]->fill(eta,et*weight);
      if(evcut[3]) _histETHighQb[bin[3]]->fill(eta,et*weight);
    }
    // fill histograms for the average quantities
    if(evcut[1]||evcut[3]) {
      _histAverETCentral->fill(q2,etcent*weight,weight);
      _histAverETFrag   ->fill(q2,etfrag*weight,weight);
    }
  }

  void H1_2000_S4129130::init() {
    // histograms and weight vectors for low  Q2 a
    _histETLowQa.reserve(17);
    _weightETLowQa.reserve(17);
    for (size_t ix = 0; ix < 17; ++ix) {
      string title = "Transverse energy flow for ";
      if (ix == 0)    title += "$\\langle x\rangle = 0.08\\cdot 10^{-3}, \\langle Q^2\rangle =  3.2 \\text{ GeV}^2$";
      else if (ix == 1) title += "$\\langle x\rangle = 0.14\\cdot 10^{-3}, \\langle Q^2\rangle =  3.8 \\text{ GeV}^2$";
      else if (ix == 2) title += "$\\langle x\rangle = 0.26\\cdot 10^{-3}, \\langle Q^2\rangle =  3.9 \\text{ GeV}^2$";
      else if (ix == 3) title += "$\\langle x\rangle = 0.57\\cdot 10^{-3}, \\langle Q^2\rangle =  4.2 \\text{ GeV}^2$";
      else if (ix == 4) title += "$\\langle x\rangle = 0.16\\cdot 10^{-3}, \\langle Q^2\rangle =  6.3 \\text{ GeV}^2$";
      else if (ix == 5) title += "$\\langle x\rangle = 0.27\\cdot 10^{-3}, \\langle Q^2\rangle =  7.0 \\text{ GeV}^2$";
      else if (ix == 6) title += "$\\langle x\rangle = 0.50\\cdot 10^{-3}, \\langle Q^2\rangle =  7.0 \\text{ GeV}^2$";
      else if (ix == 7) title += "$\\langle x\rangle = 1.10\\cdot 10^{-3}, \\langle Q^2\rangle =  7.3 \\text{ GeV}^2$";
      else if (ix == 8) title += "$\\langle x\rangle = 0.36\\cdot 10^{-3}, \\langle Q^2\rangle = 13.1 \\text{ GeV}^2$";
      else if (ix == 9) title += "$\\langle x\rangle = 0.63\\cdot 10^{-3}, \\langle Q^2\rangle = 14.1 \\text{ GeV}^2$";
      else if (ix ==10) title += "$\\langle x\rangle = 1.10\\cdot 10^{-3}, \\langle Q^2\rangle = 14.1 \\text{ GeV}^2$";
      else if (ix ==11) title += "$\\langle x\rangle = 2.30\\cdot 10^{-3}, \\langle Q^2\rangle = 14.9 \\text{ GeV}^2$";
      else if (ix ==12) title += "$\\langle x\rangle = 0.93\\cdot 10^{-3}, \\langle Q^2\rangle = 28.8 \\text{ GeV}^2$";
      else if (ix ==13) title += "$\\langle x\rangle = 2.10\\cdot 10^{-3}, \\langle Q^2\rangle = 31.2 \\text{ GeV}^2$";
      else if (ix ==14) title += "$\\langle x\rangle = 4.70\\cdot 10^{-3}, \\langle Q^2\rangle = 33.2 \\text{ GeV}^2$";
      else if (ix ==15) title += "$\\langle x\rangle = 2.00\\cdot 10^{-3}, \\langle Q^2\rangle = 59.4 \\text{ GeV}^2$";
      else if (ix ==16) title += "$\\langle x\rangle = 7.00\\cdot 10^{-3}, \\langle Q^2\rangle = 70.2 \\text{ GeV}^2$";
      _histETLowQa.push_back(bookHistogram1D(ix+1,  1, 1, title));
      _weightETLowQa.push_back(0.);
    }
    // histograms and weight vectors for high Q2 a
    _histETHighQa.reserve(7);
    _weightETHighQa.reserve(7);
    for (size_t ix = 0; ix < 7; ++ix) {
      string title = "Transverse energy flow for ";
      if (ix == 0)      title += "$\\langle x\\rangle = 0.0043, \\langle Q^2\\rangle =  175 \\text{ GeV}^2$";
      else if (ix == 1) title += "$\\langle x\\rangle = 0.01,   \\langle Q^2\\rangle =  253 \\text{ GeV}^2$";
      else if (ix == 2) title += "$\\langle x\\rangle = 0.026,  \\langle Q^2\\rangle =  283 \\text{ GeV}^2$";
      else if (ix == 3) title += "$\\langle x\\rangle = 0.012,  \\langle Q^2\\rangle =  511 \\text{ GeV}^2$";
      else if (ix == 4) title += "$\\langle x\\rangle = 0.026,  \\langle Q^2\\rangle =  617 \\text{ GeV}^2$";
      else if (ix == 5) title += "$\\langle x\\rangle = 0.076,  \\langle Q^2\\rangle =  682 \\text{ GeV}^2$";
      else if (ix == 6) title += "$\\langle x\\rangle = 0.11,   \\langle Q^2\\rangle = 2200 \\text{ GeV}^2$";
      _histETHighQa.push_back(bookHistogram1D(ix+18,  1, 1, title));
      _weightETHighQa.push_back(0.);
    }
    // histograms and weight vectors for low  Q2 b
    _histETLowQb.reserve(5);
    _weightETLowQb.reserve(5);
    for (size_t ix = 0; ix < 5; ++ix) {
      string title = "Transverse energy flow for ";
      if (ix == 0)      title += "$\\langle Q^2\\rangle = 2.5-5  \\text{ GeV}^2$";
      else if (ix == 1) title += "$\\langle Q^2\\rangle = 5-10   \\text{ GeV}^2$";
      else if (ix == 2) title += "$\\langle Q^2\\rangle = 10-20  \\text{ GeV}^2$";
      else if (ix == 3) title += "$\\langle Q^2\\rangle = 20-50  \\text{ GeV}^2$";
      else if (ix == 4) title += "$\\langle Q^2\\rangle = 50-100 \\text{ GeV}^2$";
      _histETLowQb.push_back(bookHistogram1D(ix+25,  1, 1, title));
      _weightETLowQb.push_back(0.);
    }
    // histograms and weight vectors for high
    _histETHighQb.reserve(3);
    _weightETHighQb.reserve(3);
    for (size_t ix = 0; ix < 3; ++ix) {
      string title = "Transverse energy flow for ";
      if (ix == 0)      title += "$\\langle Q^2\\rangle = 100-220 \\text{ GeV}^2";
      else if (ix == 1) title += "$\\langle Q^2\\rangle = 220-400 \\text{ GeV}^2";
      else if (ix == 2) title += "$\\langle Q^2\\rangle > 400 \\text{ GeV}^2";
      _histETHighQb.push_back(bookHistogram1D(ix+30,  1, 1, title));
      _weightETHighQb.push_back(0.);
    }
    // Histograms for the averages
    _histAverETCentral = bookProfile1D(33,  1, 1, 
				       "Average $E_\\perp$ in the central region");
    _histAverETFrag    = bookProfile1D(34,  1, 1, 
				       "Average $E_\\perp$ in the forward region");
  }

  // Finalize
  void H1_2000_S4129130::finalize() { 
    // normalization of the et distributions
    for(unsigned int ix=0;ix<17;++ix)
      scale(_histETLowQa[ix],1./_weightETLowQa[ix]);
    for(unsigned int ix=0;ix<7;++ix)
      scale(_histETHighQa[ix],1./_weightETHighQa[ix]);
    for(unsigned int ix=0;ix<5;++ix)
      scale(_histETLowQb[ix],1./_weightETLowQb[ix]);
    for(unsigned int ix=0;ix<3;++ix)
      scale(_histETHighQb[ix],1./_weightETHighQb[ix]);
  }

}
#include "Rivet/AnalysisLoader.hh"
