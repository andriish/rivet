// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/H1_1994_S2919893.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"

namespace Rivet {

  void H1_1994_S2919893::analyze(const Event& event) {
    const FinalState & fs = applyProjection<FinalState>(event, "FS");
    const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");
    const DISLepton    & dl = applyProjection<DISLepton>(event,"Lepton");
    // require outgoing lepton same type as incoming
    if(dl.in().getPdgId()!=dl.out().getPdgId()) vetoEvent(event);
    // get the DIS kinematics
    double x  = dk.x();
    double w2 = dk.W2();
    double w = sqrt(w2);
    // whether electron in + (false) or - (true) z
    bool order = dl.in().getMomentum().z()<0.;
    // momentum of the scattered lepton
    FourMomentum leptonMom = dl.out().getMomentum();
    double ptel = pT(leptonMom);
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
      double th = beamAngle(p->getMomentum(),order);
      if(th>4.4&&th<15.) efwd += p->getMomentum().E();
    }
    // apply the cuts
    // lepton energy and angle, w2 and forward energy
    bool cut = enel>14. && thel>157. && thel<172.5 && 
               w2>=3000. && efwd>0.5;
    // veto event if doesn't pass cuts
    if(!cut) vetoEvent(event);
    // weight of the event
    const double weight = event.weight();
    // weights for x<1e-3 and x>1e-3
    if(x<1e-3) _wEnergy.first  += weight;
    else       _wEnergy.second += weight;
    // boost to hadronicCM
    const LorentzTransform hcmboost = dk.boostHCM();
    // loop over the particles
    long ncharged(0);
    for(ParticleVector::const_iterator p = particles.begin();
	p!= particles.end(); ++p ) {
      double th = beamAngle(p->getMomentum(),order);
      // boost momentum to lab
      const FourMomentum hcmMom = hcmboost.transform(p->getMomentum());
      // angular cut
      if(th<=4.4) continue;
      // energy flow histogram
      double et = fabs(Et(hcmMom));
      double eta = -hcmMom.pseudorapidity(); 
      if(x<1e-3) _histEnergyFlowLowX ->fill(eta,et*weight);
      else       _histEnergyFlowHighX->fill(eta,et*weight);
      if(PID::threeCharge(p->getPdgId())!=0) {
	if(w> 50.&&w<=200.) {
	  double xf= -2.*hcmMom.z()/w;
	  double pt2 = pT2(hcmMom);
	  if(w> 50.&&w<=100.)      _histSpectraW77 ->fill(xf,weight); 
	  else if(w>100.&&w<=150.) _histSpectraW122->fill(xf,weight);
	  else if(w>150.&&w<=200.) _histSpectraW169->fill(xf,weight);
	  _histSpectraW117->fill(xf,weight);
	  _histPT2->fill(xf,pt2*weight,weight);
	  ++ncharged;
	}
      }
      // energy-energy correlation
      if(th<=8.) continue;
      double phi1 = p->getMomentum().azimuthalAngle(ZERO_2PI);
      double eta1 = p->getMomentum().pseudorapidity();
      double et1 = fabs(Et(p->getMomentum()));
      for(ParticleVector::const_iterator p2 = p+1;
	  p2!=particles.end();++p2) {
	double th2 = beamAngle(p2->getMomentum(),order);
	if(th2<=8.) continue;
	double phi2 = p2->getMomentum().azimuthalAngle(ZERO_2PI);
	double deltaphi = phi1-phi2;
	if(fabs(deltaphi)>pi) deltaphi=fabs(fabs(deltaphi)-2.*pi);
	double eta2 = p2->getMomentum().pseudorapidity();
	double omega = sqrt(sqr(eta1-eta2)+sqr(deltaphi));
	double et2 = fabs(Et(p2->getMomentum()));
	double wt = et1*et2/sqr(ptel)*weight;
	if(x<1e-3) _histEECLowX ->fill(omega,wt);
	else       _histEECHighX->fill(omega,wt);
      }

    }
    // factors for normalization
    if(w>50.&&w<=200.) {
      if(w>50.&&w<=100.)  {
	_w77.first  += ncharged*weight;
	_w77.second += weight;
      }
      else if(w>100.&&w<=150.) {
	_w122.first  += ncharged*weight;
	_w122.second += weight;
      }
      else {
	_w169.first  += ncharged*weight;
	_w169.second += weight;
      }
      _w117.first  += ncharged*weight;
      _w117.second += weight;
    }
  }

  void H1_1994_S2919893::init() {
    _w77  = make_pair(0.,0.);
    _w122 = make_pair(0.,0.);
    _w169 = make_pair(0.,0.);
    _w117 = make_pair(0.,0.);
    _wEnergy = make_pair(0.,0.);
    _histEnergyFlowLowX =
      bookHistogram1D(1,  1, 1, 
		      "Transverse energy flow as a function of rapidity x<1e-3");
    _histEnergyFlowHighX = 
      bookHistogram1D(1,  1, 2, 
		      "Transverse energy flow as a function of rapidity x>1e-3");
    _histEECLowX         = bookHistogram1D(2, 1, 1, "EEC for x<1e-3");
    _histEECHighX        = bookHistogram1D(2, 1, 2, "EEC for x>1e-3");
    _histSpectraW77      = bookHistogram1D(3, 1, 1, 
					   "Charged particle spectra for  50<W<100");
    _histSpectraW122     = bookHistogram1D(3, 1, 2, 
					   "Charged particle spectra for 100<W<150");
    _histSpectraW169     = bookHistogram1D(3, 1, 3, 
					   "Charged particle spectra for 150<W<200");
    _histSpectraW117     = bookHistogram1D(3, 1, 4, 
					   "Charged particle spectra for all W    ");
    _histPT2             = bookProfile1D  (4, 1, 1, 
					   "Average pT2 as function of xF");
  }

  // Finalize
  void H1_1994_S2919893::finalize() { 
    // Normalize inclusive single particle distributions to the average number 
    // of charged particles per event.
    double avgNumParts = _w77.first/_w77.second;
    normalize(_histSpectraW77, avgNumParts);
    avgNumParts = _w122.first/_w122.second;
    normalize(_histSpectraW122, avgNumParts);
    avgNumParts = _w169.first/_w169.second;
    normalize(_histSpectraW169, avgNumParts);
    avgNumParts = _w117.first/_w117.second;
    normalize(_histSpectraW117, avgNumParts);
    scale(_histEnergyFlowLowX , 1./_wEnergy.first );
    scale(_histEnergyFlowHighX, 1./_wEnergy.second);
    scale(_histEECLowX , 1./_wEnergy.first );
    scale(_histEECHighX, 1./_wEnergy.second); 
  }

}
