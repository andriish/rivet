// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Analyses/D0_2001_S4674421.hh" 

#include "AIDA/IDataPoint.h"


namespace Rivet {


  void D0_2001_S4674421::init() {	
    getLog() << Log::TRACE << "D0_2001_S4674421::init(): processing" << endl;
    _eventsFilledW = 0.0;
    _eventsFilledZ = 0.0;
    _h_dsigdpt_w = bookHistogram1D(1, 1, 1, "dsigma/dpT(W)");
    _h_dsigdpt_z = bookHistogram1D(1, 1, 2, "dsigma/dpT(Z)");

    const double binning[] = {0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 
			      25., 30., 35., 40., 50., 60., 70., 80., 100., 120., 160., 200.};
    vector<double> bins(23);
    for (int i=0; i<bins.size(); ++i)
      bins[i] = binning[i];
    
    _h_dsigdpt_scaled_z = bookHistogram1D("d01-x01-y03","dsigma/d(pT(Z)*mW/mZ)", bins);

    _h_temp = bookHistogram1D("temp", "dsigma/dpT(W) / dsigma/d(pT(Z)*mW/mZ)", bins);

  }


  void D0_2001_S4674421::isol(const ParticleVector& pvec, const VetoedFinalState& vfs, vector<double>& isolfrac, double rad) {

    isolfrac.resize(pvec.size());
    for (int i=0; i<isolfrac.size(); ++i)
      isolfrac[i] = 0.;

    int isolind = 0;
    for (ParticleVector::const_iterator p = pvec.begin(); 
	 p != pvec.end(); ++p, ++isolind) {  
      for (ParticleVector::const_iterator vfsp = vfs.particles().begin(); 
	   vfsp != vfs.particles().end(); ++vfsp) {
	if (vfsp != p  && deltaR(p->momentum(),vfsp->momentum()) < rad)
	  isolfrac[isolind] += vfsp->momentum().E();
      }
      //subtract and normalise to core/particle energy
      isolfrac[isolind] -= p->momentum().E();
      isolfrac[isolind] /= p->momentum().E();
    }
    
    return;
  }


  void D0_2001_S4674421::analyze(const Event & event) {
      const double weight = event.weight();
      
      //final state wothout neutrinos needed for isolation function
      const VetoedFinalState& vfs = applyProjection<VetoedFinalState>(event, "VFS");

      const LeadingParticlesFinalState& eeFS = applyProjection<LeadingParticlesFinalState>(event, "eeFS");
      
      //if there is a Z candidate
      if (!eeFS.isEmpty() && eeFS.particles().size()==2) {
	static size_t Zcount = 0;
	// Fill Z pT distributions
	const ParticleVector& Zees = eeFS.particles();

	//Check for isolation
	vector<double> isolfrac;
	isol(Zees, vfs, isolfrac, 0.4);
	//for (int i=0; i<isolfrac.size(); ++i) 
	  //cout << "Zees isolfrac=" << isolfrac[i] << endl; 

	ParticleVector::const_iterator p = Zees.begin();
	FourMomentum pmom = p->momentum();
	//cout << "Z 1 PDGid=" << p->getPdgId() << "   pT=" <<  pmom.pT() << endl;
	p++;
	pmom += p->momentum();
	//cout << "Z 2 PDGid=" << p->getPdgId() << "   pT=" << p->momentum().pT() << endl;
	//cout << "Zcand pT=" << pmom.pT() << "   m=" << sqrt(pmom.invariant()) << endl;
	double mass = sqrt(pmom.invariant());
	//if (isolfrac.size() == 2 && isolfrac[0]<0.15 && isolfrac[1]<0.15) {
	if (//p->momentum().pT() > 20. &&
	      mass > _mZmin && mass < _mZmax) {
	    Zcount += 1;
	    _eventsFilledZ += weight;
	    getLog() << Log::DEBUG << "Z #" << Zcount << " pmom.pT() = " << pmom.pT() << endl;
	    //cout << "Z #" << Zcount << " pmom.pT() = " << pmom.pT() << endl;
	    _h_dsigdpt_z->fill(pmom.pT(), weight);
	    _h_dsigdpt_scaled_z->fill(pmom.pT()*_mwmz, weight);
	    //}
	}
      }
      else { //no Z->ee candidate
	
	const LeadingParticlesFinalState& enuFS = applyProjection<LeadingParticlesFinalState>(event, "enuFS");
	
	const LeadingParticlesFinalState& enubFS = applyProjection<LeadingParticlesFinalState>(event, "enubFS");
	
	static size_t Wcount = 0;
	
	// Fill W pT distributions
	if (!enuFS.isEmpty() && enuFS.particles().size()==2 && enubFS.isEmpty()) {
	  const ParticleVector& Wenu = enuFS.particles();
	  ParticleVector::const_iterator p = Wenu.begin();
	  
	  ParticleVector Wel;
	  if (abs(p->getPdgId()) == 11) Wel.push_back(*p);

	  FourMomentum pmom = p->momentum();
	  //cout << "W 1 PDGid=" << p->getPdgId() << "   pT=" << pmom.pT() << endl;
	  p++;
	  pmom += p->momentum();
	  //cout << "W 2 PDGid=" << p->getPdgId() << "   pT=" << p->momentum().pT() << endl;
	  
	  if (abs(p->getPdgId()) == 11) Wel.push_back(*p);
	  //Check for isolation
	  vector<double> isolfrac;
	  isol(Wel, vfs, isolfrac, 0.4);
	  //cout << "Wel isolfrac.size()=" << isolfrac.size() << endl;
	  for (int i=0; i<isolfrac.size(); ++i) 
	    //cout << "Wel isolfrac=" << isolfrac[i] << endl; 
	
	  //if (isolfrac.size() == 1 && isolfrac[0] < 0.15) {
	  //if (p->momentum().pT() > 25.) {
	      Wcount++;
	  //cout << "W #" << Wcount << " pmom.pT() = " << pmom.pT() << endl;
	      _h_dsigdpt_w->fill(pmom.pT(), weight);
	      _eventsFilledW += weight;
	      //}
	      //}
	}
	else if (enuFS.isEmpty() && !enubFS.isEmpty() && enubFS.particles().size()==2) {
	  const ParticleVector& Wenub = enubFS.particles();
	  ParticleVector::const_iterator p = Wenub.begin();
	  
	  ParticleVector Wel;
	  if (abs(p->getPdgId()) == 11) Wel.push_back(*p);
	  
	  FourMomentum pmom = p->momentum();
	  p++;
	  pmom += p->momentum();

	  if (abs(p->getPdgId()) == 11) Wel.push_back(*p);
	  //Check for isolation
	  vector<double> isolfrac;
	  isol(Wel, vfs, isolfrac, 0.4);
	  //cout << "Wpos isolfrac.size()=" << isolfrac.size() << endl;
	  for (int i=0; i<isolfrac.size(); ++i) 
	    //cout << "Wpos isolfrac=" << isolfrac[i] << endl; 
	
	  //if (isolfrac.size() == 1 && isolfrac[0] < 0.15) {
	  //if (p->momentum().pT() > 25.) {
	      _h_dsigdpt_w->fill(pmom.pT(), weight);
	      _eventsFilledW += weight;
	      //}
	    //}
	}
	
      }    
      
      
  }


  void D0_2001_S4674421::finalize() { 
    // Get cross-section per event (i.e. per unit weight) from generator
    const double xSecPerEvent = crossSection()/picobarn / sumOfWeights();
    
    // Correct W pT distribution to W cross-section
    const double xSecW = xSecPerEvent * _eventsFilledW;
    
    // Correct Z pT distribution to Z cross-section
    const double xSecZ = xSecPerEvent * _eventsFilledZ;
    
    
    _h_temp = histogramFactory().divide(getName() + "/temp", *_h_dsigdpt_w, *_h_dsigdpt_scaled_z);
    
    
    const double wpt_integral = integral(_h_dsigdpt_w);
    normalize(_h_dsigdpt_w, xSecW);
    
     
    normalize(_h_dsigdpt_z, xSecZ);
    
    
    const double zpt_scaled_integral = integral(_h_dsigdpt_scaled_z);
    normalize(_h_dsigdpt_scaled_z, xSecZ);
    
    
    _h_temp->scale( (xSecW / wpt_integral) / (xSecZ / zpt_scaled_integral)  *  _brzee / _brwenu);
    
    std::vector<double> _x, _y, _ex, _ey;
    for ( int i = 0, N = _h_temp->axis().bins(); i < N; ++i ) {
      _x.push_back((_h_temp->axis().binLowerEdge(i) + _h_temp->axis().binUpperEdge(i))/2.0);
      _ex.push_back(_h_temp->axis().binWidth(i)/2.0);
      _y.push_back(_h_temp->binHeight(i)); 
      _ey.push_back(_h_temp->binError(i)); 
    }
    _dset_dsigpt_wz_rat = datapointsetFactory().createXY(getName() + "/d02-x01-y01", _h_temp->title(), _x, _y, _ex, _ey);
    
    
    //delete _h_temp;
    //tree.rm(getName() + "/temp");
    
    
  }

}

