// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/JetShape.hh"
#include "Rivet/Projections/Cmp.hh"
#include "HepPDT/ParticleID.hh"

using namespace inline_maths;

namespace Rivet {

  int JetShape::compare(const Projection& p) const {
    const JetShape& other = dynamic_cast<const JetShape&>(p);
    return pcmp(*_vfsproj, *other._vfsproj) 
      && (*_jetaxes == *other._jetaxes);
  }


  void JetShape::project(const Event& e) {
    Log& log = getLog();

    /*
    _momentum = LorentzVector();
    _momentum.setPx(0.0);
    _momentum.setPy(0.0);
    _momentum.setPz(0.0);
    _momentum.setE(0.0);

    */

    // Project into final state
    //const FinalState& fs = e.applyProjection(*_fsproj);

    // Project into vetoed final state if not done yet (supposed to be done by jet algorithm)
    const VetoedFinalState& vfs = e.applyProjection(*_vfsproj);
    ////const VetoedFinalState& vfs = *_vfsproj;

    //clear for each event and resize with zero vectors
    for (unsigned int i=0; i<_diffjetshapes->size(); ++i) {
      _diffjetshapes[i].clear();
    }
    _diffjetshapes->clear();

    for (unsigned int i=0; i<_jetaxes->size(); ++i) {
      vector<double> diffjetshape(_nbins,0.);
      _diffjetshapes->push_back(diffjetshape);
    }
   
    for (unsigned int i=0; i<_intjetshapes->size(); ++i) {
      _intjetshapes[i].clear();
    }
    _intjetshapes->clear();

    for (unsigned int i=0; i<_jetaxes->size(); ++i) {
      vector<double> intjetshape(_nbins,0.);
      _intjetshapes->push_back(intjetshape);
    }

    _oneminPsishape->clear();
    _oneminPsishape->resize(_jetaxes->size(),0.);
      
    
    //Determine jet shapes
    double y1, y2, eta1, eta2, phi1, phi2, drad, dradmin;
    int dradminind;
    for (ParticleVector::const_iterator p = vfs.particles().begin(); p != vfs.particles().end(); ++p) {
      //for (vector<LorentzVector>::iterator jt = _jetaxes->begin(); jt != _jetaxes->end(); ++jt) {
      for (unsigned int j=0; j<_jetaxes->size(); j++) {
	y1 = _jetaxes->at(j).rapidity();
	y2 = p->getMomentum().rapidity();
	eta1 = _jetaxes->at(j).pseudoRapidity();
	eta2 = p->getMomentum().eta();
	phi1 = _jetaxes->at(j).phi();
	phi2 = p->getMomentum().phi();
	
	if (_distscheme==snowmass) 	
	  drad = delta_rad(eta1, phi1, eta2, phi2);
	else //_distscheme = energy
	  drad = delta_rad(y1, phi1, y2, phi2);
	
	if (j==0 || drad<dradmin) {
	  dradminind=j;
	  dradmin=drad;
	}
      }

      for (int i=0; i<_nbins; ++i) {
	if (dradmin<_rmin+(i+1)*_interval) {
	  (*_intjetshapes)[dradminind][i] += p->getMomentum().perp();
	  if (dradmin>_rmin+i*_interval)
	    (*_diffjetshapes)[dradminind][i] += p->getMomentum().perp()/_interval;
	}
      }
      if (dradmin<_r1minPsi/_rmax)
	(*_oneminPsishape)[dradminind] += p->getMomentum().perp();

    }

    
    //normalize to total pT
    for (unsigned int j=0; j<_jetaxes->size(); j++) {
      if ((*_intjetshapes)[j][_nbins-1] > 0.) {
	(*_oneminPsishape)[j] = 1.-(*_oneminPsishape)[j]/(*_intjetshapes)[j][_nbins-1];
	for (int i=0; i<_nbins; ++i) {
	  (*_diffjetshapes)[j][i] /= (*_intjetshapes)[j][_nbins-1];
	  (*_intjetshapes)[j][i] /= (*_intjetshapes)[j][_nbins-1];
	}
      }
    }

    log << Log::DEBUG << "Done" << endl;
  }


}
