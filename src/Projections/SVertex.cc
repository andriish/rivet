// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Projections/SVertex.hh"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenEvent.h"
#include "HepMC/SimpleVector.h"
#include "Rivet/Cmp.hh"
#include "Rivet/RivetCLHEP.hh"

#include "Rivet/Tools/Utils.hh"
using namespace math;

namespace Rivet {

  int SVertex::compare(const Projection& p) const {
    const SVertex& other = dynamic_cast<const SVertex&>(p);
    const int fscmp = pcmp(_pvtx, other._pvtx);
    if (fscmp != 0) return fscmp;
    return (_applyVtxTrackCuts == other._applyVtxTrackCuts &&
	    _jetaxes == other._jetaxes && _detEta == other._detEta &&
	    _IPres == other._IPres && _DLS == other._DLS && 
	    _DLSres == _DLSres);

  }




  void SVertex::project(const Event& e) {
 
   Log log = getLog();

    const PVertex& pvtx = e.applyProjection(_pvtx);
    const ChargedFinalState& chfs = e.applyProjection(_chfs);
    
    const HepMC::GenVertex& gpvtx = pvtx.getPrimaryVertex();
    HepMC::FourVector pvpos = gpvtx.position();



    ///produce vector of vertices, each containing a vector of all charged 
    ///final state particles belonging to this vertex
    ParticleVector gps;
    vector<ParticleVector> chfsvtx;
    for (ParticleVector::const_iterator p = chfs.particles().begin(); p != chfs.particles().end(); ++p) {

      ///condiser only charged particles in tracker geometrical acceptance
      if (fabs(p->getMomentum().pseudoRapidity())> _detEta) continue;

      unsigned int i;
      for (i=0; i<chfsvtx.size(); ++i) {
	if (chfsvtx[i].size()>0 && 
	    chfsvtx[i].begin()->getHepMCParticle().production_vertex() 
	    == p->getHepMCParticle().production_vertex()) {
	  chfsvtx[i].push_back(*p);
	  break;
	}
      }
      if (i==chfsvtx.size()) { //particle of new vertex
	chfsvtx.push_back(gps);
	chfsvtx[chfsvtx.size()-1].push_back(*p);
      }
          
    } //end for loop over Particle Vector's

  
    ///Check if jets are tagged, by means of selected vertices 
    ///fulfilling track criteria
    _taggedjets.clear();
    for (unsigned int i=0; i<chfsvtx.size(); ++i) {
      LorentzVector vtxVisMom;
      if (_applyVtxTrackCuts(*this, chfsvtx[i], gpvtx, vtxVisMom) ) {

	for (unsigned int j=0; j<_jetaxes.size(); ++j) {
	  ///deltaR requirement between jet and visible vector sum of vertex tracks
	  if (_jetaxes[j].deltaR(vtxVisMom) < _deltaR) {
	    double dls = 
	      get2dDLS(*chfsvtx[i].begin()->getHepMCParticle().production_vertex(),
		       gpvtx, _jetaxes[j]);
	    if (dls > _DLS) _taggedjets.push_back(_jetaxes[j]);
	  }
	}
      }
    }

  }
  



  /// Return Distance of Closest Approach from track to given (primary) vertex
  double SVertex::get2dDCA(const HepMC::GenParticle& track, const HepMC::GenVertex& gvtx) {
    HepMC::ThreeVector pv3pos = gvtx.position();

    HepMC::FourVector trkvec = track;
    HepMC::ThreeVector trk3vec = trkvec;

    HepMC::ThreeVector trk3pos = track.production_vertex()->position();

    HepMC::ThreeVector diff;
    diff.set(pv3pos.x()-trk3pos.x(), pv3pos.y()-trk3pos.y(), pv3pos.z()-trk3pos.z());
    
    /// Impact Parameter in the transverse plane
    double d = fabs( trk3vec.x()*diff.y()-trk3vec.y()*diff.x() )
      / sqrt( sqr(trk3vec.x()) + sqr(trk3vec.y()) );

    return d;
    
  }



  /// Return Distance of Closest Approach from track to given (primary) vertex
  double SVertex::get3dDCA(const HepMC::GenParticle& track, const HepMC::GenVertex& gvtx) {
    HepMC::FourVector pvpos = gvtx.position();
    HepMC::ThreeVector pv3pos = pvpos;

    HepMC::FourVector trkvec = track;
    HepMC::ThreeVector trk3vec = trkvec;

    HepMC::FourVector trkpos = track.production_vertex()->position();
    HepMC::ThreeVector trk3pos = trkpos;

    HepMC::ThreeVector diff;
    diff.set(pv3pos.x()-trk3pos.x(), pv3pos.y()-trk3pos.y(), pv3pos.z()-trk3pos.z());

    /// Impact Parameter in 3 dimesnions
    double d = sqrt( sqr(trk3vec.y()*diff.z()-trk3vec.z()*diff.y())
    	    -sqr(trk3vec.x()*diff.z()-trk3vec.z()*diff.x())
    	    +sqr(trk3vec.x()*diff.y()-trk3vec.y()*diff.x()) ) 
      / trk3vec.mag();
    
    return d;
    
  }


  
  ///Return Impact Parameter Significance of givern track w.r.t. (primary) vertex
  double SVertex::get2dDCAsig(const HepMC::GenParticle& track, const HepMC::GenVertex& gvtx) {
    double d = get2dDCA(track, gvtx);

    return d/_IPres;
  }



  ///Return Impact Parameter Significance of givern track w.r.t. (primary) vertex
  double SVertex::get3dDCAsig(const HepMC::GenParticle& track, const HepMC::GenVertex& gvtx) {
    double d = get3dDCA(track, gvtx);

    return d/_IPres;
  }
  


  /// Return Decay Length Significance between two vertices in transverse plane
  double SVertex::get2dDLS(const HepMC::GenVertex& vtx1, const HepMC::GenVertex& vtx2,
			   LorentzVector& jetaxis) {
    HepMC::ThreeVector vtxpos1 = vtx1.position();
    HepMC::ThreeVector vtxpos2 = vtx2.position();
    
    HepMC::ThreeVector diff;
    diff.set(vtxpos1.x()-vtxpos2.x(), vtxpos1.y()-vtxpos2.y(), vtxpos1.z()-vtxpos2.z());

    //double l = fabs(sqr(diff.x())+sqr(diff.y()));
    double l = (jetaxis.px()*diff.x() + jetaxis.py()*diff.y() ) 
      / sqrt(sqr(jetaxis.px())+sqr(jetaxis.py()));
    
    return l/_DLSres;

  }



  /// Return 3 dimensional Decay Length Significance between vertices 
  double SVertex::get3dDLS(const HepMC::GenVertex& vtx1, const HepMC::GenVertex& vtx2,
		  LorentzVector& jetaxis) {
    HepMC::ThreeVector vtxpos1 = vtx1.position();
    HepMC::ThreeVector vtxpos2 = vtx2.position();
    
    HepMC::ThreeVector diff;
    diff.set(vtxpos1.x()-vtxpos2.x(), vtxpos1.y()-vtxpos2.y(), vtxpos1.z()-vtxpos2.z());

    double l = (jetaxis.px()*diff.x() +jetaxis.py()*diff.y() +jetaxis.pz()*diff.z()) 
      / sqrt(sqr(jetaxis.px())+sqr(jetaxis.py())+sqr(jetaxis.pz()));
    
    return l/_DLSres;





    
    return 0.;
  }
  
  



} // end namespace Rivet
