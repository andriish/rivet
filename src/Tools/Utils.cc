#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Utils.hh"


namespace Rivet {


  /// A function to get the Rivet version string
  string version() {
    return RIVET_VERSION;
  }


  ////////////////////////////////////////////////////////////////


  // Return distance of closest approach from track to given (primary) vertex position.
  double get2dClosestApproach(const HepMC::GenParticle& track, const Vector3& vtx3pos) {
    /// @todo Whoa! - implicit constructors from hell! 
    HepMC::FourVector trkvec = track;
    HepMC::ThreeVector trk3vec = trkvec;    
    HepMC::ThreeVector trk3pos = track.production_vertex()->position();
    
    Vector3 diff(vtx3pos.x()-trk3pos.x(), vtx3pos.y()-trk3pos.y(), vtx3pos.z()-trk3pos.z());
    
    // Impact parameter in the transverse plane
    const double d = fabs( trk3vec.x()*diff.y() - trk3vec.y()*diff.x() )
      / sqrt( sqr(trk3vec.x()) + sqr(trk3vec.y()) );
    return d;
  }


  // Return distance of closest approach from track to given (primary) vertex position.
  double get3dClosestApproach(const HepMC::GenParticle& track, const Vector3& vtx3pos) {    
    HepMC::FourVector trkvec = track;
    HepMC::ThreeVector trk3vec = trkvec;
    HepMC::FourVector trkpos = track.production_vertex()->position();
    HepMC::ThreeVector trk3pos = trkpos;
    Vector3 diff(vtx3pos.x()-trk3pos.x(), vtx3pos.y()-trk3pos.y(), vtx3pos.z()-trk3pos.z());
    
    // Impact parameter in 3 dimensions
    const double d = sqrt( sqr(trk3vec.y()*diff.z()-trk3vec.z()*diff.y()) - 
                           sqr(trk3vec.x()*diff.z()-trk3vec.z()*diff.x()) +
                           sqr(trk3vec.x()*diff.y()-trk3vec.y()*diff.x()) ) / trk3vec.mag();
    return d;
  }


  /// Return Decay Length Significance between two vertices in transverse plane
  double get2dDecayLength(const Vector3& vtx1, const Vector3& vtx2, const FourMomentum& jetaxis) {
    Vector3 diff = vtx1 - vtx2; 
    const double l = (jetaxis.px()*diff.x() + jetaxis.py()*diff.y() ) 
      / sqrt(sqr(jetaxis.px())+sqr(jetaxis.py()));
    return l;
  }



  /// Return 3 dimensional Decay Length Significance between vertices 
  double get3dDecayLength(const Vector3& vtx1, const Vector3& vtx2, const FourMomentum& jetaxis) {
    Vector3 diff = vtx1 - vtx2;
    const double l = (jetaxis.px()*diff.x() +jetaxis.py()*diff.y() +jetaxis.pz()*diff.z()) 
      / sqrt(sqr(jetaxis.px())+sqr(jetaxis.py())+sqr(jetaxis.pz()));
    return l;
  }

}
