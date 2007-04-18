// -*- C++ -*-
//     Jet algorithm from the Field & Stuart minimum bias study.

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/TrackJet.hh"
#include "Rivet/Projections/Cmp.hh"
#include "Rivet/RivetCLHEP.hh"
#include <algorithm>

using namespace Rivet;
using namespace CLHEP;
using namespace std;


int TrackJet::compare(const Projection& p) const {
  const TrackJet& other = dynamic_cast<const TrackJet&>(p);
  return pcmp(*_fsproj, *other._fsproj);
}


void TrackJet::project(const Event & e) {
  Log& log = getLog();
  _jets.clear();

  // Project into final state
  // NB to be true to the original, the final state projection should have 
  // specific cuts on eta, pT and require stable charged particles.
  log << Log::DEBUG << "About to apply the final state projection with eta and pt cuts" << endl;
  const FinalState& fs = e.applyProjection(*_fsproj);
  
  // Put each particle into the collection of tracks and sort (decreasing in pT)
  vector<LorentzVector> tracksvector; // Need to use a vector because you can't use std::sort with lists
  for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
    HepMC::FourVector fv = p->getMomentum();
    LorentzVector lv(fv.px(), fv.py(), fv.pz(), fv.e());
    tracksvector.push_back(lv);
  }
  // Now sort particles in pT
  sort(tracksvector.begin(), tracksvector.end(), compareVecsByPt);
  // Make it into a list
  Tracks tracks(tracksvector.begin(), tracksvector.end());

  // Now find jets using Field & Stuart criteria
  log << Log::DEBUG << "About to assign tracks into jets" << endl;
  while (!tracks.empty()) {

    Tracks::iterator t = tracks.begin();
    // Get eta and phi for this track
    const double eta = t->eta();
    const double phi = t->phi();
              
    // Make a new jet and put this seed track into it
    Jet thisjet;
    thisjet.addParticle(*t);
    tracks.erase(t);
    
    // Compare with all unassociated tracks with a smaller pT measure
    Tracks::iterator t2 = tracks.begin();
    while (t2 != tracks.end()) {
      // Compute Deta and Dphi, mapping Dphi into [0,pi]
      double Deta = eta - t2->eta();
      double Dphi = phi - t2->phi();
      if (Dphi > PI) Dphi = 2*PI - Dphi;

      // Add to this jet if eta-phi distance < 0.7 
      if (sqrt(Deta*Deta + Dphi*Dphi) <= 0.7) {
        // Move this particle into the current jet (no extra sorting needed)
        thisjet.addParticle(*t2);
        tracks.erase(t2++); // the postfix ++ is important!
      } else {
        ++t2;
      }
    }

    _jets.push_back(thisjet);
  }

  // Sort the jets by pT.
  sort(_jets.begin(), _jets.end(), compareJetsByPt);


  log << Log::INFO << "Number of jets = " << _jets.size() << endl;
  size_t njet = 1;
  for (Jets::const_iterator j = _jets.begin(); j != _jets.end(); ++j)
    log << Log::INFO << "Number of tracks in jet #" << njet++ << " = " << j->getNumParticles() << endl;



  // Need to provide a pT-weighted phi, a total pT and a total number of particles in the jet


// C --sum of PT for a given jet + numb particles in a jet -------            
// C --- phi of the jets     ------------------------------------- 
// C - also particles must satisfy -1 < eta < 1 and PT > 0.5 GeV -- 
// C ---First Calculate SUMPT and then the rest.              
                
//         DO loop6=1,Ijet
                
//            SUMPT4=0.               
//            SUMPHI4=0.
//            MAXPHI=0.
//            MINPHI=TWOPI

// C     initialisation -------------------------

/////////// ******* Get the max and min phi (per jet or over the whole event?)
                 
//            DO loop7=1,n_accepted_particles
              
//               IF(INT(celija2(4,loop7)).eq.loop6) THEN
                 
//                  SUMPT4=SUMPT4+celija2(2,loop7)                 
//                  SUMPHI4=SUMPHI4+celija2(3,loop7)*celija2(2,loop7)   
                 
//                  IF(celija2(3,loop7).gt.Maxphi) 
//      &                Maxphi=celija2(3,loop7)                  
//                  IF(celija2(3,loop7).lt.Minphi)
//      &                Minphi=celija2(3,loop7)  
                                  
//               END IF
//            END DO

//            IF((Maxphi-Minphi).gt.(3*Pi/2.)) THEN
              
//               SUMPHI4=0.
              
//               DO loop10=1,n_accepted_particles
//                  IF (int(celija2(4,loop10)).eq.loop6) THEN 

//                     IF (celija2(3,loop10).gt.Pi) THEN
//                        TPHI=celija2(3,loop10)-TWOPI
//                     ELSE
//                        TPHI=celija2(3,loop10) 
//                     END IF

//                     SUMPHI4=SUMPHI4+TPHI*celija2(2,loop10)
//                  END IF               
//               END DO 
//            END IF
                 
//            SUMPT(loop6)=SUMPT4  
//            JETPHI(loop6)=SUMPHI4/SUMPT4
           
//            IF(JETPHI(loop6).lt.0.) JETPHI(loop6)=
//      &          JETPHI(loop6)+Twopi 
//         END DO
         
// c taking the highest PT jet

//         DO loop3=1,Ijet
//            SUMPT2(loop3)=SUMPT(loop3) 
//            JETPHI2(loop3)=JETPHI(loop3)              
//         END DO
        
//         MAXSUMPT=SUMPT2(1) 
        
//         DO loop10=2,Ijet
//            IF (SUMPT2(loop10).gt.MAXSUMPT) then
//               MAXSUMPT=SUMPT2(loop10)
//            END IF
//         END DO  
        
// c pointer to the highest ET jet
//         DO loop11=1,Ijet
//            IF(ABS(SUMPT(loop11)-MAXSUMPT)
//      &          .lt.0.000001) THEN
//               i_highest_et=loop11
//            END IF
//         END DO  
//       END IF  
      
//       IF (n_accepted_particles.eq.1) THEN  
//          SUMPT(n_accepted_particles)=celija2(2,n_accepted_particles)       
//          JETPHI(n_accepted_particles)=celija2(3,n_accepted_particles)
//          i_highest_et=n_accepted_particles         
//       END IF                     
      
//       RETURN      
//       END 

		      }
