// -*- C++ -*-
//     Jet algorithm from the Field & Stuart minimum bias study.

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/TrakJet.hh"
#include "Rivet/Projections/Cmp.hh"

using namespace Rivet;

void KtJets::project(const Event & e) {

  vector<KtJet::KtLorentzVector> vecs;

  // Project into final state
  // NB to be true to the original, the final state projection should have 
  // specific cuts on eta, pT and require statble charged particles.
  const FinalState& fs = e.applyProjection(*fsproj);
  
  // Store 4 vector data about each particle into vecs
  for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
    HepMC::FourVector fv = p->getMomentum();
    // store the FourVector in the KtLorentzVector form
    KtJet::KtLorentzVector ktlv(fv.px(), fv.py(), fv.pz(), fv.e());
    vecs.push_back(ktlv);
  }

  // Now sort particles in pT

  // Now find jets using Field & Stuart criteria --------------

  //Ijet=0
  for 
              
    DO loop5=1,n_accepted_particles
               
               Deta=0.
               Dphi=0.
               
c need to include this line for not overcounting            
             
               IF(celija2(4,loop5).eq.0.) THEN 
                  
                  
                  Ijet=Ijet+1
                  
                  
                  DO loop4=loop5,n_accepted_particles
                     
                     IF (celija2(4,loop4).eq.0.) then 
                        
                        Deta=ABS(celija2(1,loop5)-celija2(1,loop4))
                        Dphi=ABS(celija2(3,loop5)-celija2(3,loop4))
                        
                        IF (Dphi.gt.Pi) then 
                           
                           Dphi=TWOPI-Dphi            
                           
                        END IF                
                        
                        IF (SQRT(Deta**2+Dphi**2).LE.(0.7)) then
                           celija2(4,loop4)=Ijet
                        END IF
                  
                     END IF

                  END DO
                  
               END IF  
              
              END DO

C --sum of PT for a given jet + numb particles in a jet -------            
C --- phi of the jets     ------------------------------------- 
C - also particles must satisfy -1 < eta < 1 and PT > 0.5 GeV -- 
C ---First Calculate SUMPT and then the rest.              
                
        DO loop6=1,Ijet
                
           SUMPT4=0.               
           SUMPHI4=0.
           MAXPHI=0.
           MINPHI=TWOPI

C     initialisation -------------------------
                 
           DO loop7=1,n_accepted_particles
              
              IF(INT(celija2(4,loop7)).eq.loop6) THEN
                 
                 SUMPT4=SUMPT4+celija2(2,loop7)                 
                 SUMPHI4=SUMPHI4+celija2(3,loop7)*celija2(2,loop7)   
                 
                 IF(celija2(3,loop7).gt.Maxphi) 
     &                Maxphi=celija2(3,loop7)                  
                 IF(celija2(3,loop7).lt.Minphi)
     &                Minphi=celija2(3,loop7)  
                                  
              END IF
           END DO

           IF((Maxphi-Minphi).gt.(3*Pi/2.)) THEN
              
              SUMPHI4=0.
              
              DO loop10=1,n_accepted_particles
                 IF (int(celija2(4,loop10)).eq.loop6) THEN 

                    IF (celija2(3,loop10).gt.Pi) THEN
                       TPHI=celija2(3,loop10)-TWOPI
                    ELSE
                       TPHI=celija2(3,loop10) 
                    END IF

                    SUMPHI4=SUMPHI4+TPHI*celija2(2,loop10)
                 END IF               
              END DO 
           END IF
                 
           SUMPT(loop6)=SUMPT4  
           JETPHI(loop6)=SUMPHI4/SUMPT4
           
           IF(JETPHI(loop6).lt.0.) JETPHI(loop6)=
     &          JETPHI(loop6)+Twopi 
        END DO
         
c taking the highest PT jet

        DO loop3=1,Ijet
           SUMPT2(loop3)=SUMPT(loop3) 
           JETPHI2(loop3)=JETPHI(loop3)              
        END DO
        
        MAXSUMPT=SUMPT2(1) 
        
        DO loop10=2,Ijet
           IF (SUMPT2(loop10).gt.MAXSUMPT) then
              MAXSUMPT=SUMPT2(loop10)
           END IF
        END DO  
        
c pointer to the highest ET jet
        DO loop11=1,Ijet
           IF(ABS(SUMPT(loop11)-MAXSUMPT)
     &          .lt.0.000001) THEN
              i_highest_et=loop11
           END IF
        END DO  
      END IF  
      
      IF (n_accepted_particles.eq.1) THEN  
         SUMPT(n_accepted_particles)=celija2(2,n_accepted_particles)       
         JETPHI(n_accepted_particles)=celija2(3,n_accepted_particles)
         i_highest_et=n_accepted_particles         
      END IF                     
      
      RETURN      
      END 

		      }
