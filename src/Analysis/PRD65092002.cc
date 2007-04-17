// -*- C++ -*-
/**
 *       Make plots for Field & Stuart paper on underlying events at CDF.
 *       
 *       Phys.Rev.D65:092002,2002
 *       FNAL-PUB 01/211-E
 *       
 */

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/PRD65092002.hh"
using namespace Rivet;

#include "AIDA/IHistogram1D.h"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;

	
        REAL PT1(50),PT2(50),PT3(49),PT1E(50),PT2E(50) 
	
        DATA PT1/0.277,0.55556,1.1667,1.722,1.9444,2.5,2.7778,
	1    3.2778,3.3889,3.611,3.8333,4.167,5.0,5.01,6.0556,5.2,
	2    4.45,5.833,6.778,7.23,7.667,7.7778,8.3889,8.4,8.61,
	3    9.16667,9.44,9.89,10.556,10.60,11.11,10.5,12.16667,
	4    13.33,13.50,12.4,13.50,13.8889,12.1667,14.44,14.8,17.167,
	5    15.0,15.5,14.78,16.05,16.667,17.778,17.80,21.667/  
	
        DATA PT2/0.118,0.529,1.024,1.47,1.8,2.0,2.235,2.19,2.238,
	1    2.34,2.47,2.706,2.588,2.46,2.706,2.235,1.588,3.12,2.69,
	2    2.812,2.705,2.69,2.705,2.70,2.702,2.69,2.529,2.705,2.94,
	3    3.0,2.94,2.824,2.577,2.647,2.88,2.90,3.235,2.929,2.589,
	4    2.94,2.92,3.118,4.18,3.118,2.647,3.0,2.824,3.118,3.88,
	5    4.47/    

        DATA PT1E/0.277,0.277,0.277,0.277,0.277,0.277,
	1    0.277,0.277,0.277,0.277,0.277,0.277,0.583,0.583,
	2    0.583,0.583,1.667,0.277,0.277,0.277,0.277,0.277,
	3    0.277,0.277,0.277,0.277,0.69,0.277,0.69,0.277,
	4    0.277,0.69,0.277,0.69,0.69,1.389,1.25,0.69,0.833,
	5    0.69,1.11,0.69,2.22,1.1,1.67,
	6    1.388,2.125,1.67,2.22,2.77/      
	

	DATA PT2E/0.0589,0.0589,0.0589,0.0589,0.0589,0.12,0.12,0.15,
	1    0.223,0.235,0.235,0.235,0.465,0.288,0.47,0.412,
	2    0.406,1.82,0.176,0.118,0.147,0.176,0.118,0.294,
	3    0.288,0.235,0.176,0.235,0.118,0.176,0.176,0.412,
	4    0.122,0.347,0.529,0.235,0.235,0.582,0.235,0.264,
	5    0.294,0.588,0.582,1.0,0.441,0.764,0.705,0.88,1.05,
	6    1.41/

	DATA PT3/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.5,10.0,11.0,12.0,
	1    13.0,14.0,15.0,15.5,16.0,17.0,19.0,20.0,21.0,22.0,
	2    23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.5,31.0,32.0,
	3    33.0,34.0,35.0,36.0,37.0,38.0,39.0,40.5,41.0,42.0,
	4    43.0,44.0,45.0,46.0,47.0,48.0,49.0,49.8/    


	   

           CALL HBOOK1(30,'Experimental data away <PTsum>',
	1	50,0.0,50.0,0.)
           CALL HBOOK1(31,'Experimental data transverse <PTsum>',
	1	50,0.0,50.0,0.)
           CALL HBOOK1(32,'Experimental data toward <PTsum>',
	1	50,0.0,50.0,0.)   
           CALL HBPROF(2,'PTjet#1 toward-',50,0.0,50.0,
	1	0.0,20.0,'')
c       CALL HBOOK2(1,'NCHG vs PTjet',50,0.0,50.0,30,0.0
c       &     ,30.0,0.)
           CALL HBPROF(22,'PTjet#1 toward soft',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(23,'PTjet#1 toward hard',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(24,'PTjet#1 away soft',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(25,'PTjet#1 away hard',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(26,'PTjet#1 transverse hard',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(3,'PTjet#1 toward+',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(4,'PTjet#1 away -',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(5,'PTjet#1 away +',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(6,'PTjet#1 transverse tot',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(7,'PTjet#1 transverse soft ',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(10,'<PTSUM> toward-',50,0.0,50.0,
	1	0.0,100.,'') 
           CALL HBPROF(11,'<PTSUM> toward tot',50,0.0,50.0,
	1	0.0,100.,'')       
           CALL HBPROF(12,'<PTSUM> away-',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(13,'<PTSUM> away tot',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(14,'<PTSUM> transverse tot',50,0.0,50.0,
	1	0.0,100.,'')   
           CALL HBPROF(15,'<PTSUM> transverse soft',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(16,'<PTSUM> toward soft',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(17,'<PTSUM> toward hard',50,0.0,50.0,
	1	0.0,100.,'')           
           CALL HBPROF(18,'<PTSUM> away soft ',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(19,'<PTSUM> away hard ',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(20,'<PTSUM> transverse hard',50,0.0,50.0,
	1	0.0,100.,'')
           CALL HBPROF(27,'total numb of partl',50,0.0,50.0,
	1	0.0,100.,'')







	   CALL HIDOPT(0,'stat')
	   CALL HBARX(0)



           CALL HCDIR('//HISTO/f01211e',' ')
	   CALL HCDIR('//PAWC/f01211e',' ')

	   

	   DO loop600=1,2000
              

	      SUMPT(loop600)=0.
	      JETPHI(loop600)=0.
	      
	   END DO

	   
	   DO loop601=1,5

	      DO loop602=1,4000
		 
		 celija2(loop601,loop602)=0.
		 
	      END DO

	   END DO
	   
	   LOC=0                       

	   
           brojdva=0
	   


	   
        ELSE IF (intflag.eq.2) THEN


C       *** Filling runs ***
	   
C       *** Change directory to our graphs ***
	   
	   
           CALL HCDIR('//PAWC/f01211e',' ')
	   
c       initialise for each successive run

c       newjetalgo2 is in the case you use soft. 
c       newjetAlgo2mion is in the case you use MIon.


	   
           CALL jetalgo2       

**********************************************************
C       ------------- ploting new Algo-data -------------
**********************************************************
           
c       each cell must satisfay next criteria:
c       PT of the cell > 0.5 and ABS(eta) < 1 this go's as well for the
c       cells
c       used to calculate the PT of the jet.

C       Toward region without particles in the jet#1 :
C       ABS(PHIpart-PHIjet(leading)) < (60deg=1.047198rad)
           
C       count over the number of charged particles 
C       entering following regions with respect to the leading jet:
c       

c*******************************************************
C       only count hits where SUMPT is different from 0 !!!!!           
c*******************************************************
	   



C       - Calculate the PTsum which is the average scalar Pt sum for
C       - each of this regions

	   


	   Pi=3.14159265    

	   IF (brojdva.ne.0) THEN

	      DO loop10=1,11
		 

		 PTSUM(loop10)=0.
		 
	      END DO
	      

	      counttot=0  
	      count1=0
	      count2=0
	      count3=0
	      count4=0
	      count5=0
	      count6=0
	      count7=0
	      count8=0
	      count9=0
	      count10=0
	      count11=0   
	      
c       print *,'broj dva',brojdva


	      
	      IF (((SUMPT(LOC).gt.0.5).and.(SUMPT(LOC).lt.50.0))
	1	   .and.(LOC.ne.0)) THEN       
		 
		 

		 
		 
		 DO loop1=1,brojdva   


		    counttot=counttot+1         

		    IF((celija2(4,loop1).ne.celija2(4,LOC)).and.
	1		 (celija2(2,loop1).gt.0.5).and.
	2		 (ABS(celija2(1,loop1)).lt.1.0).and. 
	3		 (ABS(celija2(3,loop1)-JETPHI(LOC))
	4		 .lt.1.047197551)) THEN      

		       
c       count the number of particles entering this region
		       
		       
		       count1 = count1 + 1
		       PTSUM(1)=PTSUM(1)+celija2(2,loop1)             

		    END IF
		    
		    
		    
		    IF((celija2(5,loop1).eq.1.0).and.
	1		 (celija2(2,loop1).gt.0.5).and.
	2		 (ABS(celija2(1,loop1)).lt.1.0).and.
	3		 (ABS(celija2(3,loop1)-JETPHI(LOC))
	4		 .lt.1.047197551)) THEN      

		       
c       count the number of particles entering this region
		       
		       
		       count7 = count7 + 1
		       PTSUM(7)=PTSUM(7)+celija2(2,loop1)             

		    END IF
		    

		    
		    IF((celija2(5,loop1).eq.0.0).and.
	1		 (celija2(2,loop1).gt.0.5).and.
	2		 (ABS(celija2(1,loop1)).lt.1.0).and.
	3		 (ABS(celija2(3,loop1)-JETPHI(LOC))
	4		 .lt.1.047197551)) THEN      

		       
c       count the number of particles entering this region
		       
		       
		       count8 = count8 + 1
		       PTSUM(8)=PTSUM(8)+celija2(2,loop1)             

		    END IF
		    


C       Toward region defined with particles in jet#1
C       ABS(PHIpart-PHIjet(leading)) < (60deg=1.047197551rad)
		    
		    

		    IF((celija2(2,loop1).gt.0.5).and.
	1		 (ABS(celija2(1,loop1)).lt.1.0).and.
	2		 (ABS(celija2(3,loop1)-JETPHI(LOC))
	3		 .lt.1.047197551)) THEN
		       
		       PTSUM(2)=PTSUM(2)+celija2(2,loop1) 
		       count2 = count2 + 1              

		    END IF   
		    









c       Away region defined without particles in the jet#1
c       ABS(PHIpart-PHIjet) > (120deg=2.094395)


		    
		    
		    IF((celija2(2,loop1).gt.0.5).and.
	1		 (ABS(celija2(1,loop1)).lt.1.0).and.
	2		 (ABS(celija2(3,loop1)-JETPHI(LOC)).gt.2.0943951
	2		 ).and.(ABS(celija2(3,loop1)-JETPHI(LOC)).lt.4
	3		 .1887902)) THEN


		       count3 = count3 + 1
		       PTSUM(3)=PTSUM(3)+celija2(2,loop1) 
		       
		    END IF 

		    
		    IF   ((celija2(2,loop1).gt.0.5)
	1		 .and.(celija2(5,loop1).eq.1.0).and.
	2		 (ABS(celija2(1,loop1)).lt.1.0).and.
	3		 (ABS(celija2(3,loop1)-JETPHI(LOC)).gt.2.0943951
	3		 ).and.(ABS(celija2(3,loop1)-JETPHI(LOC)).lt.4
	4		 .1887902)) THEN


		       count9 = count9 + 1
		       PTSUM(9)=PTSUM(9)+celija2(2,loop1) 
		       
		    END IF 




		    IF   ((celija2(2,loop1).gt.0.5)
	1		 .and.(celija2(5,loop1).eq.0.0).and.
	2		 (ABS(celija2(1,loop1)).lt.1.0).and.
	3		 (ABS(celija2(3,loop1)-JETPHI(LOC)).gt.2.0943951
	3		 ).and.(ABS(celija2(3,loop1)-JETPHI(LOC)).lt.4
	4		 .1887902)) THEN


		       count10 = count10 + 1
		       PTSUM(10)=PTSUM(10)+celija2(2,loop1) 
		       
		    END IF 


c       Away region              
		    
		    IF((celija2(2,loop1).gt.0.5).and.
	1		 (ABS(celija2(1,loop1)).lt.1.0).and.
	2		 (ABS(celija2(3,loop1)-JETPHI(LOC)).gt.2.094395)
	2		 .and.(ABS(celija2(3,loop1)-JETPHI(LOC)).lt.4
	3		 .1887902)) THEN
		       
		       
		       PTSUM(4)=PTSUM(4)+celija2(2,loop1)
		       count4 = count4 + 1
		       
		    END IF 
		    
		    

c       the transverse region              
		    



		    IF((celija2(2,loop1).gt.0.5)
	1		 .and.
	2		 (ABS(celija2(1,loop1)).lt.1.0)
	3		 .AND.    
	4		 (((ABS(celija2(3,loop1)-JETPHI(LOC)).gt.1
	4		 .047197551).and.(ABS(celija2(3,loop1)
	5		 -JETPHI(LOC)).lt.2.094395102)).OR
	5		 .((ABS(celija2(3,loop1)-JETPHI(LOC)).gt.4
	6		 .1887902).and.(ABS(celija2(3,loop1)-JETPHI(LOC)
	6		 ).lt.5.235987756))))  THEN
		       

		       
		       
c       print *,'cell5=all',celija2(5,loop1)
		       count5 = count5 + 1 
		       PTSUM(5)=PTSUM(5)+celija2(2,loop1)

		    END IF   


c       the transverse region   


		    IF((celija2(2,loop1).gt.0.5)
	1		 .and.(celija2(5,loop1).eq.1.0)
	2		 .and.
	3		 (ABS(celija2(1,loop1)).lt.1.0)
	4		 .AND.
	5		 (((ABS(celija2(3,loop1)-JETPHI(LOC)).gt.1
	5		 .047197551).and.(ABS(celija2(3,loop1)
	6		 -JETPHI(LOC)).lt.2.094395102)).OR
	6		 .((ABS(celija2(3,loop1)-JETPHI(LOC)).gt.4
	7		 .1887902).and.(ABS(celija2(3,loop1)-JETPHI(LOC)
	7		 ).lt.5.235987756)))) THEN
		       
		       

		       PTSUM(6)=PTSUM(6)+celija2(2,loop1)
		       count6 = count6 + 1 


		    END IF


		    IF((celija2(2,loop1).gt.0.5)
	1		 .and.(celija2(5,loop1).eq.0.)
	2		 .and.
	3		 (ABS(celija2(1,loop1)).lt.1.0)
	4		 .AND.
	5		 (((ABS(celija2(3,loop1)-JETPHI(LOC)).gt.1
	5		 .047197551).and.(ABS(celija2(3,loop1)
	6		 -JETPHI(LOC)).lt.2.094395102)).OR
	6		 .((ABS(celija2(3,loop1)-JETPHI(LOC)).gt.4
	7		 .1887902).and.(ABS(celija2(3,loop1)-JETPHI(LOC)
	7		 ).lt.5.235987756) ))) THEN
		       
		       

		       PTSUM(11)=PTSUM(11)+celija2(2,loop1)
		       count11 = count11 + 1 


		    END IF

		 END DO

	      END IF

c       this will plot particles entering the different regions
c       from above.           

	      IF (((SUMPT(LOC).gt.0.5).and.(SUMPT(LOC).lt.50.0))
	1	   .and.(LOC.ne.0).and.brojdva.ne.0) THEN       
		 


		 CALL HFILL (27,real(SUMPT(LOC)),real(counttot),wtx) 

		 CALL HFILL (22,real(SUMPT(LOC)),real(count7),wtx)
		 CALL HFILL (16,real(SUMPT(LOC)),real(PTSUM(7)),wtx)  
		 
		 CALL HFILL (23,real(SUMPT(LOC)),real(count8),wtx)
		 CALL HFILL (17,real(SUMPT(LOC)),real(PTSUM(8)),wtx)  


		 CALL HFILL (24,real(SUMPT(LOC)),real(count9),wtx)
		 CALL HFILL (18,real(SUMPT(LOC)),real(PTSUM(9)),wtx)
	1	      
		 
		 CALL HFILL (25,real(SUMPT(LOC)),real(count10),wtx)
		 CALL HFILL (19,real(SUMPT(LOC)),real(PTSUM(10)),wtx)
	1	      
		 
		 CALL HFILL(26,real(SUMPT(LOC)),real(count11),wtx)
		 CALL HFILL (20,real(SUMPT(LOC)),real(PTSUM(11)),wtx)
	1	      
		 
		 CALL HFILL (3,real(SUMPT(LOC)),real(count2),wtx)
		 CALL HFILL (11,real(SUMPT(LOC)),real(PTSUM(2)),wtx)

		 CALL HFILL (4,real(SUMPT(LOC)),real(count3),wtx)
		 CALL HFILL (12,real(SUMPT(LOC)),real(PTSUM(3)),wtx) 

		 CALL HFILL (5,real(SUMPT(LOC)),real(count4),wtx)
		 CALL HFILL (13,real(SUMPT(LOC)),real(PTSUM(4)),wtx)
	1	      

                 CALL HFILL (6,real(SUMPT(LOC)),real(count5),wtx)
                 CALL HFILL (14,real(SUMPT(LOC)),real(PTSUM(5)),wtx)

                 CALL HFILL (7,real(SUMPT(LOC)),real(count6),wtx)
                 CALL HFILL (15,real(SUMPT(LOC)),real(PTSUM(6)),wtx)

	      ENDIF

	      
	   END IF

	   
	   


C       ---This will plot leading jet PT vs Number of particles
C       in the leading jet*************************************


	   

	   
	   
	ELSE IF (intflag.eq.3) THEN
	   
C       *** Termination run ***
C       *** Change directory ***
	   CALL HCDIR('//PAWC/f01211e',' ')
C       *** Finish off histograms ***

	   CALL  hpak(30,PT1)
	   CALL  hpake(30,PT1E)
	   CALL  hpak(31,PT2)
	   CALL  hpake(31,PT2E)
	   CALL  hpak(32,PT3)  



c       ****************************************************************
c       end of histo *********************        
******************************************************************
******************************************************************


           
	   
c       IF (Xsec.eq.0) THEN
c       PRINT*,'f01211E: termination called with zero cross section'
c       PRINT*,'          cross section graph meaningless'
c       Xsec=1
c       ENDIF
c       IF (Ntot.eq.0) THEN
c       PRINT*,'f01211E: termination called with no total events'
c       PRINT*,'          cross section graph meaningless'
c       Ntot=1
c       ENDIF
	   
C       --- store Xsec and Ntot for both direct and resolved events.
	   
c       mhXsec(iproc)=Xsec
c       mhntot(iproc)=Ntot

	ELSE

c       *** End ****

	   print*,'f01211e: Please run routine with iflag set to 1,2 or 3'
	   
	ENDIF

	RETURN

	END 
     

