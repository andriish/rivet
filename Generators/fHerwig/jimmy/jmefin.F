*-- Author : JMB
      SUBROUTINE JMEFIN
C -----------------------------------------------------------------
C     Adjust HERWIG cross section in light of Multiple Scattering.
C     Must be called after (from?) HWEFIN.
C -----------------------------------------------------------------
#include "HERWIG65.INC"
#include "Herwig_i/jimmy.inc"

      WRITE (JMOUT,1)
    1 FORMAT(/10X,'MODIFIED OUTPUT ON ELEMENTARY PROCESS'/)
      IF (NLOST.EQ.0) THEN
         WRITE (6,10)
 10      FORMAT(10X,
     &        'NO MULTIPLE SCATTERS:NO MODIFICATION NECESSARY')
      ELSE
         AVWGT=AVWGT*(TOTSCAT-NLOST)/TOTSCAT
         WRITE (6,11) INT(TOTSCAT),INT(NLOST),1000.*AVWGT
         
 11      FORMAT(1P,
     &        10X,'NUMBER OF SCATTERS = ',I20/
     &        10X,'NO.OF MULT.SCATTERS= ',I20/
     &        10X,'FINAL C-S (PB)     =',G12.4)
         
      ENDIF
      END
