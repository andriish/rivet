C-----------------------------------------------------------------------
      SUBROUTINE HWHEPC(IHEP)

#include "HERWIG65.INC"

      INTEGER I,J,IHEP,JMOTH
      CHARACTER *8 DUMMY

      IF (IHEP.EQ.NHEP) RETURN ! DON'T BOTHER
      DO I=IHEP+1,NHEP
         CALL HWUIDT(1,IDHEP(I),IDHW(I),DUMMY)
         JMOTH=JMOHEP(1,I)
         DO J=JDAHEP(1,JMOTH),JDAHEP(2,JMOTH)
            CALL HWVEQU(4,VHEP(1,JMOTH),VHEP(1,J)) ! VERTEX
         ENDDO
         IF (ABS(ISTHEP(JMOHEP(1,I))).EQ.2) THEN ! CORRECT STATUS CODES
            IF (ABS(IDHEP(JMOHEP(1,I))).EQ.15) THEN ! TAU
               ISTHEP(JMOHEP(1,I))=195
            ELSE
               ISTHEP(JMOHEP(1,I))=198 ! DECAYED RESONANCE
            ENDIF
         ENDIF
         IF (ABS(IDHEP(I)).EQ.111) THEN ! PI0 - TO BE DECAYED
            ISTHEP(I)=193 
         ENDIF
         IF (ABS(IDHEP(I)).EQ.221) THEN ! ETA - TO BE DECAYED
            ISTHEP(I)=193 
         ENDIF
      ENDDO

      RETURN
      END
