CDECK  ID>, HWDSIN.
*CMZ :-        -30/09/02  14:05:28  by  Peter Richardson
*-- Author :    Peter Richardson
C-----------------------------------------------------------------------
      SUBROUTINE HWDSIN(CLSAVE)
C-----------------------------------------------------------------------
C  Subroutine to perform decays including spin correlations
C-----------------------------------------------------------------------
#include "HERWIG65.INC"
      DOUBLE PRECISION PW(5)
      INTEGER IDEC,IP,IS,IHEP,ID,IM,LHEP,MHEP,NPR,KHEP,CLSAVE(2),NTRY,
     &     ID1
      IF(IERROR.NE.0) RETURN
      NTRY = 0
      IDEC = 1
 1    NTRY = NTRY+1
C--search the decay products and decide which one to decay next
      IF(.NOT.DECSPN(IDEC)) THEN
        CALL HWDSI1(IDEC,IP)
      ELSE
        IDEC = JMOSPN(IDEC)
        GOTO 1
      ENDIF
C--first no more particles in this decay to develop so move up chain
      IF(IP.EQ.0) THEN
        IDEC = JMOSPN(IDEC)
C--reached the end of this spin chain go back to HWDHOB
        IF(IDEC.EQ.0) THEN
          NSPN = 0
          RETURN
C--otherwise keep going up the chain
        ELSE
          IF(NTRY.LE.NBTRY) THEN
            GOTO 1
          ELSE
            CALL HWWARN('HWDSIN',100,*999)
          ENDIF
        ENDIF
C--special for tau decays call spin correlation in tau decay routine
      ELSEIF(ABS(IDHEP(IDSPN(IP))).EQ.15) THEN
        CALL HWDSI3(IP)
        IF(IERROR.NE.0) RETURN
        GOTO 1
      ENDIF
C--work out where that particle is
      IHEP = IDSPN(IP)
C--if particle has daughters
 10   IF(JDAHEP(1,IHEP).NE.0) THEN
        IF(ISTHEP(IHEP).GE.141.AND.ISTHEP(IHEP).LE.144) THEN
          DO ID1=JDAHEP(1,IHEP),JDAHEP(2,IHEP)
            IF(IDHW(ID1).EQ.ID) IHEP=ID1
          ENDDO
        ELSE
          IHEP = JDAHEP(1,IHEP)
        ENDIF
      ENDIF
      IS=ISTHEP(IHEP)
      ID=IDHW(IHEP)
      NTRY = NTRY+1
      IF(NTRY.GE.NBTRY) CALL HWWARN('HWDSIN',101,*999)
      IF (.NOT.RSTAB(ID).AND.(ID.EQ.6.OR.ID.EQ.12.OR.
     & (ID.GE.203.AND.ID.LE.218).OR.ABS(IDPDG(ID)).GT.1000000).AND.
     & (IS.EQ.190.OR.(IS.GE.147.AND.IS.LE.151))) THEN
        CALL HWDHO1(IHEP,ID,IM,NPR,LHEP,MHEP)
        IF(IERROR.NE.0) RETURN
      ELSE
        GOTO 10
      ENDIF
C--perform the decay including spin correlations
      CALL HWDSI2(IHEP,IM,NPR,MHEP,LHEP,KHEP,PW)
      IF(IERROR.NE.0) RETURN
C--make the colour connections
      CALL HWDHO3(IHEP,ID,IM,NPR,MHEP,LHEP,KHEP,CLSAVE)
      IF (IERROR.NE.0) RETURN
C--perform the parton-showers
      CALL HWDHO4(IHEP,ID,IM,NPR,MHEP,LHEP,KHEP,PW)
      IF(IERROR.NE.0) RETURN
C--perform RPV colour connections
      CALL HWDHO5(IHEP,MHEP,LHEP,CLSAVE)
      IF(IERROR.NE.0) RETURN
C--continue and perform the next decay
      IDEC = IP
      IF(NTRY.LE.NBTRY) THEN
        GOTO 1
      ELSE
        CALL HWWARN('HWDSIN',102,*999)
      ENDIF
 999  END
CDECK  ID>, HWDSI1.
*CMZ :-        -30/09/02  14:05:28  by  Peter Richardson
*-- Author :    Peter Richardson
C-----------------------------------------------------------------------
      SUBROUTINE HWDSI1(IDEC,IP)
C-----------------------------------------------------------------------
C  Subroutine to check a vertex and decide which branch to treat
C-----------------------------------------------------------------------
#include "HERWIG65.INC"
      INTEGER IDEC,I,IPICK(5),IP,HWRINT,P1,P2,P3,P4,P3P,P4P,NPR,P0,P0P,
     &     P1P,P2P,IF1,IF2,P5,P5P
      DOUBLE PRECISION NORM
      DOUBLE COMPLEX RHOLP(2,2),RHOPS(2,2)
      EXTERNAL HWRINT
C--loop over the daughters and decide what to do
      IP = 0
C--if daughters of particle the same issue warning
      IF(JDASPN(1,IDEC).EQ.JDASPN(2,IDEC))
     &     CALL HWWARN('HWDSI1',100,*999)
C--loop over the decay products
      DO I=JDASPN(1,IDEC),JDASPN(2,IDEC)
        IF(.NOT.DECSPN(I)) THEN
C--first SM particles other than tau and top and stable particles
          IF(RSTAB(IDHW(IDSPN(I)))
     &    .OR.(IDHW(IDSPN(I)).LE.12.AND.ABS(IDHEP(IDSPN(I))).NE.6)
     &    .OR.(IDHW(IDSPN(I)).GE.121.AND.IDHW(IDSPN(I)).LE.132.AND.
     &          ABS(IDHEP(IDSPN(I))).NE.15)) THEN
             DECSPN(I) = .TRUE.
             RHOSPN(1,1,I) = HALF
             RHOSPN(1,2,I) = ZERO
             RHOSPN(2,1,I) = ZERO
             RHOSPN(2,2,I) = HALF
C--spinless particles
          ELSEIF(RSPIN(IDHW(IDSPN(I))).EQ.ZERO) THEN
             DECSPN(I) = .TRUE.
             RHOSPN(1,1,I) = ONE
             RHOSPN(1,2,I) = ZERO
             RHOSPN(2,1,I) = ZERO
             RHOSPN(2,2,I) = ZERO
          ELSE
C--particle which needs development
            IP = IP+1
            IPICK(IP) = I
          ENDIF
        ENDIF
      ENDDO
C--pick the particle to decay next
      IF(IP.EQ.0) THEN
        IF(JMOSPN(IDEC).EQ.0) RETURN
C--done everything compute the decay matrix and move up
        DECSPN(IDEC) = .TRUE.
        NPR = JDASPN(2,IDEC)-JDASPN(1,IDEC)+1
        DO 20 P0=1,2
        DO 20 P0P=1,2
 20     RHOSPN(P0,P0P,IDEC) = ZERO
C--two body decay
        IF(NPR.EQ.2) THEN
          DO 21 P0 =1,2
          DO 21 P0P=1,2
          DO 21 P1 =1,2
          DO 21 P1P=1,2
          DO 21 P2 =1,2
          DO 21 P2P=1,2
 21       RHOSPN(P0,P0P,IDEC) = RHOSPN(P0,P0P,IDEC)+
     &              MESPN(P0 ,P1 ,P2 ,1,NCFL(IDEC),IDEC)*
     &       DCONJG(MESPN(P0P,P1P,P2P,1,NCFL(IDEC),IDEC))*
     &       RHOSPN(P1,P1P,JDASPN(1,IDEC))*RHOSPN(P2,P2P,JDASPN(2,IDEC))
C--three body decay
        ELSEIF(NPR.EQ.3) THEN
          DO 25 P0 =1,2
          DO 25 P0P=1,2
          DO 25 P1 =1,2
          DO 25 P1P=1,2
          DO 25 P2 =1,2
          DO 25 P2P=1,2
          DO 25 P3 =1,2
          DO 25 P3P=1,2
 25       RHOSPN(P0,P0P,IDEC) = RHOSPN(P0,P0P,IDEC)+
     &           MESPN(P0 ,P1 ,P2 ,P3 ,NCFL(IDEC),IDEC)*
     &    DCONJG(MESPN(P0P,P1P,P2P,P3P,NCFL(IDEC),IDEC))*
     &    RHOSPN(P1,P1P,JDASPN(1,IDEC))*RHOSPN(P2,P2P,JDASPN(1,IDEC)+1)*
     &    RHOSPN(P3,P3P,JDASPN(2,IDEC))
C--higher
        ELSE
          CALL HWWARN('HWDSI1',500,*999)
        ENDIF
C--now normalise this
        NORM = DBLE(RHOSPN(1,1,IDEC))+DBLE(RHOSPN(2,2,IDEC))
        IF(NORM.GT.ZERO) THEN
          NORM = ONE/NORM
          DO 35 P0=1,2
          DO 35 P0P=1,2
 35       RHOSPN(P0,P0P,IDEC) = NORM*RHOSPN(P0,P0P,IDEC)
        ELSE
          CALL HWWARN('HWDSI1',101,*999)
        ENDIF
      ELSE
C--pick the particle to be decayed
        IP = IPICK(HWRINT(1,IP))
C--setup the spin density matrix for the decay
C--special for the hard process
        IF(ISTHEP(IDSPN(IDEC)).EQ.120) THEN
          NPR = JDASPN(2,IDEC)-JDASPN(1,IDEC)+1
C--set up the spin density matrices for the incoming partons
C--zero off diagonal elements
          RHOLP(2,1) = ZERO
          RHOLP(1,2) = ZERO
          RHOPS(2,1) = ZERO
          RHOPS(1,2) = ZERO
C--set up for polarized incoming beams in lepton collisons
          IF(IDHW(JMOHEP(1,IDSPN(IDEC))).GE.121.AND.
     &       IDHW(JMOHEP(1,IDSPN(IDEC))).LE.132) THEN
            RHOLP(1,1) = HALF*(ONE+EPOLN(3))
            RHOLP(2,2) = HALF*(ONE-EPOLN(3))
            RHOPS(1,1) = HALF*(ONE-PPOLN(3))
            RHOPS(2,2) = HALF*(ONE+PPOLN(3))
C--otherwise average
          ELSE
            RHOLP(1,1) = HALF
            RHOLP(2,2) = HALF
            RHOPS(1,1) = HALF
            RHOPS(2,2) = HALF
          ENDIF
C--first decay product
          IF(NPR.EQ.2) THEN
           IF(IP.EQ.JDASPN(1,IDEC)) THEN
C--if using first colour flow option
            IF(SPCOPT.EQ.1) THEN
              DO 5 P3 =1,2
              DO 5 P3P=1,2
              RHOSPN(P3,P3P,IP) = ZERO
              DO 5 IF1=1,NCFL(1)
              DO 5 IF2=1,NCFL(1)
              DO 5 P1 =1,2
              DO 5 P1P=1,2
              DO 5 P2 =1,2
              DO 5 P2P=1,2
              DO 5 P4 =1,2
              DO 5 P4P=1,2
 5            RHOSPN(P3,P3P,IP) = RHOSPN(P3,P3P,IP)+SPNCFC(IF1,IF2,1)*
     &               MESPN(P1 ,P2 ,P3 ,P4 ,IF1,1)*
     &        DCONJG(MESPN(P1P,P2P,P3P,P4P,IF2,1))*
     &        RHOLP(P1,P1P)*RHOPS(P2,P2P)*RHOSPN(P4,P4P,IP+1)
C--if using second colour flow option
            ELSEIF(SPCOPT.EQ.2) THEN
              DO 6 P3 =1,2
              DO 6 P3P=1,2
              RHOSPN(P3,P3P,IP) = ZERO
              DO 6 P1 =1,2
              DO 6 P1P=1,2
              DO 6 P2 =1,2
              DO 6 P2P=1,2
              DO 6 P4 =1,2
              DO 6 P4P=1,2
 6            RHOSPN(P3,P3P,IP) = RHOSPN(P3,P3P,IP)
     &                +SPNCFC(NCFL(1),NCFL(1),1)*
     &               MESPN(P1 ,P2 ,P3 ,P4 ,NCFL(1),1)*
     &        DCONJG(MESPN(P1P,P2P,P3P,P4P,NCFL(1),1))*
     &        RHOLP(P1,P1P)*RHOPS(P2,P2P)*RHOSPN(P4,P4P,IP+1)
C--unknown option issue warning
            ELSE
              CALL HWWARN('HWDSI1',501,*999)
            ENDIF
C--second decay product
           ELSE
            IF(SPCOPT.EQ.1) THEN
              DO 10 P4 =1,2
              DO 10 P4P=1,2
              RHOSPN(P4,P4P,IP) = (0.0D0,0.0D0)
              DO 10 IF1=1,NCFL(1)
              DO 10 IF2=1,NCFL(1)
              DO 10 P1 =1,2
              DO 10 P1P=1,2
              DO 10 P2 =1,2
              DO 10 P2P=1,2
              DO 10 P3 =1,2
              DO 10 P3P=1,2
 10           RHOSPN(P4,P4P,IP) = RHOSPN(P4,P4P,IP)+SPNCFC(IF1,IF2,1)*
     &                 MESPN(P1 ,P2 ,P3 ,P4 ,IF1,1)*
     &          DCONJG(MESPN(P1P,P2P,P3P,P4P,IF2,1))*
     &          RHOLP(P1,P1P)*RHOPS(P2,P2P)*RHOSPN(P3,P3P,IP-1)
            ELSEIF(SPCOPT.EQ.2) THEN
              DO 11 P4 =1,2
              DO 11 P4P=1,2
              RHOSPN(P4,P4P,IP) = (0.0D0,0.0D0)
              DO 11 P1 =1,2
              DO 11 P1P=1,2
              DO 11 P2 =1,2
              DO 11 P2P=1,2
              DO 11 P3 =1,2
              DO 11 P3P=1,2
 11           RHOSPN(P4,P4P,IP) = RHOSPN(P4,P4P,IP)
     &                +SPNCFC(NCFL(1),NCFL(1),1)*
     &                 MESPN(P1 ,P2 ,P3 ,P4 ,NCFL(1),1)*
     &          DCONJG(MESPN(P1P,P2P,P3P,P4P,NCFL(1),1))*
     &          RHOLP(P1,P1P)*RHOPS(P2,P2P)*RHOSPN(P3,P3P,IP-1)
            ELSE
              CALL HWWARN('HWDSI1',502,*999)
            ENDIF
           ENDIF
C--new for four body gauge boson pair processes
          ELSEIF(NPR.EQ.4) THEN
C--first particle
           IF(IP.EQ.JDASPN(1,IDEC)) THEN
             DO 41 P1 =1,2
             DO 41 P1P=1,2
             RHOSPN(P1,P1P,IP) = (0.0D0,0.0D0)
             DO 41 P3 =1,2
             DO 41 P3P=1,2
             DO 41 P5 =1,2
             DO 41 P5P=1,2
 41          RHOSPN(P1,P1P,IP) = RHOSPN(P1,P1P,IP)+
     &       MESPN(P5,P1,P3,1,1,1)*DCONJG(MESPN(P5P,P1P,P3P,1,1,1))*
     &       RHOSPN(P1,P1P,JDASPN(1,IDEC)+1)*
     &       RHOSPN(P3,P3P,JDASPN(1,IDEC)+2)*
     &       RHOSPN(P3,P3P,JDASPN(2,IDEC))
C--second particle
           ELSEIF(IP.EQ.JDASPN(1,IDEC)+1) THEN
             DO 42 P1 =1,2
             DO 42 P1P=1,2
             RHOSPN(P1,P1P,IP) = (0.0D0,0.0D0)
             DO 42 P3 =1,2
             DO 42 P3P=1,2
             DO 42 P5 =1,2
             DO 42 P5P=1,2
 42          RHOSPN(P1,P1P,IP) = RHOSPN(P1,P1P,IP)+
     &       MESPN(P5,P1,P3,1,1,1)*DCONJG(MESPN(P5P,P1P,P3P,1,1,1))*
     &       RHOSPN(P1,P1P,JDASPN(1,IDEC))*
     &       RHOSPN(P3,P3P,JDASPN(1,IDEC)+2)*
     &       RHOSPN(P3,P3P,JDASPN(2,IDEC))
C--third particle
           ELSEIF(IP.EQ.JDASPN(1,IDEC)+2) THEN
             DO 43 P3 =1,2
             DO 43 P3P=1,2
             RHOSPN(P3,P3P,IP) = (0.0D0,0.0D0)
             DO 43 P1 =1,2
             DO 43 P1P=1,2
             DO 43 P5 =1,2
             DO 43 P5P=1,2
 43          RHOSPN(P3,P3P,IP) = RHOSPN(P3,P3P,IP)+
     &       MESPN(P5,P1,P3,1,1,1)*DCONJG(MESPN(P5P,P1P,P3P,1,1,1))*
     &       RHOSPN(P1,P1P,JDASPN(1,IDEC))*
     &       RHOSPN(P1,P1P,JDASPN(1,IDEC)+1)*
     &       RHOSPN(P3,P3P,JDASPN(2,IDEC))
C--fourth particle
           ELSEIF(IP.EQ.JDASPN(2,IDEC)) THEN
             DO 44 P3 =1,2
             DO 44 P3P=1,2
             RHOSPN(P3,P3P,IP) = (0.0D0,0.0D0)
             DO 44 P1 =1,2
             DO 44 P1P=1,2
             DO 44 P5 =1,2
             DO 44 P5P=1,2
 44          RHOSPN(P3,P3P,IP) = RHOSPN(P3,P3P,IP)+
     &       MESPN(P5,P1,P3,1,1,1)*DCONJG(MESPN(P5P,P1P,P3P,1,1,1))*
     &       RHOSPN(P1,P1P,JDASPN(1,IDEC))*
     &       RHOSPN(P1,P1P,JDASPN(1,IDEC)+1)*
     &       RHOSPN(P3,P3P,JDASPN(1,IDEC)+2)
C--unrecognized issue warning
           ELSE
             CALL HWWARN('(HWDSI1)',509,*999)
           ENDIF
C--unrecognized issue warning
          ELSE
            CALL HWWARN('(HWDSI1)',508,*999)
          ENDIF
        ELSE
          NPR = JDASPN(2,IDEC)-JDASPN(1,IDEC)+1
          DO 50 P1 =1,2
          DO 50 P1P=1,2
 50       RHOSPN(P1,P1P,IP) = (0.0D0,0.0D0)
C--set-up matrix for 2-body decay
          IF(NPR.EQ.2) THEN
            IF(NCFL(IDEC).NE.1) CALL HWWARN('HWDSI1',503,*999)
            IF(IP.EQ.JDASPN(1,IDEC)) THEN
              DO 60 P0 =1,2
              DO 60 P0P=1,2
              DO 60 P1 =1,2
              DO 60 P1P=1,2
              DO 60 P2 =1,2
              DO 60 P2P=1,2
 60           RHOSPN(P1,P1P,IP) = RHOSPN(P1,P1P,IP)+RHOSPN(P0,P0P,IDEC)*
     &               MESPN(P0 ,P1 ,P2 ,1,1,IDEC)*
     &        DCONJG(MESPN(P0P,P1P,P2P,1,1,IDEC))*
     &        RHOSPN(P2,P2P,JDASPN(2,IDEC))
            ELSE
              DO 70 P0 =1,2
              DO 70 P0P=1,2
              DO 70 P1 =1,2
              DO 70 P1P=1,2
              DO 70 P2 =1,2
              DO 70 P2P=1,2
 70           RHOSPN(P2,P2P,IP) = RHOSPN(P2,P2P,IP)+RHOSPN(P0,P0P,IDEC)*
     &               MESPN(P0 ,P1 ,P2 ,1,1,IDEC)*
     &        DCONJG(MESPN(P0P,P1P,P2P,1,1,IDEC))*
     &        RHOSPN(P1,P1P,JDASPN(1,IDEC))
            ENDIF
C--set-up matrix for 3-body decay
          ELSEIF(NPR.EQ.3) THEN
            IF(SPCOPT.NE.2.AND.NCFL(IDEC).NE.1)
     &        CALL HWWARN('HWDSI1',504,*999)
C--first particle
            IF(IP.EQ.JDASPN(1,IDEC)) THEN
              DO 100 P0 =1,2
              DO 100 P0P=1,2
              DO 100 P1 =1,2
              DO 100 P1P=1,2
              DO 100 P2 =1,2
              DO 100 P2P=1,2
              DO 100 P3 =1,2
              DO 100 P3P=1,2
 100          RHOSPN(P1,P1P,IP) = RHOSPN(P1,P1P,IP)+RHOSPN(P0,P0P,IDEC)*
     &               MESPN(P0 ,P1 ,P2 ,P3 ,NCFL(IDEC),IDEC)*
     &        DCONJG(MESPN(P0P,P1P,P2P,P3P,NCFL(IDEC),IDEC))*
     &        RHOSPN(P2,P2P,JDASPN(1,IDEC)+1)*
     &        RHOSPN(P3,P3P,JDASPN(2,IDEC))
C--second particle
            ELSEIF(IP.EQ.JDASPN(1,IDEC)+1) THEN
              DO 105 P0 =1,2
              DO 105 P0P=1,2
              DO 105 P1 =1,2
              DO 105 P1P=1,2
              DO 105 P2 =1,2
              DO 105 P2P=1,2
              DO 105 P3 =1,2
              DO 105 P3P=1,2
 105          RHOSPN(P2,P2P,IP) = RHOSPN(P2,P2P,IP)+RHOSPN(P0,P0P,IDEC)*
     &               MESPN(P0 ,P1 ,P2 ,P3 ,NCFL(IDEC),IDEC)*
     &        DCONJG(MESPN(P0P,P1P,P2P,P3P,NCFL(IDEC),IDEC))*
     &        RHOSPN(P1,P1P,JDASPN(1,IDEC))*
     &        RHOSPN(P3,P3P,JDASPN(2,IDEC))
C--third particle
            ELSEIF(IP.EQ.JDASPN(2,IDEC)) THEN
              DO 110 P0 =1,2
              DO 110 P0P=1,2
              DO 110 P1 =1,2
              DO 110 P1P=1,2
              DO 110 P2 =1,2
              DO 110 P2P=1,2
              DO 110 P3 =1,2
              DO 110 P3P=1,2
 110          RHOSPN(P3,P3P,IP) = RHOSPN(P3,P3P,IP)+RHOSPN(P0,P0P,IDEC)*
     &               MESPN(P0 ,P1 ,P2 ,P3 ,NCFL(IDEC),IDEC)*
     &        DCONJG(MESPN(P0P,P1P,P2P,P3P,NCFL(IDEC),IDEC))*
     &        RHOSPN(P1,P1P,JDASPN(1,IDEC))*
     &        RHOSPN(P2,P2P,JDASPN(1,IDEC)+1)
C--unrecognized
            ELSE
              CALL HWWARN('HWDSI1',102,*999)
            ENDIF
          ELSEIF(NPR.EQ.4) THEN
C--first particle
            IF(IP.EQ.JDASPN(1,IDEC)) THEN
              DO 151 P1 =1,2
              DO 151 P1P=1,2
              RHOSPN(P1,P1P,IP) = (0.0D0,0.0D0)
              DO 151 P2 =1,2
              DO 151 P2P=1,2
              DO 151 P3 =1,2
              DO 151 P3P=1,2
              DO 151 P4 =1,2
              DO 151 P4P=1,2
 151          RHOSPN(P1,P1P,IP) = RHOSPN(P1,P1P,IP)+
     &                  MESPN(P1 ,P2 ,P3 ,P4 ,1,1)*
     &           DCONJG(MESPN(P1P,P2P,P3P,P4P,1,1))*
     &           RHOSPN(P2,P2P,JDASPN(1,IDEC)+1)*
     &           RHOSPN(P3,P3P,JDASPN(1,IDEC)+2)*
     &           RHOSPN(P4,P4P,JDASPN(2,IDEC))
C--second particle
            ELSEIF(IP.EQ.JDASPN(1,IDEC)+1) THEN
              DO 152 P2 =1,2
              DO 152 P2P=1,2
              RHOSPN(P2,P2P,IP) = (0.0D0,0.0D0)
              DO 152 P1 =1,2
              DO 152 P1P=1,2
              DO 152 P3 =1,2
              DO 152 P3P=1,2
              DO 152 P4 =1,2
              DO 152 P4P=1,2
 152             RHOSPN(P2,P2P,IP) = RHOSPN(P2,P2P,IP)+
     &                  MESPN(P1 ,P2 ,P3 ,P4 ,1,1)*
     &           DCONJG(MESPN(P1P,P2P,P3P,P4P,1,1))*
     &           RHOSPN(P1,P1P,JDASPN(1,IDEC))*
     &           RHOSPN(P3,P3P,JDASPN(1,IDEC)+2)*
     &           RHOSPN(P4,P4P,JDASPN(2,IDEC))
C--third particle
            ELSEIF(IP.EQ.JDASPN(1,IDEC)+2) THEN
              DO 153 P3 =1,2
              DO 153 P3P=1,2
              RHOSPN(P3,P3P,IP) = (0.0D0,0.0D0)
              DO 153 P1 =1,2
              DO 153 P1P=1,2
              DO 153 P2 =1,2
              DO 153 P2P=1,2
              DO 153 P4 =1,2
              DO 153 P4P=1,2
 153          RHOSPN(P3,P3P,IP) = RHOSPN(P3,P3P,IP)+
     &                  MESPN(P1 ,P2 ,P3 ,P4 ,1,1)*
     &           DCONJG(MESPN(P1P,P2P,P3P,P4P,1,1))*
     &           RHOSPN(P1,P1P,JDASPN(1,IDEC))*
     &           RHOSPN(P2,P2P,JDASPN(1,IDEC)+1)*
     &           RHOSPN(P4,P4P,JDASPN(2,IDEC))
C--fourth particle
            ELSEIF(IP.EQ.JDASPN(2,IDEC)) THEN
              DO 154 P4 =1,2
              DO 154 P4P=1,2
              RHOSPN(P4,P4P,IP) = (0.0D0,0.0D0)
              DO 154 P1 =1,2
              DO 154 P1P=1,2
              DO 154 P2 =1,2
              DO 154 P2P=1,2
              DO 154 P3 =1,2
              DO 154 P3P=1,2
 154          RHOSPN(P4,P4P,IP) = RHOSPN(P4,P4P,IP)+
     &                  MESPN(P1 ,P2 ,P3 ,P4 ,1,1)*
     &           DCONJG(MESPN(P1P,P2P,P3P,P4P,1,1))*
     &           RHOSPN(P1,P1P,JDASPN(1,IDEC))*
     &           RHOSPN(P2,P2P,JDASPN(1,IDEC)+1)*
     &           RHOSPN(P3,P3P,JDASPN(1,IDEC)+2)
            ELSE
              CALL HWWARN('HWDSI1',505,*999)
            ENDIF
          ELSE
            CALL HWWARN('HWDSI1',506,*999)
          ENDIF
        ENDIF
C--normalise the spin density matrix
        NORM = DBLE(RHOSPN(1,1,IP))+DBLE(RHOSPN(2,2,IP))
        IF(NORM.GT.ZERO) THEN
          NORM = ONE/NORM
          DO 15 P3=1,2
          DO 15 P3P=1,2
 15       RHOSPN(P3,P3P,IP) = NORM*RHOSPN(P3,P3P,IP)
        ELSE
          CALL HWWARN('HWDSI1',107,*999)
        ENDIF
      ENDIF
 999  END
CDECK  ID>, HWDSI2.
*CMZ :-        -30/09/02  14:05:28  by  Peter Richardson
*-- Author :    Peter Richardson
C-----------------------------------------------------------------------
      SUBROUTINE HWDSI2(IHEP,IM,NPR,MHEP,LHEP,KHEP,PW)
C-----------------------------------------------------------------------
C   Subroutine to perform the second part of the heavy object decays
C   IE generate the kinematics for the decay
C   including spin correlations
C   was part of HWDHOB
C-----------------------------------------------------------------------
#include "HERWIG65.INC"
      DOUBLE PRECISION HWRGEN,PW(5),HWDPWT,HWDWWT,PCM,HWUPCM
      INTEGER IHEP,IM,KHEP,LHEP,MHEP,NPR,ISN,RHEP
      EXTERNAL HWRGEN,HWDPWT,HWDWWT,HWUPCM
      IF (IERROR.NE.0) RETURN
      ISN = ISNHEP(IHEP)
      IF (NPR.EQ.2) THEN
C Two body decay: LHEP -> MHEP + NHEP
        IF(NME(IM).GT.20000.AND.NME(IM).LT.30000) THEN
C--generate a two body decay to a gauge boson as a three body decay
          CALL HWDSM3(2,IHEP,MHEP,NHEP,0,NME(IM)-20000,
     &                 RHOSPN(1,1,ISN),ISN)
C--two body decay
        ELSEIF(NME(IM).GT.30000.AND.NME(IM).LT.40000) THEN
          CALL HWDSM2(IHEP,MHEP,NHEP,NME(IM)-30000,
     &          RHOSPN(1,1,ISN),ISN)
C--otherwise issue warning
C--change by PR 9/30/02 to issue non-terminal warning and continue
        ELSE
          CALL HWWARN('HWDSI2',1,*999)
          PCM=HWUPCM(PHEP(5,IHEP),PHEP(5,MHEP),PHEP(5,NHEP))
          CALL HWDTWO(PHEP(1,IHEP),PHEP(1,MHEP),
     &                PHEP(1,NHEP),PCM,TWO,.FALSE.)
          DECSPN(ISN) = .TRUE.
          IF(RSPIN(IDHW(IHEP)).EQ.ZERO) THEN
            RHOSPN(1,1,ISN) = ONE
            RHOSPN(1,2,ISN) = ZERO
            RHOSPN(2,1,ISN) = ZERO
            RHOSPN(2,2,ISN) = ZERO
          ELSE
            RHOSPN(1,1,ISN) = HALF
            RHOSPN(1,2,ISN) = ZERO
            RHOSPN(2,1,ISN) = ZERO
            RHOSPN(2,2,ISN) = HALF
          ENDIF
        ENDIF
      ELSEIF (NPR.EQ.3) THEN
C Three body decay: LHEP -> KHEP + MHEP + NHEP
        KHEP=MHEP
        MHEP=MHEP+1
C Provisional colour self-connection of KHEP
        JMOHEP(2,KHEP)=KHEP
        JDAHEP(2,KHEP)=KHEP
C--if old codes issue warning
        IF (NME(IM).EQ.100.OR.NME(IM).EQ.200.OR.NME(IM).EQ.300) THEN
          CALL HWWARN('HWDSI2',502,*999)
C--three body spin matrix element
        ELSEIF(NME(IM).GE.10000.AND.NME(IM).LT.20000) THEN
          CALL HWDSM3(3,IHEP,MHEP,KHEP,NHEP,NME(IM)-10000,
     &            RHOSPN(1,1,ISN),ISN)
C--special for top decay
          IF(ABS(IDHEP(IHEP)).EQ.6) THEN
            CALL HWVSUM(4,PHEP(1,KHEP),PHEP(1,MHEP),PW)
            CALL HWUMAS(PW)
          ENDIF
C--unknown issue warning
        ELSE
          CALL HWWARN('HWDSI2',2,*999)
C Three body phase space decay
          CALL HWDTHR(PHEP(1,IHEP),PHEP(1,MHEP),
     &                PHEP(1,KHEP),PHEP(1,NHEP),HWDPWT)
          CALL HWVEQU(4,VHEP(1,KHEP),VHEP(1,MHEP))
          DECSPN(ISN) = .TRUE.
          IF(RSPIN(IDHW(IHEP)).EQ.ZERO) THEN
            RHOSPN(1,1,ISN) = ONE
            RHOSPN(1,2,ISN) = ZERO
            RHOSPN(2,1,ISN) = ZERO
            RHOSPN(2,2,ISN) = ZERO
          ELSE
            RHOSPN(1,1,ISN) = HALF
            RHOSPN(1,2,ISN) = ZERO
            RHOSPN(2,1,ISN) = ZERO
            RHOSPN(2,2,ISN) = HALF
          ENDIF
        ENDIF
      ELSEIF(NPR.EQ.4) THEN
        CALL HWWARN('HWDSI2',3,*999)
C Four body decay: LHEP -> KHEP + RHEP + MHEP + NHEP
        KHEP = MHEP
        RHEP = MHEP+1
        MHEP = MHEP+2
        ISTHEP(NHEP) = 114
C Provisional colour connections of KHEP and RHEP
        JMOHEP(2,KHEP)=RHEP
        JDAHEP(2,KHEP)=RHEP
        JMOHEP(2,RHEP)=KHEP
        JDAHEP(2,RHEP)=KHEP
C Four body phase space decay
        CALL HWDFOR(PHEP(1,IHEP),PHEP(1,KHEP),PHEP(1,RHEP),
     &                PHEP(1,MHEP),PHEP(1,NHEP))
        CALL HWVEQU(4,VHEP(1,KHEP),VHEP(1,RHEP))
        CALL HWVEQU(4,VHEP(1,KHEP),VHEP(1,MHEP))
        DECSPN(ISN) = .TRUE.
        IF(RSPIN(IDHW(IHEP)).EQ.ZERO) THEN
          RHOSPN(1,1,ISN) = ONE
          RHOSPN(1,2,ISN) = ZERO
          RHOSPN(2,1,ISN) = ZERO
          RHOSPN(2,2,ISN) = ZERO
        ELSE
          RHOSPN(1,1,ISN) = HALF
          RHOSPN(1,2,ISN) = ZERO
          RHOSPN(2,1,ISN) = ZERO
          RHOSPN(2,2,ISN) = HALF
        ENDIF
      ELSE
        CALL HWWARN('HWDSI2',100,*999)
      ENDIF
 999  END
CDECK  ID>, HWDSI3.
*CMZ :-        -30/09/02  14:05:28  by  Peter Richardson
*-- Author :    Peter Richardson
C-----------------------------------------------------------------------
      SUBROUTINE HWDSI3(IP)
C-----------------------------------------------------------------------
C     Subroutine to handle spin correlations in tau decays
C     averages spin if not using TAUOLA
C     if using TAUOLA selects the spin and uses TAUOLA to perform the
C     decay
C-----------------------------------------------------------------------
#include "HERWIG65.INC"
      INTEGER IP,IHEP,ID1,ID,NTRY
      DOUBLE PRECISION PPOL,HWRGEN,POL
      EXTERNAL HWRGEN
C--if HERWIG is performing tau decays average over spins and return
C--spin averaged tau decay will be done later
      IF(TAUDEC.EQ.'HERWIG') THEN
        DECSPN(IP) = .TRUE.
        RHOSPN(1,1,IP) = HALF
        RHOSPN(2,1,IP) = ZERO
        RHOSPN(1,2,IP) = ZERO
        RHOSPN(2,2,IP) = HALF
C--if using tauola select the polarization for the decay
      ELSEIF(TAUDEC.EQ.'TAUOLA') THEN
C--work out where that particle is
        IHEP = IDSPN(IP)
        NTRY = 0
 10     ID   = IDHW(IHEP)
        IF(JDAHEP(1,IHEP).NE.0) THEN
          IF(ISTHEP(IHEP).GE.141.AND.ISTHEP(IHEP).LE.144) THEN
            DO ID1=JDAHEP(1,IHEP),JDAHEP(2,IHEP)
              IF(IDHW(ID1).EQ.ID) IHEP=ID1
            ENDDO
          ELSE
            IHEP = JDAHEP(1,IHEP)
          ENDIF
          NTRY = NTRY+1
          IF(NTRY.LT.NBTRY) THEN
            GOTO 10
          ELSE
            CALL HWWARN('HWDSI3',100,*999)
          ENDIF
        ENDIF
C--select the tau polarization
        PPOL = DBLE(RHOSPN(1,1,IP))
        IF(PPOL.GE.HWRGEN(0)) THEN
          POL = 1.0D0
          RHOSPN(1,1,IP) =  ONE
          RHOSPN(2,1,IP) = ZERO
          RHOSPN(1,2,IP) = ZERO
          RHOSPN(2,2,IP) = ZERO
        ELSE
          POL =-1.0D0
          RHOSPN(1,1,IP) = ZERO
          RHOSPN(2,1,IP) = ZERO
          RHOSPN(1,2,IP) = ZERO
          RHOSPN(2,2,IP) =  ONE
        ENDIF
C--decay the particle
        CALL HWDTAU(1,IHEP,POL)
        DECSPN(IP) = .TRUE.
      ELSE
        CALL HWWARN('HWDSI3',500,*999)
      ENDIF
 999  END
CDECK  ID>, HWDSM2.
*CMZ :-        -09/04/02  13:46:07  by  Peter Richardson
*-- Author :    Peter Richardson
C-----------------------------------------------------------------------
      SUBROUTINE HWDSM2(ID,IOUT1,IOUT2,IMODE,RHOIN,IDSPIN)
C-----------------------------------------------------------------------
C  Subroutine to calculate the two body matrix element for spin
C  correlations
C-----------------------------------------------------------------------
#include "HERWIG65.INC"
      INTEGER IOUT1,IOUT2,IMODE,IDSPIN,ID,I,J,IDP(3),P0,P1,P2,O(2),P0P,
     &     NTRY
      DOUBLE PRECISION XMASS,PLAB,PRW,PCM,PREF(5),P(5,3),PM(5,3),PCMA,
     &     HWUPCM,MA(3),MA2(3),HWULDO,PP,HWVDOT,N(3),EPS,PRE,PHS,A(2),
     &     WGT,WTMAX,HWRGEN
      DOUBLE COMPLEX RHOIN(2,2),S,D,ME(2,2,2),F1(2,2,8),F0(2,2,8),
     &     F2M(2,2,8),F1M(2,2,8),F1F(2,2,8),F2(2,2,8,8),F0B(2,2,8,8)
      COMMON/HWHEWP/XMASS(10),PLAB(5,10),PRW(5,2),PCM(5,10)
      DATA PREF/1.0D0,0.0D0,0.0D0,1.0D0,0.0D0/
      DATA O/2,1/
      COMMON/HWHEWS/S(8,8,2),D(8,8)
      PARAMETER(EPS=1D-20)
      EXTERNAL HWUPCM,HWULDO,HWVDOT,HWRGEN
C--first setup if this is the start of a new spin chain
      IF(NSPN.EQ.0) THEN
C--zero the elements of the array
        CALL HWVZRI(  NMXHEP,ISNHEP)
        CALL HWVZRI(  NMXSPN,JMOSPN)
        CALL HWVZRI(2*NMXSPN,JDASPN)
        CALL HWVZRI(  NMXSPN, IDSPN)
        NSPN = NSPN+1
        JMOSPN(NSPN) = 0
        IDSPN (NSPN) = ID
        DECSPN(NSPN) = .FALSE.
        IF(RSPIN(IDHW(ID)).EQ.ZERO) THEN
          RHOSPN(1,1,NSPN) = ONE
          RHOSPN(2,1,NSPN) = ZERO
          RHOSPN(1,2,NSPN) = ZERO
          RHOSPN(2,2,NSPN) = ZERO
        ELSE
          RHOSPN(1,1,NSPN) = HALF
          RHOSPN(2,1,NSPN) = ZERO
          RHOSPN(1,2,NSPN) = ZERO
          RHOSPN(2,2,NSPN) = HALF
        ENDIF
        ISNHEP(ID)    = NSPN
      ENDIF
C--MA is mass for this decay (OFF-SHELL)
C--generate the momenta for a two body decay
      P(5,1) = PHEP(5,   ID)
      P(5,2) = PHEP(5,IOUT1)
      P(5,3) = PHEP(5,IOUT2)
      IDP(1) = IDHW(ID)
      IDP(2) = IDHW(IOUT1)
      IDP(3) = IDHW(IOUT2)
      DO 1 I=1,3
      MA(I)  = P(5,I)
 1    MA2(I) = MA(I)**2
      PCMA   = HWUPCM(P(5,1),P(5,2),P(5,3))
C--setup the couplings
      DO 2 I=1,2
 2    A(I) = A2MODE(I,IMODE)
C--phase space factor
      PHS = PCMA/MA2(1)/8.0D0/PIFAC
C--maximum weight
      WTMAX = WT2MAX(IMODE)
      NTRY = 0
 1000 NTRY = NTRY+1
      CALL HWVEQU(5,PHEP(1,ID),P(1,1))
      CALL HWDTWO(P(1,1),P(1,2),P(1,3),PCMA,2.0D0,.TRUE.)
      DO 3 I=1,3
C--compute the references vectors
C--not important if SM particle which can't have spin measured
C--ie anything other the top and tau
C--also not important if particle is approx massless
C--first the SM particles other than top and tau
      IF(IDP(I).LT.400.AND.(IDP(I).NE.6.AND.IDP(I).NE.12
     &                .AND.IDP(I).NE.125.AND.IDP(I).NE.131)) THEN
        CALL HWVEQU(5,PREF,PLAB(1,I+3))
C--all other particles
      ELSE
        PP = SQRT(HWVDOT(3,P(1,I),P(1,I)))
        CALL HWVSCA(3,ONE/PP,P(1,I),N)
        PLAB(4,I+3) = HALF*(P(4,I)-PP)
        PP = HALF*(PP-MA(I)-PP**2/(MA(I)+P(4,I)))
        CALL HWVSCA(3,PP,N,PLAB(1,I+3))
        CALL HWUMAS(PLAB(1,I+3))
        PP = HWVDOT(3,PLAB(1,I+3),PLAB(1,I+3))
C--fix to avoid problems if approx massless due to energy
        IF(PP.LT.EPS) CALL HWVEQU(5,PREF,PLAB(1,I+3))
      ENDIF
C--now the massless vectors
      PP = HALF*P(5,I)**2/HWULDO(PLAB(1,I+3),P(1,I))
      DO 4 J=1,4
 4    PLAB(J,I) = P(J,I)-PP*PLAB(J,I+3)
 3    CALL HWUMAS(PLAB(1,I))
C--change order of momenta for call to HE code
      DO 5 I=1,3
      PM(1,I) = P(3,I)
      PM(2,I) = P(1,I)
      PM(3,I) = P(2,I)
      PM(4,I) = P(4,I)
 5    PM(5,I) = P(5,I)
      DO 6 I=1,6
      PCM(1,I)=PLAB(3,I)
      PCM(2,I)=PLAB(1,I)
      PCM(3,I)=PLAB(2,I)
      PCM(4,I)=PLAB(4,I)
 6    PCM(5,I)=PLAB(5,I)
C--compute the S functions
      CALL HWHEW2(6,PCM(1,1),S(1,1,2),S(1,1,1),D)
      DO 7 I=1,6
      DO 7 J=1,6
      S(I,J,2) = -S(I,J,2)
 7    D(I,J)   = TWO*D(I,J)
C--now compute the F functions needed
      CALL HWH2F2(6,F1 ,5,PM(1,2), MA(2))
      CALL HWH2F2(6,F0 ,4,PM(1,1), MA(1))
      CALL HWH2F2(6,F1M,5,PM(1,2),-MA(2))
      CALL HWH2F2(6,F2M,6,PM(1,3),-MA(3))
      CALL HWH2F1(6,F1F,5,PM(1,2), MA(2))
      CALL HWH2F3(6,F2   ,PM(1,3),ZERO  )
      CALL HWH2F3(6,F0B  ,PM(1,1),ZERO  )
C--now compute the diagrams
C--fermion --> fermion scalar
      IF(I2DRTP(IMODE).EQ.1) THEN
        PRE = HWULDO(PM(1,1),PCM(1,4))*HWULDO(PM(1,2),PCM(1,5))
        PRE = HALF/SQRT(PRE)
        DO 10 P0=1,2
        DO 10 P1=1,2
        ME(P0,P1,2) = (0.0D0,0.0D0)
 10     ME(P0,P1,1) = PRE*( A(O(P1))*S(5,2,O(P1))*F0(  P1 ,O(P0),2)
     &                     +A(  P1 )*MA(2)*       F0(O(P1),O(P0),5))
C--fermion --> scalar fermion   diagrams
      ELSEIF(I2DRTP(IMODE).EQ.2) THEN
        PRE = HWULDO(PM(1,1),PCM(1,4))*HWULDO(PM(1,3),PCM(1,6))
        PRE = HALF/SQRT(PRE)
        DO 20 P0=1,2
        DO 20 P2=1,2
        ME(P0,2,P2) = (0.0D0,0.0D0)
 20     ME(P0,1,P2) = PRE*( A(O(P2))*S(6,3,O(P2))*F0(  P2 ,O(P0),3)
     &                     +A(  P2 )*MA(3)*       F0(O(P2),O(P0),6))
C--fermion --> scalar antifermion
      ELSEIF(I2DRTP(IMODE).EQ.3) THEN
        PRE = HWULDO(PM(1,1),PCM(1,4))*HWULDO(PM(1,3),PCM(1,6))
        PRE =-HALF/SQRT(PRE)
        DO 30 P0=1,2
        DO 30 P2=1,2
        ME(P0,2,P2) = (0.0D0,0.0D0)
 30     ME(P0,1,P2) = PRE*( A(  P0 )*S(4,1,P0)*F2M(O(P0),O(P2),1)
     &                     -A(O(P0))*MA(1)    *F2M(  P0 ,O(P2),4))
C--fermion --> fermion gauge boson
      ELSEIF(I2DRTP(IMODE).EQ.4) THEN
        PRE = HWULDO(PM(1,1),PCM(1,4))*HWULDO(PM(1,2),PCM(1,5))*
     &        HWULDO(PM(1,3),PCM(1,6))
        PRE = HALF/SQRT(PRE)
        DO 40 P0=1,2
        DO 40 P1=1,2
        ME(P0,P1,1) =-PRE*A(1)*F1F(O(P1),2,3)*S(3,6,2)*F0(1,O(P0),3)
 40     ME(P0,P1,2) = PRE*     F1F(O(P1),1,3)*S(3,6,1)*F0(2,O(P0),3)
C--scalar  --> fermion antifermion
      ELSEIF(I2DRTP(IMODE).EQ.5) THEN
        PRE = HWULDO(PM(1,2),PCM(1,5))*HWULDO(PM(1,3),PCM(1,6))
        PRE =-HALF/SQRT(PRE)
        DO 50 P1=1,2
        DO 50 P2=1,2
        ME(2,P1,P2) = (0.0D0,0.0D0)
 50     ME(1,P1,P2) = PRE*( A(O(P1))*S(5,2,O(P1))*F2M(  P1 ,O(P2),2)
     &                     +A(  P1 )*MA(2)*       F2M(O(P1),O(P2),5))
C--scalar --> fermion fermion
      ELSEIF(I2DRTP(IMODE).EQ.6) THEN
        PRE = HWULDO(PM(1,2),PCM(1,5))*HWULDO(PM(1,3),PCM(1,6))
        PRE = HALF/SQRT(PRE)
        DO 60 P1=1,2
        DO 60 P2=1,2
        ME(2,P1,P2) = (0.0D0,0.0D0)
 60     ME(1,P1,P2) = PRE*( A(O(P2))*S(6,3,O(P2))*F1M(  P2 ,P1,3)
     &                     +A(  P2 )*MA(3)*       F1M(O(P2),P1,6))
C--fermion --> fermion pion
      ELSEIF(I2DRTP(IMODE).EQ.7) THEN
        PRE = HWULDO(PM(1,1),PCM(1,4))*HWULDO(PM(1,2),PCM(1,5))
        PRE = 0.25D0/SQRT(PRE)/RMASS(198)**2
        DO 70 P0=1,2
        DO 70 P1=1,2
        ME(P0,P1,2) = (0.0D0,0.0D0)
 70     ME(P0,P1,1) =PRE*(
     &              MA(1)*A(O(P0))*( S(5,2,O(P1))*F2(  P1 ,O(P0),2,4)
     &                                     +MA(2)*F2(O(P1),O(P0),5,4))
     &            +A(P0)*S(1,4,P0)*( S(5,2,O(P1))*F2(  P1 ,  P0 ,2,1)
     &                                     +MA(2)*F2(O(P1),  P0 ,5,1)))
C--scalar  --> antifermion fermion
      ELSEIF(I2DRTP(IMODE).EQ.8) THEN
        PRE = HWULDO(PM(1,2),PCM(1,5))*HWULDO(PM(1,3),PCM(1,6))
        PRE =-HALF/SQRT(PRE)
        DO 80 P1=1,2
        DO 80 P2=1,2
        ME(2,P1,P2) = (0.0D0,0.0D0)
 80     ME(1,P1,P2) = PRE*( A(O(P2))*S(6,3,O(P2))*F1M(  P2 ,O(P1),3)
     &                     +A(  P2 )*MA(3)*       F1M(O(P2),O(P1),6))
C--neutralino --> gravitino photon
      ELSEIF(I2DRTP(IMODE).EQ.9) THEN
        PRE = TWO*HWULDO(PM(1,1),PCM(1,4))*HWULDO(PM(1,3),PCM(1,6))
        PRE = TWO/SQRT(PRE)
        DO 90 P1=1,2
        DO 90 P2=1,2
        ME(P1,P2,O(P2)) = (0.0D0,0.0D0)
 90     ME(P1,P2,  P2 ) = PRE*S(2,3,P2)*S(3,6,O(P2))*
     &       S(3,2,P2)*F0(O(P2),P1,2)
C--neutralino --> gravitino scalar
      ELSEIF(I2DRTP(IMODE).EQ.10) THEN
        PRE = TWO*HWULDO(PM(1,1),PCM(1,4))
        PRE = ONE/SQRT(PRE)
        DO 100 P1=1,2
        DO 100 P2=1,2
        ME(P1,P2,2) = (0.0D0,0.0D0)
 100    ME(P1,P2,1) = PRE*F2(P2,1,2,2)*F0(1,O(P1),2)
C--sfermion --> fermion gravitino
      ELSEIF(I2DRTP(IMODE).EQ.11) THEN
        PRE = TWO*HWULDO(PM(1,2),PCM(1,5))
        PRE = ONE/SQRT(PRE)
        DO 110 P1=1,2
        DO 110 P2=1,2
        ME(2,P1,P2) = (0.0D0,0.0D0)
 110    ME(1,P1,P2) = PRE*A(O(P2))*F1M(O(P1),P2,3)*F0B(P2,P2,3,3) 
C--antisfermion --> antifermion gravitino
      ELSEIF(I2DRTP(IMODE).EQ.12) THEN
        PRE = TWO*HWULDO(PM(1,2),PCM(1,5))
        PRE = ONE/SQRT(PRE)
        DO 120 P1=1,2
        DO 120 P2=1,2
        ME(2,P1,P2) = (0.0D0,0.0D0)
 120    ME(1,P1,P2) = PRE*A(O(P2))*F0B(P2,P2,3,3)*F1(P2,O(P1),3)
C--scalar --> antifermion antifermion
      ELSEIF(I2DRTP(IMODE).EQ.13) THEN
        PRE = HWULDO(PM(1,2),PCM(1,5))*HWULDO(PM(1,3),PCM(1,6))
        PRE = HALF/SQRT(PRE)
        DO 130 P1=1,2
        DO 130 P2=1,2
        ME(2,P1,P2) = (0.0D0,0.0D0)
 130    ME(1,P1,P2) = PRE*( A(  P1 )*S(5,2,  P1 )*F2M(O(P1),O(P2),2)
     &                     +A(O(P1))*MA(2)       *F2M(  P1 ,O(P2),5))
C--antifermion --> scalar antifermion
      ELSEIF(I2DRTP(IMODE).EQ.14) THEN
        PRE = HWULDO(PM(1,1),PCM(1,4))*HWULDO(PM(1,3),PCM(1,6))
        PRE = HALF/SQRT(PRE)
        DO 140 P0=1,2
        DO 140 P2=1,2
        ME(P0,2,P2) = (0.0D0,0.0D0)
 140    ME(P0,1,P2) = PRE*( A(O(P0))*S(4,1,O(P0))*F2M(  P0 ,O(P2),1)
     &                     -A(  P0 )*MA(1)       *F2M(O(P0),O(P2),4))
C--unrecognized type of diagram
      ELSE
        CALL HWWARN('HWDSM2',500,*999)
      ENDIF
C--now compute the weight
      WGT = ZERO
      DO 500 P0 =1,2
      DO 500 P0P=1,2
      DO 500 P1 =1,2
      DO 500 P2 =1,2
 500  WGT = WGT+PHS*P2MODE(IMODE)*ME(P0,P1,P2)*DCONJG(ME(P0P,P1,P2))*
     &          RHOIN(P0,P0P)
      IF(I2DRTP(IMODE).EQ.5.OR.I2DRTP(IMODE).EQ.6.OR.
     &   I2DRTP(IMODE).EQ.8.OR.I2DRTP(IMODE).EQ.13) GOTO 300
C--issue warning if greater than maximum
      IF(WGT.GT.WTMAX) THEN
        CALL HWWARN('HWDSM2',1,*200)
C        WRITE(LHWOUT,2000) RNAME(IDK(ID2PRT(IMODE))),
        WRITE(6,2000) RNAME(IDK(ID2PRT(IMODE))),
     &   RNAME(IDKPRD(1,ID2PRT(IMODE))),RNAME(IDKPRD(2,ID2PRT(IMODE))),
     &   WTMAX,1.1D0*WGT
        WT2MAX(IMODE) = 1.1D0*WGT
        WTMAX         = WT2MAX(IMODE)
      ENDIF
 200  IF(HWRGEN(0)*WTMAX.GT.WGT.AND.NTRY.LT.NSNTRY) GOTO 1000
      IF(NTRY.GE.NSNTRY) CALL HWWARN('HWDSM2',100,*999)
C--now enter the momenta in the common block
 300  CALL HWVEQU(5,P(1,2),PHEP(1,IOUT1))
      CALL HWVEQU(5,P(1,3),PHEP(1,IOUT2))
C--set up the spin information
C--setup for all decays
      JMOSPN(NSPN+1) = IDSPIN
      JMOSPN(NSPN+2) = IDSPIN
      JDASPN(1,IDSPIN) = NSPN+1
      JDASPN(2,IDSPIN) = NSPN+2
      IDSPN(NSPN+1) = IOUT1
      IDSPN(NSPN+2) = IOUT2
      DO 11 I=1,2
      DECSPN(NSPN+I) = .FALSE.
      DO 11 J=1,2
 11   JDASPN(I,NSPN+J) = 0
      ISNHEP(IOUT1) = NSPN+1
      ISNHEP(IOUT2) = NSPN+2
      DO 12 I=1,2
        IF(RSPIN(IDHW(IDSPN(NSPN+I))).EQ.ZERO) THEN
          RHOSPN(1,1,NSPN+I) = ONE
          RHOSPN(2,1,NSPN+I) = ZERO
          RHOSPN(1,2,NSPN+I) = ZERO
          RHOSPN(2,2,NSPN+I) = ZERO
        ELSE
          RHOSPN(1,1,NSPN+I) = HALF
          RHOSPN(2,1,NSPN+I) = ZERO
          RHOSPN(1,2,NSPN+I) = ZERO
          RHOSPN(2,2,NSPN+I) = HALF
        ENDIF
 12   CONTINUE
      NSPN = NSPN+2
C--now enter the matrix element
      DO 150 P0=1,2
      DO 150 P1=1,2
      DO 150 P2=1,2
      MESPN(P0,P1,P2,2,1,IDSPIN) = (0.0D0,0.0D0)
 150  MESPN(P0,P1,P2,1,1,IDSPIN) = ME(P0,P1,P2)
      SPNCFC(1,1,IDSPIN) = ONE
      NCFL(IDSPIN) = 1
      RETURN
C--format statements
 2000 FORMAT(/'WEIGHT FOR DECAY ',A8,' --> ',A8,' ',A8, 'EXCEEDS MAX',
     &       /10X,'    MAXIMUM WEIGHT =',1PG24.16,
     &       /10X,'NEW MAXIMUM WEIGHT =',1PG24.16)
 999  END
CDECK  ID>, HWDSM3.
*CMZ :-        -09/04/02  13:46:07  by  Peter Richardson
*-- Author :    Peter Richardson
C-----------------------------------------------------------------------
      SUBROUTINE HWDSM3(NPR,ID,IOUT1,IOUT2,IOUT3,IMODE,RHOIN,IDSPIN)
C-----------------------------------------------------------------------
C     Master subroutine for three body SUSY and spin ME's
C     Uses HWD3ME to generate the momenta etc
C-----------------------------------------------------------------------
#include "HERWIG65.INC"
      DOUBLE COMPLEX F0(2,2,8),F1(2,2,8),F1M(2,2,8),F3(2,2,8),
     &     F0M(2,2,8),F2(2,2,8),RHOIN(2,2),F01(2,2,8,8)
      DOUBLE PRECISION A,B,MS,MWD,MA,MB,MA2,MB2,M342,M232,M242,MR,
     &  P(5,4),PZ(5),HWRGEN,CV,CA,BR,PM(5,4),CFTHRE(NCFMAX,NCFMAX)
      INTEGER ID,IDP(4+NDIAGR),NPR,ITYPE,I,IB,ID1,ID2,IDSPIN,
     &     DRTYPE(NDIAGR),IOUT(3),IMODE,IOUT1,IOUT2,IOUT3,J,NCTHRE,
     &     DRCF(NDIAGR)
      COMMON/HWD3BY/F0,F0M,F1M,F1,F2,F3,F01,A(2,NDIAGR),B(2,NDIAGR),
     &     MS(NDIAGR),MWD(NDIAGR),MR(NDIAGR),MA(4),MA2(4),MB(4),MB2(4),
     &     M342,M232,M242,P,PM,CFTHRE,IDP,DRTYPE,NCTHRE,DRCF
      EXTERNAL HWRGEN
      SAVE PZ,IOUT,ITYPE,ID1,ID2
C--calculate the matrix element for a three body decay
      IF(NPR.EQ.3) THEN
C--set up the decay products, if a SUSY decay the SUSY particle
C--must be the first decay product
        IF(ABS(IDHEP(IOUT1)).GT.1000000) THEN
          IOUT(1) = IOUT1
          IOUT(2) = IOUT2
          IOUT(3) = IOUT3
        ELSEIF(ABS(IDHEP(IOUT2)).GT.1000000) THEN
          IOUT(1) = IOUT2
          IOUT(2) = IOUT1
          IOUT(3) = IOUT3
        ELSEIF(ABS(IDHEP(IOUT3)).GT.1000000) THEN
          IOUT(1) = IOUT3
          IOUT(2) = IOUT1
          IOUT(3) = IOUT3
C--special for top decay (bottom must be first)
        ELSEIF(ABS(IDHEP(ID)).EQ.6) THEN
          IOUT(1) = IOUT3
          IOUT(2) = IOUT1
          IOUT(3) = IOUT2
        ELSE
          IOUT(1) = IOUT2
          IOUT(2) = IOUT1
          IOUT(3) = IOUT3
        ENDIF
C--fermion must be second and antifermion third
        IF(IDHEP(IOUT(2)).LT.0.AND.
     &    (ABS(IDHEP(IOUT(1))).GT.1000000.OR.ABS(IDHEP(ID)).EQ.6)) THEN
          I = IOUT(2)
          IOUT(2) = IOUT(3)
          IOUT(3) = I
        ENDIF
C--setup the OFF SHELL MASSES
        MA(1) = PHEP(5,ID)
        DO 1 I=1,3
 1      MA(I+1) = PHEP(5,IOUT(I))
        DO 2 I=1,4
 2      MA2(I) = MA(I)**2
C--call to ME code
        CALL HWD3ME(ID,0,IMODE,RHOIN,IDSPIN)
        IF(IERROR.NE.0) RETURN
C--juggle the momenta for the RPV BV gluino if needed
        IF(SPCOPT.EQ.2.AND.N3NCFL(IMODE).EQ.3) THEN
          IF(NCFL(IDSPIN).EQ.2) THEN
            IOUT(1) = IOUT1 
            IOUT(2) = IOUT2 
            IOUT(3) = IOUT3 
          ELSEIF(NCFL(IDSPIN).EQ.3) THEN
            IOUT(1) = IOUT3
            IOUT(2) = IOUT2
            IOUT(3) = IOUT1
          ENDIF
          DO I=1,3
            IDHW(IOUT(I)) = IDP(I+1)
          ENDDO
        ENDIF
C--copy momenta into event record
        DO 3 I=1,3
 3      CALL HWVEQU(5,P(1,1+I),PHEP(1,IOUT(I)))
C--enter the spin information in the common block
        IF(SYSPIN) THEN
C--set up if start of new spin chain
          IF(NSPN.EQ.0) THEN
C--zero the elements
            CALL HWVZRI(  NMXHEP,ISNHEP)
            CALL HWVZRI(  NMXSPN,JMOSPN)
            CALL HWVZRI(2*NMXSPN,JDASPN)
            CALL HWVZRI(  NMXSPN, IDSPN)
            NSPN = NSPN+1
            JMOSPN(NSPN) = 0
            IDSPN (NSPN) = ID
            DECSPN(NSPN) = .FALSE.
C--set up spin density matrix for particle
            IF(RSPIN(IDHW(ID)).EQ.ZERO) THEN
              RHOSPN(1,1,NSPN) = ONE
              RHOSPN(2,1,NSPN) = ZERO
              RHOSPN(1,2,NSPN) = ZERO
              RHOSPN(2,2,NSPN) = ZERO
            ELSE
              RHOSPN(1,1,NSPN) = HALF
              RHOSPN(2,1,NSPN) = ZERO
              RHOSPN(1,2,NSPN) = ZERO
              RHOSPN(2,2,NSPN) = HALF
            ENDIF
            ISNHEP(ID)    = NSPN
          ENDIF
C--enter the decay products
          JDASPN(1,IDSPIN) = NSPN+1
          JDASPN(2,IDSPIN) = NSPN+3
          DO 7 I=1,3
          JMOSPN(NSPN+I  ) = IDSPIN
          IDSPN (NSPN+I  ) = IOUT(I)
          DECSPN(NSPN+I  ) = .FALSE.
          ISNHEP(IOUT(I) ) = NSPN+I
          IF(RSPIN(IDHW(IOUT(I))).EQ.ZERO) THEN
            RHOSPN(1,1,NSPN+I) = ONE
            RHOSPN(2,1,NSPN+I) = ZERO
            RHOSPN(1,2,NSPN+I) = ZERO
            RHOSPN(2,2,NSPN+I) = ZERO
          ELSE
            RHOSPN(1,1,NSPN+I) = HALF
            RHOSPN(2,1,NSPN+I) = ZERO
            RHOSPN(1,2,NSPN+I) = ZERO
            RHOSPN(2,2,NSPN+I) = HALF
          ENDIF
          DO 7 J=1,2
 7        JDASPN(J,NSPN+I) = 0
          NSPN = NSPN+3
        ENDIF
C--select the decay mode and generate the decay for a two body mode
      ELSEIF(NPR.EQ.2) THEN
        IF(IDHW(IOUT2).GE.198.AND.IDHW(IOUT2).LE.200) THEN
          IB = IDHW(IOUT2)
          IOUT(1) = IOUT1
          IOUT(2) = IOUT2
        ELSEIF(IDHW(IOUT1).GE.198.AND.IDHW(IOUT1).LE.200) THEN
          IB = IDHW(IOUT1)
          IOUT(1) = IOUT2
          IOUT(2) = IOUT1
        ELSE
          CALL HWWARN('HWDSM3',501,*999)
        ENDIF
C--setup the off shell masses and particle ids for me code
        MA(1) = PHEP(5,ID)
        MA(2) = PHEP(5,IOUT(1))
        CALL HWDBOZ(IB,ID1,ID2,CV,CA,BR,0)
        ITYPE = ID1
        IF(IB.EQ.199) ITYPE = ITYPE+1
        IF(ITYPE.GT.120) ITYPE = ITYPE-114
        IF(IB.NE.200) ITYPE = ITYPE/2
C--generate momenta of decay products
        CALL HWD3ME(ID,ITYPE,IMODE,RHOIN,IDSPIN)
        CALL HWVEQU(5,P(1,2),PHEP(1,IOUT(1)))
        CALL HWVSUM(4,P(1,3),P(1,4),PZ)
        CALL HWUMAS(PZ)
        CALL HWVEQU(5,PZ,PHEP(1,IOUT(2)))
C--enter the spin information in the common block if starting new chain
        IF(SYSPIN.AND.NSPN.EQ.0) THEN
C--zero elements of common block
          CALL HWVZRI(  NMXHEP,ISNHEP)
          CALL HWVZRI(  NMXSPN,JMOSPN)
          CALL HWVZRI(2*NMXSPN,JDASPN)
          CALL HWVZRI(  NMXSPN, IDSPN)
          NSPN = NSPN+1
          JMOSPN(NSPN) = 0
          IDSPN (NSPN) = ID
          DECSPN(NSPN) = .FALSE.
          IF(RSPIN(IDHW(ID)).EQ.ZERO) THEN
            RHOSPN(1,1,NSPN) = ONE
            RHOSPN(2,1,NSPN) = ZERO
            RHOSPN(1,2,NSPN) = ZERO
            RHOSPN(2,2,NSPN) = ZERO
          ELSE
            RHOSPN(1,1,NSPN) = HALF
            RHOSPN(2,1,NSPN) = ZERO
            RHOSPN(1,2,NSPN) = ZERO
            RHOSPN(2,2,NSPN) = HALF
          ENDIF
          ISNHEP(ID)    = NSPN
        ENDIF
        IF(SYSPIN) THEN
          IDSPN (NSPN+1  ) = IOUT(1)
          ISNHEP(IOUT(1))  = NSPN+1
        ENDIF
C--put the boson decay products into the event record for a two body mode
      ELSEIF(NPR.EQ.-1) THEN
        IOUT(1) = JDAHEP(1,IOUT(2))
        IOUT(2) = NHEP+1
        IOUT(3) = NHEP+2
C--set up the status of the particles
        ISTHEP(IOUT(1)) = 195
        JDAHEP(1,IOUT(1)) = NHEP+1
        JDAHEP(2,IOUT(1)) = NHEP+2
C--find the ID's of the particles
        IF(IDHW(IOUT(1)).EQ.200) THEN
          ID1 = ITYPE
          IF(ITYPE.GT.6) ID1 = ID1+114
          ID2 = ID1+6
        ELSE
          ID1 = 2*ITYPE-1
          IF(ITYPE.GT.3) ID1 = ID1+114
          ID2 = ID1+7
          IF(IDHW(IOUT(1)).EQ.198) THEN
            I   = ID1+6
            ID1 = ID2-6
            ID2 = I
          ENDIF
        ENDIF
C--put id's of decay products into the event record
        IDHW(NHEP+1)  = ID1
        IDHW(NHEP+2)  = ID2
        IDHEP(NHEP+1) = IDPDG(ID1)
        IDHEP(NHEP+2) = IDPDG(ID2)
C--boost decay products momenta to rest frame of boson
        CALL HWULOF(PZ,P(1,3),P(1,3))
        CALL HWULOF(PZ,P(1,4),P(1,4))
C--boost back to lab using new boson
        CALL HWULOB(PHEP(1,IOUT(1)),P(1,3),PHEP(1,NHEP+1))
        CALL HWULOB(PHEP(1,IOUT(1)),P(1,4),PHEP(1,NHEP+2))
C--setup for decay to quarks
        IF(ID1.LE.12) THEN
          ISTHEP(NHEP+1) = 113
          ISTHEP(NHEP+2) = 114
          JMOHEP(2,NHEP+1) = NHEP+2
          JDAHEP(2,NHEP+1) = NHEP+2
          JMOHEP(2,NHEP+2) = NHEP+1
          JDAHEP(2,NHEP+2) = NHEP+1
          JMOHEP(1,NHEP+1) = IOUT(1)
          JMOHEP(1,NHEP+2) = IOUT(1)
C--setup for decay to leptons
        ELSE
          ISTHEP(NHEP+1) = 193
          ISTHEP(NHEP+2) = 193
          JMOHEP(1,NHEP+1) = IOUT(1)
          JMOHEP(1,NHEP+2) = IOUT(1)
          JMOHEP(2,NHEP+1) = JMOHEP(1,IOUT(1))
          JMOHEP(2,NHEP+2) = JMOHEP(1,IOUT(1))
          JDAHEP(1,NHEP+1) = 0
          JDAHEP(1,NHEP+2) = 0
          JDAHEP(2,NHEP+1) = 0
          JDAHEP(2,NHEP+2) = 0
        ENDIF
        NHEP=NHEP+2
C--finish entering the spin information in the common block
        IF(SYSPIN) THEN
          JDASPN(1,IDSPIN) = NSPN+1
          JDASPN(2,IDSPIN) = NSPN+3
          DO 6 I=1,3
          JMOSPN(NSPN+I  ) = IDSPIN
          DECSPN(NSPN+I  ) = .FALSE.
          IF(RSPIN(IDHW(IOUT(I))).EQ.ZERO) THEN
            RHOSPN(1,1,NSPN+I) = ONE
            RHOSPN(2,1,NSPN+I) = ZERO
            RHOSPN(1,2,NSPN+I) = ZERO
            RHOSPN(2,2,NSPN+I) = ZERO
          ELSE
            RHOSPN(1,1,NSPN+I) = HALF
            RHOSPN(2,1,NSPN+I) = ZERO
            RHOSPN(1,2,NSPN+I) = ZERO
            RHOSPN(2,2,NSPN+I) = HALF
          ENDIF
          DO 6 J=1,2
 6        JDASPN(J,NSPN+I) =0
          NSPN = NSPN+3
          IDSPN (NSPN-1) = NHEP-1
          IDSPN (NSPN  ) = NHEP
          ISNHEP(NHEP-1) = NSPN-1
          ISNHEP(NHEP  ) = NSPN
        ENDIF
C--perform the parton shower for the decay products of the gauge boson
        IF(ID1.LE.12) CALL HWBGEN
C--error issue warning
      ELSE
        CALL HWWARN('HWDSM3',500,*999)
      ENDIF
 999  END
CDECK  ID>, HWDSM4.
*CMZ :-        -11/10/01  14:03:42  by  Peter Richardson
*-- Author :    Peter Richardson
C-----------------------------------------------------------------------
      SUBROUTINE HWDSM4(IOPT,ID,IOUT1,IOUT2,IMODE)
C-----------------------------------------------------------------------
C     Subroutine to perform the four body decays
C     IOPT = 1 select decay mode and generate momenta
C     IOPT = 2 enter first decays and perform parton shower
C-----------------------------------------------------------------------
#include "HERWIG65.INC"
      INTEGER IOPT,ID,IOUT1,IOUT2,IB(2),I,IDF(4),ITYPE(2),IMODE,
     &     IDP(4+NDIAGR),ID1,ID2,J
      DOUBLE PRECISION CV,CA,A,B,MS,MWD,MR,M,M2,P(5,5),PW(5,2),BR
      COMMON/HWD4BY/A(2),B(2),MS(2),MWD(2),MR(2),M(5),M2(5),P,IDP
      SAVE PW,ITYPE
C--generate the decay
      IF(IOPT.EQ.1) THEN
        IB(1) = IDHW(IOUT1)
        IB(2) = IDHW(IOUT2)
C--select the decays of the bosons
        DO 1 I=1,2
        CALL HWDBOZ(IB(I),IDF(2*I-1),IDF(2*I),CV,CA,BR,1)
        ITYPE(I) = IDF(2*I-1)
        IF(IB(I).EQ.199) ITYPE(I)    = ITYPE(I)+1
        IF(ITYPE(I).GT.120) ITYPE(I) = ITYPE(I)-114
 1      IF(IB(I).NE.200) ITYPE(I)    = ITYPE(I)/2
C--generate the momenta of the decay products
        CALL HWD4ME(ID,ITYPE(1),ITYPE(2),IMODE)
        DO 2 I=1,2
        CALL HWVSUM(4,P(1,2*I),P(1,2*I+1),PW(1,I))
 2      CALL HWUMAS(PW(1,I))
        CALL HWVEQU(5,PW(1,1),PHEP(1,IOUT1))
        CALL HWVEQU(5,PW(1,2),PHEP(1,IOUT2))
        IF(SYSPIN) THEN
          IDSPN(1)     = JDAHEP(1,ID)
          DECSPN(1)    = .FALSE.
          ISNHEP(JDAHEP(1,ID)) = 1
          JDASPN(1,1)  = 2
          JDASPN(2,1)  = 5
          DO 4 I=2,5
          DECSPN(I) = .FALSE.
 4        JMOSPN(I) = 1
        ENDIF
      ELSEIF(IOPT.EQ.2) THEN
        IB(1) = JDAHEP(1,IOUT1)
        IB(2) = JDAHEP(1,IOUT2)
        DO 3 I=1,2
          ISTHEP(IB(I)) = 195
          JDAHEP(1,IB(I)) = NHEP+1
          JDAHEP(2,IB(I)) = NHEP+2
C--find the ID's of the particles
          IF(IDHW(IB(I)).EQ.200) THEN
            ID1 = ITYPE(I)
            IF(ITYPE(I).GT.6) ID1 = ID1+114
            ID2 = ID1+6
          ELSE
            ID1 = 2*ITYPE(I)-1
            IF(ITYPE(I).GT.3) ID1 = ID1+114
            ID2 = ID1+7
            IF(IDHW(IB(I)).EQ.198) THEN
              J   = ID1+6
              ID1 = ID2-6
              ID2 = J
            ENDIF
          ENDIF
C--put id's of decay products into the event record
          IDHW(NHEP+1)  = ID1
          IDHW(NHEP+2)  = ID2
          IDHEP(NHEP+1) = IDPDG(ID1)
          IDHEP(NHEP+2) = IDPDG(ID2)
C--boost decay products momenta to rest frame of boson
          CALL HWULOF(PW(1,I),P(1,2*I  ),P(1,2*I  ))
          CALL HWULOF(PW(1,I),P(1,2*I+1),P(1,2*I+1))
C--boost back to lab using new boson
          CALL HWULOB(PHEP(1,IB(I)),P(1,2*I  ),PHEP(1,NHEP+1))
          CALL HWULOB(PHEP(1,IB(I)),P(1,2*I+1),PHEP(1,NHEP+2))
C--setup for decay to quarks
          IF(ID1.LE.12) THEN
            ISTHEP(NHEP+1) = 113
            ISTHEP(NHEP+2) = 114
            JMOHEP(2,NHEP+1) = NHEP+2
            JDAHEP(2,NHEP+1) = NHEP+2
            JMOHEP(2,NHEP+2) = NHEP+1
            JDAHEP(2,NHEP+2) = NHEP+1
            JMOHEP(1,NHEP+1) = IB(I)
            JMOHEP(1,NHEP+2) = IB(I)
C--setup for decay to leptons
          ELSE
            ISTHEP(NHEP+1) = 193
            ISTHEP(NHEP+2) = 193
            JMOHEP(1,NHEP+1) = IB(I)
            JMOHEP(1,NHEP+2) = IB(I)
            JMOHEP(2,NHEP+1) = JMOHEP(1,IB(I))
            JMOHEP(2,NHEP+2) = JMOHEP(1,IB(I))
          ENDIF
C--enter the information in the spin common block
          IF(SYSPIN) THEN
            IDSPN(2*I  ) = NHEP+1
            IDSPN(2*I+1) = NHEP+2
            ISNHEP(NHEP+1) = 2*I
            ISNHEP(NHEP+2) = 2*I+1
          ENDIF
          NHEP=NHEP+2
C--perform the parton shower for the decay products of the gauge boson
          IF(ID1.LE.12) CALL HWBGEN
 3      CONTINUE
      ENDIF
 999  END
CDECK  ID>, HWDTAU.
*CMZ :-        -17/10/01  09:42:21  by  Peter Richardson
*-- Author :    Peter Richardson
C Modified by Borut Kersevan to handle external processes and neutral decays
C-----------------------------------------------------------------------
      SUBROUTINE HWDTAU(IOPT,IHEP,POL)
C-----------------------------------------------------------------------
C     HERWIG-TAUOLA interface to perform tau decays using TAUOLA rather
C     than HERWIG
C     IOPT = 0 initialises
C     IOPT = 1 performs decay
C     IOPT = 2 write outs final TAUOLA information
C-----------------------------------------------------------------------
#include "HERWIG65.INC"
      INTEGER IOPT,IHEP,ID,ITAU,I,IMO,NHEPPO
      DOUBLE PRECISION POL,PLAB(4)
      REAL POL1(4)
      CHARACTER *8 DUMMY
      DATA PLAB/0.0D0,0.0D0,0.0D0,1.0D0/
C--common block for PHOTOS
      LOGICAL QEDRAD
      COMMON /PHOQED/ QEDRAD(NMXHEP)
C--common blocks for TAUOLA
      INTEGER NP1,NP2
      COMMON /TAUPOS/ NP1, NP2
      DOUBLE PRECISION Q1(4),Q2(4),P1(4),P2(4),P3(4),P4(4)
      COMMON / MOMDEC / Q1,Q2,P1,P2,P3,P4

C(BPK)-->
      INTEGER KEYPOL
      COMMON /POLVAL/KEYPOL
      SAVE /POLVAL/

C TAUOLA PARAMETER COMMON
      INTEGER IFPOL,IFRADCOR,IFDMOD,IFJAK1,IFJAK2
      REAL FFPHX
      COMMON /FFTAU/IFPOL,IFRADCOR,IFDMOD,IFJAK1,IFJAK2,FFPHX

      IFPHOT=0
     
C(BPK)--<
C--initialisation
      IF(IOPT.EQ.-1) THEN
C(BPK)     DONE ELSEWHERE
      ELSEIF(IOPT.EQ.1) THEN
C(BPK)--> 
         IF (POL.EQ.0.AND.KEYPOL.EQ.1) THEN ! ASSUME WE CAN DO IT BETTER
            ISTHEP(IHEP)=1 ! MAKE IT STABLE
            RETURN              
         ENDIF

         IF (IFDMOD.EQ.1) THEN  ! LEP-LEP
            CALL TAUBRS(1)
         ELSEIF (IFDMOD.EQ.2) THEN ! HAD-HAD
            CALL TAUBRS(2)
         ELSEIF (IFDMOD.EQ.3) THEN ! LEP-HAD
            PRINT *,
     # 'HWDTAU: LEP-HAD MODE NOT POSSIBLE IN INTERNAL HERWIG PROCESSES'
            STOP
         ELSEIF (IFDMOD.EQ.4) THEN ! TAU->PI NU
            CALL TAUBRS(4)
         ELSE                   ! INCLUSIVE DECAY
            CALL TAUBRS(0)
         ENDIF
C(BPK)--<
      
        ISTHEP(IHEP)=195
        ID = IDHW(IHEP)
        IMO = IHEP
 1      IMO = JMOHEP(1,IMO)
        IF(IDHW(IMO).EQ.ID) GOTO 1
C--id of tau for tauola
        IF(ID.EQ.125) THEN ! tau-
          ITAU = 2
          NP1 = IHEP
          NP2 = IHEP
          DO I=1,4
             P2(I) = PHEP(I,IHEP)
             Q1(I) = PLAB(I)
             P1(I)=Q1(I)-P2(I)
          ENDDO
        ELSEIF(ID.EQ.131) THEN ! tau+ 
          ITAU = 1
          NP1 = IHEP
          NP2 = IHEP
          DO I=1,4
             P1(I) = PHEP(I,IHEP)
             Q1(I) = PLAB(I)
             P2(I)=Q1(I)-P1(I)
          ENDDO
        ELSE
          CALL HWWARN('HWDTAU',501,*999)
        ENDIF
C--set up the tau polarization
        POL1(1) = 0.
        POL1(2) = 0.
        POL1(3) = REAL(POL)
        POL1(4) = 0.

C(BPK)--> 
C--tau momentum
c        DO I=1,4
c          P1(I) = PHEP(I,IHEP)
c          P2(I) = PHEP(I,IHEP)
C--we measure tau spins in lab frame
c          Q1(I) = PLAB(I)
c        ENDDO
C(BPK)--<

C--perform the decay and generate QED radiation if needed
        NHEPPO=NHEP
        CALL DEXAY(ITAU,POL1)
        IF(IFPHOT.EQ.1) THEN
          IF(ID.EQ.1) THEN
            CALL PHOTOS(NP1)
          ELSE
            CALL PHOTOS(NP2)
          ENDIF
        ENDIF
        IF(NHEPPO.NE.NHEP) THEN
           CALL HWHEPC(NHEPPO)
c          DO 2 I=NHEPPO+1,NHEP
c          CALL HWVEQU(4,VHEP(1,IHEP),VHEP(1,I))
c 2        CALL HWUIDT(1,IDHEP(I),IDHW(I),DUMMY)
        ENDIF
C--write out info at end
      ELSEIF(IOPT.EQ.2) THEN
C(BPK) DONE ELSEWHERE
C--otherwise issue warning
      ELSE
        CALL HWWARN('HWDTAU',500,*999)
      ENDIF
 999  END
CDECK  ID>, HWDTHR.
*CMZ :-        -26/04/91  14.55.44  by  Federico Carminati
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE HWDTHR(P0,P1,P2,P3,WEIGHT)
C-----------------------------------------------------------------------
C     GENERATES THREE-BODY DECAY 0->1+2+3 DISTRIBUTED
C     ACCORDING TO PHASE SPACE * WEIGHT
C-----------------------------------------------------------------------
      DOUBLE PRECISION HWRGEN,HWRUNI,A,B,C,D,AA,BB,CC,DD,EE,FF,PP,QQ,WW,
     & RR,PCM1,PC23,WEIGHT,P0(5),P1(5),P2(5),P3(5),P23(5),TWO
      EXTERNAL HWRGEN,HWRUNI,WEIGHT
      PARAMETER (TWO=2.D0)
      A=P0(5)+P1(5)
      B=P0(5)-P1(5)
      C=P2(5)+P3(5)
      IF (B.LT.C) CALL HWWARN('HWDTHR',100,*999)
      D=ABS(P2(5)-P3(5))
      AA=A*A
      BB=B*B
      CC=C*C
      DD=D*D
      EE=(B-C)*(A-D)
      A=0.5*(AA+BB)
      B=0.5*(CC+DD)
      C=4./(A-B)**2
C
C  CHOOSE MASS OF SUBSYSTEM 23 WITH PRESCRIBED DISTRIBUTION
C
   10 FF=HWRUNI(0,BB,CC)
      PP=(AA-FF)*(BB-FF)
      QQ=(CC-FF)*(DD-FF)
      WW=WEIGHT(FF,A,B,C)**2
      RR=EE*FF*HWRGEN(0)
      IF (PP*QQ*WW.LT.RR*RR) GOTO 10
C
C  FF IS MASS SQUARED OF SUBSYSTEM 23.
C
C  DO 2-BODY DECAYS 0->1+23, 23->2+3
C
      P23(5)=SQRT(FF)
      PCM1=SQRT(PP)*0.5/P0(5)
      PC23=SQRT(QQ)*0.5/P23(5)
      CALL HWDTWO(P0,P1,P23,PCM1,TWO,.TRUE.)
      CALL HWDTWO(P23,P2,P3,PC23,TWO,.TRUE.)
  999 END
CDECK  ID>, HWDTOP.
*CMZ :-        -09/12/92  11.03.46  by  Bryan Webber
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE HWDTOP(DECAY)
C-----------------------------------------------------------------------
C     DECIDES WHETHER TO DO TOP QUARK DECAY BEFORE HADRONIZATION
C-----------------------------------------------------------------------
#include "HERWIG65.INC"
      LOGICAL DECAY
      DECAY=RMASS(6).GT.130D0
      END
