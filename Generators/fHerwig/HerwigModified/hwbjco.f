CDECK  ID>, HWBJCO.
*CMZ :-        -30/09/02  09.19.58  by  Peter Richardson
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE HWBJCO
C-----------------------------------------------------------------------
C     COMBINES JETS WITH REQUIRED KINEMATICS
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      DOUBLE PRECISION HWULDO,EPS,PTX,PTY,PF,PTINF,PTCON,CN,CP,SP,PP0,
     & PM0,ET0,DET,ECM,EMJ,EMP,EMS,DMS,ES,DPF,ALF,AL(2),ET(2),PP(2),
     & PT(3),PA(5),PB(5),PC(5),PQ(5),PR(5),PS(5),RR(3,3),RS(3,3),ETC,
     & PJ(NMXJET),PM(NMXJET),PBR(5),RBR(3,3),DISP(4),PLAB(5)
      INTEGER LJET,IJ1,IST,IP,ICM,IP1,IP2,NP,IHEP,MHEP,JP,KP,LP,KHEP,
     & JHEP,NE,IJT,IEND(2),IJET(NMXJET),IPAR(NMXJET)
      LOGICAL AZCOR,JETRAD,DISPRO,DISLOW
      EXTERNAL HWULDO
      PARAMETER (EPS=1.D-4)
      IF (IERROR.NE.0) RETURN
      AZCOR=AZSOFT.OR.AZSPIN
      LJET=131
  10  IJET(1)=1
  20  IJ1=IJET(1)
      DO 40 IHEP=IJ1,NHEP
      IST=ISTHEP(IHEP)
      IF (IST.EQ.137.OR.IST.EQ.138) IST=133
      IF (IST.EQ.LJET) THEN
C---FOUND AN UNBOOSTED JET - FIND PARTNERS
        IP=JMOHEP(1,IHEP)
        ICM=JMOHEP(1,IP)
        DISPRO=IPRO/10.EQ.9.AND.IDHW(ICM).EQ.15
        DISLOW=DISPRO.AND.JDAHEP(1,ICM).EQ.JDAHEP(2,ICM)-1
        IF (IST.EQ.131) THEN
          IP1=JMOHEP(1,ICM)
          IP2=JMOHEP(2,ICM)
        ELSE
          IP1=JDAHEP(1,ICM)
          IP2=JDAHEP(2,ICM)
        ENDIF
        IF (IP1.NE.IP) CALL HWWARN('HWBJCO',100,*999)
        NP=0
        DO 30 JHEP=IP1,IP2
        NP=NP+1
        IPAR(NP)=JHEP
  30    IJET(NP)=JDAHEP(1,JHEP)
        GOTO 50
      ENDIF
  40  CONTINUE
C---NO MORE JETS?
      IF (LJET.EQ.131) THEN
        LJET=133
        GOTO 10
      ENDIF
      RETURN
  50  IF (LJET.EQ.131) THEN
C---SPACELIKE JETS: FIND SPACELIKE PARTONS
        IF (NP.NE.2) CALL HWWARN('HWBJCO',103,*999)
C---special for DIS: FIND BOOST AND ROTATION FROM LAB TO BREIT FRAME
        IF (DISPRO.AND.BREIT) THEN
          IP=2
          IF (JDAHEP(1,IP).NE.0) IP=JDAHEP(1,IP)
          CALL HWVDIF(4,PHEP(1,JMOHEP(1,ICM)),PHEP(1,JDAHEP(1,ICM)),PB)
          CALL HWUMAS(PB)
C---IF Q**2<10**-2, SOMETHING MUST HAVE ALREADY GONE WRONG
          IF (PB(5)**2.LT.1.D-2) CALL HWWARN('HWBJCO',102,*999)
          CALL HWVSCA(4,PB(5)**2/HWULDO(PHEP(1,IP),PB),PHEP(1,IP),PBR)
          CALL HWVSUM(4,PB,PBR,PBR)
          CALL HWUMAS(PBR)
          CALL HWULOF(PBR,PB,PB)
          CALL HWUROT(PB,ONE,ZERO,RBR)
        ENDIF
        PTX=0.
        PTY=0.
        PF=1.D0
        DO 90 IP=1,2
        MHEP=IJET(IP)
        IF (JDAHEP(1,MHEP).EQ.0) THEN
C---SPECIAL FOR NON-PARTON JETS
          IHEP=MHEP
          GOTO 70
        ELSE
          IST=134+IP
          DO 60 IHEP=MHEP,NHEP
  60      IF (ISTHEP(IHEP).EQ.IST) GOTO 70
C---COULDN'T FIND SPACELIKE PARTON
          CALL HWWARN('HWBJCO',101,*999)
        ENDIF
  70    CALL HWVSCA(3,PF,PHEP(1,IHEP),PS)
        IF (PTINT(3,IP).GT.ZERO) THEN
C---ADD INTRINSIC PT
          PT(1)=PTINT(1,IP)
          PT(2)=PTINT(2,IP)
          PT(3)=0.
          CALL HWUROT(PS, ONE,ZERO,RS)
          CALL HWUROB(RS,PT,PT)
          CALL HWVSUM(3,PS,PT,PS)
        ENDIF
        JP=IJET(IP)+1
        IF (AZCOR.AND.JP.LE.NHEP.AND.IDHW(JP).EQ.17) THEN
C---ALIGN CONE WITH INTERFERING PARTON
          CALL HWUROT(PS, ONE,ZERO,RS)
          CALL HWUROF(RS,PHEP(1,JP),PR)
          PTCON=PR(1)**2+PR(2)**2
          KP=JMOHEP(2,JP)
          IF (KP.EQ.0) THEN
            CALL HWWARN('HWBJCO',1,*999)
            PTINF=0.
          ELSE
            CALL HWVEQU(4,PHEP(1,KP),PB)
            IF (DISPRO.AND.BREIT) THEN
              CALL HWULOF(PBR,PB,PB)
              CALL HWUROF(RBR,PB,PB)
            ENDIF
            PTINF=PB(1)**2+PB(2)**2
            IF (PTINF.LT.EPS) THEN
C---COLLINEAR JETS: ALIGN CONES
              KP=JDAHEP(1,KP)+1
              IF (ISTHEP(KP).EQ.100.AND.(ISTHEP(KP-1)+9)/10.EQ.14) THEN
                CALL HWVEQU(4,PHEP(1,KP),PB)
                IF (DISPRO.AND.BREIT) THEN
                  CALL HWULOF(PBR,PB,PB)
                  CALL HWUROF(RBR,PB,PB)
                ENDIF
                PTINF=PB(1)**2+PB(2)**2
              ELSE
                PTINF=0.
              ENDIF
            ENDIF
          ENDIF
          IF (PTCON.NE.ZERO.AND.PTINF.NE.ZERO) THEN
            CN=1./SQRT(PTINF*PTCON)
            CP=CN*(PR(1)*PB(1)+PR(2)*PB(2))
            SP=CN*(PR(1)*PB(2)-PR(2)*PB(1))
          ELSE
            CALL HWRAZM( ONE,CP,SP)
          ENDIF
        ELSE
          CALL HWRAZM( ONE,CP,SP)
        ENDIF
C---ROTATE SO SPACELIKE IS ALONG AXIS (APART FROM INTRINSIC PT)
        CALL HWUROT(PS,CP,SP,RS)
        IHEP=IJET(IP)
        KHEP=JDAHEP(2,IHEP)
        IF (KHEP.LT.IHEP) KHEP=IHEP
        IEND(IP)=KHEP
        DO 80 JHEP=IHEP,KHEP
        CALL HWUROF(RS,PHEP(1,JHEP),PHEP(1,JHEP))
  80    CALL HWUROF(RS,VHEP(1,JHEP),VHEP(1,JHEP))
        PP(IP)=PHEP(4,IHEP)+PF*PHEP(3,IHEP)
        ET(IP)=PHEP(1,IHEP)**2+PHEP(2,IHEP)**2-PHEP(5,IHEP)**2
C---REDEFINE HARD CM
        PTX=PTX+PHEP(1,IHEP)
        PTY=PTY+PHEP(2,IHEP)
  90    PF=-PF
        PHEP(1,ICM)=PTX
        PHEP(2,ICM)=PTY
C---special for DIS: keep lepton momenta fixed
        IF (DISPRO) THEN
          IP1=JMOHEP(1,ICM)
          IP2=JDAHEP(1,ICM)
          IJT=IJET(1)
C---IJT will be used to store lepton momentum transfer
          CALL HWVDIF(4,PHEP(1,IP1),PHEP(1,IP2),PHEP(1,IJT))
          CALL HWUMAS(PHEP(1,IJT))
          IF (IDHEP(IP1).EQ.IDHEP(IP2)) THEN
            IDHW(IJT)=200
          ELSEIF (IDHEP(IP1).LT.IDHEP(IP2)) THEN
            IDHW(IJT)=199
          ELSE
            IDHW(IJT)=198
          ENDIF
          IDHEP(IJT)=IDPDG(IDHW(IJT))
          ISTHEP(IJT)=3
C---calculate boost for struck parton
C   PC is momentum of outgoing parton(s)
          IP2=JDAHEP(2,ICM)
          IF (.NOT.DISLOW) THEN
C---FOR heavy QQbar PQ and PC are old and new QQbar momenta
            CALL HWVSUM(4,PHEP(1,IP2-1),PHEP(1,IP2),PQ)
            CALL HWUMAS(PQ)
            PC(5)=PQ(5)
          ELSE
            PC(5)=PHEP(5,JDAHEP(1,IP2))
          ENDIF
          CALL HWVSUM(2,PHEP(1,IJT),PHEP(1,IJET(2)),PC)
          ET(1)=ET(2)
C---USE BREIT FRAME BOSON MOMENTUM IF NECESSARY
          IF (BREIT) THEN
            ET(2)=ET(1)+PC(5)**2+PHEP(5,IJET(2))**2
            PM0=PHEP(5,IJT)
            PP0=-PM0
          ELSE
            ET(2)=PC(1)**2+PC(2)**2+PC(5)**2
            PP0=PHEP(4,IJT)+PHEP(3,IJT)
            PM0=PHEP(4,IJT)-PHEP(3,IJT)
          ENDIF
          ET0=(PP0*PM0)+ET(1)-ET(2)
          DET=ET0**2-4.*(PP0*PM0)*ET(1)
          IF (DET.LT.ZERO) THEN
            FROST=.TRUE.
            RETURN
          ENDIF
          ALF=(SQRT(DET)-ET0)/(2.*PP0*PP(2))
          PB(1)=0.
          PB(2)=0.
          PB(5)=2.D0
          PB(3)=ALF-(1./ALF)
          PB(4)=ALF+(1./ALF)
          DO 100 IHEP=IJET(2),IEND(2)
          CALL HWULOF(PB,PHEP(1,IHEP),PHEP(1,IHEP))
          CALL HWULF4(PB,VHEP(1,IHEP),VHEP(1,IHEP))
C---BOOST FROM BREIT FRAME IF NECESSARY
          IF (BREIT) THEN
            CALL HWUROB(RBR,PHEP(1,IHEP),PHEP(1,IHEP))
            CALL HWULOB(PBR,PHEP(1,IHEP),PHEP(1,IHEP))
            CALL HWUROB(RBR,VHEP(1,IHEP),VHEP(1,IHEP))
            CALL HWULB4(PBR,VHEP(1,IHEP),VHEP(1,IHEP))
          ENDIF
  100     ISTHEP(IHEP)=ISTHEP(IHEP)+10
          CALL HWVDIF(4,VHEP(1,IPAR(2)),VHEP(1,IJET(2)),DISP)
          DO 110 IHEP=IJET(2),IEND(2)
  110     CALL HWVSUM(4,DISP,VHEP(1,IHEP),VHEP(1,IHEP))
          IF (IEND(2).GT.IJET(2)+1) ISTHEP(IJET(2)+1)=100
          CALL HWVSUM(4,PHEP(1,IJT),PHEP(1,IJET(2)),PC)
          CALL HWVSUM(4,PHEP(1,IP1),PHEP(1,IJET(2)),PHEP(1,ICM))
          CALL HWUMAS(PHEP(1,ICM))
        ELSEIF (IPRO/10.EQ.5) THEN
C Special to preserve photon momentum
           ETC=PTX**2+PTY**2+PHEP(5,ICM)**2
           ET0=ETC+ET(1)-ET(2)
           DET=ET0**2-4.*ETC*ET(1)
           IF (DET.LT.ZERO) THEN
              FROST=.TRUE.
              RETURN
           ENDIF
           ALF=(SQRT(DET)+ET0-2.*ET(1))/(2.*PP(1)*PP(2))
           PB(1)=0.
           PB(2)=0.
           PB(3)=ALF-1./ALF
           PB(4)=ALF+1./ALF
           PB(5)=2.
           IJT=IJET(2)
           DO 120 IHEP=IJT,IEND(2)
           CALL HWULOF(PB,PHEP(1,IHEP),PHEP(1,IHEP))
           CALL HWULF4(PB,VHEP(1,IHEP),VHEP(1,IHEP))
  120      ISTHEP(IHEP)=ISTHEP(IHEP)+10
           CALL HWVDIF(4,VHEP(1,IPAR(2)),VHEP(1,IJT),DISP)
           DO 130 IHEP=IJT,IEND(2)
  130      CALL HWVSUM(4,DISP,VHEP(1,IHEP),VHEP(1,IHEP))
           IF (IEND(2).GT.IJT+1) ISTHEP(IJT+1)=100
           ISTHEP(IJET(1))=ISTHEP(IJET(1))+10
           CALL HWVSUM(2,PHEP(3,IPAR(1)),PHEP(3,IJT),PHEP(3,ICM))
        ELSE
C--change to preserve either long mom or rapidity rather than long mom
C--by PR and BRW 30/9/02
          IF (PRESPL) THEN
C--PRESERVE LONG MOM OF CMF
            PHEP(4,ICM)=
     &            SQRT(PTX**2+PTY**2+PHEP(3,ICM)**2+PHEP(5,ICM)**2)
          ELSE
C--PRESERVE RAPIDITY OF CMF
            DET=SQRT(ONE+(PTX**2+PTY**2)/(PHEP(4,ICM)**2
     &                -PHEP(3,ICM)**2))
            CALL HWVSCA(2,DET,PHEP(3,ICM),PHEP(3,ICM))
          ENDIF
C---NOW BOOST TO REQUIRED Q**2 AND X-F
          PP0=PHEP(4,ICM)+PHEP(3,ICM)
          PM0=PHEP(4,ICM)-PHEP(3,ICM)
          ET0=(PP0*PM0)+ET(1)-ET(2)
          DET=ET0**2-4.*(PP0*PM0)*ET(1)
          IF (DET.LT.ZERO) THEN
            FROST=.TRUE.
            RETURN
          ENDIF
          DET=SQRT(DET)+ET0
          AL(1)= 2.*PM0*PP(1)/DET
          AL(2)=(PM0/PP(2))*(1.-2.*ET(1)/DET)
          PB(1)=0.
          PB(2)=0.
          PB(5)=2.
          DO 160 IP=1,2
          PB(3)=AL(IP)-(1./AL(IP))
          PB(4)=AL(IP)+(1./AL(IP))
          IJT=IJET(IP)
          DO 140 IHEP=IJT,IEND(IP)
          CALL HWULOF(PB,PHEP(1,IHEP),PHEP(1,IHEP))
          CALL HWULF4(PB,VHEP(1,IHEP),VHEP(1,IHEP))
  140     ISTHEP(IHEP)=ISTHEP(IHEP)+10
          CALL HWVDIF(4,VHEP(1,IPAR(IP)),VHEP(1,IJT),DISP)
          DO 150 IHEP=IJT,IEND(IP)
  150     CALL HWVSUM(4,DISP,VHEP(1,IHEP),VHEP(1,IHEP))
          IF (IEND(IP).GT.IJT+1) THEN
            ISTHEP(IJT+1)=100
          ELSEIF (IEND(IP).EQ.IJT) THEN
C---NON-PARTON JET
            ISTHEP(IJT)=3
          ENDIF
  160     CONTINUE
        ENDIF
        ISTHEP(ICM)=120
      ELSE
C---TIMELIKE JETS
C---SPECIAL CASE: IF HARD PROCESS IS W/Z DECAY, PERFORM KINEMATIC
C   RECONSTRUCTION IN ITS REST FRAME INSTEAD OF THE LAB FRAME
        IF (IDHW(ICM).GE.198.AND.IDHW(ICM).LE.200.AND.WZRFR) THEN
          CALL HWVEQU(5,PHEP(1,ICM),PLAB)
          CALL HWULOF(PLAB,PHEP(1,ICM),PHEP(1,ICM))
          CALL HWULF4(PLAB,VHEP(1,ICM),VHEP(1,ICM))
          DO 165 IP=1,NP
            CALL HWULOF(PLAB,PHEP(1,IPAR(IP)),PHEP(1,IPAR(IP)))
            CALL HWULF4(PLAB,VHEP(1,IPAR(IP)),VHEP(1,IPAR(IP)))
 165      CONTINUE
        ENDIF
C   special for DIS: preserve outgoing lepton momentum
        IF (DISPRO) THEN
          CALL HWVEQU(5,PHEP(1,IPAR(1)),PHEP(1,IJET(1)))
          ISTHEP(IJET(1))=1
          LP=2
        ELSE
          CALL HWVEQU(5,PHEP(1,ICM),PC)
C--- PQ AND PC ARE OLD AND NEW PARTON CM
          CALL HWVSUM(4,PHEP(1,IPAR(1)),PHEP(1,IPAR(2)),PQ)
          PQ(5)=PHEP(5,ICM)
          IF (NP.GT.2) THEN
            DO 170 KP=3,NP
  170       CALL HWVSUM(4,PHEP(1,IPAR(KP)),PQ,PQ)
          ENDIF
          LP=1
        ENDIF
        IF (.NOT.DISLOW) THEN
C---FIND JET CM MOMENTA
          ECM=PQ(5)
          EMS=0.
          JETRAD=.FALSE.
          DO 180 KP=LP,NP
          EMJ=PHEP(5,IJET(KP))
          EMP=PHEP(5,IPAR(KP))
          JETRAD=JETRAD.OR.EMJ.NE.EMP
          EMS=EMS+EMJ
          PM(KP)= EMJ**2
C---N.B. ROUNDING ERRORS HERE AT HIGH ENERGIES
          PJ(KP)=(HWULDO(PHEP(1,IPAR(KP)),PQ)/ECM)**2-EMP**2
          IF (PJ(KP).LE.ZERO) CALL HWWARN('HWBJCO',104,*999)
  180     CONTINUE
          PF=1.
          IF (JETRAD) THEN
C---JETS DID RADIATE
            IF (EMS.GE.ECM) THEN
              FROST=.TRUE.
              GOTO 240
            ENDIF
            DO 200 NE=1,NETRY
            EMS=-ECM
            DMS=0.
            DO 190 KP=LP,NP
            ES=SQRT(PF*PJ(KP)+PM(KP))
            EMS=EMS+ES
  190       DMS=DMS+PJ(KP)/ES
            DPF=2.*EMS/DMS
            IF (DPF.GT.PF) DPF=0.9*PF
            PF=PF-DPF
  200       IF (ABS(DPF).LT.EPS) GOTO 210
            CALL HWWARN('HWBJCO',105,*999)
          ENDIF
  210     CONTINUE
        ENDIF
C---BOOST PC AND PQ TO BREIT FRAME IF NECESSARY
        IF (DISPRO.AND.BREIT) THEN
          CALL HWULOF(PBR,PC,PC)
          CALL HWUROF(RBR,PC,PC)
          IF (.NOT.DISLOW) THEN
            CALL HWULOF(PBR,PQ,PQ)
            CALL HWUROF(RBR,PQ,PQ)
          ENDIF
        ENDIF
        DO 230 IP=LP,NP
C---FIND CM ROTATION FOR JET IP
        IF (.NOT.DISLOW) THEN
          CALL HWVEQU(4,PHEP(1,IPAR(IP)),PR)
          IF (DISPRO.AND.BREIT) THEN
            CALL HWULOF(PBR,PR,PR)
            CALL HWUROF(RBR,PR,PR)
          ENDIF
C--Modified by MHS 17/08/05 to do unboost in 2 stages (trans,long)
         PA(1)=PQ(1)
         PA(2)=PQ(2)
         PA(3)=ZERO
         PA(5)=SQRT(PQ(3)**2+PQ(5)**2)
         PA(4)=PQ(4)
         CALL HWULOF(PA,PR,PR)
         PA(1)=ZERO
         PA(2)=ZERO
         PA(3)=PQ(3)
         PA(4)=PA(5)
         PA(5)=PQ(5)
         CALL HWULOF(PA,PR,PR)
C--End mod
          CALL HWUROT(PR, ONE,ZERO,RR)
          PR(1)=ZERO
          PR(2)=ZERO
          PR(3)=SQRT(PF*PJ(IP))
          PR(4)=SQRT(PF*PJ(IP)+PM(IP))
          PR(5)=PHEP(5,IJET(IP))
          CALL HWUROB(RR,PR,PR)
C--Modified by BRW 25/10/02 to do boost in 2 stages (long,trans)
          PA(1)=ZERO
          PA(2)=ZERO
          PA(3)=PC(3)
          PA(5)=PC(5)
          PA(4)=SQRT(PA(3)**2+PA(5)**2)
          CALL HWULOB(PA,PR,PR)
          PA(1)=PC(1)
          PA(2)=PC(2)
          PA(3)=ZERO
          PA(5)=PA(4)
          PA(4)=PC(4)
          CALL HWULOB(PA,PR,PR)
C--End mod
        ELSE
          CALL HWVEQU(5,PC,PR)
        ENDIF
C---NOW PR IS LAB/BREIT MOMENTUM OF JET IP
        KP=IJET(IP)+1
        IF (AZCOR.AND.KP.LE.NHEP.AND.IDHW(KP).EQ.17) THEN
C---ALIGN CONE WITH INTERFERING PARTON
          CALL HWUROT(PR, ONE,ZERO,RS)
          JP=JMOHEP(2,KP)
          IF (JP.EQ.0) THEN
            CALL HWWARN('HWBJCO',2,*999)
            PTINF=0.
          ELSE
            CALL HWVEQU(4,PHEP(1,JP),PS)
            IF (DISPRO.AND.BREIT) THEN
              CALL HWULOF(PBR,PS,PS)
              CALL HWUROF(RBR,PS,PS)
            ENDIF
            CALL HWUROF(RS,PS,PS)
            PTINF=PS(1)**2+PS(2)**2
            IF (PTINF.LT.EPS) THEN
C---COLLINEAR JETS: ALIGN CONES
              JP=JDAHEP(1,JP)+1
              IF (ISTHEP(JP).EQ.100.AND.(ISTHEP(JP-1)+9)/10.EQ.14) THEN
                CALL HWVEQU(4,PHEP(1,JP),PS)
                IF (DISPRO.AND.BREIT) THEN
                  CALL HWULOF(PBR,PS,PS)
                  CALL HWUROF(RBR,PS,PS)
                ENDIF
                CALL HWUROF(RS,PS,PS)
                PTINF=PS(1)**2+PS(2)**2
              ELSE
                PTINF=0.
              ENDIF
            ENDIF
          ENDIF
          CALL HWVEQU(4,PHEP(1,KP),PB)
          IF (DISPRO.AND.BREIT) THEN
            CALL HWULOF(PBR,PB,PB)
            CALL HWUROF(RBR,PB,PB)
          ENDIF
          PTCON=PB(1)**2+PB(2)**2
          IF (PTCON.NE.ZERO.AND.PTINF.NE.ZERO) THEN
            CN=1./SQRT(PTINF*PTCON)
            CP=CN*(PS(1)*PB(1)+PS(2)*PB(2))
            SP=CN*(PS(1)*PB(2)-PS(2)*PB(1))
          ELSE
            CALL HWRAZM( ONE,CP,SP)
          ENDIF
        ELSE
          CALL HWRAZM( ONE,CP,SP)
        ENDIF
        CALL HWUROT(PR,CP,SP,RS)
C---FIND BOOST FOR JET IP
        ALF=(PHEP(3,IJET(IP))+PHEP(4,IJET(IP)))/
     &      (PR(4)+SQRT((PR(4)+PR(5))*(PR(4)-PR(5))))
        PB(1)=0.
        PB(2)=0.
        PB(3)=ALF-(1./ALF)
        PB(4)=ALF+(1./ALF)
        PB(5)=2.
        IHEP=IJET(IP)
        KHEP=JDAHEP(2,IHEP)
        IF (KHEP.LT.IHEP) KHEP=IHEP
        DO 220 JHEP=IHEP,KHEP
        CALL HWULOF(PB,PHEP(1,JHEP),PHEP(1,JHEP))
        CALL HWUROB(RS,PHEP(1,JHEP),PHEP(1,JHEP))
        CALL HWULF4(PB,VHEP(1,JHEP),VHEP(1,JHEP))
        CALL HWUROB(RS,VHEP(1,JHEP),VHEP(1,JHEP))
C---BOOST FROM BREIT FRAME IF NECESSARY
        IF (DISPRO.AND.BREIT) THEN
          CALL HWUROB(RBR,PHEP(1,JHEP),PHEP(1,JHEP))
          CALL HWULOB(PBR,PHEP(1,JHEP),PHEP(1,JHEP))
          CALL HWUROB(RBR,VHEP(1,JHEP),VHEP(1,JHEP))
          CALL HWULB4(PBR,VHEP(1,JHEP),VHEP(1,JHEP))
        ENDIF
        CALL HWVSUM(4,VHEP(1,JHEP),VHEP(1,IPAR(IP)),VHEP(1,JHEP))
C--MHS FIX 07/03/05 FOR VERTEX POSITION OF LONG LIVED NON-PARTON JETS
        IF (KHEP.EQ.IHEP.AND.(IDHW(JHEP).GE.121.AND.IDHW(JHEP).LE.132
     $       .OR.IDHW(JHEP).EQ.59))
     $       CALL HWVSUM(4,VTXPIP,VHEP(1,JHEP),VHEP(1,JHEP))
C--END FIX
  220   ISTHEP(JHEP)=ISTHEP(JHEP)+10
        IF (KHEP.GT.IHEP+1) THEN
          ISTHEP(IHEP+1)=100
        ELSEIF (KHEP.EQ.IHEP) THEN
C---NON-PARTON JET
          ISTHEP(IHEP)=190
        ENDIF
  230   CONTINUE
        IF (ISTHEP(ICM).EQ.110) ISTHEP(ICM)=120
C---SPECIAL CASE: FOR W/Z DECAY BOOST BACK TO THE LAB FRAME
 240    IF (IDHW(ICM).GE.198.AND.IDHW(ICM).LE.200.AND.WZRFR) THEN
          CALL HWULOB(PLAB,PHEP(1,ICM),PHEP(1,ICM))
          CALL HWULB4(PLAB,VHEP(1,ICM),VHEP(1,ICM))
          DO 260 IP=1,NP
            CALL HWULOB(PLAB,PHEP(1,IPAR(IP)),PHEP(1,IPAR(IP)))
            CALL HWULB4(PLAB,VHEP(1,IPAR(IP)),VHEP(1,IPAR(IP)))
            CALL HWULOB(PLAB,PHEP(1,IJET(IP)),PHEP(1,IJET(IP)))
C--MHS FIX 07/03/05 - DO NOT REBOOST PRIMARY VERTEX
            IF (ISTHEP(IJET(IP)).EQ.190)
     $           CALL HWVDIF(4,VHEP(1,IJET(IP)),VTXPIP,VHEP(1,IJET(IP)))
            CALL HWULB4(PLAB,VHEP(1,IJET(IP)),VHEP(1,IJET(IP)))
            IF (ISTHEP(IJET(IP)).EQ.190)
     $           CALL HWVSUM(4,VHEP(1,IJET(IP)),VTXPIP,VHEP(1,IJET(IP)))
C---END FIX
            IF (JDAHEP(1,IJET(IP)).GT.0) THEN
              IF (JDAHEP(2,IJET(IP)).GT.JDAHEP(1,IJET(IP))) THEN
                CALL HWULOB(PLAB,PHEP(1,IJET(IP)+1),PHEP(1,IJET(IP)+1))
                CALL HWULB4(PLAB,VHEP(1,IJET(IP)+1),VHEP(1,IJET(IP)+1))
              ENDIF
              DO 250 IHEP=JDAHEP(1,IJET(IP)),JDAHEP(2,IJET(IP))
                CALL HWULOB(PLAB,PHEP(1,IHEP),PHEP(1,IHEP))
                CALL HWULB4(PLAB,VHEP(1,IHEP),VHEP(1,IHEP))
 250          CONTINUE
            ENDIF
 260      CONTINUE
        ENDIF
        IF (FROST) RETURN
      ENDIF
      GOTO 20
  999 END
