      SUBROUTINE LINMIN(XPARAM,ALPHA,PVECT,NVAR,FUNCT,OKF,IC, DOTT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      DIMENSION XPARAM(NVAR),PVECT(NVAR)
      COMMON /GRAVEC/ COSINE
      COMMON /NUMCAL/ NUMCAL
C*********************************************************************
C
C  LINMIN DOES A LINE MINIMISATION.
C
C  ON INPUT:  XPARAM = STARTING COORDINATE OF SEARCH.
C             ALPHA  = STEP SIZE FOR INITIATING SEARCH.
C             PVECT  = DIRECTION OF SEARCH.
C             NVAR   = NUMBER OF VARIABLES IN XPARAM.
C             FUNCT  = INITIAL VALUE OF THE FUNCTION TO BE MINIMIZED.
C             ISOK   = NOT IMPORTANT.
C             COSINE = COSINE OF ANGLE OF CURRENT AND PREVIOUS GRADIENT.
C
C  ON OUTPUT: XPARAM = COORDINATE OF MINIMUM OF FUNCTI0N.
C             ALPHA  = NEW STEP SIZE, USED IN NEXT CALL OF LINMIN.
C             FUNCT  = FINAL, MINIMUM VALUE OF THE FUNCTION.
C             OKF    = TRUE IF LINMIN IMPROVED FUNCT, FALSE OTHERWISE.
C
C**********************************************************************
      COMMON /KEYWRD/ KEYWRD
C
C  THE FOLLOWING COMMON IS USED TO FIND OUT IF A NON-VARIATIONALLY
C  OPTIMIZED WAVE-FUNCTION IS BEING USED.
C
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,
     2                NCLOSE,NOPEN,NDUMY,FRACT
      CHARACTER KEYWRD*241
      DIMENSION PHI(3), VT(4)
      DIMENSION XSTOR(MAXPAR), XPAREF(MAXPAR)
      INTEGER LEFT,RIGHT,CENTER
      LOGICAL PRINT,OKF, HALFE, DIIS
      SAVE ICALCN, PRINT, HALFE, XMAXM, I
      SAVE MAXLIN,  YMAXST
      DATA ICALCN /0/
      IF (ICALCN.NE.NUMCAL) THEN
         HALFE =(INDEX(KEYWRD,'C.I.') .NE. 0 .OR. NCLOSE.NE.NOPEN)
         IF(INDEX(KEYWRD,'GNORM') .NE. 0)
     1DROP=DROP*MIN(READA(KEYWRD,INDEX(KEYWRD,'GNORM')),1.D0)
         XMAXM  = 0.4D0
         DELTA2 = 0.001D0
         IF(INDEX(KEYWRD,'NOTH') .EQ. 0) THEN
            DELTA1 = 0.5D0
         ELSE
            DELTA1 = 0.1D0
         ENDIF
         ALPHA  = 1.D0
         MAXLIN = 15
         IF(NVAR.EQ.1)THEN
            PVECT(1)=0.01D0
            DROP=0.01D0
            ALPHA=1.D0
            DELTA1 = 0.00005D0
            DELTA2 = 0.00001D0
            IF(INDEX(KEYWRD,'PREC') .NE. 0) DELTA1=0.0000005
            MAXLIN=30
         ENDIF
         COSINE=99.99D0
C
         YMAXST  = 0.4D0
         PRINT=(INDEX(KEYWRD,'LINMIN') .NE. 0)
         ICALCN=NUMCAL
      ENDIF
      DO 10 I=1,NVAR
         XPAREF(I)=XPARAM(I)
   10 CONTINUE         
      XMAXM=0.D0
      DO 20 I=1,NVAR
         PABS=ABS(PVECT(I))
   20 XMAXM=MAX(XMAXM,PABS)
      XMAXM=YMAXST/XMAXM
      IF(NVAR.EQ.1)
     1CALL COMPFG(XPARAM, .TRUE., FUNCT,.TRUE.,GRAD,.FALSE.)
      FIN=FUNCT
      SSQLST=FUNCT
      DIIS=IC.EQ.1.AND.NVAR.GT.1
      PHI(1)=FUNCT
      ALPHA=1.D0
      VT(1)=0.0D00
      VT(2)=ALPHA
      IF (VT(2).GT.XMAXM) VT(2)=XMAXM
      FMAX=FUNCT
      FMIN=FUNCT
      ALPHA=VT(2)
      DO 30 I=1,NVAR

   30 XPARAM(I)=XPAREF(I)+ALPHA*PVECT(I)
      CALL COMPFG(XPARAM, .TRUE., PHI(2),.TRUE.,GRAD,.FALSE.)
      IF(PHI(2).GT.FMAX) FMAX=PHI(2)
      IF(PHI(2).LT.FMIN) FMIN=PHI(2)
      CALL EXCHNG (PHI(2),SQSTOR,ENERGY,ESTOR,XPARAM,XSTOR,
     1ALPHA,ALFS,NVAR)
      IF(DIIS)GOTO 190
      IF(NVAR.GT.1)THEN
C
C   CALCULATE A NEW ALPHA BASED ON THIEL'S FORMULA
C
         ALPHA=-ALPHA**2*DOTT/(2.D0*(PHI(2)-SSQLST-ALPHA*DOTT))
         IF(ALPHA.GT.2.D0)ALPHA=2.D0
      ELSE
         IF(PHI(2).LT.PHI(1))THEN
            ALPHA=2*ALPHA
         ELSE
            ALPHA=-ALPHA
         ENDIF
      ENDIF
C#      IF(PRINT)WRITE(6,'(3(A,F12.6))')' ESTIMATED DROP:',DOTT*0.5D0,
C#     1'  ACTUAL: ',PHI(2)-SSQLST, '  PREDICTED ALPHA',ALPHA
      OKF=OKF.OR.PHI(2).LT.SSQLST
      IF(DELTA1.GT.0.3D0)THEN
C
C  THIEL'S TESTS # 18 AND 19
C
         IF(OKF.AND.ALPHA.LT.2.D0)GOTO 190
      ENDIF
      VT(3)=ALPHA
      IF (VT(3).LE.1.D0) THEN
         LEFT=3
         CENTER=1
         RIGHT=2
      ELSE
         LEFT=1
         CENTER=2
         RIGHT=3
      ENDIF
      DO 40 I=1,NVAR
   40 XPARAM(I)=XPAREF(I)+ALPHA*PVECT(I)
      CALL COMPFG (XPARAM, .TRUE., FUNCT,.TRUE.,GRAD,.FALSE.)
      IF(FUNCT.GT.FMAX) FMAX=FUNCT
      IF(FUNCT.LT.FMIN) FMIN=FUNCT
      IF (FUNCT.LT.SQSTOR) CALL EXCHNG (FUNCT,SQSTOR,ENERGY,
     1ESTOR,XPARAM,XSTOR,ALPHA,ALFS,NVAR)
      OKF=(OKF.OR.FUNCT.LT.FIN)
      PHI(3)=FUNCT
      IF (PRINT)WRITE (6,50) VT(1),PHI(1),PHI(1)-FIN,
     1                        VT(2),PHI(2),PHI(2)-FIN,
     2                        VT(3),PHI(3),PHI(3)-FIN
   50 FORMAT ( ' ---QLINMN ',/5X, 'LEFT   ...',F17.8,2F17.11/5X,
     1 'CENTER ...',F17.8,2F17.11,/5X, 'RIGHT  ...',F17.8,2F17.11,/)
      DO 180 ICTR=3,MAXLIN
         ALPHA=VT(2)-VT(3)
         BETA=VT(3)-VT(1)
         GAMMA=VT(1)-VT(2)
         IF(ABS(ALPHA*BETA*GAMMA) .GT. 1.D-4)THEN
            ALPHA=-(PHI(1)*ALPHA+PHI(2)*BETA+PHI(3)*GAMMA)/(ALPHA*BETA*G
     1AMM   A)
         ELSE
C
C   FINISH BECAUSE TWO POINTS CALCULATED ARE VERY CLOSE TOGETHER
C
            GOTO 190
         ENDIF
         BETA=((PHI(1)-PHI(2))/GAMMA)-ALPHA*(VT(1)+VT(2))
         IF (ALPHA) 60,60,90
   60    IF (PHI(RIGHT).GT.PHI(LEFT)) GO TO 70
         ALPHA=3.0D00*VT(RIGHT)-2.0D00*VT(CENTER)
         GO TO 80
   70    ALPHA=3.0D00*VT(LEFT)-2.0D00*VT(CENTER)
   80    S=ALPHA-ALPOLD
         IF (ABS(S).GT.XMAXM) S=SIGN(XMAXM,S)*(1+0.01*(XMAXM/S))
         ALPHA=S+ALPOLD
         GO TO 100
   90    ALPHA=-BETA/(2.0D00*ALPHA)
         S=ALPHA-ALPOLD
         XXM=2.0D00*XMAXM
         IF (ABS(S).GT.XXM) S=SIGN(XXM,S)*(1+0.01*(XXM/S))
         ALPHA=S+ALPOLD
  100    CONTINUE
C
C   FINISH IF CALCULATED POINT IS NEAR TO POINT ALREADY CALCULATED
C
         DO 110 I=1,3
  110    IF (ABS(ALPHA-VT(I)).LT.DELTA1*(1.D0+VT(I)).AND.OKF) GOTO 190
         DO 120 I=1,NVAR
  120    XPARAM(I)=XPAREF(I)+ALPHA*PVECT(I)
         FUNOLD=FUNCT
         CALL COMPFG (XPARAM, .TRUE., FUNCT,.TRUE.,GRAD,.FALSE.)
         IF(FUNCT.GT.FMAX) FMAX=FUNCT
         IF(FUNCT.LT.FMIN) FMIN=FUNCT
         IF (FUNCT.LT.SQSTOR) CALL EXCHNG (FUNCT,SQSTOR,ENERGY,ESTOR,
     1   XPARAM,XSTOR,ALPHA,ALFS,NVAR)
         OKF=OKF .OR. (FUNCT.LT.FIN)
         IF (PRINT) WRITE(6,130) VT(LEFT),PHI(LEFT), PHI(LEFT)-FIN,
     1                           VT(CENTER),PHI(CENTER),PHI(CENTER)-FIN,
     2                           VT(RIGHT),PHI(RIGHT),PHI(RIGHT)-FIN,
     3                           ALPHA,FUNCT,FUNCT-FIN
  130    FORMAT (5X,'LEFT    ...',F17.8,2F17.11,/5X,'CENTER  ...',
     1F17.8,2F17.11,/5X,'RIGHT   ...',F17.8,2F17.11,/5X,
     2 'NEW     ...',F17.8,2F17.11,/)
C
C TEST TO EXIT FROM LINMIN IF NOT DROPPING IN VALUE OF FUNCTION FAST.
C
         IF(ABS(FUNOLD-FUNCT) .LT. DELTA2 .AND. OKF) GOTO 190
         ALPOLD=ALPHA
         IF ((ALPHA.GT.VT(RIGHT)).OR.(ALPHA.GT.VT(CENTER)
     1        .AND.FUNCT.LT.PHI(CENTER)).OR.(ALPHA.GT.VT(LEFT)
     2        .AND.ALPHA.LT.VT(CENTER).AND.FUNCT.GT.PHI(CENTER)))
     3         GOTO 140
         VT(RIGHT)=ALPHA
         PHI(RIGHT)=FUNCT
         GO TO 150
  140    VT(LEFT)=ALPHA
         PHI(LEFT)=FUNCT
  150    IF (VT(CENTER).LT.VT(RIGHT)) GO TO 160
         I=CENTER
         CENTER=RIGHT
         RIGHT=I
  160    IF (VT(LEFT).LT.VT(CENTER)) GO TO 170
         I=LEFT
         LEFT=CENTER
         CENTER=I
  170    IF (VT(CENTER).LT.VT(RIGHT)) GO TO 180
         I=CENTER
         CENTER=RIGHT
         RIGHT=I
  180 CONTINUE
  190 CONTINUE
C
C  IC=1 IF THE LAST POINT CALCULATED WAS THE BEST POINT, IC=2 OTHERWISE
C
      IC=2
      IF(ABS(ESTOR-ENERGY).LT.1.D-12)IC=1
      CALL EXCHNG (SQSTOR,FUNCT,ESTOR,ENERGY,XSTOR,XPARAM,
     1             ALFS,ALPHA,NVAR)
      OKF = (FUNCT.LT.SSQLST.OR.DIIS)
      IF (FUNCT.GE.SSQLST) RETURN
      IF (ALPHA) 200,220,220
  200 ALPHA=-ALPHA
      DO 210 I=1,NVAR
  210 PVECT(I)=-PVECT(I)
  220 CONTINUE
      RETURN
C
C
      END
