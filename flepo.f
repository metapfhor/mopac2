      SUBROUTINE FLEPO (XPARAM,NVAR,FUNCT1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      DIMENSION XPARAM(MAXPAR)
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,
     2                NCLOSE,NOPEN,NDUMY,FRACT
      COMMON /KEYWRD/ KEYWRD
      COMMON /NUMSCF/ NSCF
      COMMON /LAST  / LAST
      COMMON /GRAVEC/ COSINE
      COMMON /PATH  / LATOM,LPARAM,REACT(200)
      COMMON /GRADNT/ GRAD(MAXPAR),GNORM
      COMMON /MESAGE/ IFLEPO,ISCF
      COMMON /TIME  / TIME0
      COMMON /FMATRX/ HESINV(MAXPAR**2+MAXPAR*3+1), IDUMY(4)
      COMMON /SCFTYP/ EMIN, LIMSCF
      COMMON /TIMDMP/ TLEFT, TDUMP
      COMMON /NUMCAL/ NUMCAL
C     PATAS
      COMMON /XXXXXX/ ICOMPF
C     PATAS
C       LAURENT MODIFICATION
      COMMON /LIN / PI,PJ,PK
      INTEGER PI,PJ,PK
      COMMON /GENRAL/ COORD(3,NUMATM), COLD(3,NUMATM*3), GOLD(MAXPAR)
      DOUBLE PRECISION LSTGRN,GMIN
      LOGICAL STBL,FAR
      SAVE LSTGRN,STBL,GMIN
      DATA GMIN,STBL,FAR /999999999.9D0,.FALSE.,.FALSE./
C       END LAURENT

      CHARACTER*241 KEYWRD
C
C
C     *
C     THIS SUBROUTINE ATTEMPTS TO MINIMIZE A REAL-VALUED FUNCTION OF
C     THE N-COMPONENT REAL VECTOR XPARAM ACCORDING TO THE
C     BFGS FORMULA. RELEVANT REFERENCES ARE
C
C     BROYDEN, C.G., JOURNAL OF THE INSTITUTE FOR MATHEMATICS AND
C                     APPLICATIONS, VOL. 6 PP 222-231, 1970.
C     FLETCHER, R., COMPUTER JOURNAL, VOL. 13, PP 317-322, 1970.
C
C     GOLDFARB, D. MATHEMATICS OF COMPUTATION, VOL. 24, PP 23-26, 1970.
C
C     SHANNO, D.F. MATHEMATICS OF COMPUTATION, VOL. 24, PP 647-656
C                    1970.
C
C   SEE ALSO SUMMARY IN
C
C    HEAD, J.D.; AND ZERNER, M.C., CHEMICAL PHYSICS LETTERS, VOL. 122,
C          264 (1985).
C    SHANNO, D.F., J. OF OPTIMIZATION THEORY AND APPLICATIONS
C          VOL.46, NO 1 PP 87-94 1985.
C     *
C     THE FUNCTION CAN ALSO BE MINIMIZED USING THE
C     DAVIDON-FLETCHER-POWELL ALGORITHM (COMPUTER JOURNAL, VOL. 6,
C     P. 163).
C
C     THE USER MUST SUPPLY THE SUBROUTINE
C     COMPFG(XPARAM,.TRUE.,FUNCT,.TRUE.,GRAD,LGRAD)
C     WHICH COMPUTES FUNCTION VALUES  FUNCT AT GIVEN VALUES FOR THE
C     VARIABLES XPARAM, AND THE GRADIENT GRAD IF LGRAD=.TRUE.
C     THE MINIMIZATION PROCEEDS BY A SEQUENCE OF ONE-DIMENSIONAL
C     MINIMIZATIONS.  THESE ARE CARRIED OUT WITHOUT GRADIENT COMPUTATION
C     BY THE SUBROUTINE LINMIN, WHICH SOLVES THE SUBPROBLEM OF
C     MINIMIZING THE FUNCTION FUNCT ALONG THE LINE XPARAM+ALPHA*PVECT,
C     WHERE XPARAM
C     IS THE VECTOR OF CURRENT VARIABLE VALUES,  ALPHA IS A SCALAR
C     VARIABLE, AND  PVECT  IS A SEARCH-DIRECTION VECTOR PROVIDED BY THE
C     BFGS OR DAVIDON-FLETCHER-POWELL ALGORITHM.  EACH ITERATION STEP CA
C     OUT BY FLEPO PROCEEDS BY LETTING LINMIN FIND A VALUE FOR ALPHA
C     WHICH MINIMIZES  FUNCT  ALONG  XPARAM+ALPHA*PVECT, BY
C     UPDATING THE VECTOR  XPARAM  BY THE AMOUNT ALPHA*PVECT, AND
C     FINALLY BY GENERATING A NEW VECTOR  PVECT.  UNDER
C     CERTAIN RESTRICTIONS (POWELL, J.INST.MATHS.APPLICS.(1971),
C     V.7,21-36)  A SEQUENCE OF FUNCT VALUES CONVERGING TO SOME
C     LOCAL MINIMUM VALUE AND A SEQUENCE OF
C     XPARAM VECTORS CONVERGING TO THE CORRESPONDING MINIMUM POINT
C     ARE PRODUCED.
C                          CONVERGENCE TESTS.
C
C     HERBERTS TEST: THE ESTIMATED DISTANCE FROM THE CURRENT POINT
C                    POINT TO THE MINIMUM IS LESS THAN DELHOF.
C
C                    "HERBERTS TEST SATISFIED - GEOMETRY OPTIMIZED"
C
C     GRADIENT TEST: THE GRADIENT NORM HAS BECOME LESS THAN TOLERG
C                    TIMES THE SQUARE ROOT OF THE NUMBER OF VARIABLES.
C
C                    "TEST ON GRADIENT SATISFIED".
C
C     XPARAM TEST:  THE RELATIVE CHANGE IN XPARAM, MEASURED BY ITS NORM,
C                   OVER ANY TWO SUCCESSIVE ITERATION STEPS DROPS BELOW
C                   TOLERX.
C
C                    "TEST ON XPARAM SATISFIED".
C
C     FUNCTION TEST: THE CALCULATED VALUE OF THE HEAT OF FORMATION
C                    BETWEEN ANY TWO CYCLES IS WITHIN TOLERF OF
C                    EACH OTHER.
C
C                    "HEAT OF FORMATION TEST SATISFIED"
C
C     FOR THE GRADIENT, FUNCTION, AND XPARAM TESTS A FURTHER CONDITION,
C     THAT NO INDIVIDUAL COMPONENT OF THE GRADIENT IS GREATER
C     THAN TOLERG, MUST BE SATISFIED, IN WHICH CASE THE
C     CALCULATION EXITS WITH THE MESSAGE
C
C                     "PETERS TEST SATISFIED"
C
C     WILL BE PRINTED, AND FUNCT AND XPARAM WILL CONTAIN THE LAST
C     FUNCTION VALUE CUM VARIABLE VALUES REACHED.
C
C
C     THE BROYDEN-FLETCHER-GOLDFARB-SHANNO AND DAVIDON-FLETCHER-POWELL
C     ALGORITHMS CHOOSE SEARCH DIRECTIONS
C     ON THE BASIS OF LOCAL PROPERTIES OF THE FUNCTION.  A MATRIX  H,
C     WHICH IN FLEPO IS PRESET WITH THE IDENTITY, IS MAINTAINED AND
C     UPDATED AT EACH ITERATION STEP.  THE MATRIX DESCRIBES A LOCAL
C     METRIC ON THE SURFACE OF FUNCTION VALUES ABOVE THE POINT XPARAM.
C     THE SEARCH-DIRECTION VECTOR  PVECT  IS SIMPLY A TRANSFORMATION
C     OF THE GRADIENT  GRAD  BY THE MATRIX H.
C
      DIMENSION XVAR(MAXPAR), GVAR(MAXPAR), XD(MAXPAR), GD(MAXPAR),
     1GLAST(MAXPAR), XLAST(MAXPAR), GG(MAXPAR), PVECT(MAXPAR)
      DIMENSION MDFP(9),XDFP(9), XTEMP(MAXPAR), GTEMP(MAXPAR)
      SAVE  ICALCN
      SAVE  RST, TDEL, SFACT, DELL, EINC, IGG1, DEL
      SAVE  RESTRT, GEOOK, DFP, CONST
      SAVE  SADDLE, MINPRT, ROOTV, PRINT, DELHOF, TOLERF, TOLERG
      SAVE  TOLERX, DROP, FREPF, IHDIM, CNCADD, ABSMIN, ITRY1
      SAVE  OKF, JCYC, LNSTOP, IREPET, ALPHA, PNORM, JNRST, CYCMX
      SAVE  COS, NCOUNT, RESFIL, MDFP, TX1, TX2, TLAST
      SAVE  TOTIME
      LOGICAL OKF, PRINT,  RESTRT, MINPRT, SADDLE, GEOOK, LOG
     1        ,RESFIL, LGRAD, DFP, LDIIS, THIEL, DIISOK, FRST,
     2         LIMSCF
      EQUIVALENCE (MDFP(1),JCYC  ),(MDFP(2),JNRST),(MDFP(3),NCOUNT),
     1            (MDFP(4),LNSTOP),(XDFP(1),ALPHA),(XDFP(2),COS   ),
     2            (XDFP(3),PNORM ),(XDFP(4),DROP ),(XDFP(5),DEL   ),
     3            (XDFP(6),FREPF ),(XDFP(7),CYCMX),(XDFP(8),TOTIME)
      DATA ICALCN /0/
C
C   START OF ONCE-ONLY SECTION
C
      EMIN=0.D0
      IF (ICALCN.NE.NUMCAL) THEN
C
C   THE FOLLOWING CONSTANTS SHOULD BE SET BY THE USER.
C
         RST   = 0.05D0
         IPRT  = 6
         TDEL  = 0.06D0
         NRST  = 30
         SFACT = 1.5
         DELL  = 0.01D0
         EINC  = 0.3D0
         IGG1  = 3
         DEL=DELL
C
C    THESE CONSTANTS SHOULD BE SET BY THE PROGRAM.
C
         RESTRT = INDEX(KEYWRD,'RESTAR').NE.0
         THIEL  = INDEX(KEYWRD,'NOTHIE').EQ.0
         GEOOK  = INDEX(KEYWRD,'GEO-OK').NE.0
         LOG    = INDEX(KEYWRD,'NOLOG').EQ.0
         LDIIS  = INDEX(KEYWRD,'NODIIS').EQ.0
         SADDLE = INDEX(KEYWRD,'SADDLE').NE.0
         MINPRT = .NOT.SADDLE
C         CONST=1.D0
         CONST=1.D-4
C
C      THE DAVIDON-FLETCHER-POWELL METHOD IS NOT RECOMMENDED
C      BUT CAN BE INVOKED BY USING THE KEYWORD 'DFP'
C
         DFP=INDEX(KEYWRD,'DFP').NE.0
C
C  ORDER OF PRECISION:   'GNORM' TAKES PRECEDENCE OVER 'FORCE', WHICH
C                        TAKES PRECEDENCE OVER 'PRECISE'.

C         TOLERG=1.D0 Laurent Modification
         TOLERG=1.D-1
         IF(INDEX(KEYWRD,'PREC') .NE. 0) TOLERG=TOLERG/5.0
         IF (INDEX(KEYWRD,'FORCE') .NE. 0) TOLERG = TOLERG/10.0
C
C      READ IN THE GRADIENT-NORM LIMIT, IF SPECIFIED
C
         IF(INDEX(KEYWRD,'GNORM=').NE.0) THEN
            ROOTV=1.D0
            CONST=1.D-20
            TOLERG=READA(KEYWRD,INDEX(KEYWRD,'GNORM='))
            IF(INDEX(KEYWRD,' LET').EQ.0.AND.TOLERG.LT.1.D-2)THEN
               WRITE(6,'(/,A)')'  GNORM HAS BEEN SET TOO LOW, RESET TO 0
     1.01'
               TOLERG=1.D-2
            ENDIF
         ELSE
            ROOTV=SQRT(NVAR+1.D-5)
         ENDIF
C         LAURENT MODIFICATION
         TOLERX =  0.0001D0*CONST
C         TOLERX =  0.00001D0*CONST
         DELHOF = 0.0010D0*CONST
         TOLERF = 0.002D0*CONST
         TOLRG  = TOLERG
C
C  MINOR BOOK-KEEPING
C
         TLAST=TLEFT
         TX2=SECOND()
         TLEFT=TLEFT-TX2+TIME0
         PRINT  = (INDEX(KEYWRD,'FLEPO').NE.0)
C
C   THE FOLLOWING CONSTANTS SHOULD BE SET TO SOME ARBITARY LARGE VALUE.
C
         DROP  = 1.D15
         FREPF = 1.D15
C
C     AND FINALLY, THE FOLLOWING CONSTANTS ARE CALCULATED.
C
         IHDIM=(NVAR*(NVAR+1))/2
         CNCADD=1.0D00/ROOTV
         IF (CNCADD.GT.0.15D00) CNCADD=0.15D00
         ICALCN=NUMCAL
         IF (RESTRT) THEN
            JNRST=1
            MDFP(9)=0
            CALL DFPSAV(TOTIME,XPARAM,GD,XLAST,FUNCT1,MDFP,XDFP)
            I=TOTIME/1000000.D0
            TOTIME=TOTIME-I*1000000.D0
            TIME0=TIME0-TOTIME
            NSCF=MDFP(5)
            WRITE(IPRT,'(//10X,''TOTAL TIME USED SO FAR:'',
     1    F13.2,'' SECONDS'')')TOTIME
            IF(INDEX(KEYWRD,'1SCF') .NE. 0) THEN
               LAST=1
               LGRAD= INDEX(KEYWRD,'GRAD').NE.0
               CALL COMPFG (XPARAM,.TRUE.,FUNCT1,.TRUE.,GRAD,LGRAD)
               IFLEPO=13
               EMIN=0.D0
               RETURN
            ENDIF
         ENDIF

C       LAURENT MODIFICATION: WRITE OUT INITIAL GEOMETRY
C        CALL WRTXYZ()
C       END LAURENT

C
C   END OF ONCE-ONLY SETUP
C
      ENDIF

C
C     FIRST, WE INITIALIZE THE VARIABLES.
C
      DIISOK=.FALSE.
      IRESET=0
      ABSMIN=1.D6
      FRST=.TRUE.
      ITRY1=0
      JCYC=0
      LNSTOP=1
      IREPET=1
      LIMSCF=.TRUE.
      ALPHA = 1.0D00
      PNORM=1.0D00
      JNRST=0
      CYCMX=0.D0
      COS=0.0D00
      TOTIME=0.D0
      NCOUNT=1
      IF( SADDLE) THEN
*
*   WE DON'T NEED HIGH PRECISION DURING A SADDLE-POINT CALCULATION.
*
         IF(NVAR.GT.0)GNORM=SQRT(DOT(GRAD,GRAD,NVAR))-3.D0
         IF(GNORM.GT.10.D0)GNORM =10.D0
         IF(GNORM.GT.1.D0) TOLERG=TOLRG*GNORM
         WRITE(IPRT,'('' GRADIENT CRITERION IN FLEPO ='',F12.5)')TOLE
     1RG
      ENDIF
      IF(NVAR.EQ.1) THEN
         PVECT(1)=0.01D0
         ALPHA=1.D0
         GOTO 300
      ENDIF
      TOTIME=0.D0
C
C CALCULATE THE VALUE OF THE FUNCTION -> FUNCT1, AND GRADIENTS -> GRAD.
C NORMAL SET-UP OF FUNCT1 AND GRAD, DONE ONCE ONLY.
C
      ICOMPF=1
      CALL COMPFG (XPARAM,.TRUE.,FUNCT1,.TRUE.,GRAD,.TRUE.)
      CALL SCOPY (NVAR,GRAD,1,GD,1)
      IF (NVAR.NE.0) THEN
         GNORM=SQRT(DOT(GRAD,GRAD,NVAR))
         GNORMR=GNORM
         IF (LNSTOP.NE.1.AND.COS.GT.RST.AND.(JNRST.LT.NRST.OR..NOT.DFP)
     1     .AND.RESTRT)THEN
            CALL SCOPY (NVAR,GD,1,GLAST,1)
         ELSE
            CALL SCOPY (NVAR,GRAD,1,GLAST,1)
         ENDIF
      ENDIF
      IF(GNORM.LT.TOLERG.OR.NVAR.EQ.0) THEN
         IFLEPO=2
         IF(RESTRT) THEN
            CALL COMPFG (XPARAM,.TRUE.,FUNCT1,.TRUE.,GRAD,.TRUE.)
         ELSE
      ICOMPF=1+ICOMPF
            CALL COMPFG (XPARAM,.TRUE.,FUNCT1,.TRUE.,GRAD,.FALSE.)
         ENDIF
         EMIN=0.D0
         RETURN
      ENDIF
      TX1 =  SECOND()
      TLEFT=TLEFT-TX1+TX2
C     *
C     START OF EACH ITERATION CYCLE ...
C     *
C
    5 GNORM=SQRT(DOT(GRAD,GRAD,NVAR))
      IF(GNORMR.LT.1.D-10)GNORMR=GNORM
      GOTO 30
   10 CONTINUE
      IF(COS .LT. RST) THEN
         DO 20 I=1,NVAR
   20    GD(I)=0.5D0
      ENDIF
   30 CONTINUE
      JCYC=JCYC+1
      JNRST=JNRST+1
      I80=0
   40 CONTINUE
      IF (I80.EQ.1.OR.
     1LNSTOP.EQ.1.OR.COS.LE.RST.OR.(JNRST.GE.NRST.AND.DFP))THEN
C
C     *
C     RESTART SECTION
C     *
C
         DO 50 I=1,NVAR
C
C  MAKE THE FIRST STEP A WEAK FUNCTION OF THE GRADIENT
C
            STEP=ABS(GRAD(I))*0.0002D0
            STEP=MAX(0.01D0,MIN(0.04D0,STEP))
C#         XD(I)=XPARAM(I)-SIGN(STEP,GRAD(I))
            XD(I)=XPARAM(I)-SIGN(DEL,GRAD(I))
   50    CONTINUE
C#      WRITE(6,'(10F8.3)')(XD(I)-XPARAM(I),I=1,NVAR)
C
C THIS CALL OF COMPFG IS USED TO CALCULATE THE SECOND-ORDER MATRIX IN H
C IF THE NEW POINT HAPPENS TO IMPROVE THE RESULT, THEN IT IS KEPT.
C OTHERWISE IT IS SCRAPPED, BUT STILL THE SECOND-ORDER MATRIX IS O.K.
C
C#      WRITE(6,*)' RESET HESSIAN'
         CALL COMPFG (XD,.TRUE.,FUNCT2,.TRUE.,GD,.TRUE.)
         IF(.NOT. GEOOK .AND. SQRT(DOT(GD,GD,NVAR))/GNORM.GT.10.
     1 AND.GNORM.GT.20.AND.JCYC.GT.2)THEN
C
C  THE GEOMETRY IS BADLY SPECIFIED IN THAT MINOR CHANGES IN INTERNAL
C  COORDINATES LEAD TO LARGE CHANGES IN CARTESIAN COORDINATES, AND THESE
C  LARGE CHANGES ARE BETWEEN PAIRS OF ATOMS THAT ARE CHEMICALLY BONDED
C  TOGETHER.
            WRITE(IPRT,'('' GRADIENTS OF OLD GEOMETRY, GNORM='',F13.6)')
     1              GNORM
            WRITE(IPRT,'(6F12.6)')(GRAD(I),I=1,NVAR)
            GDNORM=SQRT(DOT(GD,GD,NVAR))
            WRITE(IPRT,'('' GRADIENTS OF NEW GEOMETRY, GNORM='',F13.6)')
     1              GDNORM
            WRITE(IPRT,'(6F12.6)')(GD(I),I=1,NVAR)
            WRITE(IPRT,'(///20X,''CALCULATION ABANDONED AT THIS POINT!''
     1)')
            WRITE(IPRT,'(//10X,'' SMALL CHANGES IN INTERNAL COORDINATES
     1ARE   '',/10X,'' CAUSING A LARGE CHANGE IN THE DISTANCE BETWEEN'',
     2/   10X,'' CHEMICALLY-BOUND ATOMS. THE GEOMETRY OPTIMIZATION'',/
     3   10X,'' PROCEDURE WOULD LIKELY PRODUCE INCORRECT RESULTS'')')
            CALL GEOUT(1)
            STOP
         ENDIF
         NCOUNT=NCOUNT+1
         DO 60 I=1,IHDIM
   60    HESINV(I)=0.0D00
         II=0
         DO 90 I=1,NVAR
            II=II+I
            DELTAG=GRAD(I)-GD(I)
            DELTAX=XPARAM(I)-XD(I)
            IF (ABS(DELTAG).LT.1.D-12) GO TO 70
            GGD=ABS(GRAD(I))
            IF (FUNCT2.LT.FUNCT1) GGD=ABS(GD(I))
            HESINV(II)=DELTAX/DELTAG
            IF (HESINV(II).LT.0.0D00.AND.GGD.LT.1.D-12) GO TO 70
            IF (HESINV(II).LT.0.0D00) HESINV(II)=TDEL/GGD
            GO TO 80
   70       HESINV(II)=0.01D00
   80       CONTINUE
            IF (GGD.LT.1.D-12) GGD=1.D-12
            PMSTEP=ABS(0.1D0/GGD)
            IF (HESINV(II).GT.PMSTEP) HESINV(II)=PMSTEP
   90    CONTINUE
         JNRST=0
         IF(JCYC.LT.2)COSINE=1.D0
         IF(FUNCT2 .GE. FUNCT1) THEN
            IF(PRINT)WRITE (IPRT,100) FUNCT1,FUNCT2
  100       FORMAT (' FUNCTION VALUE=',F13.7,
     1           '  WILL NOT BE REPLACED BY VALUE=',F13.7,/10X,
     2           'CALCULATED BY RESTART PROCEDURE',/)
            COSINE=1.D0
         ELSE
            IF( PRINT ) WRITE (IPRT,110) FUNCT1,FUNCT2
  110       FORMAT (' FUNCTION VALUE=',F13.7,
     1           ' IS BEING REPLACED BY VALUE=',F13.7,/10X,
     2           ' FOUND IN RESTART PROCEDURE',/,6X,'THE CORRESPONDING'
     3           ' X VALUES AND GRADIENTS ARE ALSO BEING REPLACED',/)
            FUNCT1=FUNCT2
            CALL SCOPY (NVAR,XD,1,XPARAM,1)
            CALL SCOPY (NVAR,GD,1,GRAD  ,1)
            GNORM=SQRT(DOT(GRAD,GRAD,NVAR))
            IF(GNORMR.LT.1.D-10)GNORMR=GNORM
         ENDIF
      ELSE
C
C     *
C     UPDATE VARIABLE-METRIC MATRIX
C     *
C
         DO 120 I=1,NVAR
            XVAR(I)=XPARAM(I)-XLAST(I)
  120    GVAR(I)=GRAD(I)-GLAST(I)
         CALL SUPDOT(GG,HESINV,GVAR,NVAR,1)
         YHY=DOT(GG,GVAR,NVAR)
         SY =DOT(XVAR,GVAR,NVAR)
         K=0
C
C    UPDATE ACCORDING TO DAVIDON-FLETCHER-POWELL
C
         IF(DFP)THEN
            DO 130 I=1,NVAR
               XVARI=XVAR(I)/SY
               GGI=GG(I)/YHY
               DO 130 J=1,I
                  K=K+1
  130       HESINV(K)=HESINV(K)+XVAR(J)*XVARI-GG(J)*GGI
C
C     UPDATE USING THE BFGS FORMALISM
C
         ELSE
            YHY=1.0D0 + YHY/SY
            DO 140 I=1,NVAR
               XVARI=XVAR(I)/SY
               GGI=GG(I)/SY
               DO 140 J=1,I
                  K=K+1
  140       HESINV(K)=HESINV(K)-GG(J)*XVARI-XVAR(J)*GGI + YHY*XVAR(J)*XV
     1ARI
         ENDIF
      ENDIF
C#      DO 191 I=1,IHDIM
C#  191 HTEMP(I)=HESINV(I)
C#      CALL HQRII(HTEMP, NVAR, NVAR, XTEMP, VECTS)
C#      J=0
C#      DO 193 I=1,NVAR
C#      IF(XTEMP(I).LT.0.0D0)THEN
C#      J=J+1
C#      XTEMP(I)=0.00002D0
C#      ENDIF
C#  193 CONTINUE
C#      IF(J.NE.0)THEN
C#      DO 194 I=1,IHDIM
C#  194 HTEMP(I)=HESINV(I)
C#      CALL HREFM(NVAR,VECTS,XTEMP,HESINV)
C#      WRITE(6,*)' ORIGINAL HESSIAN'
C#      CALL VECPRT(HTEMP,NVAR)
C#      WRITE(6,*)' REFORMED HESSIAN'
C#      CALL VECPRT(HESINV,NVAR)
C#      ENDIF
C#      WRITE(6,*)' EIGENVALUES OF HESSIAN MATRIX'
C#      WRITE(6,'(1X,5G12.6)')(6.951D-3/XTEMP(I),I=1,NVAR)
C
C     *
C     ESTABLISH NEW SEARCH DIRECTION
C     *
      PNLAST=PNORM
C#      call vecprt(hesinv,nvar)

      CALL SUPDOT(PVECT,HESINV,GRAD,NVAR,1)
      PNORM=SQRT(DOT(PVECT,PVECT,NVAR))
      IF(PNORM.GT.1.5D0*PNLAST)THEN
C
C  TRIM PVECT BACK
C
         DO 150 I=1,NVAR
  150    PVECT(I)=PVECT(I)*1.5D0*PNLAST/PNORM
         PNORM=1.5D0*PNLAST
      ENDIF
      DOTT=-DOT(PVECT,GRAD,NVAR)
      DO 160 I=1,NVAR
  160 PVECT(I)=-PVECT(I)
      COS=-DOTT/(PNORM*GNORM)
      IF (JNRST.EQ.0) GO TO 190
      IF (COS.LE.CNCADD.AND.DROP.GT.1.0D00) GO TO 170
      IF (COS.LE.RST) GO TO 170
      GO TO 190
  170 CONTINUE
C#      K=0
C#      DO 222 I=1,NVAR
C#      DO 223 J=1,I-1
C#      K=K+1
C#  223 HESINV(K)=HESINV(K)*0.75D0
C#      K=K+1
C#  222 HESINV(K)=HESINV(K)+0.005D0
C#      GOTO 241
      PNORM=PNLAST
      IF( PRINT )WRITE (IPRT,180) COS
  180 FORMAT (//,5X, 'SINCE COS=',F9.3,5X,'THE PROGRAM WILL GO TO RE',
     1'START SECTION',/)
      I80=1
      GO TO 40
  190 CONTINUE
      IF ( PRINT ) WRITE (IPRT,200) JCYC,FUNCT1
  200 FORMAT (1H , 'AT THE BEGINNING OF CYCLE',I5, '  THE FUNCTION VA
     1LUE IS ',F13.6/, '  THE CURRENT POINT IS ...')
      IF(PRINT)WRITE (IPRT,210) GNORM,COS
  210 FORMAT ( '  GRADIENT NORM = ',F10.4/,'  ANGLE COSINE =',F10.4)
      IF( PRINT )THEN
         WRITE (6,220)
  220    FORMAT ('  THE CURRENT POINT IS ...')
         NTO6=(NVAR-1)/6+1
         IINC1=-5
         DO 270 I=1,NTO6
            WRITE (6,'(/)')
            IINC1=IINC1+6
            IINC2=MIN(IINC1+5,NVAR)
            WRITE (6,230) (J,J=IINC1,IINC2)
            WRITE (6,240) (XPARAM(J),J=IINC1,IINC2)
            WRITE (6,250) (GRAD(J),J=IINC1,IINC2)
            WRITE (6,260) (PVECT(J),J=IINC1,IINC2)
  230       FORMAT (1H ,3X,  1HI,9X,I3,9(8X,I3))
  240       FORMAT (1H ,1X, 'XPARAM(I)',1X,F9.4,2X,9(F9.4,2X))
  250       FORMAT (1H ,1X, 'GRAD  (I)',F10.4,1X,9(F10.4,1X))
  260       FORMAT (1H ,1X, 'PVECT (I)',2X,F10.6,1X,9(F10.6,1X))
  270    CONTINUE
      ENDIF
      LNSTOP=0
      ALPHA=ALPHA*PNLAST/PNORM
      CALL SCOPY (NVAR,GRAD,  1,GLAST,1)
      CALL SCOPY (NVAR,XPARAM,1,XLAST,1)
      IF (JNRST.EQ.0) ALPHA=1.0D00
      DROP=ABS(ALPHA*DOTT)
      IF(PRINT)WRITE (IPRT,280) DROP
  280 FORMAT (1H , 13H -ALPHA.P.G =,F18.6,/)
      IF (JNRST.NE.0.AND.DROP.LT.DELHOF)THEN
C
C   HERBERT'S TEST: THE PREDICTED DROP IN ENERGY IS LESS THAN DELHOF
C   IF PASSED, CALL COMPFG TO GET A GOOD SET OF EIGENVECTORS, THEN EXIT
C
         IF(MINPRT)WRITE (IPRT,290)
  290    FORMAT(//,10X,'HERBERTS TEST SATISFIED - GEOMETRY OPTIMIZED')
C
C   FLEPO IS ENDING PROPERLY. THIS IS IMMEDIATELY BEFORE THE RETURN.
C
         LAST=1
         CALL COMPFG (XPARAM,.TRUE.,FUNCT,.TRUE.,GRAD,.FALSE.)
         IFLEPO=3
         TIME0=TIME0-TOTIME
         EMIN=0.D0
         RETURN
      ENDIF
      BETA =ALPHA
      SMVAL=FUNCT1
      DROPN=-ABS(DROP/ALPHA)
C
C    UPDATE GEOMETRY USING THE G-DIIS PROCEDURE
C
      IF(DIISOK) THEN
         OKF=.TRUE.
         IC=1
      ELSE
         OKF=.FALSE.
         IC=2
      ENDIF
  300 CALL LINMIN(XPARAM,ALPHA,PVECT,NVAR,FUNCT1,OKF,IC,DROPN)





      IF(NVAR.EQ.1)THEN
         WRITE(6,'('' ONLY ONE VARIABLE, THEREFORE ENERGY A MINIMUM'')')
         LAST=1
         LGRAD=(INDEX(KEYWRD,'GRAD').NE.0)
         CALL COMPFG (XPARAM,.TRUE.,FUNCT,.TRUE.,GRAD,LGRAD)
         IFLEPO=14
         EMIN=0.D0
         RETURN
      ENDIF
C   WE WANT ACCURATE DERIVATIVES AT THIS POINT
C
C   LINMIN DOES NOT GENERATE ANY DERIVATIVES, THEREFORE COMPFG MUST BE
C   CALLED TO END THE LINE SEARCH
C
C  IF THE DERIVATIVES ARE TO BE CALCULATED USING FULL SCF'S, THEN CHECK
C  WHETHER TO DO FULL SCF'S (CRITERION FROM FLEPO: GRAD IS NULL).
C
      IF(IRESET.GT.10.OR.GNORM.LT.40.D0.AND.GNORM/GNORMR.LT.0.33D0)THEN
         IRESET=0
         GNORMR=0.D0
         DO 310 I=1,NVAR
  310    GRAD(I)=0.D0
      ENDIF
      IRESET=IRESET+1
C
C
C     RESTORE TO STANDARD VALUE BEFORE COMPUTING THE GRADIENT
      IF(THIEL)THEN
         CALL COMPFG (XPARAM, IC.NE.1,  SCRAP, .TRUE. ,GRAD,.TRUE.)
      ELSE
         CALL COMPFG (XPARAM, .TRUE.,  FUNCT1, .TRUE. ,GRAD,.TRUE.)
      ENDIF
      IF(LDIIS) THEN
C
C  UPDATE GEOMETRY AND GRADIENT AFTER MAKING A STEP USING LINMIN
C
         DO 320 I=1,NVAR
            XTEMP(I)=XPARAM(I)
  320    GTEMP(I)=GRAD(I)
!         CALL DIIS(XTEMP, XPARAM, GTEMP, GRAD, HDIIS, FUNCT1,HESINV, NVA
!     1R, FRST)
         IF(HDIIS.LT.FUNCT1.AND.
     1SQRT(DOT(GTEMP,GTEMP,NVAR)) .LT. SQRT(DOT(GRAD,GRAD,NVAR)))THEN
            DO 330 I=1,NVAR
               XPARAM(I)=XTEMP(I)
  330       GRAD(I)=GTEMP(I)
            DIISOK=.TRUE.
         ELSE
            DIISOK=.FALSE.
         ENDIF
      ENDIF
      GNORM=SQRT(DOT(GRAD,GRAD,NVAR))
      IF(GNORMR.LT.1.D-10)GNORMR=GNORM
      NCOUNT=NCOUNT+1
      IF ( .NOT. OKF) THEN
         LNSTOP = 1
         IF(MINPRT)WRITE (IPRT,'(/,20X, ''NO POINT LOWER IN ENERGY '',
     1    ''THAN THE STARTING POINT '',/,20X,''COULD BE FOUND '',
     2    ''IN THE LINE MINIMIZATION'')')
         FUNCT1=SMVAL
         ALPHA=BETA
         CALL SCOPY (NVAR,GLAST,1,GRAD  ,1)
         CALL SCOPY (NVAR,XLAST,1,XPARAM,1)
C       LAURENT MODIFACTION: WE DONT WANT TO END SO QUICKLY
         IF (JNRST.EQ.0)THEN

            WRITE (IPRT,340)
  340       FORMAT (1H ,//,20X, 'SINCE COS WAS JUST RESET,THE SEARCH',
     1        ' IS BEING ENDED')
C
C           FLEPO IS ENDING BADLY. THIS IS IMMEDIATELY BEFORE THE RETURN
C
            LAST=1
            CALL COMPFG (XPARAM, .TRUE., FUNCT, .TRUE. ,GRAD,.TRUE.)
            IFLEPO=4
            TIME0=TIME0-TOTIME
            EMIN=0.D0
            RETURN
         ENDIF
C       END LAURENT
         IF(PRINT)WRITE (IPRT,350)
  350    FORMAT (1H ,20X, 'COS WILL BE RESET AND ANOTHER '
     1    ,'ATTEMPT MADE')
         COS=0.0D00
         GO TO 470
      ENDIF
      XN=SQRT(DOT(XPARAM,XPARAM,NVAR))
      TX=ABS(ALPHA*PNORM)
      IF (XN.NE.0.0D00) TX=TX/XN
      TF=ABS(SMVAL-FUNCT1)
C       IF(ABSMIN-SMVAL.LT.1.D-7)THEN       Laurent Modification: this def of stationary makes more sense
      IF(ABS(ABSMIN-SMVAL).LT.1.D-7)THEN
         ITRY1=ITRY1+1
         IF(ITRY1.GT.10)THEN
            WRITE(6,'(//,'' HEAT OF FORMATION IS ESSENTIALLY STATIONARY'
     1')')
            GOTO 460
         ENDIF
      ELSE
         ITRY1=0
         ABSMIN=SMVAL
      ENDIF
      IF (PRINT) WRITE (6,360) NCOUNT,COS,TX*XN,ALPHA,-DROP,-TF,GNORM
  360 FORMAT (/,'           NUMBER OF COUNTS =',I6,
     1'         COS    =',F11.4,/,
     2        '  ABSOLUTE  CHANGE IN X     =',F13.6,
     3'  ALPHA  =',F11.4,/,
     4        '  PREDICTED CHANGE IN F     =  ',G11.4,
     5'  ACTUAL =  ',G11.4,/,
     6        '  GRADIENT NORM             =  ',G11.4,//)
      IF (TX.LE.TOLERX) THEN
         IF(MINPRT) WRITE (IPRT,370)
  370    FORMAT (' TEST ON X SATISFIED')
         GOTO 400
      ENDIF
      IF (TF.LE.TOLERF) THEN
C#         WRITE(6,*)TF,TOLERF
         IF(MINPRT) WRITE (IPRT,380)
  380    FORMAT (' HEAT OF FORMATION TEST SATISFIED')
         GOTO 400
      ENDIF
      IF (GNORM.LE.TOLERG*ROOTV) THEN
         IF(MINPRT) WRITE (IPRT,390)
  390    FORMAT (' TEST ON GRADIENT SATISFIED')
         GOTO 400
      ENDIF
      GO TO 470
  400 DO 440 I=1,NVAR
         IF (ABS(GRAD(I)).GT.TOLERG)THEN
            IREPET=IREPET+1
            IF (IREPET.GT.1) GO TO 410
            FREPF=FUNCT1
            COS=0.0D00
  410       IF(MINPRT) WRITE (IPRT,420)TOLERG
  420       FORMAT (20X,'HOWEVER, A COMPONENT OF GRADIENT IS ',
     1     'LARGER THAN',F6.2 ,/)
            IF (ABS(FUNCT1-FREPF).GT.EINC) IREPET=0
            IF (IREPET.GT.IGG1) THEN
               WRITE (IPRT,430)IGG1,EINC
  430          FORMAT (10X,' THERE HAVE BEEN',I2,' ATTEMPTS TO REDUCE TH
     1E ',' GRADIENT.',/10X,' DURING THESE ATTEMPTS THE ENERGY DROPPED',
     2' BY LESS THAN',F4.1,' KCAL/MOLE',/
     310X,' FURTHER CALCULATION IS NOT JUSTIFIED AT THIS TIME.',/
     410X,' TO CONTINUE, START AGAIN WITH THE WORD "PRECISE"' )
               LAST=1
               CALL COMPFG (XPARAM,.TRUE.,FUNCT,.TRUE.,GRAD,.FALSE.)
               IFLEPO=8
               TIME0=TIME0-TOTIME
               EMIN=0.D0
               RETURN
            ELSE
               GOTO 470
            ENDIF
         ENDIF
  440 CONTINUE
      IF(MINPRT) WRITE (IPRT,450)
  450 FORMAT ( 23H PETERS TEST SATISFIED )
  460 LAST=1
      CALL COMPFG (XPARAM,.TRUE.,FUNCT,.TRUE.,GRAD,.FALSE.)
      IFLEPO=6
      TIME0=TIME0-TOTIME
      EMIN=0.D0
      RETURN
C
C   ALL TESTS HAVE FAILED, WE NEED TO DO ANOTHER CYCLE.
C
  470 CONTINUE
      BSMVF=ABS(SMVAL-FUNCT1)
      IF (BSMVF.GT.10.D00) COS = 0.0D00
      DEL=0.002D00
      IF (BSMVF.GT.1.0D00) DEL=DELL/2.0D00
      IF (BSMVF.GT.5.0D00) DEL=DELL
      TX2 = SECOND ()
      TCYCLE=TX2-TX1
      TX1=TX2
C       LAURENT MODIFICATION: WRITEOUT THE GEOMETRY
      CALL WRTXYZ(JCYC)
C       END LAURENT
C END OF ITERATION LOOP, EVERYTHING IS STILL O.K. SO GO TO
C NEXT ITERATION, IF THERE IS ENOUGH TIME LEFT.
C
      IF(TCYCLE.LT.100000.D0)CYCMX=MAX(CYCMX,TCYCLE)
      TLEFT=TLEFT-TCYCLE
      IF(TLEFT.LT.0)TLEFT=-0.1D0
      IF(TCYCLE.GT.1.D5)TCYCLE=0.D0
      IF(TLAST-TLEFT.GT.TDUMP)THEN
         TOTIM=TOTIME   +   SECOND()-TIME0
         TLAST=TLEFT
         MDFP(9)=2
         RESFIL=.TRUE.
         MDFP(5)=NSCF
         CALL DFPSAV(TOTIM,XPARAM,GD,XLAST,FUNCT1,MDFP,XDFP)
      ENDIF
      IF(RESFIL)THEN
         IF(MINPRT) WRITE(6,480)MIN(TLEFT,9999999.9D0),
     1MIN(GNORM,999999.999D0),FUNCT1
  480    FORMAT(' RESTART FILE WRITTEN,   TIME LEFT:',F9.1,
     1' GRAD.:',F10.3,' HEAT:',G13.7)
         RESFIL=.FALSE.
      ELSE

         IF(MINPRT) WRITE(6,490)JCYC,MIN(TCYCLE,9999.99D0),
     1MIN(TLEFT,9999999.9D0),MIN(GNORM,999999.999D0),FUNCT1
         IF(LOG) WRITE(11,490)JCYC,MIN(TCYCLE,9999.99D0),
     1MIN(TLEFT,9999999.9D0),MIN(GNORM,999999.999D0),FUNCT1
  490    FORMAT(' CYCLE:',I4,' TIME:',F7.2,' TIME LEFT:',F9.1,
     1' GRAD.:',F10.3,' HEAT:',G13.7)
      ENDIF
      IF (TLEFT.GT.SFACT*CYCMX) GO TO 10
      WRITE(IPRT,500)
  500 FORMAT (20X, 42HTHERE IS NOT ENOUGH TIME FOR ANOTHER CYCLE,/,30X,
     118HNOW GOING TO FINAL)
      TOTIM=TOTIME   +   SECOND()-TIME0
      MDFP(9)=1
      MDFP(5)=NSCF
      CALL DFPSAV(TOTIM,XPARAM,GD,XLAST,FUNCT1,MDFP,XDFP)
      IFLEPO=-1
      RETURN
C
C
      END
