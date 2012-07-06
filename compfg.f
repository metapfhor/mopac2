      SUBROUTINE COMPFG(XPARAM,INT,ESCF,FULSCF,GRAD,LGRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      DIMENSION XPARAM(MAXPAR),GRAD(MAXPAR)
      LOGICAL LGRAD, FULSCF
      COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR),IDUMY,DUMY(MAXPAR)
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),LOCDEP(MAXPAR)
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /ATHEAT/ ATHEAT
      COMMON /WMATRX/ WJ(N2ELEC), WK(N2ELEC)
      COMMON /ENUCLR/ ENUCLR
      COMMON /NATYPE/ NZTYPE(107),MTYPE(30),LTYPE
      COMMON /ELECT / ELECT
      PARAMETER (MDUMY=MAXPAR**2-MPACK)
      COMMON /SCRACH/ RXYZ(MPACK), XDUMY(MDUMY)
      COMMON /HMATRX/ H(MPACK)
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1                NA(NUMATM), NB(NUMATM), NC(NUMATM)
      COMMON /ERRFN  / ERRFN(MAXPAR), AICORR(MAXPAR)
      COMMON /VECTOR/ C(MORB2),EIGS(MAXORB),CBETA(MORB2),EIGB(MAXORB)
      COMMON /LAST  / LAST
      COMMON /NUMCAL/ NUMCAL
      COMMON /SCFTYP/ EMIN, LIMSCF
      COMMON /MOLMEC/ HTYPE(4),NHCO(4,20),NNHCO,ITYPE
     1       /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     2                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,
     3                NCLOSE,NOPEN,NDUMY,FRACT
C COSMO change A. Klamt
      LOGICAL ISEPS, USEPS , UPDA
      COMMON /ISEPS/  ISEPS, USEPS, UPDA
C end of COSMO change
C***********************************************************************
C
C   COMPFG CALCULATES (A) THE HEAT OF FORMATION OF THE SYSTEM, AND
C                     (B) THE GRADIENTS, IF LGRAD IS .TRUE.
C
C   ON INPUT  XPARAM = ARRAY OF PARAMETERS TO BE USED IN INTERNAL COORDS
C             LGRAD  = .TRUE. IF GRADIENTS ARE NEEDED, .FALSE. OTHERWISE
C             INT    = .TRUE. IF HEAT OF FORMATION IS TO BE CALCULATED
C             FULSCF = .TRUE. IF FULL SCF TO BE DONE, .FALSE. OTHERWISE.
C
C   ON OUTPUT ESCF  = HEAT OF FORMATION.
C             GRAD   = ARRAY OF GRADIENTS, IF LGRAD = .TRUE.
C
C***********************************************************************
      COMMON /KEYWRD/KEYWRD
      CHARACTER*241 KEYWRD
      LOGICAL DEBUG, INT, PRINT, ANALYT, LARGE, USEDCI,
     1FORCE, TIMES, AIDER
      DIMENSION COORD(3,NUMATM), W(N2ELEC), DEGREE(3), XPAREF(MAXPAR)
     1,DELTAP(NMECI**2) ,DELTA(NMECI*MAXORB)
      SAVE DEGREE, PRINT, DEBUG
      EQUIVALENCE (W,WJ)
      DATA ICALCN /0/
C                 MNDO     AM1      PM3      MINDO/
      IF (ICALCN.NE.NUMCAL) THEN
         ICALCN=NUMCAL
         HTYPE(1)=6.1737D0
         HTYPE(2)=3.3191D0
         HTYPE(3)=7.1853D0
         HTYPE(4)=1.7712D0
         LTYPE=0
         DO 30 I=1,NUMAT
            IF(NAT(I).LT.99)THEN
               DO 10 J=1,LTYPE
   10          IF(NAT(I).EQ.MTYPE(J)) GOTO 20
               LTYPE=LTYPE+1
               MTYPE(LTYPE)=NAT(I)
               NZTYPE(NAT(I))=LTYPE
C
C       LTYPE = NUMBER OF TYPES OF REAL ATOM PRESENT
C       MTYPE = TYPES OF REAL ATOMS PRESENT
               J=LTYPE
   20          CONTINUE
            ENDIF
   30    CONTINUE
         AIDER=(INDEX(KEYWRD,'AIDER').NE.0)
         TIMES=(INDEX(KEYWRD,'TIMES').NE.0)
         ANALYT=(INDEX(KEYWRD,'ANALYT').NE.0)
         IF(INT.AND.ANALYT)CALL SETUPG
         DEGREE(1)=1.D0
         IF(INDEX(KEYWRD,' XYZ').NE.0)THEN
            DEGREE(2)=1.D0
         ELSE
            DEGREE(2)=180.D0/3.141592652589D0
         ENDIF
         DEGREE(3)=DEGREE(2)
         USEDCI=(NCLOSE.NE.NOPEN.AND.FRACT.NE.2.D0.AND.FRACT.NE.0.D0
     1         .OR.(INDEX(KEYWRD,'C.I.').NE.0))
         FORCE=(INDEX(KEYWRD,'FORCE').NE.0)
         LARGE=(INDEX(KEYWRD,'LARGE') .NE. 0)
         PRINT=(INDEX(KEYWRD,'COMPFG') .NE. 0)
         DEBUG=(INDEX(KEYWRD,'DEBUG') .NE. 0 .AND. PRINT)
         EMIN=0.D0
         DO 40 I=1,NVAR
   40    XPAREF(I)=XPARAM(I)
      ENDIF
C
C SET UP COORDINATES FOR CURRENT CALCULATION
C
C       PLACE THE NEW VALUES OF THE VARIABLES IN THE ARRAY GEO.
C       MAKE CHANGES IN THE GEOMETRY.
      DO 50 I=1,NVAR
         K=LOC(1,I)
         L=LOC(2,I)
   50 GEO(L,K)=XPARAM(I)
C      IMPOSE THE SYMMETRY CONDITIONS + COMPUTE THE DEPENDENT-PARAMETERS
      IF(NDEP.NE.0) CALL SYMTRY
C      NOW COMPUTE THE ATOMIC COORDINATES.
      IF( DEBUG ) THEN
         IF( LARGE ) THEN
            K=NATOMS
         ELSE
            K=MIN(5,NATOMS)
         ENDIF
         WRITE(6,FMT='('' INTERNAL COORDS'',/100(/,3F12.6))')
     1            ((GEO(J,I)*DEGREE(J),J=1,3),I=1,K)
      END IF
      CALL GMETRY(GEO,COORD)
      IF( DEBUG ) THEN
         IF( LARGE ) THEN
            K=NUMAT
         ELSE
            K=MIN(5,NUMAT)
         ENDIF
         WRITE(6,FMT='('' CARTESIAN COORDS'',/100(/,3F16.9))')
     1            ((COORD(J,I),J=1,3),I=1,K)
      ENDIF
      IF(INT.AND.ANALYT)REWIND 2
C COSMO change A. Klamt
      IF (.NOT. USEPS) THEN
C end of COSMO change
      IF(TIMES)CALL TIMER('BEFORE HCORE')
      IF(INT)CALL HCORE(COORD, H, W, WJ, WK, ENUCLR)
      IF(TIMES)CALL TIMER('AFTER HCORE')
C
C COMPUTE THE HEAT OF FORMATION.
C
      IF(NORBS.GT.0.AND.NELECS.GT.0) THEN
         IF(TIMES)CALL TIMER('BEFORE ITER')
         IF(INT) CALL ITER(H, W, WJ, WK, ELECT, FULSCF,.TRUE.)
         IF(TIMES)CALL TIMER('AFTER ITER')
      ELSE
         ELECT=0.D0
      ENDIF
      ESCF=(ELECT+ENUCLR)*23.061D0+ATHEAT
      IF(ESCF.LT.EMIN.OR.EMIN.EQ.0.D0)EMIN=ESCF
      DO 60 I=1,NNHCO
         CALL DIHED(COORD,NHCO(1,I),NHCO(2,I),NHCO(3,I),NHCO(4,I),ANGLE)
         ESCF=ESCF+HTYPE(ITYPE)*SIN(ANGLE)**2
   60 CONTINUE
C COSMO change A. Klamt 18.7.91
      ENDIF
      IF (ISEPS) THEN
C The following routine constructs the dielectric screening surface
           CALL CONSTS (COORD,.TRUE.)
C The following routine constructs dielectric response matrix CCMAT
        CALL BTOC (COORD)
C A. Klamt 18.7.91
        USEPS = .TRUE.
        IF(TIMES) CALL TIMER('BEFORE HCORE')
        IF(INT) CALL HCORE(COORD, H, W, WJ, WK, ENUCLR)
        IF(TIMES) CALL TIMER('AFTER HCORE')
C
C COMPUTE THE HEAT OF FORMATION.
C
        IF(NORBS.GT.0.AND.NELECS.GT.0) THEN
          IF(TIMES) CALL TIMER('BEFORE ITER')
          IF(INT) CALL ITER(H, W, WJ, WK, ELECT, FULSCF,.TRUE.)
          IF(TIMES) CALL TIMER('AFTER ITER')
        ELSE
          ELECT=0.D0
        ENDIF
        ESCF=(ELECT+ENUCLR)*23.060542301389D0+ATHEAT
        IF(ESCF.LT.EMIN.OR.EMIN.EQ.0.D0) EMIN=ESCF
        DO 61 I=1,NNHCO
         CALL DIHED(COORD,NHCO(1,I),NHCO(2,I),NHCO(3,I),NHCO(4,I),ANGLE)
         ESCF=ESCF+HTYPE(ITYPE)*SIN(ANGLE)**2
   61   CONTINUE
      ENDIF
C end of COSMO change
C       LAURENT MODIFICATION: ADDED LINEAR PENALTY ENERGY
        CALL PENENERGY(COORD,ESCF)
C       END LAURENT


C
C FIND DERIVATIVES IF DESIRED
C
      IF(LGRAD) THEN
         IF(TIMES)CALL TIMER('BEFORE DERIV')
         CALL DERIV(GEO,GRAD)
         IF(TIMES)CALL TIMER('AFTER DERIV')
      ENDIF
      IF(AIDER)THEN
C
C  ADD IN AB INITIO CORRECTION
C
         DO 70 I=1,NVAR
   70    ESCF=ESCF+(XPARAM(I)-XPAREF(I))*AICORR(I)
      ENDIF
      IF(INT.AND.PRINT)
     1WRITE(6,'(/10X,'' HEAT OF FORMATION'',G30.17)')ESCF
      IF(PRINT.AND.LGRAD)
     1   WRITE(6,FMT='('' GRADIENT       '',8F8.2,(/10F8.2))')
     2                (GRAD(I),I=1,NVAR)
C
C REFORM DENSITY MATRIX, IF A C.I. DONE AND EITHER THE LAST SCF OR A
C FORCE CALCULATION
C
      IF(USEDCI.AND. (LAST.EQ.1 .OR. FORCE))
     1CALL MECIP(C,NORBS,DELTAP,DELTA)
      RETURN
      END
