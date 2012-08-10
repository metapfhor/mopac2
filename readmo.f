      SUBROUTINE READMO
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INCLUDE 'SIZES'
C
C MODULE TO READ IN GEOMETRY FILE, OUTPUT IT TO THE USER,
C AND CHECK THE DATA TO SEE IF IT IS REASONABLE.
C EXIT IF NECESSARY.
C
C
C
C  ON EXIT NATOMS    = NUMBER OF ATOMS PLUS DUMMY ATOMS (IF ANY).
C          KEYWRD    = KEYWORDS TO CONTROL CALCULATION
C          KOMENT    = COMMENT CARD
C          TITLE     = TITLE CARD
C          LABELS    = ARRAY OF ATOMIC LABELS INCLUDING DUMMY ATOMS.
C          GEO       = ARRAY OF INTERNAL COORDINATES.
C          LOPT      = FLAGS FOR OPTIMIZATION OF MOLECULE
C          NA        = ARRAY OF LABELS OF ATOMS, BOND LENGTHS.
C          NB        = ARRAY OF LABELS OF ATOMS, BOND ANGLES.
C          NC        = ARRAY OF LABELS OF ATOMS, DIHEDRAL ANGLES.
C          LATOM     = LABEL OF ATOM OF REACTION COORDINATE.
C          LPARAM    = RC: 1 FOR LENGTH, 2 FOR ANGLE, AND 3 FOR DIHEDRAL
C          REACT(200)= REACTION COORDINATE PARAMETERS
C          LOC(1,I)  = LABEL OF ATOM TO BE OPTIMIZED.
C          LOC(2,I)  = 1 FOR LENGTH, 2 FOR ANGLE, AND 3 FOR DIHEDRAL.
C          NVAR      = NUMBER OF PARAMETERS TO BE OPTIMIZED.
C          XPARAM    = STARTING VALUE OF PARAMETERS TO BE OPTIMIZED.
C
************************************************************************
C *** INPUT THE TRIAL GEOMETRY  \IE.  KGEOM=0\
C   LABEL(I) = THE ATOMIC NUMBER OF ATOM\I\.
C            = 99, THEN THE I-TH ATOM IS A DUMMY ATOM USED ONLY TO
C              SIMPLIFY THE DEFINITION OF THE MOLECULAR GEOMETRY.
C   GEO(1,I) = THE INTERNUCLEAR SEPARATION \IN ANGSTROMS\ BETWEEN ATOMS
C              NA(I) AND (I).
C   GEO(2,I) = THE ANGLE NB(I):NA(I):(I) INPUT IN DEGREES; STORED IN
C              RADIANS.
C   GEO(3,I) = THE ANGLE BETWEEN THE VECTORS NC(I):NB(I) AND NA(I):(I)
C              INPUT IN DEGREES - STORED IN RADIANS.
C  LOPT(J,I) = -1 IF GEO(J,I) IS THE REACTION COORDINATE.
C            = +1 IF GEO(J,I) IS A PARAMETER TO BE OPTIMIZED
C            =  0 OTHERWISE.
C *** NOTE:    MUCH OF THIS DATA IS NOT INCLUDED FOR THE FIRST 3 ATOMS.
C     ATOM1  INPUT LABELS(1) ONLY.
C     ATOM2  INPUT LABELS(2) AND GEO(1,2) SEPARATION BETWEEN ATOMS 1+2
C     ATOM3  INPUT LABELS(3), GEO(1,3)    SEPARATION BETWEEN ATOMS 2+3
C              AND GEO(2,3)              ANGLE ATOM1 : ATOM2 : ATOM3
C
************************************************************************
C
      DIMENSION LOPT(3,NUMATM)
      CHARACTER KEYWRD*241, KOMENT*81, TITLE*81, LINE*80, BANNER*80
      CHARACTER KEYS(80)*1, SPACE*1, SPACE2*2, CH*1, CH2*2
      CHARACTER ELEMNT*2, IDATE*24, GETNAM*80, NAME*4
      COMMON /KEYWRD/ KEYWRD
      COMMON /TITLES/ KOMENT,TITLE
      COMMON /GEOVAR/ NVAR, LOC(2,MAXPAR), IDUMY, XPARAM(MAXPAR)
      COMMON /PATH  / LATOM,LPARAM,REACT(200)
      COMMON /MESH  / LATOM1, LPARA1, LATOM2, LPARA2
      COMMON /ELEMTS/ ELEMNT(107)
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     * NMIDLE(NUMATM),NLAST(NUMATM),NORBS,NELECS,NALPHA,NBETA,
     * NCLOSE,NOPEN,NDUMY,FRACT
      COMMON /OKMANY/ ISOK
      PARAMETER (MXDIM=MAXORB+NUMATM)
      COMMON /SYMRES/ TRANS,RTR,SIG,NAME,NAMO(MXDIM),INDX(MXDIM),
     * ISTA(2)
      COMMON /ISTOPE/ AMS(107)
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /NUMCAL/ NUMCAL
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /GEOSYM/ NDEP, LOCPAR(MAXPAR), IDEPFN(MAXPAR),
     1                      LOCDEP(MAXPAR)
C     PATAS
      COMMON /GRIDD/ XP(3),YP(3),ZP(3),DS
      COMMON /CONNOL/ SCALE,SCINCR,DEN,NSURF
      COMMON /MSTSOL/ EPS,DMP,MC,ICOMP,IFIELD,ICAV
      COMMON /MSTSUR/ OMEGA,RD,RET,FRO,DR,NDIV,NESF,ICENT
      COMMON /MSTCAV/ TABS,VMOL,DMOL,TCE,STEN,DSTEN,CMF
      COMMON /POLY/ XE(500),YE(500),ZE(500),RE(500),
     1SSFE(500),
     1IPLOCH(1500),AS(1500),STOT,VOL,NSF,NC1(NUMATM)
      COMMON /FACTOR/ FACTOR
      COMMON /CLASES/ ICLASS(NUMATM)
C       Laurent Modification
      COMMON /AXES / XHAT(3),YHAT(3),ZHAT(3),OFF(3),ATOT

      COMMON /PERMUTE /PR,PRT
            INTEGER PR(NUMATM),PRT(NUMATM)
      DOUBLE PRECISION ATOT(3,3)
      COMMON /SCANR / ISCAN,NSCAN,STEPS,STARTS,LIMS,NDIM,IVALS,REDSCN
      INTEGER ISCAN,NSCAN,IVALS(3*NUMATM),NDIM
      DOUBLE PRECISION STEPS(3*NUMATM),STARTS(3*NUMATM),LIMS(3*NUMATM)
      LOGICAL REDSCN

C       Laurent end
C**********************************************************************
C* SHIHAO'S MODIFICATION START
C* Added:
      COMMON /MOLCONST/ CTYPE,ITORS(4),CVALUE
C* SHIHAO'S MODIFICATION END
C**********************************************************************
C       LAURENT MODIFICATION

C       END LAURENT
C     PATAS

      LOGICAL INT, AIGEO, ISOK
      SAVE SPACE, SPACE2, IREACT, INT
      DIMENSION COORD(3,NUMATM),VALUE(40)
      EQUIVALENCE (KEYS(1),KEYWRD)
      DATA SPACE, SPACE2/' ','  '/
C     PATAS
      DATA ANTOAU,ONE/1.8897626D0,1.0D0/
C     PATAS
      CONVTR=2.D0*ASIN(1.D0)/180.D0
      AIGEO=.FALSE.
   10 CONTINUE
C
C      LAURENT MODIFICATION
      IF(.NOT.REDSCN)THEN
C       END LAURENT
      CALL GETTXT
      IF(INDEX(KEYWRD,'ECHO').NE.0)THEN
         REWIND 5
         ISOK=.FALSE.
         DO 50 I=1,1000
            READ(5,'(A)',END=60)KEYWRD
            DO 20 J=80,2,-1
   20       IF(KEYWRD(J:J).NE.' ')GOTO 30
            J=1
   30       DO 40 K=1,J
   40       IF(ICHAR(KEYWRD(K:K)).LT.32)KEYWRD(K:K)='*'
            WRITE(6,'(1X,A)')KEYWRD(1:J)
   50    CONTINUE
   60    CONTINUE
         REWIND 5
         CALL GETTXT
      ENDIF
      IF(INDEX(KEYWRD,'ECHO').NE.0)WRITE(6,'(''1'')')
      IF(KEYWRD(1:1) .NE. SPACE) THEN
         CH=KEYWRD(1:1)
         KEYWRD(1:1)=SPACE
         DO 70 I=2,239
            CH2=KEYWRD(I:I)
            KEYWRD(I:I)=CH
            CH=CH2
            IF(KEYWRD(I+1:I+2) .EQ. SPACE2) THEN
               KEYWRD(I+1:I+1)=CH
               GOTO 80
            ENDIF
   70    CONTINUE
         CH2=KEYWRD(240:240)
         KEYWRD(240:240)=CH
         KEYWRD(241:241)=CH2
   80    CONTINUE
      ENDIF
      IF(KOMENT(1:1) .NE. SPACE) THEN
         CH=KOMENT(1:1)
         KOMENT(1:1)=SPACE
         DO 90 I=2,79
            CH2=KOMENT(I:I)
            KOMENT(I:I)=CH
            CH=CH2
            IF(KOMENT(I+1:I+2) .EQ. SPACE2) THEN
               KOMENT(I+1:I+1)=CH
               GOTO 100
            ENDIF
   90    CONTINUE
         CH2=KOMENT(80:80)
         KOMENT(80:80)=CH
         KOMENT(81:81)=CH2
  100    CONTINUE
      ENDIF
      IF(TITLE(1:1) .NE. SPACE) THEN
         CH=TITLE(1:1)
         TITLE(1:1)=SPACE
         DO 110 I=2,79
            CH2=TITLE(I:I)
            TITLE(I:I)=CH
            CH=CH2
            IF(TITLE(I+1:I+2) .EQ. SPACE2) THEN
               TITLE(I+1:I+1)=CH
               GOTO 120
            ENDIF
  110    CONTINUE
         CH2=TITLE(80:80)
         TITLE(80:80)=CH
         TITLE(81:81)=CH2
  120    CONTINUE
      ENDIF
      DO 121 I=1,200
  121 REACT(I)=0.D0
      LATOM=0
      LPARAM=0
      IF(INDEX(KEYWRD,'OLDGEO').EQ.0) THEN
         NVAR=0
         NDEP=0
         IF(AIGEO.OR.INDEX(KEYWRD,'AIGIN').NE.0)THEN
            CALL GETGEG(5,LABELS,GEO,     NA,NB,NC,AMS,NATOMS,INT)
         IF(NVAR.EQ.0)THEN
         DO 122 J=1,3
         DO 122 I=1,NATOMS
  122    LOPT(J,I)=0
         ENDIF
         ELSE
            CALL GETGEO(5,LABELS,GEO,LOPT,NA,NB,NC,AMS,NATOMS,INT)
            IF(NATOMS.LT.0)THEN
               REWIND 5
               IF(NUMCAL.NE.1)THEN
                  WRITE(6,'(//,A)')'   GAUSSIAN INPUT REQUIRES STAND-ALO
     1NE JOB'
                  WRITE(6,'(/,A)')'   OR KEYWORD "AIGIN"'
                  STOP
               ENDIF
               AIGEO=.TRUE.
               GOTO 10
            ENDIF
         ENDIF
         IF(NATOMS.EQ.0)STOP
      ELSE
      DEGREE=90.D0/ASIN(1.D0)
      IF(NA(1).EQ.99)THEN
      DO 128 I=1,NATOMS
      DO 128 J=1,3
      LOPT(J,I)=1
  128 COORD(J,I)=GEO(J,I)
      LOPT(1,1)=0
      LOPT(2,1)=0
      LOPT(3,1)=0
      LOPT(2,2)=0
      LOPT(3,2)=0
      LOPT(3,3)=0
      CALL XYZINT(COORD,NATOMS,NA,NB,NC,DEGREE,GEO)
      ELSE 
      DO 130 I=1,NATOMS
      DO 130 J=2,3
 130  GEO(J,I)=GEO(J,I)*DEGREE
      ENDIF
      ENDIF
C      LAURENT MODIFICATION
      ENDIF !!!!REDCSN
C       END LAURENT

C       LAURENT: DONE READING IN AT THIS POINT

      IF(INDEX(KEYWRD,'FORCE').NE.0 .AND. LABELS(NATOMS).EQ.107) THEN
      DO 131 I=1,NA(NATOMS)
      IF(LABELS(I).EQ.99)THEN
      WRITE(6,'(A)')' NO DUMMY ATOMS ALLOWED BEFORE TRANSLATION'
      WRITE(6,'(A)')' ATOM IN A FORCE CALCULATION'
      STOP
      ENDIF
  131 CONTINUE
      ENDIF
C
C
C OUTPUT FILE TO UNIT 6
C
C    WRITE HEADER
      IDATE=' '
      CALL fdate(IDATE)
C       LAURENT MODIFICATION
      IF(.NOT.REDSCN)THEN
C       END LAURENT
      WRITE(6,'(1X,15(''*****''),''****'')')
C
C     CHANGE THE FOLLOWING LINE TO SUIT LOCAL ENVIRONMENT, IF DESIRED
C
      BANNER=' ** MOPAC (PUBLIC DOMAIN) FOR DEVELOPMENT USE '//
     1'ONLY.  NOT FOR PRODUCTION WORK  **'
      BANNER=' **                      MOPAC (PUBLIC DOMAIN)'//
     1'                                **'
      WRITE(6,'(A)')BANNER
C
C    THE BANNER DOES NOT APPEAR ANYWHERE ELSE.
C
      WRITE(6,'(1X,79(''*''))')
      LINE='   MNDO'
      IF(INDEX(KEYWRD,'MINDO') .NE. 0) LINE='MINDO/3'
      IF(INDEX(KEYWRD,'AM1')   .NE. 0) LINE='    AM1'
      IF(INDEX(KEYWRD,'PM3')   .NE. 0) LINE='    PM3'
C* SHIHAO MODIFICATION START ****************************************
C* Added:
      IF(INDEX(KEYWRD,'PDG')   .NE. 0) LINE='    PDG'
      IF(INDEX(KEYWRD,'MDG')   .NE. 0) LINE='    MDG'
C* SHIHAO MODIFICATION END ******************************************
      WRITE(6,'(/29X,A,'' CALCULATION RESULTS'',28X,///1X,15(''*****'')
     1,''****'' )')LINE(:7)
      WRITE(6,'('' *'',10X,''MOPAC:  VERSION '',F5.2,
     115X,''CALC''''D. '',A)') VERSON, IDATE
C
C CONVERT ANGLES TO RADIANS
!      DO 140 J=2,3
!C$DOIT VBEST
!         DO 140 I=1,NATOMS
!            GEO(J,I) = GEO(J,I) * CONVTR
!  140 CONTINUE
C
C CHECK DATA
C
C       LAURENT MODIFICATION
      ENDIF
C       END LAURENT
      NA(1)=0
      NB(1)=0
      NC(1)=0
      DO 150 I=1,NATOMS
         IF (LABELS(I) .LE. 0 ) THEN
            WRITE(6,'('' ATOMIC NUMBER OF '',I3,'' ?'')') LABELS(I)
            IF(I.EQ.1) THEN
               WRITE(6,'(A)')' THIS WAS THE FIRST ATOM'
            ELSE
               WRITE(6,'(A)')'    GEOMETRY UP TO, BUT NOT INCLUDING, THE
     1'//' FAULTY ATOM'
               NATOMS=I-1
               CALL GEOUT(6)
            ENDIF
            STOP
         ENDIF
         IF (  NA(I).GE.I.OR. NB(I).GE.I.OR. NC(I).GE.I
     1  .OR. (NA(I).EQ.NB(I))   .AND. I.GT.1
     2  .OR. (NA(I).EQ.NC(I).OR.NB(I).EQ.NC(I))  .AND. I.GT.2
     3  .OR.  NA(I)*NB(I)*NC(I).EQ.0  .AND. I.GT.3) THEN
            WRITE(6,'('' ATOM NUMBER '',I3,'' IS ILL-DEFINED'')') I
            IF(I.EQ.1)STOP
            WRITE(6,'(/,''  GEOMETRY READ IN'',/)')
            CALL GEOUT(6)
            STOP
         ENDIF
  150 CONTINUE
C
C WRITE KEYWORDS BACK TO USER AS FEEDBACK
C       LAURENT MODIFICATION
      IF(.NOT.REDSCN)THEN
C       END LAURENT
      CALL WRTKEY(KEYWRD)
      WRITE(6,'(1X,14(''*****''),''*'',I3.3,''BY'',I3.3)')MAXHEV,MAXLIT
C       LAURENT MODIFICATION
      ELSE
         CALL RDKEY(KEYWRD)
      ENDIF
C       END LAURENT
C
C FILL IN GEO MATRIX IF NEEDED
      IF(INDEX(KEYWRD,'OLDGEO').EQ.0.AND.INDEX(KEYWRD,'SYM') .NE. 0
     1.AND. NDEP.EQ.0) CALL GETSYM
      IF(NDEP.NE.0) CALL SYMTRY
C
C INITIALIZE FLAGS FOR OPTIMIZE AND PATH
      IFLAG = 0
      LATOM = 0
      NUMAT=0
      IF(NVAR.NE.0)THEN
         NUMAT=NATOMS
      ELSE
         DO 180 I=1,NATOMS
            IF(LABELS(I).NE.99.AND.LABELS(I).NE.107)NUMAT=NUMAT+1
            DO 180 J=1,3
               IF (LOPT(J,I) ) 160, 180, 170
C    FLAG FOR PATH
  160          CONVRT=1.D0
               IF ( IFLAG .NE. 0 ) THEN
                  IF(INDEX(KEYWRD,'STEP1').NE.0)THEN
                     LPARA1=LPARAM
                     LATOM1=LATOM
                     LPARA2=J
                     LATOM2=I
                     LATOM=0
                     IFLAG=0
                     GOTO 180
                  ELSE
                     WRITE(6,'('' ONLY ONE REACTION COORDINATE PERMITTED
     1'')')
                     STOP
                  ENDIF
               ENDIF
               LATOM  = I
               LPARAM = J
               IF(J.GT.1) CONVRT=0.01745329252D00
               REACT(1)  = GEO(J,I)
               IREACT=1
               IFLAG = 1
               GO TO 180
C    FLAG FOR OPTIMIZE
  170          NVAR = NVAR + 1
               LOC(1,NVAR) = I
               LOC(2,NVAR) = J
               XPARAM(NVAR)   = GEO(J,I)
  180    CONTINUE
      ENDIF

C READ IN PATH VALUES
      IF(IFLAG.EQ.0) GO TO 221
      IF(INDEX(KEYWRD,'NLLSQ').NE.0)THEN
         WRITE(6,'(A)')' NLLSQ USED WITH REACTION PATH; '//
     1'THIS OPTION IS NOT ALLOWED'
         STOP
      ENDIF
      IF(INDEX(KEYWRD,'SIGMA').NE.0)THEN
         WRITE(6,'(A)')' SIGMA USED WITH REACTION PATH; '//
     1'THIS OPTION IS NOT ALLOWED'
         STOP
      ENDIF
      IF(INDEX(KEYWRD,'STEP')+INDEX(KEYWRD,'POINTS').NE.0)THEN
         STEP=READA(KEYWRD,INDEX(KEYWRD,'STEP=')+5)
         NPTS=READA(KEYWRD,INDEX(KEYWRD,'POINT=')+6)
         IF(NPTS.GT.200)THEN
            WRITE(6,'(///,''    ONLY TWO HUNDRED POINTS ALLOWED IN REACT
     1'',''ION COORDINATE'')')
            STOP
         ENDIF
         IF(LPARAM.EQ.1.AND.STEP.LE.0)THEN
            WRITE(6,'(///,''    STEP FOR BOND LENGTH SHOULD BE SET POSIT
     1IVE '',''TO PREVENT TWO ATOMS COLLAPSE'')')
            STOP
         ENDIF
         GO TO 221
      ENDIF
  190 READ(5,'(A)',END=210) LINE
      CALL NUCHAR(LINE,VALUE,NREACT)
      IF(NREACT.EQ.0)GOTO 210
      DO 200 I=1,NREACT
         IJ=IREACT+I
         IF(IJ.GT.200)THEN
            WRITE(6,'(///,''    ONLY TWO HUNDRED POINTS ALLOWED IN REACT
     1ION'','' COORDINATE'')')
            STOP
         ENDIF
         REACT(IJ)=VALUE(I)*CONVRT
         IF(ABS(REACT(IJ)-REACT(IJ-1)).LT.1.D-5)THEN
            DUM1 = REACT(IJ)/CONVRT
            DUM2 = REACT(IJ-1)/CONVRT
            WRITE(6,'(///,'' TWO ADJACENT POINTS ARE IDENTICAL:  '',
     1 F7.3,2X,F7.3,/,'' THIS IS NOT ALLOWED IN A PATH CALCULATION'')')
     2 DUM1,DUM2
            STOP
         ENDIF
  200 CONTINUE
      IREACT=IREACT+NREACT
      GO TO 190
  210 CONTINUE
      DEGREE=1.D0
      IF(LPARAM.GT.1)DEGREE=90.D0/ASIN(1.D0)
      IF(IREACT.LE.1) THEN
         WRITE(6,'(//10X,'' NO POINTS SUPPLIED FOR REACTION PATH'')')
         WRITE(6,'(//10X,'' GEOMETRY AS READ IN IS AS FOLLOWS'')')
         CALL GEOUT(1)
         STOP
      ELSE
         WRITE(6,'(//10X,'' POINTS ON REACTION COORDINATE'')')
         WRITE(6,'(10X,8F8.2)')(REACT(I)*DEGREE,I=1,IREACT)
      ENDIF
      IEND=IREACT+1
      REACT(IEND)=-1.D12
C     PATAS
C
C     READ IN CUBIC GRID OR CONNOLLY SURFACE FOR M.E.P. COMPUTATION
C
  221 IF (INDEX(KEYWRD,'MEP=').NE.0) THEN
      I=(READA(KEYWRD,INDEX(KEYWRD,'MEP=')))
        IF (I.EQ.1) THEN
C     FILL COMMON /GRIDD/
        DO 222 J=1,3
           READ(5,'(A)') LINE
           CALL NUCHAR (LINE,VALUE,NVALUE)
           XP(J)=VALUE(1)*ANTOAU
           YP(J)=VALUE(2)*ANTOAU
           ZP(J)=VALUE(3)*ANTOAU
  222   CONTINUE
        READ(5,'(A)') LINE
        CALL NUCHAR (LINE,VALUE,NVALUE)
        DS=VALUE(1)*ANTOAU
        ELSE
C     FILL COMMON /CONNOL/
        READ(5,'(A)') LINE
        CALL NUCHAR (LINE,VALUE,NVALUE)
        SCALE=VALUE(1)*ONE
        SCINCR=VALUE(2)*ONE
        DEN=VALUE(3)*ONE
        NSURF=VALUE(4)
        ENDIF
      ENDIF
C
C     READ IN DATA FOR MST SOLVATION MODEL
C
      IF (INDEX(KEYWRD,'TOM').NE.0) THEN
C
C     FILL COMMON /MSTSOL/
C
      READ(5,'(A)') LINE
      CALL NUCHAR (LINE,VALUE,NVALUE)
      EPS=VALUE(1)*ONE
      DMP=VALUE(2)
      MC=VALUE(3)
      ICOMP=VALUE(4)
      IFIELD=VALUE(5)
      ICAV=VALUE(6)
C
C     FILL COMMON /MSTCAV/
C
      READ(5,'(A)') LINE
      CALL NUCHAR (LINE,VALUE,NVALUE)
      IF (ICAV.EQ.0) GO TO 303
      IF (ICAV.EQ.1.OR.ICAV.EQ.3) THEN
      TABS=VALUE(1)*ONE
      VMOL=VALUE(2)*ONE
      DMOL=VALUE(3)*ONE
      TCE=VALUE(4)*ONE
      ENDIF
      IF (ICAV.EQ.3) THEN
      STEN=VALUE(5)*ONE
      DSTEN=VALUE(6)*ONE
      CMF=VALUE(7)*ONE
      ENDIF
      IF (ICAV.EQ.2) THEN
      STEN=VALUE(1)*ONE
      DSTEN=VALUE(2)*ONE
      CMF=VALUE(3)*ONE
      ENDIF
C     FILL COMMON /MSTSUR/
  303 READ(5,'(A)') LINE
      CALL NUCHAR (LINE,VALUE,NVALUE)
      OMEGA=VALUE(1)*ONE
      RD=VALUE(2)*ONE
      RET=VALUE(3)*ONE
      FRO=VALUE(4)*ONE
      DR=VALUE(5)*ONE
      NDIV=VALUE(6)
      NESF=VALUE(7)
      ICENT=VALUE(8)
C
C     FILL COMMON /POLY/
C
      DO 300 I=1,NESF
      READ(5,'(A)') LINE
      CALL NUCHAR (LINE,VALUE,NVALUE)
      IF (ICENT.EQ.0) THEN
      XE(I)=VALUE(1)*ONE
      YE(I)=VALUE(2)*ONE
      ZE(I)=VALUE(3)*ONE
      RE(I)=VALUE(4)*ONE
      ICLASS(I)=VALUE(5)
      ELSE
      NC1(I)=VALUE(1)
      RE(I)=VALUE(2)*ONE
      ICLASS(I)=VALUE(3)
      ENDIF
  300 CONTINUE
C
C     FILL COMMON /FACTOR/
C
      READ(5,'(A)') LINE
      CALL NUCHAR (LINE,VALUE,NVALUE)
      FACTOR=VALUE(1)*ONE
      ENDIF

C* SHIHAO MODIFICATION START ****************************************
C* Added:
      IF(INDEX(KEYWRD,'CONST').NE.0)THEN
        READ(5,'(A)') LINE
        PI=3.1415926536D0
        CALL NUCHAR (LINE,VALUE,NVALUE)
        IF(NVALUE.EQ.5) THEN
          ITORS(1)=VALUE(1)
          ITORS(2)=VALUE(2)
          ITORS(3)=VALUE(3)
          ITORS(4)=VALUE(4)
          CVALUE=VALUE(5)/180.0*PI
          CTYPE=3
        ENDIF
      ENDIF
C* SHIHAO MODIFICATION END ******************************************
C     PATAS
C
C OUTPUT GEOMETRY AS FEEDBACK
C
  220 CALL WRTTXT(6)
      IF(INDEX(KEYWRD,'NOLOG').EQ.0)THEN
         OPEN(UNIT=11, FORM='FORMATTED', STATUS='UNKNOWN',
     +FILE=GETNAM('FOR011'))
         CALL WRTTXT(11)
      ENDIF
      CALL GEOUT(1)
      CALL GMETRY(GEO,COORD)
C
C  IF A POLYMER, EXPAND TO MERS
C
      IF(INDEX(KEYWRD,' MERS').NE.0)CALL MAKPOL(COORD)
      IF (INDEX(KEYWRD,'NOXYZ') .EQ. 0) THEN
         IF(INDEX(KEYWRD,'0SCF').NE.0)THEN
C
C  WRITE OUT CARTESIAN COORDINATES FOR USE AS A DATA-SET
C
            WRITE(6,'(A)')'   GEOMETRY IN CARTESIAN COORDINATE FORMAT'
            CALL WRTTXT(6)
            J=0
            DO 230 I=1,NATOMS
               IF(LABELS(I).NE.99)THEN
                  J=J+1
                  WRITE(6,'(2X,A,3(F19.13,I3))')
     1    ELEMNT(LABELS(I)),(COORD(K,J),1,K=1,3)
               ENDIF
  230       CONTINUE
         ELSE
            WRITE(6,'(//10X,''CARTESIAN COORDINATES '',/)')
            WRITE(6,'(4X,''NO.'',7X,''ATOM'',15X,''X'',
     1  15X,''Y'',15X,''Z'',/)')
            L=0
            DO 240 I=1,NATOMS
               IF(LABELS(PRT(I)) .EQ. 99.OR.LABELS(PRT(I)).EQ.107)THEN
                GOTO 240
               ENDIF
               L=L+1
               WRITE(6,'(I6,8X,A2,4X,3F16.10)')
C       Laurent Modification: Added coordinate backtransform
C        1  L,ELEMNT(LABELS(I)),(COORD(J,L),J=1,3)
     1  L,ELEMNT(LABELS(PRT(I))),ATOT(1,1)*COORD(1,PRT(I))+
     2ATOT(1,2)*COORD(2,PRT(I))+ATOT(1,3)*COORD(3,PRT(I))+OFF(1)
     3,ATOT(2,1)*COORD(1,PRT(I))+ ATOT(2,2)*COORD(2,PRT(I))+
     4 ATOT(2,3)*COORD(3,PRT(I))+OFF(2),ATOT(3,1)*COORD(1,PRT(I))+
     5ATOT(3,2)*COORD(2,PRT(I))+ATOT(3,3)*COORD(3,PRT(I))+OFF(3)
C       Laurent End



  240       CONTINUE
         ENDIF
      ENDIF
      IF(.NOT.REDSCN)THEN
      CALL SYMTRZ(COORD,C,NORBS,NORBS,.FALSE.,.FALSE.)

      WRITE(6,'(//''     MOLECULAR POINT GROUP   :   '',A4)') NAME
      ENDIF
      IF(   INDEX(KEYWRD,' XYZ') .NE. 0 )THEN
         IF( NVAR .NE. 0 .AND.
     1 INT.AND.(NDEP .NE. 0 .OR.  NVAR.LT.3*NUMAT-6)) THEN
            IF(NDEP.NE.0)
     1WRITE(6,'(//10X,'' INTERNAL COORDINATES READ IN, AND SYMMETRY''
     2,/10X,'' SPECIFIED, BUT CALCULATION TO BE RUN IN CARTESIAN ''
     3,''COORDINATES'')')
            IF(NVAR.LT.3*NUMAT-6)
     1WRITE(6,'(//10X,'' INTERNAL COORDINATES READ IN, AND'',
     2'' CALCULATION '',/10X,''TO BE RUN IN CARTESIAN COORDINATES, '',
     3/10X,''BUT NOT ALL COORDINATES MARKED FOR OPTIMISATION'')')
            WRITE(6,'(//10X,'' THIS INVOLVES A LOGICALLLY ABSURD CHOICE'
     1',/10X,'' SO THE CALCULATION IS TERMINATED AT THIS POINT'')')
            STOP
         ENDIF
         SUMX=0.D0
         SUMY=0.D0
         SUMZ=0.D0
         DO 250 J=1,NUMAT
            SUMX=SUMX+COORD(1,J)
            SUMY=SUMY+COORD(2,J)
  250    SUMZ=SUMZ+COORD(3,J)
         SUMX=SUMX/NUMAT
         SUMY=SUMY/NUMAT
         SUMZ=SUMZ/NUMAT
C       Laurent Modification: this recentering needs to be taken into account
        OFF(1)=OFF(1)+ATOT(1,1)*SUMX+ATOT(1,2)*SUMY+ATOT(1,3)*SUMZ
        OFF(2)=OFF(2)+ATOT(2,1)*SUMX+ATOT(2,2)*SUMY+ATOT(2,3)*SUMZ
        OFF(3)=OFF(3)+ATOT(3,1)*SUMX+ATOT(3,2)*SUMY+ATOT(3,3)*SUMZ
C       Laurent end
         DO 260 J=1,NUMAT
            GEO(1,J)=COORD(1,J)-SUMX
            GEO(2,J)=COORD(2,J)-SUMY
  260    GEO(3,J)=COORD(3,J)-SUMZ
         NA(1)=99
         J=0
         NVAR=0
         DO 280 I=1,NATOMS
            IF(LABELS(I).NE.99)THEN
               J=J+1
               IF(J.EQ.1)THEN
                  K=0
               ELSEIF(J.LT.4)THEN
                  K=MIN(3,I-1)
               ELSE
                  K=3
               ENDIF
               DO 270 L=1,K
                  NVAR=NVAR+1
                  LOC(1,NVAR)=J
                  LOC(2,NVAR)=L
  270          XPARAM(NVAR)=GEO(L,J)
               LABELS(J)=LABELS(I)
            ENDIF
  280    CONTINUE
         NATOMS=NUMAT
      ELSE
         IF(NVAR.EQ.0) RETURN
         IF(.NOT. INT.AND.(NDEP .NE. 0 .OR.NVAR.LT.3*NUMAT-6)
     1    .AND.INDEX(KEYWRD,'ALTCON').EQ.0) THEN
            IF(NDEP.NE.0)
     1WRITE(6,'(//10X,'' CARTESIAN COORDINATES READ IN, AND SYMMETRY''
     2,/10X,'' SPECIFIED, BUT CALCULATION TO BE RUN IN INTERNAL ''
     3,''COORDINATES'')')

            IF(NVAR.LT.3*NUMAT-6)
     1WRITE(6,'(//10X,'' CARTESIAN COORDINATES READ IN, AND'',
     2'' CALCULATION '',/10X,''TO BE RUN IN INTERNAL COORDINATES, '',
     3/10X,''BUT NOT ALL COORDINATES MARKED FOR OPTIMISATION'')')
            WRITE(6,'(//10X,''MOPAC, BY DEFAULT, USES INTERNAL COORDINAT
     1ES'',/10X,''TO SPECIFY CARTESIAN COORDINATES USE KEY-WORD :XYZ:'')
     2')
            WRITE(6,'(10X,''YOUR CURRENT CHOICE OF KEY-WORDS INVOLVES''
     1,'' A LOGICALLLY'',/10X,''ABSURD CHOICE SO THE CALCULATION IS'',
     2'' TERMINATED AT THIS POINT'')')
            STOP
         ENDIF
      ENDIF
      RETURN
      END
