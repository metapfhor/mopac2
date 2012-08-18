      SUBROUTINE XYZINT(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      DIMENSION XYZ(3,NUMATM), NA(NUMATM), NB(NUMATM), NC(NUMATM),
     1 GEO(3,NUMATM)
***********************************************************************
*
* XYZINT WORKS OUT THE INTERNAL COORDINATES OF A MOLECULE.
*        THE "RULES" FOR THE CONNECTIVITY ARE AS FOLLOWS:
*        ATOM I IS DEFINED AS BEING AT A DISTANCE FROM THE NEAREST
*        ATOM J, ATOM J ALREADY HAVING BEEN DEFINED.
*        ATOM I MAKES AN ANGLE WITH ATOM J AND THE ATOM K, WHICH HAS
*        ALREADY BEEN DEFINED, AND IS THE NEAREST ATOM TO J
*        ATOM I MAKES A DIHEDRAL ANGLE WITH ATOMS J, K, AND L. L HAVING
*        BEEN DEFINED AND IS THE NEAREST ATOM TO K, AND J, K AND L
*        HAVE A CONTAINED ANGLE IN THE RANGE 15 TO 165 DEGREES,
*        IF POSSIBLE.
*
*        IF(NA(2).EQ.-1 OR -2 THEN THE ORIGINAL CONNECTIVITY IS USED.
*
*        NOTE THAT GEO AND XYZ MUST NOT BE THE SAME IN THE CALL.
*
*   ON INPUT XYZ    = CARTESIAN ARRAY OF NUMAT ATOMS
*            DEGREE = 1 IF ANGLES ARE TO BE IN RADIANS
*            DEGREE = 57.29578 IF ANGLES ARE TO BE IN DEGREES
*
***********************************************************************
      COMMON /GEOOK/ IGEOOK
      COMMON /NUMCAL/ NUMCAL
C             Laurent Modification: added
      COMMON /AXES / XHAT(3),YHAT(3),ZHAT(3),OFF(3),ATOT
      COMMON /ALTCON / TRANSL, BANGL, DIHDR, ATOMS, ICONXN, APPLIED,
     1                  VALS, NVALS,NTRANS,NANGL,NDIHD,NATOM
      INTEGER TRANSL(5*NUMATM),NTRANS,BANGL(6*NUMATM),NANGL,
     1         DIHDR(7*NUMATM),NDIHD,ATOMS(NUMATM*NUMATM),
     2         NATOM,ICONXN(6,NUMATM),NVALS
      LOGICAL APPLIED
      DOUBLE PRECISION VALS(3*NUMATM)
      COMMON /LIN / PI,PJ,PK
      COMMON /PERMUTE /PR,PRT
      INTEGER PR(NUMATM),PRT(NUMATM),PI,PJ,PK
      LOGICAL PRMTD,FIRST
      DOUBLE PRECISION ATOT(3,3)
      DOUBlE PRECISION DX, DY, DZ
     1  XYZINIT(3,NUMATM)

      COMMON /KEYWRD/ KEYWRD
      CHARACTER KEYWRD*241
      SAVE PRMTD,FIRST
      DATA PRMTD,FIRST /.FALSE.,.TRUE./
C       Laurent End
      DATA ICALCN/0/
      IGEOOK=99

C       Laurent Modification: Recenter on the first atom


      IF(INDEX(KEYWRD,'ALTCON').NE.0.AND..NOT.APPLIED
     1 .AND.(NATOM.NE.0.OR.PI.NE.0))THEN
        IF(.NOT.PRMTD)THEN
            CALL PERATMS(XYZ)
            PRMTD=.TRUE.
        ENDIF
        CALL AALTCON(XYZ)
        APPLIED=.TRUE.
      ENDIF

      IF(FIRST)THEN
      DX=XYZ(1,1)
      DY=XYZ(2,1)
      DZ=XYZ(3,1)
      OFF(1)=OFF(1)+ATOT(1,1)*DX+ATOT(1,2)*DY+ATOT(1,3)*DZ
      OFF(2)=OFF(2)+ATOT(2,1)*DX+ATOT(2,2)*DY+ATOT(2,3)*DZ
      OFF(3)=OFF(3)+ATOT(3,1)*DX+ATOT(3,2)*DY+ATOT(3,3)*DZ
!      OFF(1)=OFF(1)+ATOT(1,1)*DX+ATOT(2,1)*DY+ATOT(3,1)*DZ
!      OFF(2)=OFF(2)+ATOT(1,2)*DX+ATOT(2,2)*DY+ATOT(3,2)*DZ
!      OFF(3)=OFF(3)+ATOT(1,3)*DX+ATOT(2,3)*DY+ATOT(3,3)*DZ




      DO 40 I=1,NUMAT
        XYZ(1,I)=XYZ(1,I)-DX
        XYZ(2,I)=XYZ(2,I)-DY
        XYZ(3,I)=XYZ(3,I)-DZ
   40 CONTINUE
      FIRST=.FALSE.
      ENDIF
C       Laurent End
      IF(.NOT.(ICALCN.NE.NUMCAL).AND.NA(2).EQ.-1 .OR. NA(2).EQ.-2)THEN
         NA(2)=1
         DO 10 I=2,NUMAT
            J=NA(I)
            IF(I.GT.3)CALL DIHED(XYZ,I,J,NB(I),NC(I),GEO(3,I))
            IF(I.GT.2)CALL BANGLE(XYZ,I,J,NB(I),GEO(2,I))
            GEO(1,I)=SQRT((XYZ(1,I)-XYZ(1,J))**2+
     1          (XYZ(2,I)-XYZ(2,J))**2+
     2          (XYZ(3,I)-XYZ(3,J))**2)
   10    CONTINUE
      ELSE
         K=1
         IF(NA(2).EQ.-1)ICALCN=NUMCAL
         DO 30 I=2,NUMAT
         IF(ICONXN(1,I).EQ.0)THEN
            NA(I)=1
            NB(I)=2
            NC(I)=3
            IM1=I-1
            IF(IM1.EQ.0)GOTO 30
            SUM=1.D30
            DO 20 J=1,IM1
               R=(XYZ(1,I)-XYZ(1,J))**2+
     1          (XYZ(2,I)-XYZ(2,J))**2+
     2          (XYZ(3,I)-XYZ(3,J))**2
               IF(R.LT.SUM.AND.NA(J).NE.J.AND.NB(J).NE.J) THEN
                  SUM=R
                  K=J
               ENDIF
   20       CONTINUE
C
C   ATOM I IS NEAREST TO ATOM K
C

        ELSE
        K=ICONXN(1,I)
        ENDIF
            NA(I)=K
            IF(I.GT.2)THEN
                IF(ICONXN(2,I).EQ.0)THEN
                        IF(NA(I).EQ.I-2)THEN
                            NB(I)=I-1
                        ELSE
                            NB(I)=I-2
                        ENDIF
                ELSE
                        NB(I)=ICONXN(2,I)
                ENDIF
            ENDIF
            IF(I.GT.3)THEN
                IF(ICONXN(3,I).EQ.0)THEN
                            IF(NA(I).EQ.I-3.OR.NB(I).EQ.I-3)THEN
                                IF(NA(I).EQ.I-3)THEN
                                    IF(NB(I).EQ.I-2)THEN
                                          NC(I)=I-1
                                    ELSE
                                        NC(I)=I-2
                                    ENDIF
                                ELSE
                                    IF(NA(I).EQ.I-2)THEN
                                         NC(I)=I-1
                                    ELSE
                                         NC(I)=I-2
                                    ENDIF
                                ENDIF
                            ELSE
                                NC(I)=I-3
                            ENDIF
                ELSE
                    NC(I)=ICONXN(3,I)
                ENDIF
            ENDIF
C
C   FIND ANY ATOM TO RELATE TO NA(I)
C
   30    CONTINUE
      ENDIF
      NA(1)=0
      NB(1)=0
      NC(1)=0
      NB(2)=0
      NC(2)=0
      NC(3)=0

      CALL XYZGEO(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
      RETURN
      END
      SUBROUTINE XYZGEO(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      DIMENSION XYZ(3,NUMATM), NA(NUMATM), NB(NUMATM), NC(NUMATM),
     1  GEO(3,NUMATM)
***********************************************************************
*
*   XYZGEO CONVERTS COORDINATES FROM CARTESIAN TO INTERNAL.
*
*     ON INPUT XYZ  = ARRAY OF CARTESIAN COORDINATES
*              NUMAT= NUMBER OF ATOMS
*              NA   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY DISTANCE
*              NB   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY ANGLE
*              NC   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY DIHEDRAL
*
*    ON OUTPUT GEO  = INTERNAL COORDINATES IN ANGSTROMS, RADIANS,
*                     AND RADIANS
*
***********************************************************************
C       Laurent
      COMMON /AXES / XHAT(3),YHAT(3),ZHAT(3),OFF(3),ATOT
      COMMON /ALTCON / TRANSL, BANGL, DIHDR, ATOMS, ICONXN, APPLIED,
     1                  VALS, NVALS,NTRANS,NANGL,NDIHD,NATOM
      INTEGER TRANSL(5*NUMATM),NTRANS,BANGL(6*NUMATM),NANGL,
     1         DIHDR(7*NUMATM),NDIHD,ATOMS(NUMATM*NUMATM),
     2         NATOM,ICONXN(6,NUMATM),NVALS
      LOGICAL APPLIED
      DOUBLE PRECISION VALS(3*NUMATM)
      COMMON /PERMUTE /PR,PRT
      INTEGER PR(NUMATM),PRT(NUMATM)
      COMMON /GENRAL/ COORD(3,NUMATM), COLD(3,NUMATM*3), GOLD(MAXPAR),
     1 XPARAM(MAXPAR)
      DOUBLE PRECISION XHP(3),YHP(3),ZHP(3),XHT(3),YHT(3),ZHT(3),LNP
     1 ,ATMP(3,3),XY,XZ,YZ,ZZ,ATOT(3,3)
      COMMON /KEYWRD/ KEYWRD
      CHARACTER KEYWRD*241
C       /Laurent

C       Laurent: This loop is all about defining NC(I), in order to give good bond angles
C                It must be bypassed in order to apply constraints on the dihedral
      DO 30 I=2,NUMAT
         J=NA(I)
         K=NB(I)
         L=NC(I)
         IF(I.LT.3) GOTO 30
         II=I
         CALL BANGLE(XYZ,II,J,K,GEO(2,I))
         GEO(2,I)=GEO(2,I)*DEGREE
         IF(I.LT.4) GOTO 30
C
C   MAKE SURE DIHEDRAL IS MEANINGLFUL
C

         CALL BANGLE(XYZ,J,K,L,ANGL)
         TOL=0.2617994D0
         IF((ANGL.GT.3.1415926D0-TOL.OR.ANGL.LT.TOL).AND.ICONXN(3,I)
     1    .EQ.0)THEN
C
C  ANGLE IS UNSATISFACTORY, LET'S SEARCH FOR ANOTHER ATOM FOR
C  DEFINING THE DIHEDRAL.
   10       SUM=100.D0
            DO 20 I1=1,II-1
               R=(XYZ(1,I1)-XYZ(1,K))**2+
     1          (XYZ(2,I1)-XYZ(2,K))**2+
     2          (XYZ(3,I1)-XYZ(3,K))**2
               IF(R.LT.SUM.AND.I1.NE.J.AND.I1.NE.K) THEN
                  CALL BANGLE(XYZ,J,K,I1,ANGL)
                  IF(ANGL.LT.3.1415926D0-TOL.AND.ANGL.GT.TOL)THEN
                     SUM=R
                     L=I1
                     NC(II)=L
                  ENDIF
               ENDIF
   20       CONTINUE
            IF(SUM.GT.99.D0.AND.TOL.GT.0.1D0)THEN
C
C ANYTHING WITHIN 5 DEGREES?
C
               TOL=0.087266D0
               GOTO 10
            ENDIF
         ENDIF
         CALL DIHED(XYZ,I,NA(I),NB(I),NC(I),GEO(3,I))
         GEO(3,I)=GEO(3,I)*DEGREE
   30 GEO(1,I)= SQRT((XYZ(1,I)-XYZ(1,NA(I)))**2+
     1                   (XYZ(2,I)-XYZ(2,NA(I)))**2+
     2                   (XYZ(3,I)-XYZ(3,NA(I)))**2)

C       Geometric Constraints should be imposed here



C       Laurent: calculate the input XYZ components for our basis
C       This may need to be reworked for cases where the first two molecules are colinear
C       Translation (should be done at the end in input coords)



C       New Axes

      XHP=XYZ(:,2)
      LNP = SQRT((XHP(1))**2+(XHP(2))**2+(XHP(3))**2)
      XHP(1)=XHP(1)/LNP
      XHP(2)=XHP(2)/LNP
      XHP(3)=XHP(3)/LNP

C       remove component along the new xhat from the coords of the third to get the second axis
      YHP=XYZ(:,3)
      LNP = dot_product(YHP,XHP)
      YHP(1)=YHP(1)-LNP*XHP(1)
      YHP(2)=YHP(2)-LNP*XHP(2)
      YHP(3)=YHP(3)-LNP*XHP(3)
C       make it unitary
      LNP = SQRT((YHP(1))**2+(YHP(2))**2+(YHP(3))**2)
      YHP(1)=YHP(1)/LNP
      YHP(2)=YHP(2)/LNP
      YHP(3)=YHP(3)/LNP
      XY=dot_product(XHP,YHP)

C       find the third axis from the fact that k=ixj
      ZHP(1)=XHP(2)*YHP(3)-XHP(3)*YHP(2)
      ZHP(2)=XHP(3)*YHP(1)-XHP(1)*YHP(3)
      ZHP(3)=XHP(1)*YHP(2)-XHP(2)*YHP(1)
      XZ=dot_product(XHP,ZHP)
      YZ=dot_product(YHP,ZHP)
      ZZ=dot_product(ZHP,ZHP)
C       Update our axes


      ATMP(:,1)=XHP
      ATMP(:,2)=YHP
      ATMP(:,3)=ZHP

      ATOT = MATMUL(ATMP,ATOT)

      IF(.NOT.APPLIED)THEN
        DO I=1,NUMATM
            GEO(2,I)=GEO(2,I)/DEGREE
            GEO(3,I)=GEO(3,I)/DEGREE
        END DO

        CALL GMETRY(GEO,COORD)

        DO I=1,NUMATM
            GEO(2,I)=GEO(2,I)*DEGREE
            GEO(3,I)=GEO(3,I)*DEGREE
        END DO
      ENDIF
C       Laurent End

      GEO(1,1)=0.D0
      GEO(2,1)=0.D0
      GEO(3,1)=0.D0
      GEO(2,2)=0.D0
      GEO(3,2)=0.D0
      GEO(3,3)=0.D0
      RETURN
      END
      SUBROUTINE BANGLE(XYZ,I,J,K,ANGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,*)
*********************************************************************
*
* BANGLE CALCULATES THE ANGLE BETWEEN ATOMS I,J, AND K. THE
*        CARTESIAN COORDINATES ARE IN XYZ.
*IVAL,
*********************************************************************
      D2IJ = (XYZ(1,I)-XYZ(1,J))**2+
     1       (XYZ(2,I)-XYZ(2,J))**2+
     2       (XYZ(3,I)-XYZ(3,J))**2
      D2JK = (XYZ(1,J)-XYZ(1,K))**2+
     1       (XYZ(2,J)-XYZ(2,K))**2+
     2       (XYZ(3,J)-XYZ(3,K))**2
      D2IK = (XYZ(1,I)-XYZ(1,K))**2+
     1       (XYZ(2,I)-XYZ(2,K))**2+
     2       (XYZ(3,I)-XYZ(3,K))**2
      XY = SQRT(D2IJ*D2JK)
      TEMP = 0.5D0 * (D2IJ+D2JK-D2IK) / XY
      IF (TEMP .GT. 1.0D0) TEMP=1.0D0
      IF (TEMP .LT. -1.0D0) TEMP=-1.0D0
      ANGLE = ACOS( TEMP )
      RETURN
      END
      SUBROUTINE DIHED(XYZ,I,J,K,L,ANGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,*)
*********************************************************************
*
*      DIHED CALCULATES THE DIHEDRAL ANGLE BETWEEN ATOMS I, J, K,
*            AND L.  THE CARTESIAN COORDINATES OF THESE ATOMS
*            ARE IN ARRAY XYZ.
*
*     DIHED IS A MODIFIED VERSION OF A SUBROUTINE OF THE SAME NAME
*           WHICH WAS WRITTEN BY DR. W. THEIL IN 1973.
*
*********************************************************************
      XI1=XYZ(1,I)-XYZ(1,K)
      XJ1=XYZ(1,J)-XYZ(1,K)
      XL1=XYZ(1,L)-XYZ(1,K)
      YI1=XYZ(2,I)-XYZ(2,K)
      YJ1=XYZ(2,J)-XYZ(2,K)
      YL1=XYZ(2,L)-XYZ(2,K)
      ZI1=XYZ(3,I)-XYZ(3,K)
      ZJ1=XYZ(3,J)-XYZ(3,K)
      ZL1=XYZ(3,L)-XYZ(3,K)
C      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
      DIST= SQRT(XJ1**2+YJ1**2+ZJ1**2)
      COSA=ZJ1/DIST
      IF(COSA.GT.1.0D0) COSA=1.0D0
      IF(COSA.LT.-1.0D0) COSA=-1.0D0
      DDD=1.0D0-COSA**2
      IF(DDD.LE.0.0) GO TO 10
      YXDIST=DIST* SQRT(DDD)
      IF(YXDIST.GT.1.0D-6) GO TO 20
   10 CONTINUE
      XI2=XI1
      XL2=XL1
      YI2=YI1
      YL2=YL1
      COSTH=COSA
      SINTH=0.D0
      GO TO 30
   20 COSPH=YJ1/YXDIST
      SINPH=XJ1/YXDIST
      XI2=XI1*COSPH-YI1*SINPH
      XL2=XL1*COSPH-YL1*SINPH
      YI2=XI1*SINPH+YI1*COSPH
      YJ2=XJ1*SINPH+YJ1*COSPH
      YL2=XL1*SINPH+YL1*COSPH
C      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
      COSTH=COSA
      SINTH=YJ2/DIST
   30 CONTINUE
      YI3=YI2*COSTH-ZI1*SINTH
      YL3=YL2*COSTH-ZL1*SINTH
      CALL DANG(XL2,YL3,XI2,YI3,ANGLE)
      IF (ANGLE .LT. 0.) ANGLE=4.0D0* ASIN(1.0D00)+ANGLE
      IF (ANGLE .GE. 6.2831853D0 ) ANGLE=0.D0
      RETURN
      END
      SUBROUTINE DANG(A1,A2,B1,B2,RCOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
**********************************************************************
*
*    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),
*          AND (B1,B2).  THE RESULT IS PUT IN RCOS.
*
**********************************************************************
      ZERO=1.0D-6
      IF( ABS(A1).LT.ZERO.AND. ABS(A2).LT.ZERO) GO TO 10
      IF( ABS(B1).LT.ZERO.AND. ABS(B2).LT.ZERO) GO TO 10
      ANORM=1.0D0/ SQRT(A1**2+A2**2)
      BNORM=1.0D0/ SQRT(B1**2+B2**2)
      A1=A1*ANORM
      A2=A2*ANORM
      B1=B1*BNORM
      B2=B2*BNORM
      SINTH=(A1*B2)-(A2*B1)
      COSTH=A1*B1+A2*B2
      IF(COSTH.GT.1.0D0) COSTH=1.0D0
      IF(COSTH.LT.-1.0D0) COSTH=-1.0D0
      RCOS= ACOS(COSTH)
      IF( ABS(RCOS).LT.4.0D-4) GO TO 10
      IF(SINTH.GT.0.D0) RCOS=4.0D0* ASIN(1.0D00)-RCOS
      RCOS=-RCOS
      RETURN
   10 RCOS=0.0D0
      RETURN
      END
