**************************************************************************
*
*       PACKAGE FOR APPLYING ALTERNATIVE GEOMETRIC CONSTRAINTS
*       PRODUCT OF NSERC SUMMER RESEARCH BY LAURENT MACKAY
*       USE "ALTCON" KEYWORD TO ACCESS THIS
*
*************************************************************************


      SUBROUTINE AALTCON(XYZ)
      INCLUDE 'SIZES'
      DOUBLE PRECISION XYZ(3,NUMATM),DEGREE
      COMMON /LIN / PI,PJ,PK
      COMMON /PERMUTE /PR,PRT
      INTEGER PI,PJ,PK,PR(NUMATM),PRT(NUMATM)
      COMMON /ALTCON / TRANSL, BANGL, DIHDR, ATOMS, ICONXN, APPLIED,
     1                  VALS, NVALS,NTRANS,NANGL,NDIHD,NATOM
      INTEGER TRANSL(5*NUMATM),NTRANS,BANGL(6*NUMATM),NANGL,
     1         DIHDR(7*NUMATM),NDIHD,ATOMS(NUMATM*NUMATM),
     2         NATOM,ICONXN(6,NUMATM),NVALS,TATM(1)
      LOGICAL APPLIED
      DOUBLE PRECISION VALS(3*NUMATM),ANGL
      INTEGER I,L
      INTEGER :: GETLENGTH

      DEGREE=90.D0/ASIN(1.D0)
      I=0



   10 IF((I+5).LE.NTRANS)THEN
      CALL SETRADI(XYZ,TRANSL(I+1),TRANSL(I+2),VALS(TRANSL(I+3)),
     1 ATOMS(TRANSL(I+4):TRANSL(I+5)),TRANSL(I+5)-TRANSL(I+4)+1)
      I=I+5
      GOTO 10
      ELSE
      GOTO 20
      ENDIF
   20 I=0
   30 IF((I+6).LE.NANGL)THEN
      CALL SETBANG(XYZ,BANGL(I+1),BANGL(I+2),BANGL(I+3),
     1 VALS(BANGL(I+4))/DEGREE,ATOMS(BANGL(I+5):BANGL(I+6)),
     2 BANGL(I+6)-BANGL(I+5)+1)
      I=I+6
      GOTO 30
      ELSE
      GOTO 40
      END IF
   40 I=0
   50 IF((I+7).LE.NDIHD)THEN

      CALL SETDIHD(XYZ,DIHDR(I+1),DIHDR(I+2),DIHDR(I+3),DIHDR(I+4),
     1 VALS(DIHDR(I+5))/DEGREE,ATOMS(DIHDR(I+6):DIHDR(I+7)),
     2 DIHDR(I+7)-DIHDR(I+6)+1)
      I=I+7
      GOTO 50
      END IF

      IF(PI.NE.O)THEN
          CALL BANGLE(XYZ,PRT(PI),PRT(PJ),PRT(PK),ANGL)
          IF(ANGL.NE.180/DEGREE)THEN
              IF(ANGL.LT.170/DEGREE)THEN
                  TATM(1)=PK
           CALL SETBANG(XYZ,PRT(PI),PRT(PJ),PRT(PK),180/DEGREE,TATM,1)
              ELSE
                  TATM(1)=PJ
           CALL SETBANG(XYZ,PRT(PI),PRT(PJ),PRT(PK),180/DEGREE,TATM,1)
              ENDIF
          ENDIF
      ENDIF


      END

      SUBROUTINE SETRADI(XYZ,I,J,R,ATMS,NATM)
      INCLUDE 'SIZES'
      COMMON /PERMUTE /PR,PRT
      INTEGER PR(NUMATM),PRT(NUMATM)
      DOUBLE PRECISION DEL(3),SC,R,XYZ(3,*)
      INTEGER NATM,ATMS(*)
      INTEGER I,J,II,JJ
**********************************************
*       SETS THE DISTANCE BETWEEN TWO ATOMS TO THE SPECIFIED LENGTH
*       THE SECOND ATOM WILL MOVE ALONG THE SEPARATION VECTOR
**********************************************
      SC=SQRT((XYZ(1,I)-XYZ(1,J))**2+(XYZ(2,I)-XYZ(2,J))**2+
     1 (XYZ(3,I)-XYZ(3,J))**2)
      SC=(R/SC-1)
      DEL(1)=SC*(XYZ(1,J)-XYZ(1,I))
      DEL(2)=SC*(XYZ(2,J)-XYZ(2,I))
      DEL(3)=SC*(XYZ(3,J)-XYZ(3,I))
      DO 10 II=1,NATM
        JJ=PRT(ATMS(II))
        XYZ(1,JJ)=XYZ(1,JJ)+DEL(1)
        XYZ(2,JJ)=XYZ(2,JJ)+DEL(2)
        XYZ(3,JJ)=XYZ(3,JJ)+DEL(3)
   10  CONTINUE
      RETURN
      END
      SUBROUTINE SETBANG(XYZ,I,J,K,ANGL,ATMS,NATM)
      INCLUDE 'SIZES'
      COMMON /PERMUTE /PR,PRT
      INTEGER PR(NUMATM),PRT(NUMATM)
      DOUBLE PRECISION XP(3),YP(3),ZP(3),OP(3),JII(3),R,LN,C,S,DEL,
     1 ANGL,XYZ(3,NUMATM),FOO1, FOO2
      INTEGER NATM,ATMS(*)
      INTEGER I,J,K,II,JJ
**********************************************
*       SETS THE BOND ANGLE BETWEEN THREE ATOMS
*       THIS MAY FAIL IF THE ATOMS ARE QUASI-COLINEAR
**********************************************
      YP(1)=XYZ(1,K)-XYZ(1,J)
      YP(2)=XYZ(2,K)-XYZ(2,J)
      YP(3)=XYZ(3,K)-XYZ(3,J)
      LN=SQRT(YP(1)**2+YP(2)**2+YP(3)**2)
      YP(1)=YP(1)/LN
      YP(2)=YP(2)/LN
      YP(3)=YP(3)/LN
      XP(1)=XYZ(1,I)-XYZ(1,J)
      XP(2)=XYZ(2,I)-XYZ(2,J)
      XP(3)=XYZ(3,I)-XYZ(3,J)
      LN=SQRT(XP(1)**2+XP(2)**2+XP(3)**2)
      XP(1)=XP(1)/LN
      XP(2)=XP(2)/LN
      XP(3)=XP(3)/LN
      LN=DOT_PRODUCT(XP,YP)
      YP(1)=YP(1)-LN*XP(1)
      YP(2)=YP(2)-LN*XP(2)
      YP(3)=YP(3)-LN*XP(3)
      LN=SQRT(YP(1)**2+YP(2)**2+YP(3)**2)
      IF(LN.LT.1.D-4)THEN
      YP(1)=-XP(2)
      YP(2)=XP(1)
      YP(3)=0
      LN=SQRT(YP(1)**2+YP(2)**2+YP(3)**2)
      ENDIF
      YP(1)=YP(1)/LN
      YP(2)=YP(2)/LN
      YP(3)=YP(3)/LN
      CALL BANGLE(XYZ,K,J,I,DEL)
      DEL=ANGL-DEL
      ZP(1)=XP(2)*YP(3)-XP(3)*YP(2)
      ZP(2)=XP(3)*YP(1)-XP(1)*YP(3)
      ZP(3)=XP(1)*YP(2)-XP(2)*YP(1)
      LN=SQRT(DOT_PRODUCT(ZP,ZP))
      ZP(1)=ZP(1)/LN
      ZP(2)=ZP(2)/LN
      ZP(3)=ZP(3)/LN
      DO 10 II=1,NATM
        JJ=PRT(ATMS(II))
        IF(ATMS(II).NE.I)THEN
            JII(1)=XYZ(1,JJ)-XYZ(1,J)
            JII(2)=XYZ(2,JJ)-XYZ(2,J)
            JII(3)=XYZ(3,JJ)-XYZ(3,J)
            LN=SQRT(JII(1)**2+JII(2)**2+JII(3)**2)
            IF(LN.GT.0)
     1       LN=ATAN2(DOT_PRODUCT(JII,YP)/LN,DOT_PRODUCT(JII,XP)/LN)
            LN=LN+DEL
            C=COS(LN)
            S=SIN(LN)
            LN=DOT_PRODUCT(JII,ZP)
            OP(1)=XYZ(1,J)+ZP(1)*LN
            OP(2)=XYZ(2,J)+ZP(2)*LN
            OP(3)=XYZ(3,J)+ZP(3)*LN
            R=SQRT((XYZ(1,JJ)-OP(1))**2+
     2   (XYZ(2,JJ)-OP(2))**2+
     3   (XYZ(3,JJ)-OP(3))**2)
            XYZ(1,JJ)=OP(1)+R*(C*XP(1)+S*YP(1))
            XYZ(2,JJ)=OP(2)+R*(C*XP(2)+S*YP(2))
            XYZ(3,JJ)=OP(3)+R*(C*XP(3)+S*YP(3))
        ENDIF
   10  CONTINUE
      RETURN
      END
      SUBROUTINE SETDIHD(XYZ,I,J,K,L,ANGL,ATMS,NATM)
      INCLUDE 'SIZES'
      COMMON /PERMUTE /PR,PRT
      INTEGER PR(NUMATM),PRT(NUMATM)
      DOUBLE PRECISION OP(3),XP(3),YP(3),ZP(3),JII(3),R,LN,C,S,DEL
     1 ,ANGL,XYZ(3,NUMATM)
      INTEGER NATM, ATMS(*),I,J,K,L,II,JJ
**********************************************************************
*       SETS THE DIHERDAL ANGLE OF THE IJK PLANE WRT JKL TO THAT SPECIFIED
*       WILL FAIL IF THE PLANES ARE DEFINED BY A QUASI-COLINEAR SET
**********************************************************************
      C=COS(ANGL)
      S=SIN(ANGL)
      ZP(1)=XYZ(1,K)-XYZ(1,J)
      ZP(2)=XYZ(2,K)-XYZ(2,J)
      ZP(3)=XYZ(3,K)-XYZ(3,J)
      LN=SQRT(ZP(1)**2+ZP(2)**2+ZP(3)**2)
      ZP(1)=ZP(1)/LN
      ZP(2)=ZP(2)/LN
      ZP(3)=ZP(3)/LN
      XP(1)=XYZ(1,I)-XYZ(1,J)
      XP(2)=XYZ(2,I)-XYZ(2,J)
      XP(3)=XYZ(3,I)-XYZ(3,J)
      LN=DOT_PRODUCT(XP,ZP)
      XP(1)=XP(1)-LN*ZP(1)
      XP(2)=XP(2)-LN*ZP(2)
      XP(3)=XP(3)-LN*ZP(3)
      LN=SQRT(XP(1)**2+XP(2)**2+XP(3)**2)
      XP(1)=XP(1)/LN
      XP(2)=XP(2)/LN
      XP(3)=XP(3)/LN
      YP(1)=ZP(2)*XP(3)-ZP(3)*XP(2)
      YP(2)=ZP(3)*XP(1)-ZP(1)*XP(3)
      YP(3)=ZP(1)*XP(2)-ZP(2)*XP(1)
      CALL DIHED(XYZ,I,J,K,L,DEL)
      DEL=ANGL-DEL
      DO 10 II=1,NATM
        JJ=PRT(ATMS(II))
        JII(1)=XYZ(1,JJ)-XYZ(1,J)
        JII(2)=XYZ(2,JJ)-XYZ(2,J)
        JII(3)=XYZ(3,JJ)-XYZ(3,J)
        LN=DOT_PRODUCT(JII,ZP)
        OP(1)=XYZ(1,J)+ZP(1)*LN
        OP(2)=XYZ(2,J)+ZP(2)*LN
        OP(3)=XYZ(3,J)+ZP(3)*LN
        R=SQRT((OP(1)-XYZ(1,JJ))**2+(OP(2)-XYZ(2,JJ))**2
     1   +(OP(3)-XYZ(3,JJ))**2)
        CALL DIHED(XYZ,I,J,K,JJ,LN)
        LN=LN+DEL
        C=COS(LN)
        S=SIN(LN)
        XYZ(1,JJ)=OP(1)+R*(C*XP(1)+S*YP(1))
        XYZ(2,JJ)=OP(2)+R*(C*XP(2)+S*YP(2))
        XYZ(3,JJ)=OP(3)+R*(C*XP(3)+S*YP(3))
   10  CONTINUE
      RETURN
      END

