**************************************************************************
*
*       PACKAGE FOR APPLYING ALTERNATIVE GEOMETRIC CONSTRAINTS
*       PRODUCT OF NSERC SUMMER RESEARCH BY LAURENT MACKAY
*       USE "ALTCON" KEYWORD TO ACCESS THIS
*
*************************************************************************


      SUBROUTINE AALTCON(XYZ,DEGREE)
      INCLUDE 'SIZES'
      DOUBLE PRECISION XYZ(3,NUMATM),DEGREE
      COMMON /ALTCON / TRLB, ROTB, ROTD, ATMS, ICONXN, APPLIED,
     1                  VALS, NVALS
      INTEGER TRLB, ROTB, ROTD, ATMS, ICONXN(6,NUMATM),NVALS
      LOGICAL APPLIED
      DOUBLE PRECISION VALS(3*NUMATM)
      INTEGER I, ATOMS(*),TRANSL(*),BANGL(*),DIHDR(*),L
      INTEGER :: GETLENGTH
      POINTER (A,ATOMS),(T,TRANSL),(B,BANGL),(D,DIHDR)
      L=GETLENGTH(ROTD)
      A=ATMS
      T=TRLB
      B=ROTB
      D=ROTD
      I=0



   10 IF((I+5).LE.GETLENGTH(TRLB))THEN
      CALL SETRADI(XYZ,TRANSL(I+1),TRANSL(I+2),VALS(TRANSL(I+3)),
     1 ATOMS(TRANSL(I+4):TRANSL(I+5)),TRANSL(I+5)-TRANSL(I+4)+1)
      I=I+5
      GOTO 10
      ELSE
      GOTO 20
      ENDIF
   20 I=0
   30 IF((I+6).LE.GETLENGTH(ROTB))THEN
      CALL SETBANG(XYZ,BANGL(I+1),BANGL(I+2),BANGL(I+3),
     1 VALS(BANGL(I+4))/DEGREE,ATOMS(BANGL(I+5):BANGL(I+6)),
     2 BANGL(I+6)-BANGL(I+5)+1)
      I=I+6
      GOTO 30
      ELSE
      GOTO 40
      END IF
   40 I=0
   50 IF((I+7).LE.GETLENGTH(ROTD))THEN

      CALL SETDIHD(XYZ,DIHDR(I+1),DIHDR(I+2),DIHDR(I+3),DIHDR(I+4),
     1 VALS(DIHDR(I+5))/DEGREE,ATOMS(DIHDR(I+6):DIHDR(I+7)),
     2 DIHDR(I+7)-DIHDR(I+6)+1)
      I=I+7
      GOTO 50
      END IF


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
            LN=ATAN2(DOT_PRODUCT(JII,YP)/LN,DOT_PRODUCT(JII,XP)/LN)
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

      SUBROUTINE SETLINPEN(XYZ)
      INCLUDE 'SIZES'
      COMMON /LINPEN / PAXIS,POFF,PI,PJ,PK,KF,PENRGY,PDIST,POW,KI
      COMMON /AXES / XHAT(3),YHAT(3),ZHAT(3),OFF(3),ATOT
      INTEGER I,J,K,PI,PJ,PK
      DOUBLE PRECISION PAXIS(3),POFF(3),XYZ(3,NUMATM),LN,KF,ATOT(3,3),
     1                   TMP(3),XHAT,YHAT,ZHAT,OFF,PENRGY,PDIST,POW,KI
      LOGICAL DET

      IF(PI.EQ.0)RETURN

      POFF(1:3)=XYZ(1:3,PI)
      PAXIS(1)=XYZ(1,PJ)-POFF(1)
      PAXIS(2)=XYZ(2,PJ)-POFF(2)
      PAXIS(3)=XYZ(3,PJ)-POFF(3)

      LN=SQRT(PAXIS(1)**2+PAXIS(2)**2+PAXIS(3)**2)
      PAXIS(1)=PAXIS(1)/LN
      PAXIS(2)=PAXIS(2)/LN
      PAXIS(3)=PAXIS(3)/LN



      END

      SUBROUTINE PENENERGY(XYZ,ESCF)
      INCLUDE 'SIZES'
      COMMON /LINPEN / PAXIS,POFF,PI,PJ,PK,KF,PENRGY,PDIST,POW,KI
      INTEGER I,J,PI,PJ,PK
      DOUBLE PRECISION PAXIS(3),POFF(3),XYZ(3,NUMATM),LN,SEP(3),
     1 ESCF,KF,PENRGY,PDIST,POW,TMP,KI
      LOGICAL DET

      PENRGY=0
      IF(PI.EQ.0)RETURN

      J=1
   10 SELECT CASE (J)
      CASE (1)
        I=PI
      CASE (2)
        I=PJ
      CASE (3)
        I=PK
      END SELECT
      IF(J.LT.4)THEN
          J=J+1
          SEP(1)=XYZ(1,I)-POFF(1)
          SEP(2)=XYZ(2,I)-POFF(2)
          SEP(3)=XYZ(3,I)-POFF(3)
          LN=DOT_PRODUCT(PAXIS,SEP)
          SEP(1)=POFF(1)+LN*PAXIS(1)
          SEP(2)=POFF(2)+LN*PAXIS(2)
          SEP(3)=POFF(3)+LN*PAXIS(3)
          SEP(1)=XYZ(1,I)-SEP(1)
          SEP(2)=XYZ(2,I)-SEP(2)
          SEP(3)=XYZ(3,I)-SEP(3)
          LN=((SEP(1)**2+SEP(2)**2+SEP(3)**2)**((POW+1.D0)/2.D0))
          PENRGY=PENRGY+(LN*KF/(POW+1.D0))
          GOTO 10
      ENDIF
        ESCF=ESCF+PENRGY

      END

      SUBROUTINE PENFORCE(XYZ,DXYZ)
      INCLUDE 'SIZES'
      COMMON /LINPEN / PAXIS,POFF,PI,PJ,PK,KF,PENRGY,PDIST,POW,KI
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /GRADNT/ GRAD(MAXPAR),GNORM
      COMMON /STFU/ DAXIS
      INTEGER I,J,PI,PJ,PK,NA,NB,NC
      DOUBLE PRECISION PAXIS(3),POFF(3),XYZ(3,NUMATM),LN,SEP(3,3),
     1                   DXYZ(3,9*NUMATM),KF,PENRGY,PDIST,POW,GRAD,
     2                   GNORM,DAXIS,DTOT,KI
!      IF(PDIST.NE.0)THEN
!        POW=-REAL(INT(DLOG10(PDIST)))+1.D0
!        POW=MAX(1.D0,POW)
!      ELSE
!        POW=1.D0
!      ENDIF

      PDIST=0
      DTOT=0
      DAXIS=0
      IF(PI.EQ.0)RETURN

      J=1
   10 SELECT CASE (J)
      CASE (1)
        I=PI
      CASE (2)
        I=PJ
      CASE (3)
        I=PK
      END SELECT
      IF(J.LT.4)THEN

          SEP(1,J)=XYZ(1,I)-POFF(1)
          SEP(2,J)=XYZ(2,I)-POFF(2)
          SEP(3,J)=XYZ(3,I)-POFF(3)
          LN=DOT_PRODUCT(PAXIS,SEP(:,J))
          SEP(1,J)=POFF(1)+LN*PAXIS(1)
          SEP(2,J)=POFF(2)+LN*PAXIS(2)
          SEP(3,J)=POFF(3)+LN*PAXIS(3)
          SEP(1,J)=XYZ(1,I)-SEP(1,J)
          SEP(2,J)=XYZ(2,I)-SEP(2,J)
          SEP(3,J)=XYZ(3,I)-SEP(3,J)

          LN=SEP(1,J)**2+SEP(2,J)**2+SEP(3,J)**2

          PDIST=PDIST+LN
          DAXIS=DAXIS+ABS(DOT_PRODUCT(DXYZ(:,I),PAXIS))

!          DXYZ(1,I)=0
!          DXYZ(2,I)=0
!          DXYZ(3,I)=0
          J=J+1
          GOTO 10
      ENDIF
      PDIST=SQRT(PDIST/3.D0)
      DAXIS=SQRT(DAXIS/3.D0)
!      IF(.NOT.DET.AND.GNORM.GT.0)THEN
!        KF=1.D1**(FLOOR(LOG10(GNORM))+1)
!        DET=.TRUE.
!      ENDIF
!      IF(PDIST.LT.1.D-1)THEN
!          IF(DAXIS.LT.10.D0)THEN
!              KF=1.D3
!          ELSE
!              KF=1.D1
!          ENDIF
!      ELSE
!          KF=1.D3
!      ENDIF

      J=1
   20 SELECT CASE (J)
      CASE (1)
        I=PI
      CASE (2)
        I=PJ
      CASE (3)
        I=PK
      END SELECT
      IF(J.LT.4)THEN
          DXYZ(1,I)=DXYZ(1,I)+KF*SEP(1,J)*(ABS(SEP(1,J))**(POW-1.D0))
          DXYZ(2,I)=DXYZ(2,I)+KF*SEP(2,J)*(ABS(SEP(2,J))**(POW-1.D0))
          DXYZ(3,I)=DXYZ(3,I)+KF*SEP(3,J)*(ABS(SEP(3,J))**(POW-1.D0))
          J=J+1
          GOTO 20
      ENDIF
!      CALL SETLINPEN(XYZ)
      END

