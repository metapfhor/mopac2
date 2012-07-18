**************************************************************************
*
*       PACKAGE FOR READING IN ALTERNATIVE GEOMETRIC CONSTRAINTS
*       PRODUCT OF NSERC SUMMER RESEARCH BY LAURENT MACKAY
*       USE "ALTCON" KEYWORD TO ACCESS THIS
*
*************************************************************************

      SUBROUTINE RALTCON(IREAD,LOPT)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      COMMON /ALTCON / TRLB, ROTB, ROTD, ATMS, ICONXN, APPLIED,
     1                  VALS, NVALS
      INTEGER TRLB, ROTB, ROTD, ATMS, ICONXN(6,NUMATM),NVALS
      LOGICAL APPLIED
      DOUBLE PRECISION VALS(3*NUMATM)
      LOGICAL LEADSP
      CHARACTER (LEN=MAXCHAR) :: LINE
      CHARACTER (LEN=MAXCHAR) :: CHUNK
      CHARACTER (LEN=MAXCHAR) :: RN
      DOUBLE PRECISION TMPC(8),XYZ(3,NUMATM)
      INTEGER II,JJ,EMPTY(1)
      INTEGER TATMS
      DIMENSION ISTART(MAXCHAR/2), LOPT(3,NUMATM)
      CHARACTER SPACE*1,COMMA*1,LSQB*1,RSQB*1,LCRB*1,RCRB*1
      DATA COMMA,SPACE,LSQB,RSQB,LCRB,RCRB,II,JJ,NVALS
     1 /',',' ','[',']','{','}',0,0,0/

      TRLB=NEWARR(5)
      ROTB=NEWARR(6)
      ROTD=NEWARR(7)
      ATMS=NEWARR(0)

      CALL BLDPERM(IREAD)

   10 READ(IREAD,'(A)',END=140)LINE
      IF(LINE.EQ.' ') GOTO 130

      LEADSP=.TRUE.
      NVALUE=0
      DO 20 I=1,MAXCHAR
         IF (LEADSP.AND.LINE(I:I).NE.SPACE) THEN
            NVALUE=NVALUE+1
            ISTART(NVALUE)=I
         END IF
         LEADSP=(LINE(I:I).EQ.SPACE)
   20 CONTINUE

      II=0
      DO WHILE(II.LE.MAXCHAR)
        II=II+1
        IF(LINE(II:II).EQ.LSQB)THEN
            JJ=II
            DO WHILE(JJ.LE.MAXCHAR)
            JJ=JJ+1
                IF(LINE(JJ:JJ).EQ.COMMA)THEN
                    LINE(JJ:JJ)=';'
                ELSEIF(LINE(JJ:JJ).EQ.RSQB)THEN
                    II=JJ+1
                    JJ=MAXCHAR
                ENDIF
           END DO
        ENDIF
      END DO

      TMPC(7)=0
      DO I=1,6
       TMPC(I)=0
       TMPC(I)=READN(LINE,ISTART(I))
       IF(TMPC(7).EQ.0.AND.TMPC(I).EQ.0)TMPC(7)=I
      END DO


      IF(TMPC(7).EQ.3) THEN
       GOTO 40
      ELSEIF(TMPC(7).EQ.4) THEN
       GOTO 50
      ELSEIF(TMPC(7).EQ.5) THEN
       GOTO 60
      ENDIF

C     FIRST INFO ON THIS LINE IS A TRANSLATION
   40 LINE=LINE(ISTART(2)-1:)
      DO I=1,MAXCHAR/2
           IF(LINE.EQ.' ')EXIT
          CALL SPLIT(LINE,',',CHUNK)
          LEADSP=.TRUE.
          NVALUE=0
          DO II=1,MAXCHAR/2
            ISTART(II)=0
          END DO
          DO II=1,MAXCHAR
             IF (LEADSP.AND.CHUNK(II:II).NE.SPACE) THEN
                NVALUE=NVALUE+1
                ISTART(NVALUE)=II
             END IF
             LEADSP=(CHUNK(II:II).EQ.SPACE)
          END DO

          TMPC(7)=0
          DO  II=1,5
           TMPC(II+1)=0
           TMPC(II+1)=READN(CHUNK,ISTART(II))
           IF(TMPC(7).EQ.0.AND.TMPC(II+1).EQ.0)TMPC(7)=II
          END DO

          IF(TMPC(7).EQ.2) THEN
            IF(ISTART(4).EQ.0)THEN
                RN(:)=''
            ELSE
                RN(1:)=CHUNK(ISTART(4):)
           ENDIF
            CALL ADDTRLB(TMPC,RN,INDEX(CHUNK,'F').NE.0,LOPT)
          ELSEIF(TMPC(7).EQ.3) THEN
            IF(ISTART(5).EQ.0)THEN
                RN(:)=''
            ELSE
                RN(1:)=CHUNK(ISTART(5):)
          ENDIF
            CALL ADDANGLE(TMPC,RN,INDEX(CHUNK,'F').NE.0,LOPT)

          ELSEIF(TMPC(7).EQ.4) THEN
             IF(ISTART(6).EQ.0)THEN
                RN(:)=''
            ELSE
                RN(1:)=CHUNK(ISTART(6):)
          ENDIF
            CALL ADDDIHDR(TMPC,RN,INDEX(CHUNK,'F').NE.0,LOPT)
          ENDIF



      END DO


      GOTO 10

C     FIRST INFO ON THIS LINE IS A BOND ANGLE
   50 LINE=LINE(ISTART(3)-1:)
      DO I=1,MAXCHAR/2
           IF(LINE.EQ.' ')EXIT
          CALL SPLIT(LINE,',',CHUNK)
          LEADSP=.TRUE.
          NVALUE=0
          DO II=1,MAXCHAR/2
            ISTART(II)=0
          END DO
          DO II=1,MAXCHAR
             IF (LEADSP.AND.CHUNK(II:II).NE.SPACE) THEN
                NVALUE=NVALUE+1
                ISTART(NVALUE)=II
             END IF
             LEADSP=(CHUNK(II:II).EQ.SPACE)
          END DO

          TMPC(7)=0
          DO  II=1,4
           TMPC(II+2)=0
           TMPC(II+2)=READN(CHUNK,ISTART(II))
           IF(TMPC(7).EQ.0.AND.TMPC(II+2).EQ.0)TMPC(7)=II
          END DO

          IF(TMPC(7).EQ.2) THEN
            IF(ISTART(4).EQ.0)THEN
                RN(:)=''
            ELSE
                RN(1:)=CHUNK(ISTART(4):)
          ENDIF
            CALL ADDANGLE(TMPC,RN,
     1       INDEX(CHUNK,'F').NE.0,LOPT)
          ELSEIF(TMPC(7).EQ.3) THEN
            IF(ISTART(5).EQ.0)THEN
                RN(:)=''
            ELSE
                RN(1:)=CHUNK(ISTART(5):)
           ENDIF
            CALL ADDDIHDR(TMPC,RN,INDEX(CHUNK,'F').NE.0,LOPT)


          ENDIF



      END DO


      GOTO 10

C     FIRST INFO ON THIS LINE IS A DIHEDRAL ANGLE
   60 LINE=LINE(ISTART(4)-1:)
      DO I=1,MAXCHAR/2
           IF(LINE.EQ.' ')EXIT
          CALL SPLIT(LINE,',',CHUNK)
          LEADSP=.TRUE.
          NVALUE=0
          DO II=1,MAXCHAR/2
            ISTART(II)=0
          END DO
          DO II=1,MAXCHAR
             IF (LEADSP.AND.CHUNK(II:II).NE.SPACE) THEN
                NVALUE=NVALUE+1
                ISTART(NVALUE)=II
             END IF
             LEADSP=(CHUNK(II:II).EQ.SPACE)
          END DO

          TMPC(7)=0
          DO  II=1,3
           TMPC(II+3)=0
           TMPC(II+3)=READN(CHUNK,ISTART(II))
           IF(TMPC(7).EQ.0.AND.TMPC(II+3).EQ.0)TMPC(7)=II
          ENDDO

          IF(TMPC(7).EQ.2) THEN
          IF(ISTART(4).EQ.0)THEN
            RN(:)=''
            ELSE
            RN(1:)=CHUNK(ISTART(4):)
          ENDIF
            CALL ADDDIHDR(TMPC,RN,INDEX(CHUNK,'F').NE.0,LOPT)
          ENDIF
      END DO


      GOTO 10

C       End Laurent
  130 RETURN
  140 BACKSPACE(IREAD)
      RETURN
      END

      SUBROUTINE SPLIT(STR,DELIMS,BEFORE)

        ! ROUTINE FINDS THE FIRST INSTANCE OF A CHARACTER FROM 'DELIMS' IN THE
        ! THE STRING 'STR'. THE CHARACTERS BEFORE THE FOUND DELIMITER ARE
        ! OUTPUT IN 'BEFORE'. THE CHARACTERS AFTER THE FOUND DELIMITER ARE
        ! OUTPUT IN 'STR'. THE OPTIONAL OUTPUT CHARACTER 'SEP' CONTAINS THE
        ! FOUND DELIMITER. A DELIMITER IN 'STR' IS TREATED LIKE AN ORDINARY
        ! CHARACTER IF IT IS PRECEDED BY A BACKSLASH (\). IF THE BACKSLASH
        ! CHARACTER IS DESIRED IN 'STR', THEN PRECEDE IT WITH ANOTHER BACKSLASH.

        !Modified version of algorithm written by Benthien, George
        !http://gbenthien.net/strings/index.html

      CHARACTER(LEN=*) :: STR,DELIMS,BEFORE
      LOGICAL :: PRES
      CHARACTER :: CH,CHA

      STR=ADJUSTL(STR)
      CALL COMPACT(STR)
      LENSTR=LEN_TRIM(STR)
      IF(LENSTR.EQ.0) RETURN        ! STRING STR IS EMPTY
      K=0
      IBSL=0                        ! BACKSLASH INITIALLY INACTIVE
      BEFORE=' '
      DO I=1,LENSTR
         CH=STR(I:I)
         IF(IBSL.EQ.1) THEN          ! BACKSLASH ACTIVE
            K=K+1
            BEFORE(K:K)=CH
            IBSL=0
            CYCLE
         END IF
         IF(CH.EQ.'\') THEN          ! BACKSLASH WITH BACKSLASH INACTIVE
            K=K+1
            BEFORE(K:K)=CH
            IBSL=1
            CYCLE
         END IF
         IPOS=INDEX(DELIMS,CH)
         IF(IPOS.EQ. 0) THEN          ! CHARACTER IS NOT A DELIMITER
            K=K+1
            BEFORE(K:K)=CH
            CYCLE
         END IF
         IF(CH.NE.' ') THEN          ! CHARACTER IS A DELIMITER THAT IS NOT A SPACE
            STR=STR(I+1:)
            EXIT
         END IF
         CHA=STR(I+1:I+1)            ! CHARACTER IS A SPACE DELIMITER
         IPOSA=INDEX(DELIMS,CHA)
         IF(IPOSA.GT.0) THEN          ! NEXT CHARACTER IS A DELIMITER
            STR=STR(I+2:)
            EXIT
         ELSE
            STR=STR(I+1:)
            EXIT
         END IF
      END DO
      IF(I.GE.LENSTR) STR=''
      STR=ADJUSTL(STR)              ! REMOVE INITIAL SPACES


      END

      SUBROUTINE COMPACT(STR)

      ! CONVERTS MULTIPLE SPACES AND TABS TO SINGLE SPACES; DELETES CONTROL CHARACTERS;
      ! REMOVES INITIAL SPACES.

      CHARACTER(LEN=*):: STR
      CHARACTER(LEN=1):: CH
      CHARACTER(LEN=LEN_TRIM(STR)):: OUTSTR

      STR=ADJUSTL(STR)
      LENSTR=LEN_TRIM(STR)
      OUTSTR=' '
      ISP=0
      K=0

      DO I=1,LENSTR
        CH=STR(I:I)
        ICH=IACHAR(CH)

        SELECT CASE(ICH)

          CASE(9,32)     ! SPACE OR TAB CHARACTER
            IF(ISP.EQ.0) THEN
              K=K+1
              OUTSTR(K:K)=' '
            END IF
            ISP=1

          CASE(33:)      ! NOT A SPACE, QUOTE, OR CONTROL CHARACTER
            K=K+1
            OUTSTR(K:K)=CH
            ISP=0

        END SELECT

      END DO

      STR=ADJUSTL(OUTSTR)

      END

      SUBROUTINE ADDRANGE(A,RN,I)
      INTEGER :: STRLEN
      CALL ADDRANGEL(A,RN,I,STRLEN(RN))
      END

      SUBROUTINE ADDRANGEL(A,RN,I,S)
      INTEGER I,LENSTR,A,RANGESIZE,J,K,L,M,S
      CHARACTER RN*(*),TMP*(S),CHUNK*(S)
      LOGICAL :: NINARR
      DOUBLE PRECISION :: READA
      TMP=RN(2:S-1)

      CALL ADDVAL(A,I)
      DO II=0,S
        IF(TMP.EQ.' ')EXIT
        CALL SPLIT(TMP,";",CHUNK)
        IF(INDEX(CHUNK,'-').NE.0)THEN
             J=INT(READA(CHUNK,1))
             K=INT(READA(CHUNK,INDEX(CHUNK,'-')+1))
             L=MAX(J,K)
             K=MIN(J,K)
             DO J=K,L
                IF(NINARR(A,J))THEN
                  CALL ADDVAL(A,J)
                ENDIF
             END DO
        ELSE
            J=INT(READA(CHUNK,1))
            IF(NINARR(A,J))THEN
               CALL ADDVAL(A,J)
            ENDIF
        ENDIF
      END DO
      END

      INTEGER FUNCTION STRLEN(STR)
      CHARACTER STR*(*)
      INTEGER I
      DO I=1,LEN(STR)
      IF(STR(I:I).EQ.'')THEN
          STRLEN=I-1
          RETURN
      ENDIF
      END DO
      END

      LOGICAL FUNCTION NINARR(A,I)
      INTEGER A,I,N
      INTEGER :: GETLENGTH,GETVAL
      NINARR=.TRUE.
      N=GETLENGTH(A)
      DO J=1,N
        IF(GETVAL(A,J).EQ.I)THEN
            NINARR=.FALSE.
            RETURN
        ENDIF
      END DO
      RETURN
      END


      SUBROUTINE ADDTRLB(CNX,RN,F,LOPT)
      INCLUDE 'SIZES'
      COMMON /ALTCON / TRLB, ROTB, ROTD, ATMS, ICONXN, APPLIED,
     1                  VALS, NVALS
      COMMON /PERMUTE / PR,PRT
      INTEGER PR(NUMATM),PRT(NUMATM)
      INTEGER TRLB, ROTB, ROTD, ATMS, ICONXN(6,NUMATM),NVALS
      LOGICAL APPLIED
      DOUBLE PRECISION VALS(3*NUMATM)
      INTEGER I,J, CONX(3),LOPT(3,NUMATM)
      INTEGER :: GETLENGTH
      CHARACTER (LEN=MAXCHAR) :: RN
      DOUBLE PRECISION CNX(4)
      LOGICAL F
      CONX(1)=PRT(INT(CNX(1)))
      CONX(2)=PRT(INT(CNX(2)))
      CALL ADDVAL(TRLB,CONX(1))
      CALL ADDVAL(TRLB,CONX(2))
      NVALS=NVALS+1
      VALS(NVALS)=CNX(4)
      CALL ADDVAL(TRLB,NVALS)

      CALL ADDVAL(TRLB,GETLENGTH(ATMS)+1)
      CALL ADDRANGE(ATMS,RN,CONX(2))
      CALL ADDVAL(TRLB,GETLENGTH(ATMS))

      IF(F)THEN
      IF(CONX(2).LT.CONX(1))THEN
        CONX(3)=CONX(1)
        CONX(1)=CONX(2)
        CONX(2)=CONX(3)
      ENDIF

      IF(ICONXN(1,CONX(2)).EQ.0.OR.ICONXN(1,CONX(2)).EQ.CONX(1))THEN
            ICONXN(1,CONX(2))=CONX(1)
            ICONXN(4,CONX(2))=1
            LOPT(1,CONX(2))=0
      ELSE
      GOTO 330
      ENDIF
      RETURN
  330 WRITE(6,'(A,I3.1,A4,I3.1,A,I3.1,I3.1,I3.1,A)')
     1 'Impossible Gemoetric Constraints:',CONX(1), '->',
     2 CONX(2) ,' (',PRT(ICONXN(1,CONX(2))),PRT(ICONXN(2,CONX(2))),
     3 PRT(ICONXN(3,CONX(2))),' )  '//RN
      STOP
      RETURN
      ENDIF
      END

      SUBROUTINE ADDANGLE(CNX,RN,F,LOPT)
      INCLUDE 'SIZES'
      COMMON /ALTCON / TRLB, ROTB, ROTD, ATMS, ICONXN, APPLIED,
     1                  VALS, NVALS
      COMMON /PERMUTE / PR,PRT
      INTEGER PR(NUMATM),PRT(NUMATM)
      INTEGER TRLB, ROTB, ROTD, ATMS, ICONXN(6,NUMATM),NVALS
      LOGICAL APPLIED
      DOUBLE PRECISION VALS(3*NUMATM)
      INTEGER I,J, CONX(4), PERMUTE,LOPT(3,NUMATM)
      INTEGER :: GETLENGTH
      DOUBLE PRECISION CNX(5)
      LOGICAL F

      CONX(1)=PRT(INT(CNX(1)))
      CONX(2)=PRT(INT(CNX(2)))
      CONX(3)=PRT(INT(CNX(3)))



      CALL ADDVAL(ROTB,CONX(1))
      CALL ADDVAL(ROTB,CONX(2))
      CALL ADDVAL(ROTB,CONX(3))

      NVALS=NVALS+1
      CALL ADDVAL(ROTB,NVALS)
      VALS(NVALS)=CNX(5)
      CALL ADDVAL(ROTB,GETLENGTH(ATMS)+1)
      CALL ADDRANGE(ATMS,RN,CONX(3))
      CALL ADDVAL(ROTB,GETLENGTH(ATMS))

      IF(F)THEN
      IF((ICONXN(2,CONX(3)).EQ.0.OR.ICONXN(2,CONX(3)).EQ.CONX(1))
     1.AND.(ICONXN(1,CONX(3)).EQ.0.OR.ICONXN(1,CONX(3)).EQ.CONX(2)))THEN
        ICONXN(2,CONX(3))=CONX(1)
        ICONXN(1,CONX(3))=CONX(2)
        ICONXN(4,CONX(3))=1
        IF(ABS(CNX(5)-180.D0).LT.1D-3)THEN
            LOPT(3,CONX(3))=0
        ENDIF
        ICONXN(5,CONX(3))=1
        LOPT(2,CONX(3))=0


        RETURN
      ENDIF

      PERMUTE=0
   10 IF(PERMUTE.LE.2)THEN
        PERMUTE=PERMUTE+1
        CONX(4)=CONX(1)
        CONX(1)=CONX(2)
        CONX(2)=CONX(3)
        CONX(3)=CONX(4)
      ELSE
        GOTO 20
      ENDIF

      IF((ICONXN(2,CONX(3)).EQ.0.OR.ICONXN(2,CONX(3)).EQ.CONX(1))
     1.AND.(ICONXN(1,CONX(3)).EQ.0.OR.ICONXN(1,CONX(3)).EQ.CONX(2)))THEN
        ICONXN(2,CONX(3))=CONX(1)
        ICONXN(1,CONX(3))=CONX(2)
        ICONXN(5,CONX(3))=1

            LOPT(2,CONX(3))=0
            RETURN
      ELSE
        GOTO 10
      ENDIF
      RETURN
   20 WRITE(6,'(A)')'Impossible Gemoetric Constraints:', RN
      STOP
      RETURN
      ENDIF
      END

      SUBROUTINE ADDDIHDR(CNX,RN,F,LOPT)
      INCLUDE 'SIZES'
      COMMON /ALTCON / TRLB, ROTB, ROTD, ATMS, ICONXN, APPLIED,
     1                  VALS, NVALS
      COMMON /PERMUTE / PR,PRT
      INTEGER PR(NUMATM),PRT(NUMATM)
      INTEGER TRLB, ROTB, ROTD, ATMS, ICONXN(6,NUMATM),NVALS
      LOGICAL APPLIED
      DOUBLE PRECISION VALS(3*NUMATM)
      INTEGER I,J, CONX(5), PERMUTE,LOPT(3,NUMATM)
      INTEGER :: GETLENGTH
      DOUBLE PRECISION CNX(6)
      LOGICAL F
      CONX(1)=PRT(INT(CNX(1)))
      CONX(2)=PRT(INT(CNX(2)))
      CONX(3)=PRT(INT(CNX(3)))
      CONX(4)=PRT(INT(CNX(4)))
      CALL ADDVAL(ROTD,CONX(1))
      CALL ADDVAL(ROTD,CONX(2))
      CALL ADDVAL(ROTD,CONX(3))
      CALL ADDVAL(ROTD,CONX(4))
      NVALS=NVALS+1
      CALL ADDVAL(ROTD,NVALS)

      CALL ADDVAL(ROTD,GETLENGTH(ATMS)+1)
      VALS(NVALS)=CNX(6)

      CALL ADDRANGE(ATMS,RN,CONX(4))

      CALL ADDVAL(ROTD,GETLENGTH(ATMS))
      IF(F)THEN
      IF((ICONXN(3,CONX(4)).EQ.0.OR.ICONXN(3,CONX(4)).EQ.CONX(1))
     1.AND.(ICONXN(2,CONX(4)).EQ.0.OR.ICONXN(2,CONX(4)).EQ.CONX(2))
     2.AND.(ICONXN(1,CONX(4)).EQ.0.OR.ICONXN(1,CONX(4)).EQ.CONX(3)))THEN
        ICONXN(3,CONX(4))=CONX(1)
        ICONXN(2,CONX(4))=CONX(2)
        ICONXN(1,CONX(4))=CONX(3)
        ICONXN(6,CONX(4))=1
        LOPT(3,CONX(4))=0
        RETURN
      ENDIF

      PERMUTE=0
   10 IF(PERMUTE.LE.3)THEN
        PERMUTE=PERMUTE+1
        CONX(5)=CONX(1)
        CONX(1)=CONX(2)
        CONX(2)=CONX(3)
        CONX(3)=CONX(4)
        CONX(4)=CONX(5)
      ELSE
        GOTO 20
      ENDIF

      IF((ICONXN(3,CONX(4)).EQ.0.OR.ICONXN(3,CONX(4)).EQ.CONX(1))
     1.AND.(ICONXN(2,CONX(4)).EQ.0.OR.ICONXN(2,CONX(4)).EQ.CONX(2))
     2.AND.(ICONXN(1,CONX(4)).EQ.0.OR.ICONXN(1,CONX(4)).EQ.CONX(3)))THEN
        ICONXN(3,CONX(4))=CONX(1)
        ICONXN(2,CONX(4))=CONX(2)
        ICONXN(1,CONX(4))=CONX(3)
        ICONXN(6,CONX(4))=1
        LOPT(3,CONX(4))=0
        RETURN
        ELSE
          GOTO 10
      ENDIF
      RETURN
   20 WRITE(6,'(A)')'Impossible Gemoetric Constraints:', RN
      STOP
      RETURN
      ENDIF
      END

      SUBROUTINE INITPEN(I,J,K)
      INCLUDE 'SIZES'
      COMMON /LIN / PI,PJ,PK
      INTEGER I,J,K,PI,PJ,PK,II
      COMMON /KEYWRD/ KEYWRD
      CHARACTER KEYWRD*241
      DOUBLE PRECISION :: READA
      DATA PI,PJ,PK,KF,POW /0,0,0,1.D1,2.D0/
      II=INDEX(KEYWRD,"KF=")
      IF(INDEX(KEYWRD,"KF=").NE.0)THEN
        KF=READA(KEYWRD,INDEX(KEYWRD,"KF=")+3)
      ENDIF
      KF=0
      IF(INDEX(KEYWRD,"PEN-POW=").NE.0)THEN
        POW=READA(KEYWRD,INDEX(KEYWRD,"PEN-POW=")+8)
      ENDIF

      IF(PI.NE.0)THEN
        WRITE(6,'(A)')'ONLY ONE LINEARIZATION ALLOWED'
      ENDIF
      IF(I.EQ.K.OR.I.EQ.J.OR.J.EQ.K)THEN
        WRITE(6,'(A,3I5.2)')'IMPOSSIBLE LINEARIZATION',I,J,K
        STOP
      ELSE
C       THINGS BECOME SOMEWHAT HECTIC IF WE ARE DEALING WITH THE FIRST TWO ATOMS
          IF(I.EQ.1.OR.J.EQ.1.OR.K.EQ.1)THEN
              IF(I.EQ.1)THEN
                PI=I
                IF(J.EQ.2)THEN
                    PJ=J
                    PK=K
                ELSE
                     PJ=K
                     PK=J
                ENDIF
              ENDIF
              IF(J.EQ.1)THEN
                  PI=J
                  IF(I.EQ.2)THEN
                      PJ=I
                      PK=K
                  ELSE
                     PJ=K
                     PK=I
                  ENDIF
              ENDIF
              IF(K.EQ.1)THEN
                  PI=K
                  IF(I.EQ.2)THEN
                      PJ=I
                      PK=J
                  ELSE
                      IF(J.EQ.2)THEN
                          PJ=J
                          PK=I
                      ELSE
                          PJ=I
                          PK=J
                      ENDIF
                  ENDIF
              ENDIF
          ELSE
              PI=I
              PJ=J
              PK=K
          ENDIF
      ENDIF



      END

      SUBROUTINE BLDPERM(IREAD)
      INCLUDE 'SIZES'
      COMMON /PERMUTE / PR,PRT
      COMMON /DEPND / DEP,MINDEP,NEXTDEP
      COMMON /LIN / PI,PJ,PK
      INTEGER PI,PJ,PK
      INTEGER PR(NUMATM),PRT(NUMATM)
      INTEGER NBKSPC, INQ(NUMATM), DEP(NUMATM)
      LOGICAL LEADSP
      CHARACTER (LEN=MAXCHAR) :: LINE
      CHARACTER (LEN=MAXCHAR) :: CHUNK
      CHARACTER (LEN=MAXCHAR) :: RN
      LOGICAL :: NINARR
      DOUBLE PRECISION :: READN
      DOUBLE PRECISION TMPC(8)
      INTEGER II,JJ,CONX(8),MINDEP,TMPDEP,TMPIND,
     1         NEXTDEP,TWO
      INTEGER TATMS
      DIMENSION ISTART(MAXCHAR/2), LOPT(3,NUMATM)
      CHARACTER SPACE*1,COMMA*1,LSQB*1,RSQB*1,LCRB*1,RCRB*1
      LOGICAL MOVE,DET
      DATA COMMA,SPACE,LSQB,RSQB,LCRB,RCRB,II,JJ,NVALS
     1 /',',' ','[',']','{','}',0,0,0/



      NBKSPC=0
   10 READ(IREAD,'(A)',END=140)LINE
      NBKSPC=NBKSPC+1
      IF(LINE.EQ.' ') GOTO 130

      LEADSP=.TRUE.
      NVALUE=0
      DO 20 I=1,MAXCHAR
         IF (LEADSP.AND.LINE(I:I).NE.SPACE) THEN
            NVALUE=NVALUE+1
            ISTART(NVALUE)=I
         END IF
         LEADSP=(LINE(I:I).EQ.SPACE)
   20 CONTINUE

      II=0
      DO WHILE(II.LE.MAXCHAR)
        II=II+1
        IF(LINE(II:II).EQ.LSQB)THEN
            JJ=II
            DO WHILE(JJ.LE.MAXCHAR)
            JJ=JJ+1
                IF(LINE(JJ:JJ).EQ.COMMA)THEN
                    LINE(JJ:JJ)=';'
                ELSEIF(LINE(JJ:JJ).EQ.RSQB)THEN
                    II=JJ+1
                    JJ=MAXCHAR
                ENDIF
           END DO
        ENDIF
      END DO

      TMPC(7)=0
      DO I=1,6
       TMPC(I)=0
       TMPC(I)=READN(LINE,ISTART(I))
       IF(TMPC(7).EQ.0.AND.TMPC(I).EQ.0)TMPC(7)=I
      END DO


      IF(TMPC(7).EQ.3) THEN
       GOTO 40
      ELSEIF(TMPC(7).EQ.4) THEN
       GOTO 50
!      ELSEIF(TMPC(7).EQ.5) THEN
!       GOTO 60
      ENDIF

C     FIRST INFO ON THIS LINE IS A TRANSLATION
   40 LINE=LINE(ISTART(2)-1:)
      DO I=1,MAXCHAR/2
           IF(LINE.EQ.' ')EXIT
          CALL SPLIT(LINE,',',CHUNK)
          LEADSP=.TRUE.
          NVALUE=0
          DO II=1,MAXCHAR/2
            ISTART(II)=0
          END DO
          DO II=1,MAXCHAR
             IF (LEADSP.AND.CHUNK(II:II).NE.SPACE) THEN
                NVALUE=NVALUE+1
                ISTART(NVALUE)=II
             END IF
             LEADSP=(CHUNK(II:II).EQ.SPACE)
          END DO

          TMPC(7)=0
          DO  II=1,5
           TMPC(II+1)=0
           TMPC(II+1)=READN(CHUNK,ISTART(II))
           IF(TMPC(7).EQ.0.AND.TMPC(II+1).EQ.0)TMPC(7)=II
          END DO

          IF(TMPC(7).EQ.2) THEN
            CONX(1)=INT(TMPC(1))
            CONX(2)=INT(TMPC(2))


            IF(DEP(CONX(2)).EQ.0)THEN
                DEP(CONX(2))=NEWARR(0)
                CALL ADDVAL(DEP(CONX(2)),CONX(1))
            ELSE
                CALL ADDVAL(DEP(CONX(2)),CONX(1))
            ENDIF

          ELSEIF(TMPC(7).EQ.3) THEN

          ELSEIF(TMPC(7).EQ.4) THEN

          ENDIF



      END DO


      GOTO 10

   50 LINE=LINE(ISTART(3)-1:)
      DO I=1,MAXCHAR/2
           IF(LINE.EQ.' ')EXIT
          CALL SPLIT(LINE,',',CHUNK)
          LEADSP=.TRUE.
          NVALUE=0
          DO II=1,MAXCHAR/2
            ISTART(II)=0
          END DO
          DO II=1,MAXCHAR
             IF (LEADSP.AND.CHUNK(II:II).NE.SPACE) THEN
                NVALUE=NVALUE+1
                ISTART(NVALUE)=II
             END IF
             LEADSP=(CHUNK(II:II).EQ.SPACE)
          END DO

          TMPC(7)=0
          DO  II=1,4
           TMPC(II+2)=0
           TMPC(II+2)=READN(CHUNK,ISTART(II))
           IF(TMPC(7).EQ.0.AND.TMPC(II+2).EQ.0)TMPC(7)=II
          END DO

          IF(TMPC(7).EQ.2) THEN
            CONX(1)=INT(TMPC(1))
            CONX(2)=INT(TMPC(2))
            CONX(3)=INT(TMPC(3))
            IF(ABS(TMPC(5)-180.D0).LT.1.D-3)
     1       CALL INITPEN(CONX(1),CONX(2),CONX(3))

          ENDIF



      END DO
      GOTO 10

  140 NBKSPC=NBKSPC+1
  130 DO I=1,NBKSPC
      BACKSPACE(IREAD)
      END DO

      DO I=1,NUMATM
          PR(I)=I
          PRT(I)=I
      END DO

      DO I=1,NUMATM
         MINDEP=0
         NEXTDEP=0
         DO J=NUMATM,1,-1
         IF(DEP(J).NE.0.AND.(.NOT.NINARR(DEP(J),I)))THEN
            NEXTDEP=MINDEP
            MINDEP=J
         ENDIF
         END DO
         MINDEP=MAX(MINDEP,NEXTDEP)
         IF(MINDEP.NE.0.AND.I.GT.MINDEP)THEN
            TMPDEP=DEP(I)
            TMPIND=PR(I)
            DO J=I,MINDEP+1,-1
                PR(J)=PR(J-1)
                DEP(J)=DEP(J-1)
            END DO
            PR(MINDEP)=TMPIND
            DEP(MINDEPS)=TMPDEP
         ENDIF
      END DO
      CALL CALCPRT

      IF(PI.NE.0)THEN
        MOVE=.TRUE.
        IF(PI.NE.1)THEN

            DO I=PI-1,1,-1
                MOVE=DEP(I).EQ.O.OR.NINARR(DEP(I),PRT(PI))
                IF(.NOT.MOVE)THEN
                    EXIT
                ENDIF
            ENDDO
            IF(MOVE)THEN
            TMPDEP=DEP(PI)

            DO I=PI,2,-1
                PR(I)=PR(I-1)
                DEP(I)=DEP(I-1)
            ENDDO
            DEP(1)=TMPDEP
            PR(1)=PI
            CALL CALCPRT()
            ENDIF
            PI=PRT(PI)
         ENDIF
!         IF(MOVE)THEN
!            TWO=MIN(PRT(PJ),PRT(PK))
!            DO I=TWO-1,2,-1
!                MOVE=DEP(I).EQ.O.OR.NINARR(DEP(I),PRT(TWO))
!                IF(.NOT.MOVE)THEN
!                    EXIT
!                ENDIF
!            ENDDO
!            IF(MOVE)THEN
!                TMPDEP=DEP(TWO)
!                TMPIND=PR(TWO)
!                DO I=TWO,3,-1
!                    PR(I)=PR(I-1)
!                    DEP(I)=DEP(I-1)
!                ENDDO
!                DEP(2)=TMPDEP
!                PR(2)=TMPIND
!                CALL CALCPRT()
!                ENDIF
!                IF(MOVE)THEN
!                    THREE=MAX(PRT(PJ),PRT(PK))
!                    DO I=THREE-1,3,-1
!                        MOVE=DEP(I).EQ.O.OR.NINARR(DEP(I),PRT(THREE))
!                        IF(.NOT.MOVE)THEN
!                            EXIT
!                        ENDIF
!                    ENDDO
!                    IF(MOVE)THEN
!                    TMPDEP=DEP(THREE)
!                    TMPIND=PR(THREE)
!                    DO I=THREE,4,-1
!                        PR(I)=PR(I-1)
!                        DEP(I)=DEP(I-1)
!                    ENDDO
!                    DEP(3)=TMPDEP
!                    PR(3)=TMPIND
!                    CALL CALCPRT()
!                ENDIF
!            ENDIF
            TMPDEP=PI
            PI=0
            CALL INITPEN(PRT(TMPDEP),PRT(PJ),PRT(PK))
!         ENDIF
      ENDIF



      RETURN

  150 WRITE(6,'(A)')'IMPOSSIBLE DEPENDENCIES'
      STOP
      END

      SUBROUTINE CALCPRT()
      INCLUDE 'SIZES'
      COMMON /PERMUTE / PR,PRT
      INTEGER PR(NUMATM),PRT(NUMATM)
      DO I=1,NUMATM
        DO J=1,NUMATM
            IF(PR(J).EQ.I)PRT(I)=J
        END DO
      END DO
      END

      SUBROUTINE PERATMS(XYZ)
      INCLUDE 'SIZES'
      COMMON /PERMUTE / PR,PRT
      COMMON /GEOKST/ NATOMS,LABELS,NA(NUMATM),NB(NUMATM),NC(NUMATM)
      INTEGER PR(NUMATM),PRT(NUMATM),LABELS(NUMATM),TMPLBL(NUMATM),
     1         DEP(NUMATM)
      DOUBLE PRECISION XYZ(3,NUMATM), TMPXYZ(3,NUMATM)

      LOGICAL MOVE
      LOGICAL :: NINARR

      DO I=1,NUMATM
      TMPXYZ(:,I)=XYZ(:,PR(I))
      TMPLBL(I)=LABELS(PR(I))
      END DO
      XYZ(:,1:NUMATM)=TMPXYZ(:,1:NUMATM)
      LABELS(:)=TMPLBL(:)

      END


      SUBROUTINE UNPERATMS(XYZ)
      INCLUDE 'SIZES'
      COMMON /PERMUTE / PR,PRT
      COMMON /GEOKST/ NATOMS,LABELS,NA(NUMATM),NB(NUMATM),NC(NUMATM)
      INTEGER PR(NUMATM),PRT(NUMATM),LABELS(NUMATM),TMPLBL(NUMATM)
      DOUBLE PRECISION XYZ(3,NUMATM), TMPXYZ(3,NUMATM)

      DO I=1,NUMATM
      TMPXYZ(:,I)=XYZ(:,PRT(I))
      TMPLBL(I)=LABELS(PRT(I))
      END DO
      XYZ(:,1:NUMATM)=TMPXYZ(:,1:NUMATM)
      LABELS(:)=TMPLBL(:)

      END
