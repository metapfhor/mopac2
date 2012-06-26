      INTEGER FUNCTION NEWARR(L)

      INTEGER L
      INTEGER X(*)
      POINTER(P,X)
      DATA DESC /2/

      CALL RESIZEL(NEWARR,L,0)
      END

       SUBROUTINE PRINTARR(A)
       INTEGER A
       INTEGER N
       INTEGER :: GETLENGTH
       CALL PRINTARRL(A,GETLENGTH(A))

      END

      SUBROUTINE PRINTARRL(A,L)
          INTEGER A,L,Y
          INTEGER X(*)
          POINTER(P,X)
          P=A
          PRINT *,A
          DO I=1,L
          Y=X(I)
          PRINT *,I,': ', X(I)
          END DO
          RETURN
      END

      SUBROUTINE SETLENGTH(A, L)
        COMMON /ARRS / DESC
        INTEGER DESC
        INTEGER L,A,K
        INTEGER X(2)
        POINTER (P,X)
        P=A-DESC*SIZEOF(A)
        X(2)=MIN(X(1),L)
      END


      INTEGER FUNCTION GETLENGTH(A)
          COMMON /ARRS / DESC
          INTEGER DESC
          INTEGER A,L
          INTEGER X(2)
          POINTER (P,X)
          P=A-DESC*SIZEOF(A)
          GETLENGTH = X(2)
          RETURN
      END FUNCTION

      INTEGER FUNCTION GETSIZE(A)
          COMMON /ARRS / DESC
          INTEGER DESC
          INTEGER A
          INTEGER X(2)
          POINTER (P,X)
          P=A-DESC*SIZEOF(A)
          GETSIZE = X(1)
          RETURN
      END FUNCTION

      SUBROUTINE INIT(A)
      INTEGER A
      INTEGER :: GETLENGTH

      CALL INITL(A,GETLENGTH(A))
      END

      SUBROUTINE INITL(A,L)
      INTEGER A,L
      INTEGER :: GETSIZE
      INTEGER X(*)
      POINTER(P,X)
      P=A
      DO I=1,L
        X(I)=I
      END DO
      END


      SUBROUTINE DBLSIZE(A)
      INTEGER A
      INTEGER :: GETSIZE
      CALL RESIZE(A,GETSIZE(A)*2)

      END

      SUBROUTINE RESIZE(A,S)
      INTEGER A,S
      INTEGER :: GETLENGTH

      CALL RESIZEL(A,S,GETLENGTH(A))

      END

      SUBROUTINE RESIZEL(A,S,L)
      COMMON /ARRS / DESC
      INTEGER DESC
      INTEGER A,S
      INTEGER X(*),Y(*)
      POINTER(P,X),(Q,Y)

      IF(L.NE.0)THEN
          IF(L.NE.S)THEN
              P=A
              A=MALLOC((S+DESC)*SIZEOF(A))
              Q=A

              Y(1)=S
              Y(2)=L
              A=A+DESC*SIZEOF(A)
              Q=A

              CALL POINTERCOPY(P,1,MIN(S,L),Q,1,MIN(S,L))
              CALL DELETE(P-DESC*SIZEOF(A))
          ENDIF
      ELSE
      A=MALLOC((S+DESC)*SIZEOF(A))
      P=A
      X(1)=S
      X(2)=0
      A=A+DESC*SIZEOF(A)
      ENDIF
C      A(1:L+2)=Q(1:L+2)
      END

      SUBROUTINE DELETE(A)
      INTEGER  A
      CALL FREE(A)
      END

      SUBROUTINE DELETEL(A)

      END

      SUBROUTINE POINTERCOPY(A,I,J,B,K,L)
      INTEGER A,B
      INTEGER I,J,K,L
      INTEGER :: GETLENGTH

      CALL POINTERCOPYL(A,I,J,ABS(J-I),B,K,L,ABS(L-K))

      END

      SUBROUTINE POINTERCOPYL(A,I,J,SA,B,K,L,SB)
      INTEGER A,B,SA,SB
      INTEGER I,J,K,L
      INTEGER X(SA),Y(SB)
      POINTER(P,X),(Q,Y)
      P=A
      Q=B
      CALL ARRAYCOPYL(X,I,J,SA,Y,K,L,SB)
      END

      SUBROUTINE ARRAYCOPYL(A,I,J,SA,B,K,L,SB)
      INTEGER SA,SB
      INTEGER I,J,K,L
      INTEGER A(SA),B(SB)
      B(K:L)=A(I:J)
      END

      SUBROUTINE ADDVAL(A,V)
      INTEGER A,V
      INTEGER :: GETLENGTH
        CALL SETVAL(A,(GETLENGTH(A)+1),V)
      END

      SUBROUTINE ADDVALS(A,V)
      INTEGER :: GETLENGTH
      INTEGER A,V
        CALL SETVALS(A,GETLENGTH(A)+1,GETLENGTH(A)+GETLENGTH(V)+1,V)
      END

      SUBROUTINE SETVAL(A,I,V)
      INTEGER J,A,I,V
      INTEGER :: GETSIZE
      IF(GETSIZE(A).LT.I)THEN
      J=0
   10   IF(I.GT.2**J)THEN
          J=J+1
          GOTO 10
      ENDIF
      CALL RESIZE(A,2**J)
      ENDIF
      CALL SETVALCHKD(A,I,V)
      END

      SUBROUTINE SETVALCHKD(A,I,V)
          INTEGER A,I,V
          INTEGER X(*)
          POINTER(P,X)
          INTEGER :: GETLENGTH
          P=A
          IF(GETLENGTH(A).LT.I)CALL SETLENGTH(A,I)
          X(I)=V
      END

      INTEGER FUNCTION GETVAL(A,I)
          INTEGER :: GETLENGTH
          INTEGER :: GETVALCHKD
          INTEGER A,I
          IF(GETLENGTH(A).GE.I)THEN
            GETVAL=GETVALCHKD(A,I)
          ELSE
            GETVAL=0
          ENDIF
      END

      INTEGER FUNCTION GETVALCHKD(A,I)
          INTEGER A,I
          INTEGER X(*)
          POINTER(P,X)
          P=A
          GETVALCHKD=X(I)
      END


      SUBROUTINE GETVALS(A,I,J,RES)
      INTEGER A,I,J,RES
      INTEGER :: GETLENGTH
      INTEGER :: GETSIZE

      IF(GETLENGTH(A).GE.I.AND.GETLENGTH(A).GE.J.AND.I.LE.J)THEN
            CALL ENSURESIZE(RES,(J-I+1))
            CALL GETVALSCHCKD(A,I,J,RES)
          ELSE
             CALL ENSURESIZE(RES,ABS(J-I)+1)
             CALL ZEROL(RES,ABS(J-I)+1)
            GETVAL=0
          ENDIF
      END

      SUBROUTINE GETVALSCHCKD(A,I,J,RES)
      INTEGER A,I,J,RES,L
      INTEGER :: GETLENGTH
      L=GETLENGTH(RES)
      CALL POINTERCOPYL(A,I,J,GETLENGTH(A),RES,1,J-I+1,L)

      END

      RECURSIVE SUBROUTINE SETVALS(A,I,J,RES)
      INTEGER A,I,J,K,RES
      INTEGER :: GETLENGTH
      INTEGER :: GETSIZE

      IF(GETSIZE(A).GE.I.AND.GETSIZE(A).GE.J.AND.I.LE.J)THEN
            CALL ENSURESIZE(RES,J-I+1)

            CALL SETVALSCHCKD(A,I,J,RES)
      ELSE
            K=0
   10   IF(I.GT.2**J)THEN
          K=K+1
          GOTO 10
      ENDIF
      CALL RESIZE(A,2**K)
      CALL SETVALS(A,1,J,K)
      ENDIF
      END

      SUBROUTINE SETVALSCHCKD(A,I,J,RES)
      INTEGER A,I,J,RES,L,M
      INTEGER :: GETLENGTH
      L=GETLENGTH(RES)
      M=GETLENGTH(A)
      CALL POINTERCOPYL(RES,1,J,L,A,I,J,M)
      CALL SETLENGTH(A,MAX(M,J))
      END

      SUBROUTINE GETVALSCHCKDL(A,I,J,RES)

      END





      SUBROUTINE ZERO(A)
      INTEGER A
      INTEGER :: GETSIZE
      CALL ZEROL(A,GETSIZE(A))

      END

      SUBROUTINE ZEROL(A,L)
      INTEGER A, X(L)
      INTEGER, INTENT(IN) :: L
      POINTER(P,X)
      P=A
      DO I=1,L
        X(I)=0
      END DO

      END


      SUBROUTINE ENSURESIZE(A,L)
      INTEGER A
      INTEGER :: GETSIZE
      IF(GETSIZE(A).LT.L)CALL RESIZE(A,L)
      END
