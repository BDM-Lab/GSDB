      SUBROUTINE SORTP(N,PERM,ORDERL,A1,A2,A3,A4,A5,A6,A7,A8)
C
C Generalised quicksort routine
C PERMutation version
C
C Structure A1,A2,A3,A4,A5,A6 of length N sorted as defined in logical
C function ORDERL. The algorithm used is quicksort.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER N, PERM(*)
      LOGICAL ORDERL
      EXTERNAL ORDERL
      INTEGER A1(*), A2(*), A3(*), A4(*), A5(*), A6(*), A7(*), A8(*)
C local
      INTEGER STORE
C begin
      IF (N.NE.0) THEN
      STORE=ALLHP(INTEG4(INT(ALOG(FLOAT(N))*4)))
      CALL SORTP2(N,PERM,ORDERL,A1,A2,A3,A4,A5,A6,A7,A8,HEAP(STORE))
      CALL FREHP(STORE,INTEG4(INT(ALOG(FLOAT(N))*4)))
      END IF
      RETURN
      END
C
      SUBROUTINE SORTP2(N,PERM,ORDER,A1,A2,A3,A4,A5,A6,A7,A8,STORE)
C
C See SORTP above.
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER N, PERM(*)
      LOGICAL ORDER
      EXTERNAL ORDER
      INTEGER A1(*), A2(*), A3(*), A4(*), A5(*), A6(*), A7(*), A8(*)
      INTEGER STORE(*)
C local
      INTEGER I, END, END1, L1, L2, START, START2, TOP, TEMP
      DOUBLE PRECISION DUMMY
      INTEGER RNDIND
      LOGICAL CONDIT, CLOOP1, CLOOP2
C begin
      DO I=1,N
      PERM(I)=I
      END DO
C
      TOP=0
      START=1
      END=N
C
3333  IF(END.EQ.START) GOTO 1111
      IF(END.NE.START+1) GOTO 2222
      IF (.NOT. ORDER(PERM(START), PERM(END),
     &                A1,A2,A3,A4,A5,A6,A7,A8)) THEN
      TEMP=PERM(START)
      PERM(START)=PERM(END)
      PERM(END)=TEMP
      END IF
1111  IF(TOP.EQ.0) GOTO 9999
      END=STORE(TOP)
      START=STORE(TOP-1)
      TOP=TOP-2
      GOTO 3333
C
2222  CONTINUE
C
C==begin=process=partition
      CALL GGUBFS(DUMMY)
      RNDIND=INT((END-START)*DUMMY)+START
      START2=START-1
      END1=END+1
      CLOOP1=.TRUE.
      DO WHILE (CLOOP1)
C
      CLOOP2=.TRUE.
      DO WHILE (CLOOP2)
      START2=START2+1
      CONDIT=ORDER(PERM(START2),PERM(RNDIND),A1,A2,A3,A4,A5,A6,A7,A8)
      IF (.NOT.CONDIT .OR. START2.EQ.END) CLOOP2=.FALSE.
      END DO
C
      CLOOP2=.TRUE.
      DO WHILE (CLOOP2)
      END1=END1-1
      CONDIT=ORDER(PERM(RNDIND),PERM(END1),A1,A2,A3,A4,A5,A6,A7,A8)
      IF (.NOT.CONDIT .OR. END1.EQ.START) CLOOP2=.FALSE.
      END DO
C
      IF (START2.LT.END1) THEN
      TEMP=PERM(START2)
      PERM(START2)=PERM(END1)
      PERM(END1)=TEMP
      END IF
      IF (START2.GE.END1) CLOOP1=.FALSE.
      END DO
C
      IF (START2.LT.RNDIND) THEN
      TEMP=PERM(START2)
      PERM(START2)=PERM(RNDIND)
      PERM(RNDIND)=TEMP
      START2=START2+1
      ELSE
      IF (RNDIND.LT.END1) THEN
      TEMP=PERM(END1)
      PERM(END1)=PERM(RNDIND)
      PERM(RNDIND)=TEMP
      END1=END1-1
      END IF
      END IF
C==end=process=partition
C
      L1=END1-START
      L2=END-START2
      IF (L1.LE.L2) THEN
C
C     L2>=L1 DO L1 FIRST
C
      TOP=TOP+2
      STORE(TOP-1)=START2
      STORE(TOP)=END
      END=END1
      ELSE
      TOP=TOP+2
      STORE(TOP-1)=START
      STORE(TOP)=END1
      START=START2
      END IF
      GOTO 3333
C
9999  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE SORT(N,EXCH,ORDERL,A1,A2,A3,A4,A5,A6,A7,A8)
C
C Generalised quicksort routine
C EXCHange version
C
C Structure A1,A2,A3,A4,A5,A6,A7,A8 of length N sorted as defined
C in logical function ORDERL. The algorithm used is quicksort.
C EXCHange of elements defined in subroutine EXCH.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INTEGER N
      EXTERNAL EXCH
      LOGICAL ORDERL
      EXTERNAL ORDERL
      INTEGER A1(*), A2(*), A3(*), A4(*), A5(*), A6(*), A7(*), A8(*)
C local
      INTEGER   STORE
C begin
      IF (N.NE.0) THEN
      STORE=ALLHP(INTEG4(INT(ALOG(FLOAT(N))*4)))
      CALL SORT2(N,EXCH,ORDERL,A1,A2,A3,A4,A5,A6,A7,A8,HEAP(STORE))
      CALL FREHP(STORE,INTEG4(INT(ALOG(FLOAT(N))*4)))
      END IF
C
      RETURN
      END
C
      SUBROUTINE SORT2(N,EXCH,ORDER,A1,A2,A3,A4,A5,A6,A7,A8,STORE)
C
C See SORT above.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER N
      EXTERNAL EXCH
      LOGICAL ORDER
      EXTERNAL ORDER
      INTEGER A1(*), A2(*), A3(*), A4(*), A5(*), A6(*), A7(*), A8(*)
      INTEGER STORE(*)
C local
      INTEGER END, END1, L1, L2, START, START2, TOP
      DOUBLE PRECISION DUMMY
      INTEGER RNDIND
      LOGICAL CONDIT, CLOOP1, CLOOP2
C begin
      TOP=0
      START=1
      END=N
C
3333  IF(END.EQ.START) GOTO 1111
      IF(END.NE.START+1) GOTO 2222
      IF (.NOT. ORDER(START,END,A1,A2,A3,A4,A5,A6,A7,A8)) THEN
      CALL EXCH(START,END,A1,A2,A3,A4,A5,A6,A7,A8)
      END IF
1111  IF(TOP.EQ.0) GOTO 9999
      END=STORE(TOP)
      START=STORE(TOP-1)
      TOP=TOP-2
      GOTO 3333
C
2222  CONTINUE
C
C process partition
      CALL GGUBFS(DUMMY)
      RNDIND=INT((END-START)*DUMMY)+START
      START2=START-1
      END1=END+1
      CLOOP1=.TRUE.
      DO WHILE (CLOOP1)
C
      CLOOP2=.TRUE.
      DO WHILE (CLOOP2)
      START2=START2+1
      CONDIT=ORDER(START2,RNDIND,A1,A2,A3,A4,A5,A6,A7,A8)
      IF (.NOT.CONDIT .OR. START2.EQ.END) CLOOP2=.FALSE.
      END DO
C
      CLOOP2=.TRUE.
      DO WHILE (CLOOP2)
      END1=END1-1
      CONDIT=ORDER(RNDIND,END1,A1,A2,A3,A4,A5,A6,A7,A8)
      IF (.NOT.CONDIT .OR. END1.EQ.START) CLOOP2=.FALSE.
      END DO
C
      IF (START2.LT.END1) THEN
      CALL EXCH(START2,END1,A1,A2,A3,A4,A5,A6,A7,A8)
      END IF
      IF (START2.GE.END1) CLOOP1=.FALSE.
      END DO
C
      IF (START2.LT.RNDIND) THEN
      CALL EXCH(START2,RNDIND,A1,A2,A3,A4,A5,A6,A7,A8)
      ELSE
      IF (RNDIND.LT.END1) THEN
      CALL EXCH(END1,RNDIND,A1,A2,A3,A4,A5,A6,A7,A8)
      END IF
      END IF
C end process partition
C
      L1=END1-START
      L2=END-START2
      IF (L1.LE.L2) THEN
C
C     L2>=L1 DO L1 FIRST
C
      TOP=TOP+2
      STORE(TOP-1)=START2
      STORE(TOP)=END
      END=END1
      ELSE
      TOP=TOP+2
      STORE(TOP-1)=START
      STORE(TOP)=END1
      START=START2
      END IF
      GOTO 3333
C
9999  CONTINUE
C
      RETURN
      END
C
      LOGICAL FUNCTION ORDER5(K,L,X1,X2,X3,X4,X5,X6,X7,WIDTH)
C
C True if X1(K) <= X1(L),...,XWIDTH(K) <= XWIDTH(L), K <= L.
C (in canonical order)
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER K, L
      INTEGER X1(*), X2(*), X3(*), X4(*), X5(*), X6(*), X7(*)
      INTEGER WIDTH
C begin
      ORDER5=(X1(K).LE.X1(L))
      IF (WIDTH.NE.1 .AND. X1(K).EQ.X1(L)) THEN
      ORDER5=(X2(K).LE.X2(L))
      IF (WIDTH.NE.2 .AND. X2(K).EQ.X2(L)) THEN
      ORDER5=(X3(K).LE.X3(L))
      IF (WIDTH.NE.3 .AND. X3(K).EQ.X3(L)) THEN
      ORDER5=(X4(K).LE.X4(L))
      IF (WIDTH.NE.4 .AND. X4(K).EQ.X4(L)) THEN
      ORDER5=(X5(K).LE.X5(L))
      IF (WIDTH.NE.5 .AND. X5(K).EQ.X5(L)) THEN
      ORDER5=(X6(K).LE.X6(L))
      IF (WIDTH.NE.6 .AND. X6(K).EQ.X6(L)) THEN
      ORDER5=(X7(K).LE.X7(L))
      IF (WIDTH.EQ.7.AND.X7(K).EQ.X7(L)) ORDER5=(K.LE.L)
      END IF
      IF (WIDTH.EQ.6.AND.X6(K).EQ.X6(L)) ORDER5=(K.LE.L)
      END IF
      IF (WIDTH.EQ.5.AND.X5(K).EQ.X5(L)) ORDER5=(K.LE.L)
      END IF
      IF (WIDTH.EQ.4.AND.X4(K).EQ.X4(L)) ORDER5=(K.LE.L)
      END IF
      IF (WIDTH.EQ.3.AND.X3(K).EQ.X3(L)) ORDER5=(K.LE.L)
      END IF
      IF (WIDTH.EQ.2.AND.X2(K).EQ.X2(L)) ORDER5=(K.LE.L)
      END IF
      IF (WIDTH.EQ.1.AND.X1(K).EQ.X1(L)) ORDER5=(K.LE.L)
      RETURN
      END
C
      SUBROUTINE EXCH5(INDEX1,INDEX2,A1,A2,A3,A4,A5,A6,A7,WIDTH)
C
C Exchanges A1(INDEX1) and A1(INDEX2),...,
C AWIDTH(INDEX1) and AWIDTH(INDEX2) .
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER INDEX1, INDEX2
      INTEGER A1(*), A2(*), A3(*), A4(*), A5(*), A6(*), A7(*)
      INTEGER WIDTH
C local
      INTEGER EXCHG
C begin
      EXCHG=A1(INDEX1)
      A1(INDEX1)=A1(INDEX2)
      A1(INDEX2)=EXCHG
      IF (WIDTH.NE.1) THEN
      EXCHG=A2(INDEX1)
      A2(INDEX1)=A2(INDEX2)
      A2(INDEX2)=EXCHG
      IF (WIDTH.NE.2) THEN
      EXCHG=A3(INDEX1)
      A3(INDEX1)=A3(INDEX2)
      A3(INDEX2)=EXCHG
      IF (WIDTH.NE.3) THEN
      EXCHG=A4(INDEX1)
      A4(INDEX1)=A4(INDEX2)
      A4(INDEX2)=EXCHG
      IF (WIDTH.NE.4) THEN
      EXCHG=A5(INDEX1)
      A5(INDEX1)=A5(INDEX2)
      A5(INDEX2)=EXCHG
      IF (WIDTH.NE.5) THEN
      EXCHG=A6(INDEX1)
      A6(INDEX1)=A6(INDEX2)
      A6(INDEX2)=EXCHG
      IF (WIDTH.NE.6) THEN
      EXCHG=A7(INDEX1)
      A7(INDEX1)=A7(INDEX2)
      A7(INDEX2)=EXCHG
      END IF
      END IF
      END IF
      END IF
      END IF
      END IF
      RETURN
      END
C
      LOGICAL FUNCTION ORDER(K,L,ARRAY,WIDTH,D1,D2,D3,D4,D5,D6)
C
C SORT order function for INTEGER's
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER WIDTH
      INTEGER K, L, ARRAY(WIDTH,*), D1, D2, D3, D4, D5, D6
C local
      INTEGER I
      LOGICAL CLOOP
C begin
      I=0
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      I=I+1
      ORDER=(ARRAY(I,K).LE.ARRAY(I,L))
      IF (ARRAY(I,K).NE.ARRAY(I,L).OR.I.GE.WIDTH) CLOOP=.FALSE.
      END DO
      IF (I.EQ.WIDTH.AND.ARRAY(I,K).EQ.ARRAY(I,L)) ORDER=K.LE.L
C
      RETURN
      END
C
      SUBROUTINE EXCH(K,L,ARRAY,WIDTH,D1,D2,D3,D4,D5,D6)
C
C SORT exchange function for INTEGER's
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER WIDTH
      INTEGER K, L, ARRAY(WIDTH,*), D1, D2, D3, D4, D5, D6
C begin
      INTEGER I, ITEMP
C local
      DO I=1,WIDTH
      ITEMP=ARRAY(I,K)
      ARRAY(I,K)=ARRAY(I,L)
      ARRAY(I,L)=ITEMP
      END DO
      RETURN
      END
C
      LOGICAL FUNCTION ORDERR(K,L,ARRAY,WIDTH,D1,D2,D3,D4,D5,D6)
C
C SORT order function for REAL's
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER WIDTH
      INTEGER K, L
      DOUBLE PRECISION ARRAY(WIDTH,*), D1, D2, D3, D4, D5, D6
C local
      INTEGER   I
      LOGICAL CLOOP
C begin
      I=0
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      I=I+1
      ORDERR=(ARRAY(I,K).LE.ARRAY(I,L))
      IF (ARRAY(I,K).NE.ARRAY(I,L).OR.I.GE.WIDTH) CLOOP=.FALSE.
      END DO
      IF (I.EQ.WIDTH.AND.ARRAY(I,K).EQ.ARRAY(I,L)) ORDERR=K.LE.L
      RETURN
      END
C
      LOGICAL FUNCTION ORDRRR(K,L,ARRAY,D1,D2,D3,D4,D5,D6,D7)
C
C SORT order function for REAL's (maintains original order in cases
C of approximate equality)
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'consta.inc'
      INTEGER K, L
      DOUBLE PRECISION ARRAY(*), D1, D2, D3, D4, D5, D6, D7
C parameter
      DOUBLE PRECISION SMALL
      PARAMETER (SMALL=R4SMAL*0.01D0)
C begin
      IF (ABS(ARRAY(K)-ARRAY(L)).LT.SMALL) THEN
      ORDRRR=K.LE.L
      ELSE
      ORDRRR=ARRAY(K).LE.ARRAY(L)
      END IF
      RETURN
      END
C
      SUBROUTINE EXCHR(K,L,ARRAY,WIDTH,D1,D2,D3,D4,D5,D6)
C
C SORT exchange routine for REAL's
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER WIDTH
      INTEGER K, L
      DOUBLE PRECISION ARRAY(WIDTH,*), D1, D2, D3, D4, D5, D6
C local
      INTEGER I
      DOUBLE PRECISION ITEMP
C begin
      DO I=1,WIDTH
      ITEMP=ARRAY(I,K)
      ARRAY(I,K)=ARRAY(I,L)
      ARRAY(I,L)=ITEMP
      END DO
      RETURN
      END
