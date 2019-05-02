      SUBROUTINE AVALUE(MAP,ARRAY,DIM,MARK)
C
C Routine maps values of "ARRAY" using the "MAP" function.
C Warning: MAP and ARRAY should be well defined, no checking
C is done. This routine assumes that ARRAY has INTEGER
C elements.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER DIM, MAP(*), ARRAY(*), MARK
C local
      INTEGER I
C begin
      DO I=1,DIM
      IF (ARRAY(I).EQ.MARK) THEN
      ARRAY(I)=MARK
      ELSE IF (ARRAY(I).GT.0) THEN
      ARRAY(I)=MAP(ARRAY(I))
      ELSE IF (ARRAY(I).LT.0) THEN
      ARRAY(I)=-MAP(-ARRAY(I))
      END IF
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE AINDR8(INVMAP,ARRAY,DIM,WORK)
C
C applies to DOUBLE PRECISION arrays
      IMPLICIT NONE
C input/output
      INTEGER DIM, INVMAP(*)
      DOUBLE PRECISION ARRAY(*), WORK(*)
C local
      INTEGER I
C begin
      DO I=1,DIM
      WORK(I)=ARRAY(INVMAP(I))
      END DO
      DO I=1,DIM
      ARRAY(I)=WORK(I)
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE AINDC8(INVMAP,ARRAY,DIM,WORK)
C
C applies to DOUBLE COMPLEX arrays
      IMPLICIT NONE
C input/output
      INTEGER DIM, INVMAP(*)
      DOUBLE COMPLEX ARRAY(*), WORK(*)
C local
      INTEGER I
C begin
      DO I=1,DIM
      WORK(I)=ARRAY(INVMAP(I))
      END DO
      DO I=1,DIM
      ARRAY(I)=WORK(I)
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE AINDX4(INVMAP,ARRAY,DIM,WORK)
C
C Same routine as AINDEX. Applies to INTEGER arrays.
      IMPLICIT NONE
C input/output
      INTEGER DIM, INVMAP(*), ARRAY(*), WORK(*)
C local
      INTEGER   I
C begin
      DO I=1,DIM
      WORK(I)=ARRAY(INVMAP(I))
      END DO
      DO I=1,DIM
      ARRAY(I)=WORK(I)
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE AINDL4(INVMAP,ARRAY,DIM,WORK)
C
C Same routine as AINDEX. Applies to LOGICAL arrays.
      IMPLICIT NONE
C input/output
      INTEGER DIM, INVMAP(*)
      LOGICAL ARRAY(*), WORK(*)
C local
      INTEGER   I
C begin
      DO I=1,DIM
      WORK(I)=ARRAY(INVMAP(I))
      END DO
      DO I=1,DIM
      ARRAY(I)=WORK(I)
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE AINDC4(INVMAP,ARRAY,DIM,WORK)
C
C Same routine as AINDEX. Applies to CHARACTER*4 arrays.
      IMPLICIT NONE
C input/output
      INTEGER DIM, INVMAP(*)
      CHARACTER*4 ARRAY(*), WORK(*)
C local
      INTEGER   I
C begin
      DO I=1,DIM
      WORK(I)=ARRAY(INVMAP(I))
      END DO
      DO I=1,DIM
      ARRAY(I)=WORK(I)
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE COPYI4(A,COPY,N)
C
C  Copies array A into COPY, applies to INTEGER arrays
      INTEGER A(*), COPY(*)
      INTEGER N, I
      IF (N.GT.0) THEN
      DO I=1,N
      COPY(I)=A(I)
      END DO
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE COPYR8(A,COPY,N)
C
C  Copies array A into COPY, applies to DOUBLE PRECISION arrays
      DOUBLE PRECISION A(*), COPY(*)
      INTEGER N, I
      IF (N.GT.0) THEN
      DO I=1,N
      COPY(I)=A(I)
      END DO
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE COPYC8(A,COPY,N)
C
C  Copies array A into COPY, applies to DOUBLE COMPLEX arrays
      DOUBLE COMPLEX A(*), COPY(*)
      INTEGER N, I
      IF (N.GT.0) THEN
      DO I=1,N
      COPY(I)=A(I)
      END DO
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE SWAPCO(DIM,X1,Y1,Z1,X2,Y2,Z2)
C
C Routine swaps coordinates X1,Y1,Z1 and X2,Y2,Z2 simultaneously
C applies to DOUBLE PRECISION arrays
C
      IMPLICIT NONE
      INTEGER DIM
      DOUBLE PRECISION  X1(*),Y1(*),Z1(*),X2(*),Y2(*),Z2(*)
C
      DOUBLE PRECISION  TEMP
      INTEGER I
C
      DO I=1,DIM
      TEMP=X1(I)
      X1(I)=X2(I)
      X2(I)=TEMP
      TEMP=Y1(I)
      Y1(I)=Y2(I)
      Y2(I)=TEMP
      TEMP=Z1(I)
      Z1(I)=Z2(I)
      Z2(I)=TEMP
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE COPYCO(DIM,X1,Y1,Z1,X2,Y2,Z2)
C
C routine copies coordinates X1,Y1,Z1 into X2,Y2,Z2
C applies to REAL arrays
C
      IMPLICIT NONE
      INTEGER DIM
      DOUBLE PRECISION  X1(*),Y1(*),Z1(*),X2(*),Y2(*),Z2(*)
C
      INTEGER I
      DO I=1,DIM
      X2(I)=X1(I)
      Y2(I)=Y1(I)
      Z2(I)=Z1(I)
      END DO
      RETURN
      END
C======================================================================
C search routines for looking up key's in (multiple) arrays
C =========================================================
C
      INTEGER FUNCTION FIND52(X1,X2,X3,X4,X5,Y1,Y2,Y3,Y4,Y5,DIM,WIDTH,
     &                        MARK)
C
C Linear search through X1(I),...,XWIDTH(I) for I=1,DIM to locate
C element Y1,...,YWIDTH. If not found FIND52 is set to MARK
C X1,...,X5 are of type INTEGER, Y1,...Y5 are of type INTEGER
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER DIM, X1(*), X2(*), X3(*), X4(*), X5(*)
      INTEGER Y1, Y2, Y3, Y4, Y5, WIDTH, MARK
C local
      INTEGER I
      LOGICAL FOUND
C begin
      FOUND=.FALSE.
      I=0
      DO WHILE (.NOT.(FOUND.OR.I.GE.DIM))
      I=I+1
      FOUND=(Y1.EQ.X1(I))
      IF (FOUND.AND.WIDTH.NE.1) THEN
      FOUND=(Y2.EQ.X2(I))
      IF (FOUND.AND.WIDTH.NE.2) THEN
      FOUND=(Y3.EQ.X3(I))
      IF (FOUND.AND.WIDTH.NE.3) THEN
      FOUND=(Y4.EQ.X4(I))
      IF (FOUND.AND.WIDTH.NE.4) THEN
      FOUND=(Y5.EQ.X5(I))
      END IF
      END IF
      END IF
      END IF
      END DO
      IF (FOUND) THEN
      FIND52=I
      ELSE
      FIND52=MARK
      END IF
      RETURN
      END
C======================================================================
      INTEGER FUNCTION SRCHWD(WORDS,NUMWRD,WORD)
C
C Searches the words array for word and returns its index
C if found and zero otherwise. Applies to INTEGER arrays.
      IMPLICIT NONE
C I/O
      INTEGER  WORDS(*), NUMWRD, WORD
C local
      INTEGER I
C begin
      SRCHWD = 0
      I = 1
      DO WHILE (SRCHWD .EQ. 0 .AND. I .LE. NUMWRD)
        IF (WORD .EQ. WORDS(I)) THEN
          SRCHWD = I
        ELSE
          I = I + 1
        END IF
      END DO
      RETURN
      END
C======================================================================
      INTEGER FUNCTION SRCHC4(WORD,ARRAY,START,STOP,MARK)
C
C Searches the words array for word and returns its index
C if found and MARK otherwise. Applies to character*4 arrays.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      CHARACTER*(*) WORD, ARRAY(*)
      INTEGER START, STOP, MARK
C begin
      SRCHC4 = START
C
      DO WHILE (SRCHC4 .LE. STOP)
        IF (WORD .EQ. ARRAY(SRCHC4)) RETURN
        SRCHC4 = SRCHC4 + 1
      END DO
C
      SRCHC4 = MARK
C
      RETURN
      END
C======================================================================
C  array initializing routines
C  ===========================
C
      SUBROUTINE FILL4(A,N,VALUE)
C
C Fills A with value, for I arrays
      IMPLICIT NONE
      INTEGER A(*), N, VALUE
C local
      INTEGER I
C begin
      IF (N.GT.0) THEN
      DO I=1,N
      A(I)=VALUE
      END DO
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE FILLR8(A,N,VALUE)
C
C Fills A with value, for R arrays
C
      IMPLICIT NONE
C
C input/ output
      DOUBLE PRECISION A(*), VALUE
      INTEGER N
C
C local
      INTEGER I
C
C begin
C
      IF (N.GT.0) THEN
        DO I=1,N
          A(I)=VALUE
        END DO
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE FILLR4(A,N,VALUE)
C
C Fills A with value, for single precision arrays
C
      IMPLICIT NONE
C
C input/ output
      REAL A(*), VALUE
      INTEGER N
C
C local
      INTEGER I
C
C begin
C
      IF (N.GT.0) THEN
        DO I=1,N
          A(I)=VALUE
        END DO
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE FILLC8(A,N,VALUE)
C
C Fills A with value, for C8 arrays
C
      IMPLICIT NONE
C
C input/ output
      DOUBLE COMPLEX A(*), VALUE
      INTEGER N
C
C local
      INTEGER I
C
C begin
C
      IF (N.GT.0) THEN
        DO I=1,N
          A(I)=VALUE
        END DO
      END IF
C
      RETURN
      END
C======================================================================
C  transformation of arrays
C
      SUBROUTINE MKFRAT(IMOVE,NATOM,FREEAT,NFREAT)
C
C Given a marker array IMOVE a new array FREEAT[1...NFREAT]
C is constructed which contains the indices of all elements
C with IMOVE(i).eq.0.
C
C Author: Robert Bruccoleri
C
      IMPLICIT NONE
C input/output
      INTEGER IMOVE(*), NATOM, FREEAT(*), NFREAT
C local
      INTEGER I
C begin
      NFREAT=0
      DO I=1,NATOM
      IF (IMOVE(I).EQ.0) THEN
      NFREAT=NFREAT+1
      FREEAT(NFREAT)=I
      END IF
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE MAKIND(INDEX,NATOM,NIND)
C
C Transforms the array INDEX into a new array which contains
C the indices of all elements (NIND) with INDEX(i).eq.1.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER INDEX(*), NATOM, NIND
C local
      INTEGER I
C begin
      NIND=0
      DO I=1,NATOM
      IF (INDEX(I).EQ.1) THEN
      NIND=NIND+1
      INDEX(NIND)=I
      END IF
      END DO
      RETURN
      END
C======================================================================
C  some vector manipulation routines
C  =================================
C
      SUBROUTINE VECADD(A,B,C,DIM)
C
C     C = A+B
C
      DOUBLE PRECISION A(*), B(*), C(*)
      INTEGER DIM, I
      DO I=1,DIM
      C(I)=A(I)+B(I)
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE VECMUL(A,B,C,DIM)
C
C     C(i) = B(i)*C(i)
C
      DOUBLE PRECISION A(*), B(*), C(*)
      INTEGER DIM, I
      DO I=1,DIM
      C(I)=A(I)*B(I)
      END DO
      RETURN
      END
