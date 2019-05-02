C=====================================
C c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
C     SUBROUTINE GIVEIS(N,NVECT,NV,A,B,IND,ROOT,VECT,IERR)
C
C     EISPACK-BASED SUBSTITUTE FOR QCPE SUBROUTINE GIVENS.
C     FINDS ALL EIGENVALUES AND SOME EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX.   AUTHOR.. C. MOLER AND D. SPANGLER, N.R.C.C., 4/1/79.
C
C     INPUT..
C     N     = ORDER OF MATRIX .
C     NVECT = NUMBER OF VECTORS DESIRED.  0 .LE. NVECT .LE. N .
C     NV    = ROW DIMENSION OF VECT .
C     A     = INPUT MATRIX, COLUMNS OF THE UPPER TRIANGLE PACKED INTO
C             LINEAR ARRAY OF DIMENSION N*(N+1)/2 .
C     B     = SCRATCH ARRAY, 9*N ELEMENTS (NOTE THIS IS MORE THAN
C             PREVIOUS VERSIONS OF GIVENS.)
C
C     OUTPUT..
C     A       DESTORYED .
C     ROOT  = ALL EIGENVALUES, ROOT(1) .LE. ... .LE. ROOT(N) .
C             (FOR OTHER ORDERINGS, SEE BELOW.)
C     VECT  = EIGENVECTORS FOR ROOT(1),..., ROOT(NVECT) .
C     IERR  = 0 IF NO ERROR DETECTED,
C           = K IF ITERATION FOR K-TH EIGENVALUE FAILED,
C           = -K IF ITERATION FOR K-TH EIGENVECTOR FAILED.
C             (FAILURES SHOULD BE VERY RARE.  CONTACT MOLER.)
C
C     CALLS MODIFIED EISPACK SUBROUTINES TRED3B, IMTQLV, TINVTB, AND
C     TRBK3B.  THE SUBROUTINES TRED3B, TINVTB, AND TRBK3B.
C     THE ORIGINAL EISPACK SUBROUTINES TRED3, TINVIT, AND TRBAK3
C     WERE MODIFIED BY THE INTRODUCTION OF TWO SUBROUTINES FROM THE
C     BLAS LIBRARY - DDOT AND DAXPY.
C     THESE WERE MOVED IN-LINE TO REDUCE SUBROUTINE CALL OVERHEAD
C       -B.A. BORGIAS
C
C         IF TINVIT FAILS TO CONVERGE, TQL2 IS CALLED
C
C         SEE EISPACK USERS GUIDE, B. T. SMITH ET AL, SPRINGER-VERLAG
C     LECTURE NOTES IN COMPUTER SCIENCE, VOL. 6, 2-ND EDITION, 1976 .
C     NOTE THAT IMTQLV AND TINVTB HAVE INTERNAL MACHINE
C     DEPENDENT CONSTANTS.
C
Cc c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
      SUBROUTINE GIVEIS(N,NVECT,NV,A,B,IND,ROOT,VECT,IERR)
      IMPLICIT NONE
      INTEGER N, NVECT, NV, IERR
      DOUBLE PRECISION A(*),B(N,8),ROOT(N),VECT(NV,*)
      DOUBLE PRECISION ONE, ZERO
      INTEGER IND(N)
      INTEGER I, J
      DATA ONE, ZERO /1.0D+00, 0.0D+00/
      CALL TRED3B(N,(N*N+N)/2,A,B(1,1),B(1,2),B(1,3))
C       IN THE LISTING I RECEIVED,THE ELEMENT B(1,9) WAS
C       REPLACED BY IND
      CALL IMTQLV(N,B(1,1),B(1,2),B(1,3),ROOT,IND,IERR,B(1,4))
C
C
      IF (IERR .NE. 0) RETURN
C
C     TO REORDER ROOTS...
C     K = N/2
C     B(1,3) = 2.0D+00
C     DO 50 I = 1, K
C        J = N+1-I
C        T = ROOT(I)
C        ROOT(I) = ROOT(J)
C        ROOT(J) = T
C 50  CONTINUE
C
      IF (NVECT .LE. 0) RETURN
      CALL TINVTB(NV,N,B(1,1),B(1,2),B(1,3),NVECT,ROOT,IND,VECT,IERR,
     +     B(1,4),B(1,5),B(1,6),B(1,7),B(1,8))
C
      IF (IERR .EQ. 0) GO TO 160
C
C      IF INVERSE ITERATION GIVES AN ERROR IN DETERMINING THE
C      EIGENVECTORS, TRY THE QL ALGORITHM IF ALL THE EIGENVECTORS
C      ARE DESIRED.
C
      IF (NVECT .NE. N) RETURN
      DO 120 I = 1, NVECT
      DO 100 J = 1, N
      VECT(I,J) = ZERO
  100 CONTINUE
      VECT(I,I) = ONE
  120 CONTINUE
      CALL TQL2 (NV,N,B(1,1),B(1,2),VECT,IERR)
C
      DO 140 I = 1, NVECT
      ROOT(I) = B(I,1)
  140 CONTINUE
      IF (IERR .NE. 0) RETURN
  160 CALL TRBK3B(NV,N,(N*N+N)/2,A,NVECT,VECT)
      RETURN
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE IMTQLV(N,D,E,E2,W,IND,IERR,RV1)
      IMPLICIT NONE
C
      INCLUDE 'machvar.inc'
      INTEGER I,J,K,L,M,N,II,MML,TAG,IERR
      DOUBLE PRECISION D(N),E(N),E2(N),W(N),RV1(N)
      DOUBLE PRECISION B,C,F,G,P,R,S
C     DOUBLE PRECISION DSQRT,DABS,DSIGN
      INTEGER IND(N)
C
C     THIS SUBROUTINE IS A VARIANT OF  IMTQL1  WHICH IS A TRANSLATION
C     OF
C     ALGOL PROCEDURE IMTQL1, NUM. MATH. 12, 377-383(1968) BY MARTIN
C     AND
C     WILKINSON, AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC TRIDIAGONAL
C     MATRIX BY THE IMPLICIT QL METHOD AND ASSOCIATES WITH THEM
C     THEIR CORRESPONDING SUBMATRIX INDICES.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C     ON OUTPUT-
C
C        D AND E ARE UNALTERED,
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO,
C
C        W CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES,
C
C        IND CONTAINS THE SUBMATRIX INDICES ASSOCIATED WITH THE
C          CORRESPONDING EIGENVALUES IN W -- 1 FOR EIGENVALUES
C          BELONGING TO THE FIRST SUBMATRIX FROM THE TOP,
C          2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS,
C
C        RV1 IS A TEMPORARY STORAGE ARRAY.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C ------------------------------------------------------------------
C
      IERR = 0
      K = 0
      TAG = 0
C
C      DO 100 I = 1, N
C      W(I) = D(I)
C      IF (I .NE. 1) RV1(I-1) = E(I)
C  REPLACE WITH VECTORIZED LOOP
      W(1) = D(1)
      DO 100 I = 2, N
      W(I) = D(I)
      RV1(I-1) = E(I)
  100 CONTINUE
C
      E2(1) = 0.0D+00
      RV1(N) = 0.0D+00
C
      DO 360 L = 1, N
      J = 0
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
  120 DO 140 M = L, N-1
C     IF (M .EQ. N) GO TO 160
      IF (DABS(RV1(M)) .LE. FPEPS * (DABS(W(M)) + DABS(W(M+1)))) GO TO
     +     160
C     ********** GUARD AGAINST UNDERFLOWED ELEMENT OF E2 **********
      IF (E2(M+1) .EQ. 0.0D+00) GO TO 180
  140 CONTINUE
      M = N
C
  160 IF (M .LE. K) GO TO 200
      IF (M .NE. N) E2(M+1) = 0.0D+00
  180 K = M
      TAG = TAG + 1
  200 P = W(L)
      IF (M .EQ. L) GO TO 280
      IF (J .EQ. 30) GO TO 380
      J = J + 1
C     ********** FORM SHIFT **********
      G = (W(L+1) - P) / (2.0D+00 * RV1(L))
      R = DSQRT(G*G+1.0D+00)
      G = W(M) - P + RV1(L) / (G + DSIGN(R,G))
      S = 1.0D+00
      C = 1.0D+00
      P = 0.0D+00
      MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
      DO 260 II = 1, MML
      I = M - II
      F = S * RV1(I)
      B = C * RV1(I)
      IF (DABS(F) .LT. DABS(G)) GO TO 220
      C = G / F
      R = DSQRT(C*C+1.0D+00)
      RV1(I+1) = F * R
      S = 1.0D+00 / R
      C = C * S
      GO TO 240
  220 S = F / G
      R = DSQRT(S*S+1.0D+00)
      RV1(I+1) = G * R
      C = 1.0D+00 / R
      S = S * C
  240 G = W(I+1) - P
      R = (W(I) - G) * S + 2.0D+00 * C * B
      P = S * R
      W(I+1) = G + P
      G = C * R - B
  260 CONTINUE
C
      W(L) = W(L) - P
      RV1(L) = G
      RV1(M) = 0.0D+00
      GO TO 120
C     ********** ORDER EIGENVALUES **********
  280 IF (L .EQ. 1) GO TO 320
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
      DO 300 II = 2, L
      I = L + 2 - II
      IF (P .GE. W(I-1)) GO TO 340
      W(I) = W(I-1)
      IND(I) = IND(I-1)
  300 CONTINUE
C
  320 I = 1
  340 W(I) = P
      IND(I) = TAG
  360 CONTINUE
C
      GO TO 400
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
  380 IERR = L
  400 RETURN
C     ********** LAST CARD OF IMTQLV **********
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE TINVTB(NM,N,D,E,E2,M,W,IND,Z,
     +                  IERR,RV1,RV2,RV3,RV4,RV6)
      IMPLICIT NONE
C
      INCLUDE 'machvar.inc'
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),Z(NM,M),
     +       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      DOUBLE PRECISION U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,NORM2,ORDER
      INTEGER IND(M)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-
C     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
C     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
C          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
C          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
C          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
C          0.0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0
C          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
C          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
C          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE,
C
C        M IS THE NUMBER OF SPECIFIED EIGENVALUES,
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER,
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
C
C     ON OUTPUT-
C
C        ALL INPUT ARRAYS ARE UNALTERED,
C
C        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
C          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS,
C
C        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C   ------------------------------------------------------------------
C
      INTEGER  IDOT, IXPY
C
      IERR = 0
      IF (M .EQ. 0) GO TO 680
      TAG = 0
      ORDER = 1.0D+00 - E2(1)
      Q = 0
C     ********** ESTABLISH AND PROCESS NEXT SUBMATRIX **********
  100 P = Q + 1
C
      DO 120 Q = P, N
      IF (Q .EQ. N) GO TO 140
      IF (E2(Q+1) .EQ. 0.0D+00) GO TO 140
  120 CONTINUE
C     ********** FIND VECTORS BY INVERSE ITERATION **********
  140 TAG = TAG + 1
      S = 0
C
      DO 660 R = 1, M
      IF (IND(R) .NE. TAG) GO TO 660
      ITS = 1
      X1 = W(R)
      IF (S .NE. 0) GO TO 220
C     ********** CHECK FOR ISOLATED ROOT **********
      XU = 1.0D+00
      IF (P .NE. Q) GO TO 160
      RV6(P) = 1.0D+00
      GO TO 600
  160 NORM = DABS(D(P))
      IP = P + 1
C
      DO 180 I = IP, Q
  180 NORM = NORM + DABS(D(I)) + DABS(E(I))
C     ********** EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW **********
      NORM2 = NORM / (Q - P + 1)
      EPS2 = 1.0D-3 * NORM
      EPS3 = FPEPS * NORM
      UK = DFLOAT(Q-P+1)
      EPS4 = UK * EPS3
      UK = EPS4 / DSQRT(UK)
      S = P
  200 GROUP = 0
      GO TO 240
C     ********** LOOK FOR CLOSE OR COINCIDENT ROOTS **********
  220 IF (DABS(X1-X0) .GE. EPS2) GO TO 200
      GROUP = GROUP + 1
      IF (ORDER * (X1 - X0) .LE. 0.0D+00) X1 = X0 + ORDER * EPS3
C     ********** ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR **********
  240 V = 0.0D+00
C
      DO 300 I = P, Q
      RV6(I) = UK
      IF (I .EQ. P) GO TO 280
      IF (DABS(E(I)) .LT. DABS(U)) GO TO 260
C     ********** WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY **********
      XU = U / E(I)
      RV4(I) = XU
      RV1(I-1) = E(I)
      RV2(I-1) = D(I) - X1
      RV3(I-1) = 0.0D+00
      IF (I .NE. Q) RV3(I-1) = E(I+1)
      U = V - XU * RV2(I-1)
      V = -XU * RV3(I-1)
      GO TO 300
  260 XU = E(I) / U
      RV4(I) = XU
      RV1(I-1) = U
      RV2(I-1) = V
      RV3(I-1) = 0.0D+00
  280 U = D(I) - X1 - XU * V
      IF (I .NE. Q) V = E(I+1)
  300 CONTINUE
C
      IF (U .EQ. 0.0D+00) U = EPS3
      RV1(Q) = U
      RV2(Q) = 0.0D+00
      RV3(Q) = 0.0D+00
C     ********** BACK SUBSTITUTION
C                FOR I=Q STEP -1 UNTIL P DO -- **********
  320 DO 340 II = P, Q
      I = P + Q - II
      RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
      V = U
      U = RV6(I)
  340 CONTINUE
C     ********** ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP **********
      IF (GROUP .EQ. 0) GO TO 400
      J = R
C
      DO 380 JJ = 1, GROUP
  360 J = J - 1
      IF (IND(J) .NE. TAG) GO TO 360
C Try replacing the function calls to ddot and daxpy with in-line loops
C     XU = DDOT(Q-P+1,RV6(P),1,Z(P,J),1)
      XU = 0.0
      do 5000 idot=0,Q-P
      XU = XU + RV6(P+IDOT)*Z(P+IDOT,J)
5000  continue
      IF(XU.NE.0.0) then
C         CALL DAXPY(Q-P+1,-XU,Z(P,J),1,RV6(P),1)
          do 5001 ixpy = 0,Q-P
          RV6(P+ixpy) = RV6(P+ixpy) - XU*Z(P+ixpy,J)
5001      continue
      ENDIF
C
  380 CONTINUE
C
  400 NORM = 0.0D+00
C
      DO 420 I = P, Q
  420 NORM = NORM + DABS(RV6(I))
C
      IF (NORM .GE. 1.0D+00) GO TO 560
C     ********** FORWARD SUBSTITUTION **********
      IF (ITS .EQ. 5) GO TO 540
      IF (NORM .NE. 0.0D+00) GO TO 440
      RV6(S) = EPS4
      S = S + 1
      IF (S .GT. Q) S = P
      GO TO 480
  440 XU = EPS4 / NORM
C
      DO 460 I = P, Q
  460 RV6(I) = RV6(I) * XU
C     ********** ELIMINATION OPERATIONS ON NEXT VECTOR
C                ITERATE **********
  480 DO 520 I = IP, Q
      U = RV6(I)
C     ********** IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                WAS PERFORMED EARLIER IN THE
C                TRIANGULARIZATION PROCESS **********
      IF (RV1(I-1) .NE. E(I)) GO TO 500
      U = RV6(I-1)
      RV6(I-1) = RV6(I)
  500 RV6(I) = U - RV4(I) * RV6(I-1)
  520 CONTINUE
C
      ITS = ITS + 1
      GO TO 320
C     ********** SET ERROR -- NON-CONVERGED EIGENVECTOR **********
  540 IERR = -R
      XU = 0.0D+00
      GO TO 600
C     ********** NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER **********
  560 U = 0.0D+00
C
C included the following change for zero eigenvalues (ATB)
      IF (ABS(X1).LT.FPEPS) THEN
      DO 5555 I = 1, N
 5555 RV6(I) = 0.0D+00
      RV6(R)=1.0D+00
      XU=1.0D+00
      ELSE
      DO 580 I = P, Q
  580 U = U + RV6(I)**2
C
      XU = 1.0D+00 / DSQRT(U)
      END IF
C
  600 DO 620 I = 1, N
  620 Z(I,R) = 0.0D+00
C
      DO 640 I = P, Q
  640 Z(I,R) = RV6(I) * XU
C
      X0 = X1
  660 CONTINUE
C
      IF (Q .LT. N) GO TO 100
  680 RETURN
C     ********** LAST CARD OF TINVIT **********
      END
C
C   ------------------------------------------------------------------
C
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
      IMPLICIT NONE
C
      INCLUDE 'machvar.inc'
      INTEGER I,J,K,L,M,N,II,L1,NM,MML,IERR
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
      DOUBLE PRECISION B,C,F,G,H,P,R,S
C     DOUBLE PRECISION SQRT,ABS,SIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1,
C
C        E HAS BEEN DESTROYED,
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C   ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 400
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0D+00
      B = 0.0D+00
      E(N) = 0.0D+00
C
      DO 300 L = 1, N
      J = 0
      H = FPEPS * (DABS(D(L)) + DABS(E(L)))
      IF (B .LT. H) B = H
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
      DO 120 M = L, N
      IF (DABS(E(M)) .LE. B) GO TO 140
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
  120 CONTINUE
C
  140 IF (M .EQ. L) GO TO 280
  160 IF (J .EQ. 30) GO TO 380
      J = J + 1
C     ********** FORM SHIFT **********
      L1 = L + 1
      G = D(L)
      P = (D(L1) - G) / (2.0D+00 * E(L))
      R = DSQRT(P*P+1.0D+00)
      D(L) = E(L) / (P + DSIGN(R,P))
      H = G - D(L)
C
      DO 180 I = L1, N
  180 D(I) = D(I) - H
C
      F = F + H
C     ********** QL TRANSFORMATION **********
      P = D(M)
      C = 1.0D+00
      S = 0.0D+00
      MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
      DO 260 II = 1, MML
      I = M - II
      G = C * E(I)
      H = C * P
      IF (DABS(P) .LT. DABS(E(I))) GO TO 200
      C = E(I) / P
      R = DSQRT(C*C+1.0D+00)
      E(I+1) = S * P * R
      S = C / R
      C = 1.0D+00 / R
      GO TO 220
  200 C = P / E(I)
      R = DSQRT(C*C+1.0D+00)
      E(I+1) = S * E(I) * R
      S = 1.0D+00 / R
      C = C * S
  220 P = C * D(I) - S * G
      D(I+1) = H + S * (C * G + S * D(I))
C     ********** FORM VECTOR **********
      DO 240 K = 1, N
      H = Z(K,I+1)
      Z(K,I+1) = S * Z(K,I) + C * H
      Z(K,I) = C * Z(K,I) - S * H
  240 CONTINUE
C
  260 CONTINUE
C
      E(L) = S * P
      D(L) = C * P
      IF (DABS(E(L)) .GT. B) GO TO 160
  280 D(L) = D(L) + F
  300 CONTINUE
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********
      DO 360 II = 2, N
      I = II - 1
      K = I
      P = D(I)
C
      DO 320 J = II, N
      IF (D(J) .GE. P) GO TO 320
      K = J
      P = D(J)
  320 CONTINUE
C
      IF (K .EQ. I) GO TO 360
      D(K) = D(I)
      D(I) = P
C
      DO 340 J = 1, N
      P = Z(J,I)
      Z(J,I) = Z(J,K)
      Z(J,K) = P
  340 CONTINUE
C
  360 CONTINUE
C
      GO TO 400
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
  380 IERR = L
  400 RETURN
C     ********** LAST CARD OF TQL2 **********
      END
C
C   ------------------------------------------------------------------
C
      SUBROUTINE TRBK3B(NM,N,NV,A,M,Z)
      IMPLICIT NONE
C
      INTEGER I,J,L,M,N,IK,IZ,NM,NV
      DOUBLE PRECISION A(NV),Z(NM,M)
      DOUBLE PRECISION H,S
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3B.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS
C          USED IN THE REDUCTION BY  TRED3B IN ITS FIRST
C          N*(N+1)/2 POSITIONS,
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT-
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C  ------------------------------------------------------------------
C
      INTEGER  IDOT, IXPY
C
      IF (M .EQ. 0) GO TO 140
      IF (N .EQ. 1) GO TO 140
C
      DO 120 I = 2, N
      L = I - 1
      IZ = (I * L) / 2
      IK = IZ + I
      H = A(IK)
      IF (H .EQ. 0.0D+00) GO TO 120
C
      DO 100 J = 1, M
C REPLACE FUNCTION CALL to  DDOT AND DAXPY with in-line CODE
C     S = -DDOT(L,A(IZ+1),1,Z(1,J),1)
      S = 0.0
      do 5002 idot = 1,L
         S = S - A(IZ+IDOT)*Z(IDOT,J)
5002  continue
C
C     ********** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW **********
      S = (S / H) / H
C
C     CALL DAXPY(L,S,A(IZ+1),1,Z(1,J),1)
      IF (S.NE. 0.0) then
         do 5003 ixpy = 1,L
            Z(IXPY,j) = Z(IXPY,j) + S*A(IZ+IXPY)
5003     continue
      ENDIF
C
  100 CONTINUE
C
  120 CONTINUE
C
  140 RETURN
C     ********** LAST CARD OF TRBAK3 **********
      END
C
C   ------------------------------------------------------------------
C
      SUBROUTINE TRED3B(N,NV,A,D,E,E2)
      IMPLICIT NONE
C
      INTEGER I,J,K,L,N,II,IZ,JK,NV
      DOUBLE PRECISION A(NV),D(N),E(N),E2(N)
      DOUBLE PRECISION DT,F,G,H,HH,SCALE
C      DOUBLE PRECISION DSQRT,DABS,DSIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
C
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.
C
C     ON OUTPUT-
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
C          TRANSFORMATIONS USED IN THE REDUCTION,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C   ------------------------------------------------------------------
C
      INTEGER JM1
C
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********
      DO 300 II = 1, N
      I = N + 1 - II
      L = I - 1
      IZ = (I * L) / 2
      H = 0.0D+00
      SCALE = 0.0D+00
      IF (L .LT. 1) GO TO 120
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********
      DO 100 K = 1, L
      IZ = IZ + 1
      D(K) = A(IZ)
      SCALE = SCALE + DABS(D(K))
  100 CONTINUE
C
      IF (SCALE .NE. 0.0D+00) GO TO 140
  120 E(I) = 0.0D+00
      E2(I) = 0.0D+00
      GO TO 280
C
  140 DO 160 K = 1, L
      D(K) = D(K) / SCALE
      H = H + D(K) * D(K)
  160 CONTINUE
C
      E2(I) = SCALE * SCALE * H
      F = D(L)
      G = -DSIGN(DSQRT(H),F)
      E(I) = SCALE * G
      H = H - F * G
      D(L) = F - G
      A(IZ) = SCALE * D(L)
      IF (L .EQ. 1) GO TO 280
      F = 0.0D+00
C
      JK = 1
      DO 220 J = 1, L
      JM1 = J - 1
      DT = D(J)
      G = 0.0D+00
C     ********** FORM ELEMENT OF A*U **********
      IF (JM1 .EQ. 0) GO TO 200
      DO 180 K = 1, JM1
      E(K) = E(K) + DT * A(JK)
      G = G + D(K) * A(JK)
      JK = JK + 1
  180 CONTINUE
  200 E(J) = G + A(JK) * DT
      JK = JK + 1
C     ********** FORM ELEMENT OF P **********
  220 CONTINUE
      F = 0.0D+00
      DO 240 J = 1, L
      E(J) = E(J) / H
      F = F + E(J) * D(J)
  240 CONTINUE
C
      HH = F / (H + H)
      JK = 0
C     ********** FORM REDUCED A **********
      DO 260 J = 1, L
      F = D(J)
      G = E(J) - HH * F
      E(J) = G
C
      DO 260 K = 1, J
      JK = JK + 1
      A(JK) = A(JK) - F * E(K) - G * D(K)
  260 CONTINUE
C
  280 D(I) = A(IZ+1)
      A(IZ+1) = SCALE * DSQRT(H)
  300 CONTINUE
C
      RETURN
C     ********** LAST CARD OF TRED3B **********
      END
C====6====1=========2=========3=========4=========5=========6=========72
      SUBROUTINE MATMUL
     @ (M,N,L,DIM1,MAT1,TRANS1,DIM2,MAT2,TRANS2,DIM3,MAT3)
C
C takes matrix product MAT3 := MAT1*MAT2
C
C MAT1: M*L matrix, transpose if TRANSA
C MAT2: L*N matrix, transpose if TRANSB
C MAT3: M*N matrix
C DIM1, DIM2, DIM3 are the leading dimensions of the space
C   allocated for MAT1, MAT2, MAT3
C
C Author: Michael Nilges
C
      IMPLICIT NONE
C
C input/ output
      INTEGER M,N,L
      INTEGER DIM1,DIM2,DIM3
      DOUBLE PRECISION MAT1(DIM1,*),MAT2(DIM2,*),MAT3(DIM3,*)
      LOGICAL TRANS1,TRANS2
C
C local
      INTEGER I,J,K
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
C begin
C
      IF (TRANS1.AND.TRANS2) THEN
      DO I=1,N
      DO J=1,M
      MAT3(J,I)=ZERO
      DO K=1,L
      MAT3(J,I)=MAT3(J,I)+MAT1(J,K)*MAT2(K,I)
      END DO
      END DO
      END DO
C
      ELSE IF (TRANS1) THEN
      DO I=1,N
      DO J=1,M
      MAT3(J,I)=ZERO
      DO K=1,L
      MAT3(J,I)=MAT3(J,I)+MAT1(K,J)*MAT2(K,I)
      END DO
      END DO
      END DO
C
      ELSE IF (TRANS2) THEN
      DO I=1,N
      DO J=1,M
      MAT3(J,I)=ZERO
      DO K=1,L
      MAT3(J,I)=MAT3(J,I)+MAT1(J,K)*MAT2(I,K)
      END DO
      END DO
      END DO
C
      ELSE
      DO I=1,N
      DO J=1,M
      MAT3(J,I)=ZERO
      DO K=1,L
      MAT3(J,I)=MAT3(J,I)+MAT1(J,K)*MAT2(K,I)
      END DO
      END DO
      END DO
C
      END IF
C
      RETURN
      END
C====6====1=========2=========3=========4=========5=========6=========72
      SUBROUTINE DIAGSQ(DIM,N,SQM,EVECT,EVAL)
C
C finds eigenvalues and eigenvectors of square symmetric matrix SQM.
C
C SQM   (input): square symmetric matrix
C N     (input): is the dimension of the matrix,
C DIM   (input): the leading dimension of the space allocated
C                for SQM and EVECT.
C EVECT (output): eigenvectors are stored in array EVECT.
C EVAL  (output): eigenvalues are stored in EVAL.
C
C the square matrix is first copied to lower triangular matrix row by
C   row by CSQLTR. thus, the original square matrix is not destroyed.
C
C this machine independent version calls EISPACK routine GIVEIS.
C
C Author: Michael Nilges
C
      IMPLICIT NONE
C
C global
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
C input/ output
      INTEGER DIM, N
      DOUBLE PRECISION SQM(DIM,*), EVECT(DIM,*), EVAL(*)
C
C local
      INTEGER IERR, NPAIR
C
C heap pointers
      INTEGER HPTRMA, HPBMAT, HPIND
C
C begin
C
      IF (N .LT. 1) RETURN
      NPAIR=(N*(N+1))/2
      HPTRMA=ALLHP(IREAL8(NPAIR))
      HPBMAT=ALLHP(IREAL8(8*N))
      HPIND =ALLHP(INTEG4(N))
C
      CALL CSQLTR(N,DIM,SQM,HEAP(HPTRMA))
      CALL GIVEIS(N,N,DIM,HEAP(HPTRMA),HEAP(HPBMAT),HEAP(HPIND),
     @            EVAL,EVECT,IERR)
C
      CALL FREHP(HPIND, INTEG4(N))
      CALL FREHP(HPBMAT,IREAL8(8*N))
      CALL FREHP(HPTRMA,IREAL8(NPAIR))
C
      IF (IERR.NE.0) THEN
        WRITE(6,'(A,I3)')
     @' %GIVEIS-ERR: diagonalization failed, IERR= ', IERR
      END IF
C
      RETURN
      END
C====6====1=========2=========3=========4=========5=========6=========72
      SUBROUTINE CSQLTR(N,DIM,SQM,TRM)
C
C copies lower triangle of symmetric sqare matrix SQM with leading
C dimension DIM and of order N into linear array TRM
C row by row.
C
C N   (input): order of square matrix
C DIM (input): leading dimension of its space allocation
C SQM (input): square matrix
C TRM (output): trianglular matrix
C
C Author: Michael Nilges
C
      IMPLICIT NONE
C
C input/ output
      INTEGER N,DIM
      DOUBLE PRECISION SQM(DIM,*),TRM(*)
C
C local
      INTEGER I,J,IIM
C
C begin
C
      DO I=1,N
        IIM=I*(I-1)/2
        DO J=1,I
          TRM(IIM+J)=SQM(I,J)
        END DO
      END DO
C
      RETURN
      END
