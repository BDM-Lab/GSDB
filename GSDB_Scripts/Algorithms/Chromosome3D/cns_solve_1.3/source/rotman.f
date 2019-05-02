C   ROTMAN.S
C   ========
C
      SUBROUTINE ROTMAN
C
C Routine computes relationships between angles
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'xcrystal.inc'
C local
      DOUBLE PRECISION ROT(3,3), ROTB(3,3), ROTM(3,3), TEMP, DET
      DOUBLE PRECISION AOP(3,3), DIST, ROTMM(3,3), ROTA(3,3)
      DOUBLE PRECISION DISTOLD
      DOUBLE COMPLEX DBCOMP
      INTEGER I, J, K, U, V, X, Y, ISYM1, ISYM2, IISYM1, IISYM2
C parameter
      DOUBLE PRECISION ZERO, ONE, BIG
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, BIG=9999.9D0)
C begin
CSGI-specific optimization
C*$* OPTIMIZE(4)
C
C defaults (set both primary and secondary matrix to identity)
      DO I=1,3
      DO J=1,3
      ROT(I,J)=ZERO
      ROTB(I,J)=ZERO
      END DO
      ROT(I,I)=ONE
      ROTB(I,I)=ONE
      END DO
C
      CALL PUSEND('ROTMAN>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('ROTMAN>')
      CALL MISCOM('ROTMAN>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-rotman')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MATR'.OR.WD(1:4).EQ.'EULE'.OR.
     &         WD(1:4).EQ.'LATT'.OR.WD(1:4).EQ.'SPHE'.OR.
     &         WD(1:4).EQ.'AXIS') THEN
C
C compute primary matrix
      CALL MATPAR(WD(1:4),ROT)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'COPY') THEN
C
C copy primary marix into secondary matrix
      DO I=1,3
      DO J=1,3
      ROTB(I,J)=ROT(I,J)
      END DO
      END DO
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SWAP') THEN
C
C swap primary and secondary marix
      DO I=1,3
      DO J=1,3
      TEMP=ROTB(I,J)
      ROTB(I,J)=ROT(I,J)
      ROT(I,J)=TEMP
      END DO
      END DO
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
C
C loop over crystallographic symmetry images
      DO ISYM1=1,XRNSYM
C
C compute symmetry operator in real space
      DO U=1,3
      DO V=1,3
      AOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      AOP(U,V)=AOP(U,V)+XRINTR(U,X)*XRSYMM(ISYM1,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
      DO I=1,3
      DO K=1,3
      ROTM(I,K)=ZERO
      DO J=1,3
      ROTM(I,K)=ROTM(I,K)+AOP(I,J)*ROT(J,K)
      END DO
      END DO
      END DO
      WRITE(6,'(2A)')
     & ' =====================================================',
     & '========================='
      WRITE(6,'(/A,I2,A)') ' [symmetry operator No. ',ISYM1,
     &  '] * [primary matrix]'
      CALL MATPRI(ROTM)
      WRITE(6,'(2A)')
     & ' =====================================================',
     & '========================='
      END DO
C
C now loop over crystallographic symmetry images
      DO ISYM1=1,XRNSYM
C
C compute symmetry operator in real space
      DO U=1,3
      DO V=1,3
      AOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      AOP(U,V)=AOP(U,V)+XRINTR(U,X)*XRSYMM(ISYM1,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
      DO I=1,3
      DO K=1,3
      ROTM(I,K)=ZERO
      DO J=1,3
      ROTM(I,K)=ROTM(I,K)+AOP(I,J)*ROTB(J,K)
      END DO
      END DO
      END DO
      WRITE(6,'(2A)')
     & ' =====================================================',
     & '========================='
      WRITE(6,'(/A,I2,A)') ' [symmetry operator No. ',ISYM1,
     &  '] * [secondary matrix]'
      CALL MATPRI(ROTM)
      WRITE(6,'(2A)')
     & ' =====================================================',
     & '========================='
      END DO
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PROD') THEN
C
C multiply primary and secondary matrix. Stores result in primary
C matrix.  Secondary matrix is untouched.
      DO I=1,3
      DO K=1,3
      ROTM(I,K)=ZERO
      DO J=1,3
      ROTM(I,K)=ROTM(I,K)+ROT(I,J)*ROTB(J,K)
      END DO
      END DO
      END DO
      DO I=1,3
      DO K=1,3
      ROT(I,K)=ROTM(I,K)
      END DO
      END DO
C=====================================================================
      ELSE IF (WD(1:4).EQ.'INVE') THEN
C
C compute the inverse of the primary matrix and stores result in the
C primary matrix
C Calculate the determinant of the matrix ROT.
      DET=ROT(1,1)*(ROT(2,2)*ROT(3,3)-ROT(2,3)*ROT(3,2))
     &   -ROT(1,2)*(ROT(2,1)*ROT(3,3)-ROT(2,3)*ROT(3,1))
     &   +ROT(1,3)*(ROT(2,1)*ROT(3,2)-ROT(2,2)*ROT(3,1))
C
C calculate ADJOINT matrix of R1
C transpose of the signed cofactor matrix = Adj[R1]
      ROTM(1,1)=ROT(2,2)*ROT(3,3)-ROT(2,3)*ROT(3,2)
      ROTM(2,1)=ROT(3,1)*ROT(2,3)-ROT(2,1)*ROT(3,3)
      ROTM(3,1)=ROT(2,1)*ROT(3,2)-ROT(3,1)*ROT(2,2)
      ROTM(1,2)=ROT(3,2)*ROT(1,3)-ROT(1,2)*ROT(3,3)
      ROTM(2,2)=ROT(1,1)*ROT(3,3)-ROT(3,1)*ROT(1,3)
      ROTM(3,2)=ROT(3,1)*ROT(1,2)-ROT(1,1)*ROT(3,2)
      ROTM(1,3)=ROT(1,2)*ROT(2,3)-ROT(1,3)*ROT(2,2)
      ROTM(2,3)=ROT(2,1)*ROT(1,3)-ROT(1,1)*ROT(2,3)
      ROTM(3,3)=ROT(1,1)*ROT(2,2)-ROT(2,1)*ROT(1,2)
C
C now calculate the inverse matrix only for non zero det
      IF(ABS(DET).GE.RSMALL) THEN
      DO I=1,3
      DO J=1,3
      ROT(I,J)=ROTM(I,J)/DET
      END DO
      END DO
      ELSE
      WRITE(6,'(A)') ' %ROTMAN-ERR: primary matrix singular.'
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DIST') THEN
C
C compute distance between primary and secondary matrix
      DISTOLD=BIG
C
C 1st loop through all crystallographic symmetry operators
      DO ISYM1=1,XRNSYM
C
C apply symmetry operator to primary matrix
      DO U=1,3
      DO V=1,3
      AOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      AOP(U,V)=AOP(U,V)+XRINTR(U,X)*XRSYMM(ISYM1,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
      DO I=1,3
      DO K=1,3
      ROTA(I,K)=ZERO
      DO J=1,3
      ROTA(I,K)=ROTA(I,K)+AOP(I,J)*ROT(J,K)
      END DO
      END DO
      END DO
C
C 2nd loop through all crystallographic symmetry operators
      DO ISYM2=1,XRNSYM
C
C apply symmetry operator to secondary matrix
      DO U=1,3
      DO V=1,3
      AOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      AOP(U,V)=AOP(U,V)+XRINTR(U,X)*XRSYMM(ISYM2,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
      DO I=1,3
      DO K=1,3
      ROTM(I,K)=ZERO
      DO J=1,3
      ROTM(I,K)=ROTM(I,K)+AOP(I,J)*ROTB(J,K)
      END DO
      END DO
      END DO
      DIST=SQRT(
     &      (ROTA(1,1)-ROTM(1,1))**2
     &     +(ROTA(1,2)-ROTM(1,2))**2
     &     +(ROTA(1,3)-ROTM(1,3))**2
     &     +(ROTA(2,1)-ROTM(2,1))**2
     &     +(ROTA(2,2)-ROTM(2,2))**2
     &     +(ROTA(2,3)-ROTM(2,3))**2
     &     +(ROTA(3,1)-ROTM(3,1))**2
     &     +(ROTA(3,2)-ROTM(3,2))**2
     &     +(ROTA(3,3)-ROTM(3,3))**2 )
      IF (DIST.LT.DISTOLD-RSMALL) THEN
      DISTOLD=DIST
      IISYM1=ISYM1
      IISYM2=ISYM2
      END IF
      END DO
      END DO
C
C now print the info about the minimal distance
C
C apply symmetry operator to primary matrix
      DO U=1,3
      DO V=1,3
      AOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      AOP(U,V)=AOP(U,V)+XRINTR(U,X)*XRSYMM(IISYM1,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
      DO I=1,3
      DO K=1,3
      ROTA(I,K)=ZERO
      DO J=1,3
      ROTA(I,K)=ROTA(I,K)+AOP(I,J)*ROT(J,K)
      END DO
      END DO
      END DO
C
C apply symmetry operator to secondary matrix
      DO U=1,3
      DO V=1,3
      AOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      AOP(U,V)=AOP(U,V)+XRINTR(U,X)*XRSYMM(IISYM2,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
      DO I=1,3
      DO K=1,3
      ROTM(I,K)=ZERO
      DO J=1,3
      ROTM(I,K)=ROTM(I,K)+AOP(I,J)*ROTB(J,K)
      END DO
      END DO
      END DO
      WRITE(6,'(2A)')
     & ' ===following symmetry operators yield the minimum distance',
     & '===================='
      WRITE(6,'(/A,I2,A)') ' Matrix A = [symmetry operator No. ',
     & IISYM1,
     &  '] * [primary matrix]'
      WRITE(6,'(A,I2,A)')
     &  ' Matrix B = [symmetry operator No. ',IISYM2,
     &  '] * [secondary matrix]'
      WRITE(6,'(/A,F10.4)')
     & ' Distance between matrices A and B=',DISTOLD
      CALL DECLAR( 'DISTANCE', 'DP', ' ', DBCOMP, DISTOLD )
C
C now compute S[1]*primary * [ S[2]*secondary ]-1
C Calculate the determinant of the matrix ROTM.
      DET=ROTM(1,1)*(ROTM(2,2)*ROTM(3,3)-ROTM(2,3)*ROTM(3,2))
     &   -ROTM(1,2)*(ROTM(2,1)*ROTM(3,3)-ROTM(2,3)*ROTM(3,1))
     &   +ROTM(1,3)*(ROTM(2,1)*ROTM(3,2)-ROTM(2,2)*ROTM(3,1))
C
C calculate ADJOINT matrix of R1
C transpose of the signed cofactor matrix = Adj[R1]
      ROTMM(1,1)=ROTM(2,2)*ROTM(3,3)-ROTM(2,3)*ROTM(3,2)
      ROTMM(2,1)=ROTM(3,1)*ROTM(2,3)-ROTM(2,1)*ROTM(3,3)
      ROTMM(3,1)=ROTM(2,1)*ROTM(3,2)-ROTM(3,1)*ROTM(2,2)
      ROTMM(1,2)=ROTM(3,2)*ROTM(1,3)-ROTM(1,2)*ROTM(3,3)
      ROTMM(2,2)=ROTM(1,1)*ROTM(3,3)-ROTM(3,1)*ROTM(1,3)
      ROTMM(3,2)=ROTM(3,1)*ROTM(1,2)-ROTM(1,1)*ROTM(3,2)
      ROTMM(1,3)=ROTM(1,2)*ROTM(2,3)-ROTM(1,3)*ROTM(2,2)
      ROTMM(2,3)=ROTM(2,1)*ROTM(1,3)-ROTM(1,1)*ROTM(2,3)
      ROTMM(3,3)=ROTM(1,1)*ROTM(2,2)-ROTM(2,1)*ROTM(1,2)
C
C now calculate the inverse matrix only for non zero det
      IF(ABS(DET).GE.RSMALL) THEN
      DO I=1,3
      DO J=1,3
      ROTMM(I,J)=ROTMM(I,J)/DET
      END DO
      END DO
      ELSE
      WRITE(6,'(A)') ' %ROTMAN-ERR: secondary matrix singular.'
      END IF
      DO I=1,3
      DO K=1,3
      ROTM(I,K)=ZERO
      DO J=1,3
      ROTM(I,K)=ROTM(I,K)+ROTA(I,J)*ROTMM(J,K)
      END DO
      END DO
      END DO
C
      WRITE(6,'(A,I2,A,I2,A)') ' Matrix A * B-1 '
      CALL MATPRI(ROTM)
      CALL MATDCL(ROTM)
      WRITE(6,'(2A)')
     & ' =====================================================',
     & '========================='
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ALLP') THEN
C
C compute distance between primary and secondary matrix
C
C 1st loop through all crystallographic symmetry operators
      DO ISYM1=1,XRNSYM
C
C apply symmetry operator to primary matrix
      DO U=1,3
      DO V=1,3
      AOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      AOP(U,V)=AOP(U,V)+XRINTR(U,X)*XRSYMM(ISYM1,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
      DO I=1,3
      DO K=1,3
      ROTA(I,K)=ZERO
      DO J=1,3
      ROTA(I,K)=ROTA(I,K)+AOP(I,J)*ROT(J,K)
      END DO
      END DO
      END DO
C
C 2nd loop through all crystallographic symmetry operators
      DO ISYM2=1,XRNSYM
C
C apply symmetry operator to secondary matrix
      DO U=1,3
      DO V=1,3
      AOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      AOP(U,V)=AOP(U,V)+XRINTR(U,X)*XRSYMM(ISYM2,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
      DO I=1,3
      DO K=1,3
      ROTM(I,K)=ZERO
      DO J=1,3
      ROTM(I,K)=ROTM(I,K)+AOP(I,J)*ROTB(J,K)
      END DO
      END DO
      END DO
      DIST=SQRT(
     &      (ROTA(1,1)-ROTM(1,1))**2
     &     +(ROTA(1,2)-ROTM(1,2))**2
     &     +(ROTA(1,3)-ROTM(1,3))**2
     &     +(ROTA(2,1)-ROTM(2,1))**2
     &     +(ROTA(2,2)-ROTM(2,2))**2
     &     +(ROTA(2,3)-ROTM(2,3))**2
     &     +(ROTA(3,1)-ROTM(3,1))**2
     &     +(ROTA(3,2)-ROTM(3,2))**2
     &     +(ROTA(3,3)-ROTM(3,3))**2 )
C
C now print the info
C
C apply symmetry operator to primary matrix
      DO U=1,3
      DO V=1,3
      AOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      AOP(U,V)=AOP(U,V)+XRINTR(U,X)*XRSYMM(ISYM1,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
      DO I=1,3
      DO K=1,3
      ROTA(I,K)=ZERO
      DO J=1,3
      ROTA(I,K)=ROTA(I,K)+AOP(I,J)*ROT(J,K)
      END DO
      END DO
      END DO
C
C apply symmetry operator to secondary matrix
      DO U=1,3
      DO V=1,3
      AOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      AOP(U,V)=AOP(U,V)+XRINTR(U,X)*XRSYMM(ISYM2,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
      DO I=1,3
      DO K=1,3
      ROTM(I,K)=ZERO
      DO J=1,3
      ROTM(I,K)=ROTM(I,K)+AOP(I,J)*ROTB(J,K)
      END DO
      END DO
      END DO
      WRITE(6,'(2A)')
     & ' =====================================================',
     & '========================='
      WRITE(6,'(/A,I2,A)') ' Matrix A = [symmetry operator No. ',
     & ISYM1,
     &  '] * [primary matrix]'
      WRITE(6,'(A,I2,A)')
     &  ' Matrix B = [symmetry operator No. ',ISYM2,
     &  '] * [secondary matrix]'
      WRITE(6,'(/A,F10.4)')
     & ' Distance between matrices A and B=',DIST
C
C now compute S[1]*primary * [ S[2]*secondary ]-1
C Calculate the determinant of the matrix ROTM.
      DET=ROTM(1,1)*(ROTM(2,2)*ROTM(3,3)-ROTM(2,3)*ROTM(3,2))
     &   -ROTM(1,2)*(ROTM(2,1)*ROTM(3,3)-ROTM(2,3)*ROTM(3,1))
     &   +ROTM(1,3)*(ROTM(2,1)*ROTM(3,2)-ROTM(2,2)*ROTM(3,1))
C
C calculate ADJOINT matrix of R1
C transpose of the signed cofactor matrix = Adj[R1]
      ROTMM(1,1)=ROTM(2,2)*ROTM(3,3)-ROTM(2,3)*ROTM(3,2)
      ROTMM(2,1)=ROTM(3,1)*ROTM(2,3)-ROTM(2,1)*ROTM(3,3)
      ROTMM(3,1)=ROTM(2,1)*ROTM(3,2)-ROTM(3,1)*ROTM(2,2)
      ROTMM(1,2)=ROTM(3,2)*ROTM(1,3)-ROTM(1,2)*ROTM(3,3)
      ROTMM(2,2)=ROTM(1,1)*ROTM(3,3)-ROTM(3,1)*ROTM(1,3)
      ROTMM(3,2)=ROTM(3,1)*ROTM(1,2)-ROTM(1,1)*ROTM(3,2)
      ROTMM(1,3)=ROTM(1,2)*ROTM(2,3)-ROTM(1,3)*ROTM(2,2)
      ROTMM(2,3)=ROTM(2,1)*ROTM(1,3)-ROTM(1,1)*ROTM(2,3)
      ROTMM(3,3)=ROTM(1,1)*ROTM(2,2)-ROTM(2,1)*ROTM(1,2)
C
C now calculate the inverse matrix only for non zero det
      IF(ABS(DET).GE.RSMALL) THEN
      DO I=1,3
      DO J=1,3
      ROTMM(I,J)=ROTMM(I,J)/DET
      END DO
      END DO
      ELSE
      WRITE(6,'(A)') ' %ROTMAN-ERR: secondary matrix singular.'
      END IF
      DO I=1,3
      DO K=1,3
      ROTM(I,K)=ZERO
      DO J=1,3
      ROTM(I,K)=ROTM(I,K)+ROTA(I,J)*ROTMM(J,K)
      END DO
      END DO
      END DO
C
      WRITE(6,'(A,I2,A,I2,A)') ' Matrix A * B-1 '
      CALL MATPRI(ROTM)
      WRITE(6,'(2A)')
     & ' =====================================================',
     & '========================='
C
      END DO
      END DO
C=====================================================================
      ELSE
      CALL CHKEND('ROTMAN>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
CSGI-specific optimization
C*$* OPTIMIZE(5)
      RETURN
      END
