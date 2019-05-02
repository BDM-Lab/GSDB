C
C ==========================================================================
C
      SUBROUTINE DPSTOR(X)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'machvar.inc'
      DOUBLE PRECISION  X
C
C FORCES THE ARGUMENT VALUE X TO BE STORED IN MEMORY LOCATION V.
C Source code taken from fftpack41/tfftpk.f
C
      DPTMPV=X
      RETURN
      END
C
C ==========================================================================
C
      DOUBLE PRECISION FUNCTION DPTRUNC(X)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'machvar.inc'
      DOUBLE PRECISION  X
C
C DPTRUNC IS A PORTABLE FORTRAN FUNCTION WHICH TRUNCATES A VALUE TO THE
C MACHINE DOUBLE PRECISION WORD SIZE, REGARDLESS OF WHETHER LONGER
C PRECISION INTERNAL REGISTERS ARE USED FOR FLOATING POINT ARITHMETIC IN
C COMPUTING THE VALUE INPUT TO DPTRUNC.  THE METHOD USED IS TO FORCE A
C STORE INTO MEMORY BY USING A COMMON BLOCK IN ANOTHER SUBROUTINE.
C Source code taken from fftpack41/tfftpk.f
C
      CALL DPSTOR(X)
      DPTRUNC=DPTMPV
      RETURN
      END
C
C ==========================================================================
C
      SUBROUTINE SETFPEPS
C
      IMPLICIT NONE
C I/O
      INCLUDE 'numbers.inc'
      INCLUDE 'machvar.inc'
C
C Determine smallest FPEPS such that both 1+FPEPS and 1-FPEPS
C are different from 1.
C
C local
      INTEGER           I
      DOUBLE PRECISION  ONEP, ONEM
      DOUBLE COMPLEX    DBCOMP
C
C external
      DOUBLE PRECISION  DPTRUNC
      EXTERNAL          DPTRUNC
C
C begin
      DBCOMP = DCMPLX(ZERO,ZERO)
      I = 1
      DO WHILE (I .LT. MXFPEPS2)
        FPEPS = TWO**(-I)
        ONEP = DPTRUNC(ONE) + DPTRUNC(FPEPS)
        ONEM = DPTRUNC(ONE) - DPTRUNC(FPEPS)
        IF (ONE .EQ. ONEP .OR. ONE .EQ. ONEM) THEN
          I = I - 1
          FPEPS = TWO**(-I)
          CALL DECLAR('FP_EPSILON', 'DP', ' ', DBCOMP, FPEPS)
          IF ( I .GT. 128 ) THEN
             WRITE (6,'(A,E10.3)')
     &       ' %SETFPEPS Machine epsilon determined to be ',FPEPS
             CALL WRNDIE(-5,'SETFPEPS',
     &       'Machine epsilon value is too small')
             CALL DIE
          END IF
          RETURN
        END IF
        I = I + 1
      END DO
C
      WRITE(6,'(A)')
     &      ' %SETFPEPS increase value of MXFPEPS2 and recompile'
      CALL WRNDIE(-5, 'SETFPEPS',
     &            'Could not determine machine epsilon')
      CALL DIE
C
      RETURN
      END
C
C ==========================================================================
C
