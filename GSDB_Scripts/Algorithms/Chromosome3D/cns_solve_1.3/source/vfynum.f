C
C============================================================================
C
      SUBROUTINE VFYNUM
C
      IMPLICIT NONE
C
      INTEGER NR, N2DIM
      INTEGER ANSWER
C
      ANSWER=393
C
      CALL MAKTHE('SCAN',NR,0,'LATT',.FALSE.,
     &     0.0D0,360.0D0,0.0D0,90.0D0,0.0D0,720.0D0,
     &     30.0D0,0,0,0,0,N2DIM,0)
C
      IF (NR.NE.ANSWER) THEN
         WRITE(6, '(1X, A)') 
     &   'Fatal Error: test of numerical function failed'
         CALL FAQDIE
      END IF
C
      RETURN
      END
C
C============================================================================
C
