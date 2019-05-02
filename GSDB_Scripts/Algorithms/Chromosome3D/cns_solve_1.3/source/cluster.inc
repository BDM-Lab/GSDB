C rsearch.inc
C ==========
C Rotation search
C
C Paul Adams 3/14/98
C
C allocated space for peaklist
      INTEGER PEAKALLOC
C
C number of peaks in list
      INTEGER NPEAK
C
C pointers
      INTEGER HPT1, HPT2, HPT3, HPRF
C
C ====================================================================
C
C     integer block
C
      COMMON /IRSH/ PEAKALLOC, NPEAK
C
C     pointer block
C
      COMMON /PRSH/ HPT1, HPT2, HPT3, HPRF
C
      SAVE /IRSH/
      SAVE /PRSH/
C
