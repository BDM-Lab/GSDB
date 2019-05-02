C
C trtraj.inc
C ==========
C
C variables for trajectory writing routines
C
      LOGICAL TRFIRS, TRQFOR, TRINIT
      INTEGER TRNUM, TRIND, LTRIN
      INTEGER TRUNIT
      CHARACTER*10 TRFORM
      DOUBLE PRECISION TRSCAL, TROFFS
C
C variables for trajectory reading routines
C
C
C pointers for READC
      INTEGER TRTEMP, TRFREEAT, TRIX, TRIY, TRIZ, TRNAT
C
C
C COMMON BLOCKS
C =============
C
C double precision block
C
      COMMON /TRTRA/  TRSCAL, TROFFS
C integer block
C
      COMMON /ITRTRA/ TRNUM, TRIND, LTRIN, TRUNIT,
     &                TRTEMP, TRFREEAT, TRIX, TRIY, TRIZ,
     &                TRNAT
C character common block
C
      COMMON /STRTRA/ TRFORM
C
C logical common block
      COMMON /LTRTRA/ TRFIRS, TRQFOR, TRINIT
C
      SAVE /TRTRA/
      SAVE /ITRTRA/
      SAVE /STRTRA/
      SAVE /LTRTRA/
