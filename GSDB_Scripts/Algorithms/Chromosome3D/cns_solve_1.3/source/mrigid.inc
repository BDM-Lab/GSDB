C data structure MRIGID
C =====================
C
C Common Block for some parameters of the rigid body minimizer
C The main purpose of the common block is to pass certain
C parameters to the FORTRAN function used by the minimization
C routine ZXCGR.
C
C max. no. groups
      INTEGER MAXNGR
C
C number of groups
      INTEGER NGROUP
C
C group pointers
      INTEGER IPOINT
C
C heap pointer for atom list
      INTEGER HPLIST
C
C print flag
      LOGICAL RPRINT
C
C translation flag
      LOGICAL QTRANS
C
C print period
      INTEGER NPRINT
C
C number of calls
      INTEGER NCALLS
C==========================================================================
C
C     integer block
C
      COMMON/MRIGI/ NGROUP, IPOINT, HPLIST, MAXNGR, NPRINT, NCALLS
C
C     logical block
C
      COMMON /LMRIG/ RPRINT, QTRANS
C
      SAVE /MRIGI/
      SAVE /LMRIG/
