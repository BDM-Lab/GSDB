C===============
C
C Cshifts.inc
C
C
C by John Kuszewski Nov 1993
C===============
      INTEGER MAXSHIFTCLASSES
      PARAMETER (MAXSHIFTCLASSES = 10)
C
C arrays that hold coupling info
C shiftassndx tells ending index of the shift_x arrays (below)
C for each class.
C shiftforces holds k1 for each class
C
C these are small enough not to bother with the heap.
C
      INTEGER SHIFTASSNDX (MAXSHIFTCLASSES)
      DOUBLE PRECISION SHIFTFORCES (MAXSHIFTCLASSES)
      CHARACTER*8 SHIFTCLASSNAMES (MAXSHIFTCLASSES)
C
C maxshifts = number of slots set aside for chemical shift
C             assignments
C nshifts = total number of coupling constants entered
C
      INTEGER MAXSHIFTS, NSHIFTS, NSHIFTCLASSES, CURSHIFTCLASS
      INTEGER CARBNATOM
C
C pointers to arrays to hold atom numbers, observed Ca and Cb shifts,
C calculated Ca and Cb shifts, and random coil Ca and Cb shifts
C
      INTEGER SHIFTIPTR, SHIFTJPTR, SHIFTKPTR, SHIFTLPTR,
     &        SHIFTMPTR, CAOBSPTR, CBOBSPTR, CACALCPTR, CBCALCPTR,
     &        EXCAPTR, EXCBPTR, EXCAERRPTR, EXCBERRPTR,
     &        RCAPTR, RCBPTR
C
C whether or not we've allocated the rcoil heap
C
      LOGICAL CARBRCALLOCATED
C
C number of entries along the phi and psi dimensions
C of the expectation arrays
C
      INTEGER PHISTEPS, PSISTEPS
C
C input modes
C
      INTEGER SHIFTMODE, NEW, UPDATE
      PARAMETER (NEW = 1)
      PARAMETER (UPDATE = 2)
C
C potential types
C
       INTEGER SHIFTPOTENTIAL, SQUARE, HARMONIC
       PARAMETER (SQUARE = 1)
       PARAMETER (HARMONIC = 2)
C
C use this to identify no-expectation slots
C This should be the same as in rama.inc
C
       DOUBLE PRECISION NOEXPECTATION
       PARAMETER (NOEXPECTATION = 1.0E20)
       DOUBLE PRECISION NOEXPTEST
       PARAMETER (NOEXPTEST = 1.0E17)
C
C parameters as set up in ETOR
C
C===>parameter MCONST truncates the reciprocals and makes sure
C===>we don't devide by zero
      DOUBLE PRECISION MCONST
      PARAMETER (MCONST=0.0001D0)
C===>parameter EPS determines the absolute value of SIN(PHI) at
C===>which we want to switch to the alternative computation for the
C===>derivative
      DOUBLE PRECISION CARBEPS
      PARAMETER (CARBEPS=0.1D0)
C
C
C common blocks
C
      COMMON /CCSHIF/ SHIFTCLASSNAMES
      COMMON /ICSHIF1/ SHIFTASSNDX, SHIFTPOTENTIAL,
     &       MAXSHIFTS, NSHIFTS, CURSHIFTCLASS, NSHIFTCLASSES,
     &       PHISTEPS, PSISTEPS, SHIFTMODE, CARBNATOM
      COMMON /ICSHIF2/ SHIFTIPTR, SHIFTJPTR, SHIFTKPTR, SHIFTLPTR,
     &                SHIFTMPTR, CAOBSPTR, CBOBSPTR, CACALCPTR,
     &                CBCALCPTR, RCAPTR, RCBPTR,
     &                EXCAPTR, EXCBPTR, EXCAERRPTR, EXCBERRPTR
      COMMON /LCSHIF/ CARBRCALLOCATED
      COMMON /DCSHIF/ SHIFTFORCES
      SAVE /CCSHIF/
      SAVE /ICSHIF1/
      SAVE /ICSHIF2/
      SAVE /DCSHIF/
      SAVE /LCSHIF/
