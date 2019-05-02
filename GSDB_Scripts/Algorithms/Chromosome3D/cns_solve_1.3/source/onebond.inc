C===============
C
C OneBond.inc
C
C by John Kuszewski May 1995
C===============
      INTEGER MAXONECLASSES
      PARAMETER (MAXONECLASSES = 10)
C
C arrays that hold coupling info
C coupassndx tells ending index of the coup_x arrays (below)
C for each class.
C coupforces holds k1 and coupswitches holds k2 for each class
C coupclasstypes holds the type of each class
C
C these are small enough not to bother with the heap.
C
      INTEGER ONEASSNDX (MAXONECLASSES)
      INTEGER ONESWITCHES(MAXONECLASSES)
      DOUBLE PRECISION ONEFORCES (MAXONECLASSES)
      CHARACTER*8 ONECLASSNAMES (MAXONECLASSES)
C
C holders for the Karplus curve coefficients
C for each class
C
      DOUBLE PRECISION ONEAS(MAXONECLASSES),
     &     ONEBS(MAXONECLASSES),
     &     ONECS(MAXONECLASSES), ONEDS(MAXONECLASSES),
     &     ONEPHASE1S(MAXONECLASSES),
     &     ONEPHASE2S(MAXONECLASSES), ONEPHASE3S(MAXONECLASSES)
C
C maxcoups = number of slots set aside for coupling constant
C            assignments
C ncoups = total number of coupling constants entered
C
      INTEGER MAXONES, NONES, NCLASSES, CURCLASS
C
C pointers to arrays to hold atom numbers, J-obs, and
C error for each coupling assignment
C
      INTEGER ONEIPTR, ONEJPTR, ONEKPTR, ONELPTR, ONEMPTR,
     &     ONENPTR, ONEPPTR, ONEQPTR, ONEJOBSPTR, ONEJERRPTR,
     &     CALCONEPTR
C
C input modes
C
      INTEGER MODE, NEW, UPDATE
      PARAMETER (NEW = 1)
      PARAMETER (UPDATE = 2)
C
C potential types
C
       INTEGER POTENTIAL, SQUARE, HARMONIC
       PARAMETER (SQUARE = 1)
       PARAMETER (HARMONIC = 2)
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
      DOUBLE PRECISION EPS
      PARAMETER (EPS=0.1D0)
C
C common blocks
C
      COMMON /CONEB/ ONECLASSNAMES
      COMMON /IONEB1/ ONEASSNDX, POTENTIAL,
     &     MAXONES, NONES, ONESWITCHES, CURCLASS, NCLASSES
      COMMON /IONEB2/ ONEIPTR, ONEJPTR, ONEKPTR, ONELPTR,
     &     ONEMPTR, ONENPTR, ONEPPTR, ONEQPTR,
     &     ONEJOBSPTR, ONEJERRPTR, CALCONEPTR, MODE
      COMMON /RONEB/ ONEFORCES, ONEAS, ONEBS, ONECS, ONEDS,
     &     ONEPHASE1S, ONEPHASE2S, ONEPHASE3S
      SAVE /CONEB/
      SAVE /IONEB1/
      SAVE /IONEB2/
      SAVE /RONEB/
