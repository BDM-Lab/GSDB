C===============
C
C vectangl.fcm
C
C by Jens Meiler / Michael Nilges Jan 1999
C===============
      INTEGER MAXVEANCLASSES
      PARAMETER (MAXVEANCLASSES = 10)
C
C arrays that hold vectorangle info
C vectassndx tells ending index of the vect_x arrays (below)
C for each class.
C vectforces holds k1 and k2 for each class
C vectclasstypes holds the type of each class
C
C these are small enough not to bother with the heap.
C
      INTEGER VEANASSNDX    (MAXVEANCLASSES)
      DOUBLE PRECISION VEANFORCES(2,MAXVEANCLASSES)
      CHARACTER*8 VEANCLASSNAMES (MAXVEANCLASSES)
      LOGICAL PRINTCLASS(MAXVEANCLASSES)
C
C maxvects = number of slots set aside for vectorangle constant
C            assignments
C nvects = total number of vectorangle constants entered
C
      INTEGER MAXVEANS, NVEANS, NCLASSES, CURCLASS
C
C pointers to arrays to hold atom numbers, J-obs, and
C error for each vectorangle assignment
C
      INTEGER VEANIPTR, VEANJPTR, VEANKPTR, VEANLPTR,
     &     VEANVOBSPTR, VEANVERRPTR, CALCVEANPTR
C
C  cross-validation
C  by Jens Meiler / Michael Nilges Jan 1999
C
C cross-validation array
      INTEGER VEANCV
C
C cross-validation test number
      INTEGER VICV
C
C input modes
C
      INTEGER MODE, NEW, UPDATE
      PARAMETER (NEW = 1)
      PARAMETER (UPDATE = 2)
C
C===>parameter MCONST truncates the reciprocals and makes sure
C===>we don't devide by zero
      DOUBLE PRECISION MCONST
      PARAMETER (MCONST=0.0001D0)
C
C common blocks
C
      COMMON /CVEAN/ VEANCLASSNAMES
      COMMON /IVEAN/ VEANASSNDX,
     &     MAXVEANS, NVEANS, CURCLASS, NCLASSES,
     &     VEANIPTR, VEANJPTR, VEANKPTR, VEANLPTR,
     &     VEANVOBSPTR, VEANVERRPTR, CALCVEANPTR, MODE,
     &     VICV, VEANCV
      COMMON /RVEAN/ VEANFORCES
      COMMON /LVEAN/ PRINTCLASS
      SAVE /CVEAN/
      SAVE /IVEAN/
      SAVE /RVEAN/
      SAVE /LVEAN/

