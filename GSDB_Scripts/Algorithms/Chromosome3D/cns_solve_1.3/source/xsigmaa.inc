C
C xsigmaa.inc
C
C used to store pointers for certain data elements
C required for the function passed to the powell minimizer
C - data must be passed through a common block
C
C Paul Adams and Axel Brunger
C
C parameters
      INTEGER XSNSTP
      DOUBLE PRECISION XSTOL, XSSTSZ
      DOUBLE PRECISION XSMAX, XSMIN
C number of minimization steps
      PARAMETER (XSNSTP=50)
C tolerance for minimizer and initial step size
      PARAMETER (XSTOL=1.0D-3, XSSTSZ=1.0D-1)
C maximum and minimum allowed values for sigmaa
      PARAMETER (XSMAX=0.999D0, XSMIN=0.001D0)
C number of active reflections and current cycle
      INTEGER XSREF, XSCYCLE
C weight for log likelihood and linear line restraint
C XSUSRW is a user supplied weight that is applied to the
C the initial weight determined from comparing gradients to give
C XSWGHT (the actual weight that will be applied)
      DOUBLE PRECISION XSWGHT, XSUSRW
C phistep for integration in SIGAHLCV code
      DOUBLE PRECISION PHISTEP
C lowest target value and corresponding cycle number
      DOUBLE PRECISION XSTARG
      INTEGER XSTCYC
C
C pointers
      INTEGER XSEOBS, XSECAL, XSINDX, XSSHLL, XSASTL2, XSWSUM
      INTEGER PSIGMAA, HPPTYPE, XGOOD, XPCALC, XPA, XPB, XPC, XPD
C
C integer common block
C
      COMMON /XSIGI/ XSREF, XSCYCLE, XSEOBS, XSECAL, XSINDX, XSSHLL,
     &               XSASTL2, XSWSUM, XSTCYC, PSIGMAA, HPPTYPE, XGOOD,
     &               XPCALC, XPA, XPB, XPC, XPD
C
C double precision common block
C
      COMMON /XSIGR/ XSWGHT, XSUSRW, XSTARG, PHISTEP
C
      SAVE /XSIGI/
      SAVE /XSIGR/
C
