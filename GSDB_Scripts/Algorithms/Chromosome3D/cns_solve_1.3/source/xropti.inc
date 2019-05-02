C data structure XROPTI
C =====================
C
C dynamic data structure for B-factor refinement
C Weight on the isotropic B-factor restraint terms.
      DOUBLE PRECISION RWGHT
C logical flag for B-factor restraints
      LOGICAL QRESTR
C minimum allowed value for B-factor
      DOUBLE PRECISION BFMIN
C maximum allowed value for B-factor
      DOUBLE PRECISION BFMAX
C minimum allowed value for Q
      DOUBLE PRECISION QFMIN
C maximum allowed value for Q
      DOUBLE PRECISION QFMAX
C
C minimum allowed value for f'
      DOUBLE PRECISION FPMIN
C maximum allowed value for f'
      DOUBLE PRECISION FPMAX
C
C minimum allowed value for f''
      DOUBLE PRECISION FDPMIN
C maximum allowed value for f''
      DOUBLE PRECISION FDPMAX
C
C Number of bond restraints on the isotropic B-factors
      INTEGER NBR
C b-factor bond restraint list (HEAP pointer)
      INTEGER IIBR, IJBR
C
C Number of angle restraints on the isotropic B-factors
      INTEGER NTR
C b-factor angle restraint list (HEAP pointer)
      INTEGER IITR, IKTR
C
C Non-crystallographic symmetry restraint list (HEAP pointer)
      INTEGER IINCSX
C
C restrictions due to symmetry for overall anisotropic B parameters
      CHARACTER*4 RESTRC
C
C lowest target value and cycle number that produced it
      INTEGER MCYCLE
      DOUBLE PRECISION MTARG
C
C temporary arrays (HEAP pointers)
      INTEGER ITEMP, IFLAG, ILIST, IPOINT, IOMAIN, IQOLD, ITYPE
      INTEGER IBSIGM, IASIGM, FPOLD, FDPOLD
      INTEGER FWMAIN, FQMAIN, FXRFP, FXRFDP
C
C number of cycles counter
      INTEGER ICYCLE
C
C overall B and k scale factor routines
      LOGICAL QKSCALE, QBSCALE, QANISO, QISO
      INTEGER PAROLD
      DOUBLE PRECISION XSCK, XSCB, XSCB11, XSCB22, XSCB33
      DOUBLE PRECISION XSCB12, XSCB13, XSCB23
      DOUBLE PRECISION XSCKO, XSCBO, XSCB11O, XSCB22O, XSCB33O
      DOUBLE PRECISION XSCB12O, XSCB13O, XSCB23O
      INTEGER OVERII, OVERPFCALC
C
C
C common block
C
C     double precision block
C
      COMMON /XROPT/  RWGHT, BFMIN, BFMAX, QFMIN, QFMAX,
     &       FPMIN, FPMAX, FDPMIN, FDPMAX, MTARG,
     &       XSCK, XSCB, XSCB11, XSCB22, XSCB33,
     &       XSCB12, XSCB13, XSCB23,
     &       XSCKO, XSCBO, XSCB11O, XSCB22O, XSCB33O,
     &       XSCB12O, XSCB13O, XSCB23O
C
C     integer block
C
      COMMON /IXROPT/
     &       NBR, IIBR, IJBR, NTR, IITR, IKTR, ITEMP, IFLAG,
     &       ILIST, IPOINT, IOMAIN, ICYCLE, IINCSX, IQOLD, ITYPE,
     &       IBSIGM, IASIGM, FPOLD, FDPOLD,
     &       FWMAIN, FQMAIN, FXRFP, FXRFDP, MCYCLE, PAROLD,
     &       OVERII, OVERPFCALC
C
C     logical block
C
      COMMON /LXROPT/ QRESTR, QKSCALE, QBSCALE, QANISO, QISO
C
C     character block
C
      COMMON /AXROPT/ RESTRC
C
      SAVE /XROPT/
      SAVE /IXROPT/
      SAVE /LXROPT/
      SAVE /AXROPT/
