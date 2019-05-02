C
C xtarget.inc
C
C Definition of target, dtarget, monitor expressions, scattering tables,
C associate atom selections and object names
C
C Author: Axel T. Brunger
C
C -------------------------------------------
C XNAMEAS must match RPNX in expression.inc!!
C -------------------------------------------
C
C update flags
C XRQCHK indicates if atomic model structure factors need to be updated
      LOGICAL XRQCHK
C
C XRUPAT indicates that there are no ASSOciate atom selections defined
C this can happend when the molecular structure database is changed
C (change of atom numbers will eleminate the ASSOciate selections).
      LOGICAL XRUPAT
C
C XRREUP is a flag that indicates that the tselection and cvselection
C expressions need to be interpreted (routine XRTEST)
      LOGICAL XRREUP
C
C associate atom selections
C XNAMEAS must match RPNX in expression.inc!!
      INTEGER MAXHPFLAG, IHPFLAG, XNAMEAS
      PARAMETER (MAXHPFLAG=9,XNAMEAS=XRNAME_SIZE)
      INTEGER HPFLAG(MAXHPFLAG)
      CHARACTER*(XNAMEAS) XASSOC(MAXHPFLAG)
C form factor atom selections
      INTEGER XRNATF, HPATOF, HPINDF
C anomalous atom flag array
      INTEGER HPANOMFLAG
C coordinate stores for "tolerance" feature
      INTEGER HPFX, HPFY, HPFZ
C derivative arrays
      INTEGER HPDX, HPDY, HPDZ, HPDT, HPDQ
C number of atoms in molecular structure file (to detect changes in NATOM)
      INTEGER XRMATO
C scale factor for xray energy function
      DOUBLE PRECISION XRSCAL
C unscaled restraint energy
      DOUBLE PRECISION XRE
C tolerance for linear approximation
      DOUBLE PRECISION XRLTOL
C scattering tables (four Gaussians (a,b), a constant (c), f', and f'')
      INTEGER XRSM, XRSN
      PARAMETER (XRSM=400)
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(XRSM), XRFDP(XRSM)
C flag for default ASELEction ( attr scatter_fdp > 0 )
      LOGICAL QASELE
C HPTSEL: pointer for selection array for target function
      INTEGER HPTSEL
C logical flag indicating the presence of a tselection selection
      LOGICAL XCVTEST
C
C command stack that defines the target expression
C ---------------------------------------------
      INTEGER TRPNMX, TRPNX
      PARAMETER (TRPNMX=800, TRPNX=XNAMEAS)
      CHARACTER*(TRPNX) TRPN(4,TRPNMX)
      INTEGER TRPNL(4,TRPNMX), TRPNN
C arpndb constains double precision/double complex constants
      DOUBLE COMPLEX TRPNDB(4,TRPNMX)
C arpnmlt is the number of arguments for functions
      INTEGER TRPNMLT(TRPNMX)
C arpntyp is the data type of the command, arpndom is the domain  type
      CHARACTER*2 TRPNTYP(TRPNMX), TRPNDOM(TRPNMX)
C arpnlev is the level of the operation
      INTEGER TRPNLEV(TRPNMX)
C adepth is the depth of the expression
      INTEGER TDEPTH
C
C command stack that defines the derivative expression
C ----------------------------------------------------
      INTEGER DRPNMX, DRPNX
      PARAMETER (DRPNMX=800, DRPNX=XNAMEAS)
      CHARACTER*(DRPNX) DRPN(4,DRPNMX,MAXHPFLAG)
      INTEGER DRPNL(4,DRPNMX,MAXHPFLAG), DRPNN(MAXHPFLAG)
C arpndb constains double precision/double complex constants
      DOUBLE COMPLEX DRPNDB(4,DRPNMX,MAXHPFLAG)
C arpnmlt is the number of arguments for functions
      INTEGER DRPNMLT(DRPNMX,MAXHPFLAG)
C arpntyp is the data type of the command, arpndom is the domain  type
      CHARACTER*2 DRPNTYP(DRPNMX,MAXHPFLAG)
      CHARACTER*2 DRPNDOM(DRPNMX,MAXHPFLAG)
C arpnlev is the level of the operation
      INTEGER DRPNLEV(DRPNMX,MAXHPFLAG)
C adepth is the depth of the expression
      INTEGER DDEPTH(MAXHPFLAG)
C
C command stack that defines the target selection expression
C ----------------------------------------------------------
      INTEGER SRPNMX, SRPNX
      PARAMETER (SRPNMX=800, SRPNX=XNAMEAS)
      CHARACTER*(SRPNX) SRPN(4,SRPNMX)
      INTEGER SRPNL(4,SRPNMX), SRPNN
C arpndb constains double precision/double complex constants
      DOUBLE COMPLEX SRPNDB(4,SRPNMX)
C arpnmlt is the number of arguments for functions
      INTEGER SRPNMLT(SRPNMX)
C arpntyp is the data type of the command, arpndom is the domain  type
      CHARACTER*2 SRPNTYP(SRPNMX), SRPNDOM(SRPNMX)
C arpnlev is the level of the operation
      INTEGER SRPNLEV(SRPNMX)
C adepth is the depth of the expression
      INTEGER SDEPTH
C
C command stack that defines the selection expression for cross-validation
C ----------------------------------------------------------
      INTEGER CRPNMX, CRPNX
      PARAMETER (CRPNMX=800, CRPNX=XNAMEAS)
      CHARACTER*(CRPNX) CRPN(4,CRPNMX)
      INTEGER CRPNL(4,CRPNMX), CRPNN
C arpndb constains double precision/double complex constants
      DOUBLE COMPLEX CRPNDB(4,CRPNMX)
C arpnmlt is the number of arguments for functions
      INTEGER CRPNMLT(CRPNMX)
C arpntyp is the data type of the command, arpndom is the domain  type
      CHARACTER*2 CRPNTYP(CRPNMX), CRPNDOM(CRPNMX)
C arpnlev is the level of the operation
      INTEGER CRPNLEV(CRPNMX)
C adepth is the depth of the expression
      INTEGER CDEPTH
C
C command stack that defines the expression for monitoring the refinement
C ----------------------------------------------------------
      INTEGER MRPNMX, MRPNX
      PARAMETER (MRPNMX=800, MRPNX=XNAMEAS)
      CHARACTER*(MRPNX) MRPN(4,MRPNMX)
      INTEGER MRPNL(4,MRPNMX), MRPNN
C arpndb constains double precision/double complex constants
      DOUBLE COMPLEX MRPNDB(4,MRPNMX)
C arpnmlt is the number of arguments for functions
      INTEGER MRPNMLT(MRPNMX)
C arpntyp is the data type of the command, arpndom is the domain  type
      CHARACTER*2 MRPNTYP(MRPNMX), MRPNDOM(MRPNMX)
C arpnlev is the level of the operation
      INTEGER MRPNLEV(MRPNMX)
C adepth is the depth of the expression
      INTEGER MDEPTH
C
C     double precision common block
C
      COMMON /DXTARG/ TRPNDB, DRPNDB, SRPNDB, CRPNDB, MRPNDB,
     &                XRSCAL, XRE, XRLTOL, XRSA, XRSB, XRSC,
     &                XRFP, XRFDP
C
C     integer common block
C
      COMMON /IXTARG/ TRPNL, TRPNN, TRPNMLT, TRPNLEV, TDEPTH,
     &                DRPNL, DRPNN, DRPNMLT, DRPNLEV, DDEPTH,
     &                SRPNL, SRPNN, SRPNMLT, SRPNLEV, SDEPTH,
     &                CRPNL, CRPNN, CRPNMLT, CRPNLEV, CDEPTH,
     &                MRPNL, MRPNN, MRPNMLT, MRPNLEV, MDEPTH,
     &                XRMATO, XRNATF,HPATOF,HPINDF,
     &                HPDX, HPDY, HPDZ, HPDT, HPDQ,
     &                HPFX, HPFY, HPFZ, HPANOMFLAG,
     &                IHPFLAG, HPFLAG, XRSN, HPTSEL
C
C     character common block
C
      COMMON /SXTARG/ TRPN, TRPNTYP, TRPNDOM,
     &                DRPN, DRPNTYP, DRPNDOM,
     &                SRPN, SRPNTYP, SRPNDOM,
     &                CRPN, CRPNTYP, CRPNDOM,
     &                MRPN, MRPNTYP, MRPNDOM,
     &                XASSOC
C
C     logical common block
C
      COMMON /LXTARG/ QASELE, XCVTEST, XRREUP, XRUPAT, XRQCHK
C
      SAVE /DXTARG/
      SAVE /IXTARG/
      SAVE /SXTARG/
      SAVE /LXTARG/
