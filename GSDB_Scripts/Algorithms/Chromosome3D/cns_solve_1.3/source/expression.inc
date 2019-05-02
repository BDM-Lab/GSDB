C
C expression.inc
C
C Data structure for command stacks that are generated
C by "expression.s" using Reverse Polish Notation.
C
C Author: Axel T. Brunger
C
C
C first rpn command stack
C -----------------------
      INTEGER RPNMX, RPNX
C XNAMEAS (xtarget.inc) must match RPNX !!
C XNAMEMX (xcrystal.inc) must match RPNX !!
      PARAMETER (RPNMX=200, RPNX=XRNAME_SIZE)
      CHARACTER*(RPNX) RPN(4,RPNMX)
      INTEGER RPNL(4,RPNMX), RPNN
C rpndb constains double precision/double complex constants
      DOUBLE COMPLEX RPNDB(4,RPNMX)
C rpnmlt is the number of arguments for functions
      INTEGER RPNMLT(RPNMX)
C rpntyp is the data type of the command
      CHARACTER*2 RPNTYP(RPNMX), RPNDOM(RPNMX)
C rpnlev is the level of the operation
      INTEGER RPNLEV(RPNMX)
C
C second rpn command stack
C ------------------------
      CHARACTER*(RPNX) RPN2(4,RPNMX)
      INTEGER RPNL2(4,RPNMX), RPNN2
C rpndb constains double precision/double complex constants
      DOUBLE COMPLEX RPNDB2(4,RPNMX)
C rpnmlt is the number of arguments for functions
      INTEGER RPNMLT2(RPNMX)
C rpntyp is the data type of the command
      CHARACTER*2 RPNTYP2(RPNMX), RPNDOM2(RPNMX)
C rpnlev is the level of the operation
      INTEGER RPNLEV2(RPNMX)
C
C third rpn command stack
C ------------------------
      CHARACTER*(RPNX) RPN3(4,RPNMX)
      INTEGER RPNL3(4,RPNMX), RPNN3
C rpndb constains double precision/double complex constants
      DOUBLE COMPLEX RPNDB3(4,RPNMX)
C rpnmlt is the number of arguments for functions
      INTEGER RPNMLT3(RPNMX)
C rpntyp is the data type of the command
      CHARACTER*2 RPNTYP3(RPNMX), RPNDOM3(RPNMX)
C rpnlev is the level of the operation
      INTEGER RPNLEV3(RPNMX)
C---------------------------------------------------------------------
C
C     double precision common block
C
      COMMON /DEXPRE/ RPNDB, RPNDB2, RPNDB3
C
C     integer common block
C
      COMMON /IEXPRE/ RPNL, RPNL2, RPNL3, RPNN, RPNN2, RPNN3,
     &                RPNMLT, RPNMLT2, RPNMLT3, RPNLEV, RPNLEV2,
     &                RPNLEV3
C
C     logical common block
C
CCCC      COMMON /LEXPRE/
C
C     character common block
C
      COMMON /SEXPRE/ RPN, RPN2, RPN3, RPNTYP, RPNTYP2, RPNTYP3,
     &                RPNDOM, RPNDOM2, RPNDOM3
C
