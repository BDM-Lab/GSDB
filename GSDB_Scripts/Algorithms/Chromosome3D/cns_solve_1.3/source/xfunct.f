C======================================================================
      SUBROUTINE XDOFUNC(OK,NAME,LENGTH,VARS,NARGS,QALL,DEPTH)
C
C Parsing routine for functions.  This simply
C checks if the current word matches a function definition.
C If a match is found the routine returns the name of the
C function and optional parameters (NAME(1), NAME(2), NAME(3),
C NAME(4)), the length of the NAME strings (LENGTH(1), LENGTH(2),
C LENGTH(3), LENGTH(4)), the parameters (VARS(1),VARS(2),VARS(3),
C VARS(4)), the number of arguments (NARGS), and the depth
C of the function (DEPTH).
C
C If QALL is true then all parameters are parsed (e.g.,
C  RESId 1:2).  If QALL is false, the parameters are not
C parsed.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      LOGICAL OK
      CHARACTER*(*) NAME(4)
      INTEGER LENGTH(4)
      DOUBLE COMPLEX VARS(4)
      INTEGER NARGS
      LOGICAL QALL
      INTEGER DEPTH
C local
      INTEGER LMAX, LPROMPT
      LOGICAL OKL, ERR
      DOUBLE PRECISION NEGW
      CHARACTER*(WDMAX) PROMPT
      INTEGER VARNUM, I
      DOUBLE PRECISION NEG, ELEM1, ELEM2
C parameter
      DOUBLE PRECISION MARK
      DOUBLE PRECISION ZERO, ONE, TWO, THREE
      PARAMETER (MARK=-9999.0D0)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
C begin
      OK=.TRUE.
      NAME(1)=WD
      NAME(2)=' '
      NAME(3)=' '
      NAME(4)=' '
      LMAX=LEN(NAME(1))
      LENGTH(1)=MIN(LMAX,WDLEN)
      LENGTH(2)=1
      LENGTH(3)=1
      LENGTH(4)=1
      VARS(1)=DCMPLX(ZERO,ZERO)
      VARS(2)=DCMPLX(ZERO,ZERO)
      VARS(3)=DCMPLX(ZERO,ZERO)
      VARS(4)=DCMPLX(ZERO,ZERO)
      DEPTH=0
C
C single-argument functions.
      IF (WD(1:WDLEN).EQ.'EXP') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'LOG') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'LOG10') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'SIN') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'COS') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'TAN') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'TANH') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'ASIN') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'ACOS') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'ATAN') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'I0') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'I1') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'I1OVERI0') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'I1OVERI0INV') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'SQRT') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'SIGN') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'STEP'.OR.WD(1:WDLEN).EQ.'HEAVY') THEN
      NARGS=1
      NAME(1)='STEP'
      LENGTH(1)=4
      ELSEIF (WD(1:WDLEN).EQ.'INT') THEN
      NAME(1)='INTEGER'
      LENGTH(1)=7
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'REAL') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'IMAG') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'ABS') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'AMPLITUDE') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'PHASE') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'CONJUGATE') THEN
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'GAUSS') THEN
      NAME(1)='GAUSSIAN'
      LENGTH(1)=8
      NARGS=1
      ELSEIF (WD(1:WDLEN).EQ.'DISTRIBUTE') THEN
      NARGS=1
C
C two-argument functions
      ELSEIF (WD(1:WDLEN).EQ.'MOD') THEN
      NARGS=2
      ELSEIF (WD(1:WDLEN).EQ.'COMPLEX') THEN
      NARGS=2
      ELSEIF (WD(1:WDLEN).EQ.'COMBINE') THEN
      NARGS=2
C
C multi-argument functions (indicated by the "-1")
      ELSEIF (WD(1:WDLEN).EQ.'MIN') THEN
      NARGS=-1
      ELSEIF (WD(1:WDLEN).EQ.'MAX') THEN
      NARGS=-1
C
C int/add/mult variable
      ELSEIF (WD(1:1).EQ.'_'.AND.WDLEN.GT.1) THEN
      IF (QALL) THEN
      NAME(1)='_VAR_'
      LENGTH(1)=5
C
      CALL COPYST(NAME(2),LMAX,LENGTH(2),WD,WDLEN)
      END IF
      NARGS=0
C
C integral, sum, and product functions
      ELSEIF (WD(1:WDLEN).EQ.'INTEGRATE'
     &    .OR.WD(1:WDLEN).EQ.'MAXIMIZE'
     &    .OR.WD(1:WDLEN).EQ.'IMAXIMIZE'
     &    .OR.WD(1:WDLEN).EQ.'ADD'
     &    .OR.WD(1:WDLEN).EQ.'MULTIPLY')  THEN
C
      IF (QALL) THEN
C
      ERR=.FALSE.
C
      CALL NEXTDO('int/add/mult=')
      IF (WD(1:WDLEN).NE.'[') THEN
      ERR=.TRUE.
      CALL DSPERR('XDOFUNC','Expecting "[".')
      END IF
C
C variable name
      CALL NEXTDO('Variable=')
      IF (WD(1:1).NE.'_') THEN
      ERR=.TRUE.
      CALL DSPERR('XDOFUNC','Variable name has to begin with "_"')
      END IF
      IF (WDLEN.LE.1) THEN
      ERR=.TRUE.
      CALL DSPERR('XDOFUNC','Incorrect variable name.')
      END IF
      CALL COPYST(NAME(2),LMAX,LENGTH(2),WD,WDLEN)
C
C start value
      CALL NEXTDO('int/add/mult=')
      IF (WD(1:WDLEN).NE.',') THEN
      ERR=.TRUE.
      CALL DSPERR('XDOFUNC','Expecting ",".')
      END IF
      CALL NEXTDX('Lower-bound=')
      VARS(1)=DCMPLX(DECODF(WD,WDLEN,OKL),ZERO)
      IF (.NOT.OKL) THEN
      CALL DSPERR('XDOFUNC','Error converting lower bound.')
      ERR=.TRUE.
      END IF
C
C stop value
      CALL NEXTDO('int/add/mult=')
      IF (WD(1:WDLEN).NE.',') THEN
      ERR=.TRUE.
      CALL DSPERR('XDOFUNC','Expecting ",".')
      END IF
      CALL NEXTDX('Upper-bound=')
      VARS(2)=DCMPLX(DECODF(WD,WDLEN,OKL),ZERO)
      IF (.NOT.OKL) THEN
      CALL DSPERR('XDOFUNC','Error converting upper bound.')
      ERR=.TRUE.
      END IF
C
C increment
      CALL NEXTDO('int/add/mult=')
      IF (WD(1:WDLEN).NE.',') THEN
      ERR=.TRUE.
      CALL DSPERR('XDOFUNC','Expecting ",".')
      END IF
      CALL NEXTDX('Step=')
      VARS(3)=DCMPLX(DECODF(WD,WDLEN,OKL),ZERO)
      IF (.NOT.OKL) THEN
      CALL DSPERR('XDOFUNC','Error converting step.')
      ERR=.TRUE.
      END IF
      IF (DBLE(VARS(3)).LE.ZERO) THEN
      CALL DSPERR('XDOFUNC','step has to be positive.')
      ERR=.TRUE.
      END IF
C
      CALL NEXTDO('int/add/mult=')
      IF (WD(1:WDLEN).NE.']') THEN
      ERR=.TRUE.
      CALL DSPERR('XDOFUNC','Expecting "]".')
      END IF
C
C
      IF (ERR) THEN
      CALL DSPERR('XDOFUNC',
     & 'Incorrect definition.  Must be [<name>,<real>,<real>,<real>]')
      END IF
      END IF
C
      NARGS=1
C
C depth of the function is 2 (for summation of the IMAXIMIZE function)
      IF (NAME(1).EQ.'IMAXIMIZE') THEN
      DEPTH=2
      ELSE
C depth of the function is 1 (for summation of all other functions)
      DEPTH=1
      END IF
C
C functions without arguments and no auxiliary fields
      ELSEIF (WD(1:WDLEN).EQ.'RANDOM') THEN
      NAME(1)='RANDOM'
      LENGTH(1)=6
      NARGS=0
C parse old syntax "()"
      IF (QALL) THEN
      CALL NEXTDO('XDO>')
      IF (WD(1:WDLEN).EQ.'(') THEN
      CALL NEXTDO('XDO>')
      ELSE
      CALL SAVEWD
      END IF
      END IF
C
      ELSEIF (WD(1:WDLEN).EQ.'I') THEN
      NARGS=0
C
C structure-factor-specific functions
      ELSEIF (WD(1:WDLEN).EQ.'CENTRIC_PHASE') THEN
      NAME(1)='CENTRIC_P'
      LENGTH(1)=9
      NARGS=0
C
      ELSEIF (WD(1:WDLEN).EQ.'S') THEN
      NARGS=0
C
      ELSEIF (WD(1:WDLEN).EQ.'H') THEN
      NARGS=0
C
      ELSEIF (WD(1:WDLEN).EQ.'K') THEN
      NARGS=0
C
      ELSEIF (WD(1:WDLEN).EQ.'L') THEN
      NARGS=0
C
      ELSEIF (WD(1:WDLEN).EQ.'D') THEN
      NARGS=0
C
      ELSEIF (WD(1:WDLEN).EQ.'CENTRIC') THEN
      NARGS=0
C
      ELSEIF (WD(1:WDLEN).EQ.'ACENTRIC') THEN
      NARGS=0
C
      ELSEIF (WD(1:WDLEN).EQ.'FRIEDEL_PAIR') THEN
      NAME(1)='FRIEDEL_PA'
      LENGTH(1)=10
      NARGS=1
C
      ELSEIF (WD(1:WDLEN).EQ.'MULT') THEN
      NARGS=0
C
      ELSEIF (WD(1:WDLEN).EQ.'EPSILON') THEN
      NAME(1)='EPSI'
      LENGTH(1)=4
      NARGS=0
C
      ELSEIF (WD(1:WDLEN).EQ.'TYPE') THEN
      NARGS=0
C
C MAP-specific functions
      ELSEIF (WD(1:WDLEN).EQ.'X') THEN
      NARGS=0
      ELSEIF (WD(1:WDLEN).EQ.'Y') THEN
      NARGS=0
      ELSEIF (WD(1:WDLEN).EQ.'Z') THEN
      NARGS=0
      ELSEIF (WD(1:WDLEN).EQ.'A') THEN
      NARGS=0
      ELSEIF (WD(1:WDLEN).EQ.'B') THEN
      NARGS=0
      ELSEIF (WD(1:WDLEN).EQ.'C') THEN
      NARGS=0
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      ELSEIF (WD(1:WDLEN).EQ.'SDEV') THEN
      NARGS=1
C=====================================================================
C #endif
C=====================================================================
      PROMPT=WD(1:WDLEN)//'>'
      LPROMPT=WDLEN+1
      IF (QALL) THEN
      VARS(1)=DCMPLX(ZERO,ZERO)
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'[') THEN
      DO WHILE (WD(1:WDLEN).NE.']')
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:1).EQ.',') THEN
      CONTINUE
      ELSEIF (WD(1:1).EQ.']') THEN
      CONTINUE
      ELSEIF (WD(1:4).EQ.'RADI') THEN
      CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'=') CALL NEXTDX(PROMPT(1:LPROMPT))
      VARS(1)=DCMPLX(DECODF(WD,WDLEN,OKL),ZERO)
      IF (.NOT.OKL) THEN
      CALL DSPERR('XDOFUNC','Real number expected. ')
      END IF
      ELSEIF (WD(1:4).EQ.'MEAN') THEN
      CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'=') CALL NEXTDX(PROMPT(1:LPROMPT))
      VARS(2)=DCMPLX(DECODF(WD,WDLEN,OKL),ZERO)
      IF (.NOT.OKL) THEN
      CALL DSPERR('XDOFUNC','Real number expected. ')
      END IF
      ELSE
      CALL DSPERR('XDOFUNC','Unknown parameter. ')
      END IF
      END DO
      ELSE
      CALL SAVEWD
      END IF
      END IF
C
      ELSEIF (WD(1:WDLEN).EQ.'GAVE') THEN
      NARGS=2
C
C Fourier transform operations
      ELSEIF (WD(1:WDLEN).EQ.'FT') THEN
      NARGS=1
C
C sum and average functions
      ELSEIF (WD(1:WDLEN).EQ.'SAVE'.OR.
     &        WD(1:WDLEN).EQ.'AVE'.OR.
     &        WD(1:WDLEN).EQ.'SUM') THEN
      NARGS=1
C
      VARS(2)=DCMPLX(ONE,ZERO)
      PROMPT=WD(1:WDLEN)//'>'
      LPROMPT=WDLEN+1
      IF (QALL) THEN
      VARS(1)=DCMPLX(ZERO,ZERO)
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'[') THEN
      DO WHILE (WD(1:WDLEN).NE.']')
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:1).EQ.',') THEN
      CONTINUE
      ELSEIF (WD(1:1).EQ.']') THEN
      CONTINUE
      ELSEIF (WD(1:WDLEN).EQ.'OVERALL') THEN
      VARS(2)=DCMPLX(-ONE,ZERO)
      ELSEIF (WD(1:WDLEN).EQ.'BINWISE') THEN
      VARS(2)=DCMPLX(ONE,ZERO)
      ELSE
      CALL DSPERR('XDOFUNC','Unknown parameter. ')
      END IF
      END DO
      ELSE
      CALL SAVEWD
      END IF
      END IF
C
C
C normalization function
      ELSEIF (WD(1:WDLEN).EQ.'NORM') THEN
      NARGS=1
C
      ELSEIF (WD(1:WDLEN).EQ.'FRIEDEL') THEN
      NARGS=1
C
      ELSEIF (WD(1:WDLEN).EQ.'SIGA'.OR.WD(1:WDLEN).EQ.'SIGACV'.OR.
     &        WD(1:WDLEN).EQ.'SCALE'.OR.
     &        WD(1:WDLEN).EQ.'SHAPE') THEN
      NARGS=2
C
      IF (WD(1:WDLEN).EQ.'SIGACV') THEN
      PROMPT=WD(1:WDLEN)//'>'
      LPROMPT=WDLEN+1
      IF (QALL) THEN
      VARS(1)=DCMPLX(ZERO,ZERO)
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'[') THEN
      DO WHILE (WD(1:WDLEN).NE.']')
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:1).EQ.',') THEN
      CONTINUE
      ELSEIF (WD(1:1).EQ.']') THEN
      CONTINUE
      ELSEIF (WD(1:4).EQ.'SIGM') THEN
      NEGW=ONE
      CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'=') CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'-') THEN
         CALL NEXTDX(PROMPT(1:LPROMPT))
         NEGW=-ONE
      END IF
      VARS(1)=NEGW*DCMPLX(DECODF(WD,WDLEN,OKL),ZERO)
      IF (.NOT.OKL) THEN
      CALL DSPERR('XDOFUNC','Real number expected. ')
      END IF
      ELSE
      CALL DSPERR('XDOFUNC','Unknown parameter. ')
      END IF
      END DO
      ELSE
      CALL SAVEWD
      END IF
      END IF
      END IF
C
      ELSEIF (WD(1:WDLEN).EQ.'CORR') THEN
      NARGS=2
C
      PROMPT=WD(1:WDLEN)//'>'
      LPROMPT=WDLEN+1
      IF (QALL) THEN
      VARS(1)=DCMPLX(ZERO,ZERO)
      VARS(2)=DCMPLX(ONE,ZERO)
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'[') THEN
      DO WHILE (WD(1:WDLEN).NE.']')
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:1).EQ.',') THEN
      CONTINUE
      ELSEIF (WD(1:1).EQ.']') THEN
      CONTINUE
      ELSEIF (WD(1:4).EQ.'MULT') THEN
      CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'=') CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'TRUE') THEN
      VARS(1)=DCMPLX(ONE,ZERO)
      ELSEIF (WD(1:WDLEN).EQ.'FALSE') THEN
      VARS(1)=DCMPLX(ZERO,ZERO)
      ELSE
      CALL DSPERR('XDOFUNC','TRUE or FALSE expected. ')
      END IF
      ELSEIF (WD(1:WDLEN).EQ.'OVERALL') THEN
      VARS(2)=DCMPLX(-ONE,ZERO)
      ELSEIF (WD(1:WDLEN).EQ.'BINWISE') THEN
      VARS(2)=DCMPLX(ONE,ZERO)
      ELSE
      CALL DSPERR('XDOFUNC','Unknown parameter. ')
      END IF
      END DO
      ELSE
      CALL SAVEWD
      END IF
      END IF
C
C R-value
      ELSEIF (WD(1:WDLEN).EQ.'RVALUE') THEN
      NARGS=2
C
      PROMPT=WD(1:WDLEN)//'>'
      LPROMPT=WDLEN+1
      IF (QALL) THEN
      VARS(1)=DCMPLX(-ONE,ZERO)
      VARS(2)=DCMPLX(ONE,ZERO)
      VARS(3)=DCMPLX(ZERO,ZERO)
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'[') THEN
      DO WHILE (WD(1:WDLEN).NE.']')
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:1).EQ.',') THEN
      CONTINUE
      ELSEIF (WD(1:1).EQ.']') THEN
      CONTINUE
      ELSEIF (WD(1:1).EQ.'K') THEN
      NEGW=ONE
      CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'=') CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'-') THEN
         CALL NEXTDX(PROMPT(1:LPROMPT))
         NEGW=-ONE
      END IF
      VARS(1)=NEGW*DCMPLX(DECODF(WD,WDLEN,OKL),ZERO)
      IF (.NOT.OKL) THEN
      CALL DSPERR('XDOFUNC','Real number expected. ')
      END IF
      ELSEIF (WD(1:WDLEN).EQ.'OVERALL') THEN
      VARS(2)=DCMPLX(-ONE,ZERO)
      ELSEIF (WD(1:WDLEN).EQ.'BINWISE') THEN
      VARS(2)=DCMPLX(ONE,ZERO)
      ELSEIF (WD(1:4).EQ.'MULT') THEN
      CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'=') CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'TRUE') THEN
      VARS(3)=DCMPLX(ONE,ZERO)
      ELSEIF (WD(1:WDLEN).EQ.'FALSE') THEN
      VARS(3)=DCMPLX(ZERO,ZERO)
      ELSE
      CALL DSPERR('XDOFUNC','TRUE or FALSE expected. ')
      END IF
      ELSE
      CALL DSPERR('XDOFUNC','Unknown parameter. ')
      END IF
      END DO
      ELSE
      CALL SAVEWD
      END IF
      END IF
C
C maximum likelihood functions
      ELSEIF (WD(1:WDLEN).EQ.'MLF'.OR.WD(1:WDLEN).EQ.'DMLF'.OR.
     &        WD(1:WDLEN).EQ.'MLFF'.OR.WD(1:WDLEN).EQ.'MLVF'.OR.
     &        WD(1:WDLEN).EQ.'MLI'.OR.WD(1:WDLEN).EQ.'DMLI' .OR.
     &        WD(1:WDLEN).EQ.'MLHL'.OR.WD(1:WDLEN).EQ.'DMLHL') THEN
      IF (WD(1:WDLEN).EQ.'MLF'.OR.WD(1:WDLEN).EQ.'DMLF'
     &.OR.WD(1:WDLEN).EQ.'MLFF'.OR.WD(1:WDLEN).EQ.'MLVF'
     &.OR.WD(1:WDLEN).EQ.'MLI'.OR.WD(1:WDLEN).EQ.'DMLI') THEN
         NARGS=5
      END IF
      IF (WD(1:WDLEN).EQ.'MLHL'.OR.WD(1:WDLEN).EQ.'DMLHL') THEN
         NARGS=10
      ENDIF
C
      IF (WD(1:WDLEN).EQ.'MLHL'.OR.WD(1:WDLEN).EQ.'DMLHL') THEN
      PROMPT=WD(1:WDLEN)//'>'
      LPROMPT=WDLEN+1
      IF (QALL) THEN
      VARS(1)=DCMPLX(ZERO,ZERO)
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'[') THEN
      DO WHILE (WD(1:WDLEN).NE.']')
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:1).EQ.',') THEN
      CONTINUE
      ELSEIF (WD(1:1).EQ.']') THEN
      CONTINUE
      ELSEIF (WD(1:3).EQ.'PHI') THEN
      CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'=') CALL NEXTDX(PROMPT(1:LPROMPT))
      VARS(1)=DCMPLX(DECODF(WD,WDLEN,OKL),ZERO)
      IF (.NOT.OKL) THEN
      CALL DSPERR('XDOFUNC','Real number expected. ')
      END IF
C
      ELSE
      CALL DSPERR('XDOFUNC','Unknown parameter. ')
      END IF
      END DO
      ELSE
      CALL SAVEWD
      END IF
      END IF
      END IF
C
C correlation functions
      ELSEIF (WD(1:WDLEN).EQ.'E2E2'.OR.WD(1:WDLEN).EQ.'DE2E2'.OR.
     &        WD(1:WDLEN).EQ.'E1E1'.OR.WD(1:WDLEN).EQ.'DE1E1'.OR.
     &        WD(1:WDLEN).EQ.'F2F2'.OR.WD(1:WDLEN).EQ.'DF2F2'.OR.
     &        WD(1:WDLEN).EQ.'F1F1'.OR.WD(1:WDLEN).EQ.'DF1F1') THEN
      NARGS=2
C
      IF (WD(1:WDLEN).EQ.'F2F2'.OR.WD(1:WDLEN).EQ.'DF2F2'.OR.
     &    WD(1:WDLEN).EQ.'F1F1'.OR.WD(1:WDLEN).EQ.'DF1F1') THEN
      PROMPT=WD(1:WDLEN)//'>'
      LPROMPT=WDLEN+1
      IF (QALL) THEN
      VARS(1)=DCMPLX(ZERO,ZERO)
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'[') THEN
      DO WHILE (WD(1:WDLEN).NE.']')
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:1).EQ.',') THEN
      CONTINUE
      ELSEIF (WD(1:1).EQ.']') THEN
      CONTINUE
      ELSEIF (WD(1:4).EQ.'MULT') THEN
      CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'=') CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'TRUE') THEN
      VARS(1)=DCMPLX(ONE,ZERO)
      ELSEIF (WD(1:WDLEN).EQ.'FALSE') THEN
      VARS(1)=DCMPLX(ZERO,ZERO)
      ELSE
      CALL DSPERR('XDOFUNC','TRUE or FALSE expected. ')
      END IF
C
      ELSE
      CALL DSPERR('XDOFUNC','Unknown parameter. ')
      END IF
      END DO
      ELSE
      CALL SAVEWD
      END IF
      END IF
      END IF
C
C residual target functions
      ELSEIF (WD(1:WDLEN).EQ.'RESI'.OR.
     &        WD(1:WDLEN).EQ.'DRESI'.OR.
     &        WD(1:WDLEN).EQ.'VECTOR'.OR.
     &        WD(1:WDLEN).EQ.'DVECTOR') THEN
      NARGS=3
C
      PROMPT=WD(1:WDLEN)//'>'
      LPROMPT=WDLEN+1
      IF (QALL) THEN
      VARS(1)=DCMPLX(ZERO,ZERO)
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'[') THEN
      DO WHILE (WD(1:WDLEN).NE.']')
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:1).EQ.',') THEN
      CONTINUE
      ELSEIF (WD(1:1).EQ.']') THEN
      CONTINUE
      ELSEIF (WD(1:1).EQ.'K') THEN
      CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'=') CALL NEXTDX(PROMPT(1:LPROMPT))
      VARS(1)=DCMPLX(DECODF(WD,WDLEN,OKL),ZERO)
      IF (.NOT.OKL) THEN
      CALL DSPERR('XDOFUNC','Real number expected. ')
      END IF
C
      ELSE
      CALL DSPERR('XDOFUNC','Unknown parameter. ')
      END IF
      END DO
      ELSE
      CALL SAVEWD
      END IF
      END IF
C
      ELSEIF (WD(1:WDLEN).EQ.'GET_FOM'.OR.
     &        WD(1:WDLEN).EQ.'GET_NORM'.OR.
     &        WD(1:WDLEN).EQ.'GET_ML'.OR.
     &        WD(1:WDLEN).EQ.'GET_DML'.OR.
     &        WD(1:WDLEN).EQ.'GET_DSML'.OR.
     &        WD(1:WDLEN).EQ.'GET_AML'.OR.
     &        WD(1:WDLEN).EQ.'GET_DAML'.OR.
     &        WD(1:WDLEN).EQ.'GET_DSAML') THEN
      IF (NAME(1).EQ.'GET_FOM'.OR.NAME(1).EQ.'GET_NORM') NARGS=4
      IF (NAME(1).EQ.'GET_ML'.OR.NAME(1).EQ.'GET_DML') NARGS=9
      IF (NAME(1).EQ.'GET_AML'.OR.NAME(1).EQ.'GET_DAML') NARGS=9
      IF (NAME(1).EQ.'GET_DSML'.OR.NAME(1).EQ.'GET_DSAML') NARGS=9
C
      PROMPT=WD(1:WDLEN)//'>'
      LPROMPT=WDLEN+1
      IF (QALL) THEN
      VARS(1)=DCMPLX(ZERO,ZERO)
      VARS(2)=DCMPLX(ZERO,ZERO)
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'[') THEN
      DO WHILE (WD(1:WDLEN).NE.']')
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:1).EQ.',') THEN
      CONTINUE
      ELSEIF (WD(1:1).EQ.']') THEN
      CONTINUE
      ELSEIF (WD(1:3).EQ.'PHI') THEN
      CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'=') CALL NEXTDX(PROMPT(1:LPROMPT))
      VARS(1)=DCMPLX(DECODF(WD,WDLEN,OKL),ZERO)
      IF (.NOT.OKL) THEN
      CALL DSPERR('XDOFUNC','Real number expected. ')
      END IF
      ELSEIF (WD(1:4).EQ.'CEN3') THEN
      CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'=') CALL NEXTDX(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'TRUE') THEN
      VARS(2)=DCMPLX(ONE,ZERO)
      ELSEIF (WD(1:WDLEN).EQ.'FALSE') THEN
      VARS(2)=DCMPLX(ZERO,ZERO)
      ELSE
      CALL DSPERR('XDOFUNC','TRUE or FALSE expected. ')
      END IF
C
      ELSE
      CALL DSPERR('XDOFUNC','Unknown parameter. ')
      END IF
      END DO
      ELSE
      CALL SAVEWD
      END IF
      END IF
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      ELSEIF (WD(1:WDLEN).EQ.'REMAP') THEN
      NARGS=1
C=====================================================================
C #endif
C=====================================================================
      PROMPT=WD(1:WDLEN)//'>'
      LPROMPT=WDLEN+1
      IF (QALL) THEN
      VARS(1)=DCMPLX(ONE,  ZERO)
      VARS(2)=DCMPLX(TWO,  ZERO)
      VARS(3)=DCMPLX(THREE,ZERO)
      VARNUM=1
      NEG=ONE
      ELEM1=MARK
      ELEM2=MARK
      CALL NEXTDO(PROMPT(1:LPROMPT))
      IF (WD(1:WDLEN).EQ.'[') THEN
      DO WHILE (WD(1:WDLEN).NE.']')
      CALL NEXTDO(PROMPT(1:LPROMPT))
C
      DO I=1,WDLEN
      IF (WD(I:I).EQ.'-') THEN
      NEG=-(ONE)
      ELSEIF (WD(I:I).EQ.'+') THEN
      NEG=ONE
      ELSEIF (WD(I:I).EQ.'H') THEN
      IF (ELEM1.EQ.MARK) THEN
         ELEM1=NEG*ONE
      ELSEIF (ELEM2.EQ.MARK) THEN
         ELEM2=NEG*ONE
      END IF
      NEG=ONE
      ELSEIF (WD(I:I).EQ.'K') THEN
      IF (ELEM1.EQ.MARK) THEN
         ELEM1=NEG*TWO
      ELSEIF (ELEM2.EQ.MARK) THEN
         ELEM2=NEG*TWO
      END IF
      NEG=ONE
      ELSEIF (WD(I:I).EQ.'L') THEN
      IF (ELEM1.EQ.MARK) THEN
         ELEM1=NEG*THREE
      ELSEIF (ELEM2.EQ.MARK) THEN
         ELEM2=NEG*THREE
      END IF
      NEG=ONE
      ELSEIF (WD(I:I).EQ.',') THEN
      VARS(VARNUM)=DCMPLX(ELEM1,ELEM2)
      VARNUM=VARNUM+1
      ELEM1=MARK
      ELEM2=MARK
      ELSEIF (WD(I:I).EQ.']') THEN
      VARS(VARNUM)=DCMPLX(ELEM1,ELEM2)
      ELEM1=MARK
      ELEM2=MARK
      IF (VARNUM.NE.3) THEN
         CALL WRNDIE(-5,'XDOFUNC','Incomplete remapping expression')
      END IF
      ELSE
      CALL DSPERR('XDOFUNC','Unknown parameter. ')
      END IF
      END DO
C
      END DO
      ELSE
      CALL SAVEWD
      END IF
      END IF
C function not found
      ELSE
      OK=.FALSE.
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE XDOTYPE(ERR,WD,WDLEN,NL,NARGS,
     &                   ATYPE,FTYPE,ADOMAIN,FDOMAIN,QHERM,
     &                   XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
C
C Type and domain definition for all functions, operations,
C and operands.
C
C id   data type
C --------------
C ST   string
C LO   logical
C DP   double precision
C DC   double complex
C
C id   domain type
C ----------------
C '  '   single element (number or string or logical)
C SF     structure factor
C MP     map
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      LOGICAL ERR
      CHARACTER*(*) WD
      INTEGER NL, WDLEN, NARGS
      CHARACTER*2 ATYPE(*), FTYPE, ADOMAIN(*), FDOMAIN
      LOGICAL QHERM
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*), XRHOTYP(*)
C local
      INTEGER I
      LOGICAL OK
C begin
C
      FTYPE='??'
      FDOMAIN='??'
C
C single-argument functions.
      IF    (
     &       WD(1:WDLEN).EQ.'EXP'
     &   .OR.WD(1:WDLEN).EQ.'LOG'
     &   .OR.WD(1:WDLEN).EQ.'LOG10'
     &   .OR.WD(1:WDLEN).EQ.'SQRT'
     &   .OR.WD(1:WDLEN).EQ.'CONJUGATE'
     &   .OR.WD(1:WDLEN).EQ.'CHS'
     &      )  THEN
      IF (ATYPE(NL).EQ.'DP') FTYPE='DP'
      IF (ATYPE(NL).EQ.'DC') FTYPE='DC'
      IF (ADOMAIN(NL).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'MP') FDOMAIN='MP'
      ELSEIF (
     &       WD(1:WDLEN).EQ.'SIN'
     &   .OR.WD(1:WDLEN).EQ.'COS'
     &   .OR.WD(1:WDLEN).EQ.'TAN'
     &   .OR.WD(1:WDLEN).EQ.'TANH'
     &   .OR.WD(1:WDLEN).EQ.'ASIN'
     &   .OR.WD(1:WDLEN).EQ.'ACOS'
     &   .OR.WD(1:WDLEN).EQ.'ATAN'
     &   .OR.WD(1:WDLEN).EQ.'I0'
     &   .OR.WD(1:WDLEN).EQ.'I1'
     &   .OR.WD(1:WDLEN).EQ.'I1OVERI0'
     &   .OR.WD(1:WDLEN).EQ.'I1OVERI0INV'
     &   .OR.WD(1:WDLEN).EQ.'SIGN'
     &   .OR.WD(1:WDLEN).EQ.'STEP'
     &   .OR.WD(1:WDLEN).EQ.'HEAVY'
     &   .OR.WD(1:WDLEN).EQ.'INTEGER'
     &   .OR.WD(1:WDLEN).EQ.'GAUSSIAN'
     &       ) THEN
C
      IF (ATYPE(NL).EQ.'DP') FTYPE='DP'
      IF (ADOMAIN(NL).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'MP') FDOMAIN='MP'
C distribute function
      ELSEIF (WD(1:WDLEN).EQ.'DISTRIBUTE') THEN
C
      IF (ATYPE(NL).EQ.'DP') FTYPE='DP'
      IF (ATYPE(NL).EQ.'DC') FTYPE='DC'
      IF (ADOMAIN(NL).EQ.'SF') FDOMAIN='SF'
C
      ELSEIF (
     &       WD(1:WDLEN).EQ.'REAL'
     &   .OR.WD(1:WDLEN).EQ.'IMAG'
     &   .OR.WD(1:WDLEN).EQ.'ABS'
     &   .OR.WD(1:WDLEN).EQ.'AMPLITUDE'
     &   .OR.WD(1:WDLEN).EQ.'PHASE'
     &                               ) THEN
      IF (ATYPE(NL).EQ.'DP') FTYPE='DP'
      IF (ATYPE(NL).EQ.'DC') FTYPE='DP'
      IF (ADOMAIN(NL).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'MP') FDOMAIN='MP'
C
      ELSEIF (WD(1:WDLEN).EQ.'NOT') THEN
      IF (ATYPE(NL).EQ.'LO') FTYPE='LO'
      IF (ADOMAIN(NL).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'MP') FDOMAIN='MP'
C
C Fourier transformation function
      ELSEIF (WD(1:WDLEN).EQ.'FT') THEN
      IF (QHERM.AND.ADOMAIN(NL).EQ.'MP'.AND.ATYPE(NL).EQ.'DP') THEN
      FTYPE='DC'
      FDOMAIN='SF'
      ELSEIF(.NOT.QHERM.AND.ADOMAIN(NL).EQ.'MP'.AND.ATYPE(NL).EQ.'DC')
     &   THEN
      FTYPE='DC'
      FDOMAIN='SF'
      ELSEIF (QHERM.AND.ADOMAIN(NL).EQ.'SF'.AND.ATYPE(NL).EQ.'DC') THEN
      FTYPE='DP'
      FDOMAIN='MP'
      ELSEIF (.NOT.QHERM.AND.ADOMAIN(NL).EQ.'SF'.AND.ATYPE(NL).EQ.'DC')
     &    THEN
      FTYPE='DC'
      FDOMAIN='MP'
      END IF
C
C special structure factor functions
      ELSEIF (WD(1:4).EQ.'NORM') THEN
      IF (ATYPE(NL).EQ.'DP') FTYPE='DP'
      IF (ADOMAIN(NL).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  ') FDOMAIN='  '
      ELSEIF (WD(1:4).EQ.'SAVE'.OR.WD(1:3).EQ.'AVE'.OR.
     &        WD(1:3).EQ.'SUM'.OR.WD(1:WDLEN).EQ.'FRIEDEL') THEN
      IF (ATYPE(NL).EQ.'DP') FTYPE='DP'
      IF (ATYPE(NL).EQ.'DC') FTYPE='DC'
      IF (ADOMAIN(NL).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  ') FDOMAIN='  '
      ELSEIF (WD(1:10).EQ.'FRIEDEL_PA') THEN
      IF (ATYPE(NL).EQ.'LO') FTYPE='LO'
      IF (ADOMAIN(NL).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  ') FDOMAIN='  '
C
C maximum likelihood functions
      ELSEIF (WD(1:WDLEN).EQ.'MLF'.OR.WD(1:WDLEN).EQ.'DMLF'.OR.
     &        WD(1:WDLEN).EQ.'MLFF'.OR.WD(1:WDLEN).EQ.'MLVF'.OR.
     &        WD(1:WDLEN).EQ.'MLI'.OR.WD(1:WDLEN).EQ.'DMLI') THEN
      IF ((ADOMAIN(NL).EQ.'SF'.OR.ADOMAIN(NL).EQ.'  ').AND.
     &    (ADOMAIN(NL-1).EQ.'SF'.OR.ADOMAIN(NL-1).EQ.'  ').AND.
     &    (ADOMAIN(NL-2).EQ.'SF'.OR.ADOMAIN(NL-2).EQ.'  ').AND.
     &    (ADOMAIN(NL-3).EQ.'SF'.OR.ADOMAIN(NL-3).EQ.'  ').AND.
     &    (ADOMAIN(NL-4).EQ.'SF'.OR.ADOMAIN(NL-4).EQ.'  ')
     &     ) FDOMAIN='SF'
      IF (WD(1:WDLEN).EQ.'DMLF'.OR.WD(1:WDLEN).EQ.'DMLI') THEN
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP'.AND.
     &    ATYPE(NL-2).EQ.'DC'.AND.ATYPE(NL-3).EQ.'DP'.AND.
     &    ATYPE(NL-4).EQ.'DP') FTYPE='DC'
      ELSE
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP'.AND.
     &    ATYPE(NL-2).EQ.'DC'.AND.ATYPE(NL-3).EQ.'DP'.AND.
     &    ATYPE(NL-4).EQ.'DP') FTYPE='DP'
      END IF
      ELSEIF (WD(1:WDLEN).EQ.'MLHL'.OR.WD(1:WDLEN).EQ.'DMLHL') THEN
      IF ((ADOMAIN(NL).EQ.'SF'.OR.ADOMAIN(NL).EQ.'  ').AND.
     &    (ADOMAIN(NL-1).EQ.'SF'.OR.ADOMAIN(NL-1).EQ.'  ').AND.
     &    (ADOMAIN(NL-2).EQ.'SF'.OR.ADOMAIN(NL-2).EQ.'  ').AND.
     &    (ADOMAIN(NL-3).EQ.'SF'.OR.ADOMAIN(NL-3).EQ.'  ').AND.
     &    (ADOMAIN(NL-4).EQ.'SF'.OR.ADOMAIN(NL-4).EQ.'  ').AND.
     &    (ADOMAIN(NL-5).EQ.'SF'.OR.ADOMAIN(NL-5).EQ.'  ').AND.
     &    (ADOMAIN(NL-6).EQ.'SF'.OR.ADOMAIN(NL-6).EQ.'  ').AND.
     &    (ADOMAIN(NL-7).EQ.'SF'.OR.ADOMAIN(NL-7).EQ.'  ').AND.
     &    (ADOMAIN(NL-8).EQ.'SF'.OR.ADOMAIN(NL-8).EQ.'  ')
     &     ) FDOMAIN='SF'
      IF (WD(1:WDLEN).EQ.'MLHL') THEN
        IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP'.AND.
     &      ATYPE(NL-2).EQ.'DP'.AND.ATYPE(NL-3).EQ.'DP'.AND.
     &      ATYPE(NL-4).EQ.'DP'.AND.ATYPE(NL-5).EQ.'DP'.AND.
     &      ATYPE(NL-6).EQ.'DC'.AND.ATYPE(NL-7).EQ.'DP'.AND.
     &      ATYPE(NL-8).EQ.'DP') FTYPE='DP'
      ELSE
        IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP'.AND.
     &      ATYPE(NL-2).EQ.'DP'.AND.ATYPE(NL-3).EQ.'DP'.AND.
     &      ATYPE(NL-4).EQ.'DP'.AND.ATYPE(NL-5).EQ.'DP'.AND.
     &      ATYPE(NL-6).EQ.'DC'.AND.ATYPE(NL-7).EQ.'DP'.AND.
     &      ATYPE(NL-8).EQ.'DP') FTYPE='DC'
      ENDIF
C
C correlation target functions
      ELSEIF (WD(1:WDLEN).EQ.'E2E2'.OR.WD(1:WDLEN).EQ.'DE2E2'.OR.
     &        WD(1:WDLEN).EQ.'E1E1'.OR.WD(1:WDLEN).EQ.'DE1E1'.OR.
     &        WD(1:WDLEN).EQ.'F2F2'.OR.WD(1:WDLEN).EQ.'DF2F2'.OR.
     &        WD(1:WDLEN).EQ.'F1F1'.OR.WD(1:WDLEN).EQ.'DF1F1') THEN
      IF ((ADOMAIN(NL).EQ.'SF'.OR.ADOMAIN(NL).EQ.'  ').AND.
     &    (ADOMAIN(NL-1).EQ.'SF'.OR.ADOMAIN(NL-1).EQ.'  ')
     &     ) FDOMAIN='SF'
      IF (WD(1:WDLEN).EQ.'DE2E2'.OR.
     &    WD(1:WDLEN).EQ.'DE1E1'.OR.
     &    WD(1:WDLEN).EQ.'DF2F2'.OR.
     &    WD(1:WDLEN).EQ.'DF1F1') THEN
      IF (ATYPE(NL).EQ.'DC'.AND.ATYPE(NL-1).EQ.'DP') FTYPE='DC'
      ELSE
      IF (ATYPE(NL).EQ.'DC'.AND.ATYPE(NL-1).EQ.'DP') FTYPE='DP'
      END IF
C
C residual target functions
      ELSEIF (WD(1:WDLEN).EQ.'RESI'.OR.
     &        WD(1:WDLEN).EQ.'DRESI'.OR.
     &        WD(1:WDLEN).EQ.'VECTOR'.OR.
     &        WD(1:WDLEN).EQ.'DVECTOR') THEN
      IF ((ADOMAIN(NL).EQ.'SF'.OR.ADOMAIN(NL).EQ.'  ').AND.
     &    (ADOMAIN(NL-1).EQ.'SF'.OR.ADOMAIN(NL-1).EQ.'  ').AND.
     &    (ADOMAIN(NL-2).EQ.'SF'.OR.ADOMAIN(NL-2).EQ.'  ')
     &     ) FDOMAIN='SF'
      IF (WD(1:WDLEN).EQ.'DRESI') THEN
      IF (ATYPE(NL).EQ.'DP'.AND.
     &    ATYPE(NL-1).EQ.'DC'.AND.
     &    ATYPE(NL-2).EQ.'DP') FTYPE='DC'
      ELSEIF (WD(1:WDLEN).EQ.'RESI') THEN
      IF (ATYPE(NL).EQ.'DP'.AND.
     &    ATYPE(NL-1).EQ.'DC'.AND.
     &    ATYPE(NL-2).EQ.'DP') FTYPE='DP'
      ELSEIF (WD(1:WDLEN).EQ.'DVECTOR') THEN
      IF (ATYPE(NL).EQ.'DP'.AND.
     &    ATYPE(NL-1).EQ.'DC'.AND.
     &    ATYPE(NL-2).EQ.'DC') FTYPE='DC'
      ELSEIF (WD(1:WDLEN).EQ.'VECTOR') THEN
      IF (ATYPE(NL).EQ.'DP'.AND.
     &    ATYPE(NL-1).EQ.'DC'.AND.
     &    ATYPE(NL-2).EQ.'DC') FTYPE='DP'
      END IF
C phase probability functions
      ELSEIF (WD(1:7).EQ.'GET_FOM') THEN
      IF ((ADOMAIN(NL).EQ.'SF'.OR.ADOMAIN(NL).EQ.'  ').AND.
     &    (ADOMAIN(NL-1).EQ.'SF'.OR.ADOMAIN(NL-1).EQ.'  ').AND.
     &    (ADOMAIN(NL-2).EQ.'SF'.OR.ADOMAIN(NL-2).EQ.'  ').AND.
     &    (ADOMAIN(NL-3).EQ.'SF'.OR.ADOMAIN(NL-3).EQ.'  ')
     &     ) FDOMAIN='SF'
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP'.AND.
     &    ATYPE(NL-2).EQ.'DP'.AND.ATYPE(NL-3).EQ.'DP') FTYPE='DC'
      ELSEIF (WD(1:8).EQ.'GET_NORM') THEN
      IF ((ADOMAIN(NL).EQ.'SF'.OR.ADOMAIN(NL).EQ.'  ').AND.
     &    (ADOMAIN(NL-1).EQ.'SF'.OR.ADOMAIN(NL-1).EQ.'  ').AND.
     &    (ADOMAIN(NL-2).EQ.'SF'.OR.ADOMAIN(NL-2).EQ.'  ').AND.
     &    (ADOMAIN(NL-3).EQ.'SF'.OR.ADOMAIN(NL-3).EQ.'  ')
     &     ) FDOMAIN='SF'
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP'.AND.
     &    ATYPE(NL-2).EQ.'DP'.AND.ATYPE(NL-3).EQ.'DP') FTYPE='DP'
C
      ELSEIF (WD(1:6).EQ.'GET_ML'.OR.WD(1:7).EQ.'GET_AML') THEN
      IF ((ADOMAIN(NL).EQ.'SF'.OR.ADOMAIN(NL).EQ.'  ').AND.
     &    (ADOMAIN(NL-1).EQ.'SF'.OR.ADOMAIN(NL-1).EQ.'  ').AND.
     &    (ADOMAIN(NL-2).EQ.'SF'.OR.ADOMAIN(NL-2).EQ.'  ').AND.
     &    (ADOMAIN(NL-3).EQ.'SF'.OR.ADOMAIN(NL-3).EQ.'  ').AND.
     &    (ADOMAIN(NL-4).EQ.'SF'.OR.ADOMAIN(NL-4).EQ.'  ').AND.
     &    (ADOMAIN(NL-5).EQ.'SF'.OR.ADOMAIN(NL-5).EQ.'  ').AND.
     &    (ADOMAIN(NL-6).EQ.'SF'.OR.ADOMAIN(NL-6).EQ.'  ').AND.
     &    (ADOMAIN(NL-7).EQ.'SF'.OR.ADOMAIN(NL-7).EQ.'  ').AND.
     &    (ADOMAIN(NL-8).EQ.'SF'.OR.ADOMAIN(NL-8).EQ.'  ')
     &     ) FDOMAIN='SF'
CCC modification ATB 4/27/08
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP'.AND.
     &    ATYPE(NL-2).EQ.'DP'.AND.ATYPE(NL-3).EQ.'DP'.AND.
     &    ATYPE(NL-4).EQ.'DP'.AND.ATYPE(NL-5).EQ.'DP'.AND.
     &    ATYPE(NL-6)(1:1).EQ.'D'.AND.ATYPE(NL-7).EQ.'DC'.AND.
     &    ATYPE(NL-8)(1:1).EQ.'D') FTYPE='DP'
C
      ELSEIF (WD(1:7).EQ.'GET_DML'.OR.WD(1:8).EQ.'GET_DAML'.OR.
     &  WD(1:8).EQ.'GET_DSML'.OR.WD(1:9).EQ.'GET_DSAML') THEN
      IF ((ADOMAIN(NL).EQ.'SF'.OR.ADOMAIN(NL).EQ.'  ').AND.
     &    (ADOMAIN(NL-1).EQ.'SF'.OR.ADOMAIN(NL-1).EQ.'  ').AND.
     &    (ADOMAIN(NL-2).EQ.'SF'.OR.ADOMAIN(NL-2).EQ.'  ').AND.
     &    (ADOMAIN(NL-3).EQ.'SF'.OR.ADOMAIN(NL-3).EQ.'  ').AND.
     &    (ADOMAIN(NL-4).EQ.'SF'.OR.ADOMAIN(NL-4).EQ.'  ').AND.
     &    (ADOMAIN(NL-5).EQ.'SF'.OR.ADOMAIN(NL-5).EQ.'  ').AND.
     &    (ADOMAIN(NL-6).EQ.'SF'.OR.ADOMAIN(NL-6).EQ.'  ').AND.
     &    (ADOMAIN(NL-7).EQ.'SF'.OR.ADOMAIN(NL-7).EQ.'  ').AND.
     &    (ADOMAIN(NL-8).EQ.'SF'.OR.ADOMAIN(NL-8).EQ.'  ')
     &     ) FDOMAIN='SF'
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP'.AND.
     &    ATYPE(NL-2).EQ.'DP'.AND.ATYPE(NL-3).EQ.'DP'.AND.
     &    ATYPE(NL-4).EQ.'DP'.AND.ATYPE(NL-5).EQ.'DP'.AND.
     &    ATYPE(NL-6)(1:1).EQ.'D'.OR.ATYPE(NL-7).EQ.'DC'.AND.
     &    ATYPE(NL-8)(1:1).EQ.'D') FTYPE='DC'
C special two-argument structure factor functions
      ELSEIF (WD(1:WDLEN).EQ.'CORR'.OR.WD(1:WDLEN).EQ.'SIGA'.OR.
     &        WD(1:WDLEN).EQ.'SIGACV') THEN
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP') FTYPE='DP'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
C
C R-value
      ELSEIF (WD(1:WDLEN).EQ.'RVALUE') THEN
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP') FTYPE='DP'
      IF (ATYPE(NL).EQ.'DC'.AND.ATYPE(NL-1).EQ.'DP') FTYPE='DP'
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DC') FTYPE='DP'
      IF (ATYPE(NL).EQ.'DC'.AND.ATYPE(NL-1).EQ.'DC') FTYPE='DP'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
C
      ELSEIF (
     &        WD(1:WDLEN).EQ.'SCALE'
     &   .OR. WD(1:WDLEN).EQ.'SHAPE'
     &        ) THEN
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP') FTYPE='DP'
      IF (ATYPE(NL).EQ.'DC'.AND.ATYPE(NL-1).EQ.'DC') FTYPE='DP'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
C
C remapping function
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      ELSEIF (WD(1:WDLEN).EQ.'REMAP') THEN
      IF (ATYPE(NL).EQ.'LO') FTYPE='LO'
      IF (ATYPE(NL).EQ.'DP') FTYPE='DP'
      IF (ATYPE(NL).EQ.'DC') FTYPE='DC'
      IF (ADOMAIN(NL).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  ') FDOMAIN='SF'
C=====================================================================
C #endif
C=====================================================================
C
C two-argument functions
      ELSEIF (WD(1:WDLEN).EQ.'MOD') THEN
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP') FTYPE='DP'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
C
      ELSEIF (
     &        WD(1:WDLEN).EQ.'COMPLEX'
     &   .OR. WD(1:WDLEN).EQ.'COMBINE'
     &        ) THEN
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP') FTYPE='DC'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
C
C integration variable
      ELSEIF (WD(1:5).EQ.'_VAR_') THEN
      FTYPE='DP'
      FDOMAIN='  '
C
C integration function
      ELSEIF (
     &        WD(1:WDLEN).EQ.'INTEGRATE'
     &    .OR.WD(1:WDLEN).EQ.'ADD'
     &    .OR.WD(1:WDLEN).EQ.'MULTIPLY'
     &        ) THEN
      IF (ATYPE(NL).EQ.'DP') FTYPE='DP'
      IF (ATYPE(NL).EQ.'DC') FTYPE='DC'
      IF (ADOMAIN(NL).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'MP') FDOMAIN='MP'
C
C maximize function
      ELSEIF (WD(1:WDLEN).EQ.'MAXIMIZE'
     &    .OR.WD(1:WDLEN).EQ.'IMAXIMIZE'
     &        ) THEN
      IF (ATYPE(NL).EQ.'DP') FTYPE='DP'
      IF (ADOMAIN(NL).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'MP') FDOMAIN='MP'
C
C multi-argument functions
      ELSEIF (
     &           WD(1:WDLEN).EQ.'MIN'
     &      .OR. WD(1:WDLEN).EQ.'MAX'
     &        ) THEN
      DO I=1,NARGS
      ERR=ERR.OR.ATYPE(NL-I+1).NE.'DP'
      END DO
      IF (.NOT.ERR) FTYPE='DP'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
C
C
C functions without arguments and no auxiliary fields
C Miller indices
      ELSEIF (
     &       WD(1:WDLEN).EQ.'H'
     &   .OR.WD(1:WDLEN).EQ.'K'
     &   .OR.WD(1:WDLEN).EQ.'L'
     &       ) THEN
      FTYPE='DP'
      FDOMAIN='SF'
C
      ELSEIF (WD(1:WDLEN).EQ.'RANDOM') THEN
      FTYPE='DP'
      FDOMAIN='  '
C
      ELSEIF (WD(1:WDLEN).EQ.'I') THEN
      FTYPE='DC'
      FDOMAIN='  '
C
C structure-factor-specific functions
      ELSEIF (WD(1:WDLEN).EQ.'CENTRIC_P') THEN
      FTYPE='DP'
      FDOMAIN='SF'
C
      ELSEIF (WD(1:WDLEN).EQ.'S') THEN
      FTYPE='DP'
      FDOMAIN='SF'
C
      ELSEIF (WD(1:WDLEN).EQ.'MULT') THEN
      FTYPE='DP'
      FDOMAIN='  '
C
      ELSEIF (WD(1:WDLEN).EQ.'EPSI') THEN
      FTYPE='DP'
      FDOMAIN='SF'
C
      ELSEIF (WD(1:WDLEN).EQ.'TYPE') THEN
      FTYPE='DP'
      FDOMAIN='SF'
C
      ELSEIF (WD(1:WDLEN).EQ.'D') THEN
      FTYPE='DP'
      FDOMAIN='SF'
C
      ELSEIF (WD(1:WDLEN).EQ.'CENTRIC') THEN
      FTYPE='LO'
      FDOMAIN='SF'
C
      ELSEIF (WD(1:WDLEN).EQ.'ACENTRIC') THEN
      FTYPE='LO'
      FDOMAIN='SF'
C
      ELSEIF (WD(1:10).EQ.'FRIEDEL_PA') THEN
      FTYPE='LO'
      FDOMAIN='SF'
C
C MAP-specific functions
      ELSEIF (WD(1:WDLEN).EQ.'X') THEN
      FTYPE='DP'
      FDOMAIN='MP'
C
      ELSEIF (WD(1:WDLEN).EQ.'Y') THEN
      FTYPE='DP'
      FDOMAIN='MP'
C
      ELSEIF (WD(1:WDLEN).EQ.'Z') THEN
      FTYPE='DP'
      FDOMAIN='MP'
C
      ELSEIF (WD(1:WDLEN).EQ.'A') THEN
      FTYPE='DP'
      FDOMAIN='MP'
C
      ELSEIF (WD(1:WDLEN).EQ.'B') THEN
      FTYPE='DP'
      FDOMAIN='MP'
C
      ELSEIF (WD(1:WDLEN).EQ.'C') THEN
      FTYPE='DP'
      FDOMAIN='MP'
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      ELSEIF (WD(1:WDLEN).EQ.'SDEV') THEN
      FTYPE='DP'
      FDOMAIN='MP'
C=====================================================================
C #endif
C=====================================================================
      ELSEIF (WD(1:WDLEN).EQ.'GAVE') THEN
      FTYPE='DP'
      FDOMAIN='MP'
C
C
C numerical operations
      ELSEIF (
     &             WD(1:WDLEN).EQ.'+'
     &         .OR.WD(1:WDLEN).EQ.'-'
     &         .OR.WD(1:WDLEN).EQ.'*'
     &         .OR.WD(1:WDLEN).EQ.'/'
     &         .OR.WD(1:WDLEN).EQ.'^'
     &       ) THEN
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP') FTYPE='DP'
      IF (ATYPE(NL).EQ.'DC'.AND.ATYPE(NL-1).EQ.'DP') FTYPE='DC'
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DC') FTYPE='DC'
      IF (ATYPE(NL).EQ.'DC'.AND.ATYPE(NL-1).EQ.'DC') FTYPE='DC'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
C
C
C logical operations
      ELSEIF (
     &         WD(1:WDLEN).EQ.'AND'
     &     .OR.WD(1:WDLEN).EQ.'OR'
     &       ) THEN
      IF (ATYPE(NL).EQ.'LO'.AND.ATYPE(NL-1).EQ.'LO') FTYPE='LO'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
C
C relational operations
      ELSEIF (
     &             WD(1:WDLEN).EQ.'='
     &         .OR.WD(1:WDLEN).EQ.'#'
     &       ) THEN
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP') FTYPE='LO'
      IF (ATYPE(NL).EQ.'DC'.AND.ATYPE(NL-1).EQ.'DP') FTYPE='LO'
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DC') FTYPE='LO'
      IF (ATYPE(NL).EQ.'DC'.AND.ATYPE(NL-1).EQ.'DC') FTYPE='LO'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
C
      ELSEIF (
     &             WD(1:WDLEN).EQ.'>'
     &         .OR.WD(1:WDLEN).EQ.'<'
     &         .OR.WD(1:WDLEN).EQ.'<='
     &         .OR.WD(1:WDLEN).EQ.'>='
     &       ) THEN
      IF (ATYPE(NL).EQ.'DP'.AND.ATYPE(NL-1).EQ.'DP') FTYPE='LO'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='  '
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'SF'.AND.ADOMAIN(NL-1).EQ.'SF') FDOMAIN='SF'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'  ') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'  '.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
      IF (ADOMAIN(NL).EQ.'MP'.AND.ADOMAIN(NL-1).EQ.'MP') FDOMAIN='MP'
C
      ELSE
C
C check reciprocal space objects
      OK=.FALSE.
      DO I=1,XSFNUM
      IF (WD(1:WDLEN).EQ.XSFNAM(I)) THEN
      OK=.TRUE.
      IF (XSFTYPE(I).EQ.'COMP') THEN
      FTYPE='DC'
      FDOMAIN='SF'
      ELSEIF (XSFTYPE(I).EQ.'REAL') THEN
      FTYPE='DP'
      FDOMAIN='SF'
CCC      ELSEIF (XSFTYPE(I).EQ.'PHAS') THEN
CCC      FTYPE='DP'
CCC      FDOMAIN='SF'
      ELSEIF (XSFTYPE(I).EQ.'INTE') THEN
      FTYPE='DP'
      FDOMAIN='SF'
      END IF
      END IF
      END DO
C
      IF (.NOT.OK) THEN
C check real space objects
      DO I=1,XRHONUM
      IF (WD(1:WDLEN).EQ.XRHONAM(I)) THEN
      OK=.TRUE.
      IF (QHERM) THEN
      FTYPE='DP'
      ELSE
      FTYPE='DC'
      END IF
      FDOMAIN='MP'
      END IF
      END DO
      END IF
C
C function not found
      IF (.NOT.OK) THEN
      ERR=.TRUE.
      CALL WRNDIE(-5,'XDOTYPE',
     & 'Internal error: function,operation,or operand type undefined.')
      END IF
      END IF
C
      IF (  FTYPE.EQ.'??'.OR.
     &    FDOMAIN.EQ.'??') THEN
      CALL DSPERR('XDOTYPE','Variable/type mismatch')
      ERR=.TRUE.
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOSPCL(RPNMX,RPNX,RPNN,RPN,RPNL,QSPEC)
C
C Routine checks if a special operation is present.
C A special operation implies that all structure factor or map
C elements are required for the operation (e.g.,
C normalization of structure factors).
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER RPNMX, RPNX
      CHARACTER*(*) RPN(4,RPNMX)
      INTEGER RPNL(4,RPNMX), RPNN
      LOGICAL QSPEC
C local
      INTEGER NN
C begin
      QSPEC=.FALSE.
      DO NN=1,RPNN
      IF (RPN(1,NN)(1:RPNL(1,NN)).EQ.'GAVE'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'SCALE'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'SHAPE'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'NORM'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'SAVE'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'AVE' .OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'SUM' .OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'FRIEDEL' .OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'FRIEDEL_PA' .OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'CORR'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'RVALUE'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'SIGA'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'SIGACV'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'SDEV'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'GET_AML'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'GET_DAML'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'GET_DSAML'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'E2E2'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'DE2E2'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'E1E1'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'DE1E1'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'F2F2'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'DF2F2'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'F1F1'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'DF1F1'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'RESI'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'DRESI'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'VECTOR'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'DVECTOR'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'REMAP'.OR.
     &    RPN(1,NN)(1:RPNL(1,NN)).EQ.'DISTRIBUTE') THEN
      QSPEC=.TRUE.
      END IF
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XDOOPER(OK,NAME,LENGTH,QALL,
     &                   XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
C
C Parsing routine for structure factor operands.  This simply
C checks if the current word matches a function definition.
C If a match is found the routine returns the name of the
C name of the operand and optional parameters
C (NAME(1), NAME(2), NAME(3), NAME(4)), and the length of the NAME
C strings (LENGTH(1), LENGTH(2), LENGTH(3), LENGTH(4)).
C
C If QALL is true then all parameters are parsed (e.g.,
C  RESId 1:2).  If QALL is false, the parameters are not
C parsed.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      LOGICAL OK
      CHARACTER*(*) NAME(4)
      INTEGER LENGTH(4)
      LOGICAL QALL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*), XRHOTYP(*)
C local
      INTEGER I
C begin
      NAME(1)=WD
      NAME(2)=' '
      NAME(3)=' '
      NAME(4)=' '
      LENGTH(1)=WDLEN
      LENGTH(2)=1
      LENGTH(3)=1
      LENGTH(4)=1
C
C reciprocal space objects
C ========================
      OK=.FALSE.
      DO I=1,XSFNUM
      IF (WD(1:WDLEN).EQ.XSFNAM(I)) THEN
      OK=.TRUE.
      END IF
      END DO
C
      IF (.NOT.OK) THEN
C
C map operands
C ============
      OK=.FALSE.
      DO I=1,XRHONUM
      IF (WD(1:WDLEN).EQ.XRHONAM(I)) THEN
      OK=.TRUE.
      END IF
      END DO
      END IF
C
      RETURN
      END
C======================================================================
