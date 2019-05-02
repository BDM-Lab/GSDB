C==================================================
      SUBROUTINE XCHVAR( K, FOUND )
      IMPLICIT NONE
C
C     check to see if variable is valid, if so get type
C
      LOGICAL FOUND
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INTEGER J, LAST, I, K
C
      FOUND = .FALSE.
C
C this 20 is okay, since the constants are only 20 chars - WLD 980105
      DO J = 1,20
         IF ( VAREQ(K)(J:J) .NE. ' ' ) LAST = J
      ENDDO
C
C     check for show
C
      I = 31
C modification ATB 5/28/08
      IF ( LAST .GE. VUNIQ(I) .AND.LAST.LE.VNAMLN) THEN
         IF ( VAREQ(K)(1:LAST) .EQ. VNAMES(I)(1:LAST) ) THEN
            VAREQ(K) = VNAMES(I)
            FOUND = .TRUE.
            VARTYP(K) = VTYP(I)
            VARNUM(K) = VNUM(I)
            RETURN
         ENDIF
      ENDIF
C
C
      RETURN
      END
      SUBROUTINE XRABS( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      INTEGER K, L
C
C     take the absolute value of the table value for index k,
C     store the result in index k-1,
C     and change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there are two possible absolute operations:
C                         abs(double complex)
C                         abs(double precision)
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( L - K .NE. 1 .OR. ( VARTYP(K) .NE. 'DP' .AND. VARTYP(K)
     &       .NE. 'DC' ) ) THEN
         CALL DSPERR('XRABS>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
      IF ( VARTYP(K) .EQ. 'DP' ) THEN
         RVAREQ(K-1) = ABS( RVAREQ(K) )
         VARTYP(K-1)   = 'DP'
      ELSE
         CVAREQ(K-1) = ABS( CVAREQ(K) )
         VARTYP(K-1)   = 'DC'
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRACOS( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      DOUBLE PRECISION RADS
      LOGICAL XRERR
      INTEGER K, L
C parameter
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C
C     take the arc cosine of the table value for index k,
C     store the result in index k-1,
C     and change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there is one possible arc cosine operations:
C                         arc cosine(double precision)
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( L - K .NE. 1 .OR. VARTYP(K) .NE. 'DP'
     &    .OR.RVAREQ(K).GT.ONE.OR.RVAREQ(K).LT.-ONE ) THEN
         CALL DSPERR('XRACOS>','Illegal argument(s).')
         XRERR = .TRUE.
         RETURN
      ENDIF
      RADS = RVAREQ(K)
      RADS = ACOS( RADS )
      CALL RD2DEG(RADS,RADS)
      RVAREQ(K-1) = RADS
      VARTYP(K-1)   = 'DP'
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRADD( K, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'symbol.inc'
C
      INTEGER K
      CHARACTER*(SVARMX) CTMP
      LOGICAL XRERR
C
C     add the table value for index k to that for index k-1, store
C     the result in index k-1, modify vareq to show it is a result,
C     and if needed change the vartyp for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there are five possible additions:
C                         double complex + double complex
C                         double complex + double precision
C                         double precision + double precision
C                         double precision + double complex
C                         string + string (append)
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( VARTYP(K-1) .EQ. 'DP' .AND. VARTYP(K) .EQ. 'DP' ) THEN
         RVAREQ(K-1) = RVAREQ(K-1) + RVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DP' .AND. VARTYP(K) .EQ. 'DC' ) THEN
         CVAREQ(K-1) = RVAREQ(K-1) + CVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DC' .AND. VARTYP(K) .EQ. 'DC' ) THEN
         CVAREQ(K-1) = CVAREQ(K-1) + CVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DC' .AND. VARTYP(K) .EQ. 'DP' ) THEN
         CVAREQ(K-1) = CVAREQ(K-1) + RVAREQ(K)
         VARTYP(K-1)   = 'DC'
      ELSEIF ( VARTYP(K-1) .EQ. 'ST' .AND. VARTYP(K) .EQ. 'ST' ) THEN
C
C        cat the strings
C
         CTMP = SVAREQ(K-1)(1:SVARLN(K-1))//SVAREQ(K)(1:SVARLN(K))
         SVAREQ(K-1) = CTMP
         SVARLN(K-1) = SVARLN(K-1)+SVARLN(K)
      ELSE
         CALL DSPERR('XRADD>','Mismatched data types.')
         XRERR = .TRUE.
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRASIN( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      DOUBLE PRECISION RADS
      INTEGER K, L
C parameter
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C
C     take the arc sine of the table value for index k,
C     store the result in index k-1,
C     and change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there is one possible arc sine operations:
C                         arc sine(double precision)
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( L - K .NE. 1 .OR. VARTYP(K) .NE. 'DP'.OR.
     &     RVAREQ(K).GT.ONE.OR.RVAREQ(K).LT.-ONE ) THEN
         CALL DSPERR('XRASIN>','Illegal argument(s).')
         XRERR = .TRUE.
         RETURN
      ENDIF
      RADS = RVAREQ(K)
      RADS = ASIN( RADS )
      CALL RD2DEG( RADS, RADS )
      RVAREQ(K-1) = RADS
      VARTYP(K-1)   = 'DP'
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRCHFX( MODE, K, FOUND )
      IMPLICIT NONE
C
C     check to see if variable is valid, if so get type
C
      LOGICAL FOUND
      CHARACTER*4 MODE
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INTEGER I, J, K, LAST
      INTEGER NUMFXS
      PARAMETER(NUMFXS=38)
      CHARACTER*(VNAMLN) FXS(NUMFXS)
      CHARACTER*1 FXTEO(NUMFXS)
      INTEGER FXUNIQ(NUMFXS)
      INTEGER ISTRT
C
CCC      FXUNIQ(1)=4
      FXUNIQ(2)=4
      FXUNIQ(3)=1
      FXUNIQ(4)=1
      FXUNIQ(5)=1
      FXUNIQ(6)=1
      FXUNIQ(7)=3
      FXUNIQ(8)=4
      FXUNIQ(9)=3
      FXUNIQ(10)=3
      FXUNIQ(11)=3
      FXUNIQ(12)=3
      FXUNIQ(13)=4
      FXUNIQ(14)=3
      FXUNIQ(15)=3
      FXUNIQ(16)=4
      FXUNIQ(17)=3
      FXUNIQ(18)=4
      FXUNIQ(19)=4
      FXUNIQ(20)=4
      FXUNIQ(21)=3
      FXUNIQ(22)=3
      FXUNIQ(23)=5
      FXUNIQ(24)=6
      FXUNIQ(25)=6
      FXUNIQ(26)=3
      FXUNIQ(27)=4
      FXUNIQ(28)=4
      FXUNIQ(29)=4
      FXUNIQ(30)=10
      FXUNIQ(31)=3
      FXUNIQ(32)=6
      FXUNIQ(33)=5
      FXUNIQ(34)=6
      FXUNIQ(35)=6
      FXUNIQ(36)=7
      FXUNIQ(37)=4
      FXUNIQ(38)=5
C
CCC      FXTEO(1)='N'
      FXTEO(2)='W'
      FXTEO(3)='K'
      FXTEO(4)='U'
      FXTEO(5)='O'
      FXTEO(6)='P'
      FXTEO(7)='A'
      FXTEO(8)='H'
      FXTEO(9)='M'
      FXTEO(10)='Z'
      FXTEO(11)='L'
      FXTEO(12)='E'
      FXTEO(13)='G'
      FXTEO(14)='S'
      FXTEO(15)='C'
      FXTEO(16)='J'
      FXTEO(17)='T'
      FXTEO(18)='F'
      FXTEO(19)='Q'
      FXTEO(20)='D'
      FXTEO(21)='I'
      FXTEO(22)='R'
      FXTEO(23)='B'
      FXTEO(24)='V'
      FXTEO(25)='X'
      FXTEO(26)='Y'
      FXTEO(27)='H'
      FXTEO(28)='a'
      FXTEO(29)='b'
      FXTEO(30)='c'
      FXTEO(31)='l'
      FXTEO(32)='s'
      FXTEO(33)='i'
      FXTEO(34)='r'
      FXTEO(35)='f'
      FXTEO(36)='g'
      FXTEO(37)='t'
      FXTEO(38)='m'
C
CCC      FXS(1)='NORM                '
      FXS(2)='MAXWELL             '
      FXS(3)='S                   '
      FXS(4)='H                   '
      FXS(5)='K                   '
      FXS(6)='L                   '
      FXS(7)='ABS                 '
      FXS(8)='HEAVY               '
      FXS(9)='MINIMUM             '
      FXS(10)='MAXIMUM             '
      FXS(11)='LOG                 '
      FXS(12)='EXP                 '
      FXS(13)='GAUSSIAN            '
      FXS(14)='SINE                '
      FXS(15)='COSINE              '
      FXS(16)='ACOSINE             '
      FXS(17)='TANGENT             '
      FXS(18)='ASINE               '
      FXS(19)='SQRT                '
      FXS(20)='SIGN                '
      FXS(21)='INTEGER             '
      FXS(22)='MOD                 '
      FXS(23)='LOG10               '
      FXS(24)='ENCODE              '
      FXS(25)='DECODE              '
      FXS(26)='RANDOM              '
      FXS(27)='STEP                '
      FXS(28)='IMOD                '
      FXS(29)='NINT                '
      FXS(30)='CAPITALIZE          '
C New string functions
      FXS(31)='LEN                 '
      FXS(32)='SUBSTR              '
      FXS(33)='INDEX               '
      FXS(34)='RINDEX              '
      FXS(35)='FORMAT              '
      FXS(36)='ADJUSTL             '
      FXS(37)='TRIM                '
      FXS(38)='MATCH               '
C
      FOUND = .FALSE.
C
C     evaluate does not recognize first six functions, skip them
C
      IF ( MODE .EQ. 'EVAL' ) THEN
         ISTRT = 7
      ELSE
C
C removed NORM function ATB 12/25/06
         ISTRT = 2
      ENDIF
C
      DO I = ISTRT,NUMFXS
C
C this 20 is okay, since the constants are only 20 chars - WLD 980105
         DO J = 1,20
            IF ( VAREQ(K)(J:J) .NE. ' ' ) LAST = J
         ENDDO
C
C modification ATB 5/28/08
         IF ( LAST .GE. FXUNIQ(I).AND.LAST.LE.VNAMLN ) THEN
            IF ( VAREQ(K)(1:LAST).EQ.FXS(I)(1:LAST)) THEN
               VAREQ(K) = FXS(I)
               FOUND = .TRUE.
               TEO(K+1) = FXTEO(I)
               RETURN
            ENDIF
         ENDIF
      ENDDO
C
      RETURN
      END
      SUBROUTINE XRCLOG( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'consta.inc'
C
      LOGICAL XRERR
      INTEGER K, L
C
C     take the common log of the table value for index k,
C     store the result in index k-1,
C     and change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     common log is only applied to real numbers:
C                         log10(double precision)
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( L - K .NE. 1 .OR. VARTYP(K) .NE. 'DP' ) THEN
         CALL DSPERR('XRCLOG>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
      IF ( RVAREQ(K) .LT. RSMALL ) THEN
         CALL DSPERR('XRCLOG>','Negative argument to common log')
         XRERR = .TRUE.
      ELSE
         RVAREQ(K-1) = LOG10( RVAREQ(K) )
         VARTYP(K-1)   = 'DP'
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRCMPX( K, MODE, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      CHARACTER*4 MODE
      DOUBLE PRECISION SCRTCH, PHASE
      INTEGER K
C
C     combine two parts of complex number separated by a comma
C
      IF (VARTYP(K-1).NE.'DP'.OR. VARTYP(K).NE.'DP') THEN
         CALL DSPERR('XRCMPX>','Illegal use of comma')
         XRERR = .TRUE.
         RETURN
      ENDIF
      CVAREQ(K-1) = DCMPLX(RVAREQ(K-1),RVAREQ(K))
      VARTYP(K-1) = 'DC'
C
C     if mode was 'AMPL' or 'PHAS' or 'SCAL' take only the proper
C     part of the literal complex number
C
      IF ( MODE .EQ. 'PHAS' ) THEN
         CALL XPHASE( CVAREQ(K-1), SCRTCH, PHASE )
         CALL RD2DEG( PHASE, RVAREQ(K-1) )
         VARTYP(K-1) = 'DP'
      ELSEIF ( MODE .EQ. 'SCAL' .OR. MODE .EQ. 'AMPL' ) THEN
         CALL XPHASE( CVAREQ(K-1), RVAREQ(K-1), SCRTCH )
         VARTYP(K-1) = 'DP'
      ENDIF
      RETURN
      END
      SUBROUTINE XRCOMA( K, MODE, XRERR )
C
      IMPLICIT NONE
C
      INTEGER K
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      CHARACTER*4 MODE
      INTEGER DESLEV, PRVLEV, KPREV
C
C     for the comma operator, we want to combine two real numbers
C     surrounded by parens. into a complex number, only if they
C     are not real arguments for a function call, store results
C     in index k-1
C
C     when we have a comma, search backwards looking for the previous
C     operator that is one paren. level below the comma.  if this is
C     found and the operator is a function then the comma is an arg.
C     separator.  if no operator is one paren. level below the comma
C     or it is not a function then the comma is part of a literal
C     complex number.
C
      KPREV = K - 1
      DESLEV = VOP(K) / 10 * 10 - 10
1     PRVLEV = VOP(KPREV) / 10 * 10
      IF ( PRVLEV  .EQ. DESLEV ) THEN
C
C        if previous opr is at des lev it could be either one
C
         IF ( (TEO(KPREV) .GE. 'A' .AND. TEO(KPREV) .LE. 'Z')
     &  .OR. (TEO(KPREV) .GE. 'a' .AND. TEO(KPREV) .LE. 'z')) THEN
C
C           if this condition then this is an argument seperator and
C           not an operator, so set it to a noop arg seperator
C
            TEO(K) = '_'
         ELSE
            CALL XRCMPX( K, MODE, XRERR )
         ENDIF
      ELSEIF ( PRVLEV .GT. DESLEV ) THEN
C
C        if prev opr is greater than desired level, keep working back
C
         KPREV = KPREV - 1
C
C        if back all the way, must be complex#
C
         IF ( KPREV .LT. 1 ) THEN
            CALL XRCMPX( K, MODE, XRERR )
         ELSE
            GOTO 1
         ENDIF
      ELSEIF ( PRVLEV .LT. DESLEV ) THEN
C
C        if prev opr is lt than desired level, the comma is complex#
C
         CALL XRCMPX( K, MODE, XRERR )
      ENDIF
C
      RETURN
      END
      SUBROUTINE XRCOS( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      DOUBLE PRECISION RADS
      INTEGER K, L
C
C     take the cosine of the table value for index k,
C     store the result in index k-1,
C     and change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there are two possible cosine operations:
C                         cosine(double complex)
C                         cosine(double precision)
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( L - K .NE. 1  .OR. ( VARTYP(K) .NE. 'DP' .AND. VARTYP(K)
     &       .NE. 'DC' ) ) THEN
         CALL DSPERR('XRCOS>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
      IF ( VARTYP(K) .EQ. 'DP' ) THEN
         CALL DG2RAD( RVAREQ(K), RADS )
         RADS = COS( RADS )
         RVAREQ(K-1) = RADS
         VARTYP(K-1)   = 'DP'
      ELSE
         CVAREQ(K-1) = COS( CVAREQ(K) )
         VARTYP(K-1)   = 'DC'
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRDIV( K, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      INTEGER K
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
C     divide the table value for index k-1 by that of index k,
C     store the result in index k-1,
C     and if needed change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there are four possible divisions:
C                         double complex / double complex
C                         double complex / double precision
C                         double precision / double precision
C                         double precision / double complex
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( VARTYP(K-1) .EQ. 'DP' .AND. VARTYP(K) .EQ. 'DP' ) THEN
         IF ( RVAREQ(K) .EQ. ZERO ) GOTO 999
         RVAREQ(K-1) = RVAREQ(K-1) / RVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DC' .AND. VARTYP(K) .EQ. 'DC' ) THEN
         IF (DBLE(CVAREQ(K)).EQ.ZERO.AND.DIMAG(CVAREQ(K)).EQ.ZERO )
     &       GOTO 999
         CVAREQ(K-1) = CVAREQ(K-1) / CVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DC' .AND. VARTYP(K) .EQ. 'DP' ) THEN
         IF ( RVAREQ(K) .EQ. ZERO ) GOTO 999
         CVAREQ(K-1) = CVAREQ(K-1) / RVAREQ(K)
      ELSE
         IF (DBLE(CVAREQ(K)).EQ.ZERO.AND.DIMAG(CVAREQ(K)).EQ.ZERO )
     &       GOTO 999
         CVAREQ(K-1) = RVAREQ(K-1) / CVAREQ(K)
         VARTYP(K-1)   = 'DC'
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
C
C     if divide by zero, set error, print message and return
C
999   XRERR = .TRUE.
      CALL DSPERR('XRDIV>','Divide by zero attempted')
      RETURN
      END
      SUBROUTINE XREQU( K, MODE, XRERR )
      IMPLICIT NONE
C
C     perform the assignment operation '='
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      INTEGER K
      CHARACTER*4 MODE
      LOGICAL XRERR
C
C     possible assignments
C                 dp = dp
C                 dc = dc
C                 dc = dp
C                 st = st
C                 lo = lo
C                 vt = dp or dc or st
C
C
      IF ( K .NE. 2 ) THEN
         XRERR = .TRUE.
         CALL DSPERR('XREQU>','Only one term allowed to left of ''=''')
         RETURN
      ENDIF
C
      IF ( VARTYP(K-1) .EQ. 'DP' .AND. VARTYP(K) .EQ. 'DP' ) THEN
         RVAREQ(K-1) = RVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DC' .AND. VARTYP(K) .EQ. 'DC' ) THEN
         CVAREQ(K-1) = CVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DC' .AND. VARTYP(K) .EQ. 'DP' ) THEN
         IF ( MODE .EQ. 'AMPL' .OR. MODE .EQ. 'SCAL' .OR.
     &        MODE .EQ. 'PHAS' ) THEN
            RVAREQ(K-1) = RVAREQ(K)
            VARTYP(K-1)   = 'DP'
         ELSE
C           for complex use the double prec. number as the amplitude
            CVAREQ(K-1) = RVAREQ(K)
         ENDIF
      ELSEIF ( VARTYP(K-1) .EQ. 'DP' .AND. VARTYP(K) .EQ. 'DC' ) THEN
         XRERR = .TRUE.
         CALL DSPERR('XREQU>',
     &'Attemped to assign double complex to double precision variable')
      ELSEIF ( VARTYP(K-1) .EQ. 'ST' .AND. VARTYP(K) .EQ. 'ST' ) THEN
         SVAREQ(K-1) = SVAREQ(K)
         SVARLN(K-1) = SVARLN(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'LO' .AND. VARTYP(K) .EQ. 'LO' ) THEN
         SVAREQ(K-1) = SVAREQ(K)
         SVARLN(K-1) = SVARLN(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'VT' ) THEN
         IF ( VARTYP(K) .EQ. 'DP' ) THEN
            RVAREQ(K-1) = RVAREQ(K)
         ELSEIF ( VARTYP(K) .EQ. 'DC' ) THEN
            CVAREQ(K-1) = CVAREQ(K)
         ELSEIF ( VARTYP(K) .EQ. 'ST' ) THEN
            SVAREQ(K-1) = SVAREQ(K)
            SVARLN(K-1) = SVARLN(K)
         ELSEIF ( VARTYP(K) .EQ. 'LO' ) THEN
            SVAREQ(K-1) = SVAREQ(K)
            SVARLN(K-1) = SVARLN(K)
         ENDIF
         VARTYP(K-1) = VARTYP(K)
      ELSE
      XRERR = .TRUE.
      CALL DSPERR('XREQU>',
     &'Variable type mismatch between both sides of equation.')
      ENDIF
C
      RETURN
      END
      SUBROUTINE XREXP( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      INTEGER K, L
C
C     take the absolute value of the table value for index k,
C     store the result in index k-1,
C     and change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there are two possible absolute operations:
C                         abs(double complex)
C                         abs(double precision)
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( L - K .NE. 1 .OR. ( VARTYP(K) .NE. 'DC' .AND. VARTYP(K) .NE.
     &      'DP' ) ) THEN
         CALL DSPERR('XREXP>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
      IF ( VARTYP(K) .EQ. 'DP' ) THEN
         RVAREQ(K-1) = EXP( RVAREQ(K) )
         VARTYP(K-1)   = 'DP'
      ELSE
         CVAREQ(K-1) = EXP( CVAREQ(K) )
         VARTYP(K-1)   = 'DC'
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRGAUS( K, L, INDEX, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'comand.inc'
C
      LOGICAL XRERR
      INTEGER K, L
      INTEGER INDEX
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
C     perform the gaussian function
C     expect one argument; std. dev.
C
C     the number of arguments can be derived from k & l
C     all args must be double precision
C
      L = L - 1
C
C     check it's type and while here make sure enough args given
C
      IF ( VARTYP(K) .NE. 'DP' .OR. L-K+1 .NE. 1 ) THEN
         CALL DSPERR('XRGAUS>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
C
      CALL GAUSS( ZERO, RVAREQ(K), RVAREQ(K-1), 1 )
C
C     remove processed entry from table
C
C*      CALL TBPACK( L )
C*      L = L - 1
C
C     must change the vartyp from 'FX' to 'DP' to reflect the
C     result of the function
C
      VARTYP(K-1) = 'DP'
      RETURN
      END
      SUBROUTINE XRRAND( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'seed.inc'
C
      LOGICAL XRERR
      INTEGER K, L
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
C     perform the random function
C     expect zero arguments
C
C     the number of arguments can be derived from k & l
C
      IF ( L - K .NE. 1 ) THEN
         CALL DSPERR('XRRAND>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
C
      CALL GGUBFS( RVAREQ(K-1) )
C
C     remove processed entry from table
C
C*      CALL TBPACK( L )
C*      L = L - 1
C
C     must change the vartyp from 'FX' to 'DP' to reflect the
C     result of the function
C
      VARTYP(K-1) = 'DP'
      RETURN
      END
      SUBROUTINE XRHEAV( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      INTEGER K, L
C
C     perform step function on the  the table value for index k,
C     store the result in index k-1, modify vareq to show it is a result,
C     and change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     perform step function only on double precision numbers.
C     if the variable value is less than or equal to zero return a zero,
C     if greater than zero return a one.
C
      IF ( L - K .NE. 1 .OR. VARTYP(K) .NE. 'DP' ) THEN
         CALL DSPERR('XRHEAV>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
      IF ( RVAREQ(K) .GT. 0.0D0 ) THEN
         RVAREQ(K-1) = 1.0D0
      ELSE
         RVAREQ(K-1) = 0.0D0
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      VARTYP(K-1) = 'DP'
      RETURN
      END
      SUBROUTINE XRINT( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      INTEGER K, L
      LOGICAL XRERR
      INTEGER TMPINT
C
C     perform truncation on the given argument.
C     store the result in index k-1,
C     the calling routine takes care of eliminating the values in
C     index k
C
C     we can:
C            int( double complex ) yielding an integer, but since we
C            don't have integer types then store it in a real
C            int( double precision ) yielding a double precision
C
      IF ( ( VARTYP(K) .NE. 'DP' .AND. VARTYP(K) .NE. 'DC' ) .OR.
     &    L - K .NE. 1 ) THEN
         CALL DSPERR('XRINT>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ELSEIF ( VARTYP(K) .EQ. 'DP' ) THEN
         TMPINT = INT( RVAREQ(K) )
         RVAREQ(K-1) = TMPINT
      ELSEIF ( VARTYP(K) .EQ. 'DC' ) THEN
         TMPINT = INT( DBLE(CVAREQ(K)) )
         RVAREQ(K-1) = TMPINT
      ENDIF
      VARTYP(K-1) = 'DP'
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRNINT( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      INTEGER K, L
      LOGICAL XRERR
      INTEGER TMPINT
C
C     perform truncation on the given argument.
C     store the result in index k-1,
C     the calling routine takes care of eliminating the values in
C     index k
C
C     we can:
C            int( double complex ) yielding an integer, but since we
C            don't have integer types then store it in a real
C            int( double precision ) yielding a double precision
C
      IF ( ( VARTYP(K) .NE. 'DP' .AND. VARTYP(K) .NE. 'DC' ) .OR.
     &    L - K .NE. 1 ) THEN
         CALL DSPERR('XRNINT>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ELSEIF ( VARTYP(K) .EQ. 'DP' ) THEN
         TMPINT = NINT( RVAREQ(K) )
         RVAREQ(K-1) = TMPINT
      ELSEIF ( VARTYP(K) .EQ. 'DC' ) THEN
         TMPINT = NINT( DBLE(CVAREQ(K)) )
         RVAREQ(K-1) = TMPINT
      ENDIF
      VARTYP(K-1) = 'DP'
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRLOG( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'consta.inc'
C
      LOGICAL XRERR
      INTEGER K, L
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
C     take the natural log of the table value for index k,
C     store the result in index k-1,
C     and change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there are two possible absolute operations:
C                         log(double complex)
C                         log(double precision)
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( L - K .NE. 1 .OR. ( VARTYP(K) .NE. 'DP' .AND. VARTYP(K)
     &      .NE. 'DC' ) ) THEN
         CALL DSPERR('XRLOG>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
      IF ( VARTYP(K) .EQ. 'DP' ) THEN
         IF ( RVAREQ(K) .LT. RSMALL ) THEN
            CALL DSPERR('XRLOG>','Illegal argument.')
            XRERR = .TRUE.
         ELSE
            RVAREQ(K-1) = LOG( RVAREQ(K) )
            VARTYP(K-1)   = 'DP'
         ENDIF
      ELSE
       IF (DBLE(CVAREQ(K)).EQ.ZERO.AND.DIMAG(CVAREQ(K)).EQ.ZERO ) THEN
            CALL DSPERR('XRLOG>','Illegal argument. LOG((0.,0.))')
            XRERR = .TRUE.
         ELSE
            CVAREQ(K-1) = LOG( CVAREQ(K) )
            VARTYP(K-1)   = 'DC'
         ENDIF
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRMAX( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      INTEGER K, L
      INTEGER I
C
C     perform the max function
C
C     the number of arguments can be derived from k & l
C     all args must be double precision
C
      L = L - 1
C
C     first argument is the min until further checking is done
C     check it's type and while here make sure enough args given
C
      IF ( VARTYP(K) .NE. 'DP' .OR. L-K+1 .LT. 2 ) THEN
         CALL DSPERR('XRMAX>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
      RVAREQ(K-1) = RVAREQ(K)
      DO I = K+1,L
C
C        check arg type
C
         IF ( VARTYP(L) .NE. 'DP'  .OR. TEO(L) .NE. '_' ) THEN
            CALL DSPERR('XRMAX>','Illegal argument(s)')
            XRERR = .TRUE.
            RETURN
         ENDIF
C
C        compare vs. established min
C
         IF ( RVAREQ(L) .GT. RVAREQ(K-1) ) RVAREQ(K-1) = RVAREQ(L)
C
C        remove processed entry from table
C
         CALL TBPACK( L )
         L = L - 1
      ENDDO
C
C     must change the vartyp from 'FX' to 'DP' to reflect the
C     result of the function
C
      VARTYP(K-1) = 'DP'
      RETURN
      END
      SUBROUTINE XRMAXW( K, L, INDEX, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
C
      LOGICAL XRERR
      INTEGER K, L, INDEX
      DOUBLE PRECISION ZERO, ARGIN
      PARAMETER (ZERO=0.0D0)
C
C     perform the maxwell function
C     expect one argument; temp
C
C     the number of arguments can be derived from k & l
C     all args must be double precision
C
      L = L - 1
C
C     check it's type and while here make sure enough args given
C
      IF ( VARTYP(K) .NE. 'DP' .OR. L-K+1 .NE. 1 ) THEN
         CALL DSPERR('XRMAXW>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
C
      IF ( ABS(AMASS(INDEX)) .LT. RSMALL ) THEN
         CALL DSPERR('XRMAXW>','Zero mass encountered')
         XRERR = .TRUE.
         RETURN
      ENDIF
C
      ARGIN = SQRT(RVAREQ(K)*KBOLTZ/AMASS(INDEX))
      CALL GAUSS(ZERO,ARGIN,RVAREQ(K-1),1)
      RVAREQ(K-1) = RVAREQ(K-1) / TIMFAC
C
C     remove processed entry from table
C
C*      CALL TBPACK( L )
C*      L = L - 1
C
C     must change the vartyp from 'FX' to 'DP' to reflect the
C     result of the function
C
      VARTYP(K-1) = 'DP'
      RETURN
      END
      SUBROUTINE XRMIN( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      INTEGER K, L
      INTEGER I
C
C     perform the min function
C
C     the number of arguments can be derived from k & l
C     all args must be double precision
C
      L = L - 1
C
C     first argument is the min until further checking is done
C     check it's type and while here make sure enough args given
C
      IF ( VARTYP(K) .NE. 'DP' .OR. L-K+1 .LT. 2 ) THEN
         CALL DSPERR('XRMIN>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
      RVAREQ(K-1) = RVAREQ(K)
      DO I = K+1,L
C
C        check arg type
C
         IF ( VARTYP(L) .NE. 'DP'  .OR. TEO(L) .NE. '_' ) THEN
            CALL DSPERR('XRMIN>','Illegal argument(s)')
            XRERR = .TRUE.
            RETURN
         ENDIF
C
C        compare vs. established min
C
         IF ( RVAREQ(L) .LT. RVAREQ(K-1) ) RVAREQ(K-1) = RVAREQ(L)
C
C        remove processed entry from table
C
         CALL TBPACK( L )
         L = L - 1
      ENDDO
C
C     must change the vartyp from 'FX' to 'DP' to reflect the
C     result of the function
C
      VARTYP(K-1) = 'DP'
      RETURN
      END
      SUBROUTINE XRMOD( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      INTEGER K, L
C
C     perform the mod function
C
C     the number of arguments can be derived from k & l
C     all args must be double precision
C
      L = L - 1
C
C     return the remainder when the first argument is divided by
C     the second argument
C
      IF ( VARTYP(K) .NE. 'DP' .OR. VARTYP(K+1) .NE. 'DP' .OR.
     &                                              L-K+1 .NE. 2 ) THEN
         CALL DSPERR('XRMOD>','Illegal argument type or count.')
         XRERR = .TRUE.
         RETURN
      ENDIF
C
      RVAREQ(K-1) = MOD(RVAREQ(K),RVAREQ(K+1))
C
      CALL TBPACK( L )
      L = L - 1
C
C     must change the vartyp from 'FX' to 'DP' to reflect the
C     result of the function
C
      VARTYP(K-1) = 'DP'
      RETURN
      END
      SUBROUTINE XRIMOD( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      INTEGER K, L
C
C     perform the imod function
C
C     the number of arguments can be derived from k & l
C     all args must be double precision
C
      L = L - 1
C
C     return the remainder when the first argument is divided by
C     the second argument
C
      IF ( VARTYP(K) .NE. 'DP' .OR. VARTYP(K+1) .NE. 'DP' .OR.
     &                                              L-K+1 .NE. 2 ) THEN
         CALL DSPERR('XRMOD>','Illegal argument type or count.')
         XRERR = .TRUE.
         RETURN
      ENDIF
C
      RVAREQ(K-1) = MOD(NINT(RVAREQ(K)),NINT(RVAREQ(K+1)))
C
      CALL TBPACK( L )
      L = L - 1
C
C     must change the vartyp from 'FX' to 'DP' to reflect the
C     result of the function
C
      VARTYP(K-1) = 'DP'
      RETURN
      END
      SUBROUTINE XRMULT( K, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      INTEGER K
C local
      INTEGER I,LN,LMAX
C begin
C
C     multiply the table value for index k with that of index k-1,
C     store the result in index k-1, modify vareq to show it is a result,
C     and if needed change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there are four possible multiplications:
C                         double complex * double complex
C                         double complex * double precision
C                         double precision * double precision
C                         double precision * double complex
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( VARTYP(K-1) .EQ. 'DP' .AND. VARTYP(K) .EQ. 'DP' ) THEN
         RVAREQ(K-1) = RVAREQ(K-1) * RVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DC' .AND. VARTYP(K) .EQ. 'DC' ) THEN
         CVAREQ(K-1) = CVAREQ(K-1) * CVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DC' .AND. VARTYP(K) .EQ. 'DP' ) THEN
         CVAREQ(K-1) = CVAREQ(K-1) * RVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DP' .AND. VARTYP(K) .EQ. 'DC' ) THEN
         CVAREQ(K-1) = RVAREQ(K-1) * CVAREQ(K)
         VARTYP(K-1)   = 'DC'
      ELSEIF ( VARTYP(K-1) .EQ. 'ST' .AND. VARTYP(K) .EQ. 'DP' ) THEN
         LMAX=LEN(SVAREQ(K-1))
         LN=SVARLN(K-1)
         SVARLN(K-1)=MIN(LN*INT(RVAREQ(K)),LMAX)
         DO I=LN+1,SVARLN(K-1),LN
            SVAREQ(K-1)(I:MIN(I+LN-1,LMAX)) = SVAREQ(K-1)(1:LN)
         END DO
      ELSE
         CALL DSPERR('XRMULT>','Illegal data types')
         XRERR = .TRUE.
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
C
      SUBROUTINE XRPOWR( K, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      INTEGER K
      LOGICAL XRERR
C
C     raise the table value for index k-1 to the power of index k,
C     store the result in index k-1, modify vareq to show it is a result,
C     and if needed change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there are
C                         double complex ^ double complex
C                         double complex ^ double precision
C                         double precision ^ double precision
C                         double precision ^ double complex
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( VARTYP(K-1) .EQ. 'DP' .AND. VARTYP(K) .EQ. 'DP' ) THEN
         IF ( RVAREQ(K-1) .LT. 0.0D0 ) THEN
            IF ( RVAREQ(K) .EQ. FLOAT(INT(RVAREQ(K)))) THEN
               RVAREQ(K-1) = RVAREQ(K-1)**(INT(RVAREQ(K)))
            ELSE
               CALL DSPERR('XRPOWR>',
     &           'Attempt to raise a negative number to a real power.')
               XRERR = .TRUE.
            ENDIF
         ELSEIF ( RVAREQ(K-1) .EQ. 0.D0.AND.RVAREQ(K) .EQ. 0.D0) THEN
            RVAREQ(K-1) = 1.0D0
         ELSE
            RVAREQ(K-1) = RVAREQ(K-1)**RVAREQ(K)
         ENDIF
      ELSEIF ( VARTYP(K-1) .EQ. 'DC' .AND. VARTYP(K) .EQ. 'DC' ) THEN
         CVAREQ(K-1) = CVAREQ(K-1)**CVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DC' .AND. VARTYP(K) .EQ. 'DP' ) THEN
         CVAREQ(K-1) = CVAREQ(K-1)**RVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DP' .AND. VARTYP(K) .EQ. 'DC' ) THEN
         CVAREQ(K-1) = RVAREQ(K-1)**CVAREQ(K)
         VARTYP(K-1)   = 'DC'
      ELSE
         CALL DSPERR('XRPOWR','Illegal combination of data types')
         XRERR = .TRUE.
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRSIGN( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      INTEGER K, L
C
C     perform transfer of sign function on the  the table value index k,
C     store the result in index k-1,
C     and change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     perform sign transfer only on double precision numbers.
C     if the variable value is less than zero return a negative one,
C     otherwise return a positive one.
C
      IF ( L - K .NE. 1 .OR. VARTYP(K) .NE. 'DP' ) THEN
         CALL DSPERR('XRSIGN>','Illegal argument(s).')
         XRERR = .TRUE.
         RETURN
      ENDIF
      IF ( RVAREQ(K) .GE. 0.0D0 ) THEN
         RVAREQ(K-1) = 1.0D0
      ELSE
         RVAREQ(K-1) = -1.0D0
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      VARTYP(K-1) = 'DP'
      RETURN
      END
      SUBROUTINE XRSINE( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      DOUBLE PRECISION RADS
      INTEGER K, L
C
C     take the sine of the table value for index k,
C     store the result in index k-1,
C     and change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there are two possible sine operations:
C                         sine(double complex)
C                         sine(double precision)
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( L - K .NE. 1 .OR. ( VARTYP(K) .NE. 'DC' .AND. VARTYP(K)
     &     .NE. 'DP' ) ) THEN
         CALL DSPERR('XRSINE>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
      IF ( VARTYP(K) .EQ. 'DP' ) THEN
         CALL DG2RAD( RVAREQ(K), RADS )
         RADS = SIN( RADS )
         RVAREQ(K-1) = RADS
         VARTYP(K-1)   = 'DP'
      ELSE
         CVAREQ(K-1) = SIN( CVAREQ(K) )
         VARTYP(K-1)   = 'DC'
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRSQRT( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      INTEGER K, L
      LOGICAL XRERR
C
C     take the square root of the table value for index k,
C     store the result in index k-1,
C     and if needed change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there are:
C                   sqrt( double complex )
C                   sqrt( double precision )
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( ( VARTYP(K) .NE. 'DP' .AND. VARTYP(K) .NE. 'DC' ) .OR.
     &    L - K .NE. 1 ) THEN
         CALL DSPERR('XRSQRT>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ELSEIF ( VARTYP(K) .EQ. 'DP' ) THEN
         IF ( RVAREQ(K) .LT. 0.0D0 ) THEN
            CALL DSPERR('XRSQRT>','Illegal negative argument')
            XRERR = .TRUE.
            RETURN
         ENDIF
         RVAREQ(K-1) = SQRT( RVAREQ(K) )
      ELSEIF ( VARTYP(K) .EQ. 'DC' ) THEN
         CVAREQ(K-1) = SQRT( CVAREQ(K) )
      ENDIF
      VARTYP(K-1) = VARTYP(K)
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRSUB( K, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'symbol.inc'
C
C
      INTEGER K
      INTEGER I, LAST1, LAST2, IND
      LOGICAL XRERR
      CHARACTER*1 CTEMP
C
C     subtract the table value for index k from that of index k-1,
C     store the result in index k-1, modify vareq to show it is a result,
C     and if needed change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there are five possible subtractions:
C                         double complex - double complex
C                         double complex - double precision
C                         double precision - double precision
C                         double precision - double complex
C                         string - string ( remove substring )
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( VARTYP(K-1) .EQ. 'DP' .AND. VARTYP(K) .EQ. 'DP' ) THEN
         RVAREQ(K-1) = RVAREQ(K-1) - RVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DP' .AND. VARTYP(K) .EQ. 'DC' ) THEN
         CVAREQ(K-1) = RVAREQ(K-1) - CVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DC' .AND. VARTYP(K) .EQ. 'DC' ) THEN
         CVAREQ(K-1) = CVAREQ(K-1) - CVAREQ(K)
      ELSEIF ( VARTYP(K-1) .EQ. 'DC' .AND. VARTYP(K) .EQ. 'DP' ) THEN
         CVAREQ(K-1) = CVAREQ(K-1) - RVAREQ(K)
         VARTYP(K-1)   = 'DC'
      ELSEIF ( VARTYP(K-1) .EQ. 'ST' .AND. VARTYP(K) .EQ. 'ST' ) THEN
C
C        remove substring if it is present
C
         LAST1 = SVARLN(K-1)
         LAST2 = SVARLN(K)
C
C        if substring is longer than string it can't be there
C
         IF ( LAST1 .GE. LAST2 ) THEN
            IND = INDEX(SVAREQ(K-1)(1:LAST1),SVAREQ(K)(1:LAST2))
            IF ( IND .GT. 0 ) THEN
               DO I = 1,XSTMAX
                  IF ( I .GE. IND ) THEN
                    IF ( I + LAST2 .GT. XSTMAX ) THEN
                       SVAREQ(K-1)(I:I) = ' '
                    ELSE
                      CTEMP=SVAREQ(K-1)(I+LAST2:I+LAST2)
                      SVAREQ(K-1)(I:I)=CTEMP
                    ENDIF
                  ENDIF
               ENDDO
               SVARLN(K-1) = SVARLN(K-1) - SVARLN(K)
            ENDIF
         ENDIF
      ELSE
         CALL DSPERR('XRSUB','Illegal combination of data types')
         XRERR = .TRUE.
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRTAN( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      DOUBLE PRECISION RADS
      INTEGER K, L
C
C     take the tangent of the table value for index k,
C     store the result in index k-1,
C     and change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there is one possible tangent operation:
C                         tan(double precision)
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( L - K .NE. 1 .OR. VARTYP(K) .NE. 'DP' ) THEN
         CALL DSPERR('XRTAN>','Illegal argument(s).')
         XRERR = .TRUE.
         RETURN
      ENDIF
      CALL DG2RAD( RVAREQ(K), RADS )
      RADS = TAN( RADS )
      RVAREQ(K-1) = RADS
      VARTYP(K-1)   = 'DP'
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRUNM( K, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      INTEGER K
C
C     apply the unary minus sign to the  the table value for index k,
C     store the result in index k-1, modify vareq to show it is a result,
C     and change the vartype for k-1.  the calling routine
C     takes care of eliminating the values in index k
C
C     there are two possible absolute operations:
C                         ~(double complex)
C                         ~(double precision)
C     let the compiler make the numeric conversion, just make
C     sure to store the result in the correct type of variable
C
      IF ( VARTYP(K-1) .NE. '  ' ) THEN
          WRITE(6,'(A)')' %CNS-XRUNM-ERR: Error applying unary minus'
          XRERR = .TRUE.
          RETURN
      ENDIF
      IF ( VARTYP(K) .EQ. 'DP' ) THEN
         RVAREQ(K-1) = -1 * RVAREQ(K)
         VARTYP(K-1)   = 'DP'
      ELSEIF ( VARTYP(K) .EQ. 'DC' ) THEN
         CVAREQ(K-1) = -1 * CVAREQ(K)
         VARTYP(K-1)   = 'DC'
      ELSE
         CALL DSPERR('XRUNM>','Error applying unary minus')
         XRERR = .TRUE.
      ENDIF
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE DG2RAD( DEG, RAD )
      IMPLICIT NONE
C
      DOUBLE PRECISION DEG, RAD
C
      INCLUDE 'consta.inc'
C
      RAD = PI * DEG / 180.0
C
      RETURN
      END
      SUBROUTINE RD2DEG( RAD, DEG )
      IMPLICIT NONE
C
      DOUBLE PRECISION DEG, RAD
C
      INCLUDE 'consta.inc'
C
      DEG = RAD * 180.0 / PI
C
      RETURN
      END
      SUBROUTINE XRENCD( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      INTEGER K, L
      INTEGER RETLEN
C
C     equiv of fortran encode, turn double precision into a string
C     store the result in index k-1,
C     and change the vartype for k-1.
C
      IF ( L - K .NE. 1 .OR. VARTYP(K) .NE. 'DP' ) THEN
         CALL DSPERR('XRENCODE>','Illegal argument(s).')
         XRERR = .TRUE.
         RETURN
      ENDIF
      CALL ENCODF( RVAREQ(K), SVAREQ(K-1), XSTMAX, RETLEN )
      VARTYP(K-1)   = 'ST'
      SVARLN(K-1)   = RETLEN
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
      SUBROUTINE XRDECD( K, L, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      DOUBLE PRECISION DECODF
      LOGICAL OK
      LOGICAL XRERR
      INTEGER K, L
      INTEGER I, ERRMSG_LEN
      CHARACTER*1 ERRMSG
C
C     equiv of fortran encode, turn double precision into a string
C     store the result in index k-1,
C     and change the vartype for k-1.
C
      IF ( L - K .NE. 1 .OR. VARTYP(K) .NE. 'ST' ) THEN
         CALL DSPERR('XRDECODE>','Illegal argument(s).')
         XRERR = .TRUE.
         RETURN
      ENDIF
C
C MODIFICATION: Use hybrid-36 decoding for 4-character strings that
C begin with a letter.
      IF (SVARLN(K).EQ.4 .AND. LGE(SVAREQ(K)(1:1),'A')) THEN
      CALL HY36DECODE(4,SVAREQ(K)(1:4),I,ERRMSG,ERRMSG_LEN)
      OK = ERRMSG_LEN.EQ.0
      IF (OK) RVAREQ(K-1) = I
      IF (.NOT.OK) CALL DSPERR('XRDECODE>','Invalid Hybrid-36 string')
      ELSE
      RVAREQ(K-1) = DECODF( SVAREQ(K), SVARLN(K), OK )
      ENDIF
C END MODIFICATION
      IF ( .NOT. OK ) THEN
         CALL DSPERR('XRDECODE>','Illegal argument(s).')
         XRERR = .TRUE.
         RETURN
      ENDIF
      VARTYP(K-1)   = 'DP'
C
Ctabtrace      VAREQ(K-1) = 'RESULT              '
      RETURN
      END
C============================================
      SUBROUTINE XRCAPIT( K, L, XRERR )
      IMPLICIT NONE
C
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'comand.inc'
      LOGICAL XRERR
      INTEGER K, L
C local
      INTEGER I
      CHARACTER*(1) SCHR
      INTEGER ICHR
C begin
C
      IF ( L - K .NE. 1 .OR. VARTYP(K) .NE. 'ST' ) THEN
         CALL DSPERR('XRCAPIT>','Illegal argument(s).')
         XRERR = .TRUE.
         RETURN
      ENDIF
C
      DO I=1,SVARLN(K)
      SCHR=SVAREQ(K)(I:I)
      ICHR=ASCIIM(ICHAR(SCHR))
      SVAREQ(K-1)(I:I)=CHAR(ICHR)
      END DO
C
      VARTYP(K-1)   = 'ST'
      SVARLN(K-1)   = SVARLN(K)
      RETURN
      END
C==================================================
      SUBROUTINE XRLEN( K, L, XRERR )
      IMPLICIT NONE
C
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'comand.inc'
      LOGICAL XRERR
      INTEGER K, L
C begin
C
      IF ( L - K .NE. 1 .OR. VARTYP(K) .NE. 'ST' ) THEN
         CALL DSPERR('XRLEN>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
C
      RVAREQ(K-1) = SVARLN(K)
      VARTYP(K-1) = 'DP'
      RETURN
      END
C==================================================
      SUBROUTINE XRSUBSTR( K, L, XRERR )
      IMPLICIT NONE
C
C substring operation (similar to FORTRAN)
C negative indices are allowed:
C   substr($string,$start,-1) means to take from $start to the end, 
C   substr($string,-2,-1) takes the last two chars.
C
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'comand.inc'
      LOGICAL XRERR
      INTEGER K, L
C local
      INTEGER ISTART, IEND
C begin
C
      IF ( L - K .NE. 3 .OR. VARTYP(K) .NE. 'ST' .OR.
     &     VARTYP(K+1) .NE. 'DP' .OR. VARTYP(K+2) .NE. 'DP') THEN
         CALL DSPERR('XRSUBSTR>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
C
      ISTART = INT(RVAREQ(K+1))
      IF (ISTART.LT.0) ISTART = SVARLN(K)+1+ISTART
      ISTART = MIN(MAX(ISTART,1),SVARLN(K)+1)
      IEND = INT(RVAREQ(K+2))
      IF (IEND.LT.0) IEND = SVARLN(K)+1+IEND
      IEND = MIN(MAX(IEND,ISTART-1),SVARLN(K))
      SVAREQ(K-1) = SVAREQ(K)(ISTART:IEND)
      SVARLN(K-1) = IEND-ISTART+1
      VARTYP(K-1) = 'ST'
      CALL TBPACK(L-1)
      CALL TBPACK(L-2)
C
      RETURN
      END
C==================================================
      SUBROUTINE XRINDEX( K, L, XRERR )
      IMPLICIT NONE
C
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'comand.inc'
      LOGICAL XRERR
      INTEGER K, L
C local
      INTEGER ISTART, IEND
C begin
C
      IF ( L - K .NE. 2 .OR. VARTYP(K) .NE. 'ST' .OR.
     &     VARTYP(K+1) .NE. 'ST') THEN
         CALL DSPERR('XRINDEX>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
C
      RVAREQ(K-1) = INDEX( SVAREQ(K)(1:SVARLN(K)),
     &                SVAREQ(K+1)(1:SVARLN(K+1)) )
      VARTYP(K-1) = 'DP'
      CALL TBPACK(L-1)
      RETURN
      END
C==================================================
      SUBROUTINE XRRINDEX( K, L, XRERR )
      IMPLICIT NONE
C
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      LOGICAL XRERR
      INTEGER K, L
C local
      INTEGER I, ISTART, IEND
C begin
C
      IF ( L - K .NE. 2 .OR. VARTYP(K) .NE. 'ST' .OR.
     &     VARTYP(K+1) .NE. 'ST') THEN
         CALL DSPERR('XRRINDEX>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
C
      RVAREQ(K-1) = RINDEX( SVAREQ(K)(1:SVARLN(K)),
     &                SVAREQ(K+1)(1:SVARLN(K+1)) )
      VARTYP(K-1) = 'DP'
      CALL TBPACK(L-1)
      RETURN
      END
C==================================================
      SUBROUTINE XRFORMAT( K, L, XRERR )
      IMPLICIT NONE
C
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'comand.inc'
      LOGICAL XRERR
      INTEGER K, L
C local
      LOGICAL OK
      CHARACTER*(SVARMX+6) CFORM
C begin
C
      IF ( L - K .NE. 2 .OR. VARTYP(K) .NE. 'ST') THEN
         CALL DSPERR('XRFORMAT>','Illegal argument')
         XRERR = .TRUE.
         RETURN
      ENDIF
C=============================================================
      CALL STRFMT(SVAREQ(K-1),LEN(SVAREQ(K-1)),SVARLN(K-1),OK,
     &            VARTYP(K+1),CVAREQ(K+1),RVAREQ(K+1),SVAREQ(K+1),
     &            SVARLN(K+1),SVAREQ(K),SVARLN(K),CFORM)
      IF(.NOT.OK) XRERR=.TRUE.
      VARTYP(K-1) = 'ST'
      CALL TBPACK(L-1)
      RETURN
      END
C==================================================
      SUBROUTINE XRADJUSTL( K, L, XRERR )
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'comand.inc'
      LOGICAL XRERR
      INTEGER K, L
C begin
      IF ( L - K .NE. 1 .OR. VARTYP(K) .NE. 'ST' ) THEN
         CALL DSPERR('XRADJUSTL>','Illegal argument')
         XRERR = .TRUE.
         RETURN
      ENDIF
      SVAREQ(K-1)=SVAREQ(K)
      SVARLN(K-1)=SVARLN(K)
C     Trim left, but keep the same length.
      CALL TRIML(SVAREQ(K-1),SVARLN(K))
      VARTYP(K-1) = 'ST'
      RETURN
      END
C==================================================
      SUBROUTINE XRTRIM( K, L, XRERR )
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'comand.inc'
      LOGICAL XRERR
      INTEGER K, L
C begin
      IF ( L - K .NE. 1 .OR. VARTYP(K) .NE. 'ST' ) THEN
         CALL DSPERR('XRTRIM>','Illegal argument')
         XRERR = .TRUE.
         RETURN
      ENDIF
      SVARLN(K-1)=SVARLN(K)
      SVAREQ(K-1)(1:SVARLN(K-1))=SVAREQ(K)(1:SVARLN(K))
      CALL TRIMM(SVAREQ(K-1),SVARLN(K-1))
      VARTYP(K-1) = 'ST'
      RETURN
      END
C==================================================
      SUBROUTINE XRMATCH( K, L, XRERR )
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'comand.inc'
      LOGICAL XRERR
      INTEGER K, L
C local
      LOGICAL Q(1)
C begin
      IF ( L - K .NE. 2 .OR. VARTYP(K) .NE. 'ST' .OR.
     &     VARTYP(K+1) .NE. 'ST') THEN
         CALL DSPERR('XRMATCH>','Illegal argument(s)')
         XRERR = .TRUE.
         RETURN
      ENDIF
      CALL EQSTWC(SVAREQ(K),SVARLN(K),SVAREQ(K+1),SVARLN(K+1),1,1,Q)
      IF (Q(1)) THEN
      SVAREQ(K-1)='TRUE'
      SVARLN(K-1)=4
      ELSE
      SVAREQ(K-1)='FALSE'
      SVARLN(K-1)=5
      END IF
      VARTYP(K-1) = 'LO'
      CALL TBPACK(L-1)
      RETURN
      END
C==================================================
