C{
      subroutine FreeStarTr(StarTrN, StarTrPi, StarTrPs, nStar)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      integer  StarTrN(*), StarTrPi(*), StarTrPs(*)
      integer  nStar
C
C Free up dynamic memory allocated for STAR integer translate
C
C Author: R.W. Grosse-Kunstleve, 1998
C
C local
      integer  i
C
C begin
      do i = 1, nStar
        if (StarTrN(i) .gt. 0) then
          call  FreHP(StarTrPi(i), integ4(StarTrN(i)))
          call cFreHP(StarTrPs(i), StarTrN(i), WORD_SIZE)
          StarTrN(i) = 0
          StarTrPi(i) = 0
          StarTrPs(i) = 0
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      SUBROUTINE XWRSTARTRP(STARTRN, STARTRPI, STARTRPS)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INTEGER  STARTRN
      INTEGER  STARTRPI
      INTEGER  STARTRPS
C
C Parser for WRITE REFLECTION STAR TRANSLATE
C
C Author: R.W. Grosse-Kunstleve, 1998
C
C local
      character  lprompt*10
      logical    found
      integer    ival
C
C begin
      lprompt = 'translate>'
C parsing
      CALL PUSEND(lprompt)
      DO WHILE (.NOT. DONE)
        CALL NEXTWD(lprompt)
        CALL MISCOM(lprompt, USED)
        IF (.NOT. USED) THEN
          IF      (WD(1:4) .EQ. 'HELP') THEN
CCC            CALL CNSHELP('cns-xray-write-reflection-star-translate')
          ELSE IF (WD(1:4) .EQ. 'END ') THEN
            CALL CHKEND(lprompt, DONE)
          ELSE
            ival = DECODI(WD, WDLEN, FOUND)
            ERROR = .NOT. FOUND
            IF (ERROR) then
              CALL DSPERR(lprompt, 'integer number expected')
            else
              CALL NEXTASS('string=')
              STARTRPI =  reahp(STARTRPI, integ4(STARTRN),
     &                                    integ4(STARTRN+1))
              STARTRPS = creahp(STARTRPS, STARTRN,
     &                                    STARTRN+1, WORD_SIZE)
              STARTRN = STARTRN + 1
              call SetI4Elem(heap(STARTRPI), STARTRN, ival)
              call XCOPYCH(CWSHEAP(STARTRPS), STARTRN, wd(1:wdlen), 1)
            end if
          END IF
        END IF
      END DO
      DONE=.FALSE.
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE XWRSTARP(STARLBL, STAROBJ, STARVAL,
     &                    STARTRN, STARTRPI, STARTRPS,
     &                    NSTAR,
     &                    XSFNUM, XSFNAM, XSFTYPE, HPSF)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      CHARACTER  STARLBL(*)*(*)
      INTEGER    STAROBJ(*)
      CHARACTER  STARVAL(*)*(*)
      INTEGER    STARTRN(*)
      INTEGER    STARTRPI(*)
      INTEGER    STARTRPS(*)
      INTEGER    NSTAR
      INTEGER    XSFNUM
      CHARACTER  XSFNAM(*)*(*), XSFTYPE(*)*(*)
      INTEGER    HPSF(*)
C
C Parser for WRITE REFLECTION STAR
C
C Author: R.W. Grosse-Kunstleve, 1998
C
C local
      LOGICAL  FORGET
C
C begin
      NSTAR = NSTAR + 1
      STARLBL(NSTAR) = ' '
      STAROBJ(NSTAR) = 0
      STARVAL(NSTAR) = ' '
      STARTRN(NSTAR) = 0
      STARTRPI(NSTAR) = 0
      STARTRPS(NSTAR) = 0
C
C parsing
      CALL PUSEND('STAR>')
      DO WHILE (.NOT. DONE)
        CALL NEXTWD('STAR>')
        CALL MISCOM('STAR>', USED)
        IF (.NOT. USED) THEN
          IF      (WD(1:4) .EQ. 'HELP') THEN
CCC            CALL CNSHELP('cns-xray-write-reflection-star')
          ELSE IF (WD(1:4) .EQ. 'LABE') THEN
            CALL NEXTASS('label=')
            STARLBL(NSTAR) = WD(1:WDLEN)
          ELSE IF (WD(1:4) .EQ. 'OBJE') THEN
            CALL GETIOBJ(STAROBJ(NSTAR), XSFNUM, XSFNAM, 'HKL',
     &                   'object=', 'reciprocal', 'STAR')
          ELSE IF (WD(1:4) .EQ. 'VALU') THEN
            CALL NEXTASS('value=')
            IF (      WD(1:4) .NE. 'AMPL'
     &          .AND. WD(1:4) .NE. 'PHAS'
     &          .AND. WD(1:4) .NE. 'REAL'
     &          .AND. WD(1:4) .NE. 'IMAG') THEN
              CALL DSPERR('STAR', 'Unknown parameter.')
            ELSE
              STARVAL(NSTAR) = WD(1:WDLEN)
            END IF
          ELSE IF (WD(1:4) .EQ. 'TRAN') THEN
            CALL XWRSTARTRP(STARTRN(NSTAR),
     &                      STARTRPI(NSTAR), STARTRPS(NSTAR))
          ELSE
            CALL CHKEND('STAR>', DONE)
          END IF
        END IF
      END DO
      DONE=.FALSE.
C
C Some consistency checks
      FORGET = .FALSE.
C
      IF (      STARLBL(NSTAR) .EQ. ' '
     &    .AND. STAROBJ(NSTAR) .EQ. 0
     &    .AND. STARVAL(NSTAR) .EQ. ' '
     &    .AND. STARTRN(NSTAR) .EQ. 0) THEN
        FORGET = .TRUE.
      ELSE IF (STARLBL(NSTAR) .EQ. ' ') THEN
        CALL WRNDIE(-5, 'STAR', 'label not specified.')
        FORGET = .TRUE.
      ELSE IF (STAROBJ(NSTAR) .EQ. 0) THEN
        CALL WRNDIE(-5, 'STAR', 'object not specified.')
        FORGET = .TRUE.
      ELSE
        IF (STAROBJ(NSTAR) .GT. 0) THEN
          IF (HPSF(STAROBJ(NSTAR)) .EQ. 0) THEN
            CALL WRNDIE(-5, 'STAR', 'object not set.')
            FORGET = .TRUE.
          END IF
        END IF
        IF (      STARVAL(NSTAR) .NE. ' '
     &      .AND. XSFTYPE(STAROBJ(NSTAR)) .NE. 'COMP') THEN
          CALL WRNDIE(10, 'STAR',
     &                'value type specified for non-complex object.')
          STARVAL(NSTAR) = ' '
        END IF
        IF (      STARTRN(NSTAR) .NE. 0
     &      .AND. XSFTYPE(STAROBJ(NSTAR)) .NE. 'INTE') THEN
          CALL WRNDIE(10, 'STAR',
     &                'translate only possible for integer objects.')
          STARTRN(NSTAR) = 0
        END IF
      END IF
C
      IF (FORGET) NSTAR = NSTAR - 1
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE DoStarTr(ival, StarTrN, StarTrI, StarTrS, StarBuf)
      IMPLICIT NONE
C I/O
      integer    ival
      integer    StarTrN
      integer    StarTrI(*)
      character  StarTrS(*)*(*)
      character  StarBuf*(*)
C
C STAR TRANSLATE lookup for ival
C
C Author: R.W. Grosse-Kunstleve, 1998
C
C local
      integer  i
C
C begin
      do i = 1, StarTrN
        if (StarTrI(i) .eq. ival) then
          StarBuf = StarTrS(i)
          return
        end if
      end do
C
      write(StarBuf, '(I80)') ival
C
      return
      end
C}
C=======================================================================
C{
      SUBROUTINE XWRSTAR(OUNIT, STARBUF, LSTARBUF,
     &                   STARLBL, STAROBJ, STARVAL,
     &                   STARTRN, STARTRPI, STARTRPS,
     &                   NSTAR,
     &                   XSFTYPE, HPSF,
     &                   NREF, XRH, XRK, XRL, QSELE)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INTEGER    OUNIT
      CHARACTER  STARBUF(*)*(*)
      INTEGER    LSTARBUF(*)
      CHARACTER  STARLBL(*)*(*)
      INTEGER    STAROBJ(*)
      CHARACTER  STARVAL(*)*(*)
      INTEGER    STARTRN(*)
      INTEGER    STARTRPI(*)
      INTEGER    STARTRPS(*)
      INTEGER    NSTAR
      CHARACTER  XSFTYPE(*)*(*)
      INTEGER    HPSF(*)
      INTEGER    NREF
      INTEGER    XRH(*), XRK(*), XRL(*)
      LOGICAL    QSELE(*)
C
C Write STAR reflection file
C
C Author: R.W. Grosse-Kunstleve, 1998
C
C local
      INTEGER           N, I, IREF, SUM, FIRST, LAST, CURR
      INTEGER           I4VAL
      DOUBLE PRECISION  R8VAL, R8TMP
      DOUBLE COMPLEX    C8VAL
C
C externals
      INTEGER   LNBLNK
      EXTERNAL  LNBLNK
C
C begin
      N = 0
      DO IREF = 1, NREF
        IF (QSELE(IREF)) N = N + 1
      END DO
      IF (N .EQ. 0) RETURN
C
      WRITE(OUNIT, '(1X, A)') 'loop_'
      DO I = 1, NSTAR
        WRITE(OUNIT, '(1X, A)') STARLBL(I)(:LNBLNK(STARLBL(I)))
      END DO
C
      DO IREF = 1, NREF
        IF (QSELE(IREF)) THEN
          DO I = 1, NSTAR
            STARBUF(I) = '?'
            IF (STAROBJ(I) .LT. 0) THEN
              IF      (STAROBJ(I) .EQ. -1) THEN
                I4VAL = XRH(IREF)
              ELSE IF (STAROBJ(I) .EQ. -2) THEN
                I4VAL = XRK(IREF)
              ELSE IF (STAROBJ(I) .EQ. -3) THEN
                I4VAL = XRL(IREF)
              ELSE
                CALL WRNDIE(-5, 'XWRSTAR', 'Fatal Coding Error')
                CALL DIE
              END IF
              WRITE(STARBUF(I), '(I80)') I4VAL
            ELSE
              IF      (XSFTYPE(STAROBJ(I)) .EQ. 'INTE') THEN
                CALL XCOPYI(I4VAL, 1, HEAP(HPSF(STAROBJ(I))), IREF)
                CALL DOSTARTR(I4VAL,
     &                        STARTRN(I),
     &                           HEAP(STARTRPI(I)),
     &                        CWSHEAP(STARTRPS(I)),
     &                        STARBUF(I))
              ELSE IF (XSFTYPE(STAROBJ(I)) .EQ. 'REAL') THEN
                CALL XCOPYR(R8VAL, 1, HEAP(HPSF(STAROBJ(I))), IREF)
                IF (R8VAL .EQ. DFLOAT(0)) THEN
                  STARBUF(I) = '0'
                ELSE
                  WRITE(STARBUF(I), '(G80.6)') R8VAL
                END IF
              ELSE IF (XSFTYPE(STAROBJ(I)) .EQ. 'COMP') THEN
                CALL XCOPY(C8VAL, 1, HEAP(HPSF(STAROBJ(I))), IREF)
                IF      (STARVAL(I) .EQ. ' '
     &              .OR. STARVAL(I) .EQ. 'AMPL') THEN
                  CALL XPHASE(C8VAL, R8VAL, R8TMP)
                ELSE IF (STARVAL(I) .EQ. 'PHAS') THEN
                  CALL XPHASE(C8VAL, R8TMP, R8VAL)
                ELSE IF (STARVAL(I) .EQ. 'REAL') THEN
                  R8VAL =  DBLE(C8VAL)
                ELSE IF (STARVAL(I) .EQ. 'IMAG') THEN
                  R8VAL = DIMAG(C8VAL)
                ELSE
                  CALL WRNDIE(-5, 'XWRSTAR', 'Fatal Coding Error')
                  CALL DIE
                END IF
                IF (R8VAL .EQ. DFLOAT(0)) THEN
                  STARBUF(I) = '0'
                ELSE
                  WRITE(STARBUF(I), '(G80.6)') R8VAL
                END IF
              ELSE
                CALL WRNDIE(-5, 'XWRSTAR', 'Fatal Coding Error')
                CALL DIE
              END IF
            END IF
            LSTARBUF(I) = LEN(STARBUF(I))
            CALL TRIML(STARBUF(I), LSTARBUF(I))
            CALL TRIMM(STARBUF(I), LSTARBUF(I))
          END DO
C
          FIRST = 1
          LAST = FIRST
          SUM = 0
          DO CURR = 1, NSTAR
            SUM = SUM + 1 + LSTARBUF(CURR)
            IF (SUM .LE. 78) LAST = CURR
            IF (SUM .GT. 78 .OR. CURR .EQ. NSTAR) THEN
              WRITE(OUNIT, '(80(1X, A))')
     &          (STARBUF(I)(:LSTARBUF(I)), I = FIRST, LAST)
              FIRST = CURR
              IF (LAST .EQ. CURR) THEN
                FIRST = FIRST + 1
                SUM = 0
              ELSE IF (CURR .EQ. NSTAR) THEN
                WRITE(OUNIT, '(1X, A)') STARBUF(CURR)(:LSTARBUF(CURR))
              ELSE
                SUM = 1 + LSTARBUF(CURR)
              END IF
              LAST = FIRST
            END IF
          END DO
        END IF
      END DO
C
      RETURN
      END
C}
