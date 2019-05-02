      SUBROUTINE XRSPECL(HPATOF,HPINDF,HPANOMFLAG,
     &     QASELE,
     &     XRNATF,XRFDP,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &     XRSYMM,XRITSY,XRTR,XRINTR)
C
C Computes multiplicity for selected atoms and stores in specified
C atom object.
C Indicates selected atoms at general positions by a "1"
C           selected atoms at special positions by multiplicity > 1
C           unselected atoms by a "0"
C
C Author: Axel T. Brunger
C =======================
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'comand.inc'
      INTEGER HPATOF, HPINDF, HPANOMFLAG
      LOGICAL QASELE
      INTEGER XRNATF
      DOUBLE PRECISION XRFDP(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3)
C local
      CHARACTER*(WDMAX) ARRAYTO
      INTEGER LARRAYTO
      INTEGER NATSELE
C pointer
      INTEGER ATSELE, HPFQS
C
C begin
C allocate space for atom selection array
      ATSELE=ALLHP(INTEG4(NATOM))
C
C default values
      NATSELE=0
      ARRAYTO=' '
C
C command parsing
      CALL PUSEND('SPECial>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('SPEcial>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-special')
C
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(HEAP(ATSELE),NATSELE,X,Y,Z,.TRUE.)
      ELSE IF (WD(1:2).EQ.'TO') THEN
      CALL NEXTST('TO=',WD)
      CALL COPYST(ARRAYTO,WDMAX,LARRAYTO,WD,WDLEN)
      ELSE
      CALL CHKEND('SPECial>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      IF (NATSELE.EQ.0) THEN
      WRITE(6,'(3A)') ' %XSPEcial-ERR: zero atoms selected.'
      CALL WRNDIE(-5,'XSPEcial','zero atoms selected.')
      ELSEIF (ARRAYTO.EQ.' ') THEN
      WRITE(6,'(3A)') ' %XSPEcial-ERR: atom object not specified.'
      CALL WRNDIE(-5,'XSPEcial','atom object missing.')
      ELSE
      HPFQS=ALLHP(IREAL8(NATOM))
      CALL XRSPEC2(HEAP(ATSELE),HEAP(HPANOMFLAG),HEAP(HPATOF),
     &                   HEAP(HPINDF),
     &                   QASELE,XRNATF,XRFDP,QHERM,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &                   XRTR,XRINTR,
     &                   ARRAYTO,HEAP(HPFQS))
      CALL FREHP(HPFQS,IREAL8(NATOM))
      END IF
C free space for atom selection array
      CALL FREHP (ATSELE,INTEG4(NATOM))
C
      RETURN
      END
C==================================================================
      SUBROUTINE XRSPEC2(ATSELE,ANOMFLAG,XRATOF,XRINDF,
     &                   QASELE,XRNATF,XRFDP,QHERM,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &                   XRTR,XRINTR,
     &                   ARRAYTO,XRFQS)
C
C See routine XRSPECL above.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
      INTEGER ATSELE(*), ANOMFLAG(*), XRATOF(*), XRINDF(*)
      LOGICAL QASELE
      INTEGER XRNATF
      DOUBLE PRECISION XRFDP(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3)
      CHARACTER*(*) ARRAYTO
      DOUBLE PRECISION XRFQS(*)
C
C local
      INTEGER XRNATO, NANOM, I
C pointer
      INTEGER HPATOM, HPINDX
C
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
C
C allocate space for temporary atom flag arrays
      HPATOM=ALLHP(INTEG4(NATOM))
      HPINDX=ALLHP(INTEG4(NATOM))
      DO I=1,NATOM
      XRFQS(I)=ONE
      END DO
C
C call associate routine in order to get info about special
C positions
      CALL XRASSOC(ATSELE,ANOMFLAG,
     &            XRATOF,XRINDF,HEAP(HPATOM),
     &            HEAP(HPINDX),XRFQS,XRNATO,NANOM,
     &            5,'FCALC',QASELE,XRNATF,XRFDP,QHERM,
     &            XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &            XRTR,XRINTR)
C
C take the inverse of the XRFQS array to get the
C multiplicities
      DO I=1,NATOM
      IF (ABS(XRFQS(I)).GT.RSMALL) THEN
      XRFQS(I)=ONE/XRFQS(I)
      END IF
      END DO
C store XRFQS array in specified atom object
      CALL CPATPR(ARRAYTO,XRFQS)
C
      CALL FREHP(HPATOM,INTEG4(NATOM))
      CALL FREHP(HPINDX,INTEG4(NATOM))
C
      RETURN
      END
