      SUBROUTINE XDECLARE(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   NRHO,IRHO,XRNREF)
C
C Routine declares real space or reciprocal space objects.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'consta.inc'
      INTEGER XNAMEMX, XSFMX, XRHOMX
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER NRHO, IRHO
      INTEGER XRNREF
C local
      CHARACTER*4 SDOMAIN, STYPE
      CHARACTER*(WDMAX) SNAME
      INTEGER SNAMEL, I
      LOGICAL OK, QTOUCH
      LOGICAL FOUND
      INTEGER IEXIST
C parameter
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C begin
C
C set defaults
      SDOMAIN=' '
      STYPE=' '
      SNAME=' '
C
      CALL PUSEND('DECLare>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('DECLare>')
      CALL MISCOM('DECLare>',USED)
      IF (.NOT.USED) THEN
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-declare')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DOMA') THEN
      CALL NEXTA4('DOMAin=',SDOMAIN)
      IF (SDOMAIN.NE.'RECI'.AND.SDOMAIN.NE.'REAL') THEN
      CALL DSPERR('XDECLARE>','unknown domain.')
      END IF
      IF (SDOMAIN.EQ.'MAP') SDOMAIN='REAL'
      IF (SDOMAIN.EQ.'STRU') SDOMAIN='RECI'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TYPE') THEN
      CALL NEXTA4('TYPE=',STYPE)
      IF (STYPE.NE.'INTE'.AND.
     &    STYPE.NE.'COMP'.AND.STYPE.NE.'REAL') THEN
      CALL DSPERR('XDECLARE>','unknown type.')
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NAME') THEN
      CALL NEXTST('NAME=',WD)
C
C check that "name" begins with a letter and that it does not
C contain any special characters
      CALL CHKWRD(WD,WDLEN,OK)
      IF (.NOT.OK) THEN
      CALL DSPERR('XDECLARE>',
     & 'Invalid name: must start with letter and no special chars.')
      END IF
      CALL COPYST(SNAME,XNAMEMX,SNAMEL,WD,WDLEN)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?') THEN
C
      IF (XSFNUM.GT.0) THEN
      WRITE(6,'(A)')
     & ' XDECLARE: list of declared reciprocal space objects.'
      WRITE(6,'(2A)')
     &  ' XDECLARE: object name      type ',
     &  '   allocated       in use (<>0) '
      DO I=1,XSFNUM
C
      QTOUCH=.FALSE.
      IF (HPSF(I).NE.0) THEN
      CALL XQTOUCH(QTOUCH,XRNREF,XSFNUM,HPSF,XSFTYPE,I)
      END IF
      IF (HPSF(I).NE.0.AND.QTOUCH) THEN
      WRITE(6,'(5A)')
     &  '           ',XSFNAM(I),'       ',XSFTYPE(I),
     &  '       yes            yes'
      ELSEIF (HPSF(I).NE.0.AND..NOT.QTOUCH) THEN
      WRITE(6,'(5A)')
     &  '           ',XSFNAM(I),'       ',XSFTYPE(I),
     &  '       yes            no'
      ELSEIF (HPSF(I).EQ.0) THEN
      WRITE(6,'(5A)')
     &  '           ',XSFNAM(I),'       ',XSFTYPE(I),
     &  '       no             no'
      END IF
      END DO
      ELSE
      WRITE(6,'(A)')
     & ' XDECLARE: no reciprocal space objects declared.'
      END IF
C
      IF (XRHONUM.GT.0) THEN
      WRITE(6,'(A)')
     & ' XDECLARE: list of declared real space objects.'
      WRITE(6,'(2A)')
     &  ' XDECLARE: object name      ',
     &  '   allocated  '
      DO I=1,XRHONUM
      IF (HPRRHO(I).NE.0) THEN
      WRITE(6,'(5A)')
     & '           ',XRHONAM(I),'   ','    ',
     & '      yes'
      ELSE
      WRITE(6,'(5A)')
     & '           ',XRHONAM(I),'   ','    ',
     & '      no'
      END IF
      END DO
      ELSE
      WRITE(6,'(A)')
     & ' XDECLARE: no real space objects declared.'
      END IF
C
C=====================================================================
      ELSE
      CALL CHKEND('DECLare>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (SDOMAIN.EQ.' '.AND.STYPE.EQ.' '.AND.SNAME.EQ.' ') THEN
C don't do anything
      CONTINUE
      ELSEIF (SDOMAIN.EQ.'RECI') THEN
C
C reciprocal space object
      IF (STYPE.EQ.' ') THEN
      WRITE(6,'(A)') ' %XDECLARE-ERR: Object type not specified.'
      CALL WRNDIE(-5,'XDECLARE','Type needs to be specified.')
      ELSEIF (SNAME.EQ.' ') THEN
      WRITE(6,'(A)') ' %XDECLARE-ERR: Object name not specified.'
      CALL WRNDIE(-5,'XDECLARE','Name needs to be specified.')
      ELSE
C
      OK=.FALSE.
C
C
C does the object already exist?
      FOUND=.FALSE.
      IEXIST=0
      DO I=1,XSFNUM
      IF (SNAME.EQ.XSFNAM(I)) THEN
      FOUND=.TRUE.
      IEXIST=I
      END IF
      END DO
C
C not found -> OK to declare
      IF (.NOT.FOUND) THEN
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(3A)') ' XDECLARE: Object ',SNAME(1:SNAMEL),
     &  ' has been declared.'
      END IF
      OK=.TRUE.
      END IF
C
C don't allow re-declarations of user-defined objects
      IF (FOUND) THEN
      WRITE(6,'(3A)')
     & ' %XDECLARE-ERR: Object ',SNAME(1:SNAMEL),
     & ' already exists -- cannot be re-declared.'
      CALL WRNDIE(-5,'XDECLARE','Object cannot be re-declared.')
      END IF
C
      IF (OK) THEN
C
C declare object
      IF (XSFNUM.GE.XSFMX) THEN
      CALL WRNDIE(-5,'XDECLARE',
     & 'exceeded XSFMX parameter --> recompile program')
      ELSE
      XSFNUM=XSFNUM+1
      XSFNAM(XSFNUM)=SNAME
      XSFTYPE(XSFNUM)=STYPE
      XSFGNAM(XSFNUM)=0
      XSFGTYP(XSFNUM)=' '
      XSFGORD(XSFNUM)=0
      HPSF(XSFNUM)=0
      END IF
      END IF
C
      END IF
C
      ELSEIF (SDOMAIN.EQ.'REAL') THEN
C
C real space object
      IF (SNAME.EQ.' ') THEN
      WRITE(6,'(A)') ' %XDECLARE-ERR: Object name not specified.'
      CALL WRNDIE(-5,'XDECLARE','Object name not specified.')
      ELSEIF (STYPE.NE.' ') THEN
      WRITE(6,'(A)')
     & ' %XDECLARE-ERR: Type cannot be specified for real space obj.'
      CALL WRNDIE(-5,'XDECLARE',
     & 'Type cannot be specified for real space obj.')
      ELSE
C
C declare real space object
C does already exist?  If so, stop program.
      OK=.TRUE.
      DO I=1,XRHONUM
      IF (SNAME.EQ.XRHONAM(I)) THEN
      WRITE(6,'(3A)')
     & ' %XDECLARE-ERR: Object ',SNAME(1:SNAMEL),
     & ' already exists. '
      CALL WRNDIE(-5,'XDECLARE','Object cannot be re-declared.')
      OK=.FALSE.
      END IF
      END DO
C
      IF (OK) THEN
      IF (XRHONUM.GE.XRHOMX) THEN
      CALL WRNDIE(-5,'XDECLARE',
     & 'exceeded XRHOMX parameter --> recompile program')
      ELSE
      XRHONUM=XRHONUM+1
      XRHONAM(XRHONUM)=SNAME
      HPRRHO(XRHONUM)=0
      HPIRHO(XRHONUM)=0
      END IF
      END IF
C
      END IF
C
      ELSE
      WRITE(6,'(3A)') ' %XDECLARE-ERR: domain name "',SDOMAIN,
     &     '" missing or unknown.'
      CALL WRNDIE(-5,'XDECLARE','domain name missing or unknown.')
      END IF
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XUDECLA(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   NRHO,IRHO,XRMREF)
C
C Routine undeclares real space or reciprocal space objects.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER XNAMEMX, XSFMX, XRHOMX
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER NRHO, IRHO
      INTEGER XRMREF
C local
      CHARACTER*4 SDOMAIN
      CHARACTER*(WDMAX) SNAME
      INTEGER I, J, SNAMEL, GRPID
      LOGICAL FOUND, QEQUAL, MATCH
C begin
C
C set defaults
      SDOMAIN=' '
      SNAME=' '
C
      CALL PUSEND('UNDEclare>')
      DO WHILE (.NOT.DONE)
      CALL NEXTQL('UNDEclare>')
      CALL MISCOM('UNDEclare>',USED)
      IF (.NOT.USED) THEN
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-undeclare')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DOMA') THEN
      CALL NEXTA4('DOMAin=',SDOMAIN)
      IF (SDOMAIN.NE.'STRU'.AND.SDOMAIN.NE.'MAP'.AND.
     &    SDOMAIN.NE.'RECI'.AND.SDOMAIN.NE.'REAL') THEN
      CALL DSPERR('XDECLARE>','unknown domain.')
      END IF
      IF (SDOMAIN.EQ.'MAP') SDOMAIN='REAL'
      IF (SDOMAIN.EQ.'STRU') SDOMAIN='RECI'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NAME') THEN
      CALL NEXTQL('NAME')
      QEQUAL=.TRUE.
      IF (WD(1:WDLEN).EQ.'=') THEN
      QEQUAL=.TRUE.
      ELSEIF (WD(1:WDLEN).EQ.'#') THEN
      QEQUAL=.FALSE.
      ELSE
      CALL DSPERR('XUDECLA','operator missing (= or #).')
      END IF
      CALL NEXTWD('NAME=')
      CALL COPYST(SNAME,XNAMEMX,SNAMEL,WD,WDLEN)
C=====================================================================
      ELSE
      CALL CHKEND('UNDEclare>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (SDOMAIN.EQ.' ') THEN
      WRITE(6,'(A)') ' %XUDECLA-ERR: Domain not specified.'
      CALL WRNDIE(-5,'XUDECLA','Domain not specified.')
      ELSEIF (SNAME.EQ.' ') THEN
      WRITE(6,'(A)') '%XUDECLA-ERR Object name not specified.'
      CALL WRNDIE(-5,'XUDECLA','Object name not specified.')
C
      ELSEIF (SDOMAIN.EQ.'RECI') THEN
C
      FOUND=.FALSE.
C
      I=1
      DO WHILE (I.LE.XSFNUM)
      CALL EQSTWC(XSFNAM(I),XNAMEMX,SNAME,SNAMEL,1,1,MATCH)
      IF ((QEQUAL.AND.MATCH).OR.(.NOT.QEQUAL.AND..NOT.MATCH)) THEN
C undeclare object
      FOUND=.TRUE.
      GRPID=XSFGNAM(I)
C
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(3A)')
     & ' Object ',XSFNAM(I),' has been undeclared.'
      END IF
C
C free-up space if allocated
      IF (HPSF(I).NE.0.AND.XSFTYPE(I).EQ.'REAL') THEN
      CALL FREHP(HPSF(I),IREAL8(XRMREF))
      ELSEIF (HPSF(I).NE.0.AND.XSFTYPE(I).EQ.'COMP') THEN
      CALL FREHP(HPSF(I),ICPLX8(XRMREF))
      ELSEIF (HPSF(I).NE.0.AND.XSFTYPE(I).EQ.'INTE') THEN
      CALL FREHP(HPSF(I),INTEG4(XRMREF))
      END IF
C
      DO J=I+1,XSFNUM
      XSFNAM(J-1)=XSFNAM(J)
      XSFTYPE(J-1)=XSFTYPE(J)
      HPSF(J-1)=HPSF(J)
      XSFGNAM(J-1)=XSFGNAM(J)
      XSFGTYP(J-1)=XSFGTYP(J)
      XSFGORD(J-1)=XSFGORD(J)
      END DO
      XSFNUM=XSFNUM-1
      I=I-1
C
C remove group definition that the object belonged to
      IF (GRPID.NE.0) THEN
      WRITE(6,'(A,I4,A)') 'XUDECLA-info Object belongs to group ',
     &     XSFGNAM(I+1), '.  The group definition is removed!'
      DO J=1,XSFNUM
      IF (XSFGNAM(J).EQ.GRPID) THEN
      XSFGNAM(J)=0
      XSFGTYP(J)=' '
      XSFGORD(J)=0
      END IF
      END DO
      END IF
      END IF
      I=I+1
      END DO
C
      IF (.NOT.FOUND) THEN
      WRITE(6,'(2A)')
     &  ' XUDELCA-info: no reciprocal space objects found matching ',
     &   SNAME(1:SNAMEL)
      END IF
C
C
      ELSEIF (SDOMAIN.EQ.'REAL') THEN
C
C does already exist?  If so, free it up and undeclare.
C
      FOUND=.FALSE.
      I=1
      DO WHILE (I.LE.XRHONUM)
      CALL EQSTWC(XRHONAM(I),XNAMEMX,SNAME,SNAMEL,1,1,MATCH)
      IF ((QEQUAL.AND.MATCH).OR.(.NOT.QEQUAL.AND..NOT.MATCH)) THEN
      FOUND=.TRUE.
C
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(3A)')
     & ' Object ',XRHONAM(I),' has been undeclared.'
      END IF
C
      IF (HPRRHO(I).NE.0) CALL FREHP(HPRRHO(I),IREAL4(NRHO))
      IF (HPIRHO(I).NE.0) CALL FREHP(HPIRHO(I),IREAL4(IRHO))
      DO J=I+1,XRHONUM
      XRHONAM(J-1)=XRHONAM(J)
      HPRRHO(J-1)=HPRRHO(J)
      HPIRHO(J-1)=HPIRHO(J)
      END DO
      XRHONUM=XRHONUM-1
      I=I-1
      END IF
      I=I+1
      END DO
C
      IF (.NOT.FOUND) THEN
      WRITE(6,'(2A)')
     &  ' XUDELCA-INFO: no real space object(s) found matching',
     &   SNAME(1:SNAMEL)
      END IF
C
      END IF
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XQUERY(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   NRHO,IRHO,XRNREF)
C
C Routine queries info about real space or reciprocal space objects.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER XNAMEMX, XSFMX, XRHOMX
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER NRHO, IRHO
      INTEGER XRNREF
C local
      CHARACTER*4 SDOMAIN
      CHARACTER*(WDMAX) SNAME
      INTEGER I, SNAMEL
      LOGICAL FOUND, QTOUCH
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
      CHARACTER*7 STYPE
C begin
C
C set defaults
      SDOMAIN=' '
      SNAME=' '
C
      CALL PUSEND('QUERy>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('QUERy>')
      CALL MISCOM('QUERy>',USED)
      IF (.NOT.USED) THEN
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-query')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DOMA') THEN
      CALL NEXTA4('DOMAin=',SDOMAIN)
      IF (SDOMAIN.NE.'STRU'.AND.SDOMAIN.NE.'MAP'.AND.
     &    SDOMAIN.NE.'RECI'.AND.SDOMAIN.NE.'REAL') THEN
      CALL DSPERR('XDECLARE>','unknown domain.')
      END IF
      IF (SDOMAIN.EQ.'MAP') SDOMAIN='REAL'
      IF (SDOMAIN.EQ.'STRU') SDOMAIN='RECI'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NAME') THEN
      CALL NEXTST('NAME=',WD)
      CALL COPYST(SNAME,XNAMEMX,SNAMEL,WD,WDLEN)
C=====================================================================
      ELSE
      CALL CHKEND('QUERy>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (SDOMAIN.EQ.' '.AND.SNAME.EQ.' ') THEN
C don't do anything
      CONTINUE
      ELSEIF (SDOMAIN.EQ.' ') THEN
      WRITE(6,'(A)') ' %XQUERY-ERR: Domain not specified.'
      CALL WRNDIE(-5,'XQUERY','Domain not specified.')
      ELSEIF (SNAME.EQ.' ') THEN
      WRITE(6,'(A)') '%XQUERY-ERR Object name not specified.'
      CALL WRNDIE(-5,'XQUERY','Object name not specified.')
C
      ELSEIF (SDOMAIN.EQ.'RECI') THEN
C
      FOUND=.FALSE.
C
      DO I=1,XSFNUM
      IF (SNAME.EQ.XSFNAM(I)) THEN
      FOUND=.TRUE.
      CALL DECLAR('OBJECT_EXIST', 'LO' ,'TRUE', DCVAL, DPVAL )
C
      IF (XSFTYPE(I).EQ.'REAL') THEN
      STYPE='REAL'
      ELSEIF (XSFTYPE(I).EQ.'COMP') THEN
      STYPE='COMPLEX'
      ELSEIF (XSFTYPE(I).EQ.'INTE') THEN
      STYPE='INTEGER'
      END IF
C
      CALL DECLAR('OBJECT_TYPE', 'ST' ,STYPE, DCVAL, DPVAL )
C
      CALL XQTOUCH(QTOUCH,XRNREF,XSFNUM,HPSF,XSFTYPE,I)
      IF (QTOUCH) THEN
      CALL DECLAR('OBJECT_USED', 'LO' ,'TRUE', DCVAL, DPVAL )
      ELSE
      CALL DECLAR('OBJECT_USED', 'LO' ,'FALSE', DCVAL, DPVAL )
      END IF
C
C group stuff
      DPVAL=XSFGNAM(I)
      CALL DECLAR('OBJECT_GROUP', 'DP' ,' ' , DCVAL, DPVAL )
C
      CALL DECLAR('OBJECT_GTYPE', 'ST' ,XSFGTYP(I), DCVAL, DPVAL )
      DPVAL=XSFGORD(I)
      CALL DECLAR('OBJECT_GORD', 'DP' ,' ', DCVAL, DPVAL )
C
      IF (WRNLEV.GE.5) THEN
      IF (QTOUCH) THEN
      WRITE(6,'(4A)')
     & ' Reciprocal space object ',SNAME(1:SNAMEL),
     & ' exists, is used, and has type ',STYPE
      IF (XSFGNAM(I).EQ.0) THEN
      WRITE(6,'(4A)')
     & ' Object does not belong to any group.'
      ELSE
      WRITE(6,'(A,I4,3A,I4)')
     & ' Object belongs to group ',XSFGNAM(I),' of group-type ',
     &   XSFGTYP(I),' order=',XSFGORD(I)
      END IF
      ELSE
      WRITE(6,'(4A)')
     & ' Reciprocal space object ',SNAME(1:SNAMEL),
     & ' exists, is unused, and has type ',STYPE
      END IF
      END IF
C
      END IF
      END DO
C
      IF (.NOT.FOUND) THEN
      CALL DECLAR('OBJECT_EXIST', 'LO' ,'FALSE', DCVAL, DPVAL )
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(3A)')
     & ' Reciprocal space object ',SNAME(1:SNAMEL),
     & ' does not exist.'
      END IF
      END IF
C
C
      ELSEIF (SDOMAIN.EQ.'REAL') THEN
C
C does already exist?  If so, free it up and re-declare.
C
      FOUND=.FALSE.
      DO I=1,XRHONUM
      IF (SNAME.EQ.XRHONAM(I)) THEN
      FOUND=.TRUE.
C
      CALL DECLAR('OBJECT_EXIST', 'LO' ,'TRUE', DCVAL, DPVAL )
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(3A)')
     & ' Real space object ',SNAME(1:SNAMEL),' exists.'
      END IF
C
      END IF
      END DO
C
      IF (.NOT.FOUND) THEN
      CALL DECLAR('OBJECT_EXIST', 'LO' ,'FALSE', DCVAL, DPVAL )
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(3A)')
     & ' Real space object ',SNAME(1:SNAMEL),' does not exist.'
      END IF
      END IF
C
      END IF
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XRERES
C
C subroutine initializes or resets the x-ray refinement list
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'consta.inc'
C local
      INTEGER I
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION TEMP
C parameters
      DOUBLE PRECISION ZERO, ONE, THREE, TEN, NINETY, HALF, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, THREE=3.0D0, TEN=10.0D0)
      PARAMETER (NINETY=90.0D0, HALF=0.5D0, TWO=2.0D0)
C begin
C
C ============
C UPDATE FLAGS
C ============
C
C when XRQCHK gets set, the X-RAY term needs to be
C recomputed (otherwise a linear interpolation is performed)
      XRQCHK=.TRUE.
C
C when XRREUP gets set, the target selection has
C to be updated.
      XRREUP=.TRUE.
C
C when XRRED is set, the reflection list has to
C be mapped into an asymmetric unit.  Multiplicities
C and reflection type are also defined.
      XRRED=.TRUE.
C
C when XRUPAT gets set, the atomic selections need
C to be re-defined (this is required, for example, when
C the molecular structure file has changed).
      XRUPAT=.TRUE.
C
C when XRMAP gets set, the asymmetric unit for maps
C needs to be defined.
      XRMAP=.TRUE.
C
C set fft method
      QFFT=.TRUE.
C look-up tables turned on
      QLOOK=.TRUE.
C tolerance
      XRLTOL=HALF
C number of reflections
      XRNREF=0
C reflection list pointers into HEAP
      XRMREF=0
      HPH=0
      HPK=0
      HPL=0
      HPMULT=0
      HPTYPE=0
      HPTSEL=0
C
C initialize user-defined reciprocal space objects
      DO I=1,XSFMX
      HPSF(I)=0
      XSFNAM(I)=' '
      XSFTYPE(I)=' '
      XSFGNAM(I)=0
      XSFGTYP(I)=' '
      XSFGORD(I)=0
      END DO
C
C atom list pointers into HEAP
      HPANOMFLAG=0
      IHPFLAG=0
      DO I=1,MAXHPFLAG
      HPFLAG(I)=0
      XASSOC(I)=' '
      END DO
      HPATOF=0
      HPINDF=0
      HPFX=0
      HPFY=0
      HPFZ=0
      HPDX=0
      HPDY=0
      HPDZ=0
      HPDT=0
      HPDQ=0
C size of maps (NRHO) and size of asymmetric unit (NMASK)
      NRHO=0
      NMASK=0
C map pointers into HEAP
      HPRHOMA=0
C
C initialize reciprocal space objects
      XSFNUM=0
C
C initialize real space objects
      DO I=1,XRHOMX
      HPRRHO(I)=0
      HPIRHO(I)=0
      XRHONAM(I)=' '
      END DO
      XRHONUM=0
C
C unit cell
      XRCELL(1)=ONE
      XRCELL(2)=ONE
      XRCELL(3)=ONE
      XRCELL(4)=NINETY
      XRCELL(5)=NINETY
      XRCELL(6)=NINETY
      CALL XRFRAC(XRCELL,XRTR,XRINTR,XRVOL)
C hermitian symmetry
      QHERM=.TRUE.
C reset symmetry operator database
      CALL XSYMRES(QHERM,XRSYTH,XRMSYM,XRNSYM,XRSYMM,
     &             XRITSY,ARPNN)
C atomic scattering parameters
      XRSN=0
C scattering atoms
      XRNATF=0
      QASELE=.FALSE.
C sum (Fo-Fc)**2 scale factor
      XRSCAL=ONE
C number of bins
      MBINS=8
      TEMP=MBINS
      CALL DECLAR( 'BIN_NUMBER', 'DP', ' ', DBCOMP, TEMP)
C bin mode
      BINMODE=1
      QBINSET=.FALSE.
C default bin resolution range: none
      XBINLOW=ZERO
      XBINHIGH=ZERO
C default map resolution: none
      MAPR=ZERO
C asymmetric unit
      ARPNN=0
C target function selection expression
      SRPNN=0
C cross-validation selection expression
      CRPNN=0
C initialize target, dtarget, and monitor expressions
      CALL XTARMOIN
C reset FFT parameters
      CALL XFFTIN
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
C initialize FMAP module
      CALL FMAP(0, NA, NB, NC, NRHO, IRHO, QHERM,
     &          0, 0, 0, 0, 0, 0, 0, 0,
     &          XRHONUM, XRHONAM, HPRRHO, HPIRHO,
     &          XRNSYM, XRSYTH)
C initialize MAPYARD module
      CALL MAPYARD(0, NA, NB, NC, NRHO, IRHO, QHERM,
     &             XRHONUM, XRHONAM, HPRRHO, HPIRHO,
     &             XRNSYM)
C initialize TSMAP module (compact fast translation search)
      CALL TSLMINI
C initialize PSEARCH module
      CALL INIPSTCOMMON(.FALSE.)
C=====================================================================
C #endif
C=====================================================================
C quick expansion/reduction mode
      XQUICK=.FALSE.
      XQEXPAN=.FALSE.
      RETURN
      END
C=====================================================================
      SUBROUTINE XRAFRE
C
C routine frees up all reciprocal space objects
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
C begin
      CALL XRAREF(0,XRMREF,XRNREF,HPTSEL,
     &                  HPH,HPK,HPL,HPMULT,HPTYPE,
     &                  XSFNUM,HPSF,XSFTYPE)
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
C free up FlagMap
      CALL FMAP(-1, NA, NB, NC, NRHO, IRHO, QHERM,
     &          0, 0, 0, 0, 0, 0, 0, 0,
     &          XRHONUM, XRHONAM, HPRRHO, HPIRHO,
     &          XRNSYM, XRSYTH)
C free up MapYard
      CALL MAPYARD(-1, NA, NB, NC, NRHO, IRHO, QHERM,
     &             XRHONUM, XRHONAM, HPRRHO, HPIRHO,
     &             XRNSYM)
C free up TSMAP arrays (compact fast translation search)
      CALL TSLMRST
C free up PSEARCH arrays
      CALL INIPSTCOMMON(.FALSE.)
C=====================================================================
C #endif
C=====================================================================
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XRAREF(NEWREF,XRMREF,XRNREF,HPTSEL,
     &                  HPH,HPK,HPL,HPMULT,HPTYPE,
     &                  XSFNUM,HPSF,XSFTYPE)
C
C expands existing reflection list or deletes
C all reciprocal space objects.
C
C Author: Axel T. Brunger
C =======================
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER NEWREF
      INTEGER XRMREF, XRNREF, HPTSEL
      INTEGER HPH, HPK, HPL, HPMULT, HPTYPE
      INTEGER XSFNUM, HPSF(*)
      CHARACTER*(*) XSFTYPE(*)
C local
      INTEGER I, NEW
C parameter
      DOUBLE PRECISION ZERO, ONE
      DOUBLE COMPLEX CZERO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
      CZERO=ZERO
C
      IF (NEWREF.EQ.0) THEN
C
C wipe everything out
C check whether space was already allocated
      IF (HPL.NE.0) CALL FREHP(HPL,INTEG4(XRMREF))
      IF (HPK.NE.0) CALL FREHP(HPK,INTEG4(XRMREF))
      IF (HPH.NE.0) CALL FREHP(HPH,INTEG4(XRMREF))
C
      IF (HPMULT.NE.0) CALL FREHP(HPMULT,INTEG4(XRMREF))
      IF (HPTYPE.NE.0) CALL FREHP(HPTYPE,INTEG4(XRMREF))
      IF (HPTSEL.NE.0) CALL FREHP(HPTSEL,INTEG4(XRMREF))
C
C free-up all objects
      DO I=1,XSFNUM
      IF (HPSF(I).NE.0.AND.XSFTYPE(I).EQ.'REAL') THEN
      CALL FREHP(HPSF(I),IREAL8(XRMREF))
      ELSEIF (HPSF(I).NE.0.AND.XSFTYPE(I).EQ.'COMP') THEN
      CALL FREHP(HPSF(I),ICPLX8(XRMREF))
      ELSEIF (HPSF(I).NE.0.AND.XSFTYPE(I).EQ.'INTE') THEN
      CALL FREHP(HPSF(I),INTEG4(XRMREF))
      END IF
      END DO
C
C
      XRMREF=0
      XRNREF=0
C
      HPH=0
      HPK=0
      HPL=0
C
      HPMULT=0
      HPTYPE=0
      HPTSEL=0
C
      DO I=1,XSFNUM
      HPSF(I)=0
      END DO
C
C
      ELSEIF (XRMREF.LT.NEWREF) THEN
C
C expand objects
C
      NEW=ALLHP(INTEG4(NEWREF))
      CALL FILL4(HEAP(NEW),NEWREF,0)
      IF (HPH.NE.0) CALL COPYI4(HEAP(HPH),HEAP(NEW),XRMREF)
      IF (HPH.NE.0) CALL FREHP(HPH,INTEG4(XRMREF))
      HPH=NEW
C
      NEW=ALLHP(INTEG4(NEWREF))
      CALL FILL4(HEAP(NEW),NEWREF,0)
      IF (HPK.NE.0) CALL COPYI4(HEAP(HPK),HEAP(NEW),XRMREF)
      IF (HPK.NE.0) CALL FREHP(HPK,INTEG4(XRMREF))
      HPK=NEW
C
      NEW=ALLHP(INTEG4(NEWREF))
      CALL FILL4(HEAP(NEW),NEWREF,0)
      IF (HPL.NE.0) CALL COPYI4(HEAP(HPL),HEAP(NEW),XRMREF)
      IF (HPL.NE.0) CALL FREHP(HPL,INTEG4(XRMREF))
      HPL=NEW
C
      NEW=ALLHP(INTEG4(NEWREF))
      CALL FILL4(HEAP(NEW),NEWREF,0)
      IF (HPMULT.NE.0) CALL COPYI4(HEAP(HPMULT),HEAP(NEW),XRMREF)
      IF (HPMULT.NE.0) CALL FREHP(HPMULT,INTEG4(XRMREF))
      HPMULT=NEW
C
      NEW=ALLHP(INTEG4(NEWREF))
      CALL FILL4(HEAP(NEW),NEWREF,0)
      IF (HPTYPE.NE.0) CALL COPYI4(HEAP(HPTYPE),HEAP(NEW),XRMREF)
      IF (HPTYPE.NE.0) CALL FREHP(HPTYPE,INTEG4(XRMREF))
      HPTYPE=NEW
C
      NEW=ALLHP(INTEG4(NEWREF))
      CALL FILL4(HEAP(NEW),NEWREF,0)
      IF (HPTSEL.NE.0) CALL COPYI4(HEAP(HPTSEL),HEAP(NEW),XRMREF)
      IF (HPTSEL.NE.0) CALL FREHP(HPTSEL,INTEG4(XRMREF))
      HPTSEL=NEW
C
C expand all reciprocal space objects
      DO I=1,XSFNUM
      IF (HPSF(I).NE.0.AND.XSFTYPE(I).EQ.'REAL') THEN
C
      NEW=ALLHP(IREAL8(NEWREF))
      CALL FILLR8(HEAP(NEW),NEWREF,ZERO)
      CALL COPYR8(HEAP(HPSF(I)),HEAP(NEW),XRMREF)
      CALL FREHP(HPSF(I),IREAL8(XRMREF))
      HPSF(I)=NEW
C
      ELSEIF (HPSF(I).NE.0.AND.XSFTYPE(I).EQ.'COMP') THEN
C
      NEW=ALLHP(ICPLX8(NEWREF))
      CALL FILLC8(HEAP(NEW),NEWREF,CZERO)
      CALL COPYC8(HEAP(HPSF(I)),HEAP(NEW),XRMREF)
      CALL FREHP(HPSF(I),ICPLX8(XRMREF))
      HPSF(I)=NEW
C
      ELSEIF (HPSF(I).NE.0.AND.XSFTYPE(I).EQ.'INTE') THEN
C
      NEW=ALLHP(INTEG4(NEWREF))
      CALL FILL4(HEAP(NEW),NEWREF,0)
      CALL COPYI4(HEAP(HPSF(I)),HEAP(NEW),XRMREF)
      CALL FREHP(HPSF(I),INTEG4(XRMREF))
      HPSF(I)=NEW
C
      END IF
      END DO
C
      XRMREF=NEWREF
C
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE XSFAL(POINT,XRMREF,MODE)
C
C Routine allocates space for reciprocal space objects.
C Routine returns heap pointers.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INTEGER POINT, XRMREF
      CHARACTER*(*) MODE
C parameter
      DOUBLE COMPLEX CZERO
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0)
C local
C begin
      CZERO=DCMPLX(ZERO,ZERO)
C
      IF (MODE.EQ.'SF'.OR.MODE.EQ.'COMP') THEN
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     &  ' XSFAL: allocating space for complex reciprocal space object.'
      END IF
      POINT=ALLHP(ICPLX8(XRMREF))
      CALL FILLC8(HEAP(POINT),XRMREF,CZERO)
      ELSEIF (MODE.EQ.'HL'.OR.MODE.EQ.'REAL') THEN
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     &   ' XSFAL: allocating space for real reciprocal space object.'
      END IF
      POINT=ALLHP(IREAL8(XRMREF))
      CALL FILLR8(HEAP(POINT),XRMREF,ZERO)
      ELSEIF (MODE.EQ.'INTE') THEN
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     & ' XSFAL: allocating space for integer reciprocal space object.'
      END IF
      POINT=ALLHP(INTEG4(XRMREF))
      CALL FILL4(HEAP(POINT),XRMREF,0)
      END IF
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XRMAPR(L)
C
C Routine (re-)allocates space for the asymmetric unit for
C density maps
C
C Author: Axel T. Brunger
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'funct.inc'
      INTEGER L
C local
      INTEGER I
C begin
C
      IF (L.EQ.0.AND.HPRHOMA.NE.0) THEN
      WRITE(6,'(A)')
     &' XRMAPR: symmetry or unitcell changed. Real space obj. deleted.'
      END IF
C free-up space
      IF (HPRHOMA.NE.0) CALL FREHP(HPRHOMA,INTEG4(NRHO))
      DO I=1,XRHONUM
      IF (HPRRHO(I).NE.0) CALL FREHP(HPRRHO(I),IREAL4(NRHO))
      IF (HPIRHO(I).NE.0) CALL FREHP(HPIRHO(I),IREAL4(IRHO))
      END DO
C
      NRHO=0
      NMASK=0
      XRMAP=.TRUE.
      HPRHOMA=0
C
      DO I=1,XRHOMX
      HPRRHO(I)=0
      HPIRHO(I)=0
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMAPAL(POINTR,POINTI,QHERM,NRHO,IRHO)
C
C Routine allocates space for a map.  Routine
C returns heap pointers POINTR, POINTI
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INTEGER POINTR, POINTI
      LOGICAL QHERM
      INTEGER NRHO, IRHO
C parameter
      REAL SZERO
      PARAMETER (SZERO=0.0)
C local
C begin
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)') ' XMAPAL: allocating space for real space object.'
      END IF
C
C
C allocate HEAP space for a map.
      IF (QHERM) THEN
      IRHO=1
      ELSE
      IRHO=NRHO
      END IF
      POINTR=ALLHP(IREAL4(NRHO))
      POINTI=ALLHP(IREAL4(IRHO))
C
C fill with zero
      CALL FILLR4(HEAP(POINTR),NRHO,SZERO)
      IF (.NOT.QHERM) THEN
      CALL FILLR4(HEAP(POINTI),NRHO,SZERO)
      END IF
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XRAATM(NATOM)
C
C routine allocates space for various atom arrays
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INTEGER NATOM
C local
      INTEGER I
C begin
C
C check whether space was already allocated
      IF (HPDQ.NE.0) CALL FREHP(HPDQ,IREAL8(XRMATO))
      IF (HPDT.NE.0) CALL FREHP(HPDT,IREAL8(XRMATO))
      IF (HPDZ.NE.0) CALL FREHP(HPDZ,IREAL8(XRMATO))
      IF (HPDY.NE.0) CALL FREHP(HPDY,IREAL8(XRMATO))
      IF (HPDX.NE.0) CALL FREHP(HPDX,IREAL8(XRMATO))
      IF (HPFZ.NE.0) CALL FREHP(HPFZ,IREAL8(XRMATO))
      IF (HPFY.NE.0) CALL FREHP(HPFY,IREAL8(XRMATO))
      IF (HPFX.NE.0) CALL FREHP(HPFX,IREAL8(XRMATO))
      IF (HPINDF.NE.0) CALL FREHP(HPINDF,INTEG4(XRMATO))
      IF (HPATOF.NE.0) CALL FREHP(HPATOF,INTEG4(XRMATO))
      DO I=1,MAXHPFLAG
      IF (HPFLAG(I).NE.0) CALL FREHP(HPFLAG(I),INTEG4(XRMATO))
      END DO
      IF (HPANOMFLAG.NE.0) CALL FREHP(HPANOMFLAG,INTEG4(XRMATO))
C
      XRMATO=0
C atomic scattering parameters
      XRSN=0
C scattering atoms
      XRNATF=0
C
      HPANOMFLAG=0
      IHPFLAG=0
      DO I=1,MAXHPFLAG
      HPFLAG(I)=0
      XASSOC(I)=' '
      END DO
      HPATOF=0
      HPINDF=0
      HPFX=0
      HPFY=0
      HPFZ=0
      HPDX=0
      HPDY=0
      HPDZ=0
      HPDT=0
      HPDQ=0
C
C now allocate new space
      IF (NATOM.GT.0) THEN
      XRMATO=NATOM
      HPANOMFLAG=ALLHP(INTEG4(XRMATO))
      HPATOF=ALLHP(INTEG4(XRMATO))
      HPINDF=ALLHP(INTEG4(XRMATO))
      HPFX=ALLHP(IREAL8(XRMATO))
      HPFY=ALLHP(IREAL8(XRMATO))
      HPFZ=ALLHP(IREAL8(XRMATO))
      HPDX=ALLHP(IREAL8(XRMATO))
      HPDY=ALLHP(IREAL8(XRMATO))
      HPDZ=ALLHP(IREAL8(XRMATO))
      HPDT=ALLHP(IREAL8(XRMATO))
      HPDQ=ALLHP(IREAL8(XRMATO))
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE XRENAME(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   NRHO,IRHO,XRMREF,XRNREF)
C
C Routine renames real space or reciprocal space objects.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'consta.inc'
      INTEGER XNAMEMX, XSFMX, XRHOMX
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER NRHO, IRHO
      INTEGER XRMREF, XRNREF
C local
      CHARACTER*4 SDOMAIN
      CHARACTER*(WDMAX) SOLD, SNEW
      INTEGER LOLD, LNEW
      INTEGER I, IOLD, INEW
      LOGICAL ERR
      DOUBLE COMPLEX CZERO
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
C set defaults
      CZERO=DCMPLX(ZERO,ZERO)
      SOLD=' '
      SNEW=' '
      SDOMAIN=' '
C
      CALL PUSEND('RENAme>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('RENAme>')
      CALL MISCOM('RENAme>',USED)
      IF (.NOT.USED) THEN
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-rename')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DOMA') THEN
      CALL NEXTA4('DOMAin=',SDOMAIN)
      IF (SDOMAIN.NE.'STRU'.AND.SDOMAIN.NE.'MAP'.AND.
     &    SDOMAIN.NE.'RECI'.AND.SDOMAIN.NE.'REAL') THEN
      CALL DSPERR('XDECLARE>','unknown domain.')
      END IF
      IF (SDOMAIN.EQ.'MAP') SDOMAIN='REAL'
      IF (SDOMAIN.EQ.'STRU') SDOMAIN='RECI'
C=====================================================================
      ELSE IF (WD(1:WDLEN).EQ.'OLD') THEN
      CALL NEXTST('OLD=',WD)
      CALL COPYST(SOLD,XNAMEMX,LOLD,WD,WDLEN)
C=====================================================================
      ELSE IF (WD(1:WDLEN).EQ.'NEW') THEN
      CALL NEXTST('NEW=',WD)
      CALL COPYST(SNEW,XNAMEMX,LNEW,WD,WDLEN)
C=====================================================================
      ELSE
      CALL CHKEND('RENAme>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (SOLD.NE.SNEW) THEN
C
      IF (SOLD.EQ.' ') THEN
      WRITE(6,'(A)') ' %XRENAME-ERR: old object name missing.'
      CALL WRNDIE(-5,'XRENAME','old object name missing.')
      ELSEIF (SDOMAIN.EQ.' ') THEN
      WRITE(6,'(A)') ' %XRENAME-ERR: domain specification missing.'
      CALL WRNDIE(-5,'XRENAME','domain spec. missing.')
      ELSEIF (SNEW.EQ.' ') THEN
      WRITE(6,'(A)') ' %XRENAME-ERR: new object name missing.'
      CALL WRNDIE(-5,'XRENAME','new object name missing.')
      ELSEIF (SDOMAIN.EQ.'RECI') THEN
C
      ERR=.FALSE.
C
      IOLD=0
      DO I=1,XSFNUM
      IF (SOLD.EQ.XSFNAM(I)) THEN
      IOLD=I
      END IF
      END DO
C
      IF (IOLD.EQ.0) THEN
      WRITE(6,'(3A)') ' %XRENAME-ERR: object ',SOLD(1:LOLD),
     &    ' not found'
      CALL WRNDIE(-5,'XRENAME','object not found.')
      ERR=.TRUE.
      END IF
C
      INEW=0
      DO I=1,XSFNUM
      IF (SNEW.EQ.XSFNAM(I)) THEN
      INEW=I
      END IF
      END DO
C
      IF (INEW.NE.0) THEN
      WRITE(6,'(3A)') ' %XRENAME-ERR: new object name ',
     &   SNEW(1:LNEW),' already used.'
      CALL WRNDIE(-5,'XRENAME','new object name already used.')
      ERR=.TRUE.
      END IF
C
      IF (.NOT.ERR) THEN
C
      IF (INEW.EQ.0) THEN
C move user-defined object into new user-defined object
C
      XSFNAM(IOLD)=SNEW
C
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(4A)') ' XDECLARE: User object ',SOLD(1:LOLD),
     & ' has been renamed to user object ',SNEW(1:LNEW)
      END IF
C
      END IF
C
      END IF
C
C real space objects
C ==================
C
      ELSEIF (SDOMAIN.EQ.'REAL') THEN
C
      ERR=.FALSE.
C
      IOLD=0
      DO I=1,XRHONUM
      IF (SOLD.EQ.XRHONAM(I)) THEN
      IOLD=I
      END IF
      END DO
C
      IF (IOLD.EQ.0) THEN
      WRITE(6,'(3A)') ' %XRENAME-ERR: object ',SOLD(1:LOLD),
     &      ' not found'
      CALL WRNDIE(-5,'XRENAME','object not found.')
      ERR=.TRUE.
      END IF
C
      INEW=0
      DO I=1,XRHONUM
      IF (SNEW.EQ.XRHONAM(I)) THEN
      INEW=I
      END IF
      END DO
C
      IF (INEW.NE.0) THEN
      WRITE(6,'(3A)') ' %XRENAME-ERR: new object name ',
     &   SNEW(1:LNEW),' already used.'
      CALL WRNDIE(-5,'XRENAME','new object name already used.')
      ERR=.TRUE.
      END IF
C
      IF (.NOT.ERR) THEN
      XRHONAM(IOLD)=SNEW
      WRITE(6,'(4A)') ' XDECLARE: Real space object ',SOLD(1:LOLD),
     & ' has been renamed to ',SNEW(1:LNEW)
      END IF
C
C
      END IF
C
      END IF
C
      RETURN
      END
C===================================================================
      SUBROUTINE XQTOUCH(QTOUCH,XRNREF,XSFNUM,HPSF,
     &                   XSFTYPE,IOBJ)
C
C Routine checks if reciprocal object has any elements non-equal
C to zero.
C
C     Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      LOGICAL QTOUCH
      INTEGER XRNREF
      INTEGER XSFNUM, HPSF(*)
      CHARACTER*(*) XSFTYPE(*)
      INTEGER IOBJ
C local
      INTEGER REFLCT, ITEMP
      DOUBLE PRECISION ONE, TEMP
      DOUBLE COMPLEX CTEMP
      PARAMETER (ONE=1.0D0)
C begin
C
C check if reciprocal space object is not equal to zero
      QTOUCH=.FALSE.
      IF (HPSF(IOBJ).NE.0) THEN
      REFLCT=1
      DO WHILE (REFLCT.LE.XRNREF.AND..NOT.QTOUCH)
C
      IF (XSFTYPE(IOBJ).EQ.'COMP') THEN
      CALL XCOPY(CTEMP,1,HEAP(HPSF(IOBJ)),REFLCT)
      QTOUCH=QTOUCH.OR.(ABS(CTEMP).GT.RSMALL)
C
      ELSEIF (XSFTYPE(IOBJ).EQ.'REAL') THEN
      CALL XCOPYR(TEMP,1,HEAP(HPSF(IOBJ)),REFLCT)
      QTOUCH=QTOUCH.OR.(ABS(TEMP).GT.RSMALL)
C
      ELSEIF (XSFTYPE(IOBJ).EQ.'INTE') THEN
      CALL XCOPYI(ITEMP,1,HEAP(HPSF(IOBJ)),REFLCT)
      QTOUCH=QTOUCH.OR.(ABS(ITEMP).GT.RSMALL)
      END IF
C
      REFLCT=REFLCT+1
      END DO
      END IF
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XDECGROUP(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF)
C
C Routine groups reciprocal space objects.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER XNAMEMX, XSFMX, XRHOMX
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
C pointer
      INTEGER GRPPTR, GRPOBJ
C
      GRPPTR=ALLHP(INTEG4(XSFMX+1))
      GRPOBJ=ALLHP(INTEG4(XSFMX))
C
      CALL XDECGROU2(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF,
     &   HEAP(GRPPTR), HEAP(GRPOBJ))
C
      CALL FREHP(GRPPTR,INTEG4(XSFMX+1))
      CALL FREHP(GRPOBJ,INTEG4(XSFMX))
C
      RETURN
      END
C================================================================
      SUBROUTINE XDECGROU2(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF,
     &   GRPPTR, GRPOBJ)
C
C Routine groups reciprocal space objects.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'consta.inc'
      INTEGER XNAMEMX, XSFMX, XRHOMX
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*), GRPPTR(*), GRPOBJ(*)
C local
      INTEGER MOBJECT, NOBJECT, LOBJECT, I, II, IOBJECT
      PARAMETER (MOBJECT=10)
      CHARACTER*(WDMAX) SOBJECT(MOBJECT)
      CHARACTER*4 STYPE
      LOGICAL FOUND
      INTEGER IORDER, GRPID
      INTEGER NGROUP, IGROUP, COUNTER
C begin
C
C set defaults
      STYPE=' '
      NOBJECT=0
C
      CALL PUSEND('GROUp>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('GROUp>')
      CALL MISCOM('GROUp>',USED)
      IF (.NOT.USED) THEN
C=====================================================================
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-group')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TYPE') THEN
      CALL NEXTA4('TYPE=',STYPE)
      IF (STYPE.NE.'HL') THEN
      CALL DSPERR('XDECGROUP>','unknown type.')
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'OBJE') THEN
      CALL NEXTST('OBJEct=',WD)
      IF (NOBJECT.LT.MOBJECT) THEN
      NOBJECT=NOBJECT+1
      CALL COPYST(SOBJECT(NOBJECT),XNAMEMX,LOBJECT,WD,WDLEN)
      ELSE
      CALL WRNDIE(-5,'XDECGROUP','Max. number of objects exceeded.')
      END IF
C=====================================================================
      ELSE IF (WD(1:1).EQ.'?') THEN
C
      CALL XGLIST(XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,
     &                  XSFGTYP,XSFGORD,NGROUP,GRPOBJ,
     &                  GRPPTR)
C
      DO IGROUP=1,NGROUP
      WRITE(6,'(A,I6,2A)') ' XDECGROUP: group ',IGROUP,
     &               ' type=',XSFGTYP(GRPOBJ(GRPPTR(IGROUP)+1))
      WRITE(6,'(A)')
     &  '        object name     domain        type     group-order'
      IORDER=0
      DO COUNTER=GRPPTR(IGROUP)+1,GRPPTR(IGROUP+1)
      IORDER=IORDER+1
C
      WRITE(6,'(5A,I6)')
     &  '           ',XSFNAM(GRPOBJ(COUNTER)),
     &  '   RECIprocal    ',XSFTYPE(GRPOBJ(COUNTER)),' ',IORDER
C
      END DO
      WRITE(6,'(A)') ' '
C
      END DO
C=====================================================================
      ELSE
      CALL CHKEND('GROUp>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (STYPE.NE.' '.OR.NOBJECT.GT.0) THEN
      IF (STYPE.EQ.' ') THEN
      WRITE(6,'(A)') '%XDECGROUP-ERR: Group type not specified.'
      CALL WRNDIE(-5,'XUDECLA','Object name not specified.')
      ELSEIF (NOBJECT.EQ.0) THEN
      WRITE(6,'(A)') '%XDECGROUP-ERR: no objects specified.'
      CALL WRNDIE(-5,'XUDECLA','no objects specified.')
      ELSEIF (STYPE.EQ.'HL'.AND.NOBJECT.NE.4) THEN
      WRITE(6,'(2A)')
     & '%XDECGROUP-ERR: four objects need to be specified',
     & ' for type HL.'
      CALL WRNDIE(-5,'XUDECLA','four objects need to be specified.')
      ELSE
C
C get unique group identifier (will be stored in XSFGNAM)
      FOUND=.TRUE.
      GRPID=0
      DO WHILE (FOUND)
      GRPID=GRPID+1
      FOUND=.FALSE.
      DO I=1,XSFNUM
      FOUND=FOUND.OR.XSFGNAM(I).EQ.GRPID
      END DO
      END DO
C
      II=1
      IORDER=0
      DO IOBJECT=1,NOBJECT
      FOUND=.FALSE.
      DO I=1,XSFNUM
      IF (SOBJECT(IOBJECT).EQ.XSFNAM(I)) THEN
      FOUND=.TRUE.
      II=I
      END IF
      END DO
      IF (FOUND.AND.XSFGNAM(II).NE.0) THEN
      WRITE(6,'(3A)') '%XDECGROUP-ERR: Object ',
     &     XSFNAM(II),' already belongs to a group '
      CALL WRNDIE(-5,'XDECGROUP','Object already belongs to group.')
      ELSEIF (FOUND.AND.STYPE.EQ.'HL'.AND.XSFTYPE(II).NE.'REAL') THEN
      WRITE(6,'(3A)') '%XDECGROUP-ERR: Object ',
     &     XSFNAM(II),' must be of type real for group type HL'
      CALL WRNDIE(-5,'XDECGROUP','Object has the wrong type.')
      ELSEIF (FOUND) THEN
      IORDER=IORDER+1
      XSFGNAM(II)=GRPID
      XSFGTYP(II)=STYPE
      XSFGORD(II)=IORDER
      ELSE
      WRITE(6,'(2A)')
     &  ' %XDECGROUP-ERR: reciprocal space object not found ',
     &  SOBJECT(IOBJECT)(1:XNAMEMX)
      CALL WRNDIE(-5,'XDECGROUP','reciprocal space object not found.')
      END IF
C
      END DO
C
      END IF
      END IF
C
      RETURN
      END
C==========================================================================
      SUBROUTINE XDECUNGR(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF)
C
C Routine un-groups reciprocal space objects.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'consta.inc'
      INTEGER XNAMEMX, XSFMX, XRHOMX
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
C local
      INTEGER I, J, II, GRPID
      LOGICAL FOUND, MATCH, QEQUAL
C parameter
C begin
C
C set defaults
      GRPID=0
C
      CALL PUSEND('UNGRoup>')
      DO WHILE (.NOT.DONE)
      CALL NEXTQL('UNGRoup>')
      CALL MISCOM('UNGRoup>',USED)
      IF (.NOT.USED) THEN
C=====================================================================
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-ungroup')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NAME') THEN
      CALL NEXTQL('NAME')
      QEQUAL=.TRUE.
      IF (WD(1:WDLEN).EQ.'=') THEN
      QEQUAL=.TRUE.
      ELSEIF (WD(1:WDLEN).EQ.'#') THEN
      QEQUAL=.FALSE.
      ELSE
      CALL DSPERR('UNGRoup','operator missing (= or #).')
      END IF
      CALL NEXTWD('NAME=')
      FOUND=.FALSE.
      I=1
      DO WHILE (I.LE.XSFNUM)
      CALL EQSTWC(XSFNAM(I),XNAMEMX,WD,WDLEN,1,1,MATCH)
      IF ((QEQUAL.AND.MATCH).OR.(.NOT.QEQUAL.AND..NOT.MATCH)) THEN
      FOUND=.TRUE.
      II=I
C
      IF (XSFGNAM(II).EQ.0) THEN
      WRITE(6,'(3A)')
     &  ' %XDECUNGR-ERR: reciprocal space object ',
     &  WD(1:WDLEN),' does not belong to any group'
      ELSE
C
C mark all objects belonging to this group
      GRPID=XSFGNAM(II)
      DO J=1,XSFNUM
      IF (XSFGNAM(J).EQ.GRPID) THEN
      XSFGNAM(J)=-1
      END IF
      END DO
      END IF
      END IF
C
      I=I+1
      END DO
C
      IF (.NOT.FOUND) THEN
      WRITE(6,'(2A)')
     &  ' XDECUNGR-INFO: no reciprocal space object(s) found matching',
     &  WD(1:WDLEN)
      END IF
C
C
C=====================================================================
      ELSE
      CALL CHKEND('UNGRoup>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (GRPID.NE.0) THEN
C remove group
      DO I=1,XSFNUM
      IF (XSFGNAM(I).EQ.-1) THEN
      XSFGNAM(I)=0
      XSFGTYP(I)=' '
      XSFGORD(I)=0
      END IF
      END DO
      END IF
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XGLIST(XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,
     &                  XSFGTYP,XSFGORD,NGROUP,GRPOBJ,
     &                  GRPPTR)
C
C routine produces list of objects in groups
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER NGROUP, GRPOBJ(*), GRPPTR(*)
C local
      INTEGER COUNTER, I, J, IORDER, IGROUP
      LOGICAL FOUND
C
C begin
C
      NGROUP=0
      COUNTER=0
      GRPPTR(1)=0
C
      DO I=1,XSFNUM
      IF (XSFGNAM(I).NE.0) THEN
C
C check that this is a new group
      FOUND=.FALSE.
      DO IGROUP=1,NGROUP
      IF (XSFGNAM(I).EQ.XSFGNAM(GRPOBJ(GRPPTR(IGROUP+1)))) THEN
      FOUND=.TRUE.
      END IF
      END DO
C
      IF (.NOT.FOUND) THEN
C
      NGROUP=NGROUP+1
C
C found a group
C
C find all objects in group, one at a time.
      IORDER=0
      FOUND=.TRUE.
C
      DO WHILE (FOUND)
      IORDER=IORDER+1
      FOUND=.FALSE.
      DO J=1,XSFNUM
      IF (XSFGORD(J).EQ.IORDER.AND.XSFGNAM(J).EQ.XSFGNAM(I)) THEN
      FOUND=.TRUE.
      COUNTER=COUNTER+1
      GRPOBJ(COUNTER)=J
      END IF
      END DO
      END DO
C
      GRPPTR(NGROUP+1)=COUNTER
C
      END IF
      END IF
      END DO
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XTARMOIN
C
C Routine initialize target, dtarget, and monitor expressions to
C default functions.
C
C Author: Axel T. Brunger
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xtarget.inc'
C local
      INTEGER I, J
C parameter
      DOUBLE PRECISION ONE, ZERO
      PARAMETER (ONE=1.0D0, ZERO=0.0D0)
C begin
C
C default target function expression:
C   use RESI(amplitude(fobs),fcalc,1)
      TRPNN=7
      TDEPTH=3
      DO I=1,TRPNN
      TRPNLEV(I)=1
      TRPNTYP(I)=' '
      TRPNDOM(I)=' '
      TRPNMLT(I)=0
      DO J=1,4
      TRPN(J,I)=' '
      TRPNL(J,I)=0
      TRPNDB(J,I)=DCMPLX(ZERO,ZERO)
      END DO
      END DO
C
      TRPNLEV(1)=1
      TRPN(1,1)='FOBS'
      TRPNL(1,1)=4
      TRPNTYP(1)='DC'
      TRPNDOM(1)='SF'
      TRPNMLT(1)=0
C
      TRPNLEV(2)=1
      TRPN(1,2)='AMPLITUDE'
      TRPNL(1,2)=9
      TRPNTYP(2)='DP'
      TRPNDOM(2)='SF'
      TRPNMLT(2)=1
C
      TRPNLEV(3)=2
      TRPN(1,3)='FCALC'
      TRPNL(1,3)=5
      TRPNTYP(3)='DC'
      TRPNDOM(3)='SF'
      TRPNMLT(3)=0
C
      TRPNLEV(4)=3
      TRPN(1,4)='FPART'
      TRPNL(1,4)=5
      TRPNTYP(4)='DC'
      TRPNDOM(4)='SF'
      TRPNMLT(4)=0
C
      TRPNLEV(5)=2
      TRPN(1,5)='+'
      TRPNL(1,5)=1
      TRPNTYP(5)='DC'
      TRPNDOM(5)='SF'
      TRPNMLT(5)=0
C
      TRPNLEV(6)=3
      TRPN(1,6)='CONS'
      TRPNL(1,6)=4
      TRPNTYP(6)='DP'
      TRPNDOM(6)='SF'
      TRPNMLT(6)=0
      TRPNDB(1,6)=DCMPLX(ONE,ZERO)
C
      TRPNLEV(7)=1
      TRPN(1,7)='RESI'
      TRPNL(1,7)=4
      TRPNTYP(7)='DP'
      TRPNDOM(7)='SF'
      TRPNMLT(7)=3
C default derivative target function expression:
C   use DRESI(amplitude(fobs),fcalc,1)
      DRPNN(1)=7
      DDEPTH(1)=3
      DO I=1,DRPNN(1)
      DRPNLEV(I,1)=1
      DRPNTYP(I,1)=' '
      DRPNDOM(I,1)=' '
      DRPNMLT(I,1)=0
      DO J=1,4
      DRPN(J,I,1)=' '
      DRPNL(J,I,1)=0
      DRPNDB(J,I,1)=DCMPLX(ZERO,ZERO)
      END DO
      END DO
C
      DRPNLEV(1,1)=1
      DRPN(1,1,1)='FOBS'
      DRPNL(1,1,1)=4
      DRPNTYP(1,1)='DC'
      DRPNDOM(1,1)='SF'
      DRPNMLT(1,1)=0
C
      DRPNLEV(2,1)=1
      DRPN(1,2,1)='AMPLITUDE'
      DRPNL(1,2,1)=9
      DRPNTYP(2,1)='DP'
      DRPNDOM(2,1)='SF'
      DRPNMLT(2,1)=1
C
      DRPNLEV(3,1)=2
      DRPN(1,3,1)='FCALC'
      DRPNL(1,3,1)=5
      DRPNTYP(3,1)='DC'
      DRPNDOM(3,1)='SF'
      DRPNMLT(3,1)=0
C
      DRPNLEV(4,1)=3
      DRPN(1,4,1)='FPART'
      DRPNL(1,4,1)=5
      DRPNTYP(4,1)='DC'
      DRPNDOM(4,1)='SF'
      DRPNMLT(4,1)=0
C
      DRPNLEV(5,1)=2
      DRPN(1,5,1)='+'
      DRPNL(1,5,1)=1
      DRPNTYP(5,1)='DC'
      DRPNDOM(5,1)='SF'
      DRPNMLT(5,1)=0
C
      DRPNLEV(6,1)=3
      DRPN(1,6,1)='CONS'
      DRPNL(1,6,1)=4
      DRPNTYP(6,1)='DP'
      DRPNDOM(6,1)='SF'
      DRPNMLT(6,1)=0
      DRPNDB(1,6,1)=DCMPLX(ONE,ZERO)
C
      DRPNLEV(7,1)=1
      DRPN(1,7,1)='DRESI'
      DRPNL(1,7,1)=5
      DRPNTYP(7,1)='DP'
      DRPNDOM(7,1)='SF'
      DRPNMLT(7,1)=3
C
C default monitor function expression:
C   use RVALUE[overall]( fobs, fcalc )
      MRPNN=5
      MDEPTH=3
      DO I=1,MRPNN
      MRPNLEV(I)=1
      MRPNTYP(I)=' '
      MRPNDOM(I)=' '
      MRPNMLT(I)=0
      DO J=1,4
      MRPN(J,I)=' '
      MRPNL(J,I)=0
      MRPNDB(J,I)=DCMPLX(ZERO,ZERO)
      END DO
      END DO
C
      MRPNLEV(1)=1
      MRPN(1,1)='FOBS'
      MRPNL(1,1)=4
      MRPNTYP(1)='DC'
      MRPNDOM(1)='SF'
      MRPNMLT(1)=0
C
      MRPNLEV(2)=2
      MRPN(1,2)='FCALC'
      MRPNL(1,2)=5
      MRPNTYP(2)='DC'
      MRPNDOM(2)='SF'
      MRPNMLT(2)=0
C
      MRPNLEV(3)=3
      MRPN(1,3)='FPART'
      MRPNL(1,3)=5
      MRPNTYP(3)='DC'
      MRPNDOM(3)='SF'
      MRPNMLT(3)=0
C
      MRPNLEV(4)=2
      MRPN(1,4)='+'
      MRPNL(1,4)=1
      MRPNTYP(4)='DC'
      MRPNDOM(4)='SF'
      MRPNMLT(4)=0
      MRPNLEV(4)=1
C
      MRPN(1,5)='RVALUE'
      MRPNL(1,5)=6
      MRPNTYP(5)='DP'
      MRPNDOM(5)='SF'
      MRPNMLT(5)=2
      MRPNDB(1,5)=DCMPLX(-ONE,ZERO)
      MRPNDB(2,5)=DCMPLX(-ONE,ZERO)
C
      RETURN
      END
C
