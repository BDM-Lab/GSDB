      SUBROUTINE READDP
C
C command parser for dynamics reader
C
C For syntax see main parsing loop "HELP"
C
C variables between subsequent calls are saved.
C
C By Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'trtraj.inc'
C local
      INTEGER MXLIST
      PARAMETER (MXLIST=30)
      LOGICAL START, QFORM
      INTEGER NUNIT, FLIST(MXLIST), UNIT, BEGIN, STOP, SKIP
      LOGICAL TDONE, TEOF
      INTEGER ISTEP
      DOUBLE PRECISION DELTA
      CHARACTER*4 HDR
      DOUBLE PRECISION DPPLAC
      DOUBLE COMPLEX DCPLAC
C
      SAVE START, QFORM
      SAVE NUNIT, FLIST, UNIT, BEGIN, STOP, SKIP
      SAVE TDONE, TEOF
      SAVE ISTEP
      SAVE DELTA
      SAVE HDR
C begin
      CALL PUSEND('READ-TRAJectory>')
      CALL NEXTWD('READ-TRAJectory>')
      IF (WD(1:4).EQ.'NEXT') THEN
      IF (.NOT.TDONE) THEN
      IF (TRIZ.NE.0) THEN
C
C subsequent call
      CALL READDC(NATOM,X,Y,Z,
     &           START,TDONE,TEOF,ERROR,
     &           NUNIT,FLIST,UNIT,BEGIN,STOP,SKIP,
     &           ISTEP,DELTA,HDR,QFORM)
      WRITE(6,'(A,I8,A)') ' READCP: trajectory number ',ISTEP,
     &              ' has been copied to the main coordinate set'
      IF (TDONE) THEN
         CALL DECLAR( 'STATUS', 'ST', 'COMPLETE', DCPLAC, DPPLAC )
      ENDIF
      ELSE
      WRITE(6,'(A)')' %READCP-ERR: coordinate database was deleted.'
      END IF
      ELSE
      WRITE(6,'(A)')' %READCP-ERR: EOF reached - coordinates unchanged'
      END IF
      CALL NEXTWD('READ-TRAJectory>')
      CALL CHKEND('READ-TRAJectory>',DONE)
      DONE=.FALSE.
      ELSE
      CALL SAVEWD
C
C initial call
C
C defaults
      CALL DYNSET('INIT',USED,BEGIN,SKIP,STOP,NUNIT,
     &                  MXLIST,FLIST,QFORM)
C
      CALL DECLAR( 'STATUS', 'ST', 'READ', DCPLAC, DPPLAC )
      START=.TRUE.
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('READ-TRAJectory>')
C
      CALL DYNSET('PARSE',USED,BEGIN,SKIP,STOP,NUNIT,
     &                  MXLIST,FLIST,QFORM)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-read-trajectory')
C
      ELSE
      CALL CHKEND('READ-TRAJectory',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
      CALL READDC(NATOM,X,Y,Z,
     &           START,TDONE,TEOF,ERROR,
     &           NUNIT,FLIST,UNIT,BEGIN,STOP,SKIP,
     &           ISTEP,DELTA,HDR,QFORM)
      IF (TDONE) THEN
         CALL DECLAR( 'STATUS', 'ST', 'COMPLETE', DCPLAC, DPPLAC )
      ENDIF
      WRITE(6,'(A,I8,A)') ' READ-TRAJectory: frame # ',ISTEP,
     &              ' has been copied to the main coordinate set.'
      END IF
      RETURN
      END
C====================================================================
      SUBROUTINE DYNSET(MODE,LUSED,BEGIN,SKIP,STOP,NUNIT,
     &                  MXLIST,FLIST,QFORM)
C
C parses the dynamics I/O specifications
C For syntax see main parsing loop "HELP"
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      CHARACTER*(*) MODE
      LOGICAL LUSED
      INTEGER BEGIN, SKIP, STOP, NUNIT, MXLIST, FLIST(*)
      LOGICAL QFORM
C local
      INTEGER UNIT
C begin
      IF (MODE.EQ.'INIT') THEN
      BEGIN=1
      SKIP=1
      STOP=999999
      QFORM=.FALSE.
      NUNIT=0
C
      ELSE IF (MODE.EQ.'PARSE') THEN
      LUSED=.TRUE.
C
      IF (WD(1:4).EQ.'BEGI') THEN
      CALL NEXTI('BEGIn=',BEGIN)
      ELSE IF (WD(1:4).EQ.'SKIP') THEN
      CALL NEXTI('SKIP=',SKIP)
      ELSE IF (WD(1:4).EQ.'STOP') THEN
      CALL NEXTI('STOP=',STOP)
      ELSE IF (WD(1:4).EQ.'ASCI') THEN
      CALL NEXTLO('ASCIi=',QFORM)
      ELSE IF (WD(1:4).EQ.'INPU') THEN
      CALL NEXTFI('INPUt=',IFILE)
      IF (QFORM) THEN
      CALL ASSFIL(IFILE,UNIT,'READ','FORMATTED',ERROR)
      ELSE
      CALL ASSFIL(IFILE,UNIT,'READ','UNFORMATTED',ERROR)
      END IF
      IF (.NOT.ERROR) THEN
      IF (NUNIT.GT.MXLIST) THEN
      WRITE(6,'(A)') ' %DYNSET-ERR: MXLIST=(maximum number of'//
     &  ' trajectory files) exceeded.'
      ELSE
      NUNIT=NUNIT+1
      END IF
      FLIST(NUNIT)=UNIT
      END IF
      ELSE
      LUSED=.FALSE.
      END IF
      END IF
C
      RETURN
      END
C====================================================================
      SUBROUTINE READDC(NATOM,X,Y,Z,
     &                  START,DONE,EOF,ERROR,
     &                  NUNIT,FLIST,UNIT,BEGIN,STOP,SKIP,
     &                  ISTEP,DELTA,HDR,QFORM)
C
C "Intelligent" trajectory reader;
C Routine reads dynamics trajectories written by WRITEC;
C if EOF encoutered tries to read next file given by list
C FLIST(1...NUNIT).
C features automatic correction of mismatches between files (if
C possible). Mismatches can occur when restarts and trajectories
C are not synchronized (e.g. after system crash's).
C Reading error and end-of-file handling.
C Condition for passing a new set to calling routine:
C      MOD(ISTEP-BEGIN,SKIP)=0 and BEGIN<=ISTEP<=STOP
C
C Variables:
C
C START (I/O) if .TRUE. start reading information on UNIT
C             if .FALSE. subsequent call to routine
C DONE (on exit) if .TRUE. reading done
C EOF (on exit) encountered EOF before reading done
C ERROR (on exit) encountered ERROR before reading done
C NATOM (on input) number of atoms (for checking purposes)
C NUNIT number of units
C FLIST(*) list of units
C UNIT (I/O) current unit number
C BEGIN (I/O) number of first coordinate set; modified if incompatible
C SKIP (I/O) number to be skipped between sets;  - " -       - " -
C STOP (I/O) number of last coordinate set;      - " -       - " -
C
C ISTEP (on exit) current step number
C DELTA (on exit) time-step size between steps
C                  => current time = ISTEP * DELTA
C HDR (on exit) header
C QFORM: logical flag indicating whether to read formatted or unformatted file
C
C X(*), Y(*), Z(*) (on exit) coordinate set at step ISTEP
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'trtraj.inc'
      INTEGER NATOM
      DOUBLE PRECISION  X(*), Y(*), Z(*)
      LOGICAL START,DONE,EOF,ERROR
      INTEGER NUNIT, FLIST(*)
      INTEGER UNIT, BEGIN, STOP, SKIP
      INTEGER ISTEP
      DOUBLE PRECISION    DELTA
      CHARACTER*4 HDR
      LOGICAL QFORM
C local
      INTEGER NFREAT
      INTEGER NSAVC, NSAVCX, NFIXEX, PASSED
      DOUBLE PRECISION DELTAX
      LOGICAL QFIRST, ENSEMB
      DOUBLE PRECISION SCALE, OFFSET
      CHARACTER*10 FORM
      CHARACTER*4 HDRX
C
      SAVE NFREAT
      SAVE NSAVC, NSAVCX, NFIXEX, PASSED
      SAVE DELTAX
      SAVE QFIRST, ENSEMB
      SAVE SCALE, OFFSET
      SAVE FORM
      SAVE HDRX
C
C begin
      IF (START) THEN
      IF (TRIZ.NE.0) CALL FREHP(TRIZ,INTEG4(TRNAT))
      IF (TRIY.NE.0) CALL FREHP(TRIY,INTEG4(TRNAT))
      IF (TRIX.NE.0) CALL FREHP(TRIX,INTEG4(TRNAT))
      IF (TRTEMP.NE.0) CALL FREHP(TRTEMP,IREAL4(TRNAT))
      IF (TRFREEAT.NE.0) CALL FREHP(TRFREEAT,INTEG4(TRNAT))
      TRFREEAT=ALLHP(INTEG4(NATOM))
      TRTEMP=ALLHP(IREAL4(NATOM))
      TRIX=ALLHP(INTEG4(NATOM))
      TRIY=ALLHP(INTEG4(NATOM))
      TRIZ=ALLHP(INTEG4(NATOM))
      TRNAT=NATOM
      END IF
      CALL READC2(NATOM,X,Y,Z,
     &            START,DONE,EOF,ERROR,
     &            NUNIT,FLIST,UNIT,BEGIN,STOP,SKIP,
     &            ISTEP,DELTA,HDR,HEAP(TRTEMP),
     &            HEAP(TRFREEAT),NFREAT,
     &            NSAVC,NSAVCX,NFIXEX,DELTAX,HDRX,PASSED,QFIRST,
     &            ENSEMB,QFORM,FORM,SCALE,OFFSET,HEAP(TRIX),
     &            HEAP(TRIY),HEAP(TRIZ))
      IF (DONE) THEN
      CALL FREHP(TRIZ,INTEG4(NATOM))
      CALL FREHP(TRIY,INTEG4(NATOM))
      CALL FREHP(TRIX,INTEG4(NATOM))
      CALL FREHP(TRTEMP,IREAL4(NATOM))
      CALL FREHP(TRFREEAT,INTEG4(NATOM))
      TRIZ=0
      TRIY=0
      TRIX=0
      TRTEMP=0
      TRFREEAT=0
      END IF
      RETURN
      END
C==================================================================
      SUBROUTINE READC2(NATOM,X,Y,Z,
     &                  START,DONE,EOF,ERROR,
     &                  NUNIT,FLIST,UNIT,BEGIN,STOP,SKIP,
     &                  ISTEP,DELTA,HDR,TEMP,FREEAT,NFREAT,
     &                  NSAVC,NSAVCX,NFIXEX,DELTAX,HDRX,PASSED,QFIRST,
     &                  ENSEMB,QFORM,FORM,SCALE,OFFSET,IX,IY,IZ)
C see READC above
C
      IMPLICIT NONE
C input/output
      INCLUDE 'consta.inc'
      INTEGER NATOM
      DOUBLE PRECISION X(*), Y(*), Z(*)
      LOGICAL START, DONE, EOF, ERROR
      INTEGER NUNIT, FLIST(*), UNIT, BEGIN, STOP, SKIP
      INTEGER ISTEP
      DOUBLE PRECISION DELTA
      REAL TEMP(*)
      CHARACTER*4 HDR
      INTEGER NFREAT, FREEAT(*)
      INTEGER NSAVC, NSAVCX, NFIXEX, PASSED
      DOUBLE PRECISION DELTAX
      CHARACTER*4 HDRX
      LOGICAL QFIRST, ENSEMB, QFORM
      CHARACTER*10 FORM
      DOUBLE PRECISION SCALE,OFFSET
      INTEGER IX(*), IY(*), IZ(*)
C local
      INTEGER NFIXED, ISTART, NATOMX
      LOGICAL NEXT, COND
C
C begin
      IF (NUNIT.LE.0) THEN
      WRITE(6,'(A)') ' %READC-ERR: NUNIt=0'
      ERROR=.TRUE.
      DONE=.TRUE.
      ELSE
C
C the following test is a kluge - I did not want to modify the
C calling sequence to include a new passed variable.
      IF (.NOT.START) ENSEMB=HDR.EQ.'ENSE'
C
      IF (START) THEN
      PASSED=0
      QFIRST=.FALSE.
      DONE=.FALSE.
      EOF=.FALSE.
      ERROR=.FALSE.
      ENSEMB=.FALSE.
      UNIT=1
      REWIND(UNIT=FLIST(UNIT))
      CALL READCH(FLIST,UNIT,NFREAT,FREEAT,ISTART,NSAVC,NFIXED,
     &           NATOM,NSAVCX,NFIXEX,NATOMX,DELTA,DELTAX,
     &           HDR,HDRX,START,QFIRST,EOF,ERROR,QFORM,FORM,
     &           SCALE,OFFSET)
      ISTEP=ISTART-NSAVC
      START=.FALSE.
      IF (.NOT. (ERROR .OR. EOF)) THEN
C
C check BEGIN, STOP, SKIP parameters
      IF (MOD(BEGIN-ISTART,NSAVC).NE.0.OR.BEGIN.LT.ISTART) THEN
      BEGIN= MAX(BEGIN,ISTART)
      BEGIN= ISTART+((BEGIN-ISTART)/NSAVC) * NSAVC
      WRITE(6,'(A)') ' READC: BEGIN not compatible. Modified'
      END IF
      IF (MOD(SKIP,NSAVC).NE.0.OR.SKIP.LE.0) THEN
      SKIP=MAX(NSAVC,SKIP)
      SKIP= (SKIP/NSAVC) * NSAVC
      WRITE(6,'(A)') ' READC: SKIP not compatible. Modified'
      END IF
      IF (MOD(STOP-BEGIN,SKIP).NE.0.OR.STOP.LT.BEGIN) THEN
      STOP=MAX(BEGIN,STOP)
      STOP=BEGIN+((STOP-BEGIN)/SKIP) * SKIP
      WRITE(6,'(A)') ' READC: STOP not compatible. Modified'
      END IF
      WRITE(6,9010)
     1' READC: NUNIt=',NUNIT,
     2' READC: BEGIn=',BEGIN,'  SKIP=',SKIP,'  STOP=',STOP,
     3' READC: TIMEstep=',DELTA*TIMFAC,' ps    header= ',HDR
9010  FORMAT(A,I4,/,A,I8,A,I8,A,I8,/,A,F12.8,A,A)
      END IF
      END IF
C
      NEXT = .FALSE.
      DO WHILE (.NOT. (NEXT .OR. EOF .OR. ERROR))
      COND = .FALSE.
      DO WHILE (.NOT. (COND .OR. EOF .OR. ERROR))
      ISTEP=ISTEP+NSAVC
      COND=ISTEP.GE.BEGIN.AND.MOD(ISTEP-BEGIN,SKIP).EQ.0
      CALL READCN(FLIST,UNIT,NATOM,NFREAT,ISTEP,FREEAT,
     &                  X,Y,Z,TEMP,QFIRST,EOF,ERROR,QFORM,FORM,SCALE,
     &                  OFFSET,IX,IY,IZ)
      END DO
      NEXT=(.NOT.EOF.AND..NOT.ERROR)
      DO WHILE (EOF .AND. .NOT.ERROR .AND. UNIT.LT.NUNIT)
C open next file and read header
      ISTEP=ISTEP-NSAVC
      UNIT=UNIT+1
      EOF=.FALSE.
      REWIND(UNIT=FLIST(UNIT))
      CALL READCH(FLIST,UNIT,NFREAT,FREEAT,ISTART,NSAVC,NFIXED,
     &           NATOM,NSAVCX,NFIXEX,NATOMX,DELTA,DELTAX,
     &           HDR,HDRX,START,QFIRST,EOF,ERROR,QFORM,FORM,SCALE,
     &           OFFSET)
      IF (.NOT.EOF.AND..NOT.ERROR) THEN
      IF (ISTART-ISTEP.LT.NSAVC.AND..NOT.ENSEMB) THEN
C do correction: read until ISTEP is reached
      WRITE(6,'(A,I4)') ' READC: mismatch correction; file-number=',UNIT
      DO WHILE (ISTEP.GT.(ISTART-NSAVC) .AND. .NOT. (EOF.OR.ERROR))
      ISTART=ISTART+NSAVC
      CALL READCN(FLIST,UNIT,NATOM,NFREAT,ISTEP,FREEAT,
     &                  X,Y,Z,TEMP,QFIRST,EOF,ERROR,QFORM,FORM,
     &                  SCALE,OFFSET,IX,IY,IZ)
      END DO
      END IF
C
      IF (ISTART-ISTEP.NE.NSAVC.AND..NOT.ENSEMB) THEN
      WRITE(6,'(A,I4)') ' READC: mismatch error; file-number=',UNIT
      WRITE(6,'(A,I8)')
     & ' READC: missing coordinate sets skipped.  Now at step',ISTART
      ISTEP=ISTART-NSAVC
      END IF
      END IF
      END DO
C
      END DO
C
      IF (NEXT) PASSED=PASSED+1
      IF(EOF)WRITE(6,9000)  ' %READC-EOF file ',UNIT,' at step',ISTEP
      IF(ERROR)WRITE(6,9000)' %READC-ERROR file ',UNIT,' at step',ISTEP
      DONE=(ISTEP.GE.STOP .OR. EOF .OR. ERROR)
      IF (DONE) THEN
      WRITE(6,'(A,I8,A)') ' READC: complete. ',PASSED,
     1              ' coordinate sets were passed to calling routine'
      END IF
      END IF
9000  FORMAT(A,I4,A,I8)
C
      RETURN
      END
C
      SUBROUTINE READCH(FLIST,UNIT,NFREAT,FREEAT,ISTART,NSAVC,NFIXED,
     &           NATOM,NSAVCX,NFIXEX,NATOMX,DELTA,DELTAX,
     &           HDR,HDRX,START,QFIRST,EOF,ERROR,QFORM,FORM,SCALE,
     &           OFFSET)
C
C reads header from trajectory file UNIT
C (see routine READC2 above)
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER FLIST(*), UNIT, NFREAT, FREEAT(*)
      INTEGER ISTART, NSAVC, NFIXED, NATOM, NSAVCX, NFIXEX, NATOMX
      DOUBLE PRECISION DELTA, DELTAX
      CHARACTER*4 HDR, HDRX
      LOGICAL START, QFIRST, EOF, ERROR
      LOGICAL QFORM
      CHARACTER*10 FORM
      DOUBLE PRECISION SCALE, OFFSET
C local
      INTEGER I
C begin
      QFIRST=.TRUE.
C
      IF (QFORM) THEN
      CALL RDTITB(FLIST(UNIT),'FORMATTED')
      READ(FLIST(UNIT),'(1X,A,1X,4I7,G14.6)',ERR=999,END=888)
     &     HDR, ISTART, NSAVC, I, NATOMX, DELTA
      READ(FLIST(UNIT),'(G14.6,G14.6,A)',ERR=999,END=888)
     &     SCALE, OFFSET, FORM
      READ(FLIST(UNIT),'(10I7)',ERR=999,END=888)
     &      NFREAT, (FREEAT(I),I=1,NFREAT)
      READ(FLIST(UNIT),'(A)',ERR=999,END=888) HDR
C
      NFIXED=NATOMX-NFREAT
C
      ELSE
      READ(FLIST(UNIT),ERR=999,END=888) HDR,I,ISTART,NSAVC,I,I,
     &                            I,I,I,NFIXED,DELTA,
     &                            I,I,I,I,I,I,I,I,I
      CALL RDTITB(FLIST(UNIT),'BINARY')
      READ(FLIST(UNIT)) NATOMX
      NFREAT=NATOM-NFIXED
      IF (NFREAT.NE.NATOMX) THEN
      READ(FLIST(UNIT),ERR=999,END=888) (FREEAT(I),I=1,NFREAT)
      ELSE
      DO I=1,NATOMX
      FREEAT(I)=I
      END DO
      END IF
      END IF
C
      IF (NATOM.NE.NATOMX) THEN
      WRITE(6,'(A)') ' READC: number of atoms do not match'
      write(6,*) natom, natomx
      ERROR=.TRUE.
      END IF
      WRITE(6,'(A,I4,A,I8)')
     & ' READC: reading file=',UNIT,'   ISTART=',ISTART
      IF (.NOT.START) THEN
      ERROR=(NSAVC.NE.NSAVCX.OR.NFIXED.NE.NFIXEX
     &       .OR.ABS(DELTA-DELTAX).GT.R4SMAL.OR.HDR.NE.HDRX)
      IF(ERROR) WRITE(6,'(A,I4)')' READC: header mismatch on file',UNIT
      END IF
      NSAVCX=NSAVC
      NFIXEX=NFIXED
      DELTAX=DELTA
      HDRX=HDR
C
      WRITE(6,'(A,I8,A)') ' READC: ',NFREAT,' atoms were "free"'
      GOTO 777
999   ERROR=.TRUE.
      GOTO 777
888   EOF=.TRUE.
777   CONTINUE
      RETURN
      END
C
      SUBROUTINE READCN(FLIST,UNIT,NATOM,NFREAT,ISTEP,FREEAT,
     &                  X,Y,Z,TEMP,QFIRST,EOF,ERROR,QFORM,FORM,SCALE,
     &                  OFFSET,IX,IY,IZ)
C
C read next coordinate set from trajectory file UNIT
C (see routine READC2 above)
      IMPLICIT NONE
C I/O
      INTEGER FLIST(*), UNIT, NATOM, NFREAT, ISTEP, FREEAT(*)
      DOUBLE PRECISION X(*), Y(*), Z(*)
      REAL TEMP(*)
      LOGICAL QFIRST, EOF, ERROR
      LOGICAL QFORM
      CHARACTER*10 FORM
      DOUBLE PRECISION SCALE, OFFSET
      INTEGER IX(*), IY(*), IZ(*)
C local
      INTEGER I
      CHARACTER*4 HDR
C begin
      IF (QFIRST) THEN
      QFIRST=.FALSE.
      IF (QFORM) THEN
      READ(FLIST(UNIT),'('//FORM//')',END=666,ERR=555)
     &      (IX(I),IY(I),IZ(I), I=1,NATOM)
      DO I=1,NATOM
      X(I)=IX(I)/SCALE-OFFSET
      Y(I)=IY(I)/SCALE-OFFSET
      Z(I)=IZ(I)/SCALE-OFFSET
      END DO
      ELSE
      READ(FLIST(UNIT),END=666,ERR=555) (TEMP(I),I=1,NATOM)
      DO I=1,NATOM
      X(I)=TEMP(I)
      END DO
      READ(FLIST(UNIT),END=666,ERR=555) (TEMP(I),I=1,NATOM)
      DO I=1,NATOM
      Y(I)=TEMP(I)
      END DO
      READ(FLIST(UNIT),END=666,ERR=555) (TEMP(I),I=1,NATOM)
      DO I=1,NATOM
      Z(I)=TEMP(I)
      END DO
      END IF
      ELSE
      IF (QFORM) THEN
      READ(FLIST(UNIT),'(A)',END=666,ERR=555) HDR
      READ(FLIST(UNIT),'('//FORM//')',END=666,ERR=555)
     &     (IX(I),IY(I),IZ(I), I=1,NFREAT)
      DO I=1,NFREAT
      X(FREEAT(I))=IX(I)/SCALE-OFFSET
      Y(FREEAT(I))=IY(I)/SCALE-OFFSET
      Z(FREEAT(I))=IZ(I)/SCALE-OFFSET
      END DO
      ELSE
      READ(FLIST(UNIT),END=666,ERR=555) (TEMP(I),I=1,NFREAT)
      DO I=1,NFREAT
      X(FREEAT(I))=TEMP(I)
      END DO
      READ(FLIST(UNIT),END=666,ERR=555) (TEMP(I),I=1,NFREAT)
      DO I=1,NFREAT
      Y(FREEAT(I))=TEMP(I)
      END DO
      READ(FLIST(UNIT),END=666,ERR=555) (TEMP(I),I=1,NFREAT)
      DO I=1,NFREAT
      Z(FREEAT(I))=TEMP(I)
      END DO
      END IF
      END IF
C
      GOTO 444
666   EOF=.TRUE.
      GOTO 444
555   ERROR=.TRUE.
444   CONTINUE
      IF (EOF) THEN
      WRITE(6,'(A,I4,A,I8)')
     &' READC: EOF on file=',UNIT,'  while trying to read step=',ISTEP
      END IF
      RETURN
      END
C====================================================================
      SUBROUTINE WRITTC(HDR,NATOM,X,Y,Z,
     &                  NFREAT,FREEAT,ISTART,QFIRST,
     &                  DELTA,NSAVC,UNIT,
     &                  QFORM,FORM,SCALE,OFFSET)
C
C Writes a set of coordinates/velocities for a single dynamics step;
C if NFREAT unequal NATOM the output is compressed using FREEAT;
C ISTART should be the number of the first set written on file UNIT;
C if the flag QFIRST is true a new file with a header is created;
C NSAVC is the coordinate set writing period;
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      CHARACTER*4 HDR
      INTEGER NATOM
      DOUBLE PRECISION  X(*), Y(*), Z(*)
      INTEGER FREEAT(*), NFREAT, ISTART
      LOGICAL QFIRST
      DOUBLE PRECISION    DELTA
      INTEGER NSAVC
      INTEGER UNIT
      LOGICAL QFORM
      CHARACTER*10 FORM
      DOUBLE PRECISION SCALE, OFFSET
C pointer
      INTEGER IX, IY, IZ
C local
      IX=ALLHP(INTEG4(NATOM))
      IY=ALLHP(INTEG4(NATOM))
      IZ=ALLHP(INTEG4(NATOM))
      CALL WRITT2(HDR,NATOM,X,Y,Z,
     &                  NFREAT,FREEAT,ISTART,QFIRST,
     &                  DELTA,NSAVC,UNIT,
     &                  QFORM,FORM,SCALE,OFFSET,HEAP(IX),HEAP(IY),
     &                  HEAP(IZ))
      CALL FREHP(IZ,INTEG4(NATOM))
      CALL FREHP(IY,INTEG4(NATOM))
      CALL FREHP(IX,INTEG4(NATOM))
      RETURN
      END
C====================================================================
      SUBROUTINE WRITT2(HDR,NATOM,X,Y,Z,
     &                  NFREAT,FREEAT,ISTART,QFIRST,
     &                  DELTA,NSAVC,UNIT,
     &                  QFORM,FORM,SCALE,OFFSET,IX,IY,IZ)
C
C See routine WRITTC above
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      CHARACTER*4 HDR
      INTEGER NATOM
      DOUBLE PRECISION  X(*), Y(*), Z(*)
      INTEGER FREEAT(*), NFREAT, ISTART
      LOGICAL QFIRST
      DOUBLE PRECISION    DELTA
      INTEGER NSAVC
      INTEGER UNIT
      LOGICAL QFORM
      CHARACTER*10 FORM
      DOUBLE PRECISION SCALE, OFFSET
      INTEGER IX(*), IY(*), IZ(*)
C local
      INTEGER I0, I
      PARAMETER (I0=0)
C begin
      IF (UNIT.LT.0) CALL WRNDIE(-3,'WRITEV','invalid logical unit')
      IF (QFIRST) THEN
      WRITE(6,'(3A,I8,A,I5,A)')
     & ' WRITEC: ',HDR,' sets starting from step ',ISTART,
     & ' are written every ',NSAVC,'  steps'
      WRITE(6,'(A,I8,A)') ' WRITEC: where ',NFREAT,' atoms are "free"'
C
      IF (QFORM) THEN
      CALL WRTITL(UNIT,2)
      WRITE(UNIT,'(1X,A,1X,4I7,G14.6)')
     &      HDR, ISTART, NSAVC, I0, NATOM, DELTA
      WRITE(UNIT,'(G14.6,G14.6,A)') SCALE, OFFSET, FORM
      WRITE(UNIT,'(10I7)') NFREAT, (FREEAT(I),I=1,NFREAT)
      WRITE(UNIT,'(A)') HDR
      ELSE
      WRITE(UNIT) HDR,I0,ISTART,NSAVC,I0,I0,
     &              I0,I0,I0,NATOM-NFREAT,DELTA,
     &              I0,I0,I0,I0,I0,I0,I0,I0,I0
      CALL WRTITL(UNIT,-1)
      WRITE(UNIT) NATOM
      IF (NFREAT.NE.NATOM) WRITE(UNIT) (FREEAT(I),I=1,NFREAT)
      END IF
      END IF
C
      IF (QFIRST) THEN
      IF (QFORM) THEN
      DO I=1,NATOM
      IX(I)=INT((X(I)+OFFSET)*SCALE)
      IY(I)=INT((Y(I)+OFFSET)*SCALE)
      IZ(I)=INT((Z(I)+OFFSET)*SCALE)
      END DO
      WRITE(UNIT,'('//FORM//')') (IX(I),IY(I),IZ(I),I=1,NATOM)
      ELSE
C write coordinates in single precision
      WRITE (UNIT) (REAL(X(I)),I=1,NATOM)
      WRITE (UNIT) (REAL(Y(I)),I=1,NATOM)
      WRITE (UNIT) (REAL(Z(I)),I=1,NATOM)
      END IF
      ELSE
      IF (QFORM) THEN
      DO I=1,NFREAT
      IX(I)=INT((X(FREEAT(I))+OFFSET)*SCALE)
      IY(I)=INT((Y(FREEAT(I))+OFFSET)*SCALE)
      IZ(I)=INT((Z(FREEAT(I))+OFFSET)*SCALE)
      END DO
      WRITE(UNIT,'(A)') HDR
      WRITE(UNIT,'('//FORM//')') (IX(I),IY(I),IZ(I),I=1,NFREAT)
      ELSE
      WRITE (UNIT) (REAL(X(FREEAT(I))),I=1,NFREAT)
      WRITE (UNIT) (REAL(Y(FREEAT(I))),I=1,NFREAT)
      WRITE (UNIT) (REAL(Z(FREEAT(I))),I=1,NFREAT)
      END IF
      END IF
C
      QFIRST=.FALSE.
C
      RETURN
      END
C
      SUBROUTINE TRJMERGE
C
C reads in specified sections of a trajectory on files
C FILEs and writes them onto the single file OUTPUT;
C optionally a orientation with respect to the main coordinate set is
C applied.
C
C For syntax see main parsing loop "HELP"
C
C By Axel T. Brunger
C
      IMPLICIT NONE
C input/ouput
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
C local
      INTEGER MXLIST
      PARAMETER (MXLIST=20)
      INTEGER FLIST(MXLIST), UNIT
      INTEGER XX, YY, ZZ, FREEAT, NFREAT, FLAGS, PAIR
      INTEGER PMASS, NPAIR, NSELCT
      INTEGER NUNIT, OUTPUT, SKIP, BEGIN, STOP, ISTEP
      DOUBLE PRECISION    DELTA
      LOGICAL QFIRST, START, ORIENT, LMASS, LWEIG
      LOGICAL LNORO, TDONE, TEOF, ENSEMB, QFORM, OQFORM
      CHARACTER*10 OFORM
      DOUBLE PRECISION OSCALE, OOFFST
      CHARACTER*4 HDR
C begin
C
C defaults
C
      CALL DYNSET('INIT',USED,BEGIN,SKIP,STOP,NUNIT,
     &                  MXLIST,FLIST,QFORM)
      ORIENT=.FALSE.
      ENSEMB=.FALSE.
      OFILE='OUTPUT'
      OQFORM=QFORM
      OOFFST=800.0D0
      OSCALE=10000.0D0
      OFORM='12Z6'
C
C parsing
      CALL PUSEND('MERGE>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('MERGE>')
      CALL DYNSET('PARSE',USED,BEGIN,SKIP,STOP,NUNIT,
     &                  MXLIST,FLIST,QFORM)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-dynamics-merge')
C
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut-file=',OFILE)
      ELSE IF (WD(1:4).EQ.'OASC') THEN
      CALL NEXTLO('OASCii=',OQFORM)
      ELSE IF (WD(1:4).EQ.'FORM') THEN
      CALL NEXTST('OFORmat=',OFORM)
      ELSE IF (WD(1:4).EQ.'SCAL') THEN
      CALL NEXTF('SCALe=',OSCALE)
      ELSE IF (WD(1:4).EQ.'OFFS') THEN
      CALL NEXTF('OFFSet=',OOFFST)
      ELSE IF (WD(1:4).EQ.'ENSE') THEN
      CALL NEXTLO('ENSEmble=',ENSEMB)
      ELSE IF (WD(1:4).EQ.'ORIE') THEN
      FLAGS=ALLHP(INTEG4(NATOM))
      PAIR=ALLHP(INTEG4(2*NATOM))
      PMASS=ALLHP(IREAL8(NATOM))
C
C defaults
      ORIENT=.TRUE.
      LMASS=.FALSE.
      LNORO=.FALSE.
      LWEIG=.FALSE.
      CALL FILL4(HEAP(FLAGS),NATOM,1)
      NSELCT=NATOM
C
C parsing
      CALL PUSEND('ORIENT>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('ORIENT>')
      IF (WD(1:4).EQ.'HELP') THEN
      CALL CNSHELP('cns-dynamics-merge-orient')
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(HEAP(FLAGS),NSELCT,X,Y,Z,.TRUE.)
      ELSE IF (WD(1:4).EQ.'MASS') THEN
      CALL NEXTLO('MASS=',LMASS)
      ELSE IF (WD(1:4).EQ.'NORO') THEN
      CALL NEXTLO('NOROtation=',LNORO)
      ELSE IF (WD(1:4).EQ.'WEIG') THEN
      CALL NEXTLO('WEIGthing=',LWEIG)
      ELSE
      CALL CHKEND('ORIENT>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      ELSE
      CALL CHKEND('MERGE>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (OQFORM) THEN
      CALL ASSFIL(OFILE,OUTPUT,'WRITE','FORMATTED',ERROR)
      ELSE
      CALL ASSFIL(OFILE,OUTPUT,'WRITE','UNFORMATTED',ERROR)
      END IF
C
      IF (ENSEMB) THEN
      WRITE(6,'(A)')
     & ' MERGE: ensemble merging of several trajectory files'
      END IF
      IF (ORIENT.AND.NSELCT.GT.0) THEN
      WRITE(6,'(A,A,/,A,L1,A,L1,A,L1)')
     &' MERGE: output-trajectory oriented with respect to main ',
     &         'coordinate set',
     &' MERGE: MASSweighting=',LMASS,'  NO-ROtation=',LNORO,
     &'  wmainWEIGhting=',LWEIG
      CALL MKPAIR(NATOM,HEAP(FLAGS),HEAP(PAIR),NPAIR,X,Y,Z)
      CALL PKMASS(AMASS,AMASS,HEAP(PMASS),HEAP(PAIR),NPAIR,LMASS,
     &            LWEIG,WMAIN)
      END IF
C
      XX=ALLHP(IREAL8(NATOM))
      YY=ALLHP(IREAL8(NATOM))
      ZZ=ALLHP(IREAL8(NATOM))
      FREEAT=ALLHP(INTEG4(NATOM))
C set up the free atom array for output
      CALL MKFRAT(IMOVE,NATOM,HEAP(FREEAT),NFREAT)
C
      QFIRST=.TRUE.
      START=.TRUE.
      TDONE=.FALSE.
      TEOF=.FALSE.
      ERROR=.FALSE.
      DO WHILE (.NOT. (TDONE.OR.TEOF.OR.ERROR))
CCC modification ATB 4/27/08
      IF (ENSEMB) THEN
      HDR='ENSE'
      ELSE
      HDR='CORD'
      END IF
      CALL READDC(NATOM,HEAP(XX),HEAP(YY),HEAP(ZZ),
     &           START,TDONE,TEOF,ERROR,
     &           NUNIT,FLIST,UNIT,BEGIN,STOP,SKIP,
     &           ISTEP,DELTA,HDR,QFORM)
      IF (ORIENT) THEN
C do orientation with respect to the main coordinate set
      CALL ROTLS1(X,Y,Z,NATOM,HEAP(XX),HEAP(YY),HEAP(ZZ),
     &     NATOM,HEAP(PAIR),NPAIR,HEAP(PMASS),.FALSE.,LNORO)
      END IF
C
      IF (.NOT.(TEOF.OR.ERROR)) THEN
      CALL WRITTC(HDR,NATOM,HEAP(XX),HEAP(YY),HEAP(ZZ),
     &                  NFREAT,HEAP(FREEAT),BEGIN,QFIRST,
     &                  DELTA,SKIP,OUTPUT,
     &                  OQFORM,OFORM,OSCALE,OOFFST)
      END IF
      END DO
C
      CALL VCLOSE(OUTPUT,'KEEP',ERROR)
C
      CALL FREHP(FREEAT,INTEG4(NATOM))
      CALL FREHP(ZZ,IREAL8(NATOM))
      CALL FREHP(YY,IREAL8(NATOM))
      CALL FREHP(XX,IREAL8(NATOM))
C
      IF (ORIENT) THEN
      CALL FREHP(PMASS,IREAL8(NATOM))
      CALL FREHP(PAIR,INTEG4(2*NATOM))
      CALL FREHP(FLAGS,INTEG4(NATOM))
      END IF
C
      RETURN
      END
C=====================================================================
      SUBROUTINE MKPAIR(NATOM,FLAGS,PAIR,NPAIR,X,Y,Z)
C
C routine makes pair list based on FLAGS; checks for initialized
C coordinates.
C
      IMPLICIT NONE
C input/output
      INCLUDE 'funct.inc'
      INTEGER NATOM
      INTEGER FLAGS(*)
      INTEGER PAIR(2,*), NPAIR
      DOUBLE PRECISION  X(*), Y(*), Z(*)
C local
      INTEGER I, NMISS
C begin
      NMISS=0
      DO I=1,NATOM
      IF (FLAGS(I).EQ.1) THEN
      IF (.NOT.INITIA(I,X,Y,Z)) THEN
      FLAGS(I)=0
      NMISS=NMISS+1
      END IF
      END IF
      END DO
      NPAIR=0
      DO I=1,NATOM
      IF (FLAGS(I).EQ.1) THEN
      NPAIR=NPAIR+1
      PAIR(1,NPAIR)=I
      PAIR(2,NPAIR)=I
      END IF
      END DO
      IF (NMISS.GT.0) THEN
      CALL WRNDIE(0,'<MKPAIR>','some selected atoms unkown')
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE TRTDYN
C
C creates a trajectory file
C
C For syntax see main parsing loop "HELP"
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/ouput
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'trtraj.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
C local
      INTEGER SKIP, BEGIN
      DOUBLE PRECISION DELTA
      LOGICAL QSKIP
C begin
C
      IF (TRINIT) THEN
C
C the initialization flag is set ===> set the defaults
      TRINIT=.FALSE.
C defaults
      IF (TRIND.NE.0) THEN
      CALL FREHP(TRIND,INTEG4(LTRIN))
      TRIND=0
      LTRIN=0
      END IF
      IF (NATOM.GT.0) THEN
      LTRIN=NATOM
      TRIND=ALLHP(INTEG4(LTRIN))
      TRNUM=NATOM
      CALL FILL4(HEAP(TRIND),NATOM,1)
      END IF
      TRFIRS=.TRUE.
C
      OFILE='OUTPUT'
      TRUNIT=6
      TRQFOR=.FALSE.
      TROFFS=800.0D0
      TRSCAL=10000.0D0
      TRFORM='12Z6'
      END IF
C parsing
      QSKIP=.FALSE.
      CALL PUSEND('WRITe-DYNAmics>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('WRITe-DYNAmics>')
      CALL MISCOM('WRITe-DYNAmics>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-write-trajectory')
C
      ELSE IF (WD(1:4).EQ.'OUTP'.AND.TRFIRS) THEN
      CALL NEXTFI('OUTPut-file=',OFILE)
      QSKIP=.FALSE.
      ELSE IF (WD(1:4).EQ.'ASCI'.AND.TRFIRS) THEN
      CALL NEXTLO('ASCIi=',TRQFOR)
      QSKIP=.FALSE.
      ELSE IF (WD(1:4).EQ.'FORM'.AND.TRFIRS) THEN
      CALL NEXTST('FORMat=',TRFORM)
      QSKIP=.FALSE.
      ELSE IF (WD(1:4).EQ.'SCAL'.AND.TRFIRS) THEN
      CALL NEXTF('SCALe=',TRSCAL)
      QSKIP=.FALSE.
      ELSE IF (WD(1:4).EQ.'OFFS'.AND.TRFIRS) THEN
      CALL NEXTF('OFFSet=',TROFFS)
      QSKIP=.FALSE.
      ELSE IF (WD(1:4).EQ.'SELE'.AND.TRFIRS) THEN
      CALL SELCTA(HEAP(TRIND),TRNUM,X,Y,Z,.TRUE.)
      QSKIP=.FALSE.
      ELSE IF (WD(1:4).EQ.'NEXT') THEN
      IF (TRFIRS) THEN
      CALL WRNDIE(-5,'WRTDYN',
     &    'initialization call to WRITe DYNAmics missing')
      ELSE
C this is for the subsequent call to WRITE DYNAMICS
      BEGIN=1
      SKIP=1
      DELTA=1./TIMFAC
      CALL WRITTC('CORD',NATOM,X,Y,Z,
     &                  TRNUM,HEAP(TRIND),BEGIN,TRFIRS,
     &                  DELTA,SKIP,TRUNIT,
     &                  TRQFOR,TRFORM,TRSCAL,TROFFS)
      QSKIP=.TRUE.
      END IF
      ELSE IF (WD(1:4).EQ.'RESE') THEN
C set the initialization flag
      TRINIT=.FALSE.
      TRUNIT=6
      IF (TRIND.NE.0) THEN
      CALL FREHP(TRIND,INTEG4(LTRIN))
      TRIND=0
      LTRIN=0
      END IF
      IF (NATOM.GT.0) THEN
      LTRIN=NATOM
      TRIND=ALLHP(INTEG4(LTRIN))
      TRNUM=NATOM
      CALL FILL4(HEAP(TRIND),NATOM,1)
      END IF
      TRFIRS=.TRUE.
C
      QSKIP=.TRUE.
      OFILE='OUTPUT'
      TRUNIT=6
      TRQFOR=.FALSE.
      TROFFS=800.0D0
      TRSCAL=10000.0D0
      TRFORM='12Z6'
      ELSE
      CALL CHKEND('WRITe-DYNAmics>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (NATOM.GT.0.AND..NOT.QSKIP) THEN
      IF (TRQFOR) THEN
      CALL ASSFIL(OFILE,TRUNIT,'WRITE','FORMATTED',ERROR)
      ELSE
      CALL ASSFIL(OFILE,TRUNIT,'WRITE','UNFORMATTED',ERROR)
      END IF
C
C
C set up the trajectory index lists
      CALL MAKIND(HEAP(TRIND),NATOM,TRNUM)
C
      IF (.NOT.ERROR) THEN
      TRFIRS=.TRUE.
      BEGIN=1
      SKIP=1
      DELTA=1./TIMFAC
      CALL WRITTC('CORD',NATOM,X,Y,Z,
     &                  TRNUM,HEAP(TRIND),BEGIN,TRFIRS,
     &                  DELTA,SKIP,TRUNIT,
     &                  TRQFOR,TRFORM,TRSCAL,TROFFS)
      END IF
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE TRTINI
C
C initializes the trajectory file writing routine
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'trtraj.inc'
C local
C begin
C
C set the initialization flag
      TRINIT=.TRUE.
C
      TRUNIT=6
      TRIND=0
C
      TRFREEAT=0
      TRTEMP=0
      TRIX=0
      TRIY=0
      TRIZ=0
      RETURN
      END
C======================================================================
      SUBROUTINE TRTFRE
C
C resets trajectory read/write variables.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'trtraj.inc'
      INCLUDE 'funct.inc'
C local
      LOGICAL ERROR
C begin
C
C reset "read dynamics"
      IF (TRIND.NE.0) CALL FREHP(TRIND,INTEG4(LTRIN))
      IF (TRUNIT.NE.6) CALL VCLOSE(TRUNIT,'KEEP',ERROR)
      IF (TRIZ.NE.0) CALL FREHP(TRIZ,INTEG4(TRNAT))
      IF (TRIY.NE.0) CALL FREHP(TRIY,INTEG4(TRNAT))
      IF (TRIX.NE.0) CALL FREHP(TRIX,INTEG4(TRNAT))
      IF (TRTEMP.NE.0) CALL FREHP(TRTEMP,IREAL4(TRNAT))
      IF (TRFREEAT.NE.0) CALL FREHP(TRFREEAT,INTEG4(TRNAT))
C
C reset "write dynamics"
      TRINIT=.TRUE.
      TRUNIT=6
      TRIND=0
      TRFIRS=.TRUE.
C
      RETURN
      END
