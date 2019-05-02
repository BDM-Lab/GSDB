      INTEGER FUNCTION GETATN(SID,RID,REN,IUP,MARK)
C
C This function returns the atom number for the atom named by its
C segid SID, resid RID, residue name REN and iupac name IUP.
C If not found GETATN is set to MARK
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      CHARACTER*4 SID, RID, REN, IUP
      INTEGER MARK
C local
      INTEGER I
      LOGICAL FOUND
C begin
      FOUND=.FALSE.
      I=0
      DO WHILE ((.NOT.FOUND).AND.(I.LT.NATOM))
      I=I+1
      FOUND=(SEGID(I).EQ.SID)
      IF (FOUND) THEN
      FOUND=(RESID(I).EQ.RID)
      IF (FOUND) THEN
      FOUND=(RES(I).EQ.REN)
      IF (FOUND) THEN
      FOUND=(TYPE(I).EQ.IUP)
      END IF
      END IF
      END IF
      END DO
      IF (.NOT.FOUND) THEN
      GETATN=MARK
      ELSE
      GETATN=I
      END IF
C
      RETURN
      END
C================================================================
      LOGICAL FUNCTION HYDROG(I)
C
C This function returns a true if the atom is deemed to be a
C hydrogen atom. all tests for hydrogen atoms should use this
C routine.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/ouput
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INTEGER I
C begin
      IF (I.LE.0 .OR. I.GT.NATOM) THEN
      HYDROG=.FALSE.
      ELSE
      HYDROG=AMASS(I).LT.3.5D0
      END IF
C
      RETURN
      END
C
      LOGICAL FUNCTION INITIA(ATOM,X,Y,Z)
C
C Function true if coordinates for ATOM defined.
C
C Author: Axel Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER ATOM
      DOUBLE PRECISION X(*), Y(*), Z(*)
C local
      DOUBLE PRECISION ANUM
      PARAMETER (ANUM=9998.9D0)
C begin
      IF (ATOM.GT.0) THEN
      INITIA=(X(ATOM).LT.ANUM .AND. Y(ATOM).LT.ANUM .AND.
     2        Z(ATOM).LT.ANUM)
      ELSE
      INITIA=.FALSE.
      END IF
      RETURN
      END
C================================================================
      SUBROUTINE ATOMID(ATOM,SID,RID,REN,AC)
C
C Given the atomnumber ATOM this routine returns the SEGID, RESID,
C RESIDUE and IUPAC atom name. If ATOM out of range question marks
C are returned.
C
C Author: Axel Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INTEGER ATOM
      CHARACTER*4 SID, RID, REN, AC
C
      IF (ATOM.GT.NATOM.OR.ATOM.LE.0) THEN
      SID='x??x'
      RID='x??x'
      REN='x??x'
      AC='x??x'
      ELSE
      SID=SEGID(ATOM)
      RID=RESID(ATOM)
      REN=RES(ATOM)
      AC=TYPE(ATOM)
      END IF
C
      RETURN
      END
C================================================================
      SUBROUTINE GGUBFS(RANDOM)
C
C uniform (0,1) random number generator
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      DOUBLE PRECISION RANDOM
      INCLUDE 'seed.inc'
C begin
      IF (SEED.LT.1.0D0)  SEED=314564.D0
      SEED=MOD(16807.0D0*SEED,D2P31M)
      RANDOM=SEED/D2P31
      RETURN
      END
C
      SUBROUTINE RDTITB(UNIT,FORM)
C
C reads title from an formatted file or an binary file.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'ctitla.inc'
      INCLUDE 'comand.inc'
      INTEGER UNIT
      CHARACTER*(*) FORM
C local
      INTEGER L, J
C begin
      IF (FORM.EQ.'FORMATTED') THEN
      READ(UNIT,'(/I8)',ERR=9999,END=8888) NTITLE
      IF (NTITLE .LE. 0) THEN
      READ( UNIT, '(A)',ERR=9999, END=8888)
      ELSE
      DO J = 1, NTITLE
      TITLE(J) = ' '
      READ (UNIT, '(A)', ERR = 9999, END = 8888) TITLE(J)
      ENDDO
      ENDIF
      ELSE IF (FORM.EQ.'BINARY') THEN
      DO J=1,MXTITL
      TITLE(J)=' '
      END DO
      READ(UNIT,ERR=9999,END=8888) NTITLE,(TITLE(J)(1:80),J=1,NTITLE)
      END IF
C
C print out current title
      DO J=1,NTITLE
      L=TITMAX
      CALL TRIMM(TITLE(J),L)
      WRITE(6,'(A)') TITLE(J)(1:L)
      END DO
      GOTO 7777
9999  ERROR=.TRUE.
      GOTO 7777
8888  EOF=.TRUE.
7777  RETURN
      END
C================================================================
      SUBROUTINE WRTITL(UNIT,IOP)
C
C This routine writes titles. It accepts an type option
C   3 - Star format
C   2 - fixed field
C   1 - Brookhaven-standard
C   0 - free field
C  -1 - binary
C The date and username will be added at the end, the
C filename will be added at the beginning of the title.
C Already existing dates or name will be superceeded.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'ctitla.inc'
      INCLUDE 'version.inc'
      INTEGER UNIT, IOP
C local
      LOGICAL QFORM, QOPEN, QWRITE
      INTEGER I, J, L, LD, LU, LT
      CHARACTER*(WORD_SIZE) FILNAM
      CHARACTER*(WORD_SIZE) USERNM
      CHARACTER*11  DATE
      CHARACTER*8  TIME
      CHARACTER*(1) CNSPTMP
      CHARACTER*1  LEAD
      LOGICAL QDATE, QVERS, COND
C externals
      INTEGER   LNBLNK
      EXTERNAL  LNBLNK
C begin
C
C put in filename
CCC modification ATB 4/27/08
      COND=.FALSE.
      IF (NTITLE.EQ.0) THEN
      COND=.TRUE.
      ELSE
      IF (INDEX(TITLE(1),'FILENAME').EQ.0) COND=.TRUE. 
      END IF
      IF (COND) THEN
CCC
      NTITLE=NTITLE+1
      IF (NTITLE.GT.MXTITL) THEN
      CALL WRNDIE(-1,'WRTITL',
     &  'MXTITL (COMAND) exceeded. --> shorten title or recompile')
      ELSE
      DO J=NTITLE,2,-1
      TITLE(J)=TITLE(J-1)
      END DO
      END IF
      END IF
      CALL VINQRE('UNIT',FILNAM,WORD_SIZE,L,QOPEN,QFORM,QWRITE,UNIT)
      TITLE(1)=' REMARKS FILENAME="'//FILNAM(1:L)//'"'
C
C check for date and version
      QDATE=.TRUE.
      QVERS=.TRUE.
      DO I=NTITLE-1,NTITLE
      IF (INDEX(TITLE(I),'DATE').GT.0) THEN
      CALL VDATE(DATE,11,LD)
      CALL VTIME(TIME,8,LT)
      CALL GETNAM(USERNM,WORD_SIZE,LU)
      TITLE(I)=' REMARKS DATE:'//DATE(1:LD)//'  '//
     &TIME(1:LT)//'       created by user: '//USERNM(1:LU)
      QDATE=.FALSE.
      ELSEIF (INDEX(TITLE(I),'VERSION').GT.0) THEN
      CNSPTMP=CNSPATCH
      IF (CNSPTMP.EQ.'0') THEN
      TITLE(I)=' REMARKS VERSION:'//CNSVERSION
      ELSE
      TITLE(I)=' REMARKS VERSION:'//CNSVERSION//CNSPATCH
      END IF
      QVERS=.FALSE.
      END IF
      END DO
C
C if date and version not added before do it now
      IF (QDATE) THEN
      NTITLE=NTITLE+1
      IF (NTITLE.GT.MXTITL) THEN
      CALL WRNDIE(-1,'WRTITL',
     &'MXTITL (COMAND) exceeded. --> shorten title or recompile')
      ELSE
      CALL VDATE(DATE,11,LD)
      CALL VTIME(TIME,8,LT)
      CALL GETNAM(USERNM,WORD_SIZE,LU)
      TITLE(NTITLE)=' REMARKS DATE:'//DATE(1:LD)//'  '//
     &TIME(1:LT)//'       created by user: '//USERNM(1:LU)
      END IF
      END IF
      IF (QVERS) THEN
      NTITLE=NTITLE+1
      IF (NTITLE.GT.MXTITL) THEN
      CALL WRNDIE(-1,'WRTITL',
     &'MXTITL (COMAND) exceeded. --> shorten title or recompile')
      ELSE
      CNSPTMP=CNSPATCH
      IF (CNSPTMP.EQ.'0') THEN
      TITLE(NTITLE)=' REMARKS VERSION:'//CNSVERSION
      ELSE
      TITLE(NTITLE)=' REMARKS VERSION:'//CNSVERSION//CNSPATCH
      END IF
      END IF
      END IF
C
      IF ((-1).EQ.(IOP)) THEN
C
C for binary output
      WRITE(UNIT,ERR=9) NTITLE,(TITLE(J)(1:80),J=1,NTITLE)
C
      ELSE IF ((0).EQ.(IOP)) THEN
C
C for card output
      DO J=1,NTITLE
      L=TITMAX
      CALL TRIMM(TITLE(J),L)
      WRITE(UNIT,'(A)',ERR=9) TITLE(J)(1:L)
      END DO
C
      ELSE IF ((1).EQ.(IOP)) THEN
C
C for Brookhaven standard format.
C note that in PDB format REMARK has to begin in the first column!
      DO J=1,NTITLE
      L=TITMAX
      CALL TRIMM(TITLE(J),L)
      IF (L.GE.9) WRITE(UNIT,'(A)',ERR=9) 'REMARK'//TITLE(J)(9:L)
      END DO
C
      ELSE IF ((2).EQ.(IOP)) THEN
C
C for fixed field output
      WRITE(UNIT,'(/I8,A)',ERR=9) NTITLE,' !NTITLE'
      WRITE(UNIT,'(A)',ERR=9) (TITLE(J),J=1,NTITLE)
C
      ELSE IF ((3).EQ.(IOP)) THEN
C
C for STAR output
        J = 10
        DO I = 1, NTITLE
          IF (TITLE(I)(1:9) .NE. ' REMARKS ') J = 1
        END DO
        LEAD = ';'
        DO I = 1, NTITLE
          L = LNBLNK(TITLE(I))
          IF (L .GE. J) THEN
            WRITE(UNIT, '(A, 1X, A)', ERR=9) LEAD, TITLE(I)(J:L)
          ELSE IF (LEAD .NE. ' ') THEN
            WRITE(UNIT, '(A)', ERR=9) LEAD
          ELSE
            WRITE(UNIT, '()', ERR=9)
          END IF
          LEAD = ' '
        END DO
        WRITE(UNIT,'(A/)',ERR=9) ';'
C
      END IF
C-error-labels----
      GOTO 99
9     ERROR=.TRUE.
99    CONTINUE
C-----------------
      IF (ERROR) WRITE(6,'(A)') ' %WRTITL-ERR: ERROR DURING WRITE'
      RETURN
      END
C================================================================
      SUBROUTINE ASSFIL(NAME,UNIT,ACCESS,FORM,ERR)
C
C Checks whether file NAME already open and then returns the unit
C otherwise it gets a free FORTRAN unit, tries to open the file
C NAME and gets the full filename specification.
C
C ACCESS : 'WRITE', 'READ', 'APPEND'
C FORM : 'FORMATTED', 'UNFORMATTED'
C NAME : <filename>, 'INPUT', 'OUTPUT'
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) NAME, ACCESS, FORM
      INTEGER UNIT
      LOGICAL ERR
C local
      INTEGER L, TMAX, TLEN, LENDUM
      PARAMETER (TMAX=WORD_SIZE)
      CHARACTER*(TMAX) T
      LOGICAL QOPEN, QFORM, QWRITE
C begin
C
      ERR=.FALSE.
      L=LEN(NAME)
      CALL TRIMM(NAME,L)
C
C inquire about the file
      LENDUM=0
      CALL VINQRE('FILE',NAME,0,LENDUM,QOPEN,QFORM,QWRITE,UNIT)
C
C if file is already open do some error checking
      IF (QOPEN) THEN
      ERR=.TRUE.
      IF (FORM.EQ.'FORMATTED'.AND..NOT.QFORM) THEN
      WRITE(6,'(2A)') ' %ASSFIL-ERR: file ',
     & 'already open but not formatted as requested'
      ELSE IF (FORM.EQ.'UNFORMATTED'.AND.QFORM) THEN
      WRITE(6,'(2A)') ' %ASSFIL-ERR: file ',
     & 'already open but not unformatted as requested'
      ELSE IF ((ACCESS.EQ.'WRITE'.OR.ACCESS.EQ.'APPEND')
     & .AND..NOT.QWRITE) THEN
      WRITE(6,'(2A)') ' %ASSFIL-ERR: file ',
     & 'already open but not for write access as requested'
      ELSE IF (ACCESS.EQ.'READ'.AND.QWRITE) THEN
      WRITE(6,'(2A)') ' %ASSFIL-ERR: file ',
     & 'already open but not for read access as requested'
      ELSE
      ERR=.FALSE.
      END IF
      ELSE
C
C file not open, get a free FORTRAN unit
      UNIT=0
      QOPEN=.TRUE.
      DO WHILE ((UNIT.LT.99).AND.(QOPEN))
      UNIT=UNIT+1
      CALL VINQRE('UNIT',T,TMAX,TLEN,QOPEN,QFORM,QWRITE,UNIT)
      END DO
      IF (QOPEN) THEN
      WRITE(6,'(A)')' %ASSFIL-ERR: no free unit available'
      ERR=.TRUE.
      ELSE
      CALL VOPEN(UNIT,NAME,FORM,ACCESS,ERR)
      IF (ERR) THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      WRITE(6,'(2A)') ' %ASSFIL-ERR: error opening file ',NAME(1:L)
      END IF
      ELSE
C
C get full file name to print message about opening a file
      CALL VINQRE('UNIT',T,TMAX,TLEN,QOPEN,QFORM,QWRITE,UNIT)
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(3A)') ' ASSFIL: file ',T(1:TLEN),' opened.'
      END IF
      END IF
      END IF
      END IF
      IF (ERR.AND.(ENDACT(ENDIND).EQ.'GO  '))
     & CALL WRNDIE(-1,'ASSFIL','Error accessing file')
      RETURN
      END
C================================================================
      SUBROUTINE WRNDIE(LEV,SCALL,MESSAG)
C
C Routine produces an error message and terminates
C program execution unless the program is in interactive
C mode.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER LEV
      CHARACTER*(*) MESSAG,SCALL
      INCLUDE 'cns.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'comand.inc'
C begin
      WRITE(6,'(4A)') ' %',SCALL,' error encountered: ',MESSAG
C
      IF (BOMLEV.LE.-5) THEN
      WRITE(6,'(A)')
     & '   (CNS is in mode: SET ABORT=OFF END)'
      ELSE IF (BOMLEV.EQ.5) THEN
      WRITE(6,'(A)')
     & '   (CNS is in mode: SET ABORT=ALL END)'
      ELSE
      WRITE(6,'(A)')
     & '   (CNS is in mode: SET ABORT=NORMal END)'
      END IF
C
      IF (LEV.LE.BOMLEV) THEN
C
      IF (.NOT.QTERM) THEN
      WRITE(6,'(A)')
     & ' *****************************************************',
     & ' ABORT mode will terminate program execution. ',
     & ' *****************************************************'
      CALL DIE
      ELSE
      WRITE(6,'(A)')
     & ' WARNING: program encountered a fatal error.',
     & '    However, in interactive mode, program execution',
     & '    will continue.  Proceed at your own risk.'
      END IF
      END IF
C
      RETURN
      END
C=====================================================================
      SUBROUTINE DIE
C
      IMPLICIT NONE
C begin
      WRITE(6,'(A)') ' Program will stop immediately.'
      CALL PRFINAL
      STOP
      END
C=====================================================================
      SUBROUTINE PRFINAL
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
C local
      INTEGER           L
      DOUBLE PRECISION  SECS
      CHARACTER*11      DATE
      CHARACTER*8       TIME
      CHARACTER*11      MAXBYTES, MAXOVERH
C begin
      CALL QRYALLOC( 3, MAXBYTES)
      CALL QRYALLOC(-3, MAXOVERH)
      WRITE(6,'(10X,A)')
     &'============================================================'
      WRITE(6,'(10X,3A)')
     & ' Maximum dynamic memory allocation: ',MAXBYTES,' bytes'
      WRITE(6,'(10X,3A)')
     & ' Maximum dynamic memory overhead:   ',MAXOVERH,' bytes'
      WRITE(6,'(10X,4A)')
     & ' Program started at: ',STRTIM,' on ',STRDAT
      CALL VTIME(TIME,8,L)
      CALL VDATE(DATE,11,L)
      WRITE(6,'(10X,4A)')
     & ' Program stopped at: ',TIME,' on ',DATE(1:L)
      CALL VCPU(SECS)
      WRITE(6,'(10X,A,F12.4,A)')
     & ' CPU time used: ',SECS,' seconds'
      WRITE(6,'(10X,A)')
     &'============================================================'
      RETURN
      END
C================================================================
      SUBROUTINE DSPCPU(TEXT)
C
C displays current CPU time including "TEXT".
C
C Author: Axel T. Brunger
C I/O
      IMPLICIT NONE
      INCLUDE 'timer.inc'
      CHARACTER*(*) TEXT
C local
      DOUBLE PRECISION SECS
C begin
      CALL VCPU(SECS)
      WRITE(6,'(A,A,F12.6,A)') TEXT,' CPU=',SECS,' sec.'
      RETURN
      END
C================================================================
      SUBROUTINE DSPERR(PROMPT,LINE)
C
C markes an error by "^^^^^" and writes an message stored in LINE
C
C Author: Axel T. Brunger
C I/O
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) PROMPT, LINE
C local
      CHARACTER*(COMMAX) ERRLYN
      INTEGER I, L, MAXERR
      PARAMETER (MAXERR=100)
C
C MAXERR is the maximum number of allowed parsing errors
C begin
      L=LEN(PROMPT)
      IF (PROMPT(L:L).EQ.'>') L=L-1
      WRITE(6,'(5A)') ' %',PROMPT(1:L),'-ERR: ',LINE,':'
      ERRLYN=' '
      DO I=CUROLD,CURSOR-1
      ERRLYN(I:I)='^'
      END DO
      ERRLYN(CURSOR-1:CURSOR-1)='^'
      WRITE(6,'(2A)') ' ',COMLYN(1:CURSOR)
      WRITE(6,'(2A)') ' ',ERRLYN(1:CURSOR-1)
      NERRPA=NERRPA+1
      IF (NERRPA.EQ.MAXERR) THEN
      CALL WRNDIE(-1,'PARSER',
     &  'Encountered too many parsing errors.')
      END IF
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE SHOWAL
C
C prints pertinent information about the status of CNS
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C local
      INTEGER I
      LOGICAL QOPEN, QFORM, QWRITE
      CHARACTER*12 FORM
      CHARACTER*6  ACCESS
C begin
      CALL SHOW
C
C inquire all open files
      WRITE(6,'(A)') ' Files currently open:'
      DO I=1,99
      CALL VINQRE('UNIT',WDT,WDTMAX,WDTLEN,QOPEN,QFORM,QWRITE,I)
      IF (QOPEN) THEN
      IF (QFORM) THEN
      FORM=' '
      ELSE
      FORM=' unformatted'
      END IF
      IF (QWRITE) THEN
      ACCESS=' write'
      ELSE
      ACCESS=' read'
      END IF
      WRITE(6,'(A,3A)') ' -> ',WDT(1:WDTLEN),ACCESS,FORM
      END IF
      END DO
      RETURN
      END
C
C======================================================================
C
