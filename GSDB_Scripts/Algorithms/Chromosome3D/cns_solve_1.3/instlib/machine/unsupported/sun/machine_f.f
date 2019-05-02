C=====================================================================
C machine-dependent system-service routines (such as TIME and DATE)
C=====================================================================
      SUBROUTINE VDATE(A,AMAX,ALEN)
C
C this routine returns the date
C +++++++++++++++++++++++++++++++++++++++
C machine dependent sun/UNIX version - works for sun sparcstation 2
C +++++++++++++++++++++++++++++++++++++++
C
      IMPLICIT NONE
C input/output
      CHARACTER*(*) A
      INTEGER AMAX, ALEN
C begin
      IF (AMAX.LT.9) THEN
      WRITE(6,'(A)') ' %VDATE-ERR: string to small'
      ELSE
      A=' '
      CALL DATE(A)
      ALEN=9
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE VTIME(A,AMAX,ALEN)
C
C this routine returns day time
C +++++++++++++++++++++++++++++++++++++++
C machine dependent sun/UNIX version
C +++++++++++++++++++++++++++++++++++++++
C
      IMPLICIT NONE
C input/output
      CHARACTER*(*) A
      INTEGER AMAX, ALEN
C begin
      IF (AMAX.LT.8) THEN
      WRITE(6,'(A)') ' %VTIME-ERR: string to small'
      ELSE
      ALEN=AMAX
      CALL TIME(A)
      CALL TRIMM(A,ALEN)
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE VCPU(SECS)
C
C routine returns system time of current job in seconds
C
C +++++++++++++++++++++++++++++++++++++++
C machine dependent sun/UNIX version
C +++++++++++++++++++++++++++++++++++++++
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      DOUBLE PRECISION SECS
C local
      EXTERNAL STIME
      INTEGER STIME
      REAL SEXS(2) 
C begin
      CALL ETIME(SEXS)
      SECS=SEXS(1)
      RETURN
      END
C=====================================================================
      SUBROUTINE VINTIM
C
C initializes timer and saves initial date and time
C
C +++++++++++++++++++++++++++++++++++++++
C machine dependent sun/UNIX version
C +++++++++++++++++++++++++++++++++++++++
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C local
      INTEGER  L
      CHARACTER*8 TIME
      CHARACTER*11 DATE
C begin
      STRCPU=0.0d0 
      CALL VTIME(TIME,8,L)
      CALL VDATE(DATE,11,L)
      STRTIM=TIME
      STRDAT=' '
      STRDAT(1:L)=DATE(1:L)
      RETURN
      END
C=====================================================================
      SUBROUTINE VCHKTT(STREAM,QTERM)
C
C checks wheter unit STREAM is assigned to an interactive terminal
C
C +++++++++++++++++++++++++++++++++++++++
C machine dependent sun/UNIX version
C +++++++++++++++++++++++++++++++++++++++
C
C By Axel Brunger, 14-APR-84
C
      IMPLICIT NONE
C input/output
      INTEGER STREAM
      LOGICAL QTERM
C local
      LOGICAL ISATTY
C begin
      QTERM = ISATTY(STREAM)
      RETURN
      END
C=====================================================================
      SUBROUTINE SYTASK
C
C parses UNIX system commands
C Author: Axel T. Brunger
C 
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C local
      INTEGER SYSTEM, ISYS
      INTEGER RECMAX, RECLEN
      PARAMETER (RECMAX=160)
      CHARACTER*(RECMAX) RECORD
C begin
C
C run it through LINSUB to make symbol substutions
      CALL LINSUB(RECORD,RECMAX,RECLEN,' ',0,.FALSE.)
      ISYS=SYSTEM(RECORD(1:RECLEN))
      IF (ISYS.EQ.127) WRITE(6,'(A)') ' Can''t execute shell'
      CURSOR=COMLEN           
      RETURN
      END
C=====================================================================
      SUBROUTINE GETNAM(USERNM,MAX,LENGTH)
C
C Gets the username associated with the current process and places
C it in USERNM.
C +++++++++++++++++++++++++++++++++++++++
C machine dependent sun/UNIX version
C +++++++++++++++++++++++++++++++++++++++
C
C Author: Axel T. Brunger
C
C input/output
      CHARACTER*(*) USERNM
      INTEGER MAX, LENGTH
C local
      CALL GETLOG(USERNM)
      LENGTH=MAX
      CALL TRIMM(USERNM,LENGTH)
      RETURN
      END
C=====================================================================
      SUBROUTINE GETSYS(USERNM,MAX,LENGTH)
C
C Gets the SYSTEM associated with the current process and places
C it in USERNM.
C +++++++++++++++++++++++++++++++++++++++
C machine dependent sun/UNIX version
C +++++++++++++++++++++++++++++++++++++++
C
      IMPLICIT NONE
C input/output
      CHARACTER*14 USERNM
      INTEGER MAX, LENGTH
C begin
      CALL COPYST(USERNM,MAX,LENGTH,'SUN/SunOS',9)
      RETURN
      END
C=====================================================================
      SUBROUTINE VPROMP(UNIT,PROMPT)
C
C Makes a prompt at an interactive terminal without c/r.
C
C +++++++++++++++++++++++++++++++++++++++
C machine dependent sun/UNIX version
C +++++++++++++++++++++++++++++++++++++++
C
      IMPLICIT NONE
C input/output
      INTEGER UNIT
      CHARACTER*(*) PROMPT
C begin
      WRITE(6,'(A,$)') PROMPT
      RETURN
      END
C*********************************************************************
C=====================================================================
C=== machine-dependent I/O routines (for file opening and closing)
C=====================================================================
      SUBROUTINE NEXTVERSION(FILE)
C
C Creates "version numbers" on sun if a "NEW" file
C  is already there.
C  Input : filename  or filename_version
C  Output : filename_version+1
C
C +++++++++++++++++
C universal version
C +++++++++++++++++
C
C Michael Nilges, 21-MAY-87
C =========================
C
      IMPLICIT NONE
C input/output
      CHARACTER*(*) FILE
C local
      INTEGER I, FILLEN
      CHARACTER*10 NUMBERS
      CHARACTER*5 TEMP
      CHARACTER V
      INTEGER VERSION,NEWVERSION,VN
      INTEGER DECODI
      LOGICAL OK
C begin
      NUMBERS='1234567890'
      FILLEN=LEN(FILE)
      CALL TRIMM(FILE,FILLEN)
C Find position of last "_"
      I=FILLEN
      DO WHILE ((I.GT.2).AND.(INDEX(NUMBERS,FILE(I:I)).GT.0))
      V=FILE(I:I)
      I=I-1
      END DO
C Create version number
      IF ((FILE(I:I).EQ.'_').AND.(INDEX(NUMBERS,V).GT.0)) THEN
        TEMP=FILE((I+1):FILLEN)
        VERSION=DECODI(TEMP,FILLEN-I,OK)
        NEWVERSION=VERSION + 1
        TEMP=' '
      IF (MOD(NEWVERSION,10).EQ.0) THEN
          CALL ENCODI(NEWVERSION,TEMP,FILLEN-I+1,VN)
          FILE((I+1):FILLEN+1)=TEMP
      ELSE
          CALL ENCODI(NEWVERSION,TEMP,FILLEN-I,VN)
          FILE((I+1):(FILLEN))=TEMP
      END IF
      ELSE
        FILE((FILLEN+1):(FILLEN+2))='_1'
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE VOPEN(UNIT,FILE,FORM,ACCESS,ERR)
C
C Opens a file.
C
C +++++++++++++++++++++++++++++++++++++++
C machine dependent sun/UNIX version
C +++++++++++++++++++++++++++++++++++++++
C
C Axel Brunger, 30-APR-85
C =======================
C modified (create "version numbers" on sun): Michael Nilges 21-MAY-87
C
      IMPLICIT NONE
C input/output    
      INCLUDE 'iochan.inc'
      INTEGER UNIT
      CHARACTER*(*) FILE, FORM, ACCESS
      CHARACTER*132  LSTFIL
      LOGICAL ERR
C local
      INTEGER I, FILLEN
      LOGICAL EXISTS
C begin
      ERR=.FALSE.
      IF (FILE.EQ.'$$INITIAL$$') THEN
C---------------------------------------------------------------
C special section for initialization
C open standard input (channel 5) shared and readonly
CC      OPEN(FORM='FORMATTED',STATUS='OLD',ACCESS='SEQUENTIAL',
CC     &     UNIT=5)
C
C open standard output (channel 6)
CC      OPEN(FORM='FORMATTED',STATUS='NEW',ACCESS='SEQUENTIAL',
CC     &     UNIT=6)
C
C Initialize list of free FORTRAN units.
C It is assumed that system input is assigned to unit 5
C and system output is assigned to unit 6 and all other
C units are deassigned.
      DO I=1,MAXUN
      IFREEU(I)=0
      END DO
      IFREEU(5)=10
      IFREEU(6)=1
      ELSE
      IF (FORM.NE.'FORMATTED'.AND.FORM.NE.'UNFORMATTED') THEN
      WRITE(6,'(3A)')
     & ' %VOPEN-ERR: unknown format qualifier "',FORM,'"'
      ERR=.TRUE.
      ELSE
C---------------------------------------------------------------
C write append access
      IF (ACCESS.EQ.'APPEND') THEN
      OPEN(FILE=FILE,FORM=FORM,STATUS='OLD',ACCESS='APPEND',ERR=9999,
     &     UNIT=UNIT)
C---------------------------------------------------------------
C read access
      ELSE IF (ACCESS.EQ.'READ') THEN
      OPEN(FILE=FILE,FORM=FORM,STATUS='OLD',ACCESS='SEQUENTIAL',
     &     ERR=9999,UNIT=UNIT)
C---------------------------------------------------------------
C write access
      ELSE IF (ACCESS.EQ.'WRITE') THEN
      INQUIRE(FILE=FILE,EXIST=EXISTS)
      IF (EXISTS) THEN
      LSTFIL=FILE
      FILLEN=LEN(FILE)
      CALL TRIMM(FILE,FILLEN)
      CALL RENAME(LSTFIL(1:FILLEN),LSTFIL(1:FILLEN)//';1') 
      DO WHILE (EXISTS)
      CALL NEXTVERSION(LSTFIL)
      INQUIRE(FILE=LSTFIL,EXIST=EXISTS)
      IF (EXISTS) THEN
      CALL RENAME(LSTFIL,LSTFIL(1:FILLEN)//';2')
      CALL RENAME(LSTFIL(1:FILLEN)//';1',LSTFIL)
      CALL RENAME(LSTFIL(1:FILLEN)//';2',LSTFIL(1:FILLEN)//';1')
      END IF
      END DO
      CALL RENAME(LSTFIL(1:FILLEN)//';1',LSTFIL)
      END IF
      OPEN(FILE=FILE,FORM=FORM,STATUS='NEW',ACCESS='SEQUENTIAL',
     &     ERR=9999,UNIT=UNIT)
C---------------------------------------------------------------
      ELSE
      WRITE(6,'(3A)')
     & ' %VOPEN-ERR: unknown access qualilifer "',ACCESS,'"'
      ERR=.TRUE.
      END IF
      GOTO 8888
9999  ERR=.TRUE.
8888  CONTINUE
C
      IF (.NOT.ERR) THEN
C
C put appropriate code in IFREEU array:
C  +10 read formatted
C  +1  write/append formatted
C  -1  write/append unformatted
C  -10 read unformatted
      IF (FORM.EQ.'FORMATTED') THEN
      IFREEU(UNIT)=1
      ELSE
      IFREEU(UNIT)=-1
      END IF
      IF (ACCESS.EQ.'READ') IFREEU(UNIT)=IFREEU(UNIT)*10
      END IF
      END IF
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE VINQRE(MODE,NAME,MAXLEN,LENGTH,QOPEN,QFORM,QWRITE,
     &                  UNIT)
C
C file inquiry by file name or FORTRAN unit
C Flag QOPEN indicates whether file or unit is "open".
C Flag QFORM indicates whether file was opened formatted.
C Flag QWRITe indicates whether file was opened write-access.
C For inquiry by unit MAXLEN has to be specified (max length of NAME)
C and LENGTH returns with the length of NAME.
C For inquiry by file the two names INPUT and OUTPUT are reserved
C for the standard input and output channels 5 and 6.
C
C +++++++++++++++++++++++++++++++++++++++
C machine dependent sun/UNIX version
C +++++++++++++++++++++++++++++++++++++++
C
      IMPLICIT NONE
C I/O
      INCLUDE 'iochan.inc'
      CHARACTER*(*) MODE, NAME
      INTEGER MAXLEN, LENGTH
      LOGICAL QOPEN, QFORM, QWRITE
C sun
      logical qexist
      INTEGER UNIT
C begin
      QOPEN=.TRUE.
      IF (MODE.EQ.'FILE'.AND.NAME.EQ.'INPUT') THEN
      UNIT=5
      ELSE IF (MODE.EQ.'FILE'.AND.NAME.EQ.'OUTPUT') THEN
      UNIT=6
      ELSE IF (MODE.EQ.'FILE') THEN
C sun
      inquire(file=name,exist=qexist)
      if (qexist) then
      INQUIRE(FILE=NAME,OPENED=QOPEN,NUMBER=UNIT)
      else
        qopen=.false.
        unit=-1
      endif
      ELSE IF (MODE.EQ.'UNIT') THEN
      INQUIRE(UNIT=UNIT,OPENED=QOPEN,NAME=NAME)
      LENGTH=MAXLEN
      CALL TRIMM(NAME,LENGTH)
      END IF
C
C if file is open then get QFORM and QWRITE flags
      IF (QOPEN) THEN
      IF ((+1).EQ.(IFREEU(UNIT))) THEN
      QFORM=.TRUE.
      QWRITE=.TRUE.
      ELSE IF ((+10).EQ.(IFREEU(UNIT))) THEN
      QFORM=.TRUE.
      QWRITE=.FALSE.
      ELSE IF ((-1).EQ.(IFREEU(UNIT))) THEN
      QFORM=.FALSE.
      QWRITE=.TRUE.
      ELSE IF ((-10).EQ.(IFREEU(UNIT))) THEN
      QFORM=.FALSE.
      QWRITE=.FALSE.
      END IF
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE XXNXTFI(WD,WDLEN,WDMAX)
C
C does machine-dependent filename expansion/modification
C +++++++++++++++++++++++++++++++++++++++
C machine dependent sgi/UNIX version
C +++++++++++++++++++++++++++++++++++++++
C
C Axel T. Brunger
C ===============
C
C I/O
      IMPLICIT NONE
      CHARACTER*(*) WD
      INTEGER WDLEN, WDMAX
C local
      INTEGER L, LL, LLL, MXPATH
      PARAMETER (MXPATH=132)
      CHARACTER*(MXPATH) PATH, TEMP
      LOGICAL GOTIT
C begin
      GOTIT=.FALSE.
      LL=0
      DO WHILE (LL.LT.WDLEN.AND..NOT.GOTIT)
      LL=LL+1
      IF (WD(LL:LL).EQ.':'.AND.LL.GT.1) THEN
C
C   name DISK:file-name is translated to the value
C   of the ENVIRONMENT variable  DISK  , defined by
C
C         setenv DISK  path-name,
C
C   DISK:file-name ------>  path-name/file-name
C
      CALL GETENV(WD(1:LL-1),PATH)
      L=MXPATH
      CALL TRIMM(PATH,L)
      GOTIT=.TRUE.
      END IF
      END DO
C
      IF (GOTIT.AND.L.EQ.0) THEN
      TEMP='invalid or unknown path-name: "'//WD(1:LL-1)//'"'
      L=MXPATH
      CALL TRIMM(TEMP,L)
      CALL DSPERR('XXNXTFI',TEMP(1:L))
      ELSE IF (GOTIT) THEN
      CALL COPYST(TEMP,MXPATH,LLL,PATH,L)
      CALL ADDST(TEMP,MXPATH,LLL,'/',1)
      IF (WDLEN-LL.GT.0) THEN
      CALL ADDST(TEMP,MXPATH,LLL,WD(LL+1:WDLEN),WDLEN-LL)
      END IF
      CALL COPYST(WD,WDMAX,WDLEN,TEMP,LLL)
      END IF
C
      RETURN
      END
C=====================================================================  
      SUBROUTINE VCLOSE(UNIT,DISPOS,ERR)
C
C closes a file/unit with disposition DISPOS and sets the
C corresponding IFREEU element to zero.
C
C DONT CLOSE UNIT 5 OR 6 !!
C
C +++++++++++++++++++++++++++++++++++++++
C machine dependent sun/UNIX version
C +++++++++++++++++++++++++++++++++++++++
C
C By Axel Brunger, 21-APR-85
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'iochan.inc'
      INCLUDE 'comand.inc'
      INTEGER UNIT
      CHARACTER*(*) DISPOS
      LOGICAL ERR
C begin
      ERR=.FALSE.
      IF (UNIT.NE.6.AND.UNIT.NE.5.AND.UNIT.NE.DUNIT) THEN
         CLOSE(UNIT=UNIT,ERR=9999)
         IFREEU(UNIT)=0
      END IF
      GOTO 8888
9999  CONTINUE
      WRITE(6,'(A)') ' %VCLOSE-ERR: error during close '
      ERR=.TRUE.
8888  CONTINUE
      RETURN
      END
C=====================================================================
      SUBROUTINE GETLIN(UNIT,COMLYN,COMMAX,COMLEN,EOF,ERROR)
C
C gets a line from UNIT
C
C +++++++++++++++++++++++++++++++++++++++
C machine dependent sun/UNIX version
C +++++++++++++++++++++++++++++++++++++++
C
C Axel Brunger, 13-JAN-84
C =======================
C
      IMPLICIT NONE
C input/output
      INTEGER UNIT
      CHARACTER*(*) COMLYN
      INTEGER COMMAX, COMLEN
      LOGICAL EOF, ERROR
C local
C begin
      READ(UNIT,'(A)',END=6666,ERR=7777) COMLYN
      COMLEN=COMMAX
      CALL TRIMM(COMLYN,COMLEN)
C
C a carriage return is treated like a blank
      COMLEN=COMLEN+1
      COMLYN(COMLEN:COMLEN)=' '
C
      GOTO 8888
6666  EOF=.TRUE.
      COMLEN=1
      COMLYN=' '
      GOTO 8888
7777  ERROR=.TRUE.
      WRITE(6,'(A)') ' %GETLIN-ERR: ERROR reading stream file'
      WRITE(6,*)     ' UNIT=',UNIT
      WRITE(6,*)     ' COMLYN=',COMLYN(1:80)
8888  CONTINUE
      RETURN
      END
C
C======================================================================
C
      INTEGER FUNCTION IREAL8(N)
C
C This routine computes lengths needed to double precision
C temporary space.
C machine-dependent sun routine
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER N
C local
      IREAL8=2*N
      RETURN
      END
C======================================================================
      INTEGER FUNCTION ICPLX8(N)
C
C This routine computes lengths needed to double complex
C temporary space.
C machine-dependent sun routine
C
C Author: Axel T. Brunger
      IMPLICIT NONE
C input/output
      INTEGER N
C local
      ICPLX8=4*N
      RETURN
      END
C======================================================================
      INTEGER FUNCTION ICPLX4(N)
C
C This routine computes lengths needed to single complex
C temporary space.
C machine-dependent sun routine
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER N
C local
      ICPLX4=2*N
      RETURN
      END
C======================================================================
      INTEGER FUNCTION IREAL4(N)
C
C This routine computes lengths needed to single precision
C temporary space. 
C machine-dependent sun routine
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C intput/output
      INTEGER N
C begin
      IREAL4=N
      RETURN
      END
C======================================================================
      INTEGER FUNCTION INTEG4(N)
C
C This routine computes lengths needed to INTEGER
C temporary space.
C machine-dependent sun routine
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER N
C begin
      INTEG4=N
      RETURN
      END
C======================================================================
      INTEGER FUNCTION ILOGIC(N)
C
C routine returns length for LOGICAL type arrays
C machine-dependent sun routine
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER N
C begin
      ILOGIC=N
      RETURN
      END
C=====================================================================
      INTEGER FUNCTION ICHAR4(N)
C
C routine returns length needed to CHARACTER*4 arrays
C should be used only in connection with CALLST,
C CINIST, CFREST. The reason for putting in these
C special CHARACTER*4 routines is that there is no equivalence
C possible between a real or integer word.
C and a CHARACTER*4 word.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER N
C local
      ICHAR4=N
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE SETBSL(BACKSL)
C
C Set variable 'BACKSL' to backslash character
C
      IMPLICIT NONE
C
      CHARACTER*1 BACKSL
C
      BACKSL = '\\'
      END
C
C ===========================================================================
C
      SUBROUTINE VHOSTNAME(HOSTNM, HNLEN)
      IMPLICIT NONE
C I/O
      CHARACTER*(*)  HOSTNM
      INTEGER        HNLEN
C
C Get host name
C
C begin
      CALL FHOSTNAME(HOSTNM, LEN(HOSTNM), HNLEN)
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE CNSHELP(STRING)
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'comand.inc'
C
      CHARACTER*(*) STRING
C
      INTEGER MXPATH
      PARAMETER (MXPATH=132)
      CHARACTER*(MXPATH) PATH
      INTEGER IUNIT, STLEN
      CHARACTER*(COMMAX) LINE
      LOGICAL ERR, QEXIST
C begin
      CALL GETENV('CNS_HELPLIB',PATH)
      STLEN=MXPATH
      CALL TRIMM(PATH,STLEN)
C
      IF (STLEN.EQ.0) THEN
      WRITE(6,'(A)')' %CNS-HELP-ERR: CNS_HELPLIB undefined'
      ELSE
      IF ( WRNLEV.GE.10 ) THEN
      WRITE(6,'(A,A)')' CNS-HELP: Help library is ',
     &                PATH(1:STLEN)
      END IF
C
      IF ( PATH(STLEN:STLEN).NE.'/' ) THEN
      IFILE=PATH(1:STLEN) // '/' // STRING(1:LEN(STRING))
      STLEN=STLEN + LEN(STRING) + 1
      ELSE
      IFILE=PATH(1:STLEN) // STRING(1:LEN(STRING))
      STLEN=STLEN + LEN(STRING)
      END IF
C
      INQUIRE(FILE=IFILE(1:STLEN),EXIST=QEXIST)
      IF (QEXIST) THEN
      CALL ASSFIL(IFILE(1:STLEN),IUNIT,'READ','FORMATTED',ERR)
      IF ( .NOT.ERR ) THEN
      WRITE (6,'(1X)')
 5    CONTINUE
      READ(IUNIT,'(A)',END=11,ERR=10) LINE
      STLEN=COMMAX
      CALL TRIMM(LINE,STLEN)
      WRITE (6,'(1X,A)') LINE(1:STLEN)
      GOTO 5
 10   WRITE (6,'(A)')
     &    ' %CNS-HELP-ERR: Error reading from help library file'
      WRITE (6,'(A,A)')
     &    ' %CNS-HELP-ERR: check ',IFILE(1:STLEN)
 11   CONTINUE
      CALL VCLOSE(IUNIT,'KEEP',ERROR)
      WRITE (6,'(1X)')
      ELSE
      WRITE (6,'(A)')
     &    ' %CNS-HELP-ERR: Error opening help library file'
      WRITE (6,'(A,A)')
     &    ' %CNS-HELP-ERR: check ',IFILE(1:STLEN)
      END IF
      ELSE
      WRITE (6,'(A)')
     &    ' %CNS-HELP-ERR: Help library file does not exist'
      WRITE (6,'(A,A)')
     &    ' %CNS-HELP-ERR: check ',IFILE(1:STLEN)
      END IF
      END IF
C
      RETURN
      END
C
C======================================================================
C
