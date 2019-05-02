      SUBROUTINE DEFPROCSET(ERR)
C
C procedure definition
C (moved into parser.F to avoid an unexplained compiler warning on HPs )
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      LOGICAL ERR
C local
      CHARACTER*(VARMAX) PARNAM
      INTEGER PNMLEN
      INTEGER RECNO,FIRSTREC,PARLIN,NPAREN,LRECNO
      LOGICAL FOUND,QPART
      INTEGER MATCHLEN,NEWLEN,SCOPE
      LOGICAL QQSUBS,QQDEF,QQCAPIT,QENDPROC
      LOGICAL CLOOP
C begin
      LRECNO = 0
      FIRSTREC = 0
C save qsubs and qdef
      QQSUBS=QSUBS
      QSUBS=.FALSE.
      QQDEF=QDEF
      QDEF=.FALSE.
      QQCAPIT=QCAPIT
      QCAPIT=.FALSE.
C remember name
      PARNAM = WD(1:WDLEN)
      PNMLEN = WDLEN
      CALL DEFNEWREC(RECNO)
      IF(RECNO.EQ.0) THEN
      ERR = .TRUE.
      ELSE
      DEFNAMTXT(RECNO)=PARNAM(1:PNMLEN)
      DEFNAMLEN(RECNO)=PNMLEN
      FIRSTREC = RECNO
C now parse the parameter block
      CALL DEFNEXTW2(5)
C check for problems
      IF((WD(1:1).NE.'(').OR.(WDLEN.NE.1)) THEN
      WRITE (6,'(A)') ' %PROCEDURE-ERR: ( expected but not found.'
      ERR=.TRUE.
      ELSE
      DEFPARTXT(RECNO) = '('
      DEFPARLEN(RECNO) = 1
      NPAREN=1
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      QNEWLN = .FALSE.
      CALL DEFNEXTW2(5)
C
C check for stream switch in a MODULE file (BAD!)
      IF((WD(1:1).EQ.'(').AND.(WDLEN.EQ.1).AND.
     & (.NOT.QQUOT)) THEN
      NPAREN=NPAREN+1
      ELSE IF((WD(1:1).EQ.')').AND.(WDLEN.EQ.1).AND.
     & (.NOT.QQUOT)) THEN
      NPAREN=NPAREN-1
      END IF
      IF(QQUOT) THEN
C quoted strings will need to have quotes restored
      WDLEN=WDLEN+2
      END IF
C
      IF(WDLEN+1.GT.COMMAX) THEN
      CALL DSPERR('DEFMACSET','parameter too long')
      WDLEN=COMMAX-1
      END IF
C
      NEWLEN=WDLEN+DEFPARLEN(RECNO)+1
      IF((NEWLEN.GT.COMMAX).OR.
     & (QNEWLN.AND.(DEFPARLEN(RECNO).GT.0))) THEN
      LRECNO=RECNO
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      DEFMULTREC(LRECNO)=RECNO
      ELSE
      ERR=.TRUE.
      END IF
      END IF
C now add word if we have room
      IF(RECNO.GT.0) THEN
C add word
      IF(QQUOT) THEN
C restore quotes
      DEFPARTXT(RECNO)((DEFPARLEN(RECNO)+2):
     & (DEFPARLEN(RECNO)+WDLEN-1))=WD(1:WDLEN-2)
      DEFPARTXT(RECNO)((DEFPARLEN(RECNO)+1):
     & (DEFPARLEN(RECNO)+1))='"'
      DEFPARTXT(RECNO)((DEFPARLEN(RECNO)+WDLEN):
     & (DEFPARLEN(RECNO)+WDLEN))='"'
      ELSE
C not quoted
      IF(NPAREN.EQ.0) THEN
C remove space before last ')' -- purely cosmetic
      IF(DEFPARLEN(RECNO).GT.0) THEN
      IF(DEFPARTXT(RECNO)(DEFPARLEN(RECNO):DEFPARLEN(RECNO))
     &   .EQ.' ') DEFPARLEN(RECNO) = DEFPARLEN(RECNO) - 1
      END IF
      END IF
      DEFPARTXT(RECNO)((DEFPARLEN(RECNO)+1):
     & (DEFPARLEN(RECNO)+WDLEN))=WD(1:WDLEN)
      END IF
C
      DEFPARLEN(RECNO)=DEFPARLEN(RECNO)+WDLEN
C add trailing space (expected by parser after last word)
      DEFPARLEN(RECNO)=DEFPARLEN(RECNO)+1
      DEFPARTXT(RECNO)(DEFPARLEN(RECNO):DEFPARLEN(RECNO))=' '
      END IF
C            if(recno.gt.0)
C
      IF (NPAREN.LE.0) CLOOP=.FALSE.
      END DO
C
C now parse the body of the procedure
      LRECNO=RECNO
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      DEFMULTREC(LRECNO)=RECNO
      ELSE
      ERR=.TRUE.
      END IF
C
      IF(RECNO.GT.0) THEN
C delimit the start of the procedure body by setting
C the unused name field to 'BODY' for this record (only)
      DEFNAMTXT(RECNO) = 'BODY'
      DEFNAMLEN(RECNO) = 4
C add remaining text on the current line
      IF(CURSOR.LT.COMLEN) THEN
      DEFPARLEN(RECNO) = COMLEN - CURSOR
      DEFPARTXT(RECNO)(1:DEFPARLEN(RECNO)) = COMLYN(CURSOR+1:COMLEN)
      END IF
C
      CURSOR = COMLEN
C
      QCAPIT = .TRUE.
      QENDPROC = .FALSE.
      DO WHILE(.NOT.QENDPROC)
      CALL NEXTW2('PROC-BODY>')
C
C copy literal commands until a line is found that starts with
C ENDProcedure
      IF(EOF) THEN
      CALL DSPERR('PROC-BODY>','EOF encountered')
      QENDPROC = .TRUE.
      ELSE IF(WD(1:4).NE.'ENDP') THEN
C
      IF(DEFPARLEN(RECNO).GT.0) THEN
      LRECNO = RECNO
      CALL DEFNEWREC(RECNO)
      DEFMULTREC(LRECNO) = RECNO
      END IF
C
      DEFPARTXT(RECNO)(1:COMLEN) = COMLYN(1:COMLEN)
      DEFPARLEN(RECNO) = COMLEN
      CURSOR = COMLEN
C
      ELSE
      QENDPROC = .TRUE.
      END IF
      END DO
C
      IF(RECNO.GT.0) THEN
C remove last line if it is not the first in the body
      IF(DEFPARLEN(RECNO).EQ.0) THEN
      IF(DEFNAMLEN(RECNO).EQ.0) THEN
      DEFMULTREC(LRECNO) = 0
      CALL DEFFREEREC(RECNO,.FALSE.)
      ELSE
      DEFPARLEN(RECNO) = 1
      DEFPARTXT(RECNO) = ' '
      END IF
      END IF
      END IF
C
      END IF
C
      IF((FIRSTREC.GT.0).AND..NOT.ERR) THEN
C print new definition of parameter
      IF(WRNLEV.GE.10) THEN
      CALL DEFGETSSC(SCOPE,DEFCURSCOPE)
      CALL DEFCDUMP(FIRSTREC,.FALSE.,SCOPE,3,' ',0)
      END IF
C kill existing parameter with this name (if it exists)
      MATCHLEN=PNMLEN
      CALL DEFGETSSC(SCOPE,DEFCURSCOPE)
      CALL DEFFINDREC(PARNAM,MATCHLEN,SCOPE,3,.TRUE.,
     & FOUND,PARLIN,QPART,.TRUE.)
C add new parameter on to top of parameter list
      CALL DEFGETSSC(SCOPE,DEFCURSCOPE)
      CALL DEFINSREC(FIRSTREC,SCOPE,3)
      END IF
      END IF
      END IF
C
      IF(ERR.AND.(FIRSTREC.GT.0)) THEN
C error condition
      CALL DEFFREEREC(FIRSTREC,.TRUE.)
      END IF
C
      QSUBS=QQSUBS
      QDEF=QQDEF
      QCAPIT=QQCAPIT
      RETURN
      END
C=================================================================
      SUBROUTINE DOECHO(PROMPT,CC,CCLEN,QPRTOO,QPROECHO)
C
C from the input to be echoed, remove parameter scope directives,
C and then echo the input
C
C set QPRTOO = .true. for echoing taking place as part of the
C              get-next-line cycle in NEXTW2
C            = .false. for echoing taking place as part of BACEND
C
C set QPROECHO = .true. when echoing of commands is necessary
C                even when in terminal mode (QTERM=.true.)
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) PROMPT
      CHARACTER*(COMMAX) CC
      INTEGER CCLEN
      LOGICAL QPRTOO,QPROECHO
C local
      CHARACTER*(COMMAX) ECHOLYN
      INTEGER ECHOLEN, I
C begin
      ECHOLYN(1:CCLEN) = CC(1:CCLEN)
      ECHOLEN = CCLEN
      I=1
      DO WHILE(I.LE.(ECHOLEN-3))
C remove scope directives
      IF((ECHOLYN(I:I+3).EQ.' &%N').OR.
     &   (ECHOLYN(I:I+3).EQ.' &%K'))THEN
      ECHOLYN(I:ECHOLEN-4)=ECHOLYN(I+4:ECHOLEN)
      ECHOLEN=ECHOLEN-4
      END IF
      I=I+1
      END DO
C    echo line and journal file (for interactive session) handling
      IF (QECHO.AND.((.NOT.QTERM).OR.QPROECHO)) THEN
      IF(.NOT.QPROECHO) THEN
      WRITE(6,'(3A)') ' ',PROMPT,ECHOLYN(1:ECHOLEN)
      ELSE
      WRITE(6,'(A)') ECHOLYN(1:ECHOLEN)
      END IF
      ELSE IF (QPRTOO.AND.QTERM.AND.
     &         (BUFIND.EQ.BUFFIL).AND.PRUNIT.GE.0) THEN
      WRITE(PRUNIT,'(A)') ECHOLYN(1:ECHOLEN)
      END IF
C
      RETURN
      END
C=================================================================
      SUBROUTINE MODULE
C
C MODULE invocation routine
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C local
      INTEGER MCUNIT,INUNIT
      CHARACTER*(COMMAX) HLDMAC
      INTEGER HLDMCL
      LOGICAL QMCADD,QINADD,QQTERM
      INTEGER MCSTRM
      LOGICAL ABORT
      LOGICAL QCLOSE,QSWTCH
      LOGICAL QQSUBS, QQDEF
C begin
C
C initialize
      ABORT=.FALSE.
      QCLOSE=.TRUE.
      QSWTCH=.TRUE.
C check for invocation from interactive mode
      QQTERM=QTERM
      QTERM = .FALSE.
C
C remember current stream information
      MCSTRM=NSTRM
      MCUNIT=ISTRM(NSTRM)
      INUNIT=ISTRM(NSTRM-1)
C turn off stream switching, symbol substitution, buffering
      QMCADD=QADD(NSTRM)
      QINADD=QADD(NSTRM-1)
      QQSUBS=QSUBS
      QSUBS=.FALSE.
      QQDEF=QDEF
      QDEF=.FALSE.
C
C if this is an inline module, unshield it for parameter definitions
      IF(DEFINLCNT(DEFCURSCOPE).NE.0) DEFINLCNT(DEFCURSCOPE) = 1
C make sure the next word is a parenthesis
      CALL DEFNEXTW2(2)
C check for problems
      IF((NSTRM.LT.2).OR.(NSTRM.NE.MCSTRM).OR.
     & (ISTRM(NSTRM).NE.MCUNIT)) THEN
      WRITE (6,'(A)') ' %MODULE-ERR: ( expected but not found.'
      ERROR=.TRUE.
      QCLOSE=.FALSE.
      QSWTCH=.FALSE.
      ELSE IF((WD(1:1).NE.'(').OR.(WDLEN.NE.1)) THEN
      WRITE (6,'(A)') ' %MODULE-ERR: ( expected but not found.'
      ERROR=.TRUE.
      ELSE
C specify scope/stream into which defines will be stored
      DEFSETSCOPE = DEFCURSCOPE
      DEFEVLSCOPE = DEFCURSCOPE
C parse parameter block
      CALL DEFMACSET(ABORT,ERROR,2)
C
C check for problems
      IF((NSTRM.LT.2).OR.(MCSTRM.NE.NSTRM).OR.
     & (ISTRM(NSTRM).NE.MCUNIT)) THEN
      WRITE (6,'(A)') ' %MODULE-ERR: bad module definiton.'
      ERROR=.TRUE.
      QCLOSE=.FALSE.
      QSWTCH=.FALSE.
      END IF
C
      IF(.NOT.ERROR) THEN
C store remaining portion of current line from MODULE file
      CALL COPYST(HLDMAC,COMMAX,HLDMCL,
     & COMLYN(CURSOR:COMLEN),COMLEN-CURSOR+1)
C
C switch to previous stream file
      NSTRM=NSTRM-1
CCC      CALL VCHKTT(ISTRM(NSTRM),QTERM)
      QTERM = QQTERM
C
C    read the old "hold" line
      IF(HLDLEN(NSTRM).GT.0) THEN
      CALL COPYST(COMLYN,COMMAX,COMLEN,HLDLYN(NSTRM),HLDLEN(NSTRM))
      ELSE
      COMLYN=' '
      COMLEN=1
      END IF
C
      IF(QMCADD) THEN
C
C    copy this special line into rotating command buffer
      BUFIND=BUFIND+1
      BUFFIL=BUFFIL+1
      BUFSTK=MOD(BUFIND-1,BUFMAX)+1
      CALL COPYST(BUFLYN(BUFSTK),COMMAX,BUFLEN(BUFSTK),COMLYN,COMLEN)
      END IF
C
      CURSOR=1
      CUROLD=1
C now check for parenthesis
      CALL DEFNEXTW2(1)
C check for problems
      IF((NSTRM.NE.(MCSTRM-1)).OR.(ISTRM(NSTRM).NE.INUNIT)) THEN
      WRITE (6,'(A)') ' %MODULE-ERR: ( expected but not found.'
      ERROR=.TRUE.
      QSWTCH=.FALSE.
      ELSE IF((WD(1:1).NE.'(').OR.(WDLEN.NE.1)) THEN
      WRITE (6,'(A)') ' %MODULE-ERR: ( expected but not found.'
      ERROR=.TRUE.
      QSWTCH=.FALSE.
      ELSE
C force variables to be expanded into the surrounding scope
      DEFEVLSCOPE = MAX(1,DEFCURSCOPE-1)
C parse invocation parameter block
      CALL DEFMACSET(ABORT,ERROR,1)
C check for problems
      IF(ERROR.OR.ABORT) QSWTCH=.FALSE.
      IF((NSTRM.NE.(MCSTRM-1)).OR.(ISTRM(NSTRM).NE.INUNIT)) THEN
      WRITE (6,'(A)') ' %MODULE-ERR: bad module invocation.'
      ERROR=.TRUE.
      QSWTCH=.FALSE.
      END IF
      IF(.NOT.ERROR) THEN
      IF(.NOT.ABORT) THEN
C
C everything is peachy - switch completely over to MODULE file
      QADD(NSTRM)=QINADD
      QCLOSE=.FALSE.
      QSWTCH=.FALSE.
      NSTRM=NSTRM+1
C
C save the current hold line
      IF((COMLEN.GE.CURSOR).AND.(COMLEN.GT.0)) THEN
      CALL COPYST(HLDLYN(NSTRM-1),COMMAX,HLDLEN(NSTRM-1),
     & COMLYN(CURSOR:COMLEN),COMLEN-CURSOR+1)
      ELSE
      HLDLYN(NSTRM-1)=' '
      HLDLEN(NSTRM-1)=1
      END IF
C
C switch the stream and initialize
      EOF=.FALSE.
      ERROR=.FALSE.
      QSAVEW=.FALSE.
      QSUBS=.TRUE.
      QFILNM=.FALSE.
      QEXPRS=.FALSE.
      COMLEN=0
      CURSOR=0
C a macro file is NEVER an interactive terminal
      QTERM=.FALSE.
C
      QADD(NSTRM)=QMCADD
C
C    read the old "hold" line
      CALL COPYST(COMLYN,COMMAX,COMLEN,HLDMAC,HLDMCL)
      IF(QADD(NSTRM)) THEN
C
C    copy this special line into rotating command buffer ( see also
C    PROCESS-ADD-COMMAND in MISCOM )
      BUFIND=BUFIND+1
      BUFFIL=BUFFIL+1
      BUFSTK=MOD(BUFIND-1,BUFMAX)+1
      CALL COPYST(BUFLYN(BUFSTK),COMMAX,BUFLEN(BUFSTK),HLDMAC,HLDMCL)
      END IF
      CURSOR=1
      CUROLD=1
C reactivate INLINE definition shielding (if appropriate)
      IF(DEFINLCNT(DEFCURSCOPE).EQ.1) DEFINLCNT(DEFCURSCOPE) = 2
      ELSE
C MODULE invocation aborted
      WRITE (6,'(A)') ' MODULE: module expansion aborted.'
C restore interactive state
      QTERM = QQTERM
      END IF
      END IF
      END IF
      END IF
      END IF
C
      IF(ERROR) THEN
C MODULE invocation errored
      CALL WRNDIE(-5,'MODULE','module expansion failed.')
      END IF
C restore parser state
      QSUBS=QQSUBS
      QDEF=QQDEF
      IF(ISTRM(NSTRM).EQ.INUNIT) QADD(MCSTRM-1)=QINADD
C
C if MODULE file needs to be closed (error or abort)...
      IF(QCLOSE) THEN
      CALL VCLOSE(MCUNIT,'KEEP',ERROR)
C destroy buffers and scope
      CALL DEFCLOSE(MCSTRM)
      CALL DEFKILSCOPE
      END IF
C
C if original stream needs to be restored (error or abort)...
      IF(QSWTCH) THEN
      NSTRM=MCSTRM-1
CCC      CALL VCHKTT(ISTRM(NSTRM),QTERM)
      QTERM = QQTERM
      CALL COPYST(COMLYN,COMMAX,COMLEN,HLDLYN(NSTRM),HLDLEN(NSTRM))
      IF (QADD(NSTRM)) THEN
C    copy this special line into rotating command buffer ( see also
C    PROCESS-ADD-COMMAND in MISCOM )
      BUFIND=BUFIND+1
      BUFFIL=BUFFIL+1
      BUFSTK=MOD(BUFIND-1,BUFMAX)+1
      CALL COPYST(BUFLYN(BUFSTK),COMMAX,BUFLEN(BUFSTK),HLDMAC,HLDMCL)
      END IF
      CURSOR=1
      CUROLD=1
      END IF
      DONE=.FALSE.
      ERROR=.FALSE.
C restore stream switching and MODULE expasion
      RETURN
      END
C
      SUBROUTINE CHKEND(PROMPT,LDONE)
C
C 1. check END, END IF and END LOOP statements
C 2. issue a warning message if LDONE=.FALSE. and WD not "END"
C 3. skip all messages when LDONE=.TRUE. (fast scanning mode
C    in MISCOM)
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/ouput
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) PROMPT
      LOGICAL LDONE
C begin
      IF (EOF) THEN
      LDONE=.TRUE.
      CALL DSPERR(PROMPT,'EOF encountered')
      ELSE IF (WD(1:4).EQ.'END ') THEN
C
C process the END IF statement
      IF (ENDKEY(ENDIND)(1:4).EQ.'IF  ') THEN
      QNEWLN=.FALSE.
      CALL NEXTWD(PROMPT)
      IF(QNEWLN) THEN
C END IF not on one line
      CALL SAVEWD
      ELSE IF (WD(1:2).EQ.'IF') THEN
      ENDIND=ENDIND-1
      ELSE
      IF (.NOT.LDONE) THEN
      CALL DSPERR(PROMPT,'END IF expected')
      ELSE
      CALL SAVEWD
      END IF
      END IF
C
C process the END LOOP statement
      ELSE IF (ENDKEY(ENDIND)(1:4).EQ.'LOOP') THEN
      CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'LOOP') THEN
      CALL NEXTWD('LOOP-label=')
      IF (WD(1:4).NE.ENDLBL(ENDIND)) THEN
      CALL DSPERR(PROMPT,'wrong LOOP-label')
      ENDIND=ENDIND-1
      ELSE IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      CALL BACEND(PROMPT)
      ELSE IF (ENDACT(ENDIND).EQ.'SKIP') THEN
      ENDIND=ENDIND-1
      END IF
      ELSE
      IF (.NOT.LDONE) THEN
      CALL DSPERR(PROMPT,'END LOOP expected')
      ELSE
      CALL SAVEWD
      END IF
      END IF
C
C now process END statements not involved in control statements
      ELSE
      LDONE=.TRUE.
      ENDIND=ENDIND-1
      END IF
C
C normally we want to have a warning message if an un-recognized
C word is encountered
      ELSE IF (.NOT.LDONE) THEN
      CALL DSPERR(PROMPT,'unrecognized command')
C
C however, for LDONE=TRUE, i.e. the call from MISCOM, nothing is done
      ELSE
      END IF
      IF (ENDIND.LT.0) THEN
      CALL DSPERR(PROMPT,'END stack underflow')
      END IF
      RETURN
      END
C
      SUBROUTINE PUSEND(KEY)
C
C inverse action to CHKEND, i.e. invoke a new "END" level
C
C Author: Axel T. Brunger
C =======================
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) KEY
C
      IF (ENDIND.GE.ENDMAX) THEN
      CALL WRNDIE(-5,'PUSEND',
     & 'exceeded ENDMAX (COMAND) parameter --> recompile program')
      ELSE
      ENDIND=ENDIND+1
      END IF
      ENDKEY(ENDIND)=KEY
C position cursor at beginning of current word "KEY"
      ENDCUR(ENDIND)=CURSOR-WDLEN
      ENDRET(ENDIND)=BUFIND
      CALL PUSACT('GO  ')
      RETURN
      END
C
      SUBROUTINE NEXTW2(PROMPT)
C
C Gets the next word from the current stream file.
C
C
C Words are defined as:
C  1. single character words
C      By default, single character words are defined as
C      "(", ")", ":", "@", "=".
C
C      In expression mode (QEXPRS) the following characters are
C      treated as single words as well: ">", "<",
C      "*", "/", "^", "+", "-" (unless +,- are part of exponents of
C      scientific notation)
C
C  2. any sequence of quoted strings or non-blank characters, which are not
C     single character words or blanks, enclosed by blanks or
C     single character words.
C     A quoted string is defined as any sequence of characters enclosed by
C     double quotes (").  The double quote itself is represented by (""").
C
C
C
C Characters within braces ("{", "}") are ignored.  Characters
C following a "!" are ignored until a carriage return is reached.  The
C "!" is converted into a blank character.
C
C By default, symbol substitutions are carried out.  Symbols are indicated by
C the "$" as the first character of the word.  No substitutions
C are performed in quoted strings.
C
C Unless characters are within quotes all characters are converted
C to upper-case.   All non-printable characters (such as tabs) are
C converted into blanks.
C
C On input:
C =========
C The flag QSUBS determines whether to carry out symbol substitutions.
C The flag QEXPS determines whether the extended set of single
C character words is active.
C The flag QFILNM determines whether we're parsing a filename.
C The flag QCAPIT determines whether capitalization is done.
C The string PROMPT is the prompt string that will appear on the
C output device.
C
C On output:
C ==========
C The parsed word is returned in string WD, with length WDLEN, the
C flag QQUOT indicates whether the parsed word is a quoted string
C or a substituted quoted string.
C
C All input/output except PROMPT is performed via common block comand.inc
C
C Author: Axel Brunger
C ====================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      CHARACTER*(*) PROMPT
C local
      CHARACTER*1 CTEMP
      INTEGER I, SCAN, ICHR, TEST
      LOGICAL COND, QWORD, QSINGL, QDELIM, QEXPAND, QREPEAT
      LOGICAL FOUND, NOSKIP, CLOOP1, CLOOP2
      CHARACTER*2 WDTYP
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
      CHARACTER*1 BACKSL
C begin
      CALL SETBSL(BACKSL)
      IF (QSAVEW) THEN
C
C this flag is set by routine SAVEWD.  It indicates that NEXTWD
C should simply return the previously parsed word.  Thus, SAVEWD
C has the effect of pushing the current word back into the input
C stream.
      QSAVEW=.FALSE.
      ELSE
C
      IF (QFILNM) THEN
C
C if we're in filename parsing mode, the only single character word
C allowed is the equal sign
      TEST=0
C
      ELSE IF (QEXPRS) THEN
C
C if we're in expression parsing mode we want to use the extended
C set of single character words
      TEST=2
C
      ELSE
C
C if we're in normal parsing mode we don't want to use the extended
C set of single character words
      TEST=1
      END IF
C
C position cursor
      CURSOR=MAX(0,CURSOR-1)
C
C get a word
      WDLEN=0
      WD=' '
      SCAN=1
      QQUOT=.FALSE.
      QWORD=.FALSE.
      QSINGL=.FALSE.
      QDELIM=.FALSE.
      CLOOP1=.TRUE.
      DO WHILE (CLOOP1)
C
C get next character
      CURSOR=CURSOR+1
      QREPEAT=.TRUE.
      DO WHILE (QREPEAT)
C
C in case we need to run through GETNEXTLINE again later, we can
C set QREPEAT to true
      QREPEAT=.FALSE.
C
      IF(CURSOR.GT.COMLEN) THEN
C if we've reached the end of the physical record we've to get the next line
C======================================================================
C GET NEXT LINE   begin
C
C set newline flag (for "END IF" vs "...END /LF IF..." )
      QNEWLN=.TRUE.
C
C    the following loop makes sure that we switch to the previous stream
C    file when EOF reached. However if EOF is encountered on the first
C    stream file (NSTRM=1) the routine exits, and leaves the EOF flag on
C    and returns the word END.
C check to see if there are any buffered command lines that were created
C upon expansion of a define parameter.
      CALL DEFGETLIN(COMLYN,COMMAX,COMLEN,QEXPAND,.FALSE.)
C
      IF(.NOT.QEXPAND) THEN
C nothing was stored up, so proceed as usual
C
      EOF=.FALSE.
      ERROR=.FALSE.
      CLOOP2=.TRUE.
      DO WHILE (CLOOP2)
C
C    make a prompt at the interactive terminal
      IF (QTERM.AND.(BUFIND.EQ.BUFFIL.OR..NOT.QADD(NSTRM))) THEN
      CALL VPROMP(ISTRM(NSTRM),PROMPT)
      END IF
C
C get next record either from the define buffer, the
C rotating command buffer or from the stream file.
C
C
C single-@ stream mode?
      IF (QADD(NSTRM)) THEN
C
C increment rotating command buffer pointer
      BUFIND=BUFIND+1
      BUFSTK=MOD(BUFIND-1,BUFMAX)+1
C
C if nothing left in the rotating command buffer:
      IF (BUFIND.GT.BUFFIL) THEN
C
C increment top of rotating command buffer indicator
      BUFFIL=BUFFIL+1
C
C get new line from buffered procedures if available
      CALL DEFGETLIN(BUFLYN(BUFSTK),COMMAX,
     & BUFLEN(BUFSTK),QEXPAND,.TRUE.)
      IF(.NOT.QEXPAND) THEN
C if not, then get from the stream
      CALL GETLIN(ISTRM(NSTRM),BUFLYN(BUFSTK),COMMAX,BUFLEN(BUFSTK),
     & EOF,ERROR)
C
CCC modification ATB 5/11/08 to reduce parsing error messages after ERROR
      EOF=EOF.OR.ERROR
CCC
      ELSE
C          force echoing in terminal mode
      IF(QTERM) CALL DOECHO(PROMPT,BUFLYN(BUFSTK),BUFLEN(BUFSTK),
     &                      .TRUE.,.TRUE.)
      END IF
C
      END IF
C            (bufind.gt.buffil)
C
C something in rotating command buffer -> copy from buffer
C copy buffer record into command record for parsing
      CALL COPYST(COMLYN,COMMAX,COMLEN,BUFLYN(BUFSTK),BUFLEN(BUFSTK))
      ELSE
C
C double-@ stream mode
C any buffered procedure lines available?
      CALL DEFGETLIN(COMLYN,COMMAX,COMLEN,QEXPAND,.TRUE.)
      IF(.NOT.QEXPAND) THEN
C if not, then get a new line from the stream
      CALL GETLIN(ISTRM(NSTRM),COMLYN,COMMAX,COMLEN,EOF,ERROR)
C
CCC modification ATB 5/11/08 to reduce parsing error messages after ERROR
      EOF=EOF.OR.ERROR
CCC
      ELSE
      IF(QTERM) CALL DOECHO(PROMPT,COMLYN,COMLEN,.TRUE.,.TRUE.)
      END IF
C
C Note: The rotating command buffer remains untouched in @@ mode.
C
      END IF
C            (qadd(nstrm))
C
C
      COND=(EOF.AND.NSTRM.GT.1)
      IF (COND) THEN
      EOF=.FALSE.
C
C insert close-scope directive at end of the current BUFLYN
      IF(QADD(NSTRM)) THEN
      BUFLYN(BUFSTK)(BUFLEN(BUFSTK)+1:(BUFLEN(BUFSTK)+4))=' &%K'
      BUFLEN(BUFSTK)=BUFLEN(BUFSTK)+4
      END IF
C destroy scope associated with this stream
      CALL DEFKILSCOPE
C
C
C try to close old stream file
      CALL VCLOSE(ISTRM(NSTRM),'KEEP',ERROR)
C
C switch to previous stream file
      NSTRM=NSTRM-1
C
C    read the old "hold" line
      CALL COPYST(COMLYN,COMMAX,COMLEN,HLDLYN(NSTRM),HLDLEN(NSTRM))
      IF (QADD(NSTRM+1)) THEN
C
C    copy this special line into rotating command buffer ( see also
C    PROCESS-ADD-COMMAND in MISCOM )
      BUFIND=BUFIND+1
      BUFFIL=BUFFIL+1
      BUFSTK=MOD(BUFIND-1,BUFMAX)+1
      CALL COPYST(BUFLYN(BUFSTK),COMMAX,BUFLEN(BUFSTK),COMLYN,COMLEN)
      END IF
C
      COND=(COMLEN.LE.1)
      CALL VCHKTT(ISTRM(NSTRM),QTERM)
      END IF
      IF (.NOT.COND.OR.ERROR) CLOOP2=.FALSE.
      END DO
      IF (EOF.OR.ERROR) GOTO 9999
C
C
      CALL DOECHO(PROMPT,COMLYN,COMLEN,.TRUE.,.FALSE.)
      END IF
C
      CURSOR=1
      CUROLD=1
C
      END IF
C            if(cursor.gt.comlen)
C
CC    WRITE (6,*) DEFCURSCOPE,':',COMLYN(1:COMLEN),'>'
C
C GET NEXT LINE   end
C===================================================================
C Define and MODULE Parameter Substitution, and Scoping Routine
C
C Examines the current cursor position and looks to see
C if the next word is a $symbol or a &parameter.
C If any sort of substitution or scoping is required, it is
C carried out on the current command line.
C
C NOTE: this routine does not substitute $symbols, it only
C adds scoping information when QDEF and QSCOPE are
C True
C
C Parameters may contain multiple words or lines, and nesting is
C allowed.  Extra lines are buffered and then loaded into COMLYN
C as needed (see DEFGETLIN above).
C
C
      IF(.NOT.QQUOT) THEN
      CALL DEFSUBPAR(COMLYN,COMMAX,COMLEN,CURSOR,.FALSE.)
      END IF
C
C in case COMLEN was shortened to less than CURSOR this loop
C will run through GET NEXT LINE again;
      IF(CURSOR.GT.COMLEN) QREPEAT = .TRUE.
C
C===================================================================
C
      END DO
C            (QREPEAT)
C
C
C store current character in CHR
      CHR=COMLYN(CURSOR:CURSOR)
C
C do upper-case conversion and comment handling if we're not inside a
C quoted string
      IF (.NOT.QQUOT) THEN
C
      IF (QFILNM.OR.(.NOT.QCAPIT)) THEN
C
C filenames are case-sensitive (UNIX systems!), therefore
C only convert nonacceptable control characters into blanks
      IF (ASCIIM(ICHAR(CHR)).EQ.ICHAR(' ')) CHR=' '
      ICHR=ICHAR(CHR)
      ELSE
C
C convert everything to upper case and convert nonacceptable control
C characters into blanks.
      ICHR=ASCIIM(ICHAR(CHR))
      CHR=CHAR(ICHR)
C
C modification: 2/24/10:
C ignore a backslash escape if it precedes a blank or invalid character.
      IF (CHR.EQ.BACKSL .AND. CURSOR.LT.COMLEN .AND.
     &  ASCIIM(ICHAR(COMLYN(CURSOR+1:CURSOR+1))).EQ.ICHAR(' ')) CHR=' '
C
      END IF
C
C now: comment handling
      NOSKIP=SCAN.EQ.1
      IF      (CHR.EQ.'{' ) THEN
      SCAN=SCAN+1
      NOSKIP=.FALSE.
      ELSE IF (CHR.EQ.'}' ) THEN
      SCAN=SCAN-1
      ELSE IF (CHR.EQ.'!'.AND.NOSKIP ) THEN
      CURSOR=COMLEN
      CHR=' '
      ELSE IF (.NOT.QWORD.AND.NOSKIP.AND.CHR.NE.' ') THEN
C
C the condition above means that we've reached the beginning of a word
      QWORD=.TRUE.
C
C curold marks the beginning of the word.  This is used by DSPERR for
C nice error messages that point to the word.
      CUROLD=CURSOR
      END IF
      END IF
C            (.not.qquot)
C
C the next condition checks that we're inside a word and but not
C within comment braces.
      IF (NOSKIP.AND.QWORD) THEN
C
C first we check for quotes
      IF (CHR.EQ.'"') THEN
C
      IF (COMLYN(MAX(1,CURSOR-1):CURSOR+1).EQ.'"""') THEN
C
C special handling for """ which produces the word "
      WDLEN=WDLEN+1
      WD(WDLEN:WDLEN)='"'
      ELSE IF (QQUOT) THEN
C
C we've reached the end of the quoted string.  This also
C delimits the word.
      QDELIM=.TRUE.
      CURSOR=CURSOR+1
      ELSE
C
C this is the beginning of a quoted string
      QQUOT=.TRUE.
      END IF
      ELSE IF (QQUOT) THEN
C
C we're inside a quoted string
      WDLEN=WDLEN+1
      WD(WDLEN:WDLEN)=CHR
C
      ELSE
C
C we're not inside a quoted string
C
C Allow backslash escape to force inclusion of next character, if in
C normal parsing mode (TEST=1) where wildcard patterns may occur,
C and not the last character of the input line.
      IF (TEST.EQ.1 .AND. CHR.EQ.BACKSL .AND. CURSOR.LT.COMLEN) THEN
      WDLEN=WDLEN+1
      WD(WDLEN:WDLEN)=CHR
      WDLEN=WDLEN+1
      CURSOR=CURSOR+1
      WD(WDLEN:WDLEN)=COMLYN(CURSOR:CURSOR)
C
C now we check for delimiters (only ' ')
      ELSE IF (CHR.EQ.' ' ) THEN
      QDELIM=.TRUE.
C
C now we check for single character words.
C Note, that single character words also act as delimiters for the
C words surrounding them.
      ELSE IF (ASCIIS(ICHR).LE.TEST) THEN
      QSINGL=.TRUE.
C
C the "+" and "-" single character words require special attention
C since they may be part of the scientific notation for constants, e.g.,
C 1.0E-10.0. We check the previous characters in the current word.
      IF (QEXPRS) THEN
      IF ((CHR.EQ.'+'.OR.CHR.EQ.'-').AND.WDLEN.GE.2) THEN
C
C check whether it is part of exponent
      CTEMP=COMLYN(CURSOR-1:CURSOR-1)
      IF (CTEMP.EQ.'E'.OR.CTEMP.EQ.'e') THEN
C
C test whether the previous characters qualify as a number
      CALL CHKNUM(WD,WDLEN-1,FOUND)
      IF (FOUND) THEN
      QSINGL=.FALSE.
C +,- is part of normal word -> append it to WD.
      WDLEN=WDLEN+1
      WD(WDLEN:WDLEN)=CHR
      END IF
      END IF
      END IF
      END IF
C
      IF (QSINGL.AND.WDLEN.EQ.0) THEN
C
C parse single character words
      WDLEN=1
      WD(1:1)=CHR
      CURSOR=CURSOR+1
      END IF
C
      ELSE
C
C character is part of normal word -> append it to WD.
      WDLEN=WDLEN+1
      WD(WDLEN:WDLEN)=CHR
      END IF
C
      END IF
      END IF
      IF (QWORD.AND.(QSINGL.OR.QDELIM)) CLOOP1=.FALSE.
C prevent word overflows caused by symbol substitutions, Joe Krahn, 2/24/10 
C (requires that LEN(WD) <= STRING_SIZE+1)
      WDLEN=MIN(LEN(WD),WDLEN)
      END DO
      END IF
C
C now make symbol substitutions (words with a "$" as the first
C character)
      IF (QSUBS.AND..NOT.QQUOT.AND.
     & WD(1:1).EQ.'$'.AND.ENDACT(ENDIND).EQ.'GO  ') THEN
C
C
C if we're in FILENAME mode we've to convert the name of the
C symbol to uppercase!
      IF (QFILNM) THEN
      DO I=1,WDLEN
      CHR=WD(I:I)
      WD(I:I)=CHAR(ASCIIM(ICHAR(CHR)))
      END DO
      END IF
C
      CALL WDSUB(WD,WDMAX,WDLEN,FOUND,WDTYP,DPVAL,DCVAL)
      IF (.NOT.FOUND) THEN
      CALL DSPERR('WDSUB','symbol not found')
      END IF
C
C set the QQUOT flag if the substituted word is a quoted string
      IF (WDTYP.EQ.'ST') QQUOT=.TRUE.
      END IF
C
C
C error handling below.  If we've reached the end-of-file of all
C streams we want to make sure that the program terminates.
      GOTO 8888
9999  CALL COPYST(WD,WDMAX,WDLEN,'END',3)
      CALL DSPERR('NEXTF','EOF or ERROR encountered on input')
8888  CONTINUE
C
      RETURN
      END
C=================================================================
      SUBROUTINE NEXTWD(PROMPT)
C
C Gets the next word from the current stream file using
C routine NEXTW2 (see above).  If the word is "@" or "@@" then
C we've to try to switch to a new stream and then resume parsing.
C
C Author: Axel T. Brunger
C Modified for MODULE handling by Warren DeLano
C =======================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) PROMPT
C local
      INTEGER UNIT, I, J
      LOGICAL QSTRAM
      LOGICAL QMCHEK
      LOGICAL QQADD
      LOGICAL QQECHO, QINLNOK, QQTERM
      INTEGER OLDNSTRM
      LOGICAL CLOOP
C begin
C
C repeat until we've reached a normal word
      QMCHEK=.FALSE.
C
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
C
      CALL NEXTW2(PROMPT)
      J = CURSOR - 1
C
      QSTRAM=.FALSE.
C
C remember if we are in interactive mode...
      QQTERM = QTERM
C
C now check for the '@' statement (streaming of files)
C also allows nested @@s
      IF (WD(1:1).EQ.'@'.AND..NOT.QQUOT) THEN
C
C set the flag that we've to re-start parsing after opening the file
C and that we need to check if this is a MODULE file
      QSTRAM=.TRUE.
C
C store the old cursor position
      I=CURSOR-1
C
C @@file is read but contents is not copied to rotating command buffer.
      IF (COMLYN(CURSOR:CURSOR).EQ.'@') THEN
      CURSOR=CURSOR+1
C
C get the filename
      IFILE=' '
      CALL NEXTFI('@@-FILE=',IFILE)
C
C if we're active assign the file and switch the stream to it
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      CALL ASSFIL(IFILE,UNIT,'READ','FORMATTED',ERROR)
      IF (.NOT.ERROR) THEN
      REWIND(UNIT=UNIT)
      CALL INSTRM(UNIT)
      QADD(NSTRM)=.FALSE.
      QMCHEK=.TRUE.
C create new scope
      CALL DEFNEWSCOPE
      END IF
      ELSE
C
C always reset the QINLINE flag even we don't switch to a new
C stream because we're not in "GO" mode.  (e.g., within FALSE conditional branch)
      QINLINE=.FALSE.
C
      END IF
      ELSE
C          ( if not @@ file )
C
C file is read once only and stored in rotating command buffer.
C
C
C get the filename
      IFILE=' '
      CALL NEXTFI('@-FILE=',IFILE)
C
C always read the file even if we're not in "GO" mode
C now assign the file and switch the stream to it
      CALL ASSFIL(IFILE,UNIT,'READ','FORMATTED',ERROR)
C
      IF (.NOT.ERROR) THEN
      REWIND(UNIT=UNIT)
      CALL INSTRM(UNIT)
      QMCHEK=.TRUE.
C create new scope
C
C
      IF (QADD(NSTRM)) THEN
C remove the '@<filename> ...' from the rotating command buffer
C and insert scope-creating directive
C
      BUFLYN(BUFSTK)(I:I+4)=' &%N '
      BUFLEN(BUFSTK)=I+4
      END IF
      CALL DEFNEWSCOPE
      ELSE
C
C always reset the QINLINE flag even if the stream file doesn't
C exist (e.g., within FALSE conditional branch)
      QINLINE=.FALSE.
      END IF
      END IF
C
      END IF
C check to see if file is a MODULE stream
      IF(QMCHEK) THEN
      QMCHEK=.FALSE.
C quietly grab first command from file
      QQECHO=QECHO
      QECHO=.FALSE.
      QQADD=QADD(NSTRM)
      QADD(NSTRM)=.FALSE.
      OLDNSTRM=NSTRM
      UNIT=ISTRM(NSTRM)
      CALL NEXTW2(PROMPT)
      QECHO=QQECHO
      QADD(OLDNSTRM)=QQADD
      IF(UNIT.EQ.ISTRM(NSTRM)) THEN
C rewind file
      REWIND(UNIT=UNIT)
      COMLEN=0
      CURSOR=0
C check for module command
      IF(WD(1:4).EQ.'MACR'.OR.WD(1:4).EQ.'MODU') THEN
      CALL NEXTW2(PROMPT)
C let module subroutine know if module was called from interactive mode
      QTERM = QQTERM
      CALL MODULE
      END IF
C
      ELSE
      CALL SAVEWD
      END IF
      END IF
      IF (.NOT.QSTRAM) CLOOP=.FALSE.
      END DO
C
C do error checking after the INLINE directive
      IF(QINLINE) THEN
      QINLNOK = .FALSE.
C special handline for inline with CALL and CAL_ commands
      IF(WDLEN.EQ.4) THEN
      IF(WD(1:3).EQ.'CAL') THEN
      QINLNOK = .TRUE.
      END IF
      END IF
      IF(.NOT.QINLNOK) THEN
      CALL WRNDIE(-5,'NEXTWD',
     & 'INLINE directive not appropriate in this context')
      QINLINE = .FALSE.
      END IF
      END IF
      RETURN
      END
C=======================================================================
      SUBROUTINE NEXTASS(PROMPT)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER  PROMPT*(*)
C
C begin
      CALL NEXTWD(PROMPT)
      IF (WD(1:4) .EQ. '=   ') CALL NEXTWD(PROMPT)
C
      RETURN
      END
C=============================================================
      SUBROUTINE NEXTFI(PROMPT,FILE)
C
C Gets a filename.  A filename is defined as a sequence of
C non-blank characters.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      CHARACTER*(*) FILE, PROMPT
C local
      DOUBLE PRECISION DOUBLE
      DOUBLE COMPLEX DUCOMP
      INTEGER L
C begin
C
C turn on the filename flag to tell the parser that we're parsing
C a filename.
      QFILNM=.TRUE.
      CALL NEXTW2(PROMPT)
      IF (WD(1:1).EQ.'=') CALL NEXTW2(PROMPT)
C we've got an equal sign, get rid of it.
      IF (WD(1:1).EQ.'?') THEN
C
C if we've got a question mark then this is a inquery
      L=LEN(FILE)
      CALL TRIMM(FILE,L)
C modification, ATB, 12/04/08
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(1X,2A)') PROMPT,FILE(1:L)
      END IF
      CALL DECLAR( 'RESULT', 'ST', FILE(1:L), DUCOMP, DOUBLE )
      ELSE
C
C do machine-dependent filename expansion/modification
      CALL SUBTILDE(WD,WDLEN,WDMAX)
      CALL XXNXTFI(WD,WDLEN,WDMAX)
      CALL COPYST(FILE,LEN(FILE),L,WD,WDLEN)
      END IF
C
C turn off the filename flag
      QFILNM=.FALSE.
      RETURN
      END
C==============================================================
      SUBROUTINE INSTRM(STREAM)
C
C subroutine initializes parsing with NEXTWD on a new
C STREAM channel; NSTRM is incremented by one;
C all necessary information is passed via the comand.inc data
C structure.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER STREAM
C local
      LOGICAL COND
C begin
      COND=NSTRM.EQ.0
      IF (.NOT.COND) COND=(STREAM.NE.ISTRM(NSTRM))
      IF (COND) THEN
      NSTRM=NSTRM+1
      IF (NSTRM.GT.MXSTRM) THEN
      CALL WRNDIE(-5,'INSTRM',
     & 'exceeded MXSTRM (COMAND) parameter --> recompile program')
      ELSE
C
C the last line read on the old stream file is saved in HLDLYN(NSTM-1)
      IF (NSTRM.GT.1) THEN
      IF((COMLEN.GE.CURSOR).AND.(COMLEN.GT.0)) THEN
      CALL COPYST(HLDLYN(NSTRM-1),COMMAX,HLDLEN(NSTRM-1),
     & COMLYN(CURSOR:COMLEN),COMLEN-CURSOR+1)
      ELSE
      HLDLYN(NSTRM-1)=' '
      HLDLEN(NSTRM-1)=1
      END IF
      END IF
C
C fixed ATB 12/4/93
      IF (NSTRM.GT.1) THEN
      QADD(NSTRM)=QADD(NSTRM-1)
      ELSE
      QADD(NSTRM)=.TRUE.
      END IF
      ISTRM(NSTRM)=STREAM
      EOF=.FALSE.
      ERROR=.FALSE.
      QSAVEW=.FALSE.
      QSUBS=.TRUE.
      QFILNM=.FALSE.
      QEXPRS=.FALSE.
      COMLEN=0
      CURSOR=0
C
C finally check whether we are reading from an interactive terminal
      CALL VCHKTT(STREAM,QTERM)
C
      END IF
      END IF
      RETURN
      END
C
      SUBROUTINE INIPRS
C
C Initializes parsing and control registers
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'ctitla.inc'
C local
      INTEGER I
      DOUBLE PRECISION DOUBLE
      DOUBLE COMPLEX DUCOMP
      CHARACTER*1 BACKSL
C begin
C
C initialize the inlining flag to false
      QINLINE = .FALSE.
C
C initialize define feature
      CALL DEFCINIT
      QDEF=.TRUE.
C
C initialize the character mapping array which maps
C lower-case to upper-case characters and maps all
C non-printable characters to spaces
      CALL SETBSL(BACKSL)
      DO I=1,MASCII
      ASCIIM(I)=ICHAR(' ')
      END DO
      ASCIIM(ICHAR(' '))=ICHAR(' ')
      ASCIIM(ICHAR('!'))=ICHAR('!')
      ASCIIM(ICHAR('"'))=ICHAR('"')
      ASCIIM(ICHAR('#'))=ICHAR('#')
      ASCIIM(ICHAR('$'))=ICHAR('$')
      ASCIIM(ICHAR('%'))=ICHAR('%')
      ASCIIM(ICHAR('&'))=ICHAR('&')
      ASCIIM(ICHAR(''''))=ICHAR('''')
      ASCIIM(ICHAR('('))=ICHAR('(')
      ASCIIM(ICHAR(')'))=ICHAR(')')
      ASCIIM(ICHAR('*'))=ICHAR('*')
      ASCIIM(ICHAR('+'))=ICHAR('+')
      ASCIIM(ICHAR(','))=ICHAR(',')
      ASCIIM(ICHAR('-'))=ICHAR('-')
      ASCIIM(ICHAR('.'))=ICHAR('.')
      ASCIIM(ICHAR('/'))=ICHAR('/')
      ASCIIM(ICHAR('0'))=ICHAR('0')
      ASCIIM(ICHAR('1'))=ICHAR('1')
      ASCIIM(ICHAR('2'))=ICHAR('2')
      ASCIIM(ICHAR('3'))=ICHAR('3')
      ASCIIM(ICHAR('4'))=ICHAR('4')
      ASCIIM(ICHAR('5'))=ICHAR('5')
      ASCIIM(ICHAR('6'))=ICHAR('6')
      ASCIIM(ICHAR('7'))=ICHAR('7')
      ASCIIM(ICHAR('8'))=ICHAR('8')
      ASCIIM(ICHAR('9'))=ICHAR('9')
      ASCIIM(ICHAR(':'))=ICHAR(':')
      ASCIIM(ICHAR(';'))=ICHAR(';')
      ASCIIM(ICHAR('<'))=ICHAR('<')
      ASCIIM(ICHAR('='))=ICHAR('=')
      ASCIIM(ICHAR('>'))=ICHAR('>')
      ASCIIM(ICHAR('?'))=ICHAR('?')
      ASCIIM(ICHAR('@'))=ICHAR('@')
      ASCIIM(ICHAR('A'))=ICHAR('A')
      ASCIIM(ICHAR('B'))=ICHAR('B')
      ASCIIM(ICHAR('C'))=ICHAR('C')
      ASCIIM(ICHAR('D'))=ICHAR('D')
      ASCIIM(ICHAR('E'))=ICHAR('E')
      ASCIIM(ICHAR('F'))=ICHAR('F')
      ASCIIM(ICHAR('G'))=ICHAR('G')
      ASCIIM(ICHAR('H'))=ICHAR('H')
      ASCIIM(ICHAR('I'))=ICHAR('I')
      ASCIIM(ICHAR('J'))=ICHAR('J')
      ASCIIM(ICHAR('K'))=ICHAR('K')
      ASCIIM(ICHAR('L'))=ICHAR('L')
      ASCIIM(ICHAR('M'))=ICHAR('M')
      ASCIIM(ICHAR('N'))=ICHAR('N')
      ASCIIM(ICHAR('O'))=ICHAR('O')
      ASCIIM(ICHAR('P'))=ICHAR('P')
      ASCIIM(ICHAR('Q'))=ICHAR('Q')
      ASCIIM(ICHAR('R'))=ICHAR('R')
      ASCIIM(ICHAR('S'))=ICHAR('S')
      ASCIIM(ICHAR('T'))=ICHAR('T')
      ASCIIM(ICHAR('U'))=ICHAR('U')
      ASCIIM(ICHAR('V'))=ICHAR('V')
      ASCIIM(ICHAR('W'))=ICHAR('W')
      ASCIIM(ICHAR('X'))=ICHAR('X')
      ASCIIM(ICHAR('Y'))=ICHAR('Y')
      ASCIIM(ICHAR('Z'))=ICHAR('Z')
      ASCIIM(ICHAR('['))=ICHAR('[')
      ASCIIM(ICHAR(BACKSL))=ICHAR(BACKSL)
      ASCIIM(ICHAR(']'))=ICHAR(']')
      ASCIIM(ICHAR('^'))=ICHAR('^')
      ASCIIM(ICHAR('_'))=ICHAR('_')
      ASCIIM(ICHAR('`'))=ICHAR('`')
      ASCIIM(ICHAR('a'))=ICHAR('A')
      ASCIIM(ICHAR('b'))=ICHAR('B')
      ASCIIM(ICHAR('c'))=ICHAR('C')
      ASCIIM(ICHAR('d'))=ICHAR('D')
      ASCIIM(ICHAR('e'))=ICHAR('E')
      ASCIIM(ICHAR('f'))=ICHAR('F')
      ASCIIM(ICHAR('g'))=ICHAR('G')
      ASCIIM(ICHAR('h'))=ICHAR('H')
      ASCIIM(ICHAR('i'))=ICHAR('I')
      ASCIIM(ICHAR('j'))=ICHAR('J')
      ASCIIM(ICHAR('k'))=ICHAR('K')
      ASCIIM(ICHAR('l'))=ICHAR('L')
      ASCIIM(ICHAR('m'))=ICHAR('M')
      ASCIIM(ICHAR('n'))=ICHAR('N')
      ASCIIM(ICHAR('o'))=ICHAR('O')
      ASCIIM(ICHAR('p'))=ICHAR('P')
      ASCIIM(ICHAR('q'))=ICHAR('Q')
      ASCIIM(ICHAR('r'))=ICHAR('R')
      ASCIIM(ICHAR('s'))=ICHAR('S')
      ASCIIM(ICHAR('t'))=ICHAR('T')
      ASCIIM(ICHAR('u'))=ICHAR('U')
      ASCIIM(ICHAR('v'))=ICHAR('V')
      ASCIIM(ICHAR('w'))=ICHAR('W')
      ASCIIM(ICHAR('x'))=ICHAR('X')
      ASCIIM(ICHAR('y'))=ICHAR('Y')
      ASCIIM(ICHAR('z'))=ICHAR('Z')
      ASCIIM(ICHAR('{'))=ICHAR('{')
      ASCIIM(ICHAR('|'))=ICHAR('|')
      ASCIIM(ICHAR('}'))=ICHAR('}')
      ASCIIM(ICHAR('~'))=ICHAR('~')
C
C initialize the character mapping array which maps
C all single-character words to <=1 (normal mode) or <=2 (expression mode)
C <=0 (filename mode)
      DO I=1,MASCII
      ASCIIS(I)=10
      END DO
      ASCIIS(ICHAR(')'))=1
      ASCIIS(ICHAR('('))=0
      ASCIIS(ICHAR(':'))=1
      ASCIIS(ICHAR('@'))=1
      ASCIIS(ICHAR('='))=0
      ASCIIS(ICHAR('>'))=2
      ASCIIS(ICHAR('<'))=2
      ASCIIS(ICHAR('#'))=2
      ASCIIS(ICHAR('/'))=2
      ASCIIS(ICHAR('*'))=2
      ASCIIS(ICHAR('^'))=2
      ASCIIS(ICHAR(','))=2
      ASCIIS(ICHAR('~'))=2
      ASCIIS(ICHAR('+'))=2
      ASCIIS(ICHAR('-'))=2
      ASCIIS(ICHAR(';'))=1
      ASCIIS(ICHAR(']'))=2
      ASCIIS(ICHAR('['))=2
      ASCIIS(ICHAR(BACKSL))=2
C
C initialize command parsing on unit 5
      NSTRM=0
      CALL INSTRM(5)
      QECHO=.TRUE.
      QCAPIT=.TRUE.
C
C "END" stack
      ENDIND=0
      DO I=1,ENDMAX
      ENDCNT(I)=0
      ENDLBL(I)='MAIN'
      END DO
C
C rotating command buffer
      BUFIND=0
      BUFFIL=0
C
C initialize format for subroutine ENCODF
      CCFORM='(1PG13.6)'
      CCLENG=13
C
C initialize error count
      NERRPA=0
C
C predefined user accessable variable declarations
      CALL DECLAR( 'KBOLTZ', 'DP', ' ', DUCOMP, KBOLTZ )
      CALL DECLAR( 'PI',     'DP', ' ', DUCOMP, PI     )
      CALL DECLAR( 'TIMFAC', 'DP', ' ', DUCOMP, TIMFAC )
      CALL GETNAM( WDT, WDTMAX, WDTLEN )
      CALL DECLAR( 'NAME',   'ST', WDT(1:WDTLEN), DUCOMP, DOUBLE )
      CALL GETSYS( WDT, WDTMAX, WDTLEN )
      CALL DECLAR( 'SYSTEM', 'ST', WDT(1:WDTLEN), DUCOMP, DOUBLE )
      CALL DECLAR( 'LOG_LEVEL', 'ST', 'QUIET', DUCOMP, DOUBLE )
C
C title buffer
      NTITLE=0
      TITMODE=0
C
C initialize output filename
      OFILE='OUTPUT'
      DFILE='OUTPUT'
C
C initalize other variables
      CURSOR=0
      WDLEN=0
      QSCOPING=.FALSE.
C
      RETURN
      END
C
      SUBROUTINE PRIEND
C
C Print information about non-terminated "END" levels.
C
C Author: Axel T. Brunger
C =======================
C I/O
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C local
      INTEGER I
C begin
      IF (ENDIND.GT.0) THEN
      WRITE(6,'(A,I4,A)')
     & ' PRIEND: ',ENDIND,' levels not terminated'
      DO I=1,ENDIND
      WRITE(6,'(A,I4,4A)')
     & '             LEVEL=',I,' KEY=',ENDKEY(I),' ACTION=',ENDACT(I)
      END DO
      ELSE IF (ENDIND.LT.0) THEN
      WRITE(6,'(A,I4)') ' %PRIEND-ERR: level underflow. LEVEL=',ENDIND
      END IF
      RETURN
      END
C
      SUBROUTINE BACEND(PROMPT)
C
C Goes back to where the current END level began using the rotating
C command buffer
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) PROMPT
C local
      INTEGER IOLD
C begin
      IOLD=BUFIND
C
C reset rotating command buffer pointer to the top where the
C current END level began.
      BUFIND=ENDRET(ENDIND)
      BUFSTK=MOD(BUFIND-1,BUFMAX)+1
      IF (IOLD-BUFIND+1.GT.BUFMAX) THEN
      WRITE(6,'(A)')
     & ' %BACEND-ERR: loop too long.  Shorten loop or recompile'
      CALL WRNDIE(-5,'BACEND',
     & 'exceeded BUFMAX (COMAND) parameter --> recompile program')
      ELSE
      CALL COPYST(COMLYN,COMMAX,COMLEN,BUFLYN(BUFSTK),BUFLEN(BUFSTK))
C
      CURSOR=ENDCUR(ENDIND)
      ENDIND=ENDIND-1
C
      CALL DOECHO(PROMPT,COMLYN,COMLEN,.FALSE.,.FALSE.)
C
      END IF
      RETURN
      END
C
      SUBROUTINE PUSACT(ACTION)
C
C Pushes the END action stack
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*4 ACTION
C begin
      IF (ENDIND.EQ.1) THEN
      ENDACT(1)=ACTION
      ELSE IF (ENDACT(ENDIND-1).EQ.'GO  ') THEN
      ENDACT(ENDIND)=ACTION
      ELSE
      ENDACT(ENDIND)='SKIP'
      END IF
      RETURN
      END
C===================================================================
      SUBROUTINE SAVEWD
C
C Routine sets a flag QSAVEW to leave the just parsed word
C in WD when calling NEXTWD again.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C begin
      QSAVEW=.TRUE.
      RETURN
      END
C================================================================
      SUBROUTINE WDSUB(LWD,LWDMAX,LWDLEN,OK,TYPE,DPVAL,DCVAL)
C
C Routine makes substitutions for word LWD beginning with a "$"
C
C I/O descriptions
C
C lwd -> string containing variable input to be substituted for
C lwdmax -> length of lwd character string.
C lwdlen -> number of characters in lwd.
C ok -> returned status of call.
C type -> type of variable that was substituted.
C dpval -> if type was double precision, value is returned here.
C dcval -> if type was double complex, value is returned here.
C (note, the string representation of the variable is always
C returned in lwd, the double complex and double precision
C values are just returned for convienence. )
C
C Authors: Mark McCallum and Axel T. Brunger
C ==========================================
C Modification: Jian-Sheng Jiang and Warren L. DeLano
C
      IMPLICIT  NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      CHARACTER*(*) LWD
      INTEGER LWDLEN, LWDMAX
      LOGICAL OK
      CHARACTER*2 TYPE
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
C local
      INTEGER I, LFORM
      CHARACTER*20 FORM, CFORM
C
C begin
C
C no error yet
      OK = .TRUE.
C
C get the FORMAT specfication if any
      FORM=' '
      I=INDEX(LWD(1:LWDLEN),'[')
      IF (I.GT.0) THEN
      IF (LWD(LWDLEN:LWDLEN).NE.']') THEN
      CALL DSPERR('WDSUB','incorrect FORMAT statement')
      OK = .FALSE.
      ELSE
      FORM=LWD(I+1:LWDLEN-1)
      LFORM=LWDLEN-1-I
      LWDLEN=I-1
C capitalize the format (needed when QCAPIT = .FALSE. )
      DO I=1,LFORM
      FORM(I:I)=CHAR(ASCIIM(ICHAR(FORM(I:I))))
      END DO
C
      END IF
      END IF
C
C substitute any and all "indices" in the symbol name
      IF (OK) THEN
      CALL SYMSYMIDX(LWD,LWDMAX,LWDLEN,.FALSE.)
      END IF
C
C substitute the LWD
      IF (OK) THEN
      CALL WDSUB2(LWD,LWDMAX,LWDLEN,OK,TYPE,DPVAL,DCVAL,
     &            FORM,CFORM,LFORM)
      END IF
C
      RETURN
      END
C
C================================================================
      SUBROUTINE WDSUB2(LWD,LWDMAX,LWDLEN,OK,TYPE,DPVAL,DCVAL,
     &                  FORM,CFORM,LFORM)
C
C See the above
C
C Authors: Mark McCallum and Axel T. Brunger
C ==========================================
C Modification: Jian-Sheng Jiang
C
      IMPLICIT  NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      CHARACTER*(*) LWD
      INTEGER LWDLEN, LWDMAX
      LOGICAL OK
      CHARACTER*2 TYPE
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
      CHARACTER*(*) FORM, CFORM
      INTEGER LFORM
C local
      CHARACTER*(VARMAX) NAME
      INTEGER MATCHLEN, NAMLEN
      LOGICAL QPART
      INTEGER I, J, L, LS, ICHR, OFFSET
      DOUBLE PRECISION SECS
      LOGICAL GOTIT, QDOUBL, QEXIST, VARFLG(NUMVARFLG), QSTRIP, QBLANK
      INTEGER SCOPE
      CHARACTER*1 MARK
      CHARACTER*1 ERRMSG
C parameters
      MARK=ACHAR(30)
C begin
C initialize
      GOTIT = .FALSE.
      OK = .TRUE.
C check if single or double $$
      QDOUBL=LWD(1:2).EQ.'$$'
      IF (QDOUBL) THEN
      LS=3
      ELSE
      LS=2
      END IF
C
C remove any embedded flags in the name...
      NAME = LWD(LS:LWDLEN)
      NAMLEN = LWDLEN - ( LS - 1 )
C
C capitalize the symbol ( needed when QCAPIT = .FALSE. )
      DO I=1,NAMLEN
      NAME(I:I)=CHAR(ASCIIM(ICHAR(NAME(I:I))))
      END DO
C
      CALL DEFGETFLG(NAME,NAMLEN,.TRUE.,OFFSET,VARFLG)
      QEXIST = VARFLG(1)
      QSTRIP = VARFLG(2)
      QBLANK = VARFLG(4)
C
      IF(QEXIST) THEN
C handling of the EXIST% directive
      SCOPE = DEFCURSCOPE
      MATCHLEN = NAMLEN
      CALL DEFFINDREC(NAME,MATCHLEN,
     &  SCOPE,1,.FALSE.,GOTIT,I,QPART,.FALSE.)
      IF(QPART) THEN
      CALL WRNDIE(-5,'SYMIDX',
     & 'can not use the value of a compound symbol as an index')
      GOTIT = .FALSE.
      END IF
C check for special symbols
      IF(.NOT.GOTIT) THEN
C remove scoping information for these special symbols
      CALL DEFREMSCO(NAME,NAMLEN)
      IF (NAME(1:4).EQ.'DATE') GOTIT=.TRUE.
      IF (NAME(1:4).EQ.'TIME') GOTIT=.TRUE.
      IF (NAME(1:3).EQ.'CPU') GOTIT=.TRUE.
      END IF
      IF (GOTIT) THEN
      CALL COPYST(LWD,LWDMAX,LWDLEN,'TRUE',4)
      ELSE
      CALL COPYST(LWD,LWDMAX,LWDLEN,'FALSE',5)
      END IF
      TYPE = 'LO'
      GOTIT=.TRUE.
C
      ELSE IF(QBLANK) THEN
C handling of the BLANK% directive
      SCOPE = DEFCURSCOPE
      MATCHLEN = NAMLEN
      CALL DEFFINDREC(NAME,MATCHLEN,
     &  SCOPE,1,.FALSE.,GOTIT,I,QPART,.FALSE.)
      IF(QPART) THEN
      CALL WRNDIE(-5,'SYMIDX',
     & 'can not use the value of a compound symbol as an index')
      GOTIT = .FALSE.
      END IF
      IF(GOTIT) THEN
C all non-strings are always non-blank
      IF(CMDTYP(I).EQ.'ST') THEN
C check to see if string is blank
      DO J=1,DEFPARLEN(I)
      IF(DEFPARTXT(I)(J:J).NE.' ') GOTIT = .FALSE.
      END DO
      END IF
      ELSE
C remove scoping information for these special symbols
      CALL DEFREMSCO(NAME,NAMLEN)
      IF (NAME(1:4).EQ.'DATE') GOTIT=.TRUE.
      IF (NAME(1:4).EQ.'TIME') GOTIT=.TRUE.
      IF (NAME(1:3).EQ.'CPU') GOTIT=.TRUE.
      END IF
C
      IF (GOTIT) THEN
      CALL COPYST(LWD,LWDMAX,LWDLEN,'TRUE',4)
      ELSE
      CALL COPYST(LWD,LWDMAX,LWDLEN,'FALSE',5)
      END IF
      TYPE = 'LO'
      GOTIT=.TRUE.
      END IF
C end of directive handling...
C
      IF (.NOT.GOTIT) THEN
C
      SCOPE = DEFCURSCOPE
      MATCHLEN = NAMLEN
      CALL DEFFINDREC(NAME,MATCHLEN,
     &  SCOPE,1,.FALSE.,GOTIT,I,QPART,.FALSE.)
C
      IF(QPART.AND.GOTIT) THEN
      CALL WRNDIE(-5,'WDSUB',
     & 'can not obtain value of a compound symbol')
      GOTIT = .FALSE.
      END IF
C
C skip over the $ when comparing
      IF(GOTIT) THEN
C
C check the variable type and copy into LWD
C=============================================================
      IF ( CMDTYP(I) .EQ. 'DP' ) THEN
      IF (FORM.EQ.' ') THEN
C use free-field output
      CALL ENCODF(DPVALU(I),LWD,LWDMAX,LWDLEN)
C MODIFICATION:
C Use hybrid-36 encoding for integer output when the format is
C exactly 'I4' or 'I5', and the value would otherwise overflow.
C ERRMSG is not read, so it is just 1 character.
      ELSE IF (FORM(1:LFORM).EQ.'I4'.AND.INT(DPVALU(I)).GT.9999) THEN
      CALL HY36ENCODE(4,INT(DPVALU(I)),LWD,ERRMSG,J)
      IF (J.NE.0) LWD='****'
      LWDLEN=4
      ELSE IF (FORM(1:LFORM).EQ.'I5'.AND.INT(DPVALU(I)).GT.99999) THEN
      CALL HY36ENCODE(5,INT(DPVALU(I)),LWD,ERRMSG,J)
      IF (J.NE.0) LWD='*****'
      LWDLEN=5
      ELSE
C use formatted output
      CFORM='('//FORM(1:LFORM)//',A1)'
      IF (INDEX(FORM(1:LFORM),'I').NE.0) THEN
      WRITE(LWD,CFORM,ERR=9999) INT(DPVALU(I)),'?'
      ELSE
      WRITE(LWD,CFORM,ERR=9999) DPVALU(I),'?'
      END IF
      LWDLEN=INDEX(LWD,'?')-1
      IF (LWDLEN.LE.0) GOTO 9999
      END IF
      DPVAL = DPVALU(I)
      TYPE = 'DP'
      GOTO 1111
9999  CALL DSPERR('WDSUB','incorrect FORMAT statement')
      OK=.FALSE.
      LWDLEN=1
      LWD='?'
1111  CONTINUE
C=============================================================
      ELSEIF ( CMDTYP(I) .EQ. 'DC' ) THEN
      IF (FORM.EQ.' ') THEN
C use free-field output
      CALL ENCODC(DCVALU(I),LWD,LWDMAX,LWDLEN)
      ELSE
C use formatted output
      CFORM='(2'//FORM(1:LFORM)//',A1)'
      WRITE(LWD,CFORM,ERR=9998) DBLE(DCVALU(I)),DIMAG(DCVALU(I)),'?'
      LWDLEN=INDEX(LWD,'?')-1
      IF (LWDLEN.LE.0) GOTO 9998
      END IF
      DCVAL = DCVALU(I)
      TYPE = 'DC'
      GOTO 1112
9998  CALL DSPERR('WDSUB','incorrect FORMAT statement')
      OK=.FALSE.
      LWDLEN=1
      LWD='?'
1112  CONTINUE
C=============================================================
      ELSEIF ( CMDTYP(I) .EQ. 'ST' ) THEN
      IF (FORM.EQ.' ') THEN
C use free-field output
      CALL COPYST(LWD,LWDMAX,LWDLEN,DEFPARTXT(I),DEFPARLEN(I))
      ELSE
C use formatted output
      CFORM='('//FORM(1:LFORM)//',A1)'
      WRITE(LWD,CFORM,ERR=9997) DEFPARTXT(I),MARK
      LWDLEN=INDEX(LWD,MARK)-1
      IF (LWDLEN.LE.0) GOTO 9997
      END IF
      TYPE = 'ST'
C
C convert into upper case if QDOUBL ($$) or QSTRIP
      IF ((QDOUBL.OR.QSTRIP).AND.QCAPIT) THEN
      DO I=1,LWDLEN
      ICHR=ASCIIM(ICHAR(LWD(I:I)))
      LWD(I:I)=CHAR(ICHR)
      END DO
      END IF
C
      GOTO 1113
9997  CALL DSPERR('WDSUB','incorrect FORMAT statement')
      OK=.FALSE.
      LWDLEN=1
      LWD='?'
1113  CONTINUE
C=============================================================
      ELSEIF ( CMDTYP(I) .EQ. 'LO' ) THEN
C logical
      CALL COPYST(LWD,LWDMAX,LWDLEN,DEFPARTXT(I),DEFPARLEN(I))
      TYPE='LO'
C=============================================================
      ELSE
      CALL DSPERR('WDSUB','Corrupt variable tables')
      OK=.FALSE.
      ENDIF
      GOTIT=.TRUE.
      END IF
C===================================================
C variable was not found in tables, try the remaining possibilities
C
      IF (.NOT.GOTIT.AND.OK) THEN
C remove scoping information from what may be special symbols
      CALL DEFREMSCO(NAME,NAMLEN)
C====================================================
C the question mark is special.  No action is taken here.  Routine
C MISCOM will interpret "$?" and print a list of the symbols
      IF ( NAME(1:1) .EQ. '?' ) THEN
      CALL DEFREMSCO(LWD,LWDLEN)
      CONTINUE
C==================================================
C print date
      ELSE IF (NAME(1:NAMLEN).EQ.'DATE') THEN
      WDT=' '
      CALL VDATE(WDT,WDTMAX,WDTLEN)
      CALL COPYST(LWD,LWDMAX,LWDLEN,WDT,WDTLEN)
      TYPE = 'ST'
C=================================================
      ELSE IF (NAME(1:NAMLEN).EQ.'CPU') THEN
      CALL VCPU(SECS)
      CALL ENCODF(SECS,LWD,LWDMAX,LWDLEN)
      TYPE = 'DP'
      DPVAL=SECS
C=================================================
      ELSE IF (NAME(1:NAMLEN).EQ.'TIME') THEN
      WDT=' '
      CALL VTIME(WDT,WDTMAX,WDTLEN)
      CALL COPYST(LWD,LWDMAX,LWDLEN,WDT,WDTLEN)
      TYPE = 'ST'
C=================================================
      ELSE IF (NAME(1:NAMLEN).EQ.'CURBYTES') THEN
      WDT=' '
      CALL QRYALLOC( 2, WDT)
      L = LEN(WDT)
      CALL TRIML(WDT, L)
      CALL COPYST(LWD,LWDMAX,LWDLEN,WDT,WDTLEN)
      TYPE = 'ST'
C=================================================
      ELSE IF (NAME(1:NAMLEN).EQ.'MAXBYTES') THEN
      WDT=' '
      CALL QRYALLOC( 3, WDT)
      L = LEN(WDT)
      CALL TRIML(WDT, L)
      CALL COPYST(LWD,LWDMAX,LWDLEN,WDT,WDTLEN)
      TYPE = 'ST'
C=================================================
      ELSE IF (NAME(1:NAMLEN).EQ.'CUROVERH') THEN
      WDT=' '
      CALL QRYALLOC(-2, WDT)
      L = LEN(WDT)
      CALL TRIML(WDT, L)
      CALL COPYST(LWD,LWDMAX,LWDLEN,WDT,WDTLEN)
      TYPE = 'ST'
C=================================================
      ELSE IF (NAME(1:NAMLEN).EQ.'MAXOVERH') THEN
      WDT=' '
      CALL QRYALLOC(-3, WDT)
      L = LEN(WDT)
      CALL TRIML(WDT, L)
      CALL COPYST(LWD,LWDMAX,LWDLEN,WDT,WDTLEN)
      TYPE = 'ST'
C=================================================
      ELSE
      OK = .FALSE.
      END IF
      END IF
      END IF
C
C set the type to ST if the string wasn't found.
      IF ( .NOT. OK ) THEN
      TYPE='ST'
      END IF
C
      RETURN
      END
C================================================================
      SUBROUTINE SYMSYMIDX(LWD,LWDMAX,LWDLEN,QSILENT)
C
C Makes substitutions for symbols in symbol name LWD
C Expects LWD to be prefaced with '$', '$$'
C
C LWD and LWDLEN are updated
C
C Author: Warren L. DeLano
C ========================
C
      IMPLICIT  NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) LWD
      INTEGER LWDLEN, LWDMAX
      LOGICAL QSILENT
C local
      CHARACTER*(VARMAX) SYMNAM
      INTEGER SYMLEN
      LOGICAL QDOUBLE
C begin
      IF (LWDLEN.GT.1) THEN
C check for double specification
      QDOUBLE = .FALSE.
      IF(LWDLEN.GT.2) THEN
      IF(LWD(1:2).EQ.'$$') THEN
      QDOUBLE = .TRUE.
      END IF
      END IF
C
C strip off leading $ or $$
      IF(QDOUBLE) THEN
      SYMNAM(1:LWDLEN-2) = LWD(3:LWDLEN)
      SYMLEN = LWDLEN - 2
      ELSE
      SYMNAM(1:LWDLEN) = LWD(2:LWDLEN)
      SYMLEN = LWDLEN - 1
      END IF
C
C make the substitution
      CALL SYMIDX(SYMNAM,VARMAX,SYMLEN,QSILENT)
C
C copy the substituted variable name back into LWD
      IF(QDOUBLE) THEN
      LWD(3:SYMLEN+2) = SYMNAM(1:SYMLEN)
      LWDLEN = SYMLEN + 2
      ELSE
      LWD(2:SYMLEN+1) = SYMNAM(1:SYMLEN)
      LWDLEN = SYMLEN + 1
      END IF
      END IF
C
      RETURN
      END
C================================================================
      SUBROUTINE SYMIDX(LWD,LWDMAX,LWDLEN,QSILENT)
C
C Makes substitutions for symbols in LWD delimited by
C by either an underscore, a period, or the end of the string
C Indexes are indicated by the "$".
C
C If LWD is a symbol or parameter name, any preceeding $ or $$
C should be removed before calling SYMIDX
C
C LWD and LWDLEN are updated
C
C Author: Warren L. DeLano and Axel T. Brunger
C ============================================
C
      IMPLICIT  NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      CHARACTER*(*) LWD
      INTEGER LWDLEN, LWDMAX
      LOGICAL QSILENT
C local
      INTEGER I,K,L,ICHR,SCOPE, OFFSET
      LOGICAL FOUND, QPART, VARFLG(NUMVARFLG)
      CHARACTER*(COMMAX) SYMNAM
      CHARACTER*(COMMAX) NAME
      INTEGER SYMLEN,OLDLEN,RGTLEN, MATCHLEN,NAMLEN
      CHARACTER*2 TYPE
      CHARACTER NXT,CCHR
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
C
C begin
C
      DPVAL=0.0D0
      DCVAL=DCMPLX(0.0D0,0.0D0)
C
      RGTLEN = 0
C find out where the directives end and the name begins
      CALL DEFGETFLG(LWD,LWDLEN,.FALSE.,OFFSET,VARFLG)
      I=OFFSET + 1
      DO WHILE(I.LT.LWDLEN)
C this iterates over all characters
      IF(LWD(I:I).EQ.'$') THEN
C maybe a symbol?
      FOUND = .FALSE.
      K = I
      DO WHILE((.NOT.FOUND).AND.(K.LT.LWDLEN))
      K = K + 1
      CCHR=LWD(K:K)
      ICHR=ASCIIM(ICHAR(CCHR))
C kooky kludge to detect '_$'
      IF(K.LT.LWDLEN) THEN
      NXT=LWD(K+1:K+1)
      ELSE
      NXT=' '
      END IF
C
      IF((ICHR.EQ.ICHAR(' ')).OR.(CCHR.EQ.'.').OR.
     & ((CCHR.EQ.'_').AND.(NXT.EQ.'$')).OR.
     & (ASCIIS(ICHAR(CCHR)).LT.3)) THEN
      K=K-1
      FOUND = .TRUE.
      END IF
      END DO
C
      RGTLEN=LWDLEN-K
C
      SYMLEN = K-I+1
      IF((SYMLEN.LE.0).OR.(SYMLEN.GT.VARMAX)) THEN
      CALL WRNDIE(-5,'SYMIDX',
     & 'Invalid symbol name (too long or non-existent).')
      ELSE
      SYMNAM = LWD(I:K)
C
      SCOPE = DEFCURSCOPE
C
C does it exist?
      IF(SYMLEN.LT.2) THEN
      FOUND=.FALSE.
      ELSE
      NAMLEN = SYMLEN - 1
      NAME = SYMNAM(2:SYMLEN)
C
C capitalize the symbol ( needed when QCAPIT = .FALSE. )
      DO L=1,NAMLEN
      NAME(L:L)=CHAR(ASCIIM(ICHAR(NAME(L:L))))
      END DO
C remove directives
      CALL DEFGETFLG(NAME,NAMLEN,.TRUE.,OFFSET,VARFLG)
      MATCHLEN = NAMLEN
C see if symbol exists
      CALL DEFFINDREC(NAME,MATCHLEN,
     &  SCOPE,1,.FALSE.,FOUND,L,QPART,.FALSE.)
      IF(QPART) THEN
      FOUND=.FALSE.
      END IF
      END IF
C
      IF(.NOT.FOUND.AND..NOT.QSILENT) THEN
      WRITE(6,'(A,A,A)') ' %SYMIDX-WRN: index-symbol ',
     & SYMNAM(1:SYMLEN),' not found.'
      ELSE
C replace the name with the symbol's value
      OLDLEN = SYMLEN
      FOUND = .TRUE.
      CALL WDSUB2(SYMNAM,COMMAX,SYMLEN,FOUND,TYPE,DPVAL,DCVAL,
     &            ' ',' ',0)
C check to make sure the value is acceptable
      IF(FOUND.AND.(TYPE.EQ.'DP')) THEN
      DO L=1,SYMLEN
      IF (INDEX('0123456789',SYMNAM(L:L)).EQ.0) FOUND=.FALSE.
      END DO
      IF(.NOT.FOUND) THEN
      CALL WRNDIE(-5,'SYMIDX',
     & 'Illegal symbol value, must be integer or character string.')
      END IF
      END IF
C
      IF(FOUND) THEN
      IF((TYPE.EQ.'DC').OR.(TYPE.EQ.'LO')) THEN
      CALL WRNDIE(-5,'SYMIDX',
     & 'Illegal symbol type, must be integer or character string.')
      ELSE IF((TYPE.EQ.'DP').AND.(DPVAL.LT.0.0)) THEN
      CALL WRNDIE(-5,'SYMIDX',
     & 'Illegal symbol value (negative not allowed).')
      ELSE IF (     (TYPE.EQ.'DP')
     &         .AND.(DPVAL-DBLE(INT(DPVAL)).GT.RSMALL)) THEN
      CALL WRNDIE(-5,'SYMIDX',
     & 'Illegal symbol value, must be integer or character string.')
      ELSE
C value checks out okay
C
C capitalize string
      DO L=1,SYMLEN
      SYMNAM(L:L)=CHAR(ASCIIM(ICHAR(SYMNAM(L:L))))
      END DO
C value checks out okay - so insert the new symbol
      IF((LWDLEN+(SYMLEN-OLDLEN)).LT.LWDMAX) THEN
      IF(RGTLEN.GT.0) THEN
      SYMNAM(SYMLEN+1:SYMLEN+RGTLEN) = LWD(K+1:LWDLEN)
      SYMLEN=SYMLEN+RGTLEN
      END IF
      LWD(I:I+SYMLEN-1) = SYMNAM(1:SYMLEN)
      LWDLEN = I + SYMLEN - 1
      ELSE
      CALL WRNDIE(-5,'SYMIDX',
     & 'Symbol too long.')
      END IF
      END IF
      END IF
C            (.not.found)
C
      END IF
      END IF
      END IF
      I = I + 1
      END DO
      RETURN
      END
C====================================================================
      SUBROUTINE LINSUB(L,LMAX,LLEN,PROMPT,LPROMP,QMODE)
C
C If LPROMP is greater zero, string gets initialized with PROMPT.
C
C Routine then appends the rest of the COMLYN line starting at
C CURSOR+1 into string L.   Routine makes symbol substitutions
C unless the symbols are in quoted strings.  Routine also
C moves CURSOR to the end of the line.
C
C QMODE = 1 normal - do substitutions
C QMODE = 2 literal - don't do substitutions
C QMODE = 3 "scan" mode - just process internal scoping directives
C
C Author: Axel T. Brunger
C =======================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) L
      INTEGER LMAX, LLEN
      CHARACTER*(*) PROMPT
      INTEGER LPROMP
      INTEGER QMODE
C local
      CHARACTER*1 BACKSL
      LOGICAL QQSUBS, QQCAPIT
      LOGICAL QQDEF
      LOGICAL CLOOP
C begin
      QQSUBS = QSUBS
      QQDEF = QDEF
      QQCAPIT = QCAPIT
      QCAPIT = .FALSE.
      IF (QMODE.GT.1) THEN
      QDEF = .FALSE.
      QSUBS = .FALSE.
      END IF
      CALL SETBSL(BACKSL)
C
C add the prompt to the output string
      IF(QMODE.LT.3) CALL COPYST(L,LMAX,LLEN,PROMPT,LPROMP)
C
C skip the first character
      CURSOR=CURSOR+1
      IF (CURSOR.LE.COMLEN) THEN
C
C go through all characters of the COMLYN line
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
C
C NEXTW2 takes care of symbol substitions and quoted strings.
C
      IF (QMODE.EQ.1) THEN
C
      IF (COMLYN(CURSOR:CURSOR).EQ.'&') THEN
C define and MODULE parameter substitutions
      CALL DEFSUBPAR(COMLYN,COMMAX,COMLEN,CURSOR,.TRUE.)
      END IF
      ELSE
C literal and/or scanning mode
      IF (COMLYN(CURSOR:CURSOR).EQ.'&') THEN
      IF (QMODE.LT.3) THEN
      CALL DEFSUBPAR(COMLYN,COMMAX,COMLEN,CURSOR,.TRUE.)
      ELSE IF(COMLEN.GT.CURSOR) THEN
      IF (COMLYN(CURSOR+1:CURSOR+1).EQ.'%') THEN
      CALL DEFSUBPAR(COMLYN,COMMAX,COMLEN,CURSOR,.TRUE.)
      END IF
      END IF
      END IF
C
      END IF
C
      IF (COMLYN(CURSOR:CURSOR).EQ.'$') THEN
C symbol substitutions
C must retrieve the format specification []
      ASCIIS(ICHAR('['))=10
      ASCIIS(ICHAR(']'))=10
      QEXPRS=.TRUE.
      CALL NEXTW2(' ')
      QEXPRS=.FALSE.
      ASCIIS(ICHAR('['))=2
      ASCIIS(ICHAR(']'))=2
      IF(QMODE.LT.3) CALL ADDST(L,LMAX,LLEN,WD,WDLEN)
C quoted strings
      ELSE IF (COMLYN(CURSOR:CURSOR).EQ.'"') THEN
      IF(QMODE.LT.3) THEN
      CALL ADDST(L,LMAX,LLEN,'"',1)
      END IF
      CALL NEXTW2(' ')
      IF(QMODE.LT.3) THEN
      CALL ADDST(L,LMAX,LLEN,WD,WDLEN)
      CALL ADDST(L,LMAX,LLEN,'"',1)
      END IF
C continuation lines
      ELSE IF (CURSOR.EQ.COMLEN-1
     &        .AND.COMLYN(CURSOR:CURSOR).EQ.BACKSL) THEN
      COMLYN(CURSOR:CURSOR) = ' '
      CURSOR=COMLEN
      ASCIIS(ICHAR('['))=10
      ASCIIS(ICHAR(']'))=10
      QEXPRS=.TRUE.
      CALL NEXTW2(' ')
      QEXPRS=.FALSE.
      ASCIIS(ICHAR('['))=2
      ASCIIS(ICHAR(']'))=2
      CURSOR=1
      ELSE
      IF(QMODE.LT.3)
     & CALL ADDST(L,LMAX,LLEN,COMLYN(CURSOR:CURSOR),1)
      CURSOR=CURSOR+1
      END IF
C
      IF (CURSOR.GE.COMLEN) CLOOP=.FALSE.
      END DO
C
      END IF
C
      QSUBS = QQSUBS
      QDEF = QQDEF
      QCAPIT = QQCAPIT
      RETURN
      END
C=============================================================
      SUBROUTINE NEXTF(PROMPT,VALUE)
C
C Routine gets next word and converts it into a floating point real
C number. In case of an ERROR in DECODF the ERROR flag in comand.inc is
C set.
C The routine also picks up a sign in front of the constant.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      CHARACTER*(*) PROMPT
      DOUBLE PRECISION VALUE
C local
      LOGICAL FOUND
      EXTERNAL DECODF
      DOUBLE PRECISION DECODF
      DOUBLE COMPLEX DCVAL
C begin
      CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'?   ') THEN
C modification, ATB, 12/04/08
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(1X,A,G12.5)') PROMPT,VALUE
      END IF
      CALL DECLAR('RESULT', 'DP' ,' ', DCVAL, VALUE )
      ELSE
      VALUE=DECODF(WD,WDLEN,FOUND)
      ERROR=.NOT.FOUND
      IF (ERROR) CALL DSPERR(PROMPT,'real number expected')
      END IF
      RETURN
      END
C
      SUBROUTINE NEXTI(PROMPT,VALUE)
C
C Routine gets next word and converts it into an integer
C number. In case of an ERROR in DECODI the ERROR flag in comand.inc is
C set.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      CHARACTER*(*) PROMPT
      INTEGER VALUE
C local
      EXTERNAL DECODI
      INTEGER DECODI
      LOGICAL FOUND
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
C begin
      CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'?   ') THEN
C modification, ATB, 12/04/08
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(1X,A,I8)') PROMPT,VALUE
      END IF
      DPVAL=VALUE
      CALL DECLAR('RESULT', 'DP' ,' ', DCVAL, DPVAL )
      ELSE
      VALUE=DECODI(WD,WDLEN,FOUND)
      ERROR=.NOT.FOUND
      IF (ERROR) CALL DSPERR(PROMPT,'integer number expected')
      END IF
      RETURN
      END
C
      SUBROUTINE NEXTIT(PROMPT,VALUE)
C
C Routine gets next word, reads it as a real and truncates it into
C an integer number. In case of an ERROR in DECODF the ERROR flag in
C comand.inc is set.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      CHARACTER*(*) PROMPT
      INTEGER VALUE
C local
      EXTERNAL DECODF
      DOUBLE PRECISION DECODF
      LOGICAL FOUND
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
C begin
      CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'?   ') THEN
C modification, ATB, 12/04/08
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(1X,A,I10)') PROMPT,VALUE
      END IF
      DPVAL=VALUE
      CALL DECLAR('RESULT', 'DP' ,' ', DCVAL, DPVAL )
      ELSE
      VALUE=INT(DECODF(WD,WDLEN,FOUND))
      ERROR=.NOT.FOUND
      IF (ERROR) CALL DSPERR(PROMPT,'integer or real number expected')
      END IF
      RETURN
      END
C
      SUBROUTINE NEXTLO(PROMPT,VALUE)
C
C Routine gets next word and converts it into a logical value.
C In case of an error the ERROR flag in comand.inc is set.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      CHARACTER*(*) PROMPT
      LOGICAL VALUE
C local
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
C begin
      CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'?   ') THEN
      IF (VALUE) THEN
C modification, ATB, 12/04/08
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(1X,A,A)') PROMPT,'TRUE {ON}'
      END IF
      CALL DECLAR('RESULT', 'LO' ,'TRUE', DCVAL, DPVAL )
      ELSE
C modification, ATB, 12/04/08
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(1X,A,A)') PROMPT,'FALSe {OFF}'
      END IF
      CALL DECLAR('RESULT', 'LO' ,'FALSE', DCVAL, DPVAL )
      END IF
      ELSE
      IF (WD(1:WDLEN).EQ.'TRUE'.OR.WD(1:WDLEN).EQ.'ON'
     & .OR.WD(1:WDLEN).EQ.'T') THEN
      VALUE=.TRUE.
      ELSE IF (WD(1:4).EQ.'FALS'.OR.WD(1:WDLEN).EQ.'OFF'
     & .OR.WD(1:WDLEN).EQ.'F') THEN
      VALUE=.FALSE.
      ELSE
      CALL DSPERR(PROMPT,
     & 'ON (True) or OFF (False) expected')
      END IF
      END IF
      RETURN
      END
C
      SUBROUTINE NEXTA4(PROMPT,VALUE)
C
C Routine gets next word and returns the first four characters.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      CHARACTER*(*) PROMPT
C local
      CHARACTER*4 VALUE
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
C begin
      CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'?   ') THEN
C modification, ATB, 12/04/08
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(1X,2A)') PROMPT,VALUE
      END IF
      CALL DECLAR('RESULT', 'ST' ,VALUE, DCVAL, DPVAL )
      ELSE
      VALUE=WD(1:4)
      END IF
      RETURN
      END
C
      SUBROUTINE NEXTQL(PROMPT)
C
C Same as NEXTWD expect that the special set of single-character
C words are used during parsing.  Symbol substitutions are performed.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) PROMPT
C begin
      QEXPRS=.TRUE.
      CALL NEXTWD(PROMPT)
      QEXPRS=.FALSE.
      RETURN
      END
C
      SUBROUTINE NEXTSL(PROMPT)
C
C Same as NEXTWD except that the special set of single-character
C words are used during parsion and no symbol substitutions are
C performed except when symbols are referenced as $$ and the symbol
C type is "string".
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) PROMPT
C local
      CHARACTER*2 TYPE
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
      INTEGER I, ICHR, OFFSET
      LOGICAL OK, VARFLG(NUMVARFLG), QUPPER
      CHARACTER*(WDMAX) TMPWD
      INTEGER TMPLEN
C begin
      QEXPRS=.TRUE.
      QSUBS=.FALSE.
      QUPPER = .FALSE.
      CALL NEXTWD(PROMPT)
C
C check if symbol is referenced as $$ and is of type string
C
      IF (WD(1:2).EQ.'$$') THEN
      CALL COPYST(WDD,WDDMAX,WDDLEN,WD,WDLEN)
      CALL WDSUB(WDD,WDDMAX,WDDLEN,OK,TYPE,DPVAL,DCVAL)
      IF (TYPE.EQ.'ST'.AND.OK) QUPPER = .TRUE.
C
C check if symbol has STRIP% directive and is of type string
      ELSE IF(WD(1:1).EQ.'$') THEN
      IF(WDLEN.GT.1) THEN
      TMPLEN = WDLEN - 1
      TMPWD(1:TMPLEN) = WD(2:WDLEN)
      CALL DEFGETFLG(TMPWD,TMPLEN,.FALSE.,OFFSET,VARFLG)
      IF(VARFLG(2)) THEN
C has STRIP% directive
      CALL COPYST(WDD,WDDMAX,WDDLEN,WD,WDLEN)
      CALL WDSUB(WDD,WDDMAX,WDDLEN,OK,TYPE,DPVAL,DCVAL)
      IF (TYPE.EQ.'ST'.AND.OK) QUPPER = .TRUE.
      END IF
      END IF
      END IF
C
      IF(QUPPER) THEN
      QQUOT=.FALSE.
C
C convert into upper case
      DO I=1,WDDLEN
      ICHR=ASCIIM(ICHAR(WDD(I:I)))
      WDD(I:I)=CHAR(ICHR)
      END DO
      CALL COPYST(WD,WDMAX,WDLEN,WDD,WDDLEN)
      END IF
C
      QSUBS=.TRUE.
      QEXPRS=.FALSE.
      RETURN
      END
C
      SUBROUTINE NEXTST(PROMPT,VALUE)
C
C Get new STRING after Blank or =
C routine gets next word, if it is an "="
C get the next word.
C
C Author Axel T. Brunger
C ========================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) PROMPT,VALUE
C local
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
C begin
      CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(1X,2A)')PROMPT,VALUE
      CALL DECLAR('RESULT', 'ST' ,WD(1:WDLEN), DCVAL, DPVAL )
      ELSE
      VALUE=WD
      END IF
      RETURN
      END
C
      SUBROUTINE NEXTEX(PROMPT,ST,STMAX,STLEN)
C
C Gets a the next stream of characters which are
C enclosed by "(", ")", and puts it into ST.
C
C Special treatment for the case that WD(1:WDLEN)
C is a string where the first character is "("
C and the last character is ")".  (ATB 8/28/94)
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) PROMPT, ST
      INTEGER STLEN, STMAX
C local
      INTEGER I
C begin
      CALL NEXTST(PROMPT,WD)
      IF (WD(1:1).NE.'(') THEN
      CALL DSPERR(PROMPT,'"(" expected')
      ELSE
      I=1
C
      CALL COPYST(ST,STMAX,STLEN,WD(1:WDLEN),WDLEN)
C
C We're all set when WD(1:WDLEN) is a string
C with the first character "(" and the last character
C ")".  Otherwise, we have to keep parsing until
C we get a ")" that closes the first "(".
C
      IF (WD(WDLEN:WDLEN).NE.')') THEN
      DO WHILE (I.NE.0)
      CALL NEXTQL(PROMPT)
      IF (WD(1:1).EQ.'(') I=I+1
      IF (WD(1:1).EQ.')') I=I-1
      CALL ADDST(ST,STMAX,STLEN,WD(1:WDLEN),WDLEN)
      END DO
C
      END IF
C
      END IF
      RETURN
      END
C
      SUBROUTINE NEXTVF(PROMPT,VECTOR)
C
C parses a 3-dimensional vector
C
C For syntax see main parsing loop "HELP"
C
C Commas made optional (ATB 8/28/94)
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      DOUBLE PRECISION VECTOR(*)
C local
      INTEGER SELCT
      CHARACTER*(*) PROMPT
C dynamic
      SELCT=ALLHP(INTEG4(NATOM))
      CALL NEXTV2(VECTOR,PROMPT,HEAP(SELCT))
      CALL FREHP(SELCT,INTEG4(NATOM))
      RETURN
      END
C
      SUBROUTINE NEXTV2(VECTOR,PROMPT,SELCT)
C
C see routine NEXTVF above
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      DOUBLE PRECISION VECTOR(*)
      CHARACTER*(*) PROMPT
      INTEGER SELCT(*)
C local
      DOUBLE PRECISION TEMP(3), TMASS
      INTEGER I, J, NSELCT, SIGN
      LOGICAL OK, CLOOP
C begin
      ERROR=.FALSE.
C
      CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(1X,2A,3G12.5,A)') PROMPT,'(',(VECTOR(I),I=1,3),')'
      ELSE
      DO I=1,3
      VECTOR(I)=0.0D0
      END DO
C
      SIGN=1
C
      IF (WD(1:1).EQ.'(') THEN
      I=0
      CALL NEXTQL(PROMPT)
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
C
C commas are optional
      IF (WD(1:4).EQ.',') THEN
C
C sign
      ELSEIF (WD(1:WDLEN).EQ.'-') THEN
      SIGN=-1
      ELSEIF (WD(1:WDLEN).EQ.'+') THEN
      SIGN=+1
C
C TAIL
      ELSEIF (WD(1:4).EQ.'TAIL'.AND.I.EQ.0) THEN
      CALL FILL4(SELCT,NATOM,0)
      NSELCT=0
      CALL SELCTA(SELCT,NSELCT,X,Y,Z,.TRUE.)
      IF (NSELCT.GT.0) THEN
      DO J=1,3
      TEMP(J)=0.0D0
      END DO
      TMASS=0.0D0
      DO J=1,NATOM
      IF (SELCT(J).EQ.1) THEN
      TEMP(1)=TEMP(1)+X(J)*AMASS(J)
      TEMP(2)=TEMP(2)+Y(J)*AMASS(J)
      TEMP(3)=TEMP(3)+Z(J)*AMASS(J)
      TMASS=TMASS+AMASS(J)
      END IF
      END DO
      DO J=1,3
      VECTOR(J)=VECTOR(J)-TEMP(J)/TMASS
      END DO
      END IF
C
C HEAD
      ELSE IF (WD(1:4).EQ.'HEAD'.AND.I.EQ.0) THEN
      CALL FILL4(SELCT,NATOM,0)
      NSELCT=0
      CALL SELCTA(SELCT,NSELCT,X,Y,Z,.TRUE.)
      IF (NSELCT.GT.0) THEN
      DO J=1,3
      TEMP(J)=0.0D0
      END DO
      TMASS=0.0D0
      DO J=1,NATOM
      IF (SELCT(J).EQ.1) THEN
      TEMP(1)=TEMP(1)+X(J)*AMASS(J)
      TEMP(2)=TEMP(2)+Y(J)*AMASS(J)
      TEMP(3)=TEMP(3)+Z(J)*AMASS(J)
      TMASS=TMASS+AMASS(J)
      END IF
      END DO
      DO J=1,3
      VECTOR(J)=VECTOR(J)+TEMP(J)/TMASS
      END DO
      END IF
C
C parse the x,y,z form of <vector>
      ELSE
      IF (I.GE.3) THEN
      ERROR=.TRUE.
      CALL DSPERR(PROMPT,'dimension of vector greater 3')
      ELSE
      I=I+1
      VECTOR(I)=SIGN*DECODF(WD,WDLEN,OK)
      SIGN=1
      ERROR=.NOT.OK
      IF (ERROR) THEN
      CALL DSPERR(PROMPT,'numerical value expected')
      END IF
      END IF
      END IF
      CALL NEXTQL(PROMPT)
      IF (WD(1:1).EQ.')'.OR.ERROR) CLOOP=.FALSE.
      END DO
C
      ELSE
      CALL DSPERR(PROMPT,'vector statement has to start with "("')
      END IF
      END IF
      RETURN
      END
C
      SUBROUTINE MISCOM(PROMPT,LUSED)
C
C Routine parses auxiliary statements
C
C For syntax see main parsing loop "HELP"
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'ctitla.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'seed.inc'
      INCLUDE 'consta.inc'
      CHARACTER*(*) PROMPT
      LOGICAL   LUSED
C local
      LOGICAL QOPEN, OK, QFORM, QWRITE, QTEMP, QEXIST, ABORT
      INTEGER UNIT, I, L, SIGNUM, WDTYPL, VLEN
      DOUBLE PRECISION TEMP
      CHARACTER*(WDMAX+4) WDTYP
      CHARACTER*6 DISPOS
      CHARACTER*12 FORM, ACCESS
      CHARACTER*4 STEMP
      CHARACTER*2 FVARTP
      CHARACTER*(WDMAX) CNSVERSION, FILEVERSION
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
      LOGICAL CLOOP
C
      CHARACTER*(COMMAX) COMCPY
      INTEGER COMCPYLEN
      CHARACTER*(WDDMAX) WDDCPY
      INTEGER WDDCPYLEN
      CHARACTER*(WDMAX) WDCPY
      INTEGER WDCPYLEN
      CHARACTER*(WDMAX+4) WDTYPCPY
      INTEGER WDTYPCPYL
      LOGICAL QDECLARE
C
C begin
      LUSED=.TRUE.
C
C initialize the symbol delaration flag for for loops
      QDECLARE = .FALSE.
C
      IF (WD(1:8).EQ.'HELP-GLO') THEN
C
      CALL CNSHELP('global')
C
      ELSE IF (WD(1:8).EQ.'HELP-DAT') THEN
C
      CALL CNSHELP('data-types')
C
      ELSE IF (WD(1:4).EQ.'SYST') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      CALL SYTASK
      END IF
C ----------------------------------------------------------------------
C process-remarks-command
C -----------------------------
      ELSE IF (WD(1:4).EQ.'REMA') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      IF (TITMODE.EQ.0) NTITLE=0
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      NTITLE=NTITLE+1
C
C leave space (5 lines) for REMARKs to be added if a file is written out
      IF (NTITLE.LE.(MXTITL-5)) THEN
C
C run it through LINSUB to make symbol substitutions
      CALL LINSUB(TITLE(NTITLE),TITMAX,L,' REMARKS ',9,1)
      IF (QTERM) WRITE(6,'(A)') TITLE(NTITLE)(1:L)
      ELSE
      NTITLE=NTITLE-1
      IF (WRNLEV.GE.15) THEN
      WRITE(6,'(A,A)')' MISCOM: exceeded maximum number of REMARKs ',
     &'- ignoring REMARK'
      END IF
      END IF
      CURSOR=COMLEN
      CALL NEXTWD(PROMPT)
      IF (.NOT.(WD(1:4).EQ.'TITL'.OR.WD(1:4).EQ.'REMA')) CLOOP=.FALSE.
      END DO
      CALL SAVEWD
      ELSE
      CALL LINSUB(' ',0,0,' ',0,3)
      CURSOR=COMLEN
      END IF
C
C ----------------------------------------------------------------------
C process-display-command
C -----------------------
      ELSE IF ((WD(1:7).EQ.'DISPLAY').OR.
     &         (WD(1:7).EQ.'LITERAL')) THEN
      IF (ENDACT(ENDIND).NE.'GO  ') THEN
      CALL LINSUB(' ',0,0,' ',0,3)
      ELSE
C
C run it through LINSUB to make symbol substutions
      IF(WD(1:7).EQ.'DISPLAY') THEN
      CALL LINSUB(DISLYN,DISMAX,DISLEN,' ',0,1)
      ELSE
      CALL LINSUB(DISLYN,DISMAX,DISLEN,' ',0,2)
      END IF
      WRITE(DUNIT,'(A)',ERR=8888) DISLYN(1:DISLEN)
      GOTO 9999
8888  CALL WRNDIE(-1,'MISCOM','error writing on display unit.')
9999  CONTINUE
      END IF
      CURSOR=COMLEN
C
C ----------------------------------------------------------------------
C process-$?-command
C ------------------
      ELSE IF (WD(1:2).EQ.'$?') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      CALL DEFVARDUMP(WD,WDLEN,3,1)
      END IF
C------------------------------------------------------------------
      ELSE IF (WD(1:6).EQ.'CRYST1') THEN
C
C ignore this info for now
      CURSOR=COMLEN
C----------------------------------------------------------------------
      ELSE IF (WD(1:2).EQ.'&?') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      CALL DEFVARDUMP(WD,WDLEN,3,2)
      END IF
C----------------------------------------------------------------------
      ELSE IF (WD(1:5).EQ.'ABORT') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      CALL WRNDIE(-5,' ','ABORT statement specified.')
      END IF
C----------------------------------------------------------------------
      ELSE IF (WD(1:6).EQ.'CHECKV') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      CALL NEXTST('VERSION=',FILEVERSION)
      CNSVERSION='$CNS_VERSION'
      WDLEN=12
      CALL WDSUB(CNSVERSION,WDMAX,WDLEN,OK,FVARTP,DPVAL,DCVAL)
      VLEN=WDMAX
      CALL TRIMM(FILEVERSION,VLEN)
C modification, ATB, 12/04/08
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(4A)')' Program version= ',CNSVERSION(1:WDLEN),
     &                   ' File version= ',FILEVERSION(1:VLEN)
      END IF
      IF (FILEVERSION.NE.CNSVERSION) THEN
      CALL WRNDIE(-5,'MISCOM',
     & 'Version numbers do not match - aborting program.')
      END IF
      END IF
C----------------------------------------------------------------------
C process-set-commands
C --------------------
      ELSE IF (WD(1:4).EQ.'SET ') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
C
C defaults in main program
C
C parsing
      CALL PUSEND('SET>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('SET>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-set')
C
C==================================================================
      ELSE IF (WD(1:4).EQ.'ECHO') THEN
      CALL NEXTLO('ECHO=',QECHO)
C==================================================================
      ELSE IF (WD(1:4).EQ.'TIMI') THEN
      IF (TIMER.GT.0) THEN
      QTEMP=.TRUE.
      ELSE
      QTEMP=.FALSE.
      END IF
      CALL NEXTLO('TIMIng=',QTEMP)
      IF (QTEMP) THEN
      TIMER=2
      ELSE
      TIMER=0
      END IF
C==================================================================
      ELSE IF (WD(1:4).EQ.'REMA') THEN
      CALL NEXTA4('REMArks_mode=',STEMP)
      IF (STEMP.EQ.'RESE') NTITLE=0
      IF (STEMP.EQ.'ACCU') TITMODE=1
      IF (STEMP.EQ.'AUTO') TITMODE=0
C==================================================================
      ELSE IF (WD(1:4).EQ.'MESS') THEN
      IF (WRNLEV.EQ.0) THEN
      STEMP='OFF '
      ELSE IF (WRNLEV.EQ.10) THEN
      STEMP='ALL '
      ELSE IF (WRNLEV.EQ.15) THEN
      STEMP='DEBU'
      ELSE
      STEMP='NORM'
      END IF
      CALL NEXTA4('MESSage=',STEMP)
      IF (STEMP.EQ.'OFF ') THEN
      WRNLEV=0
      ELSE IF (STEMP.EQ.'ALL ') THEN
      WRNLEV=10
      ELSE IF (STEMP.EQ.'DEBU') THEN
      WRNLEV=15
      ELSE
      WRNLEV=5
      END IF
C==================================================================
      ELSE IF (WD(1:4).EQ.'ABOR') THEN
      IF (BOMLEV.LE.-5) THEN
      STEMP='OFF '
      ELSE IF (BOMLEV.EQ.5) THEN
      STEMP='ALL '
      ELSE
      STEMP='NORM'
      END IF
      CALL NEXTA4('ABORt=',STEMP)
      IF (STEMP.EQ.'OFF ') THEN
      BOMLEV=-6
      ELSE IF (STEMP.EQ.'ALL ') THEN
      BOMLEV=5
      ELSE
      BOMLEV=0
      END IF
C==================================================================
      ELSE IF (WD(1:4).EQ.'INTE') THEN
      CALL NEXTLO('INTEractive-mode=',QTERM)
C==================================================================
      ELSE IF (WD(1:4).EQ.'DISP') THEN
      CALL NEXTFI('DISPlay-file=',DFILE)
      CALL ASSFIL(DFILE,DUNIT,'WRITE','FORMATTED',ERROR)
C==================================================================
      ELSE IF (WD(1:4).EQ.'PRIN') THEN
      CALL NEXTFI('PRINt-file=',OFILE)
      CALL ASSFIL(OFILE,PUNIT,'WRITE','FORMATTED',ERROR)
C==================================================================
      ELSE IF (WD(1:4).EQ.'JOUR') THEN
      CALL NEXTFI('JOUrnal-file=',OFILE)
      CALL ASSFIL(OFILE,PRUNIT,'WRITE','FORMATTED',ERROR)
C==================================================================
      ELSE IF (WD(1:4).EQ.'PREC') THEN
      CCLENG=CCLENG-7
      CALL NEXTI('PRECision=',CCLENG)
      CCLENG=CCLENG+7
      WRITE(CCFORM,'(A,I2,A,I2,A)') '(1PG',CCLENG,'.',CCLENG-7,')'
C==================================================================
      ELSE IF (WD(1:4).EQ.'SEED') THEN
      TEMP=SEED
      CALL NEXTF('SEED=',TEMP)
      SEED=TEMP
C==================================================================
C
C compatibility commands (obsolete)
      ELSE IF (WD(1:4).EQ.'TIME') THEN
      CALL NEXTI('TIMEr=',TIMER)
      ELSE IF (WD(1:4).EQ.'WRNL') THEN
      CALL NEXTI('WRNLevel=',WRNLEV)
      ELSE IF (WD(1:4).EQ.'BOMB') THEN
      CALL NEXTI('BOMBlevel=',BOMLEV)
C==================================================================
      ELSE
      CALL CHKEND('SET>',DONE)
      END IF
      END DO
      DONE=.FALSE.
      END IF
C
C ------------------------------------------------------
C inline module and procedure handling - just set a flag
C ------------------------------------------------------
      ELSE IF (WD(1:WDLEN).EQ.'INLINE') THEN
      QINLINE = .TRUE.
C ---------------------------------------------------------------------
C MODULE parsing in the rotating command buffer
C ------------------------------------------------------------
      ELSE IF (WD(1:WDLEN).EQ.'MACRO'.OR.WD(1:WDLEN).EQ.'MODULE') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
C MODULE parsing in the rotating command buffer
C make sure the next word is a parenthesis
      CALL DEFNEXTW2(2)
C check for problems
      IF((WD(1:1).NE.'(').OR.(WDLEN.NE.1)) THEN
      WRITE (6,'(A)') ' %MODULE-ERR: ( expected but not found.'
      ERROR=.TRUE.
      ELSE
C specify scope/stream into which defines will be stored
      DEFSETSCOPE = DEFCURSCOPE
      DEFEVLSCOPE = DEFCURSCOPE
C parse parameter block
      CALL DEFMACSET(ABORT,ERROR,2)
      END IF
C
      IF(.NOT.ERROR.AND..NOT.ABORT) THEN
C now check for second set of parenthesis
      CALL DEFNEXTW2(1)
C check for problems
      IF((WD(1:1).NE.'(').OR.(WDLEN.NE.1)) THEN
      WRITE (6,'(A)') ' %MODULE-ERR: ( expected but not found.'
      ERROR=.TRUE.
C parse invocation parameter block
      ELSE
      DEFEVLSCOPE = MAX(1,DEFCURSCOPE-1)
      CALL DEFMACSET(ABORT,ERROR,1)
      END IF
      END IF
C
      IF(ERROR) THEN
C MODULE invocation errored
      CALL WRNDIE(-5,'MODULE','module expansion failed.')
      END IF
      END IF
C ---------------------------------------------------------------------
      ELSE IF (WD(1:WDLEN).EQ.'DEFINE') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
C make sure the next word is a parenthesis
      CALL NEXTW2('DEFINE>')
      IF((WD(1:1).NE.'(').OR.(WDLEN.NE.1)) THEN
      WRITE (6,'(A)') ' %DEFINE-ERR: ( expected but not found.'
      ERROR=.TRUE.
      ELSE
C specify scope/stream into which defines will be stored
      CALL DEFGETSSC(DEFSETSCOPE,DEFCURSCOPE)
      DEFEVLSCOPE = DEFCURSCOPE
C parse parameter bloc
      CALL DEFMACSET(ABORT,ERROR,3)
      END IF
      END IF
C ---------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'PROC') THEN
      CALL NEXTW2('PROCEDURE-NAME=')
      IF(WD(1:1).EQ.'?') THEN
      CALL DEFVARDUMP(WD,WDLEN,2,3)
      ELSE
      CALL DEFGETSSC(DEFSETSCOPE,DEFCURSCOPE)
      CALL DEFPROCSET(ERROR)
      END IF
C ---------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'CALL') THEN
      IF((BUFIND.EQ.BUFFIL).OR.(.NOT.QADD(NSTRM))) THEN
C ------------------------------------------------------------
C procedure call (CALL) parsing NOT in rotating command buffer
C ------------------------------------------------------------
      CALL NEXTW2('PROCEDURE-NAME=')
      IF(QADD(NSTRM)) THEN
C don't buffer the remaining text on the line (it will be buffered later)
      BUFLYN(BUFSTK)(CURSOR:CURSOR) = ' '
      BUFLEN(BUFSTK)=CURSOR
      END IF
      IF(WD(1:1).EQ.'?') THEN
      DO I=1,DEFCURSCOPE
      CALL DEFCDUMP(DEFVARLIST(1+DEFCURSCOPE-I,3),
     & .TRUE.,1+DEFCURSCOPE-I,3,' ',0)
      END DO
      ELSE
      DEFSETSCOPE = DEFCURSCOPE
      DEFEVLSCOPE = DEFCURSCOPE
      CALL DEFPROCINV(ABORT,ERROR)
      END IF
      ELSE
C ----------------------------------------------------------------
C procedure call (CALL) parsing INSIDE the rotating command buffer
C ----------------------------------------------------------------
C scope creation for procedures
      CALL DEFNEWSCOPE
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      CALL NEXTW2('PROCEDURE-NAME=')
C make sure the next word is a parenthesis
      CALL DEFNEXTW2(5)
C check for problems
      IF((WD(1:1).NE.'(').OR.(WDLEN.NE.1)) THEN
      WRITE (6,'(A)') ' %PROCEDURE-ERR: ( expected but not found.'
      ERROR=.TRUE.
      ELSE
C specify scope/stream into which defines will be stored
      DEFSETSCOPE = DEFCURSCOPE
      DEFEVLSCOPE = DEFCURSCOPE
C parse parameter block
      CALL DEFMACSET(ABORT,ERROR,5)
      END IF
C
      IF(.NOT.ERROR.AND..NOT.ABORT) THEN
C now check for second set of parenthesis
      CALL DEFNEXTW2(4)
C check for problems
      IF((WD(1:1).NE.'(').OR.(WDLEN.NE.1)) THEN
      WRITE (6,'(A)') ' %PROCEDURE-ERR: ( expected but not found.'
      ERROR=.TRUE.
C parse invocation parameter block
      ELSE
      DEFEVLSCOPE = MAX(1,DEFCURSCOPE-1)
      CALL DEFMACSET(ABORT,ERROR,4)
      END IF
      END IF
C
      IF(ERROR) THEN
C procedure invocation errored
      CALL WRNDIE(-5,'PROCEDURE','procedure expansion failed.')
      END IF
      END IF
      END IF
C ---------------------------------------------------------------------
CCC      ELSE IF (WD(1:3).EQ.'GET') THEN  ! GET will be supported soon
CCC      CALL DOGET
C ---------------------------------------------------------------------
      ELSE IF (WD(1:6).EQ.'BUFFER') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      CALL NEXTW2('BUFFER-NAME=')
      IF(WD(1:1).EQ.'?') THEN
      CALL DEFVARDUMP(WD,WDLEN,2,4)
      ELSE
      DEFSETSCOPE = DEFCURSCOPE
      CALL DOBUFFER
      END IF
      END IF
C ---------------------------------------------------------------------
C process-open-command
C --------------------
      ELSE IF (WD(1:4).EQ.'OPEN') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
C
C defaults
      FORM='FORMATTED'
      ACCESS='READ'
C
C parsing
      CALL NEXTFI('FILENAME=',IFILE)
      CALL PUSEND('OPEN>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('OPEN>')
      IF (WD(1:4).EQ.'FORM') THEN
      CALL NEXTA4('FORMat=',STEMP)
      CALL COPYST(FORM,12,I,WD,WDLEN)
      ELSE IF (WD(1:4).EQ.'ACCE') THEN
      CALL NEXTA4('ACCEss=',STEMP)
      CALL COPYST(ACCESS,12,I,WD,WDLEN)
      ELSE
      CALL CHKEND('OPEN>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      CALL ASSFIL(IFILE,UNIT,ACCESS,FORM,ERROR)
      END IF
C
C ---------------------------------------------------------------------
C process-close-command
C ---------------------
      ELSE IF (WD(1:4).EQ.'CLOS') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
C
C defaults
      DISPOS='KEEP'
C
C parsing
      CALL NEXTFI('FILENAME=',IFILE)
      CALL PUSEND('CLOSE>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('CLOSE>')
      IF (WD(1:4).EQ.'DISP') THEN
      CALL NEXTA4('DISPosition=',STEMP)
      CALL COPYST(DISPOS,6,I,WD,WDLEN)
      ELSE
      CALL CHKEND('CLOSE>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      CALL VINQRE('FILE',IFILE,0,0,QOPEN,QFORM,QWRITE,UNIT)
      IF (QOPEN) THEN
C
C close the unit
      CALL VCLOSE(UNIT,DISPOS,ERROR)
C
      ELSE
      CALL WRNDIE(4,PROMPT,'File '//IFILE(1:INDEX(IFILE,'  '))//
     & ' was not open.')
      END IF
      END IF
C
C ---------------------------------------------------------------------
C process-rewind-command
C ----------------------
      ELSE IF (WD(1:4).EQ.'REWI') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
C
C rewinds specified file
      CALL NEXTFI('REWIND-FILE=',IFILE)
      CALL PUSEND('REWIND>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('REWIND>')
      CALL CHKEND('REWIND>',DONE)
      END DO
      DONE=.FALSE.
      CALL VINQRE('FILE',IFILE,0,0,QOPEN,QFORM,QWRITE,UNIT)
      IF (.NOT.QOPEN) THEN
      CALL WRNDIE(4,PROMPT,'File was not open.')
      ELSE
C
C get the full file name and rewind it
      CALL VINQRE('UNIT',WDT,WDTMAX,WDTLEN,QOPEN,QFORM,QWRITE,UNIT)
      REWIND(UNIT=UNIT)
      WRITE(6,'(2A)') ' MISCOM: rewinding file ',WDT(1:WDTLEN)
      END IF
      END IF
C
C ---------------------------------------------------------------------
C process-fileexist-command
C ------------------------
      ELSE IF (WD(1:8).EQ.'FILEEXIS') THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
C
C carries out inquire with exist for specified filename
      CALL NEXTFI('FILE-EXIST-FILE=',IFILE)
      CALL PUSEND('FILE-EXIST>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('FILE-EXIST>')
      CALL CHKEND('FILE-EXIST>',DONE)
      END DO
      DONE=.FALSE.
      QEXIST=.FALSE.
      INQUIRE(FILE=IFILE,EXIST=QEXIST)
      IF (QEXIST) THEN
         WRITE(6,'(A)') ' MISCOM: file exists'
         CALL DECLAR('RESULT', 'LO' ,'TRUE', DCVAL, DPVAL )
      ELSE
         WRITE(6,'(A)') ' MISCOM: file does not exist'
         CALL DECLAR('RESULT', 'LO' ,'FALSE', DCVAL, DPVAL )
      END IF
      END IF
C ---------------------------------------------------------------------
C process-evaluate-command
C ------------------------
      ELSEIF ( WD(1:4) .EQ. 'EVAL' ) THEN
      IF (ENDACT(ENDIND).EQ.'GO  ') CALL EVAL
C
C ---------------------------------------------------------------------
C process-if-command
C ------------------
      ELSE IF (WD(1:4).EQ.'IF  ') THEN
      CALL PUSEND('IF')
      CALL NEXTCD(PROMPT,OK)
      IF (OK) THEN
      CALL PUSACT('GO  ')
      ELSE
      CALL PUSACT('SCAN')
      END IF
      CALL NEXTWD('IF-statement=')
      IF (WD(1:4).NE.'THEN') CALL DSPERR(PROMPT,'THEN expected')
C
C ---------------------------------------------------------------------
C process-else-command
C --------------------
      ELSE IF (WD(1:6).EQ.'ELSE  ') THEN
      IF (ENDKEY(ENDIND)(1:4).NE.'IF  ') THEN
      CALL DSPERR(PROMPT,'no IF ?')
      ELSE IF (ENDACT(ENDIND).EQ.'SCAN') THEN
      CALL PUSACT('GO  ')
      ELSE
      CALL PUSACT('SKIP')
      END IF
C
C --------------------------------------------------------------------
C process-elseif-command
C ----------------------
      ELSE IF (WD(1:6).EQ.'ELSEIF') THEN
CCC
      IF(ENDACT(ENDIND).EQ.'GO  ') THEN
      CALL PUSACT('SKIP')
      END IF
CCC
      CALL NEXTCD(PROMPT,OK)
      CALL NEXTWD('ELSEIF-statement=')
      IF (WD(1:4).NE.'THEN') CALL DSPERR(PROMPT,'THEN expected')
      IF (ENDKEY(ENDIND)(1:4).NE.'IF  ') THEN
      CALL DSPERR(PROMPT,'no IF ?')
      ELSE IF (OK.AND.ENDACT(ENDIND).EQ.'SCAN') THEN
      CALL PUSACT('GO  ')
      ELSE IF (.NOT.OK.AND.ENDACT(ENDIND).EQ.'SCAN') THEN
      CALL PUSACT('SCAN')
      ELSE
      CALL PUSACT('SKIP')
      END IF
C
C ---------------------------------------------------------------------
C process-loop-command
C --------------------
      ELSE IF (WD(1:4).EQ.'LOOP'.AND.QADD(NSTRM)) THEN
      CALL PUSEND('LOOP')
      CALL PUSACT('GO  ')
      CALL NEXTA4('LOOP-LABEL=',ENDLBL(ENDIND))
C
C --------------------------------------------------------------------
C process-while-command
C ---------------------
      ELSE IF (WD(1:4).EQ.'WHIL'.AND.QADD(NSTRM)) THEN
      CALL PUSEND('LOOP')
      CALL NEXTCD(PROMPT,OK)
      IF (OK) THEN
      CALL PUSACT('GO  ')
      ELSE
      CALL PUSACT('SKIP')
      END IF
      CALL NEXTWD('WHILE-clause=')
      IF (WD(1:4).NE.'LOOP') CALL DSPERR(PROMPT,'LOOP expected')
      CALL NEXTA4('LOOP-LABEL=',ENDLBL(ENDIND))
C
C --------------------------------------------------------------------
C process-for-command
C -------------------
      ELSE IF (WD(1:4).EQ.'FOR '.AND.QADD(NSTRM)) THEN
      CALL PUSEND('LOOP')
      CALL PUSACT('GO  ')
      IF (ENDACT(ENDIND).EQ.'GO  ') ENDCNT(ENDIND)=ENDCNT(ENDIND)+1
C
C get name of loop index (symbol) and copy the symbol name into WDD
      CALL NEXTSL('FOR-symbol=')
      IF (WD(1:1).NE.'$') CALL DSPERR('NEXTCD','Symbol expected')
      CALL COPYST(WDD,WDMAX,WDDLEN,WD(2:WDLEN),WDLEN-1)
C for the declaration of the WD type - JSJ 16-DEC-94
      WDTYP=WD(2:WDLEN)//'_TYPE'
      WDTYPL=WDLEN+4
C
C strip off the IN word
      CALL NEXTWD('FOR-clause=')
      IF (WD(1:4).NE.'IN  ') CALL DSPERR(PROMPT,'IN expected')
C
C see what's next
      CALL NEXTSL('FOR-clause=')
      IF (WD(1:1).EQ.'(') THEN
C
C
C enumeration FOR statement
C -------------------------
C
      DO I=1,ENDCNT(ENDIND)
C
C enable the symbol subsitution for the filename only
C fixed by JSJ  21-DEC-94
      QSUBS=.FALSE.
      QEXPRS=.TRUE.
      CALL NEXTW2('FOR-clause')
      QSUBS=.TRUE.
      QEXPRS=.FALSE.
      CALL SAVEWD
      IF (WD(1:1).EQ.'@'.AND..NOT.QQUOT) THEN
      CALL NEXTWD('FOR-clause=')
      ELSE
      CALL NEXTSL('FOR-clause=')
      END IF
C
      END DO
C
      IF (WD(1:1).EQ.')'.AND..NOT.QQUOT) THEN
      CALL PUSACT('SKIP')
      ENDCNT(ENDIND)=0
      ELSE
C do the symbol declaration only if we're active
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
C
C check whether word is a symbol
      IF ( WD(1:1) .EQ. '$' ) THEN
      CALL WDSUB(WD,WDMAX,WDLEN,OK,FVARTP,DPVAL,DCVAL )
      IF (.NOT.OK) THEN
      CALL DSPERR('WDSUB','symbol not found')
      END IF
C
C check whether word is a quoted string.
      ELSE IF (QQUOT) THEN
      FVARTP='ST'
C
C check whether word can be interpreted as a number
      ELSE
C
C get the sign if any
      IF (WD(1:1).EQ.'+') THEN
      ENDCNT(ENDIND)=ENDCNT(ENDIND)+1
      CALL NEXTSL('FOR-clause=')
      SIGNUM=+1
      ELSE IF (WD(1:1).EQ.'-') THEN
      ENDCNT(ENDIND)=ENDCNT(ENDIND)+1
      CALL NEXTSL('FOR-clause=')
      SIGNUM=-1
      ELSE
      SIGNUM=+1
      END IF
C
      CALL CHKNUM(WD,WDLEN,OK)
      IF (OK) THEN
C
C it looks like a number, try to decode it
      FVARTP = 'DP'
      DPVAL = SIGNUM*DECODF(WD,WDLEN,OK)
C here we need a error message
      IF (.NOT.OK) CALL DSPERR(PROMPT,'Error converting number')
      ELSE
C
C we assume that word is an unquoted literal string constant
      FVARTP = 'ST'
      IF ( WRNLEV .GE. 10 ) THEN
      WRITE(6,'(3A)') ' Assuming literal string "',WD(1:WDLEN),'"'
      ENDIF
      ENDIF
      ENDIF
C
C now declare the symbol
C because of symbol scope, these declars need to be done after
C the iteration loop has been parsed - so we need to save copies
C and set a flag
      QDECLARE = .TRUE.
      CALL COPYST(WDDCPY,WDDMAX,WDDCPYLEN,WDD,WDDLEN)
      CALL COPYST(WDCPY,WDMAX,WDCPYLEN,WD,WDLEN)
      CALL COPYST(WDTYPCPY,WDMAX+4,WDTYPCPYL,WDTYP,WDTYPL)
      CALL COPYST(COMCPY,COMMAX,COMCPYLEN,COMLYN,COMLEN)
C
      IF (WRNLEV.GE.5) THEN
C
C echo the assignment
      IF (FVARTP.EQ.'ST') THEN
      WRITE(6,'(5A)')
     & ' FOR LOOP: symbol ',WDD(1:WDDLEN),' set to "',
     & WD(1:WDLEN),'" (string)'
      ELSE IF (FVARTP.EQ.'DP') THEN
      WRITE(6,'(3A,G14.6,A)')
     & ' FOR LOOP: symbol ',WDD(1:WDDLEN),' set to ',
     & DPVAL,' (real)'
      ELSE IF (FVARTP.EQ.'DC') THEN
      WRITE(6,'(3A,G14.6,A,G14.6,A)')
     & ' FOR LOOP: symbol ',WDD(1:WDDLEN),' set to (',
     & DBLE(DCVAL),',',DIMAG(DCVAL),') (complex)'
      END IF
      END IF
C
      END IF
C
C read until we've reached the closing parenthesis
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      CALL NEXTSL('FOR-clause=')
C
      IF (WD(1:1).EQ.')'.AND..NOT.QQUOT) CLOOP=.FALSE.
      END DO
C
      IF (QDECLARE) THEN
      QDECLARE = .FALSE.
C  declare the symbols since we will be now be in the proper scope
C
      CALL LDECLAR(WDDCPY(1:WDDCPYLEN), FVARTP ,
     & WDCPY(1:WDCPYLEN), DCVAL, DPVAL )
C
C declare the type of the symbol WDD - JSJ 16-DEC-94
      CALL LDECLAR(WDTYPCPY(1:WDTYPCPYL), 'ST',
     & FVARTP, DCVAL, DPVAL)
C
C declare the line/record symbol ("$for_line") - JSJ 16-DEC-94
      CALL LDECLAR( 'FOR_LINE', 'ST',
     & COMCPY(1:COMCPYLEN), DCVAL, DPVAL)
      END IF
C
      END IF
      CALL NEXTWD('FOR-clause=')
C
C FOR statement over atom ID's
C ----------------------------
      ELSE IF (WD(1:4).EQ.'ID  ') THEN
      IF (ENDACT(ENDIND).NE.'GO  ') THEN
C
C we're not active --> just parse until we reach LOOP
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      CALL NEXTWD('FOR-clause=')
      IF (WD(1:4).EQ.'LOOP') CLOOP=.FALSE.
      END DO
      ELSE IF (ENDCNT(ENDIND).EQ.1) THEN
C
C We're going through the loop the first time -> read selection.
C Selection will be stored in HEAP
C space to be used by consecutive cycles through the loop.
      ENDPTR(ENDIND)=ALLHP(INTEG4(NATOM))
      ENDNAT(ENDIND)=NATOM
      CALL SELCTA(HEAP(ENDPTR(ENDIND)),ENDSEL(ENDIND),X,Y,Z,.TRUE.)
C
C make an atom index list from the selection
      CALL MAKIND(HEAP(ENDPTR(ENDIND)),NATOM,ENDSEL(ENDIND))
C
C if zero atoms are selected, skip the loop
      IF (ENDSEL(ENDIND).EQ.0) THEN
      CALL PUSACT('SKIP')
      ENDCNT(ENDIND)=0
      CALL FREHP(ENDPTR(ENDIND),INTEG4(ENDNAT(ENDIND)))
      ELSE
C
C declare the symbol for the first selected atom
      CALL FORPUS(HEAP(ENDPTR(ENDIND)),ENDCNT(ENDIND),WDD(1:WDDLEN))
      END IF
C
C read next word (should be LOOP, test is below)
      CALL NEXTWD('FOR-clause=')
C
      ELSE IF (ENDCNT(ENDIND).GT.ENDSEL(ENDIND)) THEN
C
C we've gone through all selected atoms.  Terminate the loop and
C release the HEAP space that was used for the atom selection.
      CALL PUSACT('SKIP')
      ENDCNT(ENDIND)=0
      CALL FREHP(ENDPTR(ENDIND),INTEG4(ENDNAT(ENDIND)))
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      CALL NEXTWD('FOR-clause=')
      IF (WD(1:4).EQ.'LOOP') CLOOP=.FALSE.
      END DO
C
      ELSE
C
C we're going through the loop not the first time.  Use the
C atom selection stored on the HEAP.
C
C check that the number of atoms hasn't changed between loop
C iterations
      IF (NATOM.NE.ENDNAT(ENDIND)) THEN
      CALL WRNDIE(-1,'MISCOM',
     & 'Number of atoms has changed within FOR IN ID loop. Not allowed')
      END IF
C
C declare the symbol for the currently selected atom
      CALL FORPUS(HEAP(ENDPTR(ENDIND)),ENDCNT(ENDIND),WDD(1:WDDLEN))
C
C skip over the selection statement (we're careful here
C and check that the parenthesis aren't part of quoted
C strings).
      I=0
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      CALL NEXTWD('FOR-clause=')
      IF (WD(1:1).EQ.'('.AND..NOT.QQUOT) THEN
      I=I+1
      ELSE IF (WD(1:1).EQ.')'.AND..NOT.QQUOT) THEN
      I=I-1
      END IF
      IF (I.EQ.0) CLOOP=.FALSE.
      END DO
      CALL NEXTWD('FOR-clause=')
      END IF
C
      END IF
C
      IF (WD(1:4).NE.'LOOP') CALL DSPERR(PROMPT,'LOOP expected')
      CALL NEXTA4('LOOP-LABEL=',ENDLBL(ENDIND))
C
C --------------------------------------------------------------------
C process-exit <loop-label> -command
C ----------------------------------
      ELSE IF (WD(1:4).EQ.'EXIT') THEN
      CALL NEXTA4('LOOP-LABEL=',STEMP)
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      I=ENDIND
      DO WHILE (I.GT.1.AND.(ENDKEY(I).NE.'LOOP'.OR.ENDLBL(I).NE.STEMP))
      ENDACT(I)='SKIP'
      ENDCNT(I)=0
      I=I-1
      END DO
      IF (ENDKEY(I).EQ.'LOOP'.AND.ENDLBL(I).EQ.STEMP) THEN
      ENDACT(I)='SKIP'
      ENDCNT(I)=0
      END IF
      END IF
C
C --------------------------------------------------------------------
C skip all other commands if ENDACT(ENDIND).NE.'GO'
C -------------------------------------------
      ELSE IF (ENDACT(ENDIND).NE.'GO  '.AND..NOT.EOF) THEN
      OK=.TRUE.
      CALL CHKEND(PROMPT,OK)
C
C --------------------------------------------------------------------
C otherwise return word with LUSED=.FALSE.
C ----------------------------------------
      ELSE
      LUSED=.FALSE.
      END IF
C
      RETURN
      END
C=============================================================
      SUBROUTINE FORPUS(IFORID,N,VARNAM)
C
C Routine declares symbol (name VARNAM) to be the
C index of the Nth selected atom.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'timer.inc'
      INTEGER IFORID(*), N
      CHARACTER*(*) VARNAM
C local
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
C begin
      DPVAL=IFORID(N)
      CALL LDECLAR(VARNAM,'DP',' ',DCVAL,DPVAL)
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(3A,G14.6,A)')
     & ' FOR ID LOOP: symbol ',VARNAM,' set to ',
     & DPVAL,' (real)'
      END IF
      RETURN
      END
C=============================================================
      SUBROUTINE NEXTCD(PROMPT,OK)
C
C Parses an conditional expression:
C ( <symbol> <condition>  <word> )
C <condition> can be:  > < = <= >= #       gt lt le ge eq ne
C
C Author: Axel Brunger
C ====================
C I/O
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'funct.inc'
      CHARACTER*(*) PROMPT
      LOGICAL OK
C local
      INTEGER N, SIGNUM, SIGNUM2
      LOGICAL FOUND, TMATCH, ERR
      CHARACTER*2 TYPE1, TYPE2, RELAT
      DOUBLE PRECISION DPVAL1, DPVAL2
      DOUBLE COMPLEX DCVAL1, DCVAL2
      LOGICAL CLOOP
C begin
C
C initialize OK, ERR
      OK=.FALSE.
      ERR=.FALSE.
      IF (ENDACT(ENDIND).EQ.'GO  '.OR.ENDACT(ENDIND).EQ.'SCAN') THEN
C
C pick up the opening parenthesis
      CALL NEXTWD('CONDITION=')
      IF (WD(1:1).NE.'(') THEN
      CALL DSPERR('CONDITION','"(" expected')
      ERR=.TRUE.
      END IF
C
C get the l.h.s.
      CALL NEXTSL('CONDITION=')
      SIGNUM=+1
      TYPE1='??'
C
C
C copy into WDD
      CALL COPYST(WDD,WDMAX,WDDLEN,WD,WDLEN)
C
      IF (.NOT.QQUOT) THEN
      IF (WDD(1:1).EQ.'+') THEN
      CALL NEXTSL('CONDITION=')
      CALL COPYST(WDD,WDMAX,WDDLEN,WD,WDLEN)
      SIGNUM=+1
      ELSE IF (WDD(1:1).EQ.'-') THEN
      CALL NEXTSL('CONDITION=')
      CALL COPYST(WDD,WDMAX,WDDLEN,WD,WDLEN)
      SIGNUM=-1
      END IF
C
C substitute lhs symbol if possible
      IF (WDD(1:1).EQ.'$') THEN
      CALL WDSUB(WDD,WDMAX,WDDLEN,FOUND,TYPE1,DPVAL1,DCVAL1)
      IF ( .NOT. FOUND ) THEN
      CALL DSPERR('NEXTCD','Symbol not found')
      ERR=.TRUE.
      END IF
      END IF
C
      IF (TYPE1.EQ.'??') THEN
C is it a double precision number?
      CALL CHKNUM(WDD,WDDLEN,TMATCH)
      IF (TMATCH) THEN
      DPVAL1 = SIGNUM*DECODF( WDD, WDDLEN, TMATCH )
      END IF
      IF (TMATCH) TYPE1='DP'
      END IF
      END IF
C
C is it a logical?
      IF (WDD(1:WDDLEN).EQ.'FALSE'.OR.WDD(1:WDDLEN).EQ.'TRUE') THEN
      TYPE1='LO'
      END IF
C
C if must be a string!
      IF (TYPE1.EQ.'??') TYPE1='ST'
C
C
C now get the relational operator
      CALL NEXTSL('CONDITION=')
      RELAT(1:2)='  '
      N=MIN(2,WDLEN)
      RELAT(1:N)=WD(1:N)
C
C check if we have a composite relat. operator, such as ">=" or "<=" or "=="
      IF (RELAT.EQ.'< '.OR.RELAT.EQ.'> '.OR.RELAT.EQ.'= ') THEN
      CALL NEXTSL('CONDITION=')
      IF (WD(1:1).EQ.'=') THEN
      RELAT(2:2) = '='
      ELSE
      CALL SAVEWD
      END IF
      END IF
C
C convert substitution operators
      IF (RELAT.EQ.'EQ') RELAT = '= '
      IF (RELAT.EQ.'NE') RELAT = '# '
      IF (RELAT.EQ.'GT') RELAT = '> '
      IF (RELAT.EQ.'GE') RELAT = '>='
      IF (RELAT.EQ.'LT') RELAT = '< '
      IF (RELAT.EQ.'LE') RELAT = '<='
C
C check that we have a legal relational operator
      IF (RELAT.NE.'= '.AND.
     & RELAT.NE.'# '.AND.
     & RELAT.NE.'> '.AND.
     & RELAT.NE.'=='.AND.
     & RELAT.NE.'>='.AND.
     & RELAT.NE.'< '.AND.
     & RELAT.NE.'<=' )  THEN
      CALL DSPERR('NEXTCD',
     & 'Illegal relational operator, only = # < > <= >= == allowed')
      ERR=.TRUE.
      END IF
C
C get the r.h.s. of the relation
      CALL NEXTSL('COMPARISON=')
      SIGNUM2=+1
      TYPE2='??'
      IF (.NOT.QQUOT) THEN
      IF (WD(1:1).EQ.'+') THEN
      CALL NEXTSL('CONDITION=')
      SIGNUM2=+1
      ELSE IF (WD(1:1).EQ.'-') THEN
      CALL NEXTSL('CONDITION=')
      SIGNUM2=-1
      END IF
C
C
C substitute lhs symbol if possible
      IF (WD(1:1).EQ.'$') THEN
      CALL WDSUB(WD,WDMAX,WDLEN,FOUND,TYPE2,DPVAL2,DCVAL2)
      IF ( .NOT. FOUND ) THEN
      CALL DSPERR('NEXTCD','Symbol not found')
      ERR=.TRUE.
      END IF
      END IF
C
      IF (TYPE2.EQ.'??') THEN
C is it a double precision number?
      CALL CHKNUM(WD,WDLEN,TMATCH)
      IF (TMATCH) THEN
      DPVAL2 = SIGNUM2*DECODF( WD, WDLEN, TMATCH )
      END IF
      IF (TMATCH) TYPE2='DP'
      END IF
      END IF
C
C
C is it a logical?
      IF ((WDD(1:WDDLEN).EQ.'FALSE'.OR.WDD(1:WDDLEN).EQ.'TRUE').AND.
     &     .NOT.QQUOT) THEN
      TYPE2='LO'
      END IF
C
C it must be a string!
      IF (TYPE2.EQ.'??') TYPE2='ST'
C
      IF (TYPE1.NE.TYPE2) THEN
C
      IF (TYPE1.EQ.'LO'.OR.TYPE2.EQ.'LO') THEN
C
C always mismatches involving logicals
      CALL DSPERR('NEXTCD','Mismatched data types')
      ERR=.TRUE.
      ELSEIF (RELAT.NE.'# '.AND.RELAT.NE.'= ') THEN
C
C otherwise we allow mismatches for equal and "not equal" operations
      CALL DSPERR('NEXTCD','Mismatched data types')
      ERR=.TRUE.
      END IF
      END IF
C
      IF (RELAT.NE.'= '.AND.
     & RELAT.NE.'# '.AND.TYPE1.EQ.'DC') THEN
      CALL DSPERR('NEXTCD',
     & 'Illegal relat. operator for complex type, only = # allowed')
      ERR=.TRUE.
      END IF
C
      IF (RELAT.NE.'= '.AND.
     & RELAT.NE.'# '.AND.TYPE1.EQ.'LO') THEN
      CALL DSPERR('NEXTCD',
     & 'Illegal relat. operator for logical type, only = # allowed')
      ERR=.TRUE.
      END IF
C
      IF (.NOT.ERR) THEN
C======================================
      IF (TYPE1.NE.TYPE2.AND.RELAT.EQ.'= ') THEN
      OK=.FALSE.
      ELSEIF (TYPE1.NE.TYPE2.AND.RELAT.EQ.'# ') THEN
      OK=.TRUE.
C string comparisons
      ELSEIF ( TYPE1 .EQ. 'ST' ) THEN
      IF      (RELAT.EQ.'= ') THEN
      CALL EQSTWC(WDD,WDDLEN,WD,WDLEN,1,1,OK)
      ELSEIF (RELAT.EQ.'==') THEN
      OK = WDD(1:WDDLEN).EQ.WD(1:WDLEN)
      ELSEIF (RELAT.EQ.'# ') THEN
      CALL EQSTWC(WDD,WDDLEN,WD,WDLEN,1,1,OK)
      OK=.NOT.OK
      ELSEIF (RELAT.EQ.'> ') THEN
      OK = WDD(1:WDDLEN).GT.WD(1:WDLEN)
      ELSEIF (RELAT.EQ.'< ') THEN
      OK = WDD(1:WDDLEN).LT.WD(1:WDLEN)
      ELSEIF (RELAT.EQ.'>=') THEN
      OK = WDD(1:WDDLEN).GE.WD(1:WDLEN)
      ELSEIF (RELAT.EQ.'<=') THEN
      OK = WDD(1:WDDLEN).LE.WD(1:WDLEN)
      END IF
C======================================
C double precision comparisons
      ELSEIF ( TYPE1 .EQ. 'DP' ) THEN
      IF      (RELAT.EQ.'> ') THEN
      OK = DPVAL1 .GT. DPVAL2
      ELSEIF (RELAT.EQ.'>=') THEN
      OK = DPVAL1 .GE. DPVAL2
      ELSEIF (RELAT.EQ.'< ') THEN
      OK = DPVAL1 .LT. DPVAL2
      ELSEIF (RELAT.EQ.'<=') THEN
      OK = DPVAL1 .LE. DPVAL2
      ELSEIF (RELAT.EQ.'= '.OR.RELAT.EQ.'==') THEN
      OK = ABS(DPVAL1-DPVAL2).LT.RSMALL
      ELSEIF (RELAT.EQ.'# ') THEN
      OK = ABS(DPVAL1-DPVAL2).GT.RSMALL
      END IF
C======================================
C double complex comparisons
      ELSEIF ( TYPE1 .EQ. 'DC' ) THEN
      IF (RELAT.EQ.'= '.OR.RELAT.EQ.'==') THEN
      OK =      ABS(DBLE(DCVAL1)-DBLE(DCVAL2)).LT.RSMALL
     & .AND.ABS(DIMAG(DCVAL1)-DIMAG(DCVAL2)).LT.RSMALL
      ELSEIF (RELAT.EQ.'# ') THEN
      OK =      ABS(DBLE(DCVAL1)-DBLE(DCVAL2)).GT.RSMALL
     & .AND.ABS(DIMAG(DCVAL1)-DIMAG(DCVAL2)).GT.RSMALL
      END IF
C======================================
C logical comparisons
      ELSEIF ( TYPE1 .EQ. 'LO' ) THEN
      IF      (RELAT.EQ.'= '.OR.RELAT.EQ.'==') THEN
      CALL EQSTWC(WDD,WDDLEN,WD,WDLEN,1,1,OK)
      ELSEIF (RELAT.EQ.'# ') THEN
      CALL EQSTWC(WDD,WDDLEN,WD,WDLEN,1,1,OK)
      OK=.NOT.OK
      END IF
C======================================
      END IF
      END IF
C
C pick up the closing parenthesis
      CALL NEXTWD('CONDITION=')
      IF (WD(1:1).NE.')') CALL DSPERR('CONDITION','")" expected')
C
      IF (WRNLEV.GE.5) THEN
      IF (OK) THEN
      WRITE(6,'(A)') ' NEXTCD: condition evaluated as true'
      ELSE
      WRITE(6,'(A)') ' NEXTCD: condition evaluated as false'
      END IF
      END IF
      ELSE
C
C read until we've reached the closing parenthesis
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      CALL NEXTSL('CONDITION=')
      IF (WD(1:1).EQ.')'.AND..NOT.QQUOT) CLOOP=.FALSE.
      END DO
C
      END IF
C
      RETURN
      END
C===============================================================
      SUBROUTINE DECLAR( NAME1, TYPE, STINIT, DCINIT, DPINIT )
C
C Safe wrapper around DECSYM for interal symbol definitions:
C defines variables into all active scopes!
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) NAME1
      CHARACTER*2 TYPE
      CHARACTER*(*) STINIT
      DOUBLE PRECISION DPINIT
      DOUBLE COMPLEX DCINIT
C local
      INTEGER NAMELEN, SCOPE, MATCHLEN, RECNO
      LOGICAL FOUND, QPART
      CHARACTER*(VARMAX) NAME
C begin
      NAMELEN = LEN(NAME1)
      IF(NAMELEN.GT.VARMAX) NAMELEN = VARMAX
      FOUND = .FALSE.
      DO I=1,DEFCURSCOPE
      NAME = NAME1(1:NAMELEN)
      SCOPE = I
      MATCHLEN = NAMELEN
C look for symbol in this scope
      IF(I.GT.1) THEN
      CALL DEFFINDREC(NAME,MATCHLEN,SCOPE,1,.FALSE.,
     & FOUND,RECNO,QPART,.TRUE.)
      END IF
      IF(FOUND.OR.(I.EQ.1)) THEN
C if found, then replace
      SCOPE = I
      CALL DECSYM(NAME,NAMELEN, TYPE, STINIT, DCINIT,
     &            DPINIT, SCOPE)
      END IF
      END DO
C
      RETURN
      END
C===============================================================
      SUBROUTINE LDECLAR( NAME1, TYPE, STINIT, DCINIT, DPINIT )
C
C Safe wrapper around DECSYM for external symbol definitions:
C defines variables into the current scope (DEFCURSCOPE)
C
C NOTE: NAME1 will be modified with substituted symbol name.
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) NAME1
      CHARACTER*2 TYPE
      CHARACTER*(*) STINIT
      DOUBLE PRECISION DPINIT
      DOUBLE COMPLEX DCINIT
C local
      INTEGER NAMELEN, SCOPE
C begin
      NAMELEN = LEN(NAME1)
      IF(NAMELEN.GT.VARMAX) NAMELEN = VARMAX
      CALL DEFGETSSC(SCOPE,DEFCURSCOPE)
      CALL DECSYM(NAME1,NAMELEN, TYPE, STINIT, DCINIT,
     &            DPINIT, SCOPE)
C
      RETURN
      END
C===============================================================
      SUBROUTINE DECSYM( NAME1, NAME1LEN, TYPE, STINIT, DCINIT,
     & DPINIT, INSCOPE )
C
C Add a new variable to the symbol tables, if it doesn't
C already exist, reinitialize it if it does.
C
C I/O descriptions
C
C  name ->  the variable identifier, max of VARMAX characters, VARMAX
C           defined in symbol.inc
C  type ->  variable type, 2 characters, 'DP' for double precision,
C           'DC' for double complex, 'ST' for string, 'LO' for logical.
C  stinit -> if type is 'ST' uses this as initialization value.
C  stinit -> if type is 'LO' uses this as initialization value.
C  dcinit -> if type is 'DC' uses this as initialization value.
C  dpinit -> if type is 'DP' uses this as initialization value.
C
C Authors: Axel T. Brunger and Mark McCallum
C ==========================================
C Modification: Jian-Sheng Jiang and Warren L. DeLano
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) NAME1
      INTEGER NAME1LEN, INSCOPE
      CHARACTER*2 TYPE
      CHARACTER*(*) STINIT
      DOUBLE PRECISION DPINIT
      DOUBLE COMPLEX DCINIT
C local
      CHARACTER*(VARMAX) NAME, NAME3
      INTEGER ATIND, SCOPE
      INTEGER NAMLEN, INILEN, MATCHLEN, NAM3LEN, OFFSET
      LOGICAL EXISTS, QPART, VARFLG(NUMVARFLG)
C begin
      EXISTS = .FALSE.
C
C get length of the variable identifier string
      NAMLEN = NAME1LEN
C
C check if declaration used double dollar symbols -- strip the second one
      IF (NAMLEN.GT.0) THEN
      IF (NAME1(1:1).EQ.'$') THEN
      NAME(1:NAMLEN)=NAME1(1:NAMLEN)
      NAMLEN=NAMLEN-1
      NAME1(1:NAMLEN)=NAME(2:NAMLEN+1)
      END IF
      END IF
      CALL TRIMM(NAME1,NAMLEN)
C
      IF (NAMLEN.GT.VARMAX) THEN
      CALL WRNDIE(-5,'DECLAR',
     & 'parameter name too long (VARMAX exceeded)')
      ELSE IF (NAMLEN.LE.0) THEN
      CALL WRNDIE(-5,'DECLAR',
     & 'null parameter name not allowed')
      ELSE
C
C make a local copy of the name
      NAME(1:NAMLEN)=NAME1(1:NAMLEN)
C save a copy for comparison
      NAME3=NAME
      NAM3LEN=NAMLEN
C remove any and all directives
      CALL DEFGETFLG(NAME,NAMLEN,.TRUE.,OFFSET,VARFLG)
C substitute indices
      CALL SYMIDX(NAME,VARMAX,NAMLEN,.FALSE.)
C
C see if parameter already exists in top scope
      MATCHLEN = NAMLEN
      SCOPE = INSCOPE
      CALL DEFFINDREC(NAME,MATCHLEN,SCOPE,1,.FALSE.,
     & EXISTS,ATIND,QPART,.TRUE.)
      IF(QPART) THEN
      IF(EXISTS) THEN
      CALL DSPERR('DECLAR',
     & 'symbol already defined as compound.')
      ELSE
      CALL DSPERR('DECLAR',
     & 'symbol already defined as atomic.')
      END IF
      ATIND = 0
      END IF
C allocate a new parameter record if necessary
      IF(.NOT.EXISTS) CALL DEFNEWREC(ATIND)
      IF(ATIND.EQ.0) THEN
      CALL WRNDIE(-5,'DECLAR',
     & 'can not set value of symbol')
      ELSE
C we now have a parameter
      IF(.NOT.EXISTS) THEN
C set the name and insert if necessary
      DEFNAMTXT(ATIND)(1:NAMLEN)= NAME(1:NAMLEN)
      DEFNAMLEN(ATIND)=NAMLEN
      CALL DEFINSREC(ATIND,INSCOPE,1)
      END IF
C return updated name (for the purposes of echoing
C the substituted symbol name in evaluate)
      IF (NAMLEN.GT.0) THEN
      IF (NAM3LEN.NE.NAMLEN.OR.
     &    NAME3(1:NAMLEN).NE.NAME(1:NAMLEN)) THEN
      CALL COPYST(NAME1,VARMAX,MATCHLEN,NAME(1:NAMLEN),NAMLEN)
      END IF
      END IF
C
C assign the type of the variable
      CMDTYP(ATIND) = TYPE
C
C assign the value of the variable
      IF ( TYPE .EQ. 'DP' ) THEN
      DPVALU(ATIND) = DPINIT
C============================================================
      ELSEIF ( TYPE .EQ. 'DC' ) THEN
      DCVALU(ATIND) = DCINIT
C============================================================
      ELSEIF ( TYPE .EQ. 'ST' .OR. TYPE .EQ. 'LO' ) THEN
      INILEN=LEN(STINIT)
      CALL COPYST(DEFPARTXT(ATIND),VARMAX,
     & DEFPARLEN(ATIND),STINIT,INILEN)
C============================================================
      ELSE
      CALL DSPERR('DECLAR','Corrupt variable tables')
      END IF
      END IF
      END IF
C
      RETURN
      END
C============================================================
      SUBROUTINE SUBTILDE(WD,WDLEN,WDMAX)
C
C I/O
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      CHARACTER*(*) WD
      INTEGER WDLEN, WDMAX
C local
      INTEGER MXPATH, LL, TLEN
      CHARACTER*(WORD_SIZE) TEMP
      LOGICAL GOTIT, FINISH
C begin
      GOTIT=.FALSE.
      FINISH=.FALSE.
      LL=0
      DO WHILE (LL.LT.(WDLEN-1).AND..NOT.FINISH)
      LL=LL+1
      IF (WD(LL:LL+1).EQ.'~/') THEN
         GOTIT=.TRUE.
         FINISH=.TRUE.
         LL=LL+1
      ELSE IF (.NOT.GOTIT.AND.WD(LL:LL).NE.' ') THEN
         FINISH=.TRUE.
      END IF
      END DO
C
      IF (GOTIT) THEN
         TEMP(1:5)='HOME:'
         TLEN=5
         CALL ADDST(TEMP,WORD_SIZE,TLEN,WD(LL+1:WDLEN),
     &              WDLEN-LL)
         CALL COPYST(WD,WDMAX,WDLEN,TEMP,TLEN)
      END IF
C
      RETURN
      END
C
C============================================================
C
