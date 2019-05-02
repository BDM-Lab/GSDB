C==================================================================
      SUBROUTINE DEFGETSSC(SETSCOPE,FROMSCOPE)
C
C gets the number of the scope into which parameters and
C symbols should be defined taking into account any inline
C streams or modules, starting in FROMSCOPE
C ( if FROMSCOPE < 1 then SETSCOPE defaults to 1 )
C
C Authors: Warren L. DeLano
C =========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER SETSCOPE, FROMSCOPE
C local
      IF(FROMSCOPE.LT.1) THEN
      SETSCOPE = 1
      ELSE
      SETSCOPE = FROMSCOPE
      DO WHILE((DEFINLCNT(SETSCOPE).EQ.2).AND.(SETSCOPE.GT.1))
      SETSCOPE = SETSCOPE - 1
      END DO
      END IF
      RETURN
      END
C==================================================================
      SUBROUTINE DEFCPYLST(FIRST,NEWLIST)
C
C copies a list of records (DEFMULTREC lists only)
C
C Authors: Warren L. DeLano
C =========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER FIRST, NEWLIST
C local
      INTEGER RECNO, NEWREC, LAST
C begin
      RECNO = FIRST
      LAST = 0
      NEWLIST = 0
      DO WHILE(RECNO.GT.0)
      CALL DEFNEWREC(NEWREC)
      IF(NEWREC.GT.0) THEN
C
      IF(LAST.GT.0) THEN
      DEFMULTREC(LAST) = NEWREC
      ELSE
      NEWLIST = NEWREC
      END IF
      LAST = NEWREC
C
      DEFPARLEN(NEWREC) = DEFPARLEN(RECNO)
      DEFPARTXT(NEWREC) = DEFPARTXT(RECNO)(1:DEFPARLEN(RECNO))
      DEFNAMLEN(NEWREC) = DEFNAMLEN(RECNO)
      DEFNAMTXT(NEWREC) = DEFNAMTXT(RECNO)(1:DEFNAMLEN(RECNO))
      CMDTYP(NEWREC) = CMDTYP(RECNO)
      DCVALU(NEWREC) = DCVALU(RECNO)
      DPVALU(NEWREC) = DPVALU(RECNO)
      BODEST(NEWREC) = BODEST(RECNO)
      BOUNIT(NEWREC) = BOUNIT(RECNO)
C
      END IF
      RECNO = DEFMULTREC(RECNO)
      END DO
      RETURN
      END
C==================================================================
      SUBROUTINE DEFCHKBLK(FIRST,QBLANK)
C
C returns QBLANK = .TRUE. if the records in the list
C contain nothing but spaces
C
C Authors: Warren L. DeLano
C =========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER FIRST
      LOGICAL QBLANK
C local
      INTEGER RECNO
      INTEGER I
C begin
      QBLANK = .TRUE.
      RECNO = FIRST
      DO WHILE((RECNO.GT.0).AND.(QBLANK))
      I = 1
      DO WHILE((QBLANK).AND.(I.LE.DEFPARLEN(RECNO)))
      IF ((DEFPARTXT(RECNO)(I:I).NE.' ').AND.
     &    (DEFPARTXT(RECNO)(I:I).NE.'"')) QBLANK = .FALSE.
      I = I + 1
      END DO
      RECNO = DEFMULTREC(RECNO)
      END DO
      RETURN
      END
C=================================================================
      SUBROUTINE DEFVARDUMP(LWD,LWDLEN,OFFSET,VARTYPE)
C
C write out symbol(s), parameter(s), or procedure(s)
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(WDMAX) LWD
      INTEGER LWDLEN,OFFSET,VARTYPE
C local
      INTEGER I
      CHARACTER*(VARMAX) NAME
      INTEGER NAMLEN, RECNO, MATCHLEN, SCOPE, TMPNEXT
      LOGICAL QPART,FOUND,EVRFND
C begin
      IF(LWDLEN.LT.OFFSET) THEN
C print out the whole list
      IF(VARTYPE.EQ.1) THEN
      WRITE(6,'(A)') ' Current Symbol Table'
      WRITE(6,'(A)')
     & ' NOTE: The scope number (#) of $<name> is shown as $_#_<name> '
      ELSE IF(VARTYPE.EQ.2) THEN
      WRITE (6,'(A)') ' Current Define Parameters'
      WRITE (6,'(A)')
     & ' NOTE: The scope number (#) of &<name> is shown as &_#_<name> '
      ELSE IF(VARTYPE.EQ.3) THEN
      WRITE(6,'(A)') ' Current Procedures'
      WRITE(6,'(A)') ' NOTE: The scope number (#) of the' //
     & ' procedure is shown as _#_<name> '
      END IF
      DO I=1,DEFCURSCOPE
      CALL DEFCDUMP(DEFVARLIST(1+DEFCURSCOPE-I,VARTYPE),
     & .TRUE.,1+DEFCURSCOPE-I,VARTYPE,' ',0)
      END DO
      ELSE
C print out just a single variable (can be a compound variable)
      NAMLEN = LWDLEN - OFFSET + 1
      NAME = LWD(OFFSET:LWDLEN)
      CALL SYMIDX(NAME,VARMAX,NAMLEN,.TRUE.)
      MATCHLEN = NAMLEN
      SCOPE = DEFCURSCOPE
      EVRFND = .FALSE.
      DO WHILE (SCOPE.GT.0)
      CALL DEFFINDREC(NAME,MATCHLEN,SCOPE,VARTYPE,.FALSE.,FOUND,
     & RECNO,QPART,.FALSE.)
      IF(FOUND.OR.QPART) THEN
      EVRFND = .TRUE.
C temporarily unhook this variable from the rest of the list
      TMPNEXT = DEFNEXTREC(RECNO)
      DEFNEXTREC(RECNO) = 0
C print out this one variable
      CALL DEFCDUMP(RECNO,.TRUE.,SCOPE,VARTYPE,
     &  NAME(1:MATCHLEN),MATCHLEN)
      DEFNEXTREC(RECNO) = TMPNEXT
      END IF
      SCOPE = SCOPE - 1
      END DO
      IF(.NOT.EVRFND) THEN
      CALL DSPERR('DEFVARDUMP','variable not found')
      END IF
      END IF
C
      RETURN
      END
C=================================================================
      SUBROUTINE BFODUMP(BUFLIST,BFDEST,FILUNIT)
C
C Write the contents of the output buffer to the
C specified destination
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'ctitla.inc'
      INTEGER BFDEST,FILUNIT,BUFLIST
C local
      CHARACTER*(MXBFLN) LINE
      LOGICAL QDONE
      INTEGER LINELEN,CURREC
      LOGICAL QNLINE
      INTEGER BFUNIT
C begin
C initialize
C 1 for output, 2 for display, 3 for remarks, 4 for file
      IF(BFDEST.EQ.1) THEN
C writing to standard output
      BFUNIT = 6
      ELSE IF(BFDEST.EQ.2) THEN
C writing to the display file
      BFUNIT = DUNIT
      ELSE IF(BFDEST.EQ.3) THEN
C writing to the remarks buffer
      IF (TITMODE.EQ.0) NTITLE=0
      ELSE IF(BFDEST.EQ.4) THEN
C writing to a named file (unit # in FILUNIT)
      BFUNIT = FILUNIT
      END IF
C local copy needed
      CURREC = BUFLIST
C skip over the first record - that just contains the name
      IF(CURREC.GT.0) THEN
      CURREC = DEFMULTREC(CURREC)
      END IF
C initialize
      LINELEN = 0
      QDONE = .FALSE.
      DO WHILE(.NOT.QDONE)
      IF((CURREC.EQ.0).AND.(LINELEN.EQ.0)) THEN
C at end of the list and nothing more to write out
      QDONE = .TRUE.
      END IF
C
      IF(.NOT.QDONE) THEN
C new line flag gets set to true if we need to start a new line
      QNLINE = .FALSE.
      IF(CURREC.EQ.0) THEN
C ...at end of list, but still have stuff to write
      QNLINE = .TRUE.
      ELSE IF(DEFNAMLEN(CURREC).GT.0) THEN
C ...in middle of list, but DEFNAMLEN nonzero, signalling a new line
      QNLINE = .TRUE.
      END IF
C
      IF(QNLINE) THEN
C starting a new line, so write out what has been stored up til now
      IF(LINELEN.EQ.0) THEN
C handle case where line is blank - write one space only
      LINELEN = 1
      LINE = ' '
      END IF
      IF(BFDEST.NE.3) THEN
C nomal output to a unit
      WRITE(BFUNIT,'(A)',ERR=8888) LINE(1:LINELEN)
      GOTO 9999
8888  CALL WRNDIE(-1,'BUFFER','error writing on bufferunit.')
9999  CONTINUE
C
      ELSE
C special output to the remarks buffer
      NTITLE=NTITLE+1
C
C leave space (5 lines) for REMARKs to be added if a file is written out
      IF (NTITLE.LE.(MXTITL-5)) THEN
C
      IF(LINELEN.GT.(TITMAX - 9)) THEN
      CALL WRNDIE(-1,'BUFFER',
     & 'line too long -- exceeded TITMAX')
      LINELEN = TITMAX - 9
      END IF
      TITLE(NTITLE) = ' REMARKS '//LINE(1:LINELEN)
      ELSE
      CALL WRNDIE(-1,'BUFFER',
     & 'exceeded MXTITL (COMAND) parameter --> recompile program')
      NTITLE = NTITLE - 1
      END IF
      END IF
      LINELEN = 0
      ELSE IF(CURREC.GT.0) THEN
C normal append...continue building up the current line
      IF(DEFPARLEN(CURREC).GT.0) THEN
      LINE(LINELEN+1:LINELEN+DEFPARLEN(CURREC)) =
     & DEFPARTXT(CURREC)(1:DEFPARLEN(CURREC))
      LINELEN = LINELEN + DEFPARLEN(CURREC)
      END IF
      END IF
      IF(CURREC.GT.0) THEN
C if we haven't yet reached the end of the list, then go on to the
C next record
      CURREC = DEFMULTREC(CURREC)
      END IF
      END IF
      END DO
C
      RETURN
      END
C
C=================================================================
      SUBROUTINE DOBUFFER
C
C handle BUFFER command
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(VARMAX) BUFNAM
      INTEGER BUFNAMLEN,OUTBUFLIST
      INTEGER LASTREC,RECNO, SCOPE, FRMLIST
      INTEGER BFDEST,BFUNIT,MATCHLEN
      LOGICAL QPART, FOUND
C begin
C initialize
      BFDEST = 1
      BFUNIT = 6
C 1 for output, 2 for display, 3 for remarks, 4 for file
C
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
C
C get buffer name
      BUFNAMLEN = WDLEN
      IF(BUFNAMLEN.GT.VARMAX) BUFNAMLEN=VARMAX
      BUFNAM = WD(1:BUFNAMLEN)
C
C currently buffer names share one scope
      SCOPE = 1
      MATCHLEN = BUFNAMLEN
      CALL DEFFINDREC(BUFNAM,MATCHLEN,SCOPE,4,.FALSE.,
     & FOUND,OUTBUFLIST,QPART,.TRUE.)
      IF((.NOT.FOUND).OR.(QPART)) THEN
C create new buffer if the name is not found
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      DEFNAMTXT(RECNO) = BUFNAM(1:BUFNAMLEN)
      DEFNAMLEN(RECNO) = BUFNAMLEN
C currently buffer names share one scope
      SCOPE = 1
      CALL DEFINSREC(RECNO,SCOPE,4)
      OUTBUFLIST = RECNO
      LASTREC = RECNO
      END IF
      ELSE
C remember what we are writing to
      BFDEST = BODEST(OUTBUFLIST)
      BFUNIT = BOUNIT(OUTBUFLIST)
      END IF
C
C find the end of the list
      LASTREC = OUTBUFLIST
      IF(LASTREC.GT.0) THEN
      DO WHILE(DEFMULTREC(LASTREC).GT.0)
      LASTREC = DEFMULTREC(LASTREC)
      END DO
      END IF
C
      IF(LASTREC.GT.0) THEN
C parsing
      CALL PUSEND('BUFFER>')
      DONE = .FALSE.
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('BUFFER>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-buffer')
C ----------------------------------------
      ELSE IF (WD(1:4).EQ.'RESE') THEN
C destroy contents of the buffer
      IF(OUTBUFLIST.GT.0) THEN
C delete all but the first record in the list
      IF(DEFMULTREC(OUTBUFLIST).GT.0) THEN
      CALL DEFFREEREC(DEFMULTREC(OUTBUFLIST),.TRUE.)
      DEFMULTREC(OUTBUFLIST) = 0
      LASTREC = OUTBUFLIST
      END IF
C
      END IF
C ----------------------------------------
      ELSE IF (WD(1:2).EQ.'TO') THEN
C
      CALL NEXTWD('TO=')
      IF (WD(1:1).EQ.'=') CALL NEXTWD('TO=')
      IF(WD(1:4).EQ.'FILE') THEN
      BFDEST = 4
      CALL NEXTFI('DISPlay-file=',OFILE)
      CALL ASSFIL(OFILE,BFUNIT,'WRITE','FORMATTED',ERROR)
      ELSE IF(WD(1:4).EQ.'REMA') THEN
      BFDEST = 3
      ELSE IF(WD(1:4).EQ.'DISP') THEN
      BFDEST = 2
      ELSE IF(WD(1:4).EQ.'OUTP') THEN
      BFDEST = 1
      ELSE
      CALL WRNDIE(-5,'BUFFER',
     & 'invalid buffer output direction')
      END IF
C ----------------------------------------
      ELSE IF (WD(1:4).EQ.'FROM') THEN
C
      CALL NEXTWD('FROM=')
      IF (WD(1:1).EQ.'=') CALL NEXTWD('FROM=')
      IF(WD(1:4).EQ.'BUFF') THEN
C read contents of another buffer into this one
      CALL NEXTWD('FROM-BUFFER=')
      IF (WD(1:1).EQ.'=') CALL NEXTWD('FROM-BUFFER=')
      MATCHLEN = WDLEN
C currently buffer names share one scope
      SCOPE = 1
      CALL DEFFINDREC(WD,MATCHLEN,SCOPE,4,.FALSE.,
     & FOUND,FRMLIST,QPART,.TRUE.)
      IF(FOUND.AND.(.NOT.QPART)) THEN
C skip over record containing name and destiation information
      FRMLIST = DEFMULTREC(FRMLIST)
      IF(FRMLIST.GT.0) THEN
C insert new-line signal if needed
      IF(DEFNAMLEN(LASTREC).EQ.0) THEN
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      DEFMULTREC(LASTREC) = RECNO
      LASTREC = RECNO
C signals new line
      DEFNAMLEN(LASTREC) = 1
      END IF
      END IF
C copy the contents of the list
      CALL DEFCPYLST(FRMLIST,RECNO)
      DEFMULTREC(LASTREC) = RECNO
C find end of list
      DO WHILE(DEFMULTREC(LASTREC).GT.0)
      LASTREC = DEFMULTREC(LASTREC)
      END DO
C
      END IF
      ELSE
      CALL WRNDIE(-5,'BUFFER',
     & 'named buffer not found')
      END IF
C
CC      ELSE IF(WD(1:4).EQ.'REMA') THEN
CC FROM remarks not yet implemented
CC      ELSE IF(WD(1:4).EQ.'FILE') THEN
CC FROM a file not yet implemented
      ELSE
      CALL WRNDIE(-5,'BUFFER',
     & 'invalid buffer input source')
      END IF
C ----------------------------------------
      ELSE IF (WD(1:4).EQ.'FLUS') THEN
      IF(OUTBUFLIST.GT.0) THEN
C write contents of buffer to appropriate unit
      CALL BFODUMP(OUTBUFLIST,BFDEST,BFUNIT)
C delete all but the first record in the list
      IF(DEFMULTREC(OUTBUFLIST).GT.0) THEN
      CALL DEFFREEREC(DEFMULTREC(OUTBUFLIST),.TRUE.)
      DEFMULTREC(OUTBUFLIST) = 0
      LASTREC = OUTBUFLIST
      END IF
      END IF
C ----------------------------------------
      ELSE IF (WD(1:4).EQ.'DUMP') THEN
      IF(OUTBUFLIST.GT.0) THEN
      CALL BFODUMP(OUTBUFLIST,BFDEST,BFUNIT)
      END IF
C ----------------------------------------
      ELSE IF (WD(1:4).EQ.'DISP') THEN
C buffered DISPLAY command
      IF(DEFNAMLEN(LASTREC).EQ.0) THEN
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      DEFMULTREC(LASTREC) = RECNO
      LASTREC = RECNO
C signals new line
      DEFNAMLEN(LASTREC) = 1
      END IF
      END IF
C
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      DEFMULTREC(LASTREC) = RECNO
      LASTREC = RECNO
      CALL LINSUB(DEFPARTXT(LASTREC),COMMAX,
     & DEFPARLEN(LASTREC),' ',0,1)
      END IF
C ----------------------------------------
      ELSE IF (WD(1:4).EQ.'LITE') THEN
C buffered LITERAL command
      IF(DEFNAMLEN(LASTREC).EQ.0) THEN
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      DEFMULTREC(LASTREC) = RECNO
      LASTREC = RECNO
C signals new line
      DEFNAMLEN(LASTREC) = 1
      END IF
      END IF
C
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      DEFMULTREC(LASTREC) = RECNO
      LASTREC = RECNO
      CALL LINSUB(DEFPARTXT(LASTREC),COMMAX,
     & DEFPARLEN(LASTREC),' ',0,2)
      END IF
C ----------------------------------------
      ELSE IF (WD(1:4).EQ.'CONC') THEN
C CONCatenate command
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      DEFMULTREC(LASTREC) = RECNO
      LASTREC = RECNO
      CALL LINSUB(DEFPARTXT(LASTREC),COMMAX,
     & DEFPARLEN(LASTREC),' ',0,1)
      END IF
C ----------------------------------------
      ELSE IF (WD(1:4).EQ.'NEXT') THEN
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      DEFMULTREC(LASTREC) = RECNO
      LASTREC = RECNO
C signals new line
      DEFNAMLEN(LASTREC) = 1
      END IF
C
      ELSE
      CALL CHKEND('BUFFER>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      END IF
      END IF
      IF(OUTBUFLIST.GT.0) THEN
      BODEST(OUTBUFLIST) = BFDEST
      BOUNIT(OUTBUFLIST) = BFUNIT
      END IF
      RETURN
      END
C=================================================================
      SUBROUTINE DEFGETFLG(SUBNAM,SUBLEN,QREMOVE,OFFSET,VARFLG)
C
C reads embedded flags from within symbol and parameter names
C this simplifies the code by moving the character
C handling of these commands into just two procedures
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(VARMAX) SUBNAM
      INTEGER SUBLEN, OFFSET
      LOGICAL QREMOVE,VARFLG(NUMVARFLG)
C local
      CHARACTER*(VARMAX) CPYNAM,TMPNAM
      INTEGER CPYLEN,I
      LOGICAL QFOUND
C begin
      DO I=1,NUMVARFLG
      VARFLG(I) = .FALSE.
      END DO
C
      CPYNAM = SUBNAM(1:SUBLEN)
      CPYLEN = SUBLEN
      QFOUND=.TRUE.
      DO WHILE(QFOUND.AND.(CPYLEN.GT.0))
      QFOUND=.FALSE.
C VARFLG # 1
C existence
      IF(CPYLEN.GT.6) THEN
      IF((CPYNAM(1:6).EQ.'EXIST%').OR.
     & ((CPYNAM(1:6).EQ.'EXIST_'))) THEN
      QFOUND = .TRUE.
      VARFLG(1) = .TRUE.
      TMPNAM = CPYNAM(1:CPYLEN)
      CPYNAM(1:CPYLEN-6)=TMPNAM(7:CPYLEN)
      CPYLEN = CPYLEN - 6
      END IF
      END IF
C VARFLG # 2
C stripping of quotes and/or type information from variables
      IF(CPYLEN.GT.6) THEN
      IF(CPYNAM(1:6).EQ.'STRIP%') THEN
      QFOUND = .TRUE.
      VARFLG(2) = .TRUE.
      TMPNAM = CPYNAM(1:CPYLEN)
      CPYNAM(1:CPYLEN-6)=TMPNAM(7:CPYLEN)
      CPYLEN = CPYLEN - 6
      END IF
      END IF
C VARFLG # 3
C immediate expansion of symbols/parameters
C -- only relevant in parameter blocks
      IF(CPYLEN.GT.7) THEN
      IF(CPYNAM(1:7).EQ.'EXPAND%') THEN
      QFOUND = .TRUE.
      VARFLG(3) = .TRUE.
      TMPNAM = CPYNAM(1:CPYLEN)
      CPYNAM(1:CPYLEN-7)=TMPNAM(8:CPYLEN)
      CPYLEN = CPYLEN - 7
      END IF
      END IF
C VARFLG # 4
C blank string/parameter check
      IF(CPYLEN.GT.6) THEN
      IF(CPYNAM(1:6).EQ.'BLANK%') THEN
      QFOUND = .TRUE.
      VARFLG(4) = .TRUE.
      TMPNAM = CPYNAM(1:CPYLEN)
      CPYNAM(1:CPYLEN-6)=TMPNAM(7:CPYLEN)
      CPYLEN = CPYLEN - 6
      END IF
      END IF
C VARFLG # 5
C quote directive
      IF(CPYLEN.GT.6) THEN
      IF(CPYNAM(1:6).EQ.'QUOTE%') THEN
      QFOUND = .TRUE.
      VARFLG(5) = .TRUE.
      TMPNAM = CPYNAM(1:CPYLEN)
      CPYNAM(1:CPYLEN-6)=TMPNAM(7:CPYLEN)
      CPYLEN = CPYLEN - 6
      END IF
      END IF
C extra $ or & (such as in &EXIST%&PARAMETER)
C blank string/parameter check
      IF(CPYLEN.GT.1) THEN
      IF((CPYNAM(1:1).EQ.'&').OR.
     & (CPYNAM(1:1).EQ.'$')) THEN
      QFOUND = .TRUE.
      TMPNAM = CPYNAM(1:CPYLEN)
      CPYNAM(1:CPYLEN-1)=TMPNAM(2:CPYLEN)
      CPYLEN = CPYLEN - 1
      END IF
      END IF
C
      END DO
      OFFSET = SUBLEN - CPYLEN
      IF(QREMOVE) THEN
      SUBLEN = CPYLEN
      IF(CPYLEN.GT.0) THEN
      SUBNAM(1:SUBLEN) = CPYNAM(1:SUBLEN)
      END IF
      END IF
      RETURN
      END
C=================================================================
      SUBROUTINE DEFADDFLG(SUBNAM,SUBLEN,VARFLG)
C
C adds embedded flags onto a symbol name
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(VARMAX) SUBNAM
      INTEGER SUBLEN
      LOGICAL VARFLG(NUMVARFLG)
C local
      CHARACTER*(VARMAX) CPYNAM,TMPNAM
      INTEGER CPYLEN
      LOGICAL OK
C begin
      OK = .TRUE.
      CPYLEN = SUBLEN
      CPYNAM = SUBNAM
C VARFLG # 1
C existence
      IF(VARFLG(1)) THEN
      IF((VARMAX-CPYLEN).GE.6) THEN
      TMPNAM = CPYNAM(1:CPYLEN)
      CPYNAM = 'EXIST%'//TMPNAM(1:CPYLEN)
      CPYLEN = CPYLEN + 6
      ELSE
      OK = .FALSE.
      END IF
      END IF
C VARFLG # 2
C striping of quotes and/or type information from variables
      IF(VARFLG(2)) THEN
      IF((VARMAX-CPYLEN).GE.6) THEN
      TMPNAM = CPYNAM(1:CPYLEN)
      CPYNAM = 'STRIP%'//TMPNAM(1:CPYLEN)
      CPYLEN = CPYLEN + 6
      ELSE
      OK = .FALSE.
      END IF
      END IF
C VARFLG # 3
C immediate expansion of symbols/parameters
C -- only relevant in parameter blocks
      IF(VARFLG(3)) THEN
      IF((VARMAX-CPYLEN).GE.7) THEN
      TMPNAM = CPYNAM(1:CPYLEN)
      CPYNAM = 'EXPAND%'//TMPNAM(1:CPYLEN)
      CPYLEN = CPYLEN + 7
      ELSE
      OK = .FALSE.
      END IF
      END IF
C VARFLG # 4
C blank string/parameter check
      IF(VARFLG(4)) THEN
      IF((VARMAX-CPYLEN).GE.6) THEN
      TMPNAM = CPYNAM(1:CPYLEN)
      CPYNAM = 'BLANK%'//TMPNAM(1:CPYLEN)
      CPYLEN = CPYLEN + 6
      ELSE
      OK = .FALSE.
      END IF
      END IF
C VARFLG # 5
C quote directive
      IF(VARFLG(5)) THEN
      IF((VARMAX-CPYLEN).GE.6) THEN
      TMPNAM = CPYNAM(1:CPYLEN)
      CPYNAM = 'QUOTE%'//TMPNAM(1:CPYLEN)
      CPYLEN = CPYLEN + 6
      ELSE
      OK = .FALSE.
      END IF
      END IF
C
      IF(.NOT.OK) THEN
      CALL WRNDIE(-5,'DEFADDFLG',
     & 'string too short - some flags skipped')
      END IF
      SUBLEN = CPYLEN
      IF(CPYLEN.GT.0) THEN
      SUBNAM(1:SUBLEN) = CPYNAM(1:SUBLEN)
      END IF
      RETURN
      END
C==============================================================
      SUBROUTINE DEFPROCINV(ABORT,ERR)
C
C procedure definition
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      LOGICAL ERR,ABORT
C local
      CHARACTER*(VARMAX) PARNAM
      LOGICAL DUVARFLG(NUMVARFLG)
      INTEGER PNMLEN, MATCHLEN, RECNO, PROREC, SCOPE
      LOGICAL QPART, FOUND
C begin
      PARNAM(1:WDLEN)=WD(1:WDLEN)
      PNMLEN = WDLEN
C see if it exists
      MATCHLEN = PNMLEN
      SCOPE = DEFCURSCOPE
      CALL DEFFINDREC(PARNAM,MATCHLEN,SCOPE,3,.FALSE.,
     & FOUND,PROREC,QPART,.FALSE.)
      IF(FOUND.AND..NOT.QPART) THEN
C procedure exists
C
C allocate space for current line in the storage buffer
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
C copy the current line into the buffer entry in front of
C any existing buffer entries
      DEFNEXTREC(RECNO) = DEFBUFLIST(NSTRM)
      DEFBUFLIST(NSTRM) = RECNO
      DEFPARLEN(RECNO) = COMLEN-CURSOR+1
      DEFPARTXT(RECNO) = COMLYN(CURSOR:COMLEN)
C signals echoing and buffering of this line
      DEFNAMLEN(RECNO) = 1
C
      CALL DEFNEWREC(RECNO)
C now make a blank line for the declaration parameter list
      IF(RECNO.GT.0) THEN
C
      DEFNEXTREC(RECNO) = DEFBUFLIST(NSTRM)
      DEFBUFLIST(NSTRM) = RECNO
C signals echoing and buffering of this line
      DEFNAMLEN(RECNO) = 1
      CURSOR = 1
C
C perform the substitution of the declaration part of the procedure
      CALL DEFSUBEXPA(1,CURSOR-2,PROREC,.FALSE.,
     & .FALSE.,.FALSE.,.FALSE.,1,.FALSE.,.FALSE.,DUVARFLG)
C
C force next line to be read from buffer
      CURSOR = COMLEN
C
C create a new scope for the procedure
      CALL DEFNEWSCOPE
C
C if this is an inline scope, unshield it for parameter definitions
      IF(DEFINLCNT(DEFCURSCOPE).NE.0) DEFINLCNT(DEFCURSCOPE) = 1
C
      CALL DEFNEXTW2(5)
      IF((WD(1:1).NE.'(').OR.(WDLEN.NE.1)) THEN
      WRITE (6,'(A)') ' %PROCEDURE-ERR: ( expected but not found.'
      ERR=.TRUE.
      ELSE
C specify scope/stream into which defines will be stored
      DEFSETSCOPE = DEFCURSCOPE
      DEFEVLSCOPE = DEFCURSCOPE
C parse declaration parameter block
      CALL DEFMACSET(ABORT,ERR,5)
      END IF
C
      CALL DEFNEXTW2(4)
      IF((WD(1:1).NE.'(').OR.(WDLEN.NE.1)) THEN
      WRITE (6,'(A)') ' %PROCEDURE-ERR: ( expected but not found.'
      ERR=.TRUE.
      ELSE
C specify scope/stream into which defines will be stored
      DEFEVLSCOPE = MAX(1,DEFCURSCOPE-1)
C parse invocation parameter block
      CALL DEFMACSET(ABORT,ERR,4)
C
      IF(.NOT.ABORT.AND..NOT.ERR) THEN
C
C perform the substitution of the procedure body
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
C copy the current line into the buffer entry in front of
C any existing buffer entries
      DEFNEXTREC(RECNO) = DEFBUFLIST(NSTRM)
      DEFBUFLIST(NSTRM) = RECNO
      DEFPARLEN(RECNO) = COMLEN
      DEFPARTXT(RECNO) = COMLYN(1:COMLEN)
C signals echoing and buffering of this line
      DEFNAMLEN(RECNO) = 1
C
C perform the substitution of the procedure body
      CALL DEFSUBEXPA(1,CURSOR-2,PROREC,.FALSE.,
     & .FALSE.,.FALSE.,.FALSE.,2,.FALSE.,.FALSE.,DUVARFLG)
C reactivate INLINE definition shielding (if appropriate)
      IF(DEFINLCNT(DEFCURSCOPE).EQ.1) DEFINLCNT(DEFCURSCOPE) = 2
C now let 'er rip!
      END IF
      ELSE
      CALL DEFKILSCOPE
      END IF
      END IF
C
C force next line to be read from buffer
      CURSOR = COMLEN
C
      END IF
      END IF
      END IF
      RETURN
      END
C=================================================================
      SUBROUTINE DEFADDSCO(PARNAM,PARLEN,SCOPE,QPREFIX)
C
C add scoping information onto a parameter string name
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER SCOPE, PARLEN
      CHARACTER*(VARMAX) PARNAM
      LOGICAL QPREFIX
C local
      CHARACTER*(VARMAX) PARCPY, CCHR
C begin
      PARCPY(1:PARLEN)=PARNAM(1:PARLEN)
      IF(QPREFIX) THEN
      PARCPY(1:PARLEN)=PARNAM(1:PARLEN)
      CCHR = PARCPY(1:1)
      WRITE(PARNAM,'(2A1,I1,A1,A)') CCHR,'_',SCOPE,
     & '_',PARCPY(2:PARLEN)
      PARLEN=3+PARLEN
      ELSE
      WRITE(PARNAM,'(A1,I1,A1,A)')
     & '_',SCOPE,'_',PARCPY(1:PARLEN)
      PARLEN=3+PARLEN
      END IF
      RETURN
      END
C=================================================================
      SUBROUTINE DEFGETSCO(PARNAM,PNMLEN,SCOPECNT,QSCOPED)
C
C extract scoping information from a symbol or parameter
C check to see if the scoping is valid - and return the scope
C (but leave it in place - undisturbed)
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER SCOPECNT
      LOGICAL QSCOPED
      CHARACTER*(VARMAX) PARNAM
      INTEGER PNMLEN
C begin
      QSCOPED = .FALSE.
C check to see if the parameter is explicitly scoped
C if so, then extract the scope index from the name
      SCOPECNT = 0
C explicitly scoped parameters
      IF(PNMLEN.GT.3) THEN
      IF((PARNAM(1:1).EQ.'_').AND.(PARNAM(3:3).EQ.'_')) THEN
      SCOPECNT = ICHAR(PARNAM(2:2))-ICHAR('0')
      IF((SCOPECNT.GE.1).OR.
     & (SCOPECNT.LE.DEFCURSCOPE)) THEN
C scope index is valid - use this scope
      QSCOPED = .TRUE.
      END IF
C explicitly scoped symbols
      ELSE IF((PNMLEN.GT.4).AND.(PARNAM(1:1).EQ.'$')) THEN
      IF((PARNAM(2:2).EQ.'_').AND.(PARNAM(4:4).EQ.'_')) THEN
      SCOPECNT = ICHAR(PARNAM(3:3))-ICHAR('0')
      IF((SCOPECNT.GE.1).OR.
     & (SCOPECNT.LE.DEFCURSCOPE)) THEN
C scope index is valid - use this scope
      QSCOPED = .TRUE.
      END IF
      END IF
      END IF
      END IF
      RETURN
      END
C=================================================================
      SUBROUTINE DEFREMSCO(SUBNAM,SUBLEN)
C
C removes scoping from a symbol or parameter
C (expects symbols to be proceeded by $ and
C parameters to be proceeded by nothing)
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(VARMAX) SUBNAM
      INTEGER SUBLEN
C local
      CHARACTER*(VARMAX) CPYNAM
C begin
C for symbols
      IF(SUBLEN.GT.1) THEN
      IF((SUBNAM(1:1).EQ.'$')) THEN
C
      IF(SUBLEN.GT.4) THEN
      IF((SUBNAM(2:2).EQ.'_').AND.(SUBNAM(4:4).EQ.'_')) THEN
      CPYNAM(1:SUBLEN) = SUBNAM(1:SUBLEN)
      SUBNAM(2:SUBLEN-3) = CPYNAM(5:SUBLEN)
      SUBLEN = SUBLEN - 3
      END IF
      END IF
C
      ELSE IF(SUBLEN.GT.3) THEN
      IF((SUBNAM(1:1).EQ.'_').AND.(SUBNAM(3:3).EQ.'_')) THEN
C for parameters
      CPYNAM(1:SUBLEN) = SUBNAM(1:SUBLEN)
      SUBNAM(1:SUBLEN-3) = CPYNAM(4:SUBLEN)
      SUBLEN = SUBLEN - 3
      END IF
      END IF
      END IF
      RETURN
      END
C=================================================================
      SUBROUTINE DEFNEWSCOPE
C
C creates and initializes a new scope, adjusts DEFCURSCOPE
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C begin
      IF(DEFCURSCOPE.LT.DEFSCPMAX) THEN
      DEFCURSCOPE = DEFCURSCOPE + 1
      DEFINLCNT(DEFCURSCOPE) = 0
      ELSE
C (some left-over debugging code was removed here)
      CALL WRNDIE(-5,'DEFNEWSCOPE',
     & 'No more scopes available')
      END IF
C
      IF(QINLINE) THEN
      QINLINE = .FALSE.
      DEFINLCNT(DEFCURSCOPE) = 2
      END IF
C
      RETURN
      END
C=================================================================
      SUBROUTINE DEFKILSCOPE
C
C destroys the current scope, adjusts DEFCURSCOPE
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C local
      INTEGER I
C begin
C can't destroy root scope
      IF(DEFCURSCOPE.GT.1) THEN
C free the storage
      DO I = 1,DEFTYPMAX
      IF(DEFVARLIST(DEFCURSCOPE,I).GT.0) THEN
      CALL DEFFREEREC(DEFVARLIST(DEFCURSCOPE,I),.TRUE.)
      DEFVARLIST(DEFCURSCOPE,I) = 0
      END IF
      END DO
C reduce the scope
      DEFCURSCOPE = DEFCURSCOPE - 1
      END IF
      RETURN
      END
C=================================================================
      SUBROUTINE DEFGETLIN(CC,CCMAX,CCLEN,QEXPAND,QBUFFER)
C
C if there are buffered lines to be processed, this
C command returns QEXPAND = .TRUE. and a line of
C commands
C
C otherwise, returns QEXPAND = .FALSE.
C
C Author: Warren L. DeLano and Axel T. Brunger
C ============================================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) CC
      INTEGER CCMAX, CCLEN
      LOGICAL QEXPAND,QBUFFER
C local
      INTEGER RECNO
C begin
C
      QEXPAND = .FALSE.
      IF(DEFBUFLIST(NSTRM).NE.0) THEN
C buffered lines are available
      RECNO=DEFBUFLIST(NSTRM)
      IF((QBUFFER.AND.(DEFNAMLEN(RECNO).NE.0)).OR.
     &   (DEFNAMLEN(RECNO).EQ.0)) THEN
      QEXPAND = .TRUE.
C copy the buffer into the command line
C free up the current line, and return
      DEFBUFLIST(NSTRM)=DEFNEXTREC(RECNO)
C
C IMPORTANT CONFUSION AVOIDANCE:  please note that the
C the PARTXT and PARLEN fields of parameter list records
C are used to store buffered lines of commands, and
C the NAMLEN record is used to store an echoing flag for
C this line (non-zero => echo line)
C
      IF(DEFPARLEN(RECNO).GT.CCMAX) DEFPARLEN(RECNO)=CCMAX
      CCLEN=DEFPARLEN(RECNO)
      CC(1:CCLEN)=DEFPARTXT(RECNO)(1:CCLEN)
      CALL DEFFREEREC(RECNO,.FALSE.)
      END IF
      END IF
      RETURN
      END
C=================================================================
      SUBROUTINE DEFCLOSE(MSTRM)
C
C frees up any and all storage associated with a particular
C scope/stream
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER MSTRM
C begin
C
C free buffered lines (normally, there shouldn't be any)
      IF(DEFBUFLIST(MSTRM).GT.0) THEN
      CALL DEFFREEREC(DEFBUFLIST(MSTRM),.TRUE.)
      END IF
C
C initialize this stream for next time
      DEFBUFLIST(MSTRM)=0
      RETURN
      END
C
C=================================================================
      SUBROUTINE DEFSUBPAR(CC,CCMAX,CCLEN,CURS,QDISPLAY)
C
C substitutes only if a parameter is located at
C the current cursor position (CURS)
C
C Author: Warren L. DeLano and Axel T. Brunger
C ============================================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) CC
      INTEGER CCMAX, CCLEN
      INTEGER CURS
      LOGICAL QDISPLAY
C local
      INTEGER RECNO
      LOGICAL QDONE
      INTEGER DEFCURS
      CHARACTER*(COMMAX) CHKSTR
      INTEGER CHKLEN,ICHR
      LOGICAL QCHANGED,QSEEK
      CHARACTER CCHR
C begin
      QCHANGED = .TRUE.
      DO WHILE (QCHANGED)
      QCHANGED = .FALSE.
      CHKSTR = CC
      CHKLEN = CCLEN
C
      DEFCURS = CURS
      QDONE = .FALSE.
      QSEEK = .FALSE.
      DEFCURS = CURS
      DO WHILE (.NOT.QDONE.AND.(DEFCURS.LE.CCLEN))
C keep looping through the word until a delimiter or the
C end of the string is reached
C
C if this appears to be a parameter, then enable seeking
C for other parameters which may be indices in the parameter name
      IF (CC(DEFCURS:DEFCURS).EQ.'&') QSEEK = .TRUE.
C
      IF ((CC(DEFCURS:DEFCURS).EQ.'&').OR.
     &    ((CC(DEFCURS:DEFCURS).EQ.'$').AND.
     &     QSCOPING)) THEN
C
C check for special symbols
      IF((CCLEN-DEFCURS).GE.2) THEN
      IF(CC(DEFCURS:DEFCURS+1).EQ.'&?') THEN
C &? is handled by miscom
      QDONE = .TRUE.
      ELSE IF(CC(DEFCURS:DEFCURS+1).EQ.'&%') THEN
C interpret directives
      IF(CC(DEFCURS+2:DEFCURS+2).EQ.'K') THEN
C kill scope
      CC(DEFCURS:DEFCURS+2) = '   '
      CALL DEFKILSCOPE
      QDONE=.TRUE.
      ELSE IF(CC(DEFCURS+2:DEFCURS+2).EQ.'N') THEN
C create new scope
      CC(DEFCURS:DEFCURS+2) = '   '
      CALL DEFNEWSCOPE
      QDONE=.TRUE.
      END IF
      END IF
      END IF
C
      IF (QDEF.AND.(.NOT.QDONE).AND.
     & ((ENDACT(ENDIND).EQ.'GO  ').OR.
     &  (ENDACT(ENDIND).EQ.'SCAN'))) THEN
C
C OKAY - we've found a symbol or a parameter, we're
C going to need to do something about it
C
C allocate space for current line in the storage buffer
      CALL DEFNEWREC(RECNO)
      IF(RECNO.EQ.0) THEN
      QDONE=.TRUE.
      ELSE
C copy the current line into the buffer entry in front of
C any existing buffer entries
      DEFNEXTREC(RECNO) = DEFBUFLIST(NSTRM)
      DEFBUFLIST(NSTRM) = RECNO
      DEFPARTXT(RECNO) = CC
      DEFPARLEN(RECNO) = CCLEN
C substitute
      CALL DEFSUBAT(DEFCURS,QDISPLAY)
C copy the buffer into the command line and free up the current record
      RECNO=DEFBUFLIST(NSTRM)
      DEFBUFLIST(NSTRM)=DEFNEXTREC(RECNO)
      IF(DEFPARLEN(RECNO).GT.CCMAX) DEFPARLEN(RECNO)=CCMAX
      CCLEN=DEFPARLEN(RECNO)
      CC(1:CCLEN)=DEFPARTXT(RECNO)(1:CCLEN)
      CALL DEFFREEREC(RECNO,.FALSE.)
C
      IF(CHKLEN.NE.CCLEN) THEN
      QCHANGED = .TRUE.
      ELSE IF(CHKSTR(1:CHKLEN).NE.CC(1:CHKLEN)) THEN
      QCHANGED = .TRUE.
      END IF
      IF(QCHANGED) QDONE = .TRUE.
      END IF
C            (recno.eq.0 )
C
      END IF
C            (qdef)
      END IF
C
C eliminate conditions other than expansion
      IF(QSCOPING) QDONE = .TRUE.
      IF(.NOT.QDEF) QDONE = .TRUE.
C if we are seeking then keep advancing through word
      IF(QSEEK.AND.(.NOT.QDONE).AND.(DEFCURS.LT.CCLEN)) THEN
      IF(CC(DEFCURS:DEFCURS+1).EQ.'$$') THEN
      DEFCURS = DEFCURS + 2
      ELSE IF(CC(DEFCURS:DEFCURS+1).EQ.'&&') THEN
      DEFCURS = DEFCURS + 2
      ELSE
      DEFCURS = DEFCURS + 1
      END IF
      IF(DEFCURS.LE.CCLEN) THEN
      CCHR=CC(DEFCURS:DEFCURS)
      ICHR=ASCIIM(ICHAR(CCHR))
      IF((ASCIIS(ICHAR(CCHR)).LT.10).OR.
     & (ICHR.EQ.ICHAR(' '))) QDONE = .TRUE.
      ELSE
      QDONE = .TRUE.
      END IF
      ELSE
      QDONE = .TRUE.
      END IF
C
      END DO
C            while(.not.qdone)
C
C
      END DO
C            while(qchanged)
      RETURN
      END
C==================================================================
      SUBROUTINE DEFSUBAT(DEFCURS,QDISPLAY)
C
C Checks current buffer line for a define parameter starting at
C position DEFCURS.  If found, then calls the substitution routine.
C
C Author: Warren L. DeLano and Axel T. Brunger
C ============================================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER DEFCURS
      LOGICAL QDISPLAY
C local
      INTEGER I,J
      LOGICAL FOUND
      CHARACTER*(VARMAX) PARNAM
      CHARACTER*(VARMAX) PARNAM2, TMPNAM
      INTEGER PNMLEN, PNMLEN2, PRESCOPE
      LOGICAL QDELIM
      INTEGER RECNO, OFFSET
      INTEGER CURLIN,SCOPE
      INTEGER ICHR
      CHARACTER CCHR
      LOGICAL QEXIST, QDOUBLE, QPART, QSKIP, QSCOPED, QEXPAND
      LOGICAL QSTRIP, QBLANK, QADDFLAG, QRESULT, QADDQUOT
      INTEGER MATCHLEN
      LOGICAL VARFLG(NUMVARFLG)
      LOGICAL DUOK
      CHARACTER*2 DUTYPE
      DOUBLE PRECISION DUDPVAL
      DOUBLE COMPLEX DUDCVAL
C begin
C
C initalize PARNAM
      PARNAM=' '
C
      I=DEFCURS
      CURLIN=DEFBUFLIST(NSTRM)
C
      IF((DEFPARTXT(CURLIN)(I:I).EQ.'&').OR.
     &   ((DEFPARTXT(CURLIN)(I:I).EQ.'$').AND.
     &    QSCOPING)) THEN
C
C check if we have a special double specification: &&
      IF ((DEFPARTXT(CURLIN)(I:I+1).EQ.'&&').OR.
     & (DEFPARTXT(CURLIN)(I:I+1).EQ.'$$')) THEN
      QDOUBLE=.TRUE.
      ELSE
      QDOUBLE=.FALSE.
      END IF
C looks like we have a define parameter
C
C delimit and capitalize the word
      PNMLEN=0
      IF (QDOUBLE) THEN
      J=I+2
      ELSE
      J=I+1
      END IF
      QDELIM=.FALSE.
      DO WHILE((PNMLEN.LT.VARMAX).AND.
     & (.NOT.QDELIM).AND.(J.LE.DEFPARLEN(CURLIN)))
C
      CCHR=DEFPARTXT(CURLIN)(J:J)
      ICHR=ICHAR(CCHR)
C
C delimiters are spaces, all non-printable characters, quotes,
C and all reserved characters defined by ASCIIS
      IF ((ASCIIM(ICHR).EQ.ICHAR(' ')).OR.(CCHR.EQ.'"')) THEN
      QDELIM=.TRUE.
      ELSE IF (ASCIIS(ICHR).LE.3) THEN
      QDELIM=.TRUE.
C enable &<paramater-name>&<parameter-name> (i.e. display &symbol&format )
      ELSE IF (CCHR.EQ.'&') THEN
      QDELIM=.TRUE.
      IF(PNMLEN.GT.0) THEN
C but don't break parameters after '_' or '%'
      IF((PARNAM(PNMLEN:PNMLEN).EQ.'_').OR.
     &   (PARNAM(PNMLEN:PNMLEN).EQ.'%')) QDELIM = .FALSE.
      END IF
      END IF
      IF(.NOT.QDELIM) THEN
      PNMLEN=PNMLEN+1
      PARNAM(PNMLEN:PNMLEN)=CHAR(ASCIIM(ICHR))
      J=J+1
      END IF
      END DO
C mrt pnmlen could be 0 : parnam(0:0) is invalid (valgrind)
      if( pnmlen .gt. 0 ) then
C ignore terminal '.' (i.e. with &sym1.&sym2 ...should only check &sym1 )
        IF(PARNAM(PNMLEN:PNMLEN).EQ.'.') THEN
          PNMLEN = PNMLEN - 1
        END IF
      endif
C
C see if parameter exists
      IF(PNMLEN.GT.0) THEN
      PARNAM2=PARNAM
      PNMLEN2=PNMLEN
      QSKIP = .FALSE.
      IF(QSCOPING) THEN
C currently we are only attaching a scope to the parameter...
C
C note and remove any explicit commands in the variable name
      CALL DEFGETFLG(PARNAM2,PNMLEN2,.TRUE.,OFFSET,VARFLG)
C handle double specification
      IF(QDOUBLE) VARFLG(2) = .TRUE.
C
      QEXPAND = VARFLG(3)
C check to see if the variable has already been scoped
      CALL DEFGETSCO(PARNAM2,PNMLEN2,PRESCOPE,QSCOPED)
C
      IF(QEXPAND) THEN
C The variable contains an explicit expansion directive...
C which means that substitution should occur immediately
C
      IF(DEFPARTXT(CURLIN)(I:I).EQ.'&') THEN
C for parameters, substitution will be handled below since QSKIP = .FALSE.
C so just put the directives back on to the name for processing later
      CALL DEFADDFLG(PARNAM2,PNMLEN2,VARFLG)
C
      ELSE IF(DEFPARTXT(CURLIN)(I:I).EQ.'$') THEN
C for symbols, do immediate symbol substitution right here using WDSUB
      QSKIP = .TRUE.
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
C restore the variable flags
      CALL DEFADDFLG(PARNAM2,PNMLEN2,VARFLG)
C put $ back on for substitution
      TMPNAM = PARNAM2(1:PNMLEN2)
      PARNAM2='$'//TMPNAM(1:PNMLEN2)
      PNMLEN2=PNMLEN2 + 1
C do substitution
      CALL WDSUB(PARNAM2,VARMAX,PNMLEN2,DUOK,DUTYPE,DUDPVAL,DUDCVAL)
      IF(.NOT.DUOK) THEN
      PNMLEN2 = 0
      END IF
C add required trailing space (will be removed by DEFSUBEXPA)
      DEFPARTXT(RECNO)(1:PNMLEN2+1)=PARNAM2(1:PNMLEN2)//' '
      DEFPARLEN(RECNO)=PNMLEN2+1
      CALL DEFSUBEXPA(I,PNMLEN,RECNO,.FALSE.,.FALSE.,
     & .FALSE.,.FALSE.,0,.FALSE.,.FALSE.,VARFLG)
      CALL DEFFREEREC(RECNO,.FALSE.)
      END IF
      END IF
      ELSE
C          if (.not.qexpand)
C just add scoping information
      QSKIP = .TRUE.
      IF(.NOT.QSCOPED) THEN
C parameter has not yet been scoped
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      I=I+1
C length adjustment - since we aren't replacing the $ or &
      IF(.NOT.QDOUBLE) PNMLEN = PNMLEN - 1
C attach scope
      CALL DEFADDSCO(PARNAM2,PNMLEN2,DEFEVLSCOPE,.FALSE.)
C restore the variable flags
      CALL DEFADDFLG(PARNAM2,PNMLEN2,VARFLG)
C add required trailing space (will be removed by DEFSUBEXPA)
      DEFPARTXT(RECNO)(1:PNMLEN2+1)=PARNAM2(1:PNMLEN2)//' '
      DEFPARLEN(RECNO)=PNMLEN2+1
      CALL DEFSUBEXPA(I,PNMLEN,RECNO,.FALSE.,.TRUE.,
     & .FALSE.,.FALSE.,0,.FALSE.,.FALSE.,VARFLG)
      CALL DEFFREEREC(RECNO,.FALSE.)
      END IF
      END IF
C            (.not.qscoped)
      END IF
      END IF
C
      IF(.NOT.QSKIP) THEN
C
C If qskip = .false. then we are doing a full substitution
C which means that we first need to find the parameter and then
C expand it, making whatever accomodations are necessary to
C squeeze it in.
C
C
C read and remove any and all embedded flags
      CALL DEFGETFLG(PARNAM2,PNMLEN2,.TRUE.,OFFSET,VARFLG)
C
C indexed symbol substitutions (use silent mode - in case the
C symbol doesn't exist it won't produce a warning message)
      CALL SYMIDX(PARNAM2,COMMAX,PNMLEN2,.TRUE.)
C handle double specification
      IF(QDOUBLE) VARFLG(2) = .TRUE.
C
      QEXIST = VARFLG(1)
      QSTRIP = VARFLG(2)
      QBLANK = VARFLG(4)
      QADDQUOT = VARFLG(5)
C
C intialize
      QADDFLAG = .FALSE.
C
C see if the parameter exists and is fully substituted
      MATCHLEN=PNMLEN2
      SCOPE = DEFCURSCOPE
      CALL DEFFINDREC(PARNAM2,MATCHLEN,SCOPE,2,.FALSE.,FOUND,
     & RECNO,QPART,.FALSE.)
C
      IF((FOUND.AND..NOT.QPART).OR.QPART) THEN
C check to see if the parameter contains a symbol or parameter...
      IF((DEFPARTXT(RECNO)(1:1).EQ.'&').OR.
     &   (DEFPARTXT(RECNO)(1:1).EQ.'$')) THEN
C it does -- so set the flag to propagate directives onto that symbol or parameter
      QADDFLAG = .TRUE.
      END IF
      END IF
C
      IF ((QEXIST.OR.QBLANK).AND..NOT.QADDFLAG.AND..NOT.QPART) THEN
C handling of EXIST% and BLANK% constructs
      IF (QEXIST) QRESULT = FOUND
      IF (QBLANK) CALL DEFCHKBLK(RECNO,QRESULT)
C replace symbol with TRUE or FALSE according to QRESULT
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      IF (QRESULT) THEN
      DEFPARTXT(RECNO)(1:5) = 'TRUE '
      DEFPARLEN(RECNO) = 5
      ELSE
      DEFPARTXT(RECNO)(1:6) = 'FALSE '
      DEFPARLEN(RECNO) = 6
      END IF
      CALL DEFSUBEXPA(I,PNMLEN,RECNO,.FALSE.,.FALSE.,
     & .FALSE.,.FALSE.,0,.FALSE.,.FALSE.,VARFLG)
      CALL DEFFREEREC(RECNO,.FALSE.)
      END IF
      ELSE
C
      IF(FOUND.AND..NOT.QPART) THEN
      IF(QADDQUOT.AND.DEFPARLEN(RECNO).GT.0) THEN
      IF(DEFPARTXT(RECNO)(1:1).EQ.'"') QADDQUOT = .FALSE.
      END IF
C normal parameter substitution mode
      CALL DEFSUBEXPA(I,PNMLEN,RECNO,QDOUBLE,QDISPLAY,
     & .FALSE.,QADDFLAG,0,QSTRIP,QADDQUOT,VARFLG)
      ELSE IF(QPART.AND.(.NOT.FOUND)) THEN
C partial parameter name substitution (propagate flags if appropriate)
      CALL DEFSUBEXPA(I,MATCHLEN+OFFSET,RECNO,QDOUBLE,QDISPLAY,
     & .TRUE.,QADDFLAG,0,.FALSE.,.FALSE.,VARFLG)
      END IF
      END IF
      END IF
      END IF
      END IF
      RETURN
      END
C
C==================================================================
      SUBROUTINE DEFSUBEXPA(PNMSTA,PNMLEN,RECNO1,QDOUBLE,
     & QDISPLAY,QINNAME,QADDFLAG,PROCMODE,QSTRIP,QADDQUOT,VARFLG)
C
C expands a define parameter in the current line
C
C if the parameter contains multiple lines, doesn't fit on the
C line, or the text following the parameter doesn't fit on the line,
C the extra text will be inserted into the define buffer for
C use by DEFCGTLN
C
C In QSTRIP mode, all double quotes will be removed from the
C substituted string.
C
C In QDOUBLE mode, the extra symbol will be removed
C
C In QADDFLAG mode, any directives in VARFLG will be appended onto
C the substituted variable name
C
C In QINNAME mode, the substitution is taking place in the middle
C of a parameter name, so if the substitution is going to require
C wrapping, then substitution needs to take the preceeding
C characters as well
C
C Authors: Warren L. DeLano and Axel T. Brunger
C =============================================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER PNMSTA,PNMLEN
      INTEGER RECNO1
      LOGICAL QDOUBLE, QSTRIP
      LOGICAL QDISPLAY,QINNAME, QADDQUOT
      LOGICAL QQINNAME, QADDFLAG, VARFLG(NUMVARFLG)
      INTEGER PROCMODE
C local
      INTEGER RECNO
      LOGICAL RGTFLG,ABORT,QSAVED,QBODY
      CHARACTER*(COMMAX) RGTTXT,FLGTXT
      CHARACTER*(COMMAX) NAMSAVTXT
      INTEGER NAMSAVLEN
      INTEGER CURLEN
      INTEGER LINLEN,RGTLEN
      INTEGER NXTLIN,CURLIN,PRVLIN
      INTEGER I, II, ICHR
      INTEGER RGTNEED, QUOTLEN
      INTEGER FLGLEN
      CHARACTER CCHR
C begin
C
C initialize
      RECNO=RECNO1
      ABORT=.FALSE.
      IF(PROCMODE.EQ.2) THEN
C body processing mode - to the procedure body
C
C skip first record - it has the procedure name and ( ...
      RECNO = DEFMULTREC(RECNO)
C
      IF(RECNO.GT.0) THEN
      QBODY=.FALSE.
      DO WHILE(.NOT.QBODY)
      IF((DEFNAMLEN(RECNO).EQ.4).AND.
     &   (DEFNAMTXT(RECNO)(1:4).EQ.'BODY')) THEN
      QBODY = .TRUE.
      ELSE
      RECNO = DEFMULTREC(RECNO)
      END IF
      END DO
      END IF
C safety
      IF(RECNO.EQ.0) THEN
      RECNO = RECNO1
      END IF
      END IF
C
      CURLIN=DEFBUFLIST(NSTRM)
      NXTLIN=DEFNEXTREC(CURLIN)
      PRVLIN=0
      FLGLEN=0
      QQINNAME = QINNAME
      NAMSAVLEN = 0
C
      IF(QADDFLAG) THEN
      FLGTXT = 'N'
      FLGLEN = 1
      CALL DEFADDFLG(FLGTXT,FLGLEN,VARFLG)
      FLGLEN = FLGLEN - 1
      END IF
C
      QUOTLEN = 0
      IF(QADDQUOT) THEN
C need to save space for opening quote
      QUOTLEN = 2
      END IF
C
C how much of current line is used?
      CURLEN=PNMSTA-1
C
C how much text is there to the right of the
C parameter name (if any) ?
      IF (QDOUBLE) THEN
      RGTLEN=DEFPARLEN(CURLIN)-(PNMSTA+PNMLEN)-1
      ELSE
      RGTLEN=DEFPARLEN(CURLIN)-(PNMSTA+PNMLEN)
      END IF
      IF(RGTLEN.GT.0) THEN
      RGTFLG=.TRUE.
C copy text to right of parameter name
C into temporary storage
      IF (QDOUBLE) THEN
      RGTTXT(1:RGTLEN)=
     &  DEFPARTXT(CURLIN)(PNMSTA+PNMLEN+2:DEFPARLEN(CURLIN))
      ELSE
      RGTTXT(1:RGTLEN)=
     &  DEFPARTXT(CURLIN)(PNMSTA+PNMLEN+1:DEFPARLEN(CURLIN))
      END IF
      ELSE
      RGTFLG=.FALSE.
      END IF
C
      RGTNEED = 0
      IF (QINNAME) RGTNEED = RGTLEN
C
C include everything up to & sign
      DEFPARLEN(CURLIN)=CURLEN
C
C keep looping until the parameter has been fully sustituted
C and all the of the rest of the text has been inserted
C on to the line or in to the buffer
      DO WHILE(((RECNO.GT.0).OR.RGTFLG).AND.(.NOT.ABORT))
C
      IF(CURLIN.EQ.0) THEN
C if we are starting a new line, allocate and link after previous
C line in buffer line list
      CALL DEFNEWREC(CURLIN)
      IF(CURLIN.GT.0) THEN
C transfer flags stored in DEFNAMLEN
      DEFNAMLEN(CURLIN) = DEFNAMLEN(PRVLIN)
      DEFNEXTREC(PRVLIN)=CURLIN
      CURLEN=0
      ELSE
C bail (if buffer is full)
      ABORT=.TRUE.
      END IF
      END IF
C
      IF(.NOT.ABORT) THEN
C
      IF(RECNO.GT.0) THEN
      LINLEN=DEFPARLEN(RECNO)
      IF (LINLEN.LE.0) THEN
      CALL WRNDIE(-5,'DEFSUBEXPA',
     & 'Attempted to use an undefined parameter.')
      RECNO = 0
      END IF
      END IF
C
      IF(RECNO.GT.0) THEN
C if we (still) need to insert text from the parameter
C
C special handling needed if we are expanding to multiple lines
C in the middle of a DISPLAY command
      IF(CURLEN.EQ.0) THEN
      IF(QDISPLAY) THEN
      DEFPARTXT(CURLIN)(1:8)='DISPLAY '
      CURLEN=8
      END IF
      IF((CURLEN+LINLEN+QUOTLEN+RGTNEED+FLGLEN).GT.COMMAX) THEN
      CALL WRNDIE(-5,'DEFSUBEXPA',
     & 'Insufficient room on the DISPLAY line to expand parameter.')
      END IF
C
      IF(QQINNAME.AND.NAMSAVLEN.GT.0) THEN
C we need to insert the beginning of the symbol/parameter
      QQINNAME = .FALSE.
      DO I=1,NAMSAVLEN
      DEFPARTXT(CURLIN)(I:I) =
     & NAMSAVTXT((NAMSAVLEN-I)+1:(NAMSAVLEN-I)+1)
      END DO
      CURLEN = NAMSAVLEN
C check space
      IF((CURLEN+LINLEN+QUOTLEN+RGTNEED+FLGLEN).GT.COMMAX) THEN
      CALL WRNDIE(-5,'DEFSUBEXPA',
     & 'Insufficient room on the line to expand parameter.')
      ENDIF
      END IF
      END IF
C            (curlen.eq.0)
      IF((CURLEN+LINLEN+QUOTLEN+RGTNEED+FLGLEN).LE.COMMAX) THEN
C fits
      IF (QSTRIP.AND..NOT.QADDFLAG) THEN
      II=CURLEN+1
      DO I=1,LINLEN
      IF (DEFPARTXT(RECNO)(I:I).NE.'"') THEN
      DEFPARTXT(CURLIN)(II:II)=DEFPARTXT(RECNO)(I:I)
      II=II+1
      END IF
      END DO
      CURLEN=CURLEN+II-CURLEN-1
      ELSE
C if we are preserving flags across substitution of a symbol
C or parameter, then here is where we insert them
      IF(FLGLEN.GT.0) THEN
      DEFPARTXT(CURLIN)(CURLEN+1:CURLEN+1) = DEFPARTXT(RECNO)(1:1)
      DEFPARTXT(CURLIN)(CURLEN+2:CURLEN+FLGLEN+1) = FLGTXT(1:FLGLEN)
      CURLEN = CURLEN + FLGLEN
      DEFPARTXT(CURLIN)(CURLEN+2:CURLEN+LINLEN)=
     &  DEFPARTXT(RECNO)(2:LINLEN-1)
      CURLEN=CURLEN+LINLEN
      FLGLEN=0
      ELSE
      IF(QUOTLEN.EQ.2) THEN
      DEFPARTXT(CURLIN)(CURLEN+1:CURLEN+1) = '"'
      QUOTLEN = 1
      CURLEN = CURLEN + 1
      END IF
C
      DEFPARTXT(CURLIN)(CURLEN+1:CURLEN+LINLEN)=
     &  DEFPARTXT(RECNO)(1:LINLEN)
      CURLEN=CURLEN+LINLEN
C
      IF(DEFMULTREC(RECNO).EQ.0) THEN
      IF(QUOTLEN.EQ.1) THEN
      DEFPARTXT(CURLIN)(CURLEN+1:CURLEN+1) =
     &  DEFPARTXT(CURLIN)(CURLEN:CURLEN)
      DEFPARTXT(CURLIN)(CURLEN:CURLEN) = '"'
      CURLEN = CURLEN + 1
      QUOTLEN = 0
      END IF
      END IF
C
      END IF
      END IF
C go onto next record
      RECNO=DEFMULTREC(RECNO)
C
      IF((PROCMODE.EQ.1).AND.(RECNO.GT.0)) THEN
C in procedure declaration insertion mode (PROCMODE=1)
C stop when 'BODY' record is reached
      IF((DEFNAMLEN(RECNO).EQ.4).AND.
     &   (DEFNAMTXT(RECNO)(1:4).EQ.'BODY')) THEN
      RECNO = 0
      END IF
C
      ELSE IF((PROCMODE.EQ.2).AND.(RECNO.EQ.0)) THEN
C in body mode, add a kill-scope directive after the last command in the body
      IF((COMMAX-CURLEN).LT.5) THEN
C won't fit on current line, so wrap onto next one
      DEFPARLEN(CURLIN)=CURLEN
      PRVLIN=CURLIN
      CALL DEFNEWREC(CURLIN)
      IF(CURLIN.GT.0) THEN
C transfer flags stored in DEFNAMLEN
      DEFNAMLEN(CURLIN) = DEFNAMLEN(PRVLIN)
      DEFNEXTREC(PRVLIN)=CURLIN
      CURLEN=0
      ELSE
      ABORT=.TRUE.
      END IF
      END IF
C place kill-scope directive
      IF((CURLIN.GT.0).AND.((COMMAX-CURLEN).GE.5)) THEN
      DEFPARTXT(CURLIN)(CURLEN+1:CURLEN+4)=' &%K '
      CURLEN = CURLEN + 5
      END IF
      END IF
C
C save current record length
      IF(CURLIN.GT.0) THEN
      DEFPARLEN(CURLIN)=CURLEN
      END IF
C
      IF(RECNO.GT.0) THEN
C separate lines as they were separated in initial define definition
      PRVLIN=CURLIN
      CURLIN=0
      END IF
      ELSE
C doesn't fit
      IF(QQINNAME) THEN
C special handling to make sure the front of the parameter
C name is kept on the same line when substituting inside
C a parameter name
C This routine scans backwards on the line until a separator
C is reached - note '&' is NOT a separator -
C (stored in reverse order)
      QSAVED=.FALSE.
      NAMSAVLEN=0
      DO WHILE((.NOT.QSAVED).AND.((CURLEN-NAMSAVLEN).GT.0))
      CCHR=DEFPARTXT(CURLIN)(CURLEN-NAMSAVLEN:CURLEN-NAMSAVLEN)
      ICHR=ASCIIM(ICHAR(CCHR))
      IF((ASCIIS(ICHAR(CCHR)).LT.10).OR.(ICHR.EQ.ICHAR(' '))) THEN
      QSAVED = .TRUE.
      ELSE
      NAMSAVLEN = NAMSAVLEN + 1
      NAMSAVTXT(NAMSAVLEN:NAMSAVLEN) = CCHR
      END IF
      END DO
      CURLEN=CURLEN-NAMSAVLEN
      END IF
C
C Remember the record number to link the next line after this one
      DEFPARLEN(CURLIN)=CURLEN
      PRVLIN=CURLIN
      CURLIN=0
      END IF
C
      ELSE IF(RGTFLG) THEN
C if we need to place text that followed the parameter name
C
C remove trailing space behind parameter when
C it is not the last thing on a line
      IF((CURLEN.GT.0).AND.((CURLEN+RGTLEN-1).LE.COMMAX)) THEN
      CURLEN=CURLEN-1
      END IF
C
C now add the remaning text from the right hand side of the parameter
      IF((CURLEN+RGTLEN).LE.COMMAX) THEN
C fits
      DEFPARTXT(CURLIN)(CURLEN+1:CURLEN+RGTLEN)=RGTTXT(1:RGTLEN)
      CURLEN=CURLEN+RGTLEN
      DEFPARLEN(CURLIN)=CURLEN
      RGTFLG=.FALSE.
      ELSE
C doesn't fit, so stick it on the next line
      PRVLIN=CURLIN
      CURLIN=0
C
      END IF
      END IF
C
      END IF
      END DO
C link to previous lines (if any exist)
      DEFNEXTREC(CURLIN)=NXTLIN
C
      RETURN
      END
C===============================================================
      SUBROUTINE DEFFINDPART(PARNAM,PNMLEN,RECNO,PREVREC)
C
C find (part of) a name in a linear list - no branching considered
C returns PREVREC so that the name can be excised from the list
C later if that is desired.  RECNO = 0 if not found
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C mrt: changed length of PARNAM from VARMAX to *
C     CHARACTER*(VARMAX) PARNAM
      INTEGER PNMLEN
      CHARACTER*(*) PARNAM
      INTEGER RECNO,PREVREC
C local
      LOGICAL FOUND
C begin
      FOUND = .FALSE.
      PREVREC = 0
      DO WHILE((.NOT.FOUND).AND.(RECNO.GT.0))
      IF(DEFNAMLEN(RECNO).EQ.PNMLEN) THEN
C parameter names must have the same length...
      IF(DEFNAMTXT(RECNO)(1:PNMLEN).EQ.PARNAM(1:PNMLEN)) THEN
C ...and the same text
      FOUND=.TRUE.
      END IF
      END IF
      IF(.NOT.FOUND) THEN
      PREVREC = RECNO
      RECNO = DEFNEXTREC(RECNO)
      END IF
      END DO
      RETURN
      END
C===============================================================
      SUBROUTINE DEFFINDREC(PARNAM,PNMLEN,SCOPE,VARTYPE,KILFLG,
     & FOUND,RECNO,QPART,QLIMIT)
C
C if KILFLG is false, track down the named parameter and return
C its index if it exists.
C
C if KILFLG is true, track down the named parameter, remove it
C from the relevant data structures, and free up the associated
C records
C
C if QLIMIT = true - limits search to a single scope
C
C WARNING: returns values in PNMLEN, SCOPE, FOUND, RECNO, QPART
C          ==> use local copies if you don't want them modified!
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      CHARACTER*(VARMAX) PARNAM
      INTEGER PNMLEN
      INTEGER VARTYPE
      INTEGER SCOPE,RECNO
      LOGICAL KILFLG,FOUND,QPART,QLIMIT
C local
      LOGICAL QPARSED,QSCOPED,QENDLOOP
      INTEGER PREVREC,NEXTREC,PARENT
      INTEGER CURLIST
      CHARACTER*(VARMAX) SUBNAM
      INTEGER SC,SUBLEN
      INTEGER SCOPECNT
      INTEGER PARENTLEN
C begin
C mrt : maybe uninitialized PARENTLEN
      parentlen = 0
      parent = 0
      QPART = .FALSE.
      QENDLOOP = .FALSE.
C check to see if the parameter is explicitly scoped
C if so, then extract the scope index from the name
      CALL DEFGETSCO(PARNAM,PNMLEN,SCOPECNT,QSCOPED)
C if no explicit scope was found, then use the scope
C requested by the calling function
      IF(.NOT.QSCOPED) SCOPECNT = SCOPE
C
      FOUND = .FALSE.
      SCOPECNT = SCOPECNT + 1
      DO WHILE((.NOT.FOUND).AND.(.NOT.QPART).AND.
     &         (SCOPECNT.GT.1).AND.(.NOT.QENDLOOP))
      SCOPECNT = SCOPECNT - 1
C this loop will search each enclosing scope until
C either the name is found or all scopes are exhausted
      CURLIST = 0
      QPARSED = .FALSE.
      PARENT = 0
      SC = 0
      SUBLEN = 0
C step through the name, seeking each substring in the
C list until the substrings are exhausted
      DO WHILE((.NOT.QPARSED).AND.(SC.LT.PNMLEN))
      SC = SC + 1
C is this the end of a delimited substring?
      IF((PARNAM(SC:SC).NE.'.').AND.(SC.LT.PNMLEN)) THEN
C no, not yet
      SUBLEN = SUBLEN + 1
      SUBNAM(SUBLEN:SUBLEN) = PARNAM(SC:SC)
      ELSE
      IF(SC.EQ.PNMLEN) THEN
C yes, but don't forget the last character of the last string
      SUBLEN = SUBLEN + 1
      SUBNAM(SUBLEN:SUBLEN) = PARNAM(SC:SC)
      END IF
C
C search for this substring
C
C need to get the relevant list from the proper
C place in the hierarchy
      IF(PARENT.EQ.0) THEN
C here we are just starting the search at the root level of
C this scope
      CURLIST = DEFVARLIST(SCOPECNT,VARTYPE)
      ELSE
C here we part-way down a compound parameter name.  The list
C we need to search is the one corresponding to the prefix
C so far
      CURLIST = DEFBRANREC(PARENT)
      END IF
C
      RECNO = CURLIST
C
C look for this part of the name
      IF((PARENT.EQ.0).AND.QSCOPED) THEN
C remove scoping if present
      CALL DEFREMSCO(SUBNAM,SUBLEN)
      CALL DEFFINDPART(SUBNAM(1:SUBLEN),SUBLEN,RECNO,PREVREC)
      ELSE
      CALL DEFFINDPART(SUBNAM(1:SUBLEN),SUBLEN,RECNO,PREVREC)
      END IF
C
      IF(RECNO.EQ.0) THEN
C substring was not found - so we need to check other scopes
      QPARSED = .TRUE.
      ELSE
C substring was found, what next?
      IF(SC.LT.PNMLEN) THEN
C well, if we haven't gotten to the end of the name, then
C we need to start over with next substring
      SUBLEN = 0
      PARENT = RECNO
      PARENTLEN = SC-1
      END IF
C
      END IF
C
      END IF
C            (at end of substring?)
      END DO
C            (qparsed)
C
C coming out of the DO loop above:
C   (1) if we have found the label then
C RECNO will be non-zero and PARENT, PREVREC,
C and CURLIST will all have meaningful values
C   (2) if we haven't found the lavel, RECNO
C will be zero and so we'll need to search
C other scopes or give up.
C
      IF(RECNO.GT.0) THEN
      FOUND = .TRUE.
      ELSE
      FOUND = .FALSE.
      END IF
C
      IF(FOUND) THEN
      IF(DEFBRANREC(RECNO).NE.0) THEN
C found record, but it is a compound parameter, so set QPART
      QPART = .TRUE.
      END IF
      END IF
C
      IF(.NOT.FOUND) THEN
C
      IF(PARENT.NE.0) THEN
C in case we have a partial match - meaning that we found some
C root part of the name and nothing else
      IF(DEFBRANREC(PARENT).EQ.0) THEN
C set the partial match flag - which says that we probably need
C to do a substitution with this prefix before getting something
C that we can match completely and then expand
      QPART = .TRUE.
      END IF
      END IF
C
      IF(QLIMIT) THEN
      QENDLOOP = .TRUE.
      END IF
C
      END IF
C            (.not.found)
      END DO
C
      IF(FOUND.AND.KILFLG) THEN
C it is an error to attempt to replace
C a compound parameter with an atomic one
      IF(DEFBRANREC(RECNO).GT.0) THEN
      CALL WRNDIE(-5,'DEFFINDREC',
     & 'parameter already defined as compound.')
      ENDIF
C first find the end of this MODULE parameter
      NEXTREC=DEFNEXTREC(RECNO)
C remove this parameter's link to the next parameter
      DEFNEXTREC(RECNO) = 0
C link the previous parameter to the next parameter (or 0)
      IF(PREVREC.GT.0) THEN
      DEFNEXTREC(PREVREC) = NEXTREC
C or handle special case if this is the first record
      ELSE IF(PARENT.EQ.0) THEN
C at lowest level of hierarchy
      DEFVARLIST(SCOPECNT,VARTYPE) = NEXTREC
      ELSE
C at intermediate level
      DEFBRANREC(PARENT) = NEXTREC
      END IF
C free up the memory
      CALL DEFFREEREC(RECNO,.TRUE.)
      RECNO = 0
C the parameter has now been removed from the list and the
C storage recovered
      END IF
C            (found.and.kilflg)
C
      IF(QPART.AND.(.NOT.FOUND)) THEN
C if a partial match was found,
C return the record
      RECNO=PARENT
C and the length of the match
      PNMLEN = PARENTLEN
      END IF
C return the scope in which the parameter was found
      SCOPE = SCOPECNT
      RETURN
      END
C===============================================================
      SUBROUTINE DEFINSREC(INSREC,SCOPE,VARTYPE)
C
C inserts a new parameter into the given scope
C and creates whatever parents are required to
C support it (ie. xxxx.yyyy.cccc = 1 ; will require
C xxxx and xxxx.yyyy before the data can be linked in
C as xxxx.yyyy.cccc )
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER INSREC,SCOPE,VARTYPE
C local
      CHARACTER*(VARMAX) PARNAM
      INTEGER PNMLEN
      LOGICAL QPARSED, QSCOPED
      INTEGER PREVREC,PARENT
      INTEGER CURLIST
      CHARACTER*(VARMAX) SUBNAM
      INTEGER SC,SUBLEN
      INTEGER RECNO
      INTEGER SCOPECNT
C begin
      PARNAM = DEFNAMTXT(INSREC)(1:DEFNAMLEN(INSREC))
      PNMLEN = DEFNAMLEN(INSREC)
      CURLIST = 0
      QPARSED = .FALSE.
      PARENT = 0
      SC = 0
      SUBLEN = 0
C check to see if the parameter is explicitly scoped
C if so, then extract the scope index from the name
      CALL DEFGETSCO(PARNAM,PNMLEN,SCOPECNT,QSCOPED)
C if no explicit SCOPECNT was found, then use the SCOPECNT
C requested by the calling function
      IF(.NOT.QSCOPED) SCOPECNT = SCOPE
C shield INLINE scopes
      CALL DEFGETSSC(SCOPECNT,SCOPECNT)
C remove explicit scoping from the name
      CALL DEFREMSCO(DEFNAMTXT(INSREC),DEFNAMLEN(INSREC))
C copy updated name
      PARNAM = DEFNAMTXT(INSREC)(1:DEFNAMLEN(INSREC))
      PNMLEN = DEFNAMLEN(INSREC)
C
C ignore terminal '.'
      IF(PARNAM(PNMLEN:PNMLEN).EQ.'.') THEN
      PNMLEN = PNMLEN - 1
      END IF
C step through the name, seeking each substring in the
C list until the substrings are exhausted
      DO WHILE((.NOT.QPARSED).AND.(SC.LT.PNMLEN))
      SC = SC + 1
C is this the end of a delimited substring?
      IF(PARNAM(SC:SC).NE.'.') THEN
C no, not yet
      SUBLEN = SUBLEN + 1
      SUBNAM(SUBLEN:SUBLEN) = PARNAM(SC:SC)
      ELSE
C
C search for this substring
C
C need to get list from proper place in the hierarchy
      IF(PARENT.EQ.0) THEN
C here just starting in the this SCOPECNT
      CURLIST = DEFVARLIST(SCOPECNT,VARTYPE)
      ELSE
C here part-way down a compound name
      CURLIST = DEFBRANREC(PARENT)
      END IF
C
      RECNO = CURLIST
      CALL DEFFINDPART(SUBNAM(1:SUBLEN),SUBLEN,RECNO,PREVREC)
C
      IF(RECNO.EQ.0) THEN
C substring not found - we need to create it
      CALL DEFNEWREC(RECNO)
C
      DEFNAMTXT(RECNO)(1:SUBLEN) = SUBNAM(1:SUBLEN)
      DEFNAMLEN(RECNO) = SUBLEN
C
      IF(PARENT.EQ.0) THEN
C at top of hierarchy
      DEFNEXTREC(RECNO) = DEFVARLIST(SCOPECNT,VARTYPE)
      DEFVARLIST(SCOPECNT,VARTYPE) = RECNO
      ELSE
C somewhere mysterious in the middle of it
      DEFNEXTREC(RECNO) = DEFBRANREC(PARENT)
      DEFBRANREC(PARENT) = RECNO
      END IF
C
      PARENT = RECNO
      ELSE
      PARENT = RECNO
      END IF
C            (recno = 0?)
C check to make sure this compound parameter
C hasn't already been used as an atomic one
      IF(DEFPARLEN(RECNO).NE.0) THEN
      CALL WRNDIE(-5,'DEFINSREC',
     & 'parameter already defined as atomic.')
      END IF
C
C prepare for next substring
      SUBLEN = 0
C
      END IF
C            (at end of substring?)
      END DO
C            (qparsed)
C
C coming out of the DO loop above:
C   PARENT will either be:
C (1) zero, in which case we need to append
C onto the start of the SCOPECNT block
C or (2) nonzero and linked to a list of
C parameters with the same prefix
C
C we need to attach the new record on at this point
C
C re-names the record with only the substring corresponding
C to the last element of the compound name
      DEFNAMTXT(INSREC)(1:SUBLEN) = SUBNAM(1:SUBLEN)
      DEFNAMLEN(INSREC) = SUBLEN
      IF(PARENT.EQ.0) THEN
C place it at the top of the hierarchy
      RECNO = DEFVARLIST(SCOPECNT,VARTYPE)
      CALL DEFFINDPART(SUBNAM(1:SUBLEN),SUBLEN,RECNO,PREVREC)
      IF(RECNO.GT.0) THEN
C -found a match
      IF(PREVREC.EQ.0) THEN
C --replace first entry
      DEFNEXTREC(INSREC) = DEFNEXTREC(RECNO)
      DEFVARLIST(SCOPECNT,VARTYPE) = INSREC
      ELSE
C --replace later entry
      DEFNEXTREC(INSREC) = DEFNEXTREC(RECNO)
      DEFNEXTREC(PREVREC) = INSREC
      END IF
C --free up old entry
      DEFNEXTREC(RECNO) = 0
      CALL DEFFREEREC(RECNO,.TRUE.)
      ELSE
C -no match found, so just insert
      DEFNEXTREC(INSREC) = DEFVARLIST(SCOPECNT,VARTYPE)
      DEFVARLIST(SCOPECNT,VARTYPE) = INSREC
      END IF
C
      ELSE
C place it somewhere in the middle
      RECNO = DEFBRANREC(PARENT)
      CALL DEFFINDPART(SUBNAM(1:SUBLEN),SUBLEN,RECNO,PREVREC)
      IF(RECNO.GT.0) THEN
C -found a match
      IF(PREVREC.EQ.0) THEN
C --replace first entry
      DEFNEXTREC(INSREC) = DEFNEXTREC(RECNO)
      DEFBRANREC(PARENT) = INSREC
      ELSE
C --replace later entry
      DEFNEXTREC(INSREC) = DEFNEXTREC(RECNO)
      DEFNEXTREC(PREVREC) = INSREC
      END IF
C --free up old entry
      DEFNEXTREC(RECNO) = 0
      CALL DEFFREEREC(RECNO,.TRUE.)
      ELSE
C -no match found, so just insert
      DEFNEXTREC(INSREC) = DEFBRANREC(PARENT)
      DEFBRANREC(PARENT) = INSREC
      END IF
C
      END IF
C
      RETURN
      END
C==================================================================
      SUBROUTINE DEFCDUMP(FIRLIN,LSTFRM,SCOPE,VARTYPE,ROOT,ROOTLEN)
C
C dumps a parameter list
C (for debugging and in response to ?)
C
C Authors: Warren L. DeLano and Axel T. Brunger
C =============================================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'symbol.inc'
      INTEGER FIRLIN,VARTYPE
      LOGICAL LSTFRM
      INTEGER SCOPE
      CHARACTER*(VARMAX) ROOT
      INTEGER ROOTLEN
C local
      INTEGER RECNO,I
      LOGICAL QFIRST
      CHARACTER*32 EQUTXT
      CHARACTER*(COMMAX) PARTXT
C NAMTXT needs to be longer than a normal symbol name
      CHARACTER*(COMMAX) NAMTXT
      INTEGER EQULEN,PARLEN,NAMLEN,NEXTREC
      INTEGER STKMAX,STKPTR
      PARAMETER (STKMAX=50)
      INTEGER STACK(STKMAX)
      LOGICAL QNEWVAR
C begin
      NAMLEN=0
      IF(LSTFRM) THEN
C list form
      EQUTXT='= '
      EQULEN=2
      ELSE
C or verbose form
      EQUTXT=' set to '
      EQULEN=8
      END IF
C need a stack to save the rest of the tree
C when we go down a particular branch
      STKPTR=0
      RECNO = FIRLIN
      DO WHILE((STKPTR.NE.0).OR.(RECNO.NE.0))
      IF(RECNO.EQ.0) THEN
      RECNO = STACK(STKPTR)
      STKPTR = STKPTR - 1
      ELSE
C
      DO WHILE(DEFBRANREC(RECNO).GT.0)
C check for overflow
      IF(STKPTR.EQ.STKMAX) THEN
      CALL WRNDIE(-5,'DEFCDUMP',
     & 'Stack blown -- too many parameters.')
      ELSE
      STKPTR=STKPTR+1
      END IF
C remember this branch
      STACK(STKPTR) = RECNO
      RECNO = DEFBRANREC(RECNO)
      END DO
      END IF
C
      IF(RECNO.GT.0) THEN
      NEXTREC = DEFNEXTREC(RECNO)
C
      IF(DEFBRANREC(RECNO).EQ.0) THEN
C only list entries for non-compound parameters
C
      QNEWVAR = .TRUE.
      DO WHILE(RECNO.NE.0)
      IF(DEFPARLEN(RECNO).GT.0) THEN
C for parameters with values to come
      PARLEN=DEFPARLEN(RECNO)
      PARTXT=DEFPARTXT(RECNO)(1:PARLEN)
      ELSE
      IF(VARTYPE.NE.1) THEN
      PARTXT='<value required>'
      PARLEN=16
      ELSE
      PARTXT=' '
      PARLEN=0
      END IF
      END IF
C
      IF(QNEWVAR) THEN
      QNEWVAR = .FALSE.
C only list name on first line of parameter
C
      NAMLEN = 0
C add prefix
      IF(VARTYPE.EQ.1) THEN
      NAMTXT(1:1)='$'
      NAMLEN = 1
      ELSE IF(VARTYPE.EQ.2) THEN
      NAMTXT(1:1)='&'
      NAMLEN = 1
      END IF
C
C collect parent substrs from the stack
      QFIRST = .TRUE.
      DO I=1,STKPTR
      IF(DEFBRANREC(STACK(I)).GT.0) THEN
      IF(.NOT.QFIRST) THEN
      NAMLEN=NAMLEN+1
      NAMTXT(NAMLEN:NAMLEN)='.'
      END IF
      IF(QFIRST.AND.(ROOTLEN.GT.0)) THEN
      NAMTXT(NAMLEN+1:ROOTLEN+NAMLEN) = ROOT(1:ROOTLEN)
      NAMLEN = NAMLEN + ROOTLEN
      ELSE
      NAMTXT(NAMLEN+1:NAMLEN+DEFNAMLEN(STACK(I)))=
     & DEFNAMTXT(STACK(I))(1:DEFNAMLEN(STACK(I)))
      NAMLEN=NAMLEN+DEFNAMLEN(STACK(I))
      END IF
      QFIRST = .FALSE.
      END IF
      END DO
C
C add substring for this final parameter
      IF(.NOT.QFIRST) THEN
      NAMLEN=NAMLEN+1
      NAMTXT(NAMLEN:NAMLEN)='.'
      END IF
      IF(QFIRST.AND.(ROOTLEN.GT.0)) THEN
      NAMTXT(NAMLEN+1:ROOTLEN+NAMLEN) = ROOT(1:ROOTLEN)
      NAMLEN = NAMLEN + ROOTLEN
      ELSE
      NAMTXT(NAMLEN+1:NAMLEN+DEFNAMLEN(RECNO)) =
     & DEFNAMTXT(RECNO)
      NAMLEN=NAMLEN+DEFNAMLEN(RECNO)
      END IF
      QFIRST = .FALSE.
C
C now add scope
      IF(VARTYPE.NE.3) THEN
      CALL DEFADDSCO(NAMTXT,NAMLEN,SCOPE,.TRUE.)
      ELSE
      CALL DEFADDSCO(NAMTXT,NAMLEN,SCOPE,.FALSE.)
      END IF
C
C now, write out the name
      IF (VARTYPE.EQ.1) THEN
C symbols - need type information too
      IF (CMDTYP(RECNO).EQ.'ST') THEN
      WRITE(6,'(5A)') ' ',
     & NAMTXT(1:NAMLEN),'="',PARTXT(1:PARLEN),
     & '" (string) '
C
      ELSE IF (CMDTYP(RECNO).EQ.'DP') THEN
      WRITE(6,'(3A,G14.6,A)') ' ',
     & NAMTXT(1:NAMLEN),'=',DPVALU(RECNO),' (real)'
C
      ELSE IF (CMDTYP(RECNO).EQ.'LO') THEN
      WRITE(6,'(5A)') ' ',
     & NAMTXT(1:NAMLEN),'="',PARTXT(1:PARLEN),
     & '" (logical) '
C
      ELSE IF (CMDTYP(RECNO).EQ.'DC') THEN
      WRITE(6,'(4A,G14.6,A,G14.6,A)') ' ',
     & NAMTXT(1:NAMLEN),'=',
     & ' (',DBLE(DCVALU(RECNO)),',',
     & DIMAG(DCVALU(RECNO)),')  (complex)'
C
      ELSE
      CALL DSPERR('WDSUB','Corrupt variable tables')
      END IF
C
      ELSE IF(VARTYPE.EQ.2) THEN
C normal parameters
      WRITE(6,'(2A,A,A)') ' ',NAMTXT(1:NAMLEN),
     & EQUTXT(1:EQULEN),PARTXT(1:PARLEN)
      ELSE IF(VARTYPE.EQ.3) THEN
      WRITE(6,'(2A,A,A)') ' PROCEDURE ',NAMTXT(1:NAMLEN),
     & EQUTXT(1:EQULEN),PARTXT(1:PARLEN)
      END IF
      ELSE IF(VARTYPE.NE.1) THEN
      WRITE(6,'(A,A)') ' ',PARTXT(1:PARLEN)
C
      END IF
      RECNO = DEFMULTREC(RECNO)
      END DO
      END IF
      RECNO=NEXTREC
      END IF
      END DO
C
      RETURN
      END
C==================================================================
      SUBROUTINE DEFCINIT
C
C define initialization
C
C Authors: Warren L. DeLano and Axel T. Brunger
C =============================================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C local
      INTEGER I,J
C begin
C initialize stream flags and pointers
      DO I=1,MXSTRM
      DEFBUFLIST(I)=0
      END DO
C initialize scope information
      DO I=1,DEFSCPMAX
      DO J=1,DEFTYPMAX
      DEFVARLIST(I,J)=0
      END DO
      DEFINLCNT(I)=0
      END DO
C
C initialize linked lists
      DO I=1,DBUFMX
      DEFNEXTREC(I)=I-1
      DEFBRANREC(I)=0
      DEFMULTREC(I)=0
      END DO
C
      DEFEMPTYREC=DBUFMX
C
C initialize the current scope
      DEFCURSCOPE=1
C
      RETURN
      END
C================================================================
      SUBROUTINE DEFNEWREC(RECNO)
C
C allocate new record in the array
C
C Authors: Warren L. DeLano and Axel T. Brunger
C =============================================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER RECNO
C begin
      IF(DEFEMPTYREC.GT.0) THEN
      RECNO=DEFEMPTYREC
      DEFEMPTYREC=DEFNEXTREC(RECNO)
C the zeroing below is VERY important
      DEFNEXTREC(RECNO)=0
      DEFMULTREC(RECNO)=0
      DEFBRANREC(RECNO)=0
      DEFPARLEN(RECNO)=0
      DEFNAMLEN(RECNO)=0
C
      ELSE
      CALL WRNDIE(-5,'DEFNEWREC','define buffer exhausted.')
      RECNO=0
      END IF
      RETURN
      END
C==================================================================
      SUBROUTINE DEFFREEREC(RECNO,QFRELISTS)
C
C free up a record or lists of records
C
C if QFRELISTS is .TRUE., then frees up all the
C records linked from this one using stack-based recursion
C
C Authors: Warren L. DeLano and Axel T. Brunger
C =============================================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER RECNO
      LOGICAL QFRELISTS
C local
      INTEGER STKMAX,STKPTR
      PARAMETER (STKMAX=50)
      INTEGER STACK(STKMAX)
      INTEGER NEXTTEMP
C begin
      IF(QFRELISTS) THEN
C frees up linked records as well
C
C put first record on the stack
      STKPTR = 1
      STACK(STKPTR)=RECNO
      RECNO = 0
C
      DO WHILE((STKPTR.GT.0).OR.(RECNO.NE.0))
C pop record off of stack
      IF(RECNO.EQ.0) THEN
      RECNO = STACK(STKPTR)
      STKPTR = STKPTR - 1
      END IF
C
      IF(RECNO.GT.0) THEN
C save branch entry
      IF(DEFBRANREC(RECNO).GT.0) THEN
C check for overflow
      IF(STKPTR.GE.STKMAX) THEN
      CALL WRNDIE(-5,'DEFFREEREC',
     & 'Stack blown -- too many parameters.')
      ELSE
      STKPTR=STKPTR+1
      END IF
C remember this branch
      STACK(STKPTR) = DEFBRANREC(RECNO)
      END IF
C break link
      DEFBRANREC(RECNO) = 0
C
C save multi-line entry
      IF(DEFMULTREC(RECNO).GT.0) THEN
C check for overflow
      IF(STKPTR.GE.STKMAX) THEN
      CALL WRNDIE(-5,'DEFFREEREC',
     & 'Stack blown -- too many parameters.')
      ELSE
      STKPTR=STKPTR+1
      END IF
C remember this multi-line entry
      STACK(STKPTR) = DEFMULTREC(RECNO)
      END IF
C break link
      DEFMULTREC(RECNO) = 0
C
C link this record into the empty list
      NEXTTEMP = DEFNEXTREC(RECNO)
      DEFNEXTREC(RECNO) = DEFEMPTYREC
      DEFEMPTYREC = RECNO
C go on to next record
      RECNO = NEXTTEMP
C
      END IF
      END DO
C
      ELSE
      DEFNEXTREC(RECNO)=DEFEMPTYREC
      DEFBRANREC(RECNO)=0
      DEFMULTREC(RECNO)=0
      DEFEMPTYREC=RECNO
      END IF
C
      RETURN
      END
C=================================================================
      SUBROUTINE DEFCNTRECS
C
C this is a routine for debugging purposes only (not normally called)
C it counts records to make sure we don't have a memory leak
C
C Authors: Warren L. DeLano and Axel T. Brunger
C =============================================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C local
      INTEGER RECNO
      INTEGER STKMAX,STKPTR
      PARAMETER (STKMAX=50)
      INTEGER STACK(STKMAX)
      INTEGER I,J,CNT
C begin
      CNT = 0
      STKPTR = 0
      RECNO = 0
      DO I=1,MXSTRM
      STKPTR = STKPTR + 1
      STACK(STKPTR)=DEFBUFLIST(I)
      END DO
      DO I=1,DEFCURSCOPE
      DO J=1,DEFTYPMAX
      STKPTR = STKPTR + 1
      STACK(STKPTR)=DEFVARLIST(I,J)
      END DO
      END DO
C
      STKPTR = STKPTR + 1
      STACK(STKPTR) = DEFEMPTYREC
C
      DO I=1,STKPTR
      END DO
C
      DO WHILE((STKPTR.GT.0).OR.(RECNO.NE.0))
C pop record off of stack
      IF(RECNO.EQ.0) THEN
      RECNO = STACK(STKPTR)
      STKPTR = STKPTR - 1
      END IF
C
      IF(RECNO.GT.0) THEN
C
      IF(DEFBRANREC(RECNO).GT.0) THEN
C check for overflow
      IF(STKPTR.GE.STKMAX) THEN
      CALL WRNDIE(-5,'DEFCNTRECS',
     & 'Stack blown -- too many parameters.')
      ELSE
      STKPTR=STKPTR+1
      END IF
C remember this branch
      STACK(STKPTR) = DEFBRANREC(RECNO)
      END IF
C
      IF(DEFMULTREC(RECNO).GT.0) THEN
C check for overflow
      IF(STKPTR.GE.STKMAX) THEN
      CALL WRNDIE(-5,'DEFCNTRECS',
     & 'Stack blown -- too many parameters.')
      ELSE
      STKPTR=STKPTR+1
      END IF
C remember this multi-line entry
      STACK(STKPTR) = DEFMULTREC(RECNO)
      END IF
C
      CNT = CNT + 1
      RECNO = DEFNEXTREC(RECNO)
C
      END IF
      END DO
C
      RETURN
      END
C==============================================================
      SUBROUTINE DEFNEXTW2(DEFCON)
C
C special interface to NEXTW2 which displays the
C appropriate prompt when parsing defines, modules,
C and procedure blocks
C
C see DEFMACSET below for DEFCON codes
C
C Author: Warren L. DeLano
C ========================
      IMPLICIT NONE
C input/output
      INTEGER DEFCON
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C begin
      IF(DEFCON.EQ.1) THEN
      CALL NEXTW2('MODULE-INVOCATION>')
      ELSE IF(DEFCON.EQ.2) THEN
      CALL NEXTW2('MODULE-DECLARATION>')
      ELSE IF(DEFCON.EQ.3) THEN
      CALL NEXTW2('DEFINE>')
      ELSE IF(DEFCON.EQ.4) THEN
      CALL NEXTW2('PROC-INVOCATION>')
      ELSE IF(DEFCON.EQ.5) THEN
      CALL NEXTW2('PROC-DECLARATION>')
      END IF
      RETURN
      END
C==============================================================
      SUBROUTINE DEFMACSET(ABORT,ERR,DEFCON)
C
C define definition parser
C
C Definition Control Level
C
C DEFCON = 5 parsing a procedure declaration parameter set
C DEFCON = 4 parsing a procedure invocation parameter set
C DEFCON = 3 standard define command
C DEFCON = 2 parsing a MODULE declaration parameter set
C DEFCON = 1 parsing a MODULE invocation parameter set
C
C Also note that the global variables DEFSETSCOPE and DEFEVLSCOPE
C need to be set appropriately before calling this function since
C they will determine what scope will be used for storing the
C parameters and what scope will be assigned to any parameters
C contained in the parameters being defined
C
C Also, this routine modifies QSUBS and QDEF, so they will need
C to be saved and restored by the enclosing subroutines
C
C Authors: Warren L. DeLano and Axel T. Brunger
C =============================================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      LOGICAL ABORT
      INTEGER DEFCON
C local
      CHARACTER*(VARMAX) PARNAM
      INTEGER PNMLEN,PARLIN
      INTEGER RECNO,LRECNO,NEWLEN,FIRSTL,SCOPE
      LOGICAL FOUND,ERR,EOFFLG,QPART
      INTEGER NPAREN
      INTEGER MATCHLEN
      INTEGER CRUNIT,CRSTRM
      LOGICAL QQSUBS,QQDEF
      LOGICAL CLOOP1, CLOOP2
C begin
C save qsubs and qdef
      QQSUBS=QSUBS
      QSUBS=.FALSE.
      QQDEF=QDEF
      QDEF=.FALSE.
C initialize
      CRUNIT=ISTRM(NSTRM)
      CRSTRM=NSTRM
      ABORT=.FALSE.
      EOFFLG=.FALSE.
C instead of substituting, force scoping of parameters
      QSCOPING=.TRUE.
C
      NPAREN=1
      CLOOP1=.TRUE.
      DO WHILE (CLOOP1)
C for parsing parameter name, turn off define substitution
      QDEF=.FALSE.
C and turn off symbol substitution (call to SYMIDX below will
C handle all the substitution)
      QSUBS = .FALSE.
      CALL DEFNEXTW2(DEFCON)
C
C check for stream switch in a MODULE file (BAD!)
      IF((DEFCON.EQ.2).AND.((CRUNIT.NE.ISTRM(NSTRM)).OR.
     &                     (NSTRM.NE.CRSTRM))) THEN
      EOFFLG=.TRUE.
      ERR=.TRUE.
      ELSE IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-define')
C
      ELSE IF((WD(1:1).EQ.')').AND.(WDLEN.EQ.1).AND.
     & (.NOT.QQUOT)) THEN
      NPAREN=NPAREN-1
      ELSE IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      IF((WD(1:1).EQ.'*').AND.(WDLEN.EQ.1).AND.
     & (.NOT.QQUOT)) THEN
      ABORT=.TRUE.
      ELSE IF((WD(1:1).EQ.'?').AND.
     & (WDLEN.EQ.1).AND.(.NOT.QQUOT)) THEN
      CALL DEFVARDUMP(' ',0,1,2)
      ELSE
C
      RECNO=0
      FIRSTL=0
      IF(WD(1:1).EQ.'&') THEN
      IF(WDLEN.GT.1) WD(1:WDLEN-1)=WD(2:WDLEN)
      WDLEN=WDLEN-1
      END IF
      IF(WDLEN.EQ.0) THEN
      CALL DSPERR('DEFMACSET',
     &            'parameter name must immediately follow &')
      ELSE
      IF(WDLEN.GT.VARMAX) THEN
      CALL DSPERR('DEFMACSET','parameter name too long; truncated')
      WDLEN=VARMAX
      END IF
C
C substitute any symbols that are present in the name
      CALL SYMIDX(WD,WDMAX,WDLEN,.TRUE.)
C
      PARNAM(1:WDLEN)=WD(1:WDLEN)
      PNMLEN=WDLEN
C
      FOUND=.FALSE.
      IF((DEFCON.EQ.1).OR.(DEFCON.EQ.4)) THEN
C verify that this parameter has been declared
C (invocation mode only)
      SCOPE = DEFSETSCOPE
      CALL DEFFINDREC(PARNAM,PNMLEN,SCOPE,2,.FALSE.,FOUND,
     & PARLIN,QPART,.FALSE.)
      IF(.NOT.FOUND) THEN
      ERR=.TRUE.
      IF(DEFCON.EQ.1) THEN
      CALL WRNDIE(-5,'MODULE',
     & 'undeclared module parameter: '//PARNAM(1:PNMLEN))
      ELSE IF(DEFCON.EQ.4) THEN
      CALL WRNDIE(-5,'PROCEDURE',
     & 'undeclared procedure parameter: '//PARNAM(1:PNMLEN))
      END IF
      END IF
      END IF
      IF(FOUND.OR.((DEFCON.NE.1).AND.(DEFCON.NE.4))) THEN
C
      QNEWLN=.FALSE.
      QCAPIT=.FALSE.
C turn on define substitution (with scoping)
      QDEF=.TRUE.
C turn off symbol substitution
      QSUBS = .FALSE.
      CALL NEXTW2(PARNAM(1:PNMLEN)//'=')
      IF(WD(1:1).NE.'=') CALL SAVEWD
      DONE=.FALSE.
      CLOOP2=.TRUE.
      DO WHILE (CLOOP2)
      CALL NEXTW2(PARNAM(1:PNMLEN)//'=')
C
C
C check for stream switch (bad!)
      IF((DEFCON.EQ.2).AND.((CRUNIT.NE.ISTRM(NSTRM))
     & .OR.(NSTRM.NE.CRSTRM))) THEN
      EOFFLG=.TRUE.
      ELSE
C
      IF((WD(1:1).EQ.'(').AND.(WDLEN.EQ.1).AND.
     & (.NOT.QQUOT)) THEN
      NPAREN=NPAREN+1
      ELSE IF((WD(1:1).EQ.')').AND.(WDLEN.EQ.1).AND.
     & (.NOT.QQUOT)) THEN
      NPAREN=NPAREN-1
      END IF
C
      IF(((WD(1:1).EQ.';').AND.(WDLEN.EQ.1).AND.
     & (.NOT.QQUOT)).OR.(NPAREN.LE.0)) THEN
      DONE=.TRUE.
      IF(RECNO.EQ.0) THEN
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      FIRSTL=RECNO
      DEFNAMTXT(RECNO)(1:PNMLEN)=PARNAM(1:PNMLEN)
      DEFNAMLEN(RECNO)=PNMLEN
      IF((DEFCON.EQ.1).OR.(DEFCON.EQ.4)) THEN
C
C blank parameter acceptable during invocation
      DEFPARLEN(RECNO)=1
      DEFPARTXT(RECNO)(1:1)=' '
      ELSE
C
C blank parameter during definition means required
C parameter upon invocation (signal: set length to -1 )
      DEFPARLEN(RECNO)=(-1)
      END IF
      ELSE
      ERR=.TRUE.
      END IF
      END IF
C
      ELSE IF(.NOT.ERR) THEN
      IF(QQUOT) THEN
C quoted strings will need to have quotes restored
      WDLEN=WDLEN+2
      END IF
      IF(WDLEN+1.GT.COMMAX) THEN
      CALL DSPERR('DEFMACSET','parameter too long')
      WDLEN=COMMAX-1
      END IF
C make sure we have a place to put this word
      IF(RECNO.EQ.0) THEN
      CALL DEFNEWREC(RECNO)
      IF(RECNO.GT.0) THEN
      FIRSTL=RECNO
      DEFNAMTXT(RECNO)(1:PNMLEN)=PARNAM(1:PNMLEN)
      DEFNAMLEN(RECNO)=PNMLEN
      ELSE
      ERR=.TRUE.
      END IF
      ELSE
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
      DEFPARTXT(RECNO)((DEFPARLEN(RECNO)+1):
     & (DEFPARLEN(RECNO)+WDLEN))=WD(1:WDLEN)
      END IF
      DEFPARLEN(RECNO)=DEFPARLEN(RECNO)+WDLEN
C add trailing space (expected by parser after last word)
      DEFPARLEN(RECNO)=DEFPARLEN(RECNO)+1
      DEFPARTXT(RECNO)(DEFPARLEN(RECNO):DEFPARLEN(RECNO))=' '
      END IF
C
      END IF
      QNEWLN=.FALSE.
      END IF
      IF (DONE.OR.(NPAREN.LE.0).OR.EOFFLG) CLOOP2=.FALSE.
      END DO
      DONE=.FALSE.
      QCAPIT=.TRUE.
C put this line on the top of the list
      IF(FIRSTL.GT.0) THEN
      IF(.NOT.ERR) THEN
C print new definition of parameter
      IF(((WRNLEV.GE.5).AND.((DEFCON.EQ.1).OR.(DEFCON.EQ.4)))
     & .OR.(WRNLEV.GE.10)) THEN
      CALL DEFCDUMP(FIRSTL,.FALSE.,DEFSETSCOPE,2,' ',0)
      END IF
C kill existing parameter with this name (if it exists)
      MATCHLEN=PNMLEN
      SCOPE = DEFSETSCOPE
      CALL DEFFINDREC(PARNAM,MATCHLEN,SCOPE,2,.TRUE.,
     & FOUND,PARLIN,QPART,.TRUE.)
C add new parameter on to top of parameter list
      CALL DEFINSREC(FIRSTL,DEFSETSCOPE,2)
      ELSE
C error condition
      CALL DEFFREEREC(FIRSTL,.TRUE.)
      END IF
      END IF
      END IF
      END IF
      END IF
      ELSE
      IF((WD(1:1).EQ.'(').AND.(WDLEN.EQ.1).AND.
     &   (.NOT.QQUOT)) THEN
      NPAREN=NPAREN+1
      ELSE IF((WD(1:1).EQ.')').AND.(WDLEN.EQ.1).AND.
     &   (.NOT.QQUOT)) THEN
      NPAREN=NPAREN-1
      END IF
      END IF
      IF ((NPAREN.LE.0).OR.EOFFLG) CLOOP1=.FALSE.
      END DO
C
      IF (ENDACT(ENDIND).EQ.'GO  ') THEN
      IF(((DEFCON.EQ.1).OR.(DEFCON.EQ.4)).AND.(.NOT.ERR)
     &   .AND.(.NOT.ABORT)) THEN
C if this is an invocation, check to make sure
C all required parameters have been defined
      RECNO=DEFVARLIST(DEFSETSCOPE,2)
      DO WHILE((RECNO.GT.0))
      IF((DEFNAMLEN(RECNO).GT.0).AND.(DEFPARLEN(RECNO).LT.0)) THEN
      ERR=.TRUE.
      IF(DEFCON.EQ.1) THEN
      WRITE(6,'(A,A,A)') ' %MODULE-ERR: ',
     & 'required parameter not defined: ',
     & DEFNAMTXT(RECNO)(1:DEFNAMLEN(RECNO))
      ELSE
      WRITE(6,'(A,A,A)') ' %PROCEDURE-ERR: ',
     & 'required parameter not defined: ',
     & DEFNAMTXT(RECNO)(1:DEFNAMLEN(RECNO))
      END IF
      END IF
      RECNO=DEFNEXTREC(RECNO)
      END DO
      END IF
      END IF
C
      DONE=.FALSE.
      QSCOPING=.FALSE.
      QSUBS=QQSUBS
      QDEF=QQDEF
      RETURN
      END
