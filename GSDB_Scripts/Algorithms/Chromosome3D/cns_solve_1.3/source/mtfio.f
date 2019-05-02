      SUBROUTINE SHOW
C
C prints current structure information.
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
C begin
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     &' Status of internal molecular topology database:'
      WRITE(6,'(6(A,I10),A)')
     &' -> NATOM= ',NATOM,'(MAXA=  ',MAXA,')  NBOND= ',NBOND,
     & '(MAXB=  ',MAXB,')'
      WRITE(6,'(6(A,I10),A)')
     &' -> NTHETA=',NTHETA,'(MAXT=  ',MAXT,')  NGRP=  ',NGRP,
     &  '(MAXGRP=',MAXGRP,')'
      WRITE(6,'(6(A,I10),A)')
     &' -> NPHI=  ',NPHI,'(MAXP=  ',MAXP,')  NIMPHI=',NIMPHI,
     &'(MAXIMP=',MAXIMP,')'
      WRITE(6,'(6(A,I10),A)')
     &' -> NNB=   ',NNB,'(MAXNB= ',MAXNB,') '
      END IF
      RETURN
      END
C
C=======================================================================
C
      SUBROUTINE MTFPRS
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
C
C subroutine parses STRUcture commands.
C
C local
      INTEGER NOLD
C
C begin
      CALL PUSEND('STRUcture>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('STRUcture>')
      CALL MISCOM('STRUcture>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-structure')
C
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      CALL SHOW
      ELSE IF (WD(1:4).EQ.'PSF ') THEN
      NOLD=NATOM
      CALL PSFRDR
      CALL ATMINI(NOLD+1,NATOM)
      CALL SCRATC
      CALL SHOW
      ELSE IF (WD(1:4).EQ.'DATA') THEN
      NOLD=NATOM
      CALL MTFRDR
      CALL ATMINI(NOLD+1,NATOM)
      CALL SCRATC
      CALL SHOW
      ELSE IF (WD(1:4).EQ.'RESE') THEN
      CALL MTFRES
      CALL SCRATC
      CALL SHOW
      ELSE
      CALL CHKEND('STRUcture>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
      RETURN
      END
C
C=======================================================================
C
      SUBROUTINE MTFRES
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
C
C this routine scratches the molecular topology
C
C Author: Axel T. Brunger
C
C begin
      NATOM=0
      NBOND=0
      NTHETA=0
      NPHI=0
      NIMPHI=0
      NNB=0
      NGRP=0
C
      RETURN
      END
C
C=======================================================================
C
      SUBROUTINE STRPRIMESC(STR, ESC)
      IMPLICIT NONE
C I/O
      CHARACTER*(*) STR, ESC
C
C Replace the prime character (') by \p.
C Author: Ralf W. Grosse-Kunstleve
C
C local
      INTEGER  I, J
      CHARACTER*1  BACKSL
C
C begin
      CALL SETBSL(BACKSL)
      ESC = ' '
      J = 1
      DO I = 1, LEN(STR)
        IF (J .GT. LEN(ESC)) RETURN
        IF (STR(I:I) .EQ. '''') THEN
          ESC(J:J) = BACKSL
          J = J + 1
          IF (J .GT. LEN(ESC)) RETURN
          ESC(J:J) = 'p'
        ELSE IF (STR(I:I) .EQ. BACKSL) THEN
          ESC(J:J) = BACKSL
          J = J + 1
          IF (J .GT. LEN(ESC)) RETURN
          ESC(J:J) = BACKSL
        ELSE
          ESC(J:J) = STR(I:I)
        END IF
        J = J + 1
      END DO
C
      RETURN
      END
C
C=======================================================================
C
      SUBROUTINE STRPRIMUNESC(STR, UNESC)
      IMPLICIT NONE
C I/O
      CHARACTER*(*) STR, UNESC
C
C Replace \p by prime character (').
C Author: Ralf W. Grosse-Kunstleve
C
C local
      INTEGER      I, J
      CHARACTER*1  BACKSL
C
C begin
      CALL SETBSL(BACKSL)
      UNESC = ' '
      I = 1
      J = 1
      DO WHILE (I .LT. LEN(STR))
        IF (J .GT. LEN(UNESC)) RETURN
        IF (      STR(I:I) .EQ. BACKSL
     &      .AND. (     STR(I+1:I+1) .EQ. 'p'
     &             .OR. STR(I+1:I+1) .EQ. 'P')) THEN
          UNESC(J:J) = ''''
          I = I + 1
        ELSE IF (      STR(I:I) .EQ. BACKSL
     &           .AND. STR(I+1:I+1) .EQ. BACKSL) THEN
          UNESC(J:J) = BACKSL
          I = I + 1
        ELSE
          UNESC(J:J) = STR(I:I)
        END IF
        I = I + 1
        J = J + 1
      END DO
C
      RETURN
      END
C
C=======================================================================
C
      SUBROUTINE STRCRUNCH(STR, ESTR)
      IMPLICIT NONE
C I/O
      CHARACTER*(*) STR
      INTEGER       ESTR
C
C Remove blanks from string.
C If portions of the string are surrounded by single quotes, only
C the trailing blanks inside the quotes are removed.
C Author: Ralf W. Grosse-Kunstleve
C
C local
      INTEGER  I, J
      INTEGER  LNB, LNBQ
      LOGICAL  QUOTED
C
C begin
      ESTR = LEN(STR)
      QUOTED = .FALSE.
      LNB = 0
      LNBQ = 0
      I = 1
      DO WHILE (I .LE. ESTR)
        IF (STR(I:I) .EQ. '''') THEN
          IF (QUOTED) THEN
            IF (      STR(LNBQ:LNBQ) .EQ. ''''
     &          .AND. STR(LNBQ+1:LNBQ+1) .EQ. ' ') LNBQ = LNBQ + 1
            IF (LNBQ + 1 .NE. I) THEN
              J = ESTR - LNBQ
              CALL TRIML(STR(LNBQ+1:ESTR), J)
              ESTR = LNBQ + J
              I = LNBQ + 1
            END IF
          END IF
          QUOTED = .NOT. QUOTED
          LNBQ = I
          LNB = I
          I = I + 1
        ELSE IF (QUOTED) THEN
          IF (STR(I:I) .NE. ' ') LNBQ = I
          LNB = I
          I = I + 1
        ELSE IF (STR(I:I) .EQ. ' ') THEN
          IF (I .EQ. ESTR) THEN
            ESTR = ESTR - 1
          ELSE
            IF (I .EQ. 1) I = 0
            J = ESTR - I
            CALL TRIML(STR(I+1:ESTR), J)
            ESTR = I + J
            I = I + 1
          END IF
        ELSE
          LNB = I
          I = I + 1
        END IF
      END DO
      ESTR = LNB
      IF (ESTR .EQ. 0) ESTR = 1
C
      RETURN
      END
C
C=======================================================================
C
      SUBROUTINE MTFWRT(FILE)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'numbers.inc'
      CHARACTER*(*) FILE
C
C Write MTF file.
C Author: Ralf W. Grosse-Kunstleve
C
C local
      INTEGER       I, U, EBUF
      LOGICAL       ERROR
      CHARACTER*128 BUF
      CHARACTER*8   ESEGID, ERESID, ERES, ETYPE, EIAC
C
C begin
      CALL ASSFIL(FILE, U, 'WRITE', 'FORMATTED', ERROR)
      IF (ERROR) RETURN
C
      IF (     2 * LEN(SEGID(1)) .GT. LEN(ESEGID)
     &    .OR. 2 * LEN(RESID(1)) .GT. LEN(ERESID)
     &    .OR. 2 * LEN(RES(1))   .GT. LEN(ERES)
     &    .OR. 2 * LEN(TYPE(1))  .GT. LEN(ETYPE)
     &    .OR. 2 * LEN(IAC(1))   .GT. LEN(EIAC) .OR. ERROR) THEN
        CALL WRNDIE(-5, 'MTFWRT', 'Fatal Coding Error')
        CALL DIE
      END IF
C
      WRITE(U,'(A/)',ERR=7) 'data_cns_mtf'
C
      WRITE(U,'(A)',ERR=7) '_cns_mtf.title'
      CALL WRTITL(U, 3)
C
      WRITE(U,'(A)',ERR=7) 'loop_'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_atom.id'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_atom.segment_id'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_atom.residue_id'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_atom.residue_name'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_atom.atom_name'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_atom.chemical_type'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_atom.charge'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_atom.atom_mass'
      DO I = 1, NATOM
        CALL STRPRIMESC(SEGID(I), ESEGID)
        CALL STRPRIMESC(RESID(I), ERESID)
        CALL STRPRIMESC(RES(I), ERES)
        CALL STRPRIMESC(TYPE(I), ETYPE)
        CALL STRPRIMESC(IAC(I), EIAC)
        WRITE(BUF, '(I8,11A,2G14.6)')
     &    I,
     &    ' ''', ESEGID, ''' ''',
     &           ERESID, ''' ''',
     &           ERES,   ''' ''',
     &           ETYPE,  ''' ''',
     &           EIAC,   ''' ',
     &           CG(I),
     &           AMASS(I)
        CALL STRCRUNCH(BUF, EBUF)
        WRITE(U, '(A)', ERR=7) BUF(1:EBUF)
      END DO
      WRITE(BUF, '(I8,11A,2G14.6)')
     &  -1,
     &  ' ''', ' ', ''' ''',
     &         ' ', ''' ''',
     &         ' ', ''' ''',
     &         ' ', ''' ''',
     &         ' ', ''' ',
     &         -ONE, -ONE
      CALL STRCRUNCH(BUF, EBUF)
      WRITE(U, '(A/)', ERR=7) BUF(1:EBUF)
C
      WRITE(U,'(A)',ERR=7) 'loop_'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_bond.id[1]'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_bond.id[2]'
      DO I = 1, NBOND
        WRITE(BUF, '(I8,1X,I8)') IB(I), JB(I)
        CALL STRCRUNCH(BUF, EBUF)
        WRITE(U, '(A)', ERR=7) BUF(1:EBUF)
      END DO
      WRITE(U, '(A/)', ERR=7) '-1 -1'
C
      WRITE(U,'(A)',ERR=7) 'loop_'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_angle.id[1]'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_angle.id[2]'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_angle.id[3]'
      DO I = 1, NTHETA
        WRITE(BUF, '(I8,1X,I8,1X,I8)') IT(I), JT(I), KT(I)
        CALL STRCRUNCH(BUF, EBUF)
        WRITE(U, '(A)', ERR=7) BUF(1:EBUF)
      END DO
      WRITE(U, '(A/)', ERR=7) '-1 -1 -1'
C
      WRITE(U,'(A)',ERR=7) 'loop_'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_dihedral.id[1]'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_dihedral.id[2]'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_dihedral.id[3]'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_dihedral.id[4]'
      DO I = 1, NPHI
        WRITE(BUF, '(I8,1X,I8,1X,I8,1X,I8)') IP(I),JP(I),KP(I),LP(I)
        CALL STRCRUNCH(BUF, EBUF)
        WRITE(U, '(A)', ERR=7) BUF(1:EBUF)
      END DO
      WRITE(U, '(A/)', ERR=7) '-1 -1 -1 -1'
C
      WRITE(U,'(A)',ERR=7) 'loop_'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_improper.id[1]'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_improper.id[2]'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_improper.id[3]'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_improper.id[4]'
      DO I = 1, NIMPHI
        WRITE(BUF, '(I8,1X,I8,1X,I8,1X,I8)') IM(I),JM(I),KM(I),LM(I)
        CALL STRCRUNCH(BUF, EBUF)
        WRITE(U, '(A)', ERR=7) BUF(1:EBUF)
      END DO
      WRITE(U, '(A/)', ERR=7) '-1 -1 -1 -1'
C
      WRITE(U,'(A)',ERR=7) 'loop_'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_explicit_nonbonded_exclusion.inb'
      DO I = 1, NNB
        WRITE(BUF, '(I8)') INB(I)
        CALL STRCRUNCH(BUF, EBUF)
        WRITE(U, '(A)', ERR=7) BUF(1:EBUF)
      END DO
      WRITE(U, '(A/)', ERR=7) '-1'
C
      WRITE(U,'(A)',ERR=7) 'loop_'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_explicit_nonbonded_exclusion.iblo'
      DO I = 1, NATOM
        WRITE(BUF, '(I8)') IBLO(I)
        CALL STRCRUNCH(BUF, EBUF)
        WRITE(U, '(A)', ERR=7) BUF(1:EBUF)
      END DO
      WRITE(U, '(A/)', ERR=7) '-1'
C
      WRITE(U,'(A)',ERR=7) 'loop_'
      WRITE(U,'(A)',ERR=7) '_cns_mtf_group_linked_list.first_atom_id'
      DO I = 1, NGRP
        WRITE(BUF, '(I8)') IGPBS(I)
        CALL STRCRUNCH(BUF, EBUF)
        WRITE(U, '(A)', ERR=7) BUF(1:EBUF)
      END DO
      WRITE(U, '(A/)', ERR=7) '-1'
C
      CALL VCLOSE(U,'KEEP',ERROR)
C
      GOTO 77
7     ERROR=.TRUE.
77    CONTINUE
C
      IF (ERROR) WRITE(6,'(A)') ' %MTFWRT-ERR: ERROR DURING WRITE'
C
      RETURN
      END
C
C=======================================================================
C
      subroutine RSTARTXT(U, mStr, nStr, Str, Error, EOF)
      IMPLICIT NONE
C I/O
      integer        U
      integer        mStr, nStr
      character*(*)   Str(*)
      logical        Error, EOF
C
C Read STAR text block.
C A STAR text block starts with a ; at the beginning of a line
C end ends with the next ; at the beginning of a line.
C
C local
      character*80  buf
C
C begin
      nStr = 0
      do while (.not. EOF)
        read(U, '(A)', ERR=7, END=6) buf
        if (nStr .eq. 0) then
          if (buf(1:1) .ne. ';') goto 7
        else if (buf(1:1) .eq. ';') then
          goto 9
        end if
        if (nStr .lt. mStr) then
          nStr = nStr + 1
          Str(nStr) = ' REMARKS' // buf(2:)
        end if
      end do
      goto 9
C
 6    EOF = .true.
      goto 9
 7    Error = .true.
C
 9    continue
      return
      end
C
C=======================================================================
C
      subroutine GetSingleKW(U, keyword, Error, EOF)
      IMPLICIT NONE
C I/O
      integer        U
      character*(*)  keyword
      logical        Error, EOF
C
C Read empty lines until a line with the keyword is found.
C
C local
      integer       l
      character*80  buf
C
C begin
      do while (.not. EOF)
        read(U, '(A)', ERR=7, END=6) buf
        l = len(buf)
        call TrimL(buf, l)
        if (l .ne. 0) then
          if (buf .eq. keyword) goto 9
          goto 7
        end if
      end do
      goto 9
C
 6    EOF = .true.
      goto 9
 7    Error = .true.
C
 9    continue
      return
      end
C
C=======================================================================
C
      SUBROUTINE MTFRDR
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'ctitla.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'numbers.inc'
C
C Read MTF file.
C Author: Ralf W. Grosse-Kunstleve
C
C local
      INTEGER           U, I, I1, I2, I3, I4
      INTEGER           NATO2, NBON2, NTHET2, NPH2, NIMPH2, NN2, NGR2
      INTEGER           NATO3
      CHARACTER*80      BUF
      CHARACTER*8       LSEGID, LRESID, LRES, LTYPE, LIAC
      DOUBLE PRECISION  LCG, LAMASS
C
C begin
      IF (     2 * LEN(SEGID(1)) .GT. LEN(LSEGID)
     &    .OR. 2 * LEN(RESID(1)) .GT. LEN(LRESID)
     &    .OR. 2 * LEN(RES(1))   .GT. LEN(LRES)
     &    .OR. 2 * LEN(TYPE(1))  .GT. LEN(LTYPE)
     &    .OR. 2 * LEN(IAC(1))   .GT. LEN(LIAC) .OR. ERROR) THEN
        CALL WRNDIE(-5, 'MTFRDR', 'Fatal Coding Error')
        CALL DIE
      END IF
C
      ERROR=.FALSE.
      EOF=.FALSE.
C
      U=ISTRM(NSTRM)
C
      CALL GetSingleKW(U, '_cns_mtf.title', ERROR, EOF)
      IF (ERROR .OR. EOF) GOTO 77
      CALL RSTARTXT(U, MXTITL, NTITLE, TITLE, ERROR, EOF)
      IF (ERROR .OR. EOF) GOTO 77
      CALL RDTITB(-1, ' ')
C
      CALL GetSingleKW(U, 'loop_', ERROR, EOF)
      IF (ERROR .OR. EOF) GOTO 77
      READ(U, '(///////A)', ERR=7, END=6) BUF
      IF (BUF .NE. '_cns_mtf_atom.atom_mass') GOTO 7
      NATO2 = 0
      I1 = 0
      DO WHILE (I1 .GE. 0)
        READ(U, *, ERR=7, END=6)
     &    I1, LSEGID, LRESID, LRES, LTYPE, LIAC, LCG, LAMASS
        IF (I1 .GE. 0) THEN
          IF (NATOM + NATO2 .GE. MAXA) THEN
            CALL WRNDIE(-5, 'MTFIO',
     &        'MAXA parameter exceeded --> recompile program')
          ELSE
            NATO2 = NATO2 + 1
            I = NATOM + NATO2
            CALL STRPRIMUNESC(LSEGID, SEGID(I))
            CALL STRPRIMUNESC(LRESID, RESID(I))
            CALL STRPRIMUNESC(LRES,   RES(I))
            CALL STRPRIMUNESC(LTYPE,  TYPE(I))
            CALL STRPRIMUNESC(LIAC,   IAC(I))
            CG(I)     = LCG
            AMASS(I)  = LAMASS
            QLOOKU(I) = .TRUE.
          END IF
        END IF
      END DO
C
      CALL GetSingleKW(U, 'loop_', ERROR, EOF)
      IF (ERROR .OR. EOF) GOTO 77
      READ(U, '(/A)', ERR=7, END=6) BUF
      IF (BUF .NE. '_cns_mtf_bond.id[2]') GOTO 7
      NBON2 = 0
      I1 = 0
      DO WHILE (I1 .GE. 0)
        READ(U, *, ERR=7, END=6) I1, I2
        IF (I1 .GE. 0) THEN
          IF (NBOND + NBON2 .GE. MAXB) THEN
            CALL WRNDIE(-5, 'MTFIO',
     &        'MAXB parameter exceeded --> recompile program')
          ELSE
            NBON2 = NBON2 + 1
            I = NBOND + NBON2
            IB(I)    = I1
            JB(I)    = I2
            QACBB(I) = .TRUE.
            QACBC(I) = .TRUE.
          END IF
        END IF
      END DO
C
      CALL GetSingleKW(U, 'loop_', ERROR, EOF)
      IF (ERROR .OR. EOF) GOTO 77
      READ(U, '(//A)', ERR=7, END=6) BUF
      IF (BUF .NE. '_cns_mtf_angle.id[3]') GOTO 7
      NTHET2 = 0
      I1 = 0
      DO WHILE (I1 .GE. 0)
        READ(U, *, ERR=7, END=6) I1, I2, I3
        IF (I1 .GE. 0) THEN
          IF (NTHETA + NTHET2 .GE. MAXT) THEN
            CALL WRNDIE(-5, 'MTFIO',
     &        'MAXT parameter exceeded --> recompile program')
          ELSE
            NTHET2 = NTHET2 + 1
            I = NTHETA + NTHET2
            IT(I)     = I1
            JT(I)     = I2
            KT(I)     = I3
            QACTB(I)  = .TRUE.
            QACTC(I)  = .TRUE.
            QACTUB(I) = .TRUE.
            QACTUC(I) = .TRUE.
            ACTUB(I)  = ZERO
            ACTUC(I)  = ZERO
          END IF
        END IF
      END DO
C
      CALL GetSingleKW(U, 'loop_', ERROR, EOF)
      IF (ERROR .OR. EOF) GOTO 77
      READ(U, '(///A)', ERR=7, END=6) BUF
      IF (BUF .NE. '_cns_mtf_dihedral.id[4]') GOTO 7
      NPH2 = 0
      I1 = 0
      DO WHILE (I1 .GE. 0)
        READ(U, *, ERR=7, END=6) I1, I2, I3, I4
        IF (I1 .GE. 0) THEN
          IF (NPHI + NPH2 .GE. MAXP) THEN
            CALL WRNDIE(-5, 'MTFIO',
     &        'MAXP parameter exceeded --> recompile program')
          ELSE
            NPH2 = NPH2 + 1
            I = NPHI + NPH2
            IP(I)    = I1
            JP(I)    = I2
            KP(I)    = I3
            LP(I)    = I4
            QACPB(I) = .TRUE.
            QACPC(I) = .TRUE.
            QACPD(I) = .TRUE.
          END IF
        END IF
      END DO
C
      CALL GetSingleKW(U, 'loop_', ERROR, EOF)
      IF (ERROR .OR. EOF) GOTO 77
      READ(U, '(///A)', ERR=7, END=6) BUF
      IF (BUF .NE. '_cns_mtf_improper.id[4]') GOTO 7
      NIMPH2 = 0
      I1 = 0
      DO WHILE (I1 .GE. 0)
        READ(U, *, ERR=7, END=6) I1, I2, I3, I4
        IF (I1 .GE. 0) THEN
          IF (NIMPHI + NIMPH2 .GE. MAXIMP) THEN
            CALL WRNDIE(-5, 'MTFIO',
     &        'MAXIMP parameter exceeded --> recompile program')
          ELSE
            NIMPH2 = NIMPH2 + 1
            I = NIMPHI + NIMPH2
            IM(I)    = I1
            JM(I)    = I2
            KM(I)    = I3
            LM(I)    = I4
            QACIB(I) = .TRUE.
            QACIC(I) = .TRUE.
            QACID(I) = .TRUE.
          END IF
        END IF
      END DO
C
      CALL GetSingleKW(U, 'loop_', ERROR, EOF)
      IF (ERROR .OR. EOF) GOTO 77
      READ(U, '(A)', ERR=7, END=6) BUF
      IF (BUF .NE. '_cns_mtf_explicit_nonbonded_exclusion.inb') GOTO 7
      NN2 = 0
      I1 = 0
      DO WHILE (I1 .GE. 0)
        READ(U, *, ERR=7, END=6) I1
        IF (I1 .GE. 0) THEN
          IF (NNB + NN2 .GE. MAXNB) THEN
            CALL WRNDIE(-5, 'MTFIO',
     &        'MAXNB parameter exceeded --> recompile program')
          ELSE
            NN2 = NN2 + 1
            I = NNB + NN2
            INB(I) = I1
          END IF
        END IF
      END DO
C
      CALL GetSingleKW(U, 'loop_', ERROR, EOF)
      IF (ERROR .OR. EOF) GOTO 77
      READ(U, '(A)', ERR=7, END=6) BUF
      IF (BUF .NE. '_cns_mtf_explicit_nonbonded_exclusion.iblo') GOTO 7
      NATO3 = 0
      I1 = 0
      DO WHILE (I1 .GE. 0)
        READ(U, *, ERR=7, END=6) I1
        IF (I1 .GE. 0) THEN
          IF (NATO3 .GE. NATO2) GOTO 7
          IF (NATOM + NATO3 .GE. MAXA) THEN
            CALL WRNDIE(-5, 'MTFIO',
     &        'MAXA parameter exceeded --> recompile program')
          ELSE
            NATO3 = NATO3 + 1
            I = NATOM + NATO3
            IBLO(I) = I1
          END IF
        END IF
      END DO
C
      CALL GetSingleKW(U, 'loop_', ERROR, EOF)
      IF (ERROR .OR. EOF) GOTO 77
      READ(U, '(A)', ERR=7, END=6) BUF
      IF (BUF .NE. '_cns_mtf_group_linked_list.first_atom_id') GOTO 7
      NGR2 = 0
      I1 = 0
      DO WHILE (I1 .GE. 0)
        READ(U, *, ERR=7, END=6) I1
        IF (I1 .GE. 0) THEN
          IF (NGRP + NGR2 .GE. MAXGRP) THEN
            CALL WRNDIE(-5, 'MTFIO',
     &        'MAXGRP parameter exceeded --> recompile program')
          ELSE
            NGR2 = NGR2 + 1
            I = NGRP + NGR2
            IGPBS(I) = I1
          END IF
        END IF
      END DO
      IGPBS(NGRP+NGR2+1)=NATOM+NATO2
C
      GOTO 77
C
C Error-labels
 7    ERROR=.TRUE.
      GOTO 77
 6    EOF=.TRUE.
C
 77   CONTINUE
C
      IF (ERROR) WRITE(6,'(A)') ' %MTFRDR-ERR: error during read'
      IF (EOF)   WRITE(6,'(A)') ' %MTFRDR-ERR: EOF during read'
C
      IF (.NOT.ERROR.AND..NOT.EOF) THEN
        IF (NATOM.GT.0) THEN
C         shift the atom indices
          DO I=NBOND+1,NBOND+NBON2
            IB(I)=IB(I)+NATOM
            JB(I)=JB(I)+NATOM
          END DO
          DO I=NTHETA+1,NTHETA+NTHET2
            IT(I)=IT(I)+NATOM
            JT(I)=JT(I)+NATOM
            KT(I)=KT(I)+NATOM
          END DO
          DO I=NPHI+1,NPHI+NPH2
            IP(I)=IP(I)+NATOM
            JP(I)=JP(I)+NATOM
            KP(I)=KP(I)+NATOM
            LP(I)=LP(I)+NATOM
          END DO
          DO I=NIMPHI+1,NIMPHI+NIMPH2
            IM(I)=IM(I)+NATOM
            JM(I)=JM(I)+NATOM
            KM(I)=KM(I)+NATOM
            LM(I)=LM(I)+NATOM
          END DO
          DO I=NNB+1,NNB+NN2
            INB(I)=INB(I)+NATOM
          END DO
          DO I=NATOM+1,NATOM+NATO2
            IBLO(I)=IBLO(I)+IBLO(NATOM)
          END DO
          DO I=NGRP+1,NGRP+NGR2
            IGPBS(I)=IGPBS(I)+NATOM
          END DO
        END IF
C
C       update the max-counters
        NATOM=NATOM+NATO2
        NBOND=NBOND+NBON2
        NTHETA=NTHETA+NTHET2
        NPHI=NPHI+NPH2
        NIMPHI=NIMPHI+NIMPH2
        NNB=NNB+NN2
        NGRP=NGRP+NGR2
      END IF
C
      RETURN
      END
