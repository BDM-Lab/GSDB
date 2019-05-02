      SUBROUTINE XMREAD(
     &           NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &           MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,
     &           ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &           XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &           HPRHOMA,NRHO,IRHO,NMASK,XRMAP,XRSYGP,XRSYIV)
C
C Parser routine for reading maps from specified file.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER NA, NB, NC, NAP, NBP, NCP, NAPP, NBPP, NCPP
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9)
      INTEGER ADEPTH, ARPNMX, ARPNN, ARPNX
      CHARACTER*(*) ARPN(4,*)
      INTEGER ARPNL(4,*)
      DOUBLE COMPLEX ARPNDB(4,*)
      INTEGER ARPNMLT(*)
      DOUBLE PRECISION MAPR
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER HPRHOMA, IRHO, NRHO, NMASK
      LOGICAL XRMAP
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
C local
      INTEGER ISTO, I
      INTEGER AMIN, AMAX, BMIN, BMAX, CMIN, CMAX
      INTEGER IUNIT, ILEN
      LOGICAL QFORM, COND, ERR
      INTEGER LNA, LNB, LNC
      DOUBLE PRECISION LCELL(6)
C pointers
      INTEGER SECTON, HPRMAP, HPIMAP
C parameters
C begin
C default values
      ERR=.FALSE.
      ISTO=1
      IFILE=' '
      IUNIT=0
      QFORM=.TRUE.
C
      CALL PUSEND('READ MAP>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('READ MAP>')
      CALL MISCOM('READ MAP>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-read-map')
C
C=====================================================================
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'TO') THEN
      CALL NEXTWD('TO=')
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD('TO=')
      IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(1X,2A)') 'TO=',XRHONAM(ISTO)
      ELSE
      COND=.FALSE.
      DO I=1,XRHONUM
      IF (WD(1:WDLEN).EQ.XRHONAM(I)) THEN
      ISTO=I
      COND=.TRUE.
      END IF
      END DO
      IF (.NOT.COND) THEN
      WRITE(6,'(3A)') ' %XMREAD-ERR: real space object ',
     & WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'XMREAD','object undeclared.')
      ERR=.TRUE.
      END IF
      END IF
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'INPU') THEN
      CALL NEXTFI('INPUt=',IFILE)
      ELSE IF (WD(1:4).EQ.'FORM') THEN
      CALL NEXTLO('FORMatted=',QFORM)
C
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' ---------------------------',
     & '-----MAP-parameters---------------------------------'
      ILEN=LEN(IFILE)
      CALL TRIMM(IFILE,ILEN)
      WRITE(6,'(3A,L1,2A)') ' | INPUt-file=',IFILE(1:ILEN),
     &    ' FORMatted=',QFORM,' TO=',XRHONAM(ISTO)
      WRITE(6,'(2A)') ' ---------------------------',
     & '----------------------------------------------------'
C=====================================================================
      ELSE
      CALL CHKEND('READ MAP>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (.NOT.ERR) THEN
C
      IF (IFILE.EQ.' ') THEN
      CALL WRNDIE(-5,'XMREAD',' no input file specified.')
      ELSE
C
C define the asymmetric unit if required
      IF (XRMAP) CALL XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   HPRHOMA,NRHO,NMASK,XRSYGP,XRSYIV)
      XRMAP=.FALSE.
C
      IF (HPRRHO(ISTO).EQ.0) THEN
      CALL XMAPAL(HPRRHO(ISTO),HPIRHO(ISTO),QHERM,NRHO,IRHO)
      END IF
      HPRMAP=HPRRHO(ISTO)
      HPIMAP=HPIRHO(ISTO)
C
C get the map pointers
C
C
C read sectioning information
      CALL RMAP1(IFILE,IUNIT,ERROR,EOF,
     &           LNA,AMIN,AMAX,LNB,BMIN,BMAX,
     &           LNC,CMIN,CMAX,LCELL,.FALSE.,QFORM)
C
C check to make sure sectioning and cell parameters are compatible
      IF (LNA.NE.NA.OR.LNB.NE.NB.OR.LNC.NE.NC) THEN
      WRITE(6,'(A,3I4,A,3I4)')
     & ' XREAD: map NA,NB,NC: ',LNA,LNB,LNC,'; expected NA,NB,NC: ',
     &  NA,NB,NC
      CALL WRNDIE(-5,'XMREAD',
     &    ' sectioning of map incompatible with resolution.')
      ELSE IF (ABS(LCELL(1)-XRCELL(1)).GT.RSMALL
     &     .OR.ABS(LCELL(2)-XRCELL(2)).GT.RSMALL
     &     .OR.ABS(LCELL(3)-XRCELL(3)).GT.RSMALL
     &     .OR.ABS(LCELL(4)-XRCELL(4)).GT.RSMALL
     &     .OR.ABS(LCELL(5)-XRCELL(5)).GT.RSMALL
     &     .OR.ABS(LCELL(6)-XRCELL(6)).GT.RSMALL) THEN
      WRITE(6,'(6(A,F9.4))')
     & ' Map: a=',LCELL(1),', b=',LCELL(2),', c=',LCELL(3),
     & ', alpha=',LCELL(4),', beta=',LCELL(5),', gamma=',LCELL(6)
      WRITE(6,'(6(A,F9.4))')
     & ' Expected: a=',XRCELL(1),', b=',XRCELL(2),', c=',XRCELL(3),
     & ', alpha=',XRCELL(4),', beta=',XRCELL(5),', gamma=',XRCELL(6)
      CALL WRNDIE(-5,'XMREAD',
     &    ' unit cell parameters incompatible.')
      ELSE
C
C write sectioning information
      WRITE(6,'(A,3I4,A,3I4,A,3I4,A)')
     & ' XREADX: extend NA=(',LNA,AMIN,AMAX,') NB=(',LNB,BMIN,BMAX,
     & ') NC=(',LNC,CMIN,CMAX,') '
C
C allocate temporary heap space for z-sections
      SECTON=ALLHP(IREAL8((AMAX-AMIN+1)*(BMAX-BMIN+1)))
C
      CALL XMREAD2(NA,NB,NC,
     &            MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &            HEAP(HPRMAP),HEAP(HPRHOMA),
     &            XRCELL,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &            AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,
     &            HEAP(SECTON),QFORM,IUNIT)
C
C deallocate temporary heap space for z-sections
      CALL FREHP(SECTON,IREAL8((AMAX-AMIN+1)*(BMAX-BMIN+1)))
C
      END IF
      END IF
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMREAD2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &            RRHO,RHOMASK,XRCELL,
     &            XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &            AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,SECTON,
     &            QFORM,IUNIT)
C
C Reads density map
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION XRCELL(*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4)
      INTEGER AMIN, AMAX, BMIN, BMAX, CMIN, CMAX
      DOUBLE PRECISION SECTON(AMIN:AMAX,BMIN:BMAX)
      LOGICAL QFORM
      INTEGER IUNIT
C local
      DOUBLE PRECISION SECS, SECS2, RAVE, RSIGMA
      CHARACTER*3 MODE
      LOGICAL COND, ERROR, EOF
      INTEGER A, B, C, AS, BS, CS, KSECT, SYM
C parameters
C
C begin
      EOF=.FALSE.
      ERROR=.FALSE.
C
C initialize timer
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
C read matrix mode
      IF (QFORM) THEN
      READ(IUNIT,'(A)',END=6,ERR=7) MODE
      ELSE
      READ(IUNIT,END=6,ERR=7) MODE
      END IF
      IF (MODE.NE.'ZYX') THEN
      CALL WRNDIE(-5,'XMREAD2','error in matrix mode')
      GOTO 7
      END IF
C
C read density matrix, c is slowest ("z-sections").
      DO C=CMIN,CMAX
C
C read next section
      IF (QFORM) THEN
      READ(IUNIT,'(I8)',END=6,ERR=7) KSECT
      READ(IUNIT,'(6E12.5)',END=6,ERR=7)
     &  ((SECTON(A,B),A=AMIN,AMAX),B=BMIN,BMAX)
      ELSE
      READ(IUNIT,END=6,ERR=7) KSECT
      READ(IUNIT,END=6,ERR=7)
     &  ((SECTON(A,B),A=AMIN,AMAX),B=BMIN,BMAX)
      END IF
C
C map the section into the asymmetric unit
C
      DO B=BMIN,BMAX
      DO A=AMIN,AMAX
C
      COND=(A.GE.MAASY.AND.A.LE.NAASY.AND.
     &      B.GE.MBASY.AND.B.LE.NBASY.AND.
     &      C.GE.MCASY.AND.C.LE.NCASY)
      IF (COND) COND=RHOMASK(A,B,C).GT.0
      IF (COND) THEN
C
C just copy
      RRHO(A,B,C)=SECTON(A,B)
      ELSE
C
C apply all symmetry operators to grid point until
C the symmetry-related grid point falls into the
C asymmetric unit
C
      SYM=1
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
      AS=MOD(XRSYMM(SYM,1,1)*A+XRSYMM(SYM,1,2)*B+XRSYMM(SYM,1,3)*C
     &       +(XRSYMM(SYM,1,4)*NA)/XRSYTH+10000*NA-MAASY,NA) +MAASY
      BS=MOD(XRSYMM(SYM,2,1)*A+XRSYMM(SYM,2,2)*B+XRSYMM(SYM,2,3)*C
     &       +(XRSYMM(SYM,2,4)*NB)/XRSYTH+10000*NB-MBASY,NB) +MBASY
      CS=MOD(XRSYMM(SYM,3,1)*A+XRSYMM(SYM,3,2)*B+XRSYMM(SYM,3,3)*C
     &       +(XRSYMM(SYM,3,4)*NC)/XRSYTH+10000*NC-MCASY,NC) +MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND) COND=RHOMASK(AS,BS,CS).GT.0
      IF (COND) RRHO(AS,BS,CS)=SECTON(A,B)
      SYM=SYM+1
      END DO
C
      IF (.NOT.COND) CALL WRNDIE(-5,'XMREAD2','Internal error no. 2')
C
      END IF
C
      END DO
      END DO
      END DO
C
C
C read average and scale factor
      IF (QFORM) THEN
      READ(IUNIT,'(I8)',END=6,ERR=7) KSECT
      READ(IUNIT,'(2(E12.4,1X))',END=6,ERR=7) RAVE, RSIGMA
      ELSE
      READ(IUNIT,END=6,ERR=7) KSECT
      READ(IUNIT,END=6,ERR=7) RAVE, RSIGMA
      END IF
C
C apply average and scale factor to map
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      RRHO(A,B,C)=RRHO(A,B,C)*RSIGMA+RAVE
      END DO
      END DO
      END DO
C
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4,/,A,F10.4)')
     &        ' XMREAD2: CPU-time: expansion=',SECS2-SECS
      END IF
C
      CALL VCLOSE(IUNIT,'KEEP',ERROR)
C
C-error-labels----
      GOTO 77
7     ERROR=.TRUE.
      GOTO 77
6     EOF=.TRUE.
77    CONTINUE
C-----------------
C
      IF (ERROR) WRITE(6,'(A)') ' %XMREAD2-ERR: error during read'
      IF (EOF) WRITE(6,'(A)') ' %XMREAD2-ERR: EOF during read'
C
      RETURN
      END
C======================================================================
      SUBROUTINE RMAP1(FILE,U,ERROR,EOF,
     &           NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX,
     &           CELL,QDEBUG,QFORM)
C
C Subroutine reads header information from a map written
C by XMAPWR from file FILE.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      CHARACTER*(*) FILE
      INTEGER U
      LOGICAL ERROR, EOF
      INTEGER NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX
      DOUBLE PRECISION CELL(6)
      LOGICAL QDEBUG, QFORM
C local
      INTEGER I
C begin
      ERROR=.FALSE.
      EOF=.FALSE.
      IF (QFORM) THEN
      CALL ASSFIL(FILE,U,'READ','FORMATTED',ERROR)
      ELSE
      CALL ASSFIL(FILE,U,'READ','UNFORMATTED',ERROR)
      END IF
      IF (.NOT.ERROR) THEN
C read title
      IF (QFORM) THEN
      CALL RDTITB(U,'FORMATTED')
      ELSE
      CALL RDTITB(U,'BINARY')
      END IF
C
C read sectioning information
      IF (QFORM) THEN
      READ(U,'(9I8)',END=6,ERR=7)
     &  NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX
      ELSE
      READ(U,END=6,ERR=7)
     &  NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX
      END IF
C
C read unit cell constants in angstroms and degrees
      IF (QFORM) THEN
      READ(U,'(6E12.5)',END=6,ERR=7) (CELL(I),I=1,6)
      ELSE
      READ(U,END=6,ERR=7) (CELL(I),I=1,6)
      END IF
      IF (QDEBUG) THEN
      WRITE(6,'(9I8)')
     &  NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX
      WRITE(6,'(6E12.5)') (CELL(I),I=1,6)
      END IF
      END IF
C
C-error-labels----
      GOTO 77
7     ERROR=.TRUE.
      GOTO 77
6     EOF=.TRUE.
77    CONTINUE
C-----------------
      IF (ERROR) WRITE(6,'(A)') ' %RMAP1-ERR: error during read'
      IF (EOF) WRITE(6,'(A)') ' %RMAP2-ERR: EOF during read'
      RETURN
      END
C======================================================================
      SUBROUTINE RMAP2(U,ERROR,EOF,
     &           AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,
     &           MAP,QFORM)
C
C Subroutine reads map written by XMAPWR from file FILE.
C Header information is read by subroutine RMAP1.
C
C This routine is currently only used by RSEARCH.S.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER U
      LOGICAL ERROR, EOF
      INTEGER AMIN,AMAX,BMIN,BMAX,CMIN,CMAX
      DOUBLE PRECISION MAP(AMIN:AMAX,BMIN:BMAX,CMIN:CMAX)
      LOGICAL QFORM
C local
      INTEGER A, B, C, KSECT
      CHARACTER*3 MODE
C begin
C
C read matrix mode
      IF (QFORM) THEN
      READ(U,'(A)',END=6,ERR=7) MODE
      ELSE
      READ(U,END=6,ERR=7) MODE
      END IF
      IF (MODE.NE.'ZYX') THEN
      CALL WRNDIE(-5,'RMAP','error in matrix mode')
      GOTO 7
      END IF
C
C read density matrix, c is slowest ("z-sections").
      DO C=CMIN,CMAX
C
C read next section
      IF (QFORM) THEN
      READ(U,'(I8)',END=6,ERR=7) KSECT
      READ(U,'(6E12.5)',END=6,ERR=7)
     &  ((MAP(A,B,C),A=AMIN,AMAX),B=BMIN,BMAX)
      ELSE
      READ(U,END=6,ERR=7) KSECT
      READ(U,END=6,ERR=7)
     &  ((MAP(A,B,C),A=AMIN,AMAX),B=BMIN,BMAX)
      END IF
      END DO
      CALL VCLOSE(U,'KEEP',ERROR)
C
C-error-labels----
      GOTO 77
7     ERROR=.TRUE.
      GOTO 77
6     EOF=.TRUE.
77    CONTINUE
C-----------------
C
      IF (ERROR) WRITE(6,'(A)') ' %RMAP2-ERR: error during read'
      IF (EOF) WRITE(6,'(A)') ' %RMAP2-ERR: EOF during read'
C
      RETURN
      END
