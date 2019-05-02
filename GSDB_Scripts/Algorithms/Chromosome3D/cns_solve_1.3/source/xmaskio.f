      SUBROUTINE XMSKRD(
     &           NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &           MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,
     &           ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &           XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &           HPRHOMA,NRHO,IRHO,NMASK,XRMAP,XRSYGP,XRSYIV)
C
C Parser routine for reading masks from specified file.
C A mask is a special case of a map file for CNS format masks
C the default format is the compressed O format
C
C Author: Paul Adams
C ==================
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
      CHARACTER*(4) MSKTYPE
      INTEGER ISTO, I
      INTEGER AMIN, AMAX, BMIN, BMAX, CMIN, CMAX
      INTEGER IUNIT, ILEN
      LOGICAL COND, ERR, REMAP
      INTEGER LNA, LNB, LNC
      DOUBLE PRECISION LCELL(6)
C pointers
      INTEGER SECTON, OMASK, HPRMAP, HPIMAP
C parameters
C begin
C default values
      ERR=.FALSE.
      ISTO=1
      IFILE=' '
      IUNIT=0
      MSKTYPE='OMAS'
      REMAP=.FALSE.
C
      CALL PUSEND('READ MASK>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('READ MASK>')
      CALL MISCOM('READ MASK>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-read-mask')
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
      WRITE(6,'(3A)') ' %XMSKRD-ERR: real space object ',
     & WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'XMSKRD','object undeclared.')
      ERR=.TRUE.
      END IF
      END IF
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'INPU') THEN
      CALL NEXTFI('INPUt=',IFILE)
      ELSE IF (WD(1:4).EQ.'TYPE') THEN
      CALL NEXTA4('TYPE=',MSKTYPE)
C
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' ---------------------------',
     & '-----MASK-parameters--------------------------------'
      ILEN=LEN(IFILE)
      CALL TRIMM(IFILE,ILEN)
      WRITE(6,'(3A,A,2A)') ' | INPUt-file=',IFILE(1:ILEN),
     &    ' TYPE=',MSKTYPE,' TO=',XRHONAM(ISTO)
      WRITE(6,'(2A)') ' ---------------------------',
     & '----------------------------------------------------'
C=====================================================================
      ELSE
      CALL CHKEND('READ MASK>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (.NOT.ERR) THEN
C
      IF (IFILE.EQ.' ') THEN
      CALL WRNDIE(-5,'XMSKRD','no input file specified.')
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
C read sectioning information
      IF (MSKTYPE.EQ.'OMAS') THEN
         CALL ROMSKHDR(IFILE,IUNIT,ERROR,EOF,
     &                 LNA,AMIN,AMAX,LNB,BMIN,BMAX,
     &                 LNC,CMIN,CMAX,LCELL)
      ELSE
         CALL RMAP1(IFILE,IUNIT,ERROR,EOF,
     &              LNA,AMIN,AMAX,LNB,BMIN,BMAX,
     &              LNC,CMIN,CMAX,LCELL,.FALSE.,.TRUE.)
      END IF
C
C check to make sure sectioning and cell parameters are compatible
      IF (MSKTYPE(1:3).EQ.'CNS') THEN
      IF (LNA.NE.NA.OR.LNB.NE.NB.OR.LNC.NE.NC) THEN
      WRITE(6,'(A,3I4,A,3I4)')
     & ' XREAD: mask NA,NB,NC: ',LNA,LNB,LNC,'; expected NA,NB,NC: ',
     &  NA,NB,NC
      CALL WRNDIE(-5,'XMSKRD',
     &    ' sectioning of mask incompatible with resolution.')
      END IF
      END IF
C
      IF (ABS(LCELL(1)-XRCELL(1)).GT.RSMALL
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
     & ' XMSKRD: extent NA=(',LNA,AMIN,AMAX,') NB=(',LNB,BMIN,BMAX,
     & ') NC=(',LNC,CMIN,CMAX,') '
C
C if omask must be read at one go, else do in sections
      IF (MSKTYPE.EQ.'OMAS') THEN
      IF (LNA.NE.NA.OR.LNB.NE.NB.OR.LNC.NE.NC) THEN
         REMAP=.TRUE.
      END IF
C
C allocate temporary heap space for z-sections
      OMASK=ALLHP(IREAL4(
     &           (AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
C
      CALL XMSKORD(NA,NB,NC,LNA,LNB,LNC,
     &             MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &             HEAP(HPRMAP),HEAP(HPRHOMA),
     &             XRCELL,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &             AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,
     &             HEAP(OMASK),IUNIT,REMAP)
C
C deallocate temporary heap space for z-sections
      CALL FREHP(OMASK,IREAL4(
     &           (AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
C
      ELSE
C
C allocate temporary heap space for z-sections
      SECTON=ALLHP(IREAL4((AMAX-AMIN+1)*(BMAX-BMIN+1)))
C
      CALL XMSKRD2(NA,NB,NC,
     &             MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &             HEAP(HPRMAP),HEAP(HPRHOMA),
     &             XRCELL,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &             AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,
     &             HEAP(SECTON),.TRUE.,IUNIT)
C
C deallocate temporary heap space for z-sections
      CALL FREHP(SECTON,IREAL4((AMAX-AMIN+1)*(BMAX-BMIN+1)))
C
      END IF
C
      END IF
      END IF
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMSKRD2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &            RRHO,RHOMASK,XRCELL,
     &            XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &            AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,SECTON,
     &            QFORM,IUNIT)
C
C Reads CNS format mask
C
C Author: Paul Adams
C ==================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'xmask.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION XRCELL(*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4)
      INTEGER AMIN, AMAX, BMIN, BMAX, CMIN, CMAX
      REAL SECTON(AMIN:AMAX,BMIN:BMAX)
      LOGICAL QFORM
      INTEGER IUNIT
C local
      DOUBLE PRECISION SECS, SECS2, RAVE, RSIGMA
      CHARACTER*3 MODE
      LOGICAL COND, ERROR, EOF
      INTEGER TMASKIN, TMASKOUT
      INTEGER A, B, C, AS, BS, CS, KSECT, SYM
      INTEGER AA, BB, CC, SYMMOP, LTX, LTY, LTZ
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
      CALL WRNDIE(-5,'XMSKRD2','error in matrix mode')
      GOTO 7
      END IF
C
C zero overall statistics accumulators
      TMASKIN=0
      TMASKOUT=0
C
C intialize mask
      DO A=MAASY,NAASY
      DO B=MBASY,NBASY
      DO C=MCASY,NCASY
         RRHO(A,B,C)=1
      END DO
      END DO
      END DO
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
C just copy - do not overwrite a grid point that has been
C set as a mask point (this is to deal with maps that extend over
C more than the unit cell and may have symmetry related points
C marked as outside of the mask)
C enter the identity operator and zero translations into the grid point
      IF (SECTON(A,B).LE.0) THEN
         CALL MSKA2M(SECTON(A,B),SYMMOP,LTX,LTY,LTZ,XRMSYM)
         CALL MSKM2A(SYMMOP,LTX,LTY,LTZ,RRHO(A,B,C),XRMSYM)
         TMASKIN=TMASKIN+1
      ELSE
         IF (RRHO(A,B,C).GT.0) THEN
            RRHO(A,B,C)=SECTON(A,B)
            IF (SECTON(A,B).GT.0) TMASKOUT=TMASKOUT+1
         END IF
      END IF
      ELSE
C
C apply all symmetry operators to grid point until
C the symmetry-related grid point falls into the
C asymmetric unit
C
      SYM=1
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
C
      AA=XRSYMM(SYM,1,1)*A+XRSYMM(SYM,1,2)*B+XRSYMM(SYM,1,3)*C
     &       +(XRSYMM(SYM,1,4)*NA)/XRSYTH
      BB=XRSYMM(SYM,2,1)*A+XRSYMM(SYM,2,2)*B+XRSYMM(SYM,2,3)*C
     &       +(XRSYMM(SYM,2,4)*NB)/XRSYTH
      CC=XRSYMM(SYM,3,1)*A+XRSYMM(SYM,3,2)*B+XRSYMM(SYM,3,3)*C
     &       +(XRSYMM(SYM,3,4)*NC)/XRSYTH
C
      AS=MOD(AA+10000*NA-MAASY,NA) +MAASY
      BS=MOD(BB+10000*NB-MBASY,NB) +MBASY
      CS=MOD(CC+10000*NC-MCASY,NC) +MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND)  COND=RHOMASK(AS,BS,CS).GT.0
      SYM=SYM+1
      END DO
C
      IF (COND) THEN
      SYM=SYM-1
C
C enter the identity operator and zero translations into the grid point
C if this is negative then it must be an averaging map, store the
C remapping information needed to put the point back to its original
C place (should only occur for a molecular mask, not ASU)
      IF (SECTON(A,B).LE.0) THEN
      LTX=(AS-AA)/NA
      LTY=(BS-BB)/NB
      LTZ=(CS-CC)/NC
      CALL MSKM2A(SYM,LTX,LTY,LTZ,RRHO(AS,BS,CS),XRMSYM)
      TMASKIN=TMASKIN+1
C
C do not overwrite a grid point that has been
C set as a mask point (this is to deal with maps that extend over
C more than the unit cell and may have symmetry related points
C marked as outside of the mask)
      ELSE IF (RRHO(AS,BS,CS).GT.0) THEN
         RRHO(AS,BS,CS)=SECTON(A,B)
         IF (SECTON(A,B).GT.0) TMASKOUT=TMASKOUT+1
      END IF
C
      END IF
C
      IF (.NOT.COND)
     &    CALL WRNDIE(-5,'XMSKRD2','point could not be remapped')
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
C print info about the whole map
      WRITE(6,'(A,I8,A,I8)')
     & ' XMSKRD: total # points inside mask',TMASKIN,
     & ' total # points outside mask',TMASKOUT
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4,/,A,F10.4)')
     &        ' XMSKRD2: CPU-time: expansion=',SECS2-SECS
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
      IF (ERROR) WRITE(6,'(A)') ' %XMSKRD2-ERR: error during read'
      IF (EOF) WRITE(6,'(A)') ' %XMSKRD2-ERR: EOF during read'
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMSKORD(NA,NB,NC,LNA,LNB,LNC,
     &            MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &            RRHO,RHOMASK,XRCELL,
     &            XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &            LAMIN,LAMAX,LBMIN,LBMAX,LCMIN,LCMAX,OMASK,
     &            IUNIT,REMAP)
C
C Reads mask in compressed O format
C front end to deal with interpolation of the O mask
C onto to the current grid, if necessary
C
C Author: Paul Adams 14-2-97
C ==========================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'xmask.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER LNA, LNB, LNC
      INTEGER LAMIN, LAMAX, LBMIN, LBMAX, LCMIN, LCMAX
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION XRCELL(*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4)
      REAL OMASK(LAMIN:LAMAX,LBMIN:LBMAX,LCMIN:LCMAX)
      INTEGER IUNIT
      LOGICAL REMAP
C local
      INTEGER LENGTH, AMIN, AMAX, BMIN, BMAX, CMIN, CMAX
      INTEGER LAMIN2, LAMAX2, LBMIN2, LBMAX2, LCMIN2, LCMAX2
C pointers
      INTEGER CGMSK,OGMSK
C begin
C
C read the mask into the complete mask array
      LENGTH=(LAMAX-LAMIN+1)*(LBMAX-LBMIN+1)*(LCMAX-LCMIN+1)
      CALL MSKOCOMPR(OMASK,LENGTH,IUNIT)
C
      IF (REMAP) THEN
         WRITE (6,'(A)')' XMSKORD: remapping mask to current grid'
C
C calculate parameters for the current O mask grid
C
C calculate dimensions of a grid in current grid units
C that covers the current O mask
         CALL GRIDSZ(LNA,LNB,LNC,LAMIN,LAMAX,LBMIN,
     &               LBMAX,LCMIN,LCMAX,NA,NB,NC,
     &               AMIN,AMAX,BMIN,BMAX,CMIN,CMAX)
         WRITE(6,'(A,3I4,A,3I4,A,3I4,A)')
     &        ' XMSKORD: from extent NA=(',LNA,LAMIN,LAMAX,
     &        ') NB=(',LNB,LBMIN,LBMAX,
     &        ') NC=(',LNC,LCMIN,LCMAX,') '
         WRITE(6,'(A,3I4,A,3I4,A,3I4,A)')
     &        ' XMSKORD:   to extent NA=(',NA,AMIN,AMAX,
     &        ') NB=(',NB,BMIN,BMAX,
     &        ') NC=(',NC,CMIN,CMAX,') '
C
C calculate dimensions of a grid in O mask grid units
C that covers the grid calculated above
         CALL GRIDSZ(NA,NB,NC,AMIN,AMAX,BMIN,
     &               BMAX,CMIN,CMAX,LNA,LNB,LNC,
     &               LAMIN2,LAMAX2,LBMIN2,LBMAX2,LCMIN2,LCMAX2)
C
C expand grid by 1 in all directions
         LAMIN2=LAMIN2-1
         LAMAX2=LAMAX2+1
         LBMIN2=LBMIN2-1
         LBMAX2=LBMAX2+1
         LCMIN2=LCMIN2-1
         LCMAX2=LCMAX2+1
C
C allocate space for mask
         CGMSK=ALLHP(IREAL4(
     &   (AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
C
C allocate space for enlarged O mask
         OGMSK=ALLHP(IREAL4(
     &   (LAMAX2-LAMIN2+1)*(LBMAX2-LBMIN2+1)*(LCMAX2-LCMIN2+1)))
C
C interpolate the current O mask onto the current grid
C nb. crystal symmetry is not applied in this step
         CALL REGRID(OMASK,LNA,LNB,LNC,
     &               LAMIN,LAMAX,LBMIN,LBMAX,LCMIN,LCMAX,
     &               HEAP(OGMSK),
     &               LAMIN2,LAMAX2,LBMIN2,LBMAX2,LCMIN2,LCMAX2,
     &               HEAP(CGMSK),NA,NB,NC,
     &               AMIN,AMAX,BMIN,BMAX,CMIN,CMAX)
C
C remap mask into the asymmetric unit
         CALL XMSKORD2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &        RRHO,RHOMASK,XRCELL,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &        AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,HEAP(CGMSK),
     &        IUNIT)
C
C free up allocated space
         CALL FREHP(OGMSK,IREAL4(
     &   (LAMAX2-LAMIN2+1)*(LBMAX2-LBMIN2+1)*(LCMAX2-LCMIN2+1)))
         CALL FREHP(CGMSK,IREAL4(
     &   (AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
C
      ELSE
C
C remap mask into the asymmetric unit
         CALL XMSKORD2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &        RRHO,RHOMASK,XRCELL,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &        LAMIN,LAMAX,LBMIN,LBMAX,LCMIN,LCMAX,OMASK,
     &        IUNIT)
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMSKORD2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &            RRHO,RHOMASK,XRCELL,
     &            XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &            AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,OMASK,
     &            IUNIT)
C
C Reads mask in compressed O format
C
C Author: Paul Adams
C ==================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'xmask.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION XRCELL(*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4)
      INTEGER AMIN, AMAX, BMIN, BMAX, CMIN, CMAX
      REAL OMASK(AMIN:AMAX,BMIN:BMAX,CMIN:CMAX)
      INTEGER IUNIT
C local
      DOUBLE PRECISION SECS, SECS2, OLPPER
      LOGICAL COND, ERROR, EOF
      INTEGER TMASKIN, TMASKOUT, OVERLAP
      INTEGER A, B, C, AS, BS, CS, SYM
      INTEGER AA, BB, CC, SYMMOP, LTX, LTY, LTZ
C parameters
C
C begin
      EOF=.FALSE.
      ERROR=.FALSE.
C
C initialize timer
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
C zero overall statistics accumulators
      OVERLAP=0
      TMASKIN=0
      TMASKOUT=0
C
C intialize mask
      DO A=MAASY,NAASY
      DO B=MBASY,NBASY
      DO C=MCASY,NCASY
         RRHO(A,B,C)=9999
      END DO
      END DO
      END DO
C
C read density matrix, c is slowest ("z-sections").
      DO C=CMIN,CMAX
C
C map the section into the asymmetric unit
      DO B=BMIN,BMAX
      DO A=AMIN,AMAX
C
      COND=(A.GE.MAASY.AND.A.LE.NAASY.AND.
     &      B.GE.MBASY.AND.B.LE.NBASY.AND.
     &      C.GE.MCASY.AND.C.LE.NCASY)
      IF (COND) COND=RHOMASK(A,B,C).GT.0
      IF (COND) THEN
C
C just copy - do not overwrite a grid point that has been
C set as a mask point (this is to deal with maps that extend over
C more than the unit cell and may have symmetry related points
C marked as outside of the mask)
C enter the identity operator and zero translations into the grid point
      IF (OMASK(A,B,C).EQ.-1) THEN
         CALL MSKA2M(OMASK(A,B,C),SYMMOP,LTX,LTY,LTZ,XRMSYM)
         IF (RRHO(A,B,C).LE.0) OVERLAP=OVERLAP+1
         CALL MSKM2A(SYMMOP,LTX,LTY,LTZ,RRHO(A,B,C),XRMSYM)
         TMASKIN=TMASKIN+1
      ELSE
         IF (RRHO(A,B,C).GT.0) THEN
            RRHO(A,B,C)=OMASK(A,B,C)
            IF (OMASK(A,B,C).EQ.1) TMASKOUT=TMASKOUT+1
         END IF
      END IF
      ELSE
C
C apply all symmetry operators to grid point until
C the symmetry-related grid point falls into the
C asymmetric unit
C
      SYM=1
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
C
      AA=XRSYMM(SYM,1,1)*A+XRSYMM(SYM,1,2)*B+XRSYMM(SYM,1,3)*C
     &       +(XRSYMM(SYM,1,4)*NA)/XRSYTH
      BB=XRSYMM(SYM,2,1)*A+XRSYMM(SYM,2,2)*B+XRSYMM(SYM,2,3)*C
     &       +(XRSYMM(SYM,2,4)*NB)/XRSYTH
      CC=XRSYMM(SYM,3,1)*A+XRSYMM(SYM,3,2)*B+XRSYMM(SYM,3,3)*C
     &       +(XRSYMM(SYM,3,4)*NC)/XRSYTH
C
      AS=MOD(AA+10000*NA-MAASY,NA) +MAASY
      BS=MOD(BB+10000*NB-MBASY,NB) +MBASY
      CS=MOD(CC+10000*NC-MCASY,NC) +MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND)  COND=RHOMASK(AS,BS,CS).GT.0
      SYM=SYM+1
      END DO
C
      IF (COND) THEN
      SYM=SYM-1
C
C enter the identity operator and zero translations into the grid point
C if this is negative then it must be an averaging map, store the
C remapping information needed to put the point back to its original
C place (should only occur for a molecular mask, not ASU)
      IF (OMASK(A,B,C).EQ.-1) THEN
      LTX=(AS-AA)/NA
      LTY=(BS-BB)/NB
      LTZ=(CS-CC)/NC
      IF (RRHO(AS,BS,CS).LE.0) OVERLAP=OVERLAP+1
      CALL MSKM2A(SYM,LTX,LTY,LTZ,RRHO(AS,BS,CS),XRMSYM)
      TMASKIN=TMASKIN+1
C
C do not overwrite a grid point that has been
C set as a mask point (this is to deal with maps that extend over
C more than the unit cell and may have symmetry related points
C marked as outside of the mask)
      ELSE
         IF (RRHO(AS,BS,CS).GT.0) THEN
            RRHO(AS,BS,CS)=OMASK(A,B,C)
            IF (OMASK(A,B,C).EQ.1) TMASKOUT=TMASKOUT+1
         END IF
      END IF
C
      END IF
C
      IF (.NOT.COND)
     &    CALL WRNDIE(-5,'XMSKORD','point could not be remapped')
C
      END IF
C
      END DO
      END DO
      END DO
C
C set any untouched mask points to solvent
      DO A=MAASY,NAASY
      DO B=MBASY,NBASY
      DO C=MCASY,NCASY
         IF (RRHO(A,B,C).EQ.9999) RRHO(A,B,C)=1
      END DO
      END DO
      END DO
C
C print info about the whole map
      WRITE(6,'(A,I8,A,I8)')
     & ' XMSKORD: total # points inside mask',TMASKIN,
     & ' total # points outside mask',TMASKOUT
C
C print warnings or stop if too much overlap
      OLPPER=(OVERLAP*100.0D0)/TMASKIN
      WRITE (6,'(A,F6.1,A)')
     &   ' XMSKORD-WRN: the mask overlaps itself by ',OLPPER,'%'
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4,/,A,F10.4)')
     &        ' XMSKORD: CPU-time: expansion=',SECS2-SECS
      END IF
C
      CALL VCLOSE(IUNIT,'KEEP',ERROR)
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMSKWR(
     &           NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &           MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,QHERM,XRCELL,MAPR,
     &           XRHONUM,XRHONAM,HPRRHO,
     &           HPRHOMA,NRHO,NMASK,XRMAP,
     &           XRTR,XRINTR,XRVOL,
     &           XRMREF,XRNREF,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,
     &           HPTYPE,XRRED,XRSYGP,XRSYIV)
C
C Parser routine for writing masks to specified file.
C A mask is a special case of a map
C
C Author: Paul Adams
C ==================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'xmask.inc'
      INTEGER NA, NB, NC, NAP, NBP, NCP, NAPP, NBPP, NCPP
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9)
      DOUBLE PRECISION MAPR
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*)
      INTEGER HPRHOMA, NRHO, NMASK
      LOGICAL XRMAP
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3), XRVOL
      INTEGER XRMREF, XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE
      LOGICAL XRRED
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
C local
      CHARACTER*4 MSKTYPE, EXTEND
      DOUBLE PRECISION CUSHON
      DOUBLE PRECISION XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
      DOUBLE PRECISION TMPMAT(MXMAVG,3,4), CENTER(3)
      INTEGER AMIN, AMAX, BMIN, BMAX, CMIN, CMAX, IFROM, I, J
      INTEGER OUNIT, OLEN
      LOGICAL COND, ERR
      INTEGER LNA, LNB, LNC
      INTEGER LMAASY, LMBASY, LMCASY, LNAASY, LNBASY, LNCASY
C pointers
      INTEGER SECTON, HPRMAP, HPMASK, GRIDRT
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, FOUR=4.0D0)
      REAL SZERO
      PARAMETER (SZERO=0.0)
C begin
C default values
      ERR=.FALSE.
      IFROM=1
      OFILE=' '
      OUNIT=0
      MSKTYPE='OMAS'
      EXTEND='MOLE'
      CUSHON=TWO
      XMIN=ZERO
      XMAX=ONE
      YMIN=ZERO
      YMAX=ONE
      ZMIN=ZERO
      ZMAX=ONE
C
      CALL PUSEND('WRITe MASK>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('WRITe MASK>')
      CALL MISCOM('WRITe MASK>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-write-mask')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FROM') THEN
      CALL NEXTWD('FROM=')
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD('FROM=')
      IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(1X,2A)') 'FROM=',XRHONAM(IFROM)
      ELSE
C
      COND=.FALSE.
      DO I=1,XRHONUM
      IF (WD(1:WDLEN).EQ.XRHONAM(I)) THEN
      COND=.TRUE.
      IFROM=I
      END IF
      END DO
C
      IF (.NOT.COND) THEN
      WRITE(6,'(3A)') ' %XMSKWR-ERR: real space object ',
     & WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'XMSKWR','object undeclared.')
      ERR=.TRUE.
      END IF
      END IF
C
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
      ELSE IF (WD(1:4).EQ.'TYPE') THEN
      CALL NEXTA4('TYPE=',MSKTYPE)
      COND= MSKTYPE(1:3).EQ.'CNS'.OR.
     &      MSKTYPE(1:4).EQ.'OMAS'
      IF (.NOT.COND) CALL DSPERR('XMSKWR','unknown value for TYPE')
      ELSE IF (WD(1:4).EQ.'CUSH') THEN
      CALL NEXTF('CUSHion=',CUSHON)
      ELSE IF (WD(1:4).EQ.'EXTE') THEN
      CALL NEXTA4('EXTEnt=',EXTEND)
      COND= EXTEND.EQ.'UNIT'.OR.
     &      EXTEND.EQ.'MOLE'.OR.
     &      EXTEND.EQ.'BOX '.OR.
     &      EXTEND.EQ.'FRAC'.OR.
     &      EXTEND.EQ.'ASYM'
      IF (.NOT.COND) CALL DSPERR('XMSKWR','unknown value for EXTEnt')
      ELSE IF (WD(1:4).EQ.'XMIN') THEN
         CALL NEXTF('XMIN=',XMIN)
      ELSE IF (WD(1:4).EQ.'XMAX') THEN
         CALL NEXTF('XMAX=',XMAX)
      ELSE IF (WD(1:4).EQ.'YMIN') THEN
         CALL NEXTF('YMIN=',YMIN)
      ELSE IF (WD(1:4).EQ.'YMAX') THEN
         CALL NEXTF('YMAX=',YMAX)
      ELSE IF (WD(1:4).EQ.'ZMIN') THEN
         CALL NEXTF('ZMIN=',ZMIN)
      ELSE IF (WD(1:4).EQ.'ZMAX') THEN
         CALL NEXTF('ZMAX=',ZMAX)
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' ---------------------------',
     & '-----MASK-parameters--------------------------------'
      OLEN=LEN(OFILE)
      CALL TRIMM(OFILE,OLEN)
      WRITE(6,'(6A)') ' | OUTPut-file=',OFILE(1:OLEN),
     &    ' TYPE=',MSKTYPE,
     &    ' | FROM=',XRHONAM(IFROM)
      IF (MSKTYPE.EQ.'OMAS') THEN
      WRITE(6,'(2A)')' | EXTEnt=',EXTEND
      IF (EXTEND.EQ.'BOX '.OR.EXTEND.EQ.'FRAC') THEN
      WRITE(6,'(6(A,F9.4))')
     & ' | XMIN=',XMIN,' XMAX=',XMAX,' YMIN=',YMIN,' YMAX=',YMAX,
     & ' ZMIN=',ZMIN,' ZMAX=',ZMAX
      ELSE IF (EXTEND.EQ.'MOLE') THEN
      WRITE(6,'(A,F6.2)') ' | CUSHion=',CUSHON
      END IF
      ELSE IF (MSKTYPE(1:3).EQ.'CNS') THEN
      EXTEND='ASYM'
      WRITE(6,'(2A)')' | EXTEnt=',EXTEND
      END IF
      WRITE(6,'(2A)') ' ---------------------------',
     & '----------------------------------------------------'
C=====================================================================
      ELSE
      CALL CHKEND('WRITe MASK>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (.NOT.ERR) THEN
C
      IF (OFILE.NE.' ') THEN
      CALL ASSFIL(OFILE,OUNIT,'WRITE','FORMATTED',ERROR)
      END IF
C
C get the map pointers
      IF (HPRRHO(IFROM).EQ.0) THEN
      WRITE(6,'(3A)') ' %XMSKWR-ERR: real space object ',
     & XRHONAM(IFROM),' undefined.'
      CALL WRNDIE(-5,'XMSKWR','object undefined.')
      ERR=.TRUE.
      ELSE
      HPRMAP=HPRRHO(IFROM)
      END IF
C
      LNA=NA
      LNB=NB
      LNC=NC
      LMAASY=MAASY
      LMBASY=MBASY
      LMCASY=MCASY
      LNAASY=NAASY
      LNBASY=NBASY
      LNCASY=NCASY
      HPMASK=HPRHOMA
C
C error checking
      IF (OUNIT.EQ.0) THEN
      CALL WRNDIE(-5,'XMSKWR','no output file specified.')
      ELSE IF (.NOT.ERR) THEN
C
C determine extent of map
      IF (MSKTYPE(1:3).EQ.'CNS') THEN
C
C use asymmetric unit
      AMIN=LMAASY
      AMAX=LNAASY
      BMIN=LMBASY
      BMAX=LNBASY
      CMIN=LMCASY
      CMAX=LNCASY
C
C write sectioning information
      WRITE(6,'(A,3I4,A,3I4,A,3I4,A)')
     & ' XMSKWR: extent NA=(',LNA,AMIN,AMAX,') NB=(',LNB,BMIN,BMAX,
     & ') NC=(',LNC,CMIN,CMAX,') '
C
C allocate temporary heap space for z-sections
      SECTON=ALLHP(IREAL4((AMAX-AMIN+1)*(BMAX-BMIN+1)))
C
      CALL XMSKWR2(LNA,LNB,LNC,
     &             LMAASY,LMBASY,LMCASY,LNAASY,LNBASY,LNCASY,
     &             QHERM,HEAP(HPRMAP),HEAP(HPMASK),
     &             XRCELL,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &             AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,
     &             HEAP(SECTON),
     &             OUNIT,MSKTYPE)
      IF (OUNIT.GT.0) CALL VCLOSE(OUNIT,'KEEP',ERROR)
C
C deallocate temporary heap space for z-sections
      CALL FREHP(SECTON,IREAL4((AMAX-AMIN+1)*(BMAX-BMIN+1)))
C
      ELSE
      IF (EXTEND.EQ.'MOLE') THEN
C
      DO I=1,3
      DO J=1,3
      IF (I.EQ.J) THEN
      TMPMAT(1,I,J)=1.0D0
      ELSE
      TMPMAT(1,I,J)=0.0D0
      END IF
      END DO
      TMPMAT(1,I,4)=0.0D0
      CENTER(I)=0.0D0
      END DO
C
      GRIDRT=ALLHP(IREAL8(XRNSYM*3*4))
C
C find the limits of the mask
      CALL MSKDIM(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,XRTR,XRINTR,
     &     HEAP(HPMASK),HEAP(HPRMAP),HEAP(GRIDRT),
     &     TMPMAT,CENTER,CUSHON,AMIN,AMAX,BMIN,BMAX,CMIN,CMAX)
C
      CALL FREHP(GRIDRT,IREAL8(XRNSYM*3*4))
C
      ELSE IF (EXTEND.EQ.'ASYM') THEN
C
C use asymmetric unit
      AMIN=LMAASY
      AMAX=LNAASY
      BMIN=LMBASY
      BMAX=LNBASY
      CMIN=LMCASY
      CMAX=LNCASY
C
      ELSE IF (EXTEND.EQ.'UNIT'.OR.
     &         EXTEND.EQ.'BOX '.OR.
     &         EXTEND.EQ.'FRAC') THEN
C
      CALL XMSKMX(EXTEND,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,XRTR,XRCELL)
C
      AMIN=NINT(XMIN*LNA)
      AMAX=NINT(XMAX*LNA)
      BMIN=NINT(YMIN*LNB)
      BMAX=NINT(YMAX*LNB)
      CMIN=NINT(ZMIN*LNC)
      CMAX=NINT(ZMAX*LNC)
C
      ELSE
      CALL WRNDIE(-5,'XMSKWR','EXTEnt mode incorrectly defined.')
C
      END IF
C
C write sectioning information
      WRITE(6,'(A,3I4,A,3I4,A,3I4,A)')
     & ' XMSKWR: extent NA=(',LNA,AMIN,AMAX,') NB=(',LNB,BMIN,BMAX,
     & ') NC=(',LNC,CMIN,CMAX,') '
C
C allocate temporary heap space for mask
      SECTON=ALLHP(IREAL4(
     &          (AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
C
      CALL XMSKOWR(LNA,LNB,LNC,LMAASY,LMBASY,LMCASY,LNAASY,LNBASY,
     &             LNCASY,QHERM,HEAP(HPRMAP),HEAP(HPMASK),XRCELL,
     &             XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,AMIN,AMAX,
     &             BMIN,BMAX,CMIN,CMAX,HEAP(SECTON),OUNIT,EXTEND)
      IF (OUNIT.GT.0) CALL VCLOSE(OUNIT,'KEEP',ERROR)
C
C deallocate temporary heap space for mask
      CALL FREHP(SECTON,IREAL4(
     &          (AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
C
      END IF
C
      END IF
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE  XMSKWR2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &            QHERM,RRHO,RHOMASK,XRCELL,
     &            XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &            AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,SECTON,
     &            OUNIT,MSKTYPE)
C
C Write mask
C
C Author: Paul Adams
C ==================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'xmask.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      LOGICAL QHERM
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRCELL(*)
      INTEGER AMIN, AMAX, BMIN, BMAX, CMIN, CMAX
      REAL SECTON(AMAX-AMIN+1,BMAX-BMIN+1)
      INTEGER OUNIT
      CHARACTER*(*) MSKTYPE
C local
      INTEGER A, B, C
      INTEGER KSECT, ISECT, JSECT, I, J
      INTEGER MASKIN, MASKOUT, TMASKIN, TMASKOUT
      DOUBLE PRECISION SECS, SECS2
C parameters
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
C begin
C
C initialize timers
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
C write map headers
      IF (OUNIT.NE.0) THEN
      CALL WRTITL(OUNIT,2)
      WRITE(OUNIT,'(9I8)')
     &  NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX
      WRITE(OUNIT,'(6E12.5)') (XRCELL(I),I=1,6)
      WRITE(OUNIT,'(A)') 'ZYX'
      END IF
C
C zero overall statistics accumulators
      TMASKIN=0
      TMASKOUT=0
C
C write density matrix, c is slowest ("z-sections").
      KSECT=-1
      DO C=CMIN,CMAX
      KSECT=KSECT+1
C
C zero section statistics accumulators
      MASKIN=0
      MASKOUT=0
C
C zero points in section first
      JSECT=0
      DO B=BMIN,BMAX
      JSECT=JSECT+1
      ISECT=0
      DO A=AMIN,AMAX
      ISECT=ISECT+1
      SECTON(ISECT,JSECT)=1
      END DO
      END DO
C
C we're just writing the asymmetric unit
      JSECT=0
      DO B=BMIN,BMAX
      JSECT=JSECT+1
      ISECT=0
      DO A=AMIN,AMAX
      ISECT=ISECT+1
C
C only check if point is in ASU - even if point is outside RHOMASK
C it will be written as is (assumed to be zero)
      IF (A.GE.MAASY.AND.B.GE.MBASY.AND.C.GE.MCASY.AND.
     &    A.LE.NAASY.AND.B.LE.NBASY.AND.C.LE.NCASY) THEN
C
C keep remapping information
      SECTON(ISECT,JSECT)=RRHO(A,B,C)
      IF (RHOMASK(A,B,C).GT.0) THEN
         IF (RRHO(A,B,C).LE.0) MASKIN=MASKIN+1
         IF (RRHO(A,B,C).GT.0) MASKOUT=MASKOUT+1
      END IF
      ELSE
      CALL WRNDIE(-5,'XMSKWR2','Mask point outside ASU in ASYM mode')
      END IF
      END DO
      END DO
C
      IF (OUNIT.GT.0) THEN
C
C print info about the section
      WRITE(6,'(A,I4,A,I8,A,I8)')
     & ' XMSKWR: mask section #',KSECT,' points inside mask',MASKIN,
     & ' points outside mask',MASKOUT
C
C write the section
      WRITE(OUNIT,'(I8)') KSECT
      WRITE(OUNIT,'(1P6E12.5)')
     &  ((SECTON(I,J),I=1,ISECT),J=1,JSECT)
C
      END IF
C
C update overall statistics counters
      TMASKIN=TMASKIN+MASKIN
      TMASKOUT=TMASKOUT+MASKOUT
C
      END DO
C
C print info about the whole map
      WRITE(6,'(A,I8,A,I8)')
     & ' XMSKWR: total # points inside mask',TMASKIN,
     & ' total # points outside mask',TMASKOUT
C
C terminate output file with negative integer number
C also write map scaling parameters (always zero and one for masks).
      IF (OUNIT.GT.0) THEN
      WRITE(OUNIT,'(I8)') -9999
      WRITE(OUNIT,'(2(E12.4,1X))') ZERO,ONE
      END IF
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4,/,A,F10.4)')
     &        ' XMSKWR: CPU-time: expansion=',SECS2-SECS
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMSKMX(EXTEND,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,
     &                  XRTR,XRCELL)
C
C determine max and min values, convert into fractional coordinates
C
C Author: Axel T. Brunger and Paul Adams
C ======================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'funct.inc'
      CHARACTER*4 EXTEND
      DOUBLE PRECISION XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
      DOUBLE PRECISION XRTR(3,3), XRCELL(*)
C local
      INTEGER IAT
      DOUBLE PRECISION TEMPX, TEMPY, TEMPZ
      DOUBLE PRECISION XL(8), YL(8), ZL(8)
C parameters
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
      IF (EXTEND.EQ.'UNIT') THEN
      XMIN=ZERO
      XMAX=ONE
      YMIN=ZERO
      YMAX=ONE
      ZMIN=ZERO
      ZMAX=ONE
C
      ELSE IF (EXTEND.EQ.'FRAC') THEN
C
      CONTINUE
C
      ELSE IF (EXTEND.EQ.'BOX ') THEN
C
C set up fake coordinate list of cube edges
      XL(1)=XMIN
      YL(1)=YMIN
      ZL(1)=ZMIN
      XL(2)=XMAX
      YL(2)=YMIN
      ZL(2)=ZMIN
      XL(3)=XMIN
      YL(3)=YMAX
      ZL(3)=ZMIN
      XL(4)=XMAX
      YL(4)=YMAX
      ZL(4)=ZMIN
      XL(5)=XMIN
      YL(5)=YMIN
      ZL(5)=ZMAX
      XL(6)=XMAX
      YL(6)=YMIN
      ZL(6)=ZMAX
      XL(7)=XMIN
      YL(7)=YMAX
      ZL(7)=ZMAX
      XL(8)=XMAX
      YL(8)=YMAX
      ZL(8)=ZMAX
C
C
C fractionalize cube edge coordinates
      DO IAT=1,8
      TEMPX=XL(IAT)
      TEMPY=YL(IAT)
      TEMPZ=ZL(IAT)
      XL(IAT)=XRTR(1,1)*TEMPX +XRTR(1,2)*TEMPY +XRTR(1,3)*TEMPZ
      YL(IAT)=XRTR(2,1)*TEMPX +XRTR(2,2)*TEMPY +XRTR(2,3)*TEMPZ
      ZL(IAT)=XRTR(3,1)*TEMPX +XRTR(3,2)*TEMPY +XRTR(3,3)*TEMPZ
      END DO
C
C determine min and max
      XMIN=XL(1)
      XMAX=XL(1)
      YMIN=YL(1)
      YMAX=YL(1)
      ZMIN=ZL(1)
      ZMAX=ZL(1)
      DO IAT=2,8
      XMIN=MIN(XMIN,XL(IAT))
      XMAX=MAX(XMAX,XL(IAT))
      YMIN=MIN(YMIN,YL(IAT))
      YMAX=MAX(YMAX,YL(IAT))
      ZMIN=MIN(ZMIN,ZL(IAT))
      ZMAX=MAX(ZMAX,ZL(IAT))
      END DO
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE  XMSKOWR(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &            QHERM,RRHO,RHOMASK,XRCELL,
     &            XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &            AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,OMASK,
     &            OUNIT,EXTEND)
C
C Write mask in compressed O format
C
C Author: Paul Adams
C ==================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'xmask.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      LOGICAL QHERM
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRCELL(*)
      INTEGER AMIN, AMAX, BMIN, BMAX, CMIN, CMAX
      REAL OMASK(AMAX-AMIN+1,BMAX-BMIN+1,CMAX-CMIN+1)
      INTEGER OUNIT
      CHARACTER*(*) EXTEND
C local
      INTEGER A, B, C, SYM
      INTEGER KSECT, ISECT, JSECT, I
      INTEGER AA1, AA2, AA3, BB1, BB2, BB3
      INTEGER AS, BS, CS, ASS, BSS, CSS, ASSS, BSSS, CSSS
      INTEGER LTX, LTY, LTZ, SYMMOP
      INTEGER MASKIN, MASKOUT, TMASKIN, TMASKOUT, LENGTH
      DOUBLE PRECISION SECS, SECS2
C parameters
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
C begin
C
C initialize timers
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
C write map headers
      IF (OUNIT.NE.0) THEN
      WRITE(OUNIT,'(A)')'.MASK_INPUT'
      WRITE(OUNIT,'(A,3I6)')
     &  'ORIGIN ',AMIN,BMIN,CMIN
      WRITE(OUNIT,'(A,3I6)')
     &  'EXTENT ',AMAX-AMIN+1,BMAX-BMIN+1,CMAX-CMIN+1
      WRITE(OUNIT,'(A,3I6)')
     &  'GRID ',NA,NB,NC
      WRITE(OUNIT,'(A,6F10.5)') 'CELL',(XRCELL(I),I=1,6)
      WRITE(OUNIT,'(A)') 'COMPRESSED'
      END IF
C
C zero overall statistics accumulators
      TMASKIN=0
      TMASKOUT=0
C
C write density matrix, c is slowest ("z-sections").
      KSECT=0
      DO C=CMIN,CMAX
      KSECT=KSECT+1
C
C zero section statistics accumulators
      MASKIN=0
      MASKOUT=0
C
C zero points in section first
      JSECT=0
      DO B=BMIN,BMAX
      JSECT=JSECT+1
      ISECT=0
      DO A=AMIN,AMAX
      ISECT=ISECT+1
      OMASK(ISECT,JSECT,KSECT)=0
      END DO
      END DO
C
      DO SYM=1,XRNSYM
C
      AA1=XRSYMM(SYM,1,3)*C+(XRSYMM(SYM,1,4)*NA)/XRSYTH
      AA2=XRSYMM(SYM,2,3)*C+(XRSYMM(SYM,2,4)*NB)/XRSYTH
      AA3=XRSYMM(SYM,3,3)*C+(XRSYMM(SYM,3,4)*NC)/XRSYTH
C
      JSECT=0
      DO B=BMIN,BMAX
      BB1=AA1+XRSYMM(SYM,1,2)*B+10000*NA-MAASY
      BB2=AA2+XRSYMM(SYM,2,2)*B+10000*NB-MBASY
      BB3=AA3+XRSYMM(SYM,3,2)*B+10000*NC-MCASY
      JSECT=JSECT+1
C
      ISECT=0
      DO A=AMIN,AMAX
C compute symmetry-related indices and project into ASU brick
      AS=MOD(XRSYMM(SYM,1,1)*A+BB1,NA) +MAASY
      BS=MOD(XRSYMM(SYM,2,1)*A+BB2,NB) +MBASY
      CS=MOD(XRSYMM(SYM,3,1)*A+BB3,NC) +MCASY
      ISECT=ISECT+1
      IF (AS.GE.MAASY.AND.BS.GE.MBASY.AND.CS.GE.MCASY.AND.
     &    AS.LE.NAASY.AND.BS.LE.NBASY.AND.CS.LE.NCASY) THEN
      IF (RHOMASK(AS,BS,CS).GT.0) THEN
C
C if this is an averaging map (points in mask are negative) then
C check this is really a point around molecule and not xtal symmetry
C map the symmetry point back using the stored remapping operators
C does this point match the one we are trying to fill?
      IF (RRHO(AS,BS,CS).LE.0) THEN
         IF (EXTEND.EQ.'MOLE') THEN
            CALL MSKA2M(RRHO(AS,BS,CS),SYMMOP,LTX,LTY,LTZ,XRMSYM)
            ASS=AS-(LTX*NA)-(XRSYMM(SYMMOP,1,4)*NA)/XRSYTH
            BSS=BS-(LTY*NB)-(XRSYMM(SYMMOP,2,4)*NB)/XRSYTH
            CSS=CS-(LTZ*NC)-(XRSYMM(SYMMOP,3,4)*NC)/XRSYTH
            ASSS=XRITSY(SYMMOP,1,1)*ASS+XRITSY(SYMMOP,2,1)*BSS
     &           +XRITSY(SYMMOP,3,1)*CSS
            BSSS=XRITSY(SYMMOP,1,2)*ASS+XRITSY(SYMMOP,2,2)*BSS
     &           +XRITSY(SYMMOP,3,2)*CSS
            CSSS=XRITSY(SYMMOP,1,3)*ASS+XRITSY(SYMMOP,2,3)*BSS
     &           +XRITSY(SYMMOP,3,3)*CSS
C
C if its xtal symmetry then write 0 as point is not in mask
            IF (ASSS.EQ.A.AND.BSSS.EQ.B.AND.CSSS.EQ.C) THEN
               OMASK(ISECT,JSECT,KSECT)=1
               MASKIN=MASKIN+1
            ELSE
               OMASK(ISECT,JSECT,KSECT)=0
               MASKOUT=MASKOUT+1
            END IF
         ELSE
            OMASK(ISECT,JSECT,KSECT)=1
            MASKIN=MASKIN+1
         END IF
C
C if point is zero or one then we have to write it as is
C unless the point in the section has already been written as non-zero
C hence non-averaging masks will be expanded by xtal symmetry
      ELSE
         IF (RRHO(AS,BS,CS).LE.0) THEN
            MASKIN=MASKIN+1
            OMASK(ISECT,JSECT,KSECT)=1
         END IF
         IF (RRHO(AS,BS,CS).GT.0) THEN
            MASKOUT=MASKOUT+1
            OMASK(ISECT,JSECT,KSECT)=0
         END IF
      END IF
C
      END IF
      END IF
      END DO
      END DO
      END DO
C
      IF (OUNIT.GT.0) THEN
C print info about the section
      WRITE(6,'(A,I4,A,I8,A,I8)')
     & ' XMSKWR: mask section #',KSECT,' points inside mask',MASKIN,
     & ' points outside mask',MASKOUT
      END IF
C
C update overall statistics counters
      TMASKIN=TMASKIN+MASKIN
      TMASKOUT=TMASKOUT+MASKOUT
C
      END DO
C
C print info about the whole map
      WRITE(6,'(A,I8,A,I8)')
     & ' XMSKWR: total # points inside mask',TMASKIN,
     & ' total # points outside mask',TMASKOUT
C
C compress and write the mask
      LENGTH=(AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)
      CALL MSKOCOMPW(OMASK,LENGTH,OUNIT)
C
C terminate file with END card
      IF (OUNIT.GT.0) THEN
      WRITE(OUNIT,'(A)') 'END'
      END IF
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4,/,A,F10.4)')
     &        ' XMSKWR: CPU-time: mask write=',SECS2-SECS
      END IF
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE MSKOCOMPW(OMASK,LENGTH,OUNIT)
C
C write a compressed O mask
C
C Paul Adams and Gerard Kleywegt
C
      IMPLICIT NONE
C
      INTEGER LENGTH, OUNIT
      REAL OMASK(LENGTH)
C local
      INTEGER I, START, FINISH
      LOGICAL FOUND, ISONE, ERROR
C
      ERROR=.FALSE.
C
      START=0
      FINISH=0
      FOUND=.FALSE.
      ISONE=.FALSE.
C
      DO I=1,LENGTH
         IF (.NOT.ISONE.AND.OMASK(I).EQ.1) THEN
            FOUND=.TRUE.
            ISONE=.TRUE.
            START=I
         ELSE IF (ISONE.AND.OMASK(I).NE.1) THEN
            ISONE=.FALSE.
            FINISH=I-1
            WRITE(OUNIT,'(2I10)') START,FINISH
         END IF
      END DO
      IF (ISONE.AND.OMASK(LENGTH).EQ.1) THEN
         FINISH=LENGTH
         WRITE(OUNIT,'(2I10)') START,FINISH
      END IF
C
      IF (.NOT.FOUND) THEN
         WRITE(6,'(A)')' MSKOCOMP: zero mask points written to output'
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE MSKOCOMPR(OMASK,LENGTH,IUNIT)
C
C read a compressed O mask - assume the header info has been read already
C
C Paul Adams and Gerard Kleywegt
C
      IMPLICIT NONE
C
      INTEGER LENGTH, IUNIT
      REAL OMASK(LENGTH)
C local
      INTEGER LINMAX
      PARAMETER (LINMAX=132)
C
      INTEGER DECODI
      CHARACTER*(LINMAX) LINE
      CHARACTER*20 STRING
      INTEGER I, START, FINISH, STLEN, SLENGTH
      LOGICAL ERROR, EOF, FOUND, DONE, OK
C
      FOUND=.FALSE.
      OK=.TRUE.
      ERROR=.FALSE.
      EOF=.FALSE.
C
C set the omask array to all non-mask first
      DO I=1,LENGTH
         OMASK(I)=1
      END DO
C
C read a line from the file
      READ(IUNIT,'(A)',END=6,ERR=7) LINE
C
C trim the line
      STLEN=LINMAX
      CALL TRIMM(LINE,STLEN)
      CALL TRIML(LINE,STLEN)
C
      DONE=.FALSE.
C
      DO WHILE (.NOT.DONE)
         CALL NXTSTR(LINE,LINMAX,STRING,SLENGTH,OK)
         IF (.NOT.OK)
     &   CALL WRNDIE(-5,'MSKOCOMPR','Error reading string')
         START=DECODI(STRING,SLENGTH,OK)
         IF (.NOT.OK)
     &   CALL WRNDIE(-5,'MSKOCOMPR','Error decoding integer')
         CALL LSHFTSTR(LINE,LINMAX)
         CALL NXTSTR(LINE,LINMAX,STRING,SLENGTH,OK)
         IF (.NOT.OK)
     &   CALL WRNDIE(-5,'MSKOCOMPR','Error reading string')
         FINISH=DECODI(STRING,SLENGTH,OK)
         IF (.NOT.OK)
     &   CALL WRNDIE(-5,'MSKOCOMPR','Error decoding integer')
         IF ((START.LE.0).OR.(FINISH.GT.LENGTH)) THEN
            CALL WRNDIE(-5,'MSKOCOMPR',
     &                  'Lattice points outside mask array')
         END IF
         DO I=START,FINISH
            FOUND=.TRUE.
            OMASK(I)=-1
         END DO
         READ(IUNIT,'(A)',END=6,ERR=7) LINE
         STLEN=LINMAX
         CALL TRIMM(LINE,STLEN)
         CALL TRIML(LINE,STLEN)
         IF (LINE(1:3).EQ.'END'.OR.LINE(1:3).EQ.'end') DONE=.TRUE.
      END DO
C
      IF (.NOT.FOUND) THEN
         CALL WRNDIE(-5,'MSKOCOMP',
     &     'zero mask points read from input')
      END IF
C
C-error-labels----
      GOTO 77
7     ERROR=.TRUE.
      GOTO 77
6     EOF=.TRUE.
77    CONTINUE
C
      IF (ERROR) WRITE(6,'(A)') ' %MSKOCOMPR-ERR: error during read'
      IF (EOF) WRITE(6,'(A)') ' %MSKOCOMPR-ERR: EOF during read'
C
      RETURN
      END
C======================================================================
      SUBROUTINE ROMSKHDR(FILE,U,ERROR,EOF,
     &                    NA,AMIN,AMAX,NB,BMIN,BMAX,
     &                    NC,CMIN,CMAX,LCELL)
C
C Subroutine reads header information from a compressed O mask
C
C Paul Adams
C
      IMPLICIT NONE
C I/O
      CHARACTER*(*) FILE
      INTEGER U
      LOGICAL ERROR, EOF
      INTEGER NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX
      DOUBLE PRECISION LCELL(6)
C local
      INTEGER LINMAX
      PARAMETER (LINMAX=132)
C
      INTEGER DECODI
      DOUBLE PRECISION DECODF
      CHARACTER*(LINMAX) LINE
      CHARACTER*20 STRING
      INTEGER I, ATEMP, BTEMP, CTEMP, STLEN, LENGTH, INT(3)
      LOGICAL QORIG, QEXT, QGRID, QCELL, QCOMP, DONE, OK
C begin
      OK=.TRUE.
      ERROR=.FALSE.
      EOF=.FALSE.
      CALL ASSFIL(FILE,U,'READ','FORMATTED',ERROR)
      IF (.NOT.ERROR) THEN
C
C read stuff until we hit the COMPRESSED flag
      QORIG=.FALSE.
      QEXT=.FALSE.
      QGRID=.FALSE.
      QCELL=.FALSE.
      QCOMP=.FALSE.
C
      DONE=.FALSE.
C
      DO WHILE (.NOT.DONE)
         READ(U,'(A)',END=6,ERR=7) LINE
         STLEN=LINMAX
         CALL TRIMM(LINE,STLEN)
         CALL TRIML(LINE,STLEN)
         IF (LINE(1:6).EQ.'ORIGIN') THEN
            QORIG=.TRUE.
            DO I=1,3
               CALL LSHFTSTR(LINE,LINMAX)
               CALL NXTSTR(LINE,LINMAX,STRING,LENGTH,OK)
               IF (.NOT.OK)
     &         CALL WRNDIE(-5,'ROMSKHDR','Error reading string')
               INT(I)=DECODI(STRING,LENGTH,OK)
               IF (.NOT.OK)
     &         CALL WRNDIE(-5,'ROMSKHDR','Error decoding integer')
            END DO
            AMIN=INT(1)
            BMIN=INT(2)
            CMIN=INT(3)
         ELSE IF (LINE(1:6).EQ.'EXTENT') THEN
            QEXT=.TRUE.
            DO I=1,3
               CALL LSHFTSTR(LINE,LINMAX)
               CALL NXTSTR(LINE,LINMAX,STRING,LENGTH,OK)
               IF (.NOT.OK)
     &         CALL WRNDIE(-5,'ROMSKHDR','Error reading string')
               INT(I)=DECODI(STRING,LENGTH,OK)
               IF (.NOT.OK)
     &         CALL WRNDIE(-5,'ROMSKHDR','Error decoding integer')
            END DO
            ATEMP=INT(1)
            BTEMP=INT(2)
            CTEMP=INT(3)
         ELSE IF (LINE(1:4).EQ.'GRID') THEN
            QGRID=.TRUE.
            DO I=1,3
               CALL LSHFTSTR(LINE,LINMAX)
               CALL NXTSTR(LINE,LINMAX,STRING,LENGTH,OK)
               IF (.NOT.OK)
     &         CALL WRNDIE(-5,'ROMSKHDR','Error reading string')
               INT(I)=DECODI(STRING,LENGTH,OK)
               IF (.NOT.OK)
     &         CALL WRNDIE(-5,'ROMSKHDR','Error decoding integer')
            END DO
            NA=INT(1)
            NB=INT(2)
            NC=INT(3)
         ELSE IF (LINE(1:4).EQ.'CELL') THEN
            QCELL=.TRUE.
            DO I=1,6
               CALL LSHFTSTR(LINE,LINMAX)
               CALL NXTSTR(LINE,LINMAX,STRING,LENGTH,OK)
               IF (.NOT.OK)
     &         CALL WRNDIE(-5,'ROMSKHDR','Error reading string')
               LCELL(I)=DECODF(STRING,LENGTH,OK)
               IF (.NOT.OK)
     &         CALL WRNDIE(-5,'ROMSKHDR','Error decoding real')
            END DO
         ELSE IF (LINE(1:10).EQ.'COMPRESSED') THEN
            QCOMP=.TRUE.
            DONE=.TRUE.
         ELSE IF (LINE(1:3).EQ.'END'.OR.LINE(1:3).EQ.'end') THEN
            DONE=.TRUE.
         END IF
      END DO
C
      ELSE
         CALL WRNDIE(-5,'ROMSKHDR','Error opening mask file')
      END IF
C
C-error-labels----
      GOTO 77
 7    ERROR=.TRUE.
      GOTO 77
 6    EOF=.TRUE.
 77   CONTINUE
C
C check we got all the right stuff
      IF (.NOT.QORIG) THEN
         CALL WRNDIE(-5,'ROMSKHDR','ORIGIN undefined in mask file')
      ELSE IF (.NOT.QEXT) THEN
         CALL WRNDIE(-5,'ROMSKHDR','EXTENT undefined in mask file')
      ELSE IF (.NOT.QGRID) THEN
         CALL WRNDIE(-5,'ROMSKHDR','GRID undefined in mask file')
      ELSE IF (.NOT.QCELL) THEN
         CALL WRNDIE(-5,'ROMSKHDR','CELL undefined in mask file')
      ELSE IF (.NOT.QCOMP) THEN
         CALL WRNDIE(-5,'ROMSKHDR','data not marked as compressed')
      END IF
C
      AMAX=AMIN+ATEMP-1
      BMAX=BMIN+BTEMP-1
      CMAX=CMIN+CTEMP-1
C
      IF (ERROR) WRITE(6,'(A)') ' %ROMSKHDR-ERR: error during read'
      IF (EOF) WRITE(6,'(A)') ' %ROMSKHDR-ERR: EOF during read'
C
      RETURN
      END
C======================================================================
      SUBROUTINE LSHFTSTR(LINE,LINMAX)
C
      IMPLICIT NONE
C
      INTEGER LINMAX
      CHARACTER*(*) LINE
C local
      INTEGER POINT, SPPNT, CHPNT, I
      LOGICAL DONE
C
      POINT=1
      DONE=.FALSE.
      DO WHILE (.NOT.DONE)
         IF (LINE(POINT:POINT).EQ.' ') THEN
            SPPNT=POINT
            DONE=.TRUE.
         ELSE IF (POINT.GE.LINMAX) THEN
            SPPNT=LINMAX
            DONE=.TRUE.
         END IF
         POINT=POINT+1
      END DO
C
      POINT=SPPNT
      DONE=.FALSE.
      DO WHILE (.NOT.DONE)
         IF (LINE(POINT:POINT).NE.' ') THEN
            CHPNT=POINT
            DONE=.TRUE.
         ELSE IF (POINT.GE.LINMAX) THEN
            CHPNT=LINMAX
            DONE=.TRUE.
         END IF
         POINT=POINT+1
      END DO
C
      IF (SPPNT.LT.LINMAX) THEN
         LINE(1:(LINMAX-SPPNT))=LINE(CHPNT:LINMAX)
         DO I=LINMAX-SPPNT+1,LINMAX
            LINE(I:I)=' '
         END DO
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE NXTSTR(LINE,LINMAX,STRING,LENGTH,ERROR)
C
      IMPLICIT NONE
C
      INTEGER LENGTH, LINMAX
      CHARACTER*(*) LINE
      CHARACTER*20 STRING
      LOGICAL ERROR
C local
      INTEGER POINT, I
      LOGICAL DONE
C
      POINT=1
      DONE=.FALSE.
      DO WHILE (.NOT.DONE)
         IF (LINE(POINT:POINT).EQ.' ') THEN
            LENGTH=POINT-1
            DONE=.TRUE.
         ELSE IF (POINT.GE.LINMAX) THEN
            LENGTH=LINMAX
            DONE=.TRUE.
         END IF
         POINT=POINT+1
      END DO
C
      IF (LENGTH.LE.20) THEN
         STRING(1:LENGTH)=LINE(1:LENGTH)
         DO I=(LENGTH+1),20
            STRING(I:I)=' '
         END DO
      ELSE
         ERROR=.TRUE.
      END IF
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE GRIDSZ(LNA,LNB,LNC,LAMIN,LAMAX,LBMIN,
     &                  LBMAX,LCMIN,LCMAX,NA,NB,NC,
     &                  AMIN,AMAX,BMIN,BMAX,CMIN,CMAX)
C
      IMPLICIT NONE
C IO
      INTEGER LNA,LNB,LNC
      INTEGER LAMIN,LAMAX,LBMIN,LBMAX,LCMIN,LCMAX
      INTEGER NA,NB,NC,AMIN,AMAX,BMIN,BMAX,CMIN,CMAX
C local
      DOUBLE PRECISION FAMIN,FAMAX,FBMIN,FBMAX,FCMIN,FCMAX
      DOUBLE PRECISION TEMP
C
      FAMIN=DBLE(LAMIN)/LNA
      FAMAX=DBLE(LAMAX)/LNA
      FBMIN=DBLE(LBMIN)/LNB
      FBMAX=DBLE(LBMAX)/LNB
      FCMIN=DBLE(LCMIN)/LNC
      FCMAX=DBLE(LCMAX)/LNC
C
      AMIN=INT(NA*FAMIN)
      AMAX=INT(NA*FAMAX)
      BMIN=INT(NB*FBMIN)
      BMAX=INT(NB*FBMAX)
      CMIN=INT(NC*FCMIN)
      CMAX=INT(NC*FCMAX)
C
      TEMP=DBLE(AMIN)/NA
      IF (TEMP.GE.FAMIN) AMIN=AMIN-1
      TEMP=DBLE(AMAX)/NA
      IF (TEMP.LE.FAMAX) AMAX=AMAX+1
      TEMP=DBLE(BMIN)/NB
      IF (TEMP.GE.FBMIN) BMIN=BMIN-1
      TEMP=DBLE(BMAX)/NB
      IF (TEMP.LE.FBMAX) BMAX=BMAX+1
      TEMP=DBLE(CMIN)/NC
      IF (TEMP.GE.FCMIN) CMIN=CMIN-1
      TEMP=DBLE(CMAX)/NC
      IF (TEMP.LE.FCMAX) CMAX=CMAX+1
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE REGRID(IMASK,NA,NB,NC,AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,
     &                  MASK1,AMIN1,AMAX1,BMIN1,BMAX1,CMIN1,CMAX1,
     &                  MASK2,NA2,NB2,NC2,
     &                  AMIN2,AMAX2,BMIN2,BMAX2,CMIN2,CMAX2)
C
      IMPLICIT NONE
C IO
      INTEGER NA,NB,NC,AMIN,AMAX,BMIN,BMAX,CMIN,CMAX
      INTEGER AMIN1,AMAX1,BMIN1,BMAX1,CMIN1,CMAX1
      INTEGER NA2,NB2,NC2,AMIN2,AMAX2,BMIN2,BMAX2,CMIN2,CMAX2
      REAL IMASK(AMIN:AMAX,BMIN:BMAX,CMIN:CMAX)
      REAL MASK1(AMIN1:AMAX1,BMIN1:BMAX1,CMIN1:CMAX1)
      REAL MASK2(AMIN2:AMAX2,BMIN2:BMAX2,CMIN2:CMAX2)
C local
      INTEGER A,B,C,LAMIN,LAMAX,LBMIN,LBMAX,LCMIN,LCMAX
      DOUBLE PRECISION INTERP,FA,FB,FC,DX,DY,DZ
C
C fill MASK1 (the source for interpolation) with the input mask (IMASK)
C fill values outside IMASK with one
      DO C=CMIN1,CMAX1
         DO B=BMIN1,BMAX1
            DO A=AMIN1,AMAX1
               IF (A.GE.AMIN.AND.B.GE.BMIN.AND.C.GE.CMIN.AND.
     &             A.LE.AMAX.AND.B.LE.BMAX.AND.C.LE.CMAX) THEN
                  MASK1(A,B,C)=IMASK(A,B,C)
               ELSE
                  MASK1(A,B,C)=1
               END IF
            END DO
         END DO
      END DO
C
C loop over all grid points in MASK2 (the destination for interpolation)
      DO C=CMIN2,CMAX2
      DO B=BMIN2,BMAX2
      DO A=AMIN2,AMAX2
C
C calculate nearest grid points in MASK1
      FA=DBLE(NA)*A/NA2
      FB=DBLE(NB)*B/NB2
      FC=DBLE(NC)*C/NC2
      LAMIN=MIN(INT(FA),(INT(FA)+INT(SIGN(1.0D0,FA))))
      LAMAX=MAX(INT(FA),(INT(FA)+INT(SIGN(1.0D0,FA))))
      LBMIN=MIN(INT(FB),(INT(FB)+INT(SIGN(1.0D0,FB))))
      LBMAX=MAX(INT(FB),(INT(FB)+INT(SIGN(1.0D0,FB))))
      LCMIN=MIN(INT(FC),(INT(FC)+INT(SIGN(1.0D0,FC))))
      LCMAX=MAX(INT(FC),(INT(FC)+INT(SIGN(1.0D0,FC))))
C
C check to see if points are within the original mask grid
      IF (.NOT.(LAMIN.GE.AMIN1.AND.LBMIN.GE.BMIN1.AND.
     &          LCMIN.GE.CMIN1.AND.LAMAX.LE.AMAX1.AND.
     &          LBMAX.LE.BMAX1.AND.LCMAX.LE.CMAX1)) THEN
         WRITE(6,'(6(A,I4))')' REGRID: box ',LAMIN,' ',LAMAX,' ',
     &                 LBMIN,' ',LBMAX,' ',LCMIN,' ',LCMAX
         CALL WRNDIE(-5,'REGRID','points outside input mask')
      END IF
C
C distance from point to the grid
      DX=FA-LAMIN
      DY=FB-LBMIN
      DZ=FC-LCMIN
C
      INTERP=MASK1(LAMIN,LBMIN,LCMIN)*(1-DX)*(1-DY)*(1-DZ)+
     &       MASK1(LAMAX,LBMIN,LCMIN)*DX*(1-DY)*(1-DZ)+
     &       MASK1(LAMIN,LBMAX,LCMIN)*(1-DX)*DY*(1-DZ)+
     &       MASK1(LAMAX,LBMAX,LCMIN)*DX*DY*(1-DZ)+
     &       MASK1(LAMIN,LBMIN,LCMAX)*(1-DX)*(1-DY)*DZ+
     &       MASK1(LAMAX,LBMIN,LCMAX)*DX*(1-DY)*DZ+
     &       MASK1(LAMIN,LBMAX,LCMAX)*(1-DX)*DY*DZ+
     &       MASK1(LAMAX,LBMAX,LCMAX)*DX*DY*DZ
C
      IF (INTERP.LT.0.0D0) THEN
         MASK2(A,B,C)=-1
      ELSE
         MASK2(A,B,C)=1
      END IF
C
      END DO
      END DO
      END DO
C
      RETURN
      END
C
C======================================================================
C
