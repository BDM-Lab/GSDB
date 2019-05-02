      SUBROUTINE XPEAKPIK(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &  MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &  XRITSY,QHERM,XRCELL,XRINTR,XRTR,
     &  XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &  HPRHOMA,NRHO,IRHO,NMASK)
C
C Peak picking.
C
C Author: Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'mtf.inc'
      INTEGER NA, NB, NC, NAP, NBP, NCP, NAPP, NBPP, NCPP
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9)
      DOUBLE PRECISION XRINTR(3,3), XRTR(3,3)
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER HPRHOMA, IRHO, NRHO, NMASK
C local
      INTEGER IFROM, I, MPEAK, NPEAK, NSEL, MSEL
      INTEGER NNSELE, NRHOEX, OUNIT, LENG
      LOGICAL ERR, COND, QATOM, QSYMB, QFRAC
C pointer
      INTEGER HPRMAP, HPIMAP, EXPAND, MPACK, ISEL
      INTEGER NEXT, LHPEAK, LXPEAK, LYPEAK, LZPEAK
      INTEGER LHGRID, LXGRID, LYGRID, LZGRID, ORDR, XF, YF, ZF
C begin
C default values
      ERR=.FALSE.
      IFROM=1
      MPEAK=100
      OFILE='OUTPUT'
      OUNIT=6
      QATOM=.FALSE.
      QSYMB=.FALSE.
      QFRAC=.FALSE.
C
C make default map-selection:
      IF (HPRHOMA.EQ.0) THEN
      CALL WRNDIE(-5,'XPEAKPIK','No maps are defined.')
      MPACK=0
      ELSE
      MPACK=ALLHP(INTEG4(NRHO))
      CALL XMPPCK0(HEAP(MPACK),NNSELE,
     & MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,HEAP(HPRHOMA))
      END IF
C
C make default atom-selection:
      MSEL=NATOM
      ISEL=ALLHP(INTEG4(MSEL))
      NSEL=0
C
C parsing
      CALL PUSEND('PEAKpik>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('PEAKpik>')
      CALL MISCOM('PEAKpik>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-peakpick')
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
      WRITE(6,'(3A)') ' %XMHISTO-ERR: real space object ',
     & WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'XMHISTO','object undeclared.')
      ERR=.TRUE.
      END IF
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MPEA') THEN
      CALL NEXTI('MPEAks=',MPEAK)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ATOM') THEN
      CALL NEXTLO('ATOM=',QATOM)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SYMB') THEN
      CALL NEXTLO('SYMBols=',QSYMB)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FRAC') THEN
      CALL NEXTLO('FRACtional=',QFRAC)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
      CALL ASSFIL(OFILE,OUNIT,'WRITE','FORMATTED',ERR)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SELE'.OR.WD(1:1).EQ.'(') THEN
      CALL XMPSELE(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &      QHERM,XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &      HPRHOMA,NRHO,NMASK,
     &      XRNSYM,HEAP(MPACK),NNSELE)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PROX') THEN
      CALL SELCTA(HEAP(ISEL),NSEL,X,Y,Z,.TRUE.)
      CALL MAKIND(HEAP(ISEL),NATOM,NSEL)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' |-----------------------',
     & ' peakpik ------------------------------------'
      WRITE(6,'(3A,I6)') ' | FROM = ',XRHONAM(IFROM),' MPEAk=',
     &   MPEAK
      LENG=LEN(OFILE)
      CALL TRIMM(OFILE,LENG)
      WRITE(6,'(2A)') ' | OUTPut=',OFILE(1:LENG)
      WRITE(6,'(A,L1,A,I6)')  ' | ATOM=',QATOM,
     &   ' selected proximity atoms =',NSEL
      WRITE(6,'(2A)') ' |-----------------------',
     & '---------------------------------------------'
C=====================================================================
      ELSE
      CALL CHKEND('PEAKpik>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (MPACK.NE.0) THEN
C
      IF (.NOT.ERR.AND.HPRRHO(IFROM).EQ.0) THEN
      WRITE(6,'(3A)') ' %XPEAKPIK-ERR: real space object ',
     & XRHONAM(IFROM),' undefined.'
      CALL WRNDIE(-5,'XPEAKPIK','object undefined.')
      ELSEIF (.NOT.ERR) THEN
      HPRMAP=HPRRHO(IFROM)
      HPIMAP=HPIRHO(IFROM)
C
      IF (HPRMAP.NE.0.AND.HPRHOMA.NE.0) THEN
C
      IF (NNSELE.GT.0) THEN
C
C Expand asymmetric unit by cushion (one layer) by applying
C symmetry operations.  In this way we don't have to worry
C about surface problems during peak-picking and interpolation.
C
C allocate space for expanded map
      NRHOEX=(NAASY-MAASY+3)*(NBASY-MBASY+3)*(NCASY-MCASY+3)
C
      EXPAND=ALLHP(IREAL4(NRHOEX))
      CALL XPEAKEX(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &  MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &  XRITSY,QHERM,HEAP(HPRMAP),HEAP(EXPAND),HEAP(HPRHOMA),
     &  HEAP(MPACK))
C
C search for peaks and interpolate peaks
C
C allocate space
      NEXT=ALLHP(INTEG4(MPEAK+1))
      LHPEAK=ALLHP(IREAL8(MPEAK))
      LXPEAK=ALLHP(IREAL8(MPEAK))
      LYPEAK=ALLHP(IREAL8(MPEAK))
      LZPEAK=ALLHP(IREAL8(MPEAK))
      LHGRID=ALLHP(IREAL8(MPEAK))
      LXGRID=ALLHP(IREAL8(MPEAK))
      LYGRID=ALLHP(IREAL8(MPEAK))
      LZGRID=ALLHP(IREAL8(MPEAK))
      ORDR=ALLHP(INTEG4(MPEAK))
      CALL XPEAKFN(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,
     &      MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &      HEAP(EXPAND),HEAP(HPRHOMA),HEAP(MPACK),
     &      XRINTR,NPEAK,MPEAK,
     &      HEAP(NEXT),
     &      HEAP(LHPEAK),HEAP(LXPEAK),HEAP(LYPEAK),
     &      HEAP(LZPEAK),HEAP(LHGRID),HEAP(LXGRID),HEAP(LYGRID),
     &      HEAP(LZGRID),HEAP(ORDR))
C
C if PROXimity atoms are specified, apply all symmetry and
C lattice operations
C in order to get the equivalent peak position  closest to
C the selected atoms
C
      IF (NSEL.GT.0) THEN
      XF=ALLHP(IREAL8(NSEL))
      YF=ALLHP(IREAL8(NSEL))
      ZF=ALLHP(IREAL8(NSEL))
      CALL XPEAKSL(OUNIT,NPEAK,MPEAK,
     &      HEAP(NEXT),
     &      HEAP(LHPEAK),HEAP(LXPEAK),HEAP(LYPEAK),
     &      HEAP(LZPEAK),HEAP(LHGRID),HEAP(LXGRID),HEAP(LYGRID),
     &      HEAP(LZGRID),HEAP(ORDR),NSEL,HEAP(ISEL),HEAP(XF),HEAP(YF),
     &      HEAP(ZF),XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &      XRITSY,QHERM,XRCELL,XRINTR,XRTR,X,Y,Z)
      CALL FREHP(ZF,IREAL8(NSEL))
      CALL FREHP(YF,IREAL8(NSEL))
      CALL FREHP(XF,IREAL8(NSEL))
      END IF
C
C print peak list and store info about highest peak in symbols
C $PEAK_HEIGHT, $PEAK_X, $PEAK_Y, $PEAK_Z, $NPEAKS
      CALL XPEAKPR(OUNIT,NPEAK,MPEAK,
     &      HEAP(NEXT),
     &      HEAP(LHPEAK),HEAP(LXPEAK),HEAP(LYPEAK),
     &      HEAP(LZPEAK),HEAP(LHGRID),HEAP(LXGRID),HEAP(LYGRID),
     &      HEAP(LZGRID),HEAP(ORDR),QSYMB,QFRAC,XRTR)
C
C create dummy atoms if required
      IF (QATOM) THEN
      CALL XPEAKAT(NPEAK,MPEAK,
     &      HEAP(NEXT),
     &      HEAP(LHPEAK),HEAP(LXPEAK),HEAP(LYPEAK),
     &      HEAP(LZPEAK),HEAP(LHGRID),HEAP(LXGRID),HEAP(LYGRID),
     &      HEAP(LZGRID),HEAP(ORDR))
      END IF
C
C free-up space
      CALL FREHP(ORDR,INTEG4(MPEAK))
      CALL FREHP(LZGRID,IREAL8(MPEAK))
      CALL FREHP(LYGRID,IREAL8(MPEAK))
      CALL FREHP(LXGRID,IREAL8(MPEAK))
      CALL FREHP(LHGRID,IREAL8(MPEAK))
      CALL FREHP(LZPEAK,IREAL8(MPEAK))
      CALL FREHP(LYPEAK,IREAL8(MPEAK))
      CALL FREHP(LXPEAK,IREAL8(MPEAK))
      CALL FREHP(LHPEAK,IREAL8(MPEAK))
      CALL FREHP(NEXT,INTEG4(MPEAK+1))
C
      CALL FREHP(EXPAND,IREAL4(NRHOEX))
C
      END IF
C
      END IF
      END IF
C
C free-up space if necessary
      IF (MPACK.NE.0) CALL FREHP(MPACK,INTEG4(NRHO))
C
      END IF
C
      CALL FREHP(ISEL,INTEG4(MSEL))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XPEAKEX(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &  MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &  XRITSY,QHERM,RRHO,EXPAND,RHOMASK,MPACK)
C
C Expand asymmetric unit by cushion (one layer) by applying
C symmetry operations.  In this way we don't have to worry
C about surface problems during peak-picking and interpolation.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER NA, NB, NC, NAP, NBP, NCP, NAPP, NBPP, NCPP
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL EXPAND(MAASY-1:NAASY+1,MBASY-1:NBASY+1,MCASY-1:NCASY+1)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER MPACK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
C local
      INTEGER A, B, C, AA, BB, CC, AS, BS, CS, SYM
      LOGICAL COND
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
C initialize expanded map
      DO C=MCASY-1,NCASY+1
      DO B=MBASY-1,NBASY+1
      DO A=MAASY-1,NAASY+1
      EXPAND(A,B,C)=ZERO
      END DO
      END DO
      END DO
C
C loop over all selected asymmetric unit grid points
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).NE.0.AND.MPACK(A,B,C).EQ.1) THEN
C
C store this point in expanded map
      EXPAND(A,B,C)=RRHO(A,B,C)
C
C loop over a box around the asymmetric unit grid point
      DO CC=C-1,C+1
      DO BB=B-1,B+1
      DO AA=A-1,A+1
      COND=(AA.GE.MAASY.AND.AA.LE.NAASY.AND.
     &      BB.GE.MBASY.AND.BB.LE.NBASY.AND.
     &      CC.GE.MCASY.AND.CC.LE.NCASY)
      IF (COND) COND=RHOMASK(AA,BB,CC).GT.0
      IF (.NOT.COND.AND.EXPAND(AA,BB,CC).EQ.ZERO) THEN
C
C this grid point (AA,BB,CC) is in the neighborhood of
C the asymmetric unit grid point (A,B,C) but it is
C not in the asymmetric unit.
C Apply all symmetry operators to grid point until
C the symmetry-related grid point falls into the
C asymmetric unit
C
      SYM=1
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
      AS=MOD(XRSYMM(SYM,1,1)*AA+XRSYMM(SYM,1,2)*BB+XRSYMM(SYM,1,3)*CC
     &       +(XRSYMM(SYM,1,4)*NA)/XRSYTH+10000*NA-MAASY,NA) +MAASY
      BS=MOD(XRSYMM(SYM,2,1)*AA+XRSYMM(SYM,2,2)*BB+XRSYMM(SYM,2,3)*CC
     &       +(XRSYMM(SYM,2,4)*NB)/XRSYTH+10000*NB-MBASY,NB) +MBASY
      CS=MOD(XRSYMM(SYM,3,1)*AA+XRSYMM(SYM,3,2)*BB+XRSYMM(SYM,3,3)*CC
     &       +(XRSYMM(SYM,3,4)*NC)/XRSYTH+10000*NC-MCASY,NC) +MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND) THEN
      COND=RHOMASK(AS,BS,CS).GT.0
      IF(COND.AND.MPACK(AS,BS,CS).EQ.1) EXPAND(AA,BB,CC)=RRHO(AS,BS,CS)
      END IF
      SYM=SYM+1
      END DO
C
      IF (.NOT.COND) CALL WRNDIE(-5,'XPEAKPIK','Internal error no. 2')
C
      END IF
      END DO
      END DO
      END DO
      END IF
      END DO
      END DO
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XPEAKFN(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,
     &      MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &      EXPAND,RHOMASK,MPACK,XRINTR,NPEAK,MPEAK,NEXT,
     &      LHPEAK,LXPEAK,LYPEAK,LZPEAK,LHGRID,LXGRID,LYGRID,
     &      LZGRID,ORDR)
C
C Search for peaks in selected map elements and interpolate peaks.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER NA, NB, NC, NAP, NBP, NCP, NAPP, NBPP, NCPP
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      REAL EXPAND(MAASY-1:NAASY+1,MBASY-1:NBASY+1,MCASY-1:NCASY+1)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER MPACK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION XRINTR(3,3)
      INTEGER NPEAK, MPEAK, NEXT(0:MPEAK)
      DOUBLE PRECISION LHPEAK(*), LXPEAK(*), LYPEAK(*), LZPEAK(*)
      DOUBLE PRECISION LHGRID(*), LXGRID(*), LYGRID(*), LZGRID(*)
      INTEGER ORDR(*)
C local
      INTEGER A, B, C
      DOUBLE PRECISION FX, FY, FZ, FXX, FYY, FZZ, FXY, FYZ, FXZ
      DOUBLE PRECISION AMAT(3,3), BVEC(3), TEMPX, TEMPY, TEMPZ
      DOUBLE PRECISION RHOA, RHOB, RHOC, RHOD
      DOUBLE PRECISION HPEAK, XPEAK, YPEAK, ZPEAK
      DOUBLE PRECISION HGRID, XGRID, YGRID, ZGRID
      INTEGER NEW, I, IPREV, ELEMENT
      LOGICAL DONE
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, FOUR=4.0D0)
      PARAMETER (HALF=0.5D0)
C begin
C
C initialize linked list of peaks
      NEXT(0)=0
      NPEAK=0
C
C loop over all selected asymmetric unit grid points
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).NE.0.AND.MPACK(A,B,C).EQ.1) THEN
C
      RHOA = EXPAND(A,B,C)
C
C test: all neighbors have to have equal or smaller density values
C       and at least one neighbor has to has a smaller density value
      RHOB = MAX(EXPAND(A-1,B,C),EXPAND(A,B-1,C),EXPAND(A,B,C-1),
     &           EXPAND(A+1,B,C),EXPAND(A,B+1,C),EXPAND(A,B,C+1))
      IF(RHOA.GE.RHOB) THEN
      RHOC = MAX(EXPAND(A+1,B+1,C),EXPAND(A+1,B,C+1),EXPAND(A,B+1,C+1),
     &           EXPAND(A-1,B-1,C),EXPAND(A-1,B,C-1),EXPAND(A,B-1,C-1),
     &           EXPAND(A+1,B-1,C),EXPAND(A+1,B,C-1),EXPAND(A,B+1,C-1),
     &           EXPAND(A-1,B+1,C),EXPAND(A-1,B,C+1),EXPAND(A,B-1,C+1))
      IF(RHOA.GE.RHOC) THEN
      RHOD = MAX(EXPAND(A-1,B-1,C-1),EXPAND(A+1,B-1,C-1),
     &           EXPAND(A-1,B+1,C-1),EXPAND(A+1,B+1,C-1),
     &           EXPAND(A-1,B-1,C+1),EXPAND(A+1,B-1,C+1),
     &           EXPAND(A-1,B+1,C+1),EXPAND(A+1,B+1,C+1))
      IF(RHOA.GE.RHOD.AND.
     &   (RHOA.GT.RHOB.OR.RHOA.GT.RHOC.OR.RHOA.GT.RHOD)) THEN
C
C map value at grid point A, B, C
      HGRID=RHOA
C grid point A, B, C
      XGRID=A
      YGRID=B
      ZGRID=C
C convert into fractional coordinates
      XGRID=XGRID/NA
      YGRID=YGRID/NB
      ZGRID=ZGRID/NC
C convert into orthogonal coordinates
      TEMPX=XGRID
      TEMPY=YGRID
      TEMPZ=ZGRID
      XGRID=XRINTR(1,1)*TEMPX +XRINTR(1,2)*TEMPY +XRINTR(1,3)*TEMPZ
      YGRID=XRINTR(2,1)*TEMPX +XRINTR(2,2)*TEMPY +XRINTR(2,3)*TEMPZ
      ZGRID=XRINTR(3,1)*TEMPX +XRINTR(3,2)*TEMPY +XRINTR(3,3)*TEMPZ
C
C interpolate to find the local maximum
      FX = (EXPAND(A+1,B,C) -   EXPAND(A-1,B,C))/TWO
      FY = (EXPAND(A,B+1,C) -   EXPAND(A,B-1,C))/TWO
      FZ = (EXPAND(A,B,  C+1) - EXPAND(A,B,  C-1))/TWO
      FXX = (EXPAND(A-1,B,C) + EXPAND(A+1,B,C) - TWO*RHOA)
      FYY = (EXPAND(A,B-1,C) + EXPAND(A,B+1,C) - TWO*RHOA)
      FZZ = (EXPAND(A,B,  C-1) +  EXPAND(A,B,  C+1) - TWO*RHOA)
      FXY = ((EXPAND(A+1,B+1,C) + EXPAND(A-1,B-1,C))
     &    -  (EXPAND(A+1,B-1,C) + EXPAND(A-1,B+1,C)))/FOUR
      FXZ = ((EXPAND(A+1,B,  C+1) + EXPAND(A-1,B,  C-1))
     &    -  (EXPAND(A+1,B,  C-1) + EXPAND(A-1,B,  C+1)))/FOUR
      FYZ = ((EXPAND(A  ,B+1,C+1) + EXPAND(A  ,B-1,C-1))
     &    -  (EXPAND(A  ,B+1,C-1) + EXPAND(A  ,B-1,C+1)))/FOUR
      AMAT(1,1) = FXX
      AMAT(1,2) = FXY
      AMAT(1,3) = FXZ
      AMAT(2,1) = FXY
      AMAT(2,2) = FYY
      AMAT(2,3) = FYZ
      AMAT(3,1) = FXZ
      AMAT(3,2) = FYZ
      AMAT(3,3) = FZZ
      BVEC(1) = -FX
      BVEC(2) = -FY
      BVEC(3) = -FZ
      CALL SOLV3(AMAT,BVEC)
      IF(ABS(BVEC(1)).GE.ONE.OR.ABS(BVEC(2)).GE.ONE
     &   .OR.ABS(BVEC(3)).GE.ONE) THEN
      BVEC(1)=ZERO
      BVEC(2)=ZERO
      BVEC(3)=ZERO
      END IF
C
C interpolated peak height
      HPEAK=RHOA + FX*BVEC(1) + FY*BVEC(2) + FZ*BVEC(3)
     & +HALF*(FXX*BVEC(1)**2+FYY*BVEC(2)**2+FZZ*BVEC(3)**2)
     & +FXY*BVEC(1)*BVEC(2)+FXZ*BVEC(1)*BVEC(3)+FYZ*BVEC(2)*BVEC(3)
C
C interpolated peak position
      XPEAK=A+BVEC(1)
      YPEAK=B+BVEC(2)
      ZPEAK=C+BVEC(3)
C convert into fractional coordinates
      XPEAK=XPEAK/NA
      YPEAK=YPEAK/NB
      ZPEAK=ZPEAK/NC
C convert into orthogonal coordinates
      TEMPX=XPEAK
      TEMPY=YPEAK
      TEMPZ=ZPEAK
      XPEAK=XRINTR(1,1)*TEMPX +XRINTR(1,2)*TEMPY +XRINTR(1,3)*TEMPZ
      YPEAK=XRINTR(2,1)*TEMPX +XRINTR(2,2)*TEMPY +XRINTR(2,3)*TEMPZ
      ZPEAK=XRINTR(3,1)*TEMPX +XRINTR(3,2)*TEMPY +XRINTR(3,3)*TEMPZ
C
C try to insert this peak into linked list of peaks
C =================================================
C
C go through current list of all peaks
      I=0
      DONE=.FALSE.
      NEW=-1
C
      DO WHILE (.NOT.DONE)
C
      IPREV=I
      I=NEXT(I)
C
C Condition: Reached end of linked list. ==> Append to list.
C            Height is greater than current element ==> Insert in list.
      DONE=I.EQ.0
C
C defines order on list of peaks <====
      IF (.NOT.DONE) DONE=HPEAK.LT.LHPEAK(I)
C
      IF (DONE) THEN
C
C maximum number of peaks <===
      IF (NPEAK.LT.MPEAK) THEN
C
C allocate new element for list
      NPEAK=NPEAK+1
      NEW=NPEAK
C
      ELSEIF (IPREV.NE.0) THEN
C
C overwrite first element
      NEW=NEXT(0)
C
C start list at second element
      NEXT(0)=NEXT(NEW)
C
C special condition: first element is IPREV
      IF (IPREV.EQ.NEW) IPREV=0
C
C we don't do anything if the list is full and the
C new peak is less than all peaks already in the list
      END IF
C
      IF (NEW.GE.0) THEN
      NEXT(IPREV)=NEW
      NEXT(NEW)=I
C store peak height and position in element NEW <======
      LHPEAK(NEW)=HPEAK
      LXPEAK(NEW)=XPEAK
      LYPEAK(NEW)=YPEAK
      LZPEAK(NEW)=ZPEAK
      LHGRID(NEW)=HGRID
      LXGRID(NEW)=XGRID
      LYGRID(NEW)=YGRID
      LZGRID(NEW)=ZGRID
      END IF
C
C
      END IF
      END DO
C
      END IF
      END IF
      END IF
C
      END IF
      END DO
      END DO
      END DO
C
C define order (heighest peak first)
      ELEMENT=NEXT(0)
      DO I=1,NPEAK
      ORDR(NPEAK-I+1)=ELEMENT
      ELEMENT=NEXT(ELEMENT)
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE SOLV3(A,B)
C
C Solves the system of three linear equations, AX=B
C Returns the inverse of A in A and the solution X in B.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      DOUBLE PRECISION A(3,3), B(3)
C local
      DOUBLE PRECISION COFA1, COFA2, COFA3, COFA4, COFA5, COFA6
      DOUBLE PRECISION DETER, RDET, B1, B2, B3
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
      COFA1=A(2,2)*A(3,3)-A(2,3)*A(3,2)
      COFA2=A(2,3)*A(3,1)-A(2,1)*A(3,3)
      COFA3=A(2,1)*A(3,2)-A(2,2)*A(3,1)
      COFA4=A(1,1)*A(3,3)-A(1,3)*A(3,1)
      COFA5=A(1,2)*A(3,1)-A(1,1)*A(3,2)
      COFA6=A(1,1)*A(2,2)-A(1,2)*A(2,1)
C
      DETER=A(1,1)*COFA1+A(1,2)*COFA2+A(1,3)*COFA3
      IF (ABS(DETER).LT.RSMALL) THEN
      B(1)=ZERO
      B(2)=ZERO
      B(3)=ZERO
      ELSE
      RDET=ONE/DETER
C
      A(1,1)=COFA1*RDET
      A(1,2)=COFA2*RDET
      A(1,3)=COFA3*RDET
      A(2,2)=COFA4*RDET
      A(2,3)=COFA5*RDET
      A(3,3)=COFA6*RDET
      A(2,1)=A(1,2)
      A(3,1)=A(1,3)
      A(3,2)=A(2,3)
C
      B1=A(1,1)*B(1)+A(1,2)*B(2)+A(1,3)*B(3)
      B2=A(2,1)*B(1)+A(2,2)*B(2)+A(2,3)*B(3)
      B3=A(3,1)*B(1)+A(3,2)*B(2)+A(3,3)*B(3)
      B(1)=B1
      B(2)=B2
      B(3)=B3
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XPEAKSL(OUNIT,NPEAK,MPEAK,NEXT,
     &      LHPEAK,LXPEAK,LYPEAK,LZPEAK,LHGRID,LXGRID,LYGRID,
     &      LZGRID,ORDR,NSEL,ISEL,XF,YF,ZF,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &      XRITSY,QHERM,XRCELL,XRINTR,XRTR,X,Y,Z)
C
C Applies all symmetry and lattice operations to peak list in order
C to find the equivalent position closest to the selected atoms
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INTEGER OUNIT, NPEAK, MPEAK, NEXT(0:MPEAK)
      DOUBLE PRECISION LHPEAK(*), LXPEAK(*), LYPEAK(*), LZPEAK(*)
      DOUBLE PRECISION LHGRID(*), LXGRID(*), LYGRID(*), LZGRID(*)
      INTEGER ORDR(*), NSEL, ISEL(*)
      DOUBLE PRECISION XF(*), YF(*), ZF(*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9)
      DOUBLE PRECISION XRINTR(3,3), XRTR(3,3)
      DOUBLE PRECISION X(*), Y(*), Z(*)
C local
      INTEGER I, ISYM, IAT
      DOUBLE PRECISION XEQU, YEQU, ZEQU, XTEMP, YTEMP, ZTEMP, RTH
      DOUBLE PRECISION DX, DY, DZ, DIST2, MIND2, XFF, YFF, ZFF
      INTEGER XLATT, YLATT, ZLATT, NNSEL
C parameter
      DOUBLE PRECISION HALF, ZERO, ONE
      PARAMETER (HALF=0.5D0, ZERO=0.0D0, ONE=1.0D0)
C begin
C
      RTH=ONE/XRSYTH
C
C fractionalize the selected and known PROXimity atoms
      NNSEL=0
      DO I=1,NSEL
      IF (INITIA(ISEL(I),X,Y,Z)) THEN
      NNSEL=NNSEL+1
      XF(I)=XRTR(1,1)*X(ISEL(I)) +XRTR(1,2)*Y(ISEL(I))
     &     +XRTR(1,3)*Z(ISEL(I))
      YF(I)=XRTR(2,1)*X(ISEL(I)) +XRTR(2,2)*Y(ISEL(I))
     &     +XRTR(2,3)*Z(ISEL(I))
      ZF(I)=XRTR(3,1)*X(ISEL(I)) +XRTR(3,2)*Y(ISEL(I))
     &     +XRTR(3,3)*Z(ISEL(I))
      END IF
      END DO
C
C loop over all peaks
      DO I=1,NPEAK
C
C set minimum distance^2
      MIND2=R4BIG
C
C equivalent peak position
      XEQU=LXPEAK(I)
      YEQU=LYPEAK(I)
      ZEQU=LZPEAK(I)
C
C loop over crystallographic symmetry operators
      DO ISYM=1,XRNSYM
C
C fractionalize the peak position
      XFF=XRTR(1,1)*LXPEAK(I) +XRTR(1,2)*LYPEAK(I) +XRTR(1,3)*LZPEAK(I)
      YFF=XRTR(2,1)*LXPEAK(I) +XRTR(2,2)*LYPEAK(I) +XRTR(2,3)*LZPEAK(I)
      ZFF=XRTR(3,1)*LXPEAK(I) +XRTR(3,2)*LYPEAK(I) +XRTR(3,3)*LZPEAK(I)
C
C apply the symmetry operator
      XTEMP=XRSYMM(ISYM,1,1)*XFF
     &     +XRSYMM(ISYM,1,2)*YFF
     &     +XRSYMM(ISYM,1,3)*ZFF+XRSYMM(ISYM,1,4)*RTH
      YTEMP=XRSYMM(ISYM,2,1)*XFF
     &     +XRSYMM(ISYM,2,2)*YFF
     &     +XRSYMM(ISYM,2,3)*ZFF+XRSYMM(ISYM,2,4)*RTH
      ZTEMP=XRSYMM(ISYM,3,1)*XFF
     &     +XRSYMM(ISYM,3,2)*YFF
     &     +XRSYMM(ISYM,3,3)*ZFF+XRSYMM(ISYM,3,4)*RTH
      XFF=XTEMP
      YFF=YTEMP
      ZFF=ZTEMP
C
C loop over all selected coordinates
      DO IAT=1,NNSEL
C
C compute difference vector
      XTEMP=XFF-XF(IAT)
      YTEMP=YFF-YF(IAT)
      ZTEMP=ZFF-ZF(IAT)
C
      XLATT=ZERO
      YLATT=ZERO
      ZLATT=ZERO
C
C compute minimum image difference
      IF ((ABS(XTEMP).GE.HALF).OR.(ABS(YTEMP).GE.HALF).OR.
     &       (ABS(ZTEMP).GE.HALF)) THEN
C
      XLATT=INT(ABS(XTEMP)+HALF)*SIGN(ONE,-XTEMP)
      YLATT=INT(ABS(YTEMP)+HALF)*SIGN(ONE,-YTEMP)
      ZLATT=INT(ABS(ZTEMP)+HALF)*SIGN(ONE,-ZTEMP)
      XTEMP=XTEMP+XLATT
      YTEMP=YTEMP+YLATT
      ZTEMP=ZTEMP+ZLATT
      END IF
C
C orthogonalize distance vector
      DX=XRINTR(1,1)*XTEMP +XRINTR(1,2)*YTEMP +XRINTR(1,3)*ZTEMP
      DY=XRINTR(2,1)*XTEMP +XRINTR(2,2)*YTEMP +XRINTR(2,3)*ZTEMP
      DZ=XRINTR(3,1)*XTEMP +XRINTR(3,2)*YTEMP +XRINTR(3,3)*ZTEMP
C
C compute distance^2
      DIST2=DX*DX+DY*DY+DZ*DZ
C
C compare with minimum distance
      IF (DIST2.LT.MIND2) THEN
      MIND2=DIST2
C
C compute equivalent peak position
      XTEMP=XFF+XLATT
      YTEMP=YFF+YLATT
      ZTEMP=ZFF+ZLATT
      XEQU=XRINTR(1,1)*XTEMP +XRINTR(1,2)*YTEMP +XRINTR(1,3)*ZTEMP
      YEQU=XRINTR(2,1)*XTEMP +XRINTR(2,2)*YTEMP +XRINTR(2,3)*ZTEMP
      ZEQU=XRINTR(3,1)*XTEMP +XRINTR(3,2)*YTEMP +XRINTR(3,3)*ZTEMP
      END IF
C
      END DO
      END DO
C
C update the peak position
      LXPEAK(I)=XEQU
      LYPEAK(I)=YEQU
      LZPEAK(I)=ZEQU
C
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XPEAKPR(OUNIT,NPEAK,MPEAK,NEXT,
     &      LHPEAK,LXPEAK,LYPEAK,LZPEAK,LHGRID,LXGRID,LYGRID,
     &      LZGRID,ORDR,QSYMB,QFRAC,XRTR)
C
C Prints peak list
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER OUNIT, NPEAK, MPEAK, NEXT(0:MPEAK)
      DOUBLE PRECISION LHPEAK(*), LXPEAK(*), LYPEAK(*), LZPEAK(*)
      DOUBLE PRECISION LHGRID(*), LXGRID(*), LYGRID(*), LZGRID(*)
      INTEGER ORDR(*)
      LOGICAL QSYMB,QFRAC
      DOUBLE PRECISION XRTR(3,3)
C local
      INTEGER I
      DOUBLE PRECISION TEMP
      DOUBLE PRECISION LXPF, LYPF, LZPF
      DOUBLE PRECISION LXGF, LYGF, LZGF
C dummy
      DOUBLE PRECISION DBPREC
      DOUBLE COMPLEX DBCOMP
C begin
      IF (.NOT. QFRAC) THEN
      WRITE(OUNIT,'(2A)') ' peak no.  interpolated (x,y,z,height)',
     &   '     grid point (x,y,z,rho) (cartesian coordinates)'
      ELSE
      WRITE(OUNIT,'(2A)') ' peak no.  interpolated (x,y,z,height)',
     &   '     grid point (x,y,z,rho) (fractional coordinates)'
      ENDIF
C
      DO I=1,NPEAK
      IF (.NOT. QFRAC) THEN
C Print cartesian coordinates
      WRITE(OUNIT,'(A,I6,3F9.2,F11.4,A,3F9.2,F11.4)')
     & ' ',I,LXPEAK(ORDR(I)),LYPEAK(ORDR(I)),LZPEAK(ORDR(I)),
     &       LHPEAK(ORDR(I)),
     & ' ',  LXGRID(ORDR(I)),LYGRID(ORDR(I)),LZGRID(ORDR(I)),
     &       LHGRID(ORDR(I))
      ELSE
C Print fractional coordinates
      CALL TRCOOR(XRTR,
     &            LXPEAK(ORDR(I)),LYPEAK(ORDR(I)),LZPEAK(ORDR(I)),
     &            LXPF, LYPF, LZPF)
      CALL TRCOOR(XRTR,
     &            LXGRID(ORDR(I)),LYGRID(ORDR(I)),LZGRID(ORDR(I)),
     &            LXGF, LYGF, LZGF)
      WRITE(OUNIT,'(A,I6,3F9.4,F11.4,A,3F9.4,F11.4)')
     & ' ',I,LXPF,LYPF,LZPF,
     &       LHPEAK(ORDR(I)),
     & ' ',  LXGF,LYGF,LZGF,
     &       LHGRID(ORDR(I))
      ENDIF
C
      IF (QSYMB) THEN
      CALL ENCODI(I,WD,WDMAX,WDLEN)
      DBPREC=LHPEAK(ORDR(I))
      CALL DECLAR( 'PEAK_HEIGHT_'//WD(1:WDLEN),
     &             'DP', ' ', DBCOMP, DBPREC )
      DBPREC=LXPEAK(ORDR(I))
      CALL DECLAR( 'PEAK_X_'//WD(1:WDLEN),
     &             'DP', ' ', DBCOMP, DBPREC )
      DBPREC=LYPEAK(ORDR(I))
      CALL DECLAR( 'PEAK_Y_'//WD(1:WDLEN),
     &             'DP', ' ', DBCOMP, DBPREC )
      DBPREC=LZPEAK(ORDR(I))
      CALL DECLAR( 'PEAK_Z_'//WD(1:WDLEN),
     &             'DP', ' ', DBCOMP, DBPREC )
      END IF
C
      END DO
C
C store info about highest peak in symbols
      CALL DECLAR( 'PEAK_HEIGHT', 'DP', ' ', DBCOMP,LHPEAK(ORDR(1)))
      CALL DECLAR( 'PEAK_X', 'DP', ' ', DBCOMP,LXPEAK(ORDR(1)))
      CALL DECLAR( 'PEAK_Y', 'DP', ' ', DBCOMP,LYPEAK(ORDR(1)))
      CALL DECLAR( 'PEAK_Z', 'DP', ' ', DBCOMP,LZPEAK(ORDR(1)))
      TEMP=NPEAK
      CALL DECLAR( 'NPEAKS', 'DP', ' ', DBCOMP,TEMP)
C
      RETURN
      END
C======================================================================
      SUBROUTINE XPEAKAT(NPEAK,MPEAK,NEXT,
     &      LHPEAK,LXPEAK,LYPEAK,LZPEAK,LHGRID,LXGRID,LYGRID,
     &      LZGRID,ORDR)
C
C Creates dummy atoms.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
C
      INTEGER NPEAK, MPEAK, NEXT(0:MPEAK)
      DOUBLE PRECISION LHPEAK(*), LXPEAK(*), LYPEAK(*), LZPEAK(*)
      DOUBLE PRECISION LHGRID(*), LXGRID(*), LYGRID(*), LZGRID(*)
      INTEGER ORDR(*)
C local
      INTEGER I, LEN
C parameter
      DOUBLE PRECISION ZERO, TEN
      PARAMETER (ZERO=0.0D0, TEN=10.D0)
C begin
      IF (NGRP.EQ.0) THEN
        IGPBS(1)=0
      END IF
C
      DO I=1,NPEAK
C
      IF (NATOM.GE.MAXA) THEN
      CALL WRNDIE(-5,'PEAKPIK',
     & 'exceeded MAXA parameter --> recompile program')
      ELSE
      NATOM=NATOM+1
      END IF
C
C
C set atom properties
      IBLO(NATOM)=NNB
      IAC(NATOM)='PEAK'
      AMASS(NATOM)=TEN
      TYPE(NATOM)='PEAK'
      SEGID(NATOM)='PEAK'
      RESID(NATOM)=' '
      CALL ENCODI(I,RESID(NATOM),4,LEN)
      RES(NATOM)='PEAK'
      CG(NATOM)=ZERO
      LOOKUP(NATOM)=0
      QLOOKU(NATOM)=.TRUE.
      IMOVE(NATOM)=0
      CALL ATMINI(NATOM,NATOM)
      X(NATOM)=LXPEAK(ORDR(I))
      Y(NATOM)=LYPEAK(ORDR(I))
      Z(NATOM)=LZPEAK(ORDR(I))
      WMAIN(NATOM)=LHPEAK(ORDR(I))
      NGRP=NGRP+1
      IGPBS(NGRP+1)=NATOM
      END DO
C
C store info about highest peak in symbols
      WRITE(6,'(A,I6,A)')
     & ' PEAKpik: ',NPEAK,
     & ' atoms were appended to molecular structure.'
C
      CALL SCRATC
C
      IF (XRUPAT) THEN
      XRUPAT=.FALSE.
      CALL XRAATM(NATOM)
C
C initialize the default associate selection
C ASSociate FCALC ( not hydrogen )
C
      IHPFLAG=1
      HPFLAG(1)=ALLHP(INTEG4(XRMATO))
      XASSOC(1)='FCALC'
      CALL XRINAT(HEAP(HPFLAG(1)),HEAP(HPANOMFLAG),NATOM)
      END IF
C
      RETURN
      END
