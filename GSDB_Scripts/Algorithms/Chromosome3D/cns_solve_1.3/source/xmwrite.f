      SUBROUTINE XMAPX(MODE,
     &           NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &           MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,QHERM,XRCELL,MAPR,
     &           XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &           HPRHOMA,NRHO,IRHO,NMASK,XRMAP,
     &           XRTR,XRINTR,XRVOL,
     &           XRMREF,XRNREF,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,
     &           HPTYPE,XRRED,XRSYGP,XRSYIV)
C
C Parser routine for writing maps to specified file.
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
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'xmask.inc'
      CHARACTER*(*) MODE
      INTEGER NA, NB, NC, NAP, NBP, NCP, NAPP, NBPP, NCPP
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9)
      DOUBLE PRECISION MAPR
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER HPRHOMA, IRHO, NRHO, NMASK
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
      CHARACTER*4 EXTEND
      CHARACTER*4 MAPTYPE
      DOUBLE PRECISION XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, CUSHON
      INTEGER AMIN, AMAX, BMIN, BMAX, CMIN, CMAX, IFROM, I, J
      INTEGER OUNIT, IOUNIT, OLEN, FLAGS, NFLAGS, ILEN
      LOGICAL COND, QAUTO, QFORM, ERR
      INTEGER LNA, LNB, LNC
      INTEGER LMAASY, LMBASY, LMCASY, LNAASY, LNBASY, LNCASY
      INTEGER IMFROM
      DOUBLE PRECISION TMPMAT(MXMAVG,3,4), CENTER(3)
C pointers
      INTEGER SECTON, ISECTON, HPRMAP, HPIMAP, HPMASK, TMPMAP, ITMPMAP
      INTEGER HPOMASK, GRIDRT
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
      REAL SZERO
      PARAMETER (SZERO=0.0)
C begin
C default values
      ERR=.FALSE.
      IFROM=1
      OFILE=' '
      IFILE=' '
      OUNIT=0
      IOUNIT=0
      EXTEND='MOLE'
      MAPTYPE='CNS'
      XMIN=ZERO
      XMAX=ONE
      YMIN=ZERO
      YMAX=ONE
      ZMIN=ZERO
      ZMAX=ONE
      CUSHON=TWO
      QAUTO=.TRUE.
      QFORM=.TRUE.
      FLAGS=ALLHP(INTEG4(NATOM))
      CALL FILL4(HEAP(FLAGS),NATOM,1)
      NFLAGS=NATOM
      IMFROM=0
C
      CALL PUSEND('WRITe MAP>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('WRITe MAP>')
      CALL MISCOM('WRITe MAP>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-write-map')
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
      WRITE(6,'(3A)') ' %XMAPX-ERR: real space object ',
     & WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'XMAPX','object undeclared.')
      ERR=.TRUE.
      END IF
      END IF
C ========================================================================
      ELSE IF (WD(1:4).EQ.'MASK') THEN
      CALL NEXTWD('MASK=')
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD('FROM=')
C
      COND=.FALSE.
      DO I=1,XRHONUM
      IF (WD(1:WDLEN).EQ.XRHONAM(I)) THEN
      COND=.TRUE.
      IMFROM=I
      END IF
      END DO
      IF (.NOT.COND) THEN
      WRITE(6,'(3A)') ' %XMAPX-ERR: real space object ',
     & WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'XMAPX','object undeclared.')
      ERR=.TRUE.
      END IF
C
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
      ELSE IF (WD(1:4).EQ.'IOUT') THEN
      CALL NEXTFI('IOUTput=',IFILE)
      ELSE IF (WD(1:4).EQ.'EXTE') THEN
      CALL NEXTA4('EXTEnd=',EXTEND)
      COND= EXTEND.EQ.'UNIT'.OR.
     &      EXTEND.EQ.'MOLE'.OR.
     &      EXTEND.EQ.'MASK'.OR.
     &      EXTEND.EQ.'BOX '.OR.EXTEND.EQ.'ORTH'.OR.EXTEND.EQ.'FRAC'.OR.
     &      EXTEND.EQ.'ASYM'
      IF (.NOT.COND) CALL DSPERR('XMAP','unknown value for EXTEnd')
      ELSE IF (WD(1:4).EQ.'TYPE') THEN
      CALL NEXTA4('TYPE=',MAPTYPE)
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(HEAP(FLAGS),NFLAGS,X,Y,Z,.TRUE.)
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
      ELSE IF (WD(1:4).EQ.'CUSH') THEN
      CALL NEXTF('CUSHion=',CUSHON)
      ELSE IF (WD(1:4).EQ.'AUTO') THEN
      CALL NEXTLO('AUTOmatic-scale=',QAUTO)
      ELSE IF (WD(1:4).EQ.'FORM') THEN
      CALL NEXTLO('FORMatted=',QFORM)
C
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' ---------------------------',
     & '-----MAP-parameters---------------------------------'
      OLEN=LEN(OFILE)
      CALL TRIMM(OFILE,OLEN)
      ILEN=LEN(IFILE)
      CALL TRIMM(IFILE,ILEN)
      WRITE(6,'(5A,L1,A,L1,A,A4,/,4A)') ' | OUTPut-file=',
     &    OFILE(1:OLEN),' EXTEnd=',EXTEND,' AUTOmatic-scale=',QAUTO,
     &    ' FORMatted=',QFORM,' TYPE= ',MAPTYPE,
     &    ' | FROM=',XRHONAM(IFROM),'    IOUTput-file=',IFILE(1:ILEN)
      IF (EXTEND.EQ.'BOX '.OR.EXTEND.EQ.'ORTH'.OR.EXTEND.EQ.'FRAC') THEN
      WRITE(6,'(6(A,F9.4))')
     & ' | XMIN=',XMIN,' XMAX=',XMAX,' YMIN=',YMIN,' YMAX=',YMAX,
     & ' ZMIN=',ZMIN,' ZMAX=',ZMAX
      ELSE IF (EXTEND.EQ.'MOLE') THEN
      WRITE(6,'(A,F6.2,A,I6,A)') ' | CUSHion=',CUSHON,
     &    ', map covers ',NFLAGS,' atoms'
      ELSE IF (EXTEND.EQ.'MASK') THEN
      WRITE(6,'(A,F6.2,A,A)') ' | CUSHion=',CUSHON,
     &    ', map covers mask ',XRHONAM(IMFROM)
      END IF
      WRITE(6,'(2A)') ' ---------------------------',
     & '----------------------------------------------------'
C=====================================================================
      ELSE
      CALL CHKEND('WRITe MAP>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (.NOT.ERR) THEN
C
      IF (IFILE.NE.' '.AND.QHERM) THEN
      CALL WRNDIE(-5,'XMAP',
     & 'Cannot write imaginary density (density is real)')
      IFILE=' '
      END IF
C
      IF (QFORM) THEN
      IF (OFILE.NE.' ') THEN
      CALL ASSFIL(OFILE,OUNIT,'WRITE','FORMATTED',ERROR)
      END IF
      IF (IFILE.NE.' ') THEN
      CALL ASSFIL(IFILE,IOUNIT,'WRITE','FORMATTED',ERROR)
      END IF
      ELSE
      IF (OFILE.NE.' ') THEN
      CALL ASSFIL(OFILE,OUNIT,'WRITE','UNFORMATTED',ERROR)
      END IF
      IF (IFILE.NE.' ') THEN
      CALL ASSFIL(IFILE,IOUNIT,'WRITE','UNFORMATTED',ERROR)
      END IF
      END IF
C
      IF (QAUTO) THEN
      WRITE(6,'(A)') ' XMAPX: map will be scaled.'
      ELSE
      WRITE(6,'(A)') ' XMAPX: map will be written in [e]/A^3 units.'
      END IF
C
C get the map pointers
      IF (HPRRHO(IFROM).EQ.0) THEN
      WRITE(6,'(3A)') ' %XMDOEVA-ERR: real space object ',
     & XRHONAM(IFROM),' undefined.'
      CALL WRNDIE(-5,'XMDOEVA','object undefined.')
      ERR=.TRUE.
      ELSE
      HPRMAP=HPRRHO(IFROM)
      HPIMAP=HPIRHO(IFROM)
      END IF
C
C get the mask pointer if needed
      IF (EXTEND.EQ.'MASK') THEN
      IF (HPRRHO(IMFROM).EQ.0) THEN
      WRITE(6,'(3A)') ' %XMDOEVA-ERR: real space object ',
     & XRHONAM(IMFROM),' undefined.'
      CALL WRNDIE(-5,'XMDOEVA','object undefined.')
      ERR=.TRUE.
      ELSE
      HPOMASK=HPRRHO(IMFROM)
      END IF
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
      IF (OUNIT.EQ.0.AND.IOUNIT.EQ.0) THEN
      CALL WRNDIE(-5,'XMAPX',' no output file specified.')
      ELSE IF (EXTEND.EQ.'MOLE'.AND.NFLAGS.EQ.0) THEN
      CALL WRNDIE(-5,'XMAPX',' zero atom selection.')
      ELSE IF (.NOT.ERR) THEN
C
C determine extend of map
      IF (EXTEND.EQ.'ASYM') THEN
C
C use asymmetric unit
      AMIN=LMAASY
      AMAX=LNAASY
      BMIN=LMBASY
      BMAX=LNBASY
      CMIN=LMCASY
      CMAX=LNCASY
C
      ELSE IF (EXTEND.EQ.'MASK') THEN
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
     &     HEAP(HPMASK),HEAP(HPOMASK),HEAP(GRIDRT),
     &     TMPMAT,CENTER,CUSHON,AMIN,AMAX,BMIN,BMAX,CMIN,CMAX)
C
      CALL FREHP(GRIDRT,IREAL8(XRNSYM*3*4))
C
      ELSE
      CALL XMAPMX(EXTEND,CUSHON,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,
     &            HEAP(FLAGS),XRTR,XRCELL)
      AMIN=NINT(XMIN*LNA)
      AMAX=NINT(XMAX*LNA)
      BMIN=NINT(YMIN*LNB)
      BMAX=NINT(YMAX*LNB)
      CMIN=NINT(ZMIN*LNC)
      CMAX=NINT(ZMAX*LNC)
      END IF
C
C write sectioning information
      WRITE(6,'(A,3I4,A,3I4,A,3I4,A)')
     & ' XMAPX: extend NA=(',LNA,AMIN,AMAX,') NB=(',LNB,BMIN,BMAX,
     & ') NC=(',LNC,CMIN,CMAX,') '
      IF (MAPTYPE(1:3).EQ.'CNS') THEN
C
C allocate temporary heap space for z-sections
      SECTON=ALLHP(IREAL8((AMAX-AMIN+1)*(BMAX-BMIN+1)))
      ISECTON=ALLHP(IREAL8((AMAX-AMIN+1)*(BMAX-BMIN+1)))
C
      CALL XMAPX2(LNA,LNB,LNC,
     &            LMAASY,LMBASY,LMCASY,LNAASY,LNBASY,LNCASY,
     &            QHERM,HEAP(HPRMAP),HEAP(HPIMAP),HEAP(HPMASK),
     &            XRCELL,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &            AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,
     &            HEAP(SECTON),HEAP(ISECTON),
     &            QAUTO,QFORM,OUNIT,IOUNIT,EXTEND,MODE)
      IF (OUNIT.GT.0) CALL VCLOSE(OUNIT,'KEEP',ERROR)
      IF (IOUNIT.GT.0) CALL VCLOSE(IOUNIT,'KEEP',ERROR)
C
C
C deallocate temporary heap space for z-sections
      CALL FREHP(ISECTON,IREAL8((AMAX-AMIN+1)*(BMAX-BMIN+1)))
      CALL FREHP(SECTON,IREAL8((AMAX-AMIN+1)*(BMAX-BMIN+1)))
C
      ELSEIF (MAPTYPE(1:3).EQ.'EZD') THEN
C
C allocate temporary heap space for
      TMPMAP=ALLHP(IREAL8((AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
      ITMPMAP=ALLHP(IREAL8((AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
C
      CALL XEZDMAP(LNA,LNB,LNC,
     &            LMAASY,LMBASY,LMCASY,LNAASY,LNBASY,LNCASY,
     &            QHERM,HEAP(HPRMAP),HEAP(HPIMAP),HEAP(HPMASK),
     &            XRCELL,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &            AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,
     &            HEAP(TMPMAP),HEAP(ITMPMAP),
     &            QAUTO,QFORM,OUNIT,IOUNIT,EXTEND,MODE)
C
C deallocate temporary heap space for
      CALL FREHP(TMPMAP,
     &           IREAL8((AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
      CALL FREHP(ITMPMAP,
     &           IREAL8((AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
C
      END IF
C
      END IF
      END IF
C
      CALL FREHP(FLAGS,INTEG4(NATOM))
      RETURN
      END
C======================================================================
      SUBROUTINE  XMAPX2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &            QHERM,RRHO,IRHO,RHOMASK,XRCELL,
     &            XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &            AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,SECTON,ISECTON,
     &            QAUTO,QFORM,OUNIT,IOUNIT,EXTEND,MODE)
C
C Write density map in CNS format
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      LOGICAL QHERM
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4)
      DOUBLE PRECISION XRCELL(*)
      INTEGER AMIN, AMAX, BMIN, BMAX, CMIN, CMAX
      DOUBLE PRECISION SECTON(AMAX-AMIN+1,BMAX-BMIN+1)
      DOUBLE PRECISION ISECTON(AMAX-AMIN+1,BMAX-BMIN+1)
      LOGICAL QAUTO, QFORM
      INTEGER OUNIT, IOUNIT
      CHARACTER*(*) EXTEND, MODE
C local
      INTEGER A, B, C, NR, SYM, LSYM
      INTEGER KSECT, ISECT, JSECT, I, J
      INTEGER AA1, AA2, AA3, BB1, BB2, BB3, AS, BS, CS
      DOUBLE PRECISION RAVE, RSIGMA, IRAVE, IRSIGMA, RAVE0, RSIGMA0
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
C real part first:
      IF (OUNIT.GT.0) THEN
C get average and sigma of map
      RAVE=ZERO
      RSIGMA=ZERO
      NR=0
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).NE.0) THEN
      RAVE=RAVE+RRHO(A,B,C)/RHOMASK(A,B,C)
      RSIGMA=RSIGMA+(RRHO(A,B,C)/RHOMASK(A,B,C))**2
      NR=NR+1
      END IF
      END DO
      END DO
      END DO
C
C print overall statistics
      RAVE=RAVE/NR
      RSIGMA=SQRT(MAX(ZERO,RSIGMA/NR-RAVE**2))
      WRITE(6,'(A,F14.5,A,F14.5,A)')
     & ' XMAPX: ave. density (real) in unit cell=', RAVE,
     & ' sigma=',RSIGMA,' e/A^3'
      IF (QAUTO.AND.RSIGMA.LT.R4SMAL) THEN
      WRITE(6,'(A)')
     & ' XMAPX: sigma very small for map.  No scaling performed.'
      RSIGMA=ONE
      END IF
C
      END IF
C
C now the imaginary part:
      IF (IOUNIT.GT.0) THEN
C get average and sigma of map
      IRAVE=ZERO
      IRSIGMA=ZERO
      NR=0
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).NE.0) THEN
      IRAVE=IRAVE+IRHO(A,B,C)/RHOMASK(A,B,C)
      IRSIGMA=IRSIGMA+(IRHO(A,B,C)/RHOMASK(A,B,C))**2
      NR=NR+1
      END IF
      END DO
      END DO
      END DO
C
C print overall statistics
      IRAVE=IRAVE/NR
      IRSIGMA=SQRT(MAX(ZERO,IRSIGMA/NR-IRAVE**2))
      WRITE(6,'(A,F14.5,A,F14.5,A)')
     & ' XMAPX: ave. density (imaginary) in unit cell=', IRAVE,
     & ' sigma=',IRSIGMA,' e/A^3'
      IF (QAUTO.AND.IRSIGMA.LT.R4SMAL) THEN
      WRITE(6,'(A)')
     & ' XMAPX: sigma very small for map.  No scaling performed.'
      IRSIGMA=ONE
      END IF
C
      END IF
C
C write map headers
      IF (OUNIT.NE.0) THEN
      IF (QFORM) THEN
      CALL WRTITL(OUNIT,2)
      WRITE(OUNIT,'(9I8)')
     &  NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX
      WRITE(OUNIT,'(6E12.5)') (XRCELL(I),I=1,6)
      WRITE(OUNIT,'(A)') 'ZYX'
      ELSE
      CALL WRTITL(OUNIT,-1)
      WRITE(OUNIT)
     &  NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX
      WRITE(OUNIT) (XRCELL(I),I=1,6)
      WRITE(OUNIT) 'ZYX'
      END IF
      END IF
C
      IF (IOUNIT.NE.0) THEN
      IF (QFORM) THEN
      CALL WRTITL(IOUNIT,2)
      WRITE(IOUNIT,'(9I8)')
     &  NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX
      WRITE(IOUNIT,'(6E12.5)') (XRCELL(I),I=1,6)
      WRITE(IOUNIT,'(A)') 'ZYX'
      ELSE
      CALL WRTITL(IOUNIT,-1)
      WRITE(IOUNIT)
     &  NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX
      WRITE(IOUNIT) (XRCELL(I),I=1,6)
      WRITE(IOUNIT) 'ZYX'
      END IF
      END IF
C
C write density matrix, c is slowest ("z-sections").
C
      KSECT=-1
      DO C=CMIN,CMAX
      KSECT=KSECT+1
C
C loop over all symmetry operators unless we're just
C writing the asymmetric unit or we're in MAP mode (old syntax)
      IF (EXTEND.EQ.'ASYM') THEN
      LSYM=1
      DO I=1,AMAX-AMIN+1
      DO J=1,BMAX-BMIN+1
      IF (OUNIT.GT.0) SECTON(I,J)=ZERO
      IF (IOUNIT.GT.0) ISECTON(I,J)=ZERO
      END DO
      END DO
      ELSE
      LSYM=XRNSYM
      END IF
C
      DO SYM=1,LSYM
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
      IF (AS.LE.NAASY.AND.BS.LE.NBASY.AND.CS.LE.NCASY) THEN
      IF (RHOMASK(AS,BS,CS).GT.0) THEN
      IF (OUNIT.GT.0) SECTON(ISECT,JSECT)=RRHO(AS,BS,CS)
      IF (IOUNIT.GT.0) ISECTON(ISECT,JSECT)=IRHO(AS,BS,CS)
      END IF
      END IF
      END DO
      END DO
      END DO
C
C do the real section first
      IF (OUNIT.GT.0) THEN
C do some statistics
      RAVE0=ZERO
      RSIGMA0=ZERO
      NR=0
      DO I=1,ISECT
      DO J=1,JSECT
      NR=NR+1
      RAVE0=RAVE0+SECTON(I,J)
      RSIGMA0=RSIGMA0+SECTON(I,J)**2
      END DO
      END DO
C
C print info about the section
      RAVE0=RAVE0/NR
      RSIGMA0=SQRT(MAX(ZERO,RSIGMA0/NR-RAVE0**2))
      WRITE(6,'(A,I4,A,F14.5,A,F14.5,A)')
     & ' XMAPX: real part of section #',KSECT,' ave. dens.=',RAVE0,
     & ' sigma=',RSIGMA0,' e/A^3'
C
C scale the section
      IF (QAUTO) THEN
      DO I=1,ISECT
      DO J=1,JSECT
      SECTON(I,J)=(SECTON(I,J)-RAVE)/RSIGMA
      END DO
      END DO
      END IF
C
C write the section
      IF (QFORM) THEN
      WRITE(OUNIT,'(I8)') KSECT
      WRITE(OUNIT,'(6E12.5)')
     &  ((SECTON(I,J),I=1,ISECT),J=1,JSECT)
      ELSE
      WRITE(OUNIT) KSECT
      WRITE(OUNIT)((SECTON(I,J),I=1,ISECT),J=1,JSECT)
      END IF
C
      END IF
C
C
C now do the imaginary section
      IF (IOUNIT.GT.0) THEN
C do some statistics
      RAVE0=ZERO
      RSIGMA0=ZERO
      NR=0
      DO I=1,ISECT
      DO J=1,JSECT
      NR=NR+1
      RAVE0=RAVE0+ISECTON(I,J)
      RSIGMA0=RSIGMA0+ISECTON(I,J)**2
      END DO
      END DO
C
C print info about the section
      RAVE0=RAVE0/NR
      RSIGMA0=SQRT(MAX(ZERO,RSIGMA0/NR-RAVE0**2))
      WRITE(6,'(A,I4,A,F14.5,A,F14.5,A)')
     & ' XMAPX: imaginary part of section #',KSECT,
     & ' ave. dens.=',RAVE0,' sigma=',RSIGMA0,' e/A^3'
C
C scale the section
      IF (QAUTO) THEN
      DO I=1,ISECT
      DO J=1,JSECT
      ISECTON(I,J)=(ISECTON(I,J)-IRAVE)/IRSIGMA
      END DO
      END DO
      END IF
C
C write the section
      IF (QFORM) THEN
      WRITE(IOUNIT,'(I8)') KSECT
      WRITE(IOUNIT,'(6E12.5)')
     &  ((ISECTON(I,J),I=1,ISECT),J=1,JSECT)
      ELSE
      WRITE(IOUNIT) KSECT
      WRITE(IOUNIT)((ISECTON(I,J),I=1,ISECT),J=1,JSECT)
      END IF
C
      END IF
C
      END DO
C
C terminate output file with negative integer number
C also write map scaling parameters.
      IF (OUNIT.GT.0) THEN
      IF (QFORM) THEN
      WRITE(OUNIT,'(I8)') -9999
      IF (QAUTO) THEN
      WRITE(OUNIT,'(2(E12.4,1X))') RAVE,RSIGMA
      ELSE
      WRITE(OUNIT,'(2(E12.4,1X))') ZERO,ONE
      END IF
      ELSE
      WRITE(OUNIT) -9999
      IF (QAUTO) THEN
      WRITE(OUNIT) RAVE,RSIGMA
      ELSE
      WRITE(OUNIT) ZERO,ONE
      END IF
      END IF
      END IF
C
C terminate output file with negative integer number
C also write map scaling parameters.
      IF (IOUNIT.GT.0) THEN
      IF (QFORM) THEN
      WRITE(IOUNIT,'(I8)') -9999
      IF (QAUTO) THEN
      WRITE(IOUNIT,'(2(E12.4,1X))') IRAVE, IRSIGMA
      ELSE
      WRITE(IOUNIT,'(2(E12.4,1X))') ZERO, ONE
      END IF
      ELSE
      WRITE(IOUNIT) -9999
      IF (QAUTO) THEN
      WRITE(IOUNIT) IRAVE, IRSIGMA
      ELSE
      WRITE(IOUNIT) ZERO, ONE
      END IF
      END IF
      END IF
C
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4,/,A,F10.4)')
     &        ' XMAPX: CPU-time: expansion=',SECS2-SECS
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XEZDMAP(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &           QHERM,RRHO,IRHO,RHOMASK,XRCELL,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,TMPMAP,ITMPMAP,
     &           QAUTO,QFORM,OUNIT,IOUNIT,EXTEND,MODE)
C
C Write density map in EZD format
C
C Author: Paul Adams and Axel T. Brunger
C ======================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      LOGICAL QHERM
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL TMPMAP(*)
      REAL ITMPMAP(*)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4)
      DOUBLE PRECISION XRCELL(*)
      INTEGER AMIN, AMAX, BMIN, BMAX, CMIN, CMAX
      LOGICAL QAUTO, QFORM
      INTEGER OUNIT, IOUNIT
      CHARACTER*(*) EXTEND, MODE
C local
      INTEGER A, B, C, NR, SYM, LSYM
      INTEGER I, J, I2
      INTEGER AA1, AA2, AA3, BB1, BB2, BB3, AS, BS, CS
      INTEGER NPOINT, POINT
      DOUBLE PRECISION RAVE, RSIGMA, IRAVE, IRSIGMA
      DOUBLE PRECISION SECS, SECS2, RSCALE, ISCALE
      INTEGER LINELEN, MXLINE
      PARAMETER (MXLINE=84)
      CHARACTER*(MXLINE) LINE
C parameters
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
C begin
C
C number of points in map
      NPOINT=(AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)
C
C initialize timers
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
C real part first:
      IF (OUNIT.GT.0) THEN
C get average and sigma of map
      RAVE=ZERO
      RSIGMA=ZERO
      NR=0
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).NE.0) THEN
      RAVE=RAVE+RRHO(A,B,C)/RHOMASK(A,B,C)
      RSIGMA=RSIGMA+(RRHO(A,B,C)/RHOMASK(A,B,C))**2
      NR=NR+1
      END IF
      END DO
      END DO
      END DO
C
C print overall statistics
      RAVE=RAVE/NR
      RSIGMA=SQRT(MAX(ZERO,RSIGMA/NR-RAVE**2))
      WRITE(6,'(A,F14.5,A,F14.5,A)')
     & ' XEZDMAP: ave. density (real) in unit cell=', RAVE,
     & ' sigma=',RSIGMA,' e/A^3'
      IF (QAUTO) THEN
      IF (RSIGMA.LT.R4SMAL) THEN
      WRITE(6,'(A)')
     & ' XEZDMAP: sigma very small for map.  No scaling performed.'
      RSIGMA=ONE
      RSCALE=ONE
      ELSE
      RSCALE=ONE/RSIGMA
      END IF
      ELSE
      RSCALE=ONE
      END IF
C
      END IF
C
C now the imaginary part:
      IF (IOUNIT.GT.0) THEN
C get average and sigma of map
      IRAVE=ZERO
      IRSIGMA=ZERO
      NR=0
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).NE.0) THEN
      IRAVE=IRAVE+IRHO(A,B,C)/RHOMASK(A,B,C)
      IRSIGMA=IRSIGMA+(IRHO(A,B,C)/RHOMASK(A,B,C))**2
      NR=NR+1
      END IF
      END DO
      END DO
      END DO
C
C print overall statistics
      IRAVE=IRAVE/NR
      IRSIGMA=SQRT(MAX(ZERO,IRSIGMA/NR-IRAVE**2))
      WRITE(6,'(A,F14.5,A,F14.5,A)')
     & ' XEZDMAP: ave. density (imaginary) in unit cell=', IRAVE,
     & ' sigma=',IRSIGMA,' e/A^3'
      IF (QAUTO) THEN
      IF (IRSIGMA.LT.R4SMAL) THEN
      WRITE(6,'(A)')
     & ' XEZDMAP: sigma very small for map.  No scaling performed.'
      IRSIGMA=ONE
      ISCALE=ONE
      ELSE
      ISCALE=ONE/IRSIGMA
      END IF
      ELSE
      ISCALE=ONE
      END IF
C
      END IF
C
C write map headers
      IF (OUNIT.NE.0) THEN
      WRITE(OUNIT,'(A)')'EZD_MAP'
      WRITE(OUNIT,'(A,6F10.5)') 'CELL ',(XRCELL(I),I=1,6)
      WRITE(OUNIT,'(A,3I6)')
     &  'ORIGIN ',AMIN,BMIN,CMIN
      WRITE(OUNIT,'(A,3I6)')
     &  'EXTENT ',AMAX-AMIN+1,BMAX-BMIN+1,CMAX-CMIN+1
      WRITE(OUNIT,'(A,3I6)')
     &  'GRID ',NA,NB,NC
      WRITE(OUNIT,'(A,F10.5)') 'SCALE ',RSCALE
      WRITE(OUNIT,'(A)') 'MAP'
      END IF
C
      IF (IOUNIT.NE.0) THEN
      WRITE(IOUNIT,'(A)')'EZD_MAP'
      WRITE(IOUNIT,'(A,6F10.5)') 'CELL ',(XRCELL(I),I=1,6)
      WRITE(IOUNIT,'(A,3I6)')
     &  'ORIGIN ',AMIN,BMIN,CMIN
      WRITE(IOUNIT,'(A,3I6)')
     &  'EXTENT ',AMAX-AMIN+1,BMAX-BMIN+1,CMAX-CMIN+1
      WRITE(IOUNIT,'(A,3I6)')
     &  'GRID ',NA,NB,NC
      WRITE(IOUNIT,'(A,F10.5)') 'SCALE ',ISCALE
      WRITE(IOUNIT,'(A)') 'MAP'
      END IF
C
C loop over all symmetry operators unless we're just
C writing the asymmetric unit or we're in MAP mode (old syntax)
      IF (EXTEND.EQ.'ASYM') THEN
      LSYM=1
      DO I=1,NPOINT
      IF (OUNIT.GT.0) TMPMAP(I)=ZERO
      IF (IOUNIT.GT.0) ITMPMAP(I)=ZERO
      END DO
      ELSE
      LSYM=XRNSYM
      END IF
C
      DO SYM=1,LSYM
C
      POINT=1
      DO C=CMIN,CMAX
C
      AA1=XRSYMM(SYM,1,3)*C+(XRSYMM(SYM,1,4)*NA)/XRSYTH
      AA2=XRSYMM(SYM,2,3)*C+(XRSYMM(SYM,2,4)*NB)/XRSYTH
      AA3=XRSYMM(SYM,3,3)*C+(XRSYMM(SYM,3,4)*NC)/XRSYTH
C
      DO B=BMIN,BMAX
      BB1=AA1+XRSYMM(SYM,1,2)*B+10000*NA-MAASY
      BB2=AA2+XRSYMM(SYM,2,2)*B+10000*NB-MBASY
      BB3=AA3+XRSYMM(SYM,3,2)*B+10000*NC-MCASY
C
      DO A=AMIN,AMAX
C compute symmetry-related indices and project into ASU brick
      AS=MOD(XRSYMM(SYM,1,1)*A+BB1,NA) +MAASY
      BS=MOD(XRSYMM(SYM,2,1)*A+BB2,NB) +MBASY
      CS=MOD(XRSYMM(SYM,3,1)*A+BB3,NC) +MCASY
      IF (AS.LE.NAASY.AND.BS.LE.NBASY.AND.CS.LE.NCASY) THEN
      IF (RHOMASK(AS,BS,CS).GT.0) THEN
      IF (OUNIT.GT.0) TMPMAP(POINT)=RRHO(AS,BS,CS)
      IF (IOUNIT.GT.0) ITMPMAP(POINT)=IRHO(AS,BS,CS)
      END IF
      END IF
      POINT=POINT+1
      END DO
      END DO
      END DO
      END DO
C
C do the real section first
      IF (OUNIT.GT.0) THEN
C
C scale the map
      IF (QAUTO) THEN
         DO I=1,NPOINT
            TMPMAP(I)=(TMPMAP(I)-RAVE)/RSIGMA
         END DO
      END IF
C
C write the map
      DO I=1,NPOINT,7
         I2=MIN(NPOINT,I+7-1)
         WRITE(LINE,'(7F12.5)') (TMPMAP(J),J=I,I2)
         LINELEN=MXLINE
         CALL RMSPACE(LINE,LINELEN)
         WRITE(OUNIT,'(A)') LINE(1:LINELEN)
      END DO
C
      END IF
C
C now do the imaginary section
      IF (IOUNIT.GT.0) THEN
C
C scale the section
      IF (QAUTO) THEN
         DO I=1,NPOINT
            ITMPMAP(I)=(ITMPMAP(I)-IRAVE)/IRSIGMA
         END DO
      END IF
C
C write the map
      DO I=1,NPOINT,7
         I2=MIN(NPOINT,I+7-1)
         WRITE(LINE,'(7F12.5)') (TMPMAP(J),J=I,I2)
         LINELEN=MXLINE
         CALL RMSPACE(LINE,LINELEN)
         WRITE(IOUNIT,'(A)') LINE(1:LINELEN)
      END DO
C
      END IF
C
      IF (OUNIT.GT.0) THEN
      WRITE(OUNIT,'(A)') 'END'
      WRITE(6,'(A,I10,A,A)')' XEZDMAP: ',NPOINT,' grid points ',
     &'written to EZD map'
      END IF
C
      IF (IOUNIT.GT.0) THEN
      WRITE(IOUNIT,'(A)') 'END'
      WRITE(6,'(A,I10,A,A)')' XEZDMAP: ',NPOINT,' grid points ',
     &'written to EZD map'
      END IF
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4,/,A,F10.4)')
     &        ' XEZDMAP: CPU-time: expansion=',SECS2-SECS
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMAPMX(EXTEND,CUSHON,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,
     &                  FLAGS,XRTR,XRCELL)
C
C determine max and min values, convert into fractional coordinates
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      CHARACTER*4 EXTEND
      DOUBLE PRECISION CUSHON, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
      INTEGER FLAGS(*)
      DOUBLE PRECISION XRTR(3,3), XRCELL(*)
C local
      INTEGER IAT, COUNT
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
      ELSE IF (EXTEND.EQ.'BOX '.OR.EXTEND.EQ.'ORTH') THEN
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
      ELSE IF (EXTEND.EQ.'FRAC') THEN
C
C use fractional coordinates -- just use XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
      CONTINUE
C
      ELSE IF (EXTEND.EQ.'MOLE') THEN
C
C determine min and max
      XMIN=+9999.D0
      XMAX=-9999.D0
      YMIN=+9999.D0
      YMAX=-9999.D0
      ZMIN=+9999.D0
      ZMAX=-9999.D0
      COUNT=0
C
C go through all SELECTED atoms
      DO IAT=1,NATOM
      IF (FLAGS(IAT).EQ.1.AND.INITIA(IAT,X,Y,Z)) THEN
      COUNT=COUNT+1
C
C fractionalize coordinates of selected atoms
      TEMPX=XRTR(1,1)*X(IAT) +XRTR(1,2)*Y(IAT) +XRTR(1,3)*Z(IAT)
      TEMPY=XRTR(2,1)*X(IAT) +XRTR(2,2)*Y(IAT) +XRTR(2,3)*Z(IAT)
      TEMPZ=XRTR(3,1)*X(IAT) +XRTR(3,2)*Y(IAT) +XRTR(3,3)*Z(IAT)
      XMIN=MIN(XMIN,TEMPX)
      XMAX=MAX(XMAX,TEMPX)
      YMIN=MIN(YMIN,TEMPY)
      YMAX=MAX(YMAX,TEMPY)
      ZMIN=MIN(ZMIN,TEMPZ)
      ZMAX=MAX(ZMAX,TEMPZ)
      END IF
      END DO
      XMIN=XMIN-CUSHON/XRCELL(1)
      XMAX=XMAX+CUSHON/XRCELL(1)
      YMIN=YMIN-CUSHON/XRCELL(2)
      YMAX=YMAX+CUSHON/XRCELL(2)
      ZMIN=ZMIN-CUSHON/XRCELL(3)
      ZMAX=ZMAX+CUSHON/XRCELL(3)
      IF (COUNT.EQ.0) THEN
            CALL WRNDIE(-1,'XMAP',
     & 'zero known atoms selected for map calculation')
      XMIN=ZERO
      XMAX=ZERO
      YMIN=ZERO
      YMAX=ZERO
      ZMIN=ZERO
      ZMAX=ZERO
      END IF
C
      END IF
C
      RETURN
      END
C
C =========================================================================
C
      SUBROUTINE RMSPACE(ST,STLEN)
C
      IMPLICIT NONE
C
      INTEGER STLEN
      CHARACTER*(*) ST
C local
      INTEGER I, LENGTH
      LOGICAL ONESPACE
C begin
      I=1
      LENGTH=STLEN
      ONESPACE=.FALSE.
      DO WHILE (I.LE.LENGTH)
         IF (ST(I:I).EQ.' ') THEN
            IF (ONESPACE) THEN
               ST(I:(LENGTH-1))=ST((I+1):LENGTH)
               LENGTH=LENGTH-1
               I=I-1
            ELSE
               ONESPACE=.TRUE.
            END IF
         ELSE
            ONESPACE=.FALSE.
         END IF
         I=I+1
      END DO
C
      STLEN=LENGTH
C
      RETURN
      END
C
C =========================================================================
C
