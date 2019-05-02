! high-level OpenMP parallelization in XFFTUP: the subgrid algorithm
! is used to hand out work to the OpenMP threads
!
! Kay Diederichs 12/2006
!
! this version has minimal changes wrt the original source:
!
! a) XFGRID was changed/simplified; 
! b) OpenMP directives inserted into XFFTUP
! c) due to a), the MEMORY subcommand of the FFT command has no effect
! d) a bit of trivial parallelization in XFGRCHK
! 
      SUBROUTINE XFFTIN
C
C routine initializes parameters for Fast Fourier computation
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xfft.inc'
      INCLUDE 'comand.inc'
      INTEGER DUMMY
C parameters
      DOUBLE PRECISION ONE, THREE, TWENTY
      PARAMETER (ONE=1.0D0, THREE=3.0D0, TWENTY=20.0D0)
C begin
      GRID=ONE/THREE
      XGRIDF=1
      YGRIDF=1
      ZGRIDF=1
      MEMORY=500000
      BSCAL=TWENTY
      ELIM=7.0D0
      QMEMAUTO=.FALSE.
C
C check restrictions that are imposed by 3-d FFT algorithm
      CALL FFTPRP(DUMMY,FTPRIM,FTAVOI)
      RETURN
      END
C======================================================================
      SUBROUTINE XFFTUP(QDXYZ,QDB,START,STOP,
     &            XRH,XRK,XRL,FCALC,
     &            XRATOM,
     &            XRINDX,XRFQS,XRDX,XRDY,
     &            XRDZ,XRDT,XRDQ,XNCS,YNCS,ZNCS,DXNCS,DYNCS,DZNCS,
     &            NTBL,EXPTBL,QLOOK,XRSM,XRSA,XRSB,XRSC,XRFP,
     &            TSEL,XRNREF,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &            QHERM,XRCELL,MAPR,XRTR,XRINTR,XRVOL)
C
C Space-group general structure factor and atomic derivative
C routine.  Uses FFT method.
C
C Mode:
C    QDXYZ and QDB false: compute fcalcs.
C    QDXYZ or QDB true: compute derivatives with respect to
C                       x,y,z (QDXYZ) and/or to b,q (QDB).  The
C                       routine assumes that the derivatives of the
C                       target function are stored in FCALC.
C
C All atoms between START and STOP are used in the summation.
C This routine cannot deal with an imaginary form factor.
C The direct summation routine has to be used to compute
C the structure factor for the imaginary form factors.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xfft.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'ncs.inc'
      LOGICAL QDXYZ, QDB
      INTEGER START, STOP, XRH(*), XRK(*), XRL(*)
      INTEGER XRATOM(*), XRINDX(*)
      DOUBLE PRECISION XRFQS(*), XRDX(natom), XRDY(natom), XRDZ(natom)
      DOUBLE PRECISION XRDT(natom), XRDQ(natom) 
      DOUBLE PRECISION XNCS(natom), YNCS(natom), ZNCS(natom)
      DOUBLE PRECISION DXNCS(natom), DYNCS(natom), DZNCS(natom)
      INTEGER NTBL
      DOUBLE PRECISION EXPTBL(0:NTBL)
      LOGICAL QLOOK
      INTEGER XRSM
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(XRSM)
      INTEGER TSEL(*), XRNREF, XRNSYM, XRMSYM, XRSYTH
      DOUBLE COMPLEX FCALC(xrnref)
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(*), MAPR, XRTR(3,3), XRINTR(3,3), XRVOL
C local
      INTEGER NA, NB, NC, NAP, NBP, NCP, NAPP, NBPP, NCPP, NAPX, NBPX
      INTEGER NAPPP, NBPPP, NCPPP
      INTEGER I, J, IAT, REFLCT, INSYM, H, K, L
      LOGICAL ERROR, DONE
      DOUBLE PRECISION NABC(3), ABC(3,3), ABCINV(3,3), T(3), EXPLIM
      DOUBLE PRECISION SECS, SECS2, SECELD, SEC3RC, SECAB, SSQ
C pointer
      INTEGER RHO, CWK
C parameters
      DOUBLE PRECISION ONE, ZERO, TWO, FOUR
      PARAMETER (ONE=1.0D0, ZERO=0.0D0, TWO=2.0D0, FOUR=4.0D0)
C begin
      ERROR=.FALSE.
C
C initialize timers
      IF (TIMER.GT.0) THEN
      SECELD=ZERO
      SEC3RC=ZERO
      SECAB=ZERO
      END IF
C
C setup exponential lookup tables
      IF (QLOOK) CALL XFLOOK(NTBL,EXPTBL,EXPLIM)
C
C determine grid (Na, Nb, Nc) and sublattice
C
!$    if (.true.) then ! if OpenMP, the variable called "MEMORY" is meaningless
!$                     ! so call XFGRID right away: 
!$    CALL XFGRID(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,NAPX,NBPX,
!$   &             XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
!$   &             XRCELL,MAPR)
!$    else             ! below is the code that a non-OpenMP compiler sees
      DONE=.FALSE.
      DO WHILE (.NOT.DONE)
C
      IF (QMEMAUTO.AND.MEMORY.LT.DEFMEMA) THEN
      MEMORY=DEFMEMA
      WRITE(6,'(A,I10)')
     & ' %XFFT-AUTOmem: increasing memory allocation to ',MEMORY
      END IF
C
      CALL XFGRID(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,NAPX,NBPX,
     &             XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &             XRCELL,MAPR)
C
      IF (QMEMAUTO) THEN
      IF ((NAPP*NBPP*NCPP.GT.6).AND..NOT.QDXYZ.AND..NOT.QDB) THEN
      MEMORY=INT(MEMORY*NAPP*NBPP*NCPP/6.0D0)
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I10)')
     & ' %XFFT-AUTOmem: increasing memory allocation to ',MEMORY
      END IF
      ELSE
      DONE=.TRUE.
      END IF
      ELSE
      DONE=.TRUE.
      END IF
      END DO
!$    endif 
C
C check grid to make sure that it is fine enough
      CALL XFGRCHK(NA,NB,NC,XRNREF,.TRUE.,TSEL,XRH,XRK,XRL,ZERO)
C
      IF (WRNLEV.GE.5.AND..NOT.QDXYZ.AND..NOT.QDB) THEN
      WRITE(6,'(8(A,I4),A)')
     & ' XFFT: using grid [',NA,',',NB,',',NC,'] and sublattice [',
     & NAP,'(',NAPX,'),',NBP,'(',NBPX,'),',NCP,']'
      END IF
C
!$    if (.false.) then ! if OpenMP, then no warning message
      IF (NAPP*NBPP*NCPP.GT.12.AND.WRNLEV.GE.5.AND.
     &   .NOT.QDXYZ.AND..NOT.QDB) THEN
      WRITE(6,'(A)')
     &  ' *****************************************************'
      CALL WRNDIE(+1,'XFFT','Warning: Extreme factorization of FFT.')
      WRITE(6,'(A,I6,/,A,/,A)')
     &  ' %XFFT: Factorization of FFT=',NAPP*NBPP*NCPP,
     &  ' %XFFT: This may cause excessive CPU time usage.',
     &  ' %XFFT: Resubmit with larger XRAY FFT MEMOry specification.'
      WRITE(6,'(A,I12,A)')
     &  ' %XFFT: Recommendation: MEMOry=',
     &    INT(MEMORY*NAPP*NBPP*NCPP/12.D0),' or larger.'
      WRITE(6,'(A)')
     &  ' *****************************************************'
      END IF
!$    end if
C
C allocate HEAP space for sublattice grid
!$    if (.false.) then  ! when using OpenMP, don't call allhp here
      RHO=ALLHP(ICPLX8(NAPX*NBPX*(NCP/2+1)))
      CWK=ALLHP(ICPLX8(NAPX+2))
!$    endif
C
C compute sublattice transformation matrices
      NABC(1)=NAP
      NABC(2)=NBP
      NABC(3)=NCP
      DO I=1,3
      DO J=1,3
      ABC(I,J)=NABC(I)*XRTR(I,J)
      END DO
      END DO
      DO I=1,3
      DO J=1,3
      ABCINV(I,J)=XRINTR(I,J)/NABC(J)
      END DO
      END DO
C
C compute volume of FFT grid box
      FTVOL=XRVOL/(NA*NB*NC)
C
C compute norm of A*, B*, C* in FFT grid box
      FTAS=XRCELL(7)*NAP
      FTBS=XRCELL(8)*NBP
      FTCS=XRCELL(9)*NCP
C
C branch into appropriate part of the routine
      IF (.NOT.QDXYZ.AND..NOT.QDB) THEN
C
C ************************
C structure factor section
C ========================
C
C
C loop over all sublattices
      DO NAPPP=0,NAPP-1
      DO NBPPP=0,NBPP-1
!$omp parallel default(shared) num_threads(ncpp) 
!$omp& private(rho,cwk,ncppp,t,insym,xncs,yncs,zncs) firstprivate(timer)
!$omp& reduction(.OR.:error)
!$omp& reduction(+:fcalc) ! needs ifort 11.0; ifort 9 & 10 crash (issue 410656)
!$    RHO=ALLHP(ICPLX8(NAPX*NBPX*(NCP/2+1)))
!$    CWK=ALLHP(ICPLX8(NAPX+2))
!$omp do
      DO NCPPP=0,NCPP-1
!$    if (ncppp.lt.ncpp-1) timer = 0   ! only count times with last thread
C
C compute sublattice origin translation
      T(1)= NAPPP*XRINTR(1,1)/NA
     &     +NBPPP*XRINTR(1,2)/NB
     &     +NCPPP*XRINTR(1,3)/NC
      T(2)= NAPPP*XRINTR(2,1)/NA
     &     +NBPPP*XRINTR(2,2)/NB
     &     +NCPPP*XRINTR(2,3)/NC
      T(3)= NAPPP*XRINTR(3,1)/NA
     &     +NBPPP*XRINTR(3,2)/NB
     &     +NCPPP*XRINTR(3,3)/NC
C
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
C compute electron density in sublattice
C
C Initialize density to zero
      CALL RHOINI(NAP,NBP,NCP,NAPX,NBPX,HEAP(RHO))
C
C Loop over non-crystallographic symmetry operators to generate entire
C asymmetric unit of density from protomer coordinates.  The
C transformed coordinates are held in the arrays XNCS, YNCS, ZNCS.
      DO INSYM=1,XNNSYM
      CALL NCSXYZ(NATOM,X,Y,Z,INSYM,XNCS,YNCS,ZNCS)
      CALL FFTELS(ABC,ABCINV,T,NAPX,NBPX,NAP,NBP,NCP,START,STOP,
     &    XNCS,YNCS,
     &    ZNCS,XRFQS,HEAP(RHO),.TRUE.,.FALSE.,
     &    .FALSE.,XRINDX,XRDX,XRDY,XRDZ,XRDT,XRDQ,NTBL,EXPTBL,EXPLIM,
     &    QLOOK,XRSM,XRSA,XRSB,XRSC,XRFP,XRATOM)
      END DO
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SECELD=SECELD+SECS2-SECS
      END IF
C
      IF (.NOT.ERROR) THEN
C
C FFT electron density on sublattice
      IF (TIMER.GT.0) CALL VCPU(SECS)
      CALL FFT3RC2(NAPX,NBPX,NAP,NBP,NCP,HEAP(RHO),HEAP(CWK),ERROR)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SEC3RC=SEC3RC+SECS2-SECS
      END IF
C
      IF (.NOT.ERROR) THEN
C
C apply symmetry operators, and phase shifts and copy into A, B's
      IF (TIMER.GT.0) CALL VCPU(SECS)
      CALL FFTAB(NAPX,NBPX,NAP,NBP,NCP,HEAP(RHO),HEAP(CWK),
     &     NA,NB,NC,NAPPP,NBPPP,NCPPP,.FALSE.,
     &     XRH,XRK,XRL,FCALC
     &     ,TSEL,XRNREF,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SECAB=SECAB+SECS2-SECS
      END IF
C
      END IF
      END IF
      END DO
!$omp end do nowait
!$    CALL FREHP(CWK,ICPLX8(NAPX+2))
!$    CALL FREHP(RHO,ICPLX8(NAPX*NBPX*(NCP/2+1)))
!$omp end parallel
      END DO
      END DO
C
C
C apply B-factor correction to Fcalc
C ==================================
!$omp parallel do default(none) 
!$omp& shared(xrnref,tsel,xrh,xrk,xrl,xrtr,fcalc,bscal)
!$omp& private(reflct,h,k,l,ssq)
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).NE.0) THEN
      H=XRH(REFLCT)
      K=XRK(REFLCT)
      L=XRL(REFLCT)
      SSQ=  (XRTR(1,1)*H + XRTR(2,1)*K + XRTR(3,1)*L)**2
     &    + (XRTR(1,2)*H + XRTR(2,2)*K + XRTR(3,2)*L)**2
     &    + (XRTR(1,3)*H + XRTR(2,3)*K + XRTR(3,3)*L)**2
      FCALC(REFLCT)=FCALC(REFLCT)*DEXP(BSCAL*SSQ/FOUR)
      END IF
      END DO
C
C
      ELSE
C
C ******************
C derivative section
C ==================
C
C
C apply B-factor correction to derivatives of Fcalc
C and multiply by volume of FFT grid box (for FFT method)
C =======================================================
!$omp parallel do default(none) 
!$omp& shared(xrnref,tsel,xrh,xrk,xrl,xrtr,fcalc,bscal,ftvol)
!$omp& private(reflct,h,k,l,ssq)
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).NE.0) THEN
      H=XRH(REFLCT)
      K=XRK(REFLCT)
      L=XRL(REFLCT)
      SSQ=  (XRTR(1,1)*H + XRTR(2,1)*K + XRTR(3,1)*L)**2
     &    + (XRTR(1,2)*H + XRTR(2,2)*K + XRTR(3,2)*L)**2
     &    + (XRTR(1,3)*H + XRTR(2,3)*K + XRTR(3,3)*L)**2
C NOTE: there is a factor two which takes account of
C the symmetry f(h,k,l)=f*(-h,-k,-l)
      FCALC(REFLCT)=FCALC(REFLCT)*DEXP(BSCAL*SSQ/FOUR)*FTVOL/TWO
      END IF
      END DO
C
C loop over all sublattices
      DO NAPPP=0,NAPP-1
      DO NBPPP=0,NBPP-1
!$omp parallel default(shared) num_threads(ncpp)
!$omp& private(rho,cwk,ncppp,t,insym,dxncs,dyncs,dzncs,iat,xncs,yncs,
!$omp&         zncs) firstprivate(timer)
!$omp& reduction(+:xrdx,xrdy,xrdz,xrdt,xrdq) reduction(.OR.:error)
!$    RHO=ALLHP(ICPLX8(NAPX*NBPX*(NCP/2+1)))
!$    CWK=ALLHP(ICPLX8(NAPX+2))
!$omp do
      DO NCPPP=0,NCPP-1
!$    if (ncppp.lt.ncpp-1) timer = 0   ! only count times with last thread
C
C compute sublattice origin translation
      T(1)= NAPPP*XRINTR(1,1)/NA
     &     +NBPPP*XRINTR(1,2)/NB
     &     +NCPPP*XRINTR(1,3)/NC
      T(2)= NAPPP*XRINTR(2,1)/NA
     &     +NBPPP*XRINTR(2,2)/NB
     &     +NCPPP*XRINTR(2,3)/NC
      T(3)= NAPPP*XRINTR(3,1)/NA
     &     +NBPPP*XRINTR(3,2)/NB
     &     +NCPPP*XRINTR(3,3)/NC
C
      IF (TIMER.GT.0) CALL VCPU(SECS)
      CALL FFTAB(NAPX,NBPX,NAP,NBP,NCP,
     &     HEAP(RHO),HEAP(CWK),
     &     NA,NB,NC,NAPPP,NBPPP,NCPPP,.TRUE.,
     &     XRH,XRK,XRL,FCALC,
     &     TSEL,XRNREF,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SECAB=SECAB+SECS2-SECS
      END IF
C
C FFT electron density on sublattice
      IF (TIMER.GT.0) CALL VCPU(SECS)
      CALL FFT3CR2(NAPX,NBPX,NAP,NBP,NCP,HEAP(RHO),
     &            HEAP(CWK),ERROR)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SEC3RC=SEC3RC+SECS2-SECS
      END IF
C
      IF (.NOT.ERROR) THEN
C
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
C Loop over non-crystallographic symmetry operators. The derivatives
C are returned in the in the DXNCS arrays here so that the derivatives
C can be properly accumulated over the non-crystallographic symmetries.
      DO INSYM=1,XNNSYM
C Set protomer derivative accumulators to zero for passing to FFTELD
      DO IAT=START,STOP
      DXNCS(XRATOM(IAT))=ZERO
      DYNCS(XRATOM(IAT))=ZERO
      DZNCS(XRATOM(IAT))=ZERO
      END DO
C
      CALL NCSXYZ(NATOM,X,Y,Z,INSYM,XNCS,YNCS,ZNCS)
C compute derivatives
      CALL FFTELS(ABC,ABCINV,T,NAPX,NBPX,NAP,NBP,NCP,START,STOP,
     &            XNCS,YNCS,
     &            ZNCS,XRFQS,HEAP(RHO),
     &            .FALSE.,QDXYZ,QDB,
     &            XRINDX,DXNCS,DYNCS,DZNCS,XRDT,XRDQ,NTBL,
     &            EXPTBL,EXPLIM,QLOOK,XRSM,XRSA,XRSB,XRSC,XRFP,
     &            XRATOM)
C
      DO IAT=START,STOP
      XRDX(XRATOM(IAT))=XRDX(XRATOM(IAT))
     &  +NCSOP(INSYM,1,1)*DXNCS(XRATOM(IAT))
     &  +NCSOP(INSYM,2,1)*DYNCS(XRATOM(IAT))
     &  +NCSOP(INSYM,3,1)*DZNCS(XRATOM(IAT))
      XRDY(XRATOM(IAT))=XRDY(XRATOM(IAT))
     &  +NCSOP(INSYM,1,2)*DXNCS(XRATOM(IAT))
     &  +NCSOP(INSYM,2,2)*DYNCS(XRATOM(IAT))
     &  +NCSOP(INSYM,3,2)*DZNCS(XRATOM(IAT))
      XRDZ(XRATOM(IAT))=XRDZ(XRATOM(IAT))
     &  +NCSOP(INSYM,1,3)*DXNCS(XRATOM(IAT))
     &  +NCSOP(INSYM,2,3)*DYNCS(XRATOM(IAT))
     &  +NCSOP(INSYM,3,3)*DZNCS(XRATOM(IAT))
      END DO
      END DO
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SECELD=SECELD+SECS2-SECS
      END IF
C
      END IF
      END DO
!$omp end do nowait
!$    CALL FREHP(CWK,ICPLX8(NAPX+2))
!$    CALL FREHP(RHO,ICPLX8(NAPX*NBPX*(NCP/2+1)))
!$omp end parallel
      END DO
      END DO
      END IF
C
      IF (TIMER.GT.0) THEN
      WRITE(6,'(A,F10.4,/,A,F10.4,/,A,F10.4)')
     & ' XFFT: CPU-time:  electron-density-calc=',SECELD,
     & ' XFFT: CPU-time:  FFT of electron density=',SEC3RC,
     & ' XFFT: CPU-time:  transformation into A, B=',SECAB
      END IF
C
C free up temporary space for electron density map
!$    if (2.ne.2) then   ! when using OpenMP, don't call frehp here 
      CALL FREHP(CWK,ICPLX8(NAPX+2))
      CALL FREHP(RHO,ICPLX8(NAPX*NBPX*(NCP/2+1)))
!$    endif
C
      RETURN
      END
C======================================================================
      SUBROUTINE XFGRID(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,NAPX,NBPX,
     &             XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &             XRCELL,MAPR)
C
C determine integers NAPP, NAP, NBPP, NBP, NCPP, NCP such that
C 1) NA > a/(grid*high-resolution)
C    NB > b/(grid*high-resolution)
C    NC > c/(grid*high-resolution)
C 2) NA, NB, NC as small as possible
C 3) NAP*NAPP=NA, NBP*NBPP=NB, NCP*NCPP=NC,
C 4) NAP*NBP*NCP<MEMORY,
C 5) NAP, NBP, NCP as large as possible,
C 6) NCP is an even integer,
C 7) NAP, NBP, NCP have only prime factors allowed by FFT
C    (see parameters set in machine dependent routine FFTPRP).
C In addition we check the actual physical dimension of the 3-d
C matrix (NAPX, NBPX); to avoid data bank conflicts they may be
C increased somewhat.
C
C Currently NAPP, NBPP, NCPP are allowed only to be powers
C of two.  -- the distributed CNS source code does not enforce this! KD
C
C Author: Axel T. Brunger
C =======================
C
! OpenMP version: set NAPP=NBPP=1 and NCPP=number of threads
      IMPLICIT NONE
C I/O
      INCLUDE 'xfft.inc'
      INCLUDE 'consta.inc'
      INTEGER NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,NAPX,NBPX
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9), MAPR
C local
!$    integer omp_get_max_threads
      LOGICAL DONE
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      IF (MAPR.EQ.ZERO) THEN
      CALL WRNDIE(-5,'XFGRID','MAPResolution is not set.')
      ELSE
C
C make initial guess
      NA=INT(XRCELL(1)/(GRID*MAPR)+RSMALL)
      NB=INT(XRCELL(2)/(GRID*MAPR)+RSMALL)
      NC=INT(XRCELL(3)/(GRID*MAPR)+RSMALL)
      NAP=NA
      NBP=NB
      NCP=NC
      CALL XPRIME(NAP,1,1,FTPRIM)
      CALL XPRIME(NBP,1,1,FTPRIM)
      CALL XPRIME(NCP,1,2,FTPRIM)
      NAPP=1
      NBPP=1
!$      ncpp = max(1,min(omp_get_max_threads(),nc/2))
!$      NCP=NC/NCPP
!$      IF (NCP*NCPP.LT.NC) NCP=NCP+1
!$      CALL XPRIME(NCP,1,2,FTPRIM)
!$    if (2.ne.2) then ! an OpenMP compiler will disregard the following code
      NCPP=1
C
C loop until everything fits into memory
      DONE=.FALSE.
      DO WHILE (NAP*NBP*NCP.GT.MEMORY.AND..NOT.DONE)
      IF (NCP.GT.4) THEN
      NCPP=NCPP+1
      NCP=NC/NCPP
      IF (NCP*NCPP.LT.NC) NCP=NCP+1
C
C increase base by a factor of two because of real symmetry operation
C in z - direction
      CALL XPRIME(NCP,1,2,FTPRIM)
      ELSE IF (NBP.GT.4) THEN
      NBPP=NBPP+1
      NBP=NB/NBPP
      IF (NBP*NBPP.LT.NB) NBP=NBP+1
      CALL XPRIME(NBP,1,1,FTPRIM)
      ELSE IF (NAP.GT.4) THEN
      NAPP=NAPP+1
      NAP=NA/NAPP
      IF (NAP*NAPP.LT.NA) NAP=NAP+1
      CALL XPRIME(NAP,1,1,FTPRIM)
      ELSE
      DONE=.TRUE.
      END IF
      END DO
!$    endif
C
C compute actual NA, NB, NC
      NA=NAP*NAPP
      NB=NBP*NBPP
      NC=NCP*NCPP
C
C check physical dimension of 3-d array
      NAPX=NAP
      NBPX=NBP
      IF (FTAVOI.GT.0) THEN
      IF (MOD(NAPX,FTAVOI).EQ.0) NAPX=NAPX+1
      IF (MOD(NBPX,FTAVOI).EQ.0) NBPX=NBPX+1
      END IF
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XPRIME(N,M,BASE,PRIME)
C
C Routine increments integer N until BASE is a factor of N*M
C and N does not contain any primes greater than PRIME.
C
      IMPLICIT NONE
C I/O
      INTEGER N, M, BASE, PRIME
C local
      INTEGER NN, P
C begin
      NN=0
      N=N-1
      DO WHILE (NN.NE.1)
        N = N + 1
C
C increment N such that base is a factor of N*M
        IF (BASE .GT. 1) THEN
          NN = MOD(N * M, BASE)
          IF (NN .NE. 0) N = N + BASE - NN
        END IF
C
C divide N by integers less equal PRIME until N is itself a prime
        NN=N
        P=MIN(NN,PRIME)
        DO WHILE (P.GT.1)
          DO WHILE (MOD(NN,P).EQ.0)
            NN=NN/P
          END DO
          P=P-1
        END DO
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XFGRCHK(NA,NB,NC,XRNREF,QSEL,TSEL,XRH,XRK,XRL,FCALC)
C
C Routine checks to make sure that highest h,k,l indices
C are smaller than (na-1)/2, (nb-1)/2, (nc-1)/2, respectively
C IF QSEL is set to TRUE, the test will be done on all selected
C reflections, if QSEL is set to FALSE, the test will be done
C on all reflections with FCALC # 0.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER NA, NB, NC, XRNREF
      LOGICAL QSEL
      INTEGER TSEL(*), XRH(*), XRK(*), XRL(*)
      DOUBLE COMPLEX FCALC(*)
C local
      INTEGER IREF, HMAX, LMAX, KMAX
C begin
      HMAX=0
      KMAX=0
      LMAX=0
      IF (QSEL) THEN
!$omp parallel do default(none)
!$omp& shared(xrnref,tsel,fcalc,xrl,xrk,xrh) private(iref)
!$omp& reduction(max:hmax,kmax,lmax)
      DO IREF=1,XRNREF
      IF (TSEL(IREF).NE.0) THEN
      HMAX=MAX(HMAX,ABS(XRH(IREF)))
      KMAX=MAX(KMAX,ABS(XRK(IREF)))
      LMAX=MAX(LMAX,ABS(XRL(IREF)))
      END IF
      END DO
      ELSE
!$omp parallel do default(none)
!$omp& shared(xrnref,fcalc,xrl,xrk,xrh) private(iref)
!$omp& reduction(max:hmax,kmax,lmax)
      DO IREF=1,XRNREF
      IF (ABS(FCALC(IREF)).GT.RSMALL) THEN
      HMAX=MAX(HMAX,ABS(XRH(IREF)))
      KMAX=MAX(KMAX,ABS(XRK(IREF)))
      LMAX=MAX(LMAX,ABS(XRL(IREF)))
      END IF
      END DO
      END IF
C
      IF (HMAX.GT.(NA-1)/2) THEN
      CALL WRNDIE(-5,'XFGRCHK',
     & 'Grid in x-direction too coarse (FFT GRID or MAPRes).')
      ELSEIF (KMAX.GT.(NB-1)/2) THEN
      CALL WRNDIE(-5,'XFGRCHK',
     & 'Grid in y-direction too coarse (FFT GRID or MAPRes).')
      ELSEIF (LMAX.GT.(NC-1)/2) THEN
      CALL WRNDIE(-5,'XFGRCHK',
     & 'Grid in z-direction too coarse (FFT GRID or MAPRes).')
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE FFTAB(NAPX,NBPX,NAP,NBP,NCP,RHO,CWK,
     &                 NA,NB,NC,NAPPP,
     &                 NBPPP,NCPPP,QDERIV,XRH,XRK,XRL,FCALC,
     &                 TSEL,XRNREF,
     &                 XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY)
C
C Applies symmetry operators, and phase shifts and transforms
C RHO into A, B's or, for derivatives (QDERIV) vice versa.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'xfft.inc'
      INTEGER NAPX, NBPX, NAP, NBP, NCP
      DOUBLE COMPLEX RHO(NAPX,NBPX,*), CWK(*)
      INTEGER NA, NB, NC, NAPPP, NBPPP, NCPPP
      LOGICAL QDERIV
      INTEGER XRH(*), XRK(*), XRL(*)
      DOUBLE COMPLEX FCALC(*)
      INTEGER TSEL(*), XRNREF, XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
C local
      INTEGER H, K, L, ISYM, REFLCT, HH, KK, LL
      INTEGER LM, KM, ITEMP
      DOUBLE PRECISION RTH, RNA, RNB, RNC, TEMP, T1, T2, T3
      INTEGER S11, S12, S13, S21, S22, S23, S31, S32, S33
      LOGICAL COND
      DOUBLE COMPLEX FLOC, RLOC
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, FOUR=4.0D0)
      DOUBLE COMPLEX CIMAG
      CIMAG=DCMPLX(ZERO,ONE)
C begin
C
C initialize RHO
      IF (QDERIV) THEN
      DO L=1,NCP/2+1
      DO K=1,NBP
      DO H=1,NAP
      RHO(H,K,L)=DCMPLX(ZERO,ZERO)
      END DO
      END DO
      END DO
      END IF
C
      RTH=TWO*PI/XRSYTH
      RNA=TWO*PI/NA
      RNB=TWO*PI/NB
      RNC=TWO*PI/NC
      RNA=RNA*NAPPP
      RNB=RNB*NBPPP
      RNC=RNC*NCPPP
C
C loop over all symmetry operators
      DO ISYM=1,XRNSYM
C
      T1=XRSYMM(ISYM,1,4)*RTH
      T2=XRSYMM(ISYM,2,4)*RTH
      T3=XRSYMM(ISYM,3,4)*RTH
C
      S11=XRSYMM(ISYM,1,1)
      S21=XRSYMM(ISYM,2,1)
      S31=XRSYMM(ISYM,3,1)
      S12=XRSYMM(ISYM,1,2)
      S22=XRSYMM(ISYM,2,2)
      S32=XRSYMM(ISYM,3,2)
      S13=XRSYMM(ISYM,1,3)
      S23=XRSYMM(ISYM,2,3)
      S33=XRSYMM(ISYM,3,3)
C
C
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).NE.0) THEN
C
C apply transpose of symmetry operator to H, K, L
      HH=XRH(REFLCT)
      KK=XRK(REFLCT)
      LL=XRL(REFLCT)
      H=S11*HH + S21*KK + S31*LL
      K=S12*HH + S22*KK + S32*LL
      L=S13*HH + S23*KK + S33*LL
C
C compute phase shift of translational part of symmetry operator,
C multiply with phase shift of sublattice origin translation,
C and apply B-factor scale to minimize aliasing
      TEMP= T1*HH + T2*KK + T3*LL + H*RNA + K*RNB + L*RNC
C
C complex-multiply everything and store factor in temporary FLOC array
      FLOC=DCMPLX(DCOS(TEMP),DSIN(TEMP))
C
C project H,K,L into periodic sublattices (1...NAP,...)
      H=MOD(H+ 10000*NAP,NAP) +1
      K=MOD(K+ 10000*NBP,NBP) +1
      L=MOD(L+ 10000*NCP,NCP) +1
C
C map H,K,L into L<=NCP/2+1 hemisphere (real RHO symmetry!!)
      COND=(L.GT.NCP/2+1)
      IF (COND) H=MOD(NAP-H+1,NAP)+1
      IF (COND) K=MOD(NBP-K+1,NBP)+1
      IF (COND) L=MOD(NCP-L+1,NCP)+1
C
C
C branch according to whether derivative or structure factor is
C computed
      IF (QDERIV) THEN
C
C derivative section: multiply FCALC=(A, B) derivative with factor F
C (note that the derivative of a constant factor is conjugate, as
C are the factors for the subsequent Fourier transform)
      FLOC=FLOC*DCONJG(FCALC(REFLCT))
      IF (COND) FLOC=DCONJG(FLOC)
C
      RHO(H,K,L)=RHO(H,K,L)+FLOC
C
      ELSE
C
C structure factor section:
C gather RHO
      RLOC=RHO(H,K,L)
C
C multiply RL with factor F and accumulate in FCALC
      IF (COND) RLOC=DCONJG(RLOC)
      FCALC(REFLCT)=FCALC(REFLCT)+FLOC*RLOC*FTVOL
C
      END IF
C
      END IF
      END DO
      END DO
C
      IF (QDERIV) THEN
C
C for derivatives, the data in RHO are scrambled for L=-L.
C Note that we add RHO and RHO of the Friedel mate.  If the
C data consists of a hemisphere, only
C one of the mates is going to be non-zero except
C when H.EQ.HM and K.EQ.KM and L.EQ.LM.  In this case the formula is
C still correct since the FFT computes the derivative as if all
C Friedel mates are included in the objective function.  If the
C data consists of a full sphere, things are correct since we're
C accumulating the contributions to a particular H,K,L.
C Note also the factor 1/2 in the calling routine.
      DO L=1,NCP/2+1
      LM=MOD(NCP+1-L,NCP)+1
      IF (L.EQ.LM) THEN
      DO K=1,NBP/2+1
      KM=MOD(NBP+1-K,NBP)+1
      IF (K.EQ.KM) THEN
      ITEMP=NAP/2+1
      ELSE
      ITEMP=NAP
      END IF
      CWK(1)=RHO(1,KM,LM)
      DO H=2,ITEMP
      CWK(H)=RHO(NAP-H+2,KM,LM)
      END DO
      DO H=1,ITEMP
C
C add complex conjugate of CWK to RHO
      RHO(H,K,L)=RHO(H,K,L)+DCONJG(CWK(H))
C
C store complex conjugate of RHO in CWK
      CWK(H)=DCONJG(RHO(H,K,L))
      END DO
      DO H=2,ITEMP
      RHO(NAP-H+2,KM,LM)=CWK(H)
      END DO
      RHO(1,KM,LM)=CWK(1)
      END DO
      END IF
      END DO
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE FFTELS(ABC,ABCINV,T,NAPX,NBPX,NAP,NBP,NCP,START,STOP,
     & XL,YL,ZL,XRFQS,RHO,QSF,QDXYZ,QDB,
     & XRINDX,XRDX,XRDY,XRDZ,XRDT,XRDQ,
     & NTBL,EXPTBL,EXPLIM,QLOOK,XRSM,XRSA,XRSB,XRSC,XRFP,XRATOM)
C
C    Scalar version of electron density routine.
C
C    Computes electron density RHO or derivatives on lattice given by
C T+ i*A + j*B + k*C, where i,j,k are integers and (A,B,C)=ABCINV.
C RHO is accumulated over start...stop atoms X,Y,Z in unit cell.
C RHO is periodic, i.e.
C
C    RHO(i,j,k) = RHO(i+u*NA, j+v*NB, k+w*NC) for arbitrary u, v, w.
C
C    The individual atomic formfactors are given by
C QMAIN * exp(-b*s**2/4) * (c + sum(i=1,4) a(i) *exp(b(i) s**2/4))
C where QMAIN and b are the occupancies and temperature factors
C respectively.
C
C    On output RHO is a complex array which is related to the real
C electron r density by
C
C     RHO(I,J,K/2)  =  r(I,J,2*K-1) + i*r(I,J,2K)
C
C    This is done to save memory and operations in the FFT
C routines.
C
C
C note: the electron density array RHO must be initialized before
C calling this routine.
C
C This routine either computes structure factors (QSF),
C or derivatives (QDXYZ or QDB).  If QDXYZ is specified,
C derivatives with respect to X,Y,Z are computed and returned
C in XRDX, XRDY, XRDZ.  If QDB is specified, derivatives with
C respect to Q, B, and f' are computed and returned in XRDQ, XRDB, XRDZ.
C At present, QDXYZ and QDB are mutually exclusive, although this
C could be changed (was done to save memory).
C
C
C The purpose of this routine is to allocate space for
C temporary arrays.
C
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      DOUBLE PRECISION ABC(3,3), ABCINV(3,3), T(3)
      INTEGER NAPX, NBPX, NAP, NBP, NCP, START, STOP
      DOUBLE PRECISION XL(*), YL(*), ZL(*), XRFQS(*)
      DOUBLE COMPLEX RHO(NAPX,NBPX,*)
      LOGICAL QSF, QDXYZ, QDB
      INTEGER XRINDX(*)
      DOUBLE PRECISION XRDX(*), XRDY(*), XRDZ(*), XRDT(*), XRDQ(*)
      INTEGER NTBL
      DOUBLE PRECISION EXPTBL(0:NTBL), EXPLIM
      LOGICAL QLOOK
      INTEGER XRSM
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(XRSM)
      INTEGER XRATOM
C pointer
      INTEGER AAMOD, BBMOD
C begin
C
C
      AAMOD=ALLHP(INTEG4(NAP*2+1))
      BBMOD=ALLHP(INTEG4(NBP*2+1))
      CALL FFTES2(ABC,ABCINV,T,NAPX,NBPX,NAP,NBP,NCP,START,STOP,
     & XL,YL,ZL,XRFQS,RHO,QSF,QDXYZ,QDB,
     & XRINDX,XRDX,XRDY,XRDZ,XRDT,XRDQ,
     & NTBL,EXPTBL,EXPLIM,QLOOK,XRSM,XRSA,XRSB,XRSC,XRFP,
     & HEAP(AAMOD),HEAP(BBMOD),XRATOM)
      CALL FREHP(BBMOD,INTEG4(NBP*2+1))
      CALL FREHP(AAMOD,INTEG4(NAP*2+1))
C
      RETURN
      END
C======================================================================
      SUBROUTINE FFTES2(ABC,ABCINV,T,NAPX,NBPX,NAP,NBP,NCP,START,STOP,
     & XL,YL,ZL,XRFQS,RHO,QSF,QDXYZ,QDB,
     & XRINDX,XRDX,XRDY,XRDZ,XRDT,XRDQ,
     & NTBL,EXPTBL,EXPLIM,QLOOK,XRSM,XRSA,XRSB,XRSC,XRFP,
     & AAMOD,BBMOD,XRATOM)
C
C See routine FFTEL2 above
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'xfft.inc'
      INCLUDE 'coord.inc'
      DOUBLE PRECISION ABC(3,3), ABCINV(3,3), T(3)
      INTEGER NAPX, NBPX, NAP, NBP, NCP, START, STOP
      DOUBLE PRECISION XL(*), YL(*), ZL(*), XRFQS(*)
      DOUBLE COMPLEX RHO(NAPX,NBPX,*)
      LOGICAL QSF, QDXYZ, QDB
      INTEGER XRINDX(*)
      DOUBLE PRECISION XRDX(*), XRDY(*), XRDZ(*), XRDT(*), XRDQ(*)
      INTEGER NTBL
      DOUBLE PRECISION EXPTBL(0:NTBL), EXPLIM
      LOGICAL QLOOK
      INTEGER XRSM
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(XRSM)
      INTEGER AAMOD(-NAP:NAP), BBMOD(-NBP:NBP), XRATOM(*)
C local
      INTEGER A, B, C, AX, BX, CX, IAT, AA, BB, CC, CR
      DOUBLE PRECISION RCUT, RCUT2, DX, DY, DZ, DELX, DELY, DELZ, DEL2
      DOUBLE PRECISION DELXB, DELYB, DELZB, DELXC, DELYC, DELZC
      DOUBLE PRECISION OFFSETX, OFFSETY, OFFSETZ
      DOUBLE PRECISION AEQ1, AEQ2, AEQ3, AEQ4, AEQ5, AEEC
      DOUBLE PRECISION AE1, AE2, AE3, AE4, AE5
      DOUBLE PRECISION BE1, BE2, BE3, BE4, BE5
      DOUBLE PRECISION BEE1, BEE2, BEE3, BEE4, BEE5
      DOUBLE PRECISION ABE1, ABE2, ABE3, ABE4, ABE5, TEMP
      DOUBLE PRECISION CE1, CE2, CE3, CE4, CE5, DE1, DE2, DE3
      DOUBLE PRECISION DE4, DE5, DF
      INTEGER AMAXP, AMAXN, BMAXP, BMAXN, CMAXP, CMAXN
C parameters
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, THREE, FOUR
      PARAMETER (ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0, TWO=2.0D0)
      PARAMETER (THREE=3.0D0, FOUR=4.0D0)
C begin
C
C loop over all atoms.
      DO IAT=START, STOP
      AEQ1=XRFQS(XRATOM(IAT))*XRSA(XRINDX(IAT),1)*(SQRT(FOUR*PI/(
     & XRSB(XRINDX(IAT),1)+WMAIN(XRATOM(IAT))+BSCAL)))**3
      AEQ2=XRFQS(XRATOM(IAT))*XRSA(XRINDX(IAT),2)*(SQRT(FOUR*PI/(
     & XRSB(XRINDX(IAT),2)+WMAIN(XRATOM(IAT))+BSCAL)))**3
      AEQ3=XRFQS(XRATOM(IAT))*XRSA(XRINDX(IAT),3)*(SQRT(FOUR*PI/(
     & XRSB(XRINDX(IAT),3)+WMAIN(XRATOM(IAT))+BSCAL)))**3
      AEQ4=XRFQS(XRATOM(IAT))*XRSA(XRINDX(IAT),4)*(SQRT(FOUR*PI/(
     & XRSB(XRINDX(IAT),4)+WMAIN(XRATOM(IAT))+BSCAL)))**3
      AEQ5=XRFQS(XRATOM(IAT))*(XRSC(XRINDX(IAT))+XRFP(XRINDX(IAT)))
     &      *(SQRT(FOUR*PI/(
     & WMAIN(XRATOM(IAT))+BSCAL)))**3
      AE1=AEQ1*QMAIN(XRATOM(IAT))
      AE2=AEQ2*QMAIN(XRATOM(IAT))
      AE3=AEQ3*QMAIN(XRATOM(IAT))
      AE4=AEQ4*QMAIN(XRATOM(IAT))
      AE5=AEQ5*QMAIN(XRATOM(IAT))
      AEEC=XRFQS(XRATOM(IAT))*QMAIN(XRATOM(IAT))*(SQRT(FOUR*PI/(
     & WMAIN(XRATOM(IAT))+BSCAL)))**3
      BE1=FOUR*PI**2/(XRSB(XRINDX(IAT),1)+WMAIN(XRATOM(IAT))+BSCAL)
      BE2=FOUR*PI**2/(XRSB(XRINDX(IAT),2)+WMAIN(XRATOM(IAT))+BSCAL)
      BE3=FOUR*PI**2/(XRSB(XRINDX(IAT),3)+WMAIN(XRATOM(IAT))+BSCAL)
      BE4=FOUR*PI**2/(XRSB(XRINDX(IAT),4)+WMAIN(XRATOM(IAT))+BSCAL)
      BE5=FOUR*PI**2/(WMAIN(XRATOM(IAT))+BSCAL)
C
C determine individual cutoff for electron density calculation
      RCUT2=MAX(ZERO,
     & ( ELIM+LOG(MAX( ABS(AE1) , RSMALL )) )/BE1,
     & ( ELIM+LOG(MAX( ABS(AE2) , RSMALL )) )/BE2,
     & ( ELIM+LOG(MAX( ABS(AE3) , RSMALL )) )/BE3,
     & ( ELIM+LOG(MAX( ABS(AE4) , RSMALL )) )/BE4,
     & ( ELIM+LOG(MAX( ABS(AE5) , RSMALL )) )/BE5)
C
C compute maximum atomic radius
      RCUT=SQRT(RCUT2)
C
C precompute special coefficients for lookup-table mode
      IF (QLOOK) THEN
      BEE1=BE1*EXPLIM
      BEE2=BE2*EXPLIM
      BEE3=BE3*EXPLIM
      BEE4=BE4*EXPLIM
      BEE5=BE5*EXPLIM
      END IF
C
C precompute coefficients for positional derivatives
      IF (QDXYZ) THEN
      ABE1=AE1*BE1
      ABE2=AE2*BE2
      ABE3=AE3*BE3
      ABE4=AE4*BE4
      ABE5=AE5*BE5
      END IF
C
C precompute coefficients for B-factor derivatives
      IF (QDB) THEN
      CE1=-(THREE/TWO)*AE1/(XRSB(XRINDX(IAT),1)
     &    +WMAIN(XRATOM(IAT))+BSCAL)
      CE2=-(THREE/TWO)*AE2/(XRSB(XRINDX(IAT),2)
     &    +WMAIN(XRATOM(IAT))+BSCAL)
      CE3=-(THREE/TWO)*AE3/(XRSB(XRINDX(IAT),3)
     &    +WMAIN(XRATOM(IAT))+BSCAL)
      CE4=-(THREE/TWO)*AE4/(XRSB(XRINDX(IAT),4)
     &    +WMAIN(XRATOM(IAT))+BSCAL)
      CE5=-(THREE/TWO)*AE5/(WMAIN(XRATOM(IAT))+BSCAL)
      DE1=AE1*FOUR*PI**2/(XRSB(XRINDX(IAT),1)
     &    +WMAIN(XRATOM(IAT))+BSCAL)**2
      DE2=AE2*FOUR*PI**2/(XRSB(XRINDX(IAT),2)
     &    +WMAIN(XRATOM(IAT))+BSCAL)**2
      DE3=AE3*FOUR*PI**2/(XRSB(XRINDX(IAT),3)
     &    +WMAIN(XRATOM(IAT))+BSCAL)**2
      DE4=AE4*FOUR*PI**2/(XRSB(XRINDX(IAT),4)
     &    +WMAIN(XRATOM(IAT))+BSCAL)**2
      DE5=AE5*FOUR*PI**2/(WMAIN(XRATOM(IAT))+BSCAL)**2
      END IF
C
C determine indices of lattice point closest to atomic position x,y,z
      DX= ABC(1,1)*(XL(XRATOM(IAT))-T(1))
     &   +ABC(1,2)*(YL(XRATOM(IAT))-T(2))
     &   +ABC(1,3)*(ZL(XRATOM(IAT))-T(3))
      AX=NINT(DX)
      DY= ABC(2,1)*(XL(XRATOM(IAT))-T(1))
     &   +ABC(2,2)*(YL(XRATOM(IAT))-T(2))
     &   +ABC(2,3)*(ZL(XRATOM(IAT))-T(3))
      BX=NINT(DY)
      DZ= ABC(3,1)*(XL(XRATOM(IAT))-T(1))
     &   +ABC(3,2)*(YL(XRATOM(IAT))-T(2))
     &   +ABC(3,3)*(ZL(XRATOM(IAT))-T(3))
      CX=NINT(DZ)
C
C determine the difference between lattice point AX, AY, AZ and
C actual position
      DX=AX-DX
      DY=BX-DY
      DZ=CX-DZ
C
C convert this difference into real space
      OFFSETX=ABCINV(1,3)*DZ+ABCINV(1,2)*DY+ABCINV(1,1)*DX
      OFFSETY=ABCINV(2,3)*DZ+ABCINV(2,2)*DY+ABCINV(2,1)*DX
      OFFSETZ=ABCINV(3,3)*DZ+ABCINV(3,2)*DY+ABCINV(3,1)*DX
C
C determine AMAXP, AMAXN, BMAXP, BMAXN, CMAXP, CMAXN of a
C box isomorphous to the unitcell geometry which surrounds
C the sphere with the atomic radius.  Compute maximum box.
      AMAXP=MIN(NAP,INT(RCUT*FTAS-DX+RSMALL))
      BMAXP=MIN(NBP,INT(RCUT*FTBS-DY+RSMALL))
      CMAXP=MIN(NCP,INT(RCUT*FTCS-DZ+RSMALL))
      AMAXN=MIN(NAP,INT(RCUT*FTAS+DX+RSMALL))
      BMAXN=MIN(NBP,INT(RCUT*FTBS+DY+RSMALL))
      CMAXN=MIN(NCP,INT(RCUT*FTCS+DZ+RSMALL))
C
C precompute some quantities
      DO A=-AMAXN,+AMAXP
      AAMOD(A)=MOD(A+AX +10000*NAP,NAP) +1
      END DO
      DO B=-BMAXN,+BMAXP
      BBMOD(B)=MOD(B+BX +10000*NBP,NBP) +1
      END DO
C
      DO C=-CMAXN,CMAXP
      DELXC=OFFSETX+ABCINV(1,3)*C
      DELYC=OFFSETY+ABCINV(2,3)*C
      DELZC=OFFSETZ+ABCINV(3,3)*C
C
      CC=MOD(C+CX +10000*NCP,NCP) +1
      CR=MOD(CC,2)
      CC=(CC+1)/2
C
      DO B=-BMAXN,BMAXP
      DELXB=DELXC+ABCINV(1,2)*B
      DELYB=DELYC+ABCINV(2,2)*B
      DELZB=DELZC+ABCINV(3,2)*B
      BB=BBMOD(B)
C
      DO A=-AMAXN,AMAXP
      DELX=ABCINV(1,1)*A+DELXB
      DELY=ABCINV(2,1)*A+DELYB
      DELZ=ABCINV(3,1)*A+DELZB
      DEL2=DELX**2+DELY**2+DELZ**2
      IF (RCUT2.GE.DEL2) THEN
      AA=AAMOD(A)
C
      IF (QSF) THEN
C
C structure factor section
C ========================
      IF (QLOOK) THEN
C
C use lookup tables for exponential function
      TEMP=AE1*EXPTBL(MIN(NTBL,INT(BEE1*DEL2)))
     &    +AE2*EXPTBL(MIN(NTBL,INT(BEE2*DEL2)))
     &    +AE3*EXPTBL(MIN(NTBL,INT(BEE3*DEL2)))
     &    +AE4*EXPTBL(MIN(NTBL,INT(BEE4*DEL2)))
     &    +AE5*EXPTBL(MIN(NTBL,INT(BEE5*DEL2)))
C
      ELSE
C
C use internal exponential function
      TEMP=AE1*EXP(-BE1*DEL2)
     &    +AE2*EXP(-BE2*DEL2)
     &    +AE3*EXP(-BE3*DEL2)
     &    +AE4*EXP(-BE4*DEL2)
     &    +AE5*EXP(-BE5*DEL2)
      END IF
C
      IF (CR.EQ.1) THEN
C
C odd points
      RHO(AA,BB,CC)=DCMPLX(
     & DBLE(RHO(AA,BB,CC))+TEMP,
     & DIMAG(RHO(AA,BB,CC)))
      ELSE
C
C even points
      RHO(AA,BB,CC)=DCMPLX(
     & DBLE(RHO(AA,BB,CC)),
     & DIMAG(RHO(AA,BB,CC))+TEMP)
      END IF
C
      END IF
C
      IF (QDB.OR.QDXYZ) THEN
C
C common section for positional and B-factor derivatives
C ======================================================
C
C gather RHO
      IF (CR.EQ.1) THEN
C
C odd points
      TEMP=DBLE(RHO(AA,BB,CC))
      ELSE
C
C even points
      TEMP=DIMAG(RHO(AA,BB,CC))
      END IF
C
      IF (QDXYZ) THEN
C
C positional derivative section
C =============================
C
C Accumulate derivatives over the atoms passed to this routine.  Global
C accumulation (over non-crystallographic symmetries) is
C done in calling routine.
      IF (QLOOK) THEN
      DF=TWO*TEMP*
     & (ABE1*EXPTBL(MIN(NTBL,INT(BEE1*DEL2)))
     & +ABE2*EXPTBL(MIN(NTBL,INT(BEE2*DEL2)))
     & +ABE3*EXPTBL(MIN(NTBL,INT(BEE3*DEL2)))
     & +ABE4*EXPTBL(MIN(NTBL,INT(BEE4*DEL2)))
     & +ABE5*EXPTBL(MIN(NTBL,INT(BEE5*DEL2))))
      XRDX(XRATOM(IAT))=XRDX(XRATOM(IAT))+DF*DELX
      XRDY(XRATOM(IAT))=XRDY(XRATOM(IAT))+DF*DELY
      XRDZ(XRATOM(IAT))=XRDZ(XRATOM(IAT))+DF*DELZ
      ELSE
      DF=TWO*TEMP*
     & (ABE1*EXP(-BE1*DEL2)
     & +ABE2*EXP(-BE2*DEL2)
     & +ABE3*EXP(-BE3*DEL2)
     & +ABE4*EXP(-BE4*DEL2)
     & +ABE5*EXP(-BE5*DEL2))
      XRDX(XRATOM(IAT))=XRDX(XRATOM(IAT))+DF*DELX
      XRDY(XRATOM(IAT))=XRDY(XRATOM(IAT))+DF*DELY
      XRDZ(XRATOM(IAT))=XRDZ(XRATOM(IAT))+DF*DELZ
      END IF
      END IF
C
      IF (QDB) THEN
C
C B-factor and occupancy derivative section
C =========================================
      IF (QLOOK) THEN
C
C b-factor
      XRDT(XRATOM(IAT))=XRDT(XRATOM(IAT))+TEMP*(
     & (CE1+DE1*DEL2)*
     & EXPTBL(MIN(NTBL,INT(BEE1*DEL2)))
     & +(CE2+DE2*DEL2)*
     & EXPTBL(MIN(NTBL,INT(BEE2*DEL2)))
     & +(CE3+DE3*DEL2)*
     & EXPTBL(MIN(NTBL,INT(BEE3*DEL2)))
     & +(CE4+DE4*DEL2)*
     & EXPTBL(MIN(NTBL,INT(BEE4*DEL2)))
     & +(CE5+DE5*DEL2)*
     & EXPTBL(MIN(NTBL,INT(BEE5*DEL2))))
C
C occupancy
      XRDQ(XRATOM(IAT))=XRDQ(XRATOM(IAT))+TEMP*(
     & +AEQ1*EXPTBL(MIN(NTBL,INT(BEE1*DEL2)))
     & +AEQ2*EXPTBL(MIN(NTBL,INT(BEE2*DEL2)))
     & +AEQ3*EXPTBL(MIN(NTBL,INT(BEE3*DEL2)))
     & +AEQ4*EXPTBL(MIN(NTBL,INT(BEE4*DEL2)))
     & +AEQ5*EXPTBL(MIN(NTBL,INT(BEE5*DEL2))))
C
C f'
      XRDX(XRATOM(IAT))=XRDX(XRATOM(IAT))+TEMP*
     & AEEC*EXPTBL(MIN(NTBL,INT(BEE5*DEL2)))
C
      ELSE
C
C b-factor
      XRDT(XRATOM(IAT))=XRDT(XRATOM(IAT))+TEMP*(
     & (CE1+DE1*DEL2)*EXP(-BE1*DEL2)
     & +(CE2+DE2*DEL2)*EXP(-BE2*DEL2)
     & +(CE3+DE3*DEL2)*EXP(-BE3*DEL2)
     & +(CE4+DE4*DEL2)*EXP(-BE4*DEL2)
     & +(CE5+DE5*DEL2)*EXP(-BE5*DEL2))
C
C occupancy
      XRDQ(XRATOM(IAT))=XRDQ(XRATOM(IAT))+TEMP*(
     & +AEQ1*EXP(-BE1*DEL2)
     & +AEQ2*EXP(-BE2*DEL2)
     & +AEQ3*EXP(-BE3*DEL2)
     & +AEQ4*EXP(-BE4*DEL2)
     & +AEQ5*EXP(-BE5*DEL2))
C
C f'
      XRDX(XRATOM(IAT))=XRDX(XRATOM(IAT))+TEMP*AEEC*EXP(-BE5*DEL2)
C
      END IF
      END IF
      END IF
C
      END IF
C
      END DO
      END DO
      END DO
C
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XFLOOK(NTBL,EXPTBL,EXPLIM)
C
C routine sets up lookup tables for EXP. Used in XFFTELD and XFFTELS
C
C Author: Axel T. Brunger
C =======================
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'xfft.inc'
      INTEGER NTBL
      DOUBLE PRECISION EXPTBL(0:NTBL), EXPLIM
C local
      INTEGER I
C parameters
      DOUBLE PRECISION THREE
      PARAMETER (THREE=3.0D0)
C begin
C
C go a bit beyond ELIM...
      EXPLIM=(NTBL-1)/(ELIM+THREE)
C sin and cos
      DO I=0,NTBL
      EXPTBL(I)=EXP( -I/EXPLIM )
      END DO
C
C to lookup do I=MIN(NTBL,INT( -R*EXPLIM)), R has to be positive.
      RETURN
      END
C======================================================================
      SUBROUTINE RHOINI(NAP,NBP,NCP,NAPX,NBPX,RHO)
C
C Initialize density on sublattice to ZERO.
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C i/o
      INTEGER NAP,NBP,NCP,NAPX,NBPX,A,B,C
      DOUBLE COMPLEX RHO(NAPX,NBPX,*)
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
C initialize RHO for structure factor calc.
      DO C=1,NCP/2+1
      DO B=1,NBP
      DO A=1,NAP
      RHO(A,B,C)=DCMPLX(ZERO,ZERO)
      END DO
      END DO
      END DO
      RETURN
      END
C
C=======================================================================
C{
      subroutine mFFTcr(d, n, FTAvoi, m)
      IMPLICIT NONE
C I/O
      integer  d, n(d), FTAvoi, m(d)
C
C Determine physical dimensions m(i) of d-dimensional array suitable
C for complex-to-real FFT.
C n(i) is the number of grid points along dimension i.
C
C local
      integer  i
C
C begin
C m(1) must be even (because 2 real = 1 complex)
C m(1) must be >= max(1,n(1)+2)
C m(1)/2 must not contain the factor FTAvoi
      m(1) = (n(1) + 3) / 2
      do i = 2, d
        m(i) = n(i)
      end do
C
      if (FTAvoi .gt. 0) then
        do i = 1, d - 1
          if (mod(m(i), FTAvoi) .eq. 0) m(i) = m(i) + 1
        end do
      end if
C
      m(1) = m(1) * 2
C
      return
      end
C}
C=======================================================================
