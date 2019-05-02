C=======================================================================
C Unified Conventional and F2F2 FFT Translation Function
C
C The F2F2 FFT Translation Function code is based on:
C   J. Navaza & E. Vernoslova
C   "On the Fast Translation Functions for Molecular Replacement"
C   Acta Cryst. (1995), A51, 445-449
C
C Implementation by R.W. Grosse-Kunstleve, 1997
C=======================================================================
C{
      subroutine HmulMt(h, k, l,
     &                  mSymMx, SymMx, iSym,
     &                  hM, kM, lM, Ht, STBF)
      IMPLICIT NONE
C I/O
      integer  h, k, l
      integer  mSymMx, SymMx(mSymMx, 3, 4), iSym
      integer  hM, kM, lM, Ht, STBF
C
C Multiply hkl by symmetry matrix (rotation & translation parts)
C
C begin
      hM =   SymMx(iSym, 1, 1) * h
     &     + SymMx(iSym, 2, 1) * k
     &     + SymMx(iSym, 3, 1) * l
      kM =   SymMx(iSym, 1, 2) * h
     &     + SymMx(iSym, 2, 2) * k
     &     + SymMx(iSym, 3, 2) * l
      lM =   SymMx(iSym, 1, 3) * h
     &     + SymMx(iSym, 2, 3) * k
     &     + SymMx(iSym, 3, 3) * l
      Ht =   SymMx(iSym, 1, 4) * h
     &     + SymMx(iSym, 2, 4) * k
     &     + SymMx(iSym, 3, 4) * l
C
      Ht = mod(Ht, STBF)
      if (Ht .lt. 0) Ht = Ht + STBF
C
      return
      end
C}
C=======================================================================
C{
      subroutine DimH3D(nRef, xrh, xrk, xrl, TSel, QHerm,
     &                  h, n)
      IMPLICIT NONE
C I/O
      integer  nRef
      integer  xrh(*), xrk(*), xrl(*), TSel(*)
      logical  QHerm
      integer  h(3, 2), n(3)
C
C Compute dimensions for 3D hkl array
C
C local
      integer  iRef, j
C
C begin
      do j = 1, 3
        h(j, 2) = 0
      end do
C
      do iRef = 1, nRef
        if (TSel(iRef) .ne. 0) then
          h(1, 2) = max(h(1, 2), xrh(iRef), -xrh(iRef))
          h(2, 2) = max(h(2, 2), xrk(iRef), -xrk(iRef))
          h(3, 2) = max(h(3, 2), xrl(iRef), -xrl(iRef))
        end if
      end do
C
      do j = 1, 3
        h(j, 1) = -h(j, 2)
      end do
C
      if (QHerm) h(3, 1) = 0
C
      do j = 1, 3
        n(j) = h(j, 2) - h(j, 1) + 1
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine hklAsyPt(nRef, xrh, xrk, xrl, TSel, AsyPt, QHerm,
     &                    mSymMx, SymMx, nSymMx, STBF,
     &                    WrnLev)
      IMPLICIT NONE
C I/O
      integer  nRef, xrh(*), xrk(*), xrl(*), TSel(*), AsyPt(2, *)
      logical  QHerm
      integer  mSymMx, SymMx(mSymMx, 3, 4), nSymMx, STBF
      integer  WrnLev
C
C Compute asymmetric unit in reciprocal space for use in fcAtX.
C
C local
      integer  iRef, jRef, iSym
      integer  h, k, l, hM, kM, lM, Ht, sM, hA, kA, lA, HtA, sA
      integer  nSel, nAsy, nSym
C
C externals
      integer   CmpHKL, RecHem
      external  CmpHKL, RecHem
C
C begin
      nSel = 0
      nAsy = 0
      nSym = 0
C
      do iRef = 1, nRef
        if (TSel(iRef) .ne. 0) then
          nSel = nSel + 1
          h = xrh(iRef)
          k = xrk(iRef)
          l = xrl(iRef)
C
          if (QHerm .and. RecHem(h, k, l) .eq. -1) then
            call WrnDie(-5, 'hklAsyPt', 'Fatal Coding Error')
            call Die
          end if
C
          hA  = h
          kA  = k
          lA  = l
          HtA = 0
          sA  = 0
C
          do iSym = 1, nSymMx
            call HmulMt(h, k, l,
     &                  mSymMx, SymMx, iSym,
     &                  hM, kM, lM, Ht, STBF)
C
            if (QHerm .and. RecHem(hM, kM, lM) .eq. -1) then
              hM = -hM
              kM = -kM
              lM = -lM
              sM = -1
            else
              sM =  1
            end if
C
            if (CmpHKL(hA, kA, lA, hM, kM, lM) .gt. 0) then
              hA  = hM
              kA  = kM
              lA  = lM
              HtA = Ht
              sA  = sM
            end if
          end do
C
          AsyPt(1, iRef) = 0
          AsyPt(2, iRef) = 0
C
          if (sA .eq. 0) then
            nAsy = nAsy + 1
          else
            nSym = nSym + 1
C
                      jRef = 1
            do while (jRef .le. nRef .and. AsyPt(1, iRef) .eq. 0)
              if (      TSel(jRef) .ne. 0
     &            .and.  xrh(jRef) .eq. hA
     &            .and.  xrk(jRef) .eq. kA
     &            .and.  xrl(jRef) .eq. lA) AsyPt(1, iRef) = sA * jRef
              jRef = jRef + 1
            end do
C
            if (AsyPt(1, iRef) .eq. 0) then
              call WrnDie(-5, 'hklAsyPt',
     &                    'Corrupt Symmetry Operations.')
              call Die
            end if
C
            AsyPt(2, iRef) = -sA * HtA
          end if
        end if
      end do
C
      if (WrnLev .ge. 5) then
        write(6,'(1X,A,I7)') 'TSMAP: Total number of reflections', nRef
        write(6,'(1X,A,I7)') 'TSMAP: Selected reflections       ', nSel
        write(6,'(1X,A,I7)') 'TSMAP: Reflections in asym. unit  ', nAsy
        write(6,'(1X,A,I7)') 'TSMAP: Sym. equiv. reflections    ', nSym
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine SetLiRef3D(nRef, xrh, xrk, xrl, TSel, QHerm,
     &                      hfcg11, hfcg12,
     &                      hfcg21, hfcg22,
     &                      hfcg31, hfcg32,
     &                      LiRef3D)
      IMPLICIT NONE
C I/O
      integer  nRef, xrh(*), xrk(*), xrl(*), TSel(*)
      logical  QHerm
      integer  hfcg11, hfcg12
      integer  hfcg21, hfcg22
      integer  hfcg31, hfcg32
      integer  LiRef3D(hfcg11:hfcg12,
     &                 hfcg21:hfcg22,
     &                 hfcg31:hfcg32)
C
C Store indices for 1D list of hkl in 3D array
C
C local
      integer  iRef, h, k, l
C
C externals
      integer   RecHem
      external  RecHem
C
C begin
      do iRef = 1, nRef
        if (TSel(iRef) .ne. 0) then
          h = xrh(iRef)
          k = xrk(iRef)
          l = xrl(iRef)
C
          if (QHerm .and. RecHem(h, k, l) .eq. -1) then
            call WrnDie(-5, 'SetLiRef3D', 'Fatal Coding Error.')
            call Die
          end if
C
          LiRef3D(h, k, l) = iRef
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine GetLiRef3D(hfcg11, hfcg12,
     &                      hfcg21, hfcg22,
     &                      hfcg31, hfcg32,
     &                      LiRef3D, h, k, l, iRef)
      IMPLICIT NONE
C I/O
      integer  hfcg11, hfcg12
      integer  hfcg21, hfcg22
      integer  hfcg31, hfcg32
      integer  LiRef3D(hfcg11:hfcg12,
     &                 hfcg21:hfcg22,
     &                 hfcg31:hfcg32)
      integer  h, k, l, iRef
C
C Get index in 1D list of hkl from 3D array
C
C begin
      iRef = LiRef3D(h, k, l)
C
      return
      end
C}
C=======================================================================
C{
      subroutine SetTjRef(nRef, xrh, xrk, xrl, TSel, QHerm,
     &                    mSymMx, SymMx, nSymMx, STBF,
     &                    TjRef, LessMemory)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      integer  nRef, xrh(*), xrk(*), xrl(*), TSel(*)
      logical  QHerm
      integer  mSymMx, SymMx(mSymMx, 3, 4), nSymMx, STBF
      integer  TjRef(nRef, *)
      logical  LessMemory
C
C Compute pointers to sym. equiv. reflections for use in fcAtX.
C
C local
      integer  h3D(3, 2), n3D(3), Pn3D
      integer  iSym, iRef, jRef
      integer  hM, kM, lM, Ht, sM
C
C pointers
      integer  ptLiRef3D
C
C externals
      integer   RecHem
      external  RecHem
C
C begin
      if (.not. LessMemory) then
        call DimH3D(nRef, xrh, xrk, xrl, TSel, QHerm, h3D, n3D)
        Pn3D = n3D(1) * n3D(2) * n3D(3)
        ptLiRef3D = AllHp(integ4(Pn3D))
        call Fill4(heap(ptLiRef3D), Pn3D, 0)
        call SetLiRef3D(nRef, xrh, xrk, xrl, TSel, QHerm,
     &                  h3D(1, 1), h3D(1, 2),
     &                  h3D(2, 1), h3D(2, 2),
     &                  h3D(3, 1), h3D(3, 2),
     &                  heap(ptLiRef3D))
      end if
C
      do iSym = 1, nSymMx
        do iRef = 1, nRef
          if (TSel(iRef) .ne. 0) then
            call HmulMt(xrh(iRef), xrk(iRef), xrl(iRef),
     &                  mSymMx, SymMx, iSym,
     &                  hM, kM, lM, Ht, STBF)
C
            if (QHerm .and. RecHem(hM, kM, lM) .eq. -1) then
              hM = -hM
              kM = -kM
              lM = -lM
              sM = -1
            else
              sM =  1
            end if
C
            if (.not. LessMemory) then
              call GetLiRef3D(h3D(1, 1), h3D(1, 2),
     &                        h3D(2, 1), h3D(2, 2),
     &                        h3D(3, 1), h3D(3, 2),
     &                        heap(ptLiRef3D), hM, kM, lM, jRef)
              TjRef(iRef, iSym) = sM * jRef
            else
              TjRef(iRef, iSym) = 0
C
                        jRef = 1
              do while (jRef .le. nRef .and. TjRef(iRef, iSym) .eq. 0)
                if (      TSel(jRef) .ne. 0
     &              .and.  xrh(jRef) .eq. hM
     &              .and.  xrk(jRef) .eq. kM
     &              .and.  xrl(jRef) .eq. lM) then
                  TjRef(iRef, iSym) = sM * jRef
                end if
                jRef = jRef + 1
              end do
            end if
C
            if (TjRef(iRef, iSym) .eq. 0) then
              call WrnDie(-5, 'SetTjRef', 'Fatal Coding Error.')
              call Die
            end if
          end if
        end do
      end do
C
      if (.not. LessMemory)
     &  call FreHp(ptLiRef3D, integ4(Pn3D))
C
      return
      end
C}
C=======================================================================
C{
      subroutine fcAtX(nRef, xrh, xrk, xrl, TSel, fcalc, HAsyPt,
     &                 TjRef, SymMx, mSymMx, nSymMx, STBF,
     &                 X, SpFac, fcX)
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      integer           nRef, xrh(*), xrk(*), xrl(*), TSel(*)
      double complex    fcalc(*)
      integer           HAsyPt(2, *), TjRef(nRef, *)
      integer           mSymMx, SymMx(mSymMx, 3, 4), nSymMx, STBF
      double precision  X(3)
      integer           SpFac
      double complex    fcX(*)
C
C Compute Fcalc at grid point X (fractional coordinates).
C
C local
      integer           iSym, iRef, jRef, HMs(3)
      double precision  Ts(3), eiArg
      double complex    fHMs
C
C parameter
      double complex    CZero
      double precision  Zero
      parameter(Zero = 0.0D0)
C
C externals
      logical   IsUTr
      external  IsUTr
C
C begin
      CZero = dcmplx(Zero, Zero)
C
      call FillC8(fcX, nRef, CZero)
C
      do iSym = 1, nSymMx
        Ts(1) = SymMx(iSym, 1, 4) / dfloat(STBF)
        Ts(2) = SymMx(iSym, 2, 4) / dfloat(STBF)
        Ts(3) = SymMx(iSym, 3, 4) / dfloat(STBF)
C
        do iRef = 1, nRef
          if (      TSel(iRef) .ne. 0
     &        .and. HAsyPt(1, iRef) .eq. 0) then
                jRef = TjRef(iRef, iSym)
            if (jRef .gt. 0) then
              fHMs    =        fcalc( jRef)
               HMs(1) =          xrh( jRef)
               HMs(2) =          xrk( jRef)
               HMs(3) =          xrl( jRef)
            else
              fHMs    = dconjg(fcalc(-jRef))
               HMs(1) =         -xrh(-jRef)
               HMs(2) =         -xrk(-jRef)
               HMs(3) =         -xrl(-jRef)
            end if
C
            eiArg = (  xrh(iRef) * Ts(1) + HMs(1) * X(1)
     &               + xrk(iRef) * Ts(2) + HMs(2) * X(2)
     &               + xrl(iRef) * Ts(3) + HMs(3) * X(3)) * 2 * PI
C
            fcX(iRef) = fcX(iRef)
     &        + fHMs * dcmplx(dcos(eiArg), dsin(eiArg)) / SpFac
          end if
        end do
      end do
C
C Do "internal expand", because fcX are only computed in the
C asymmetric unit (see HAsyPt).
      do iRef = 1, nRef
        jRef = HAsyPt(1, iRef)
        if (jRef .ne. 0 .and. TSel(iRef) .ne. 0) then
C phase(H*Ms) = phase(H) - 2*pi*H*Ts
          eiArg = -2 * PI * HAsyPt(2, iRef) / STBF
          fcX(iRef) = fcX(abs(jRef)) * dcmplx(dcos(eiArg), dsin(eiArg))
          if (jRef .lt. 0) fcX(iRef) = dconjg(fcX(iRef))
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine XSpFac(X, SpFac, SymMx, mSymMx, nSymMx, STBF)
      IMPLICIT NONE
C I/O
      double precision  X(3)
      integer           SpFac
      integer           mSymMx, SymMx(mSymMx, 3, 4), nSymMx, STBF
C
C Compute "SpFac" for correct scaling of Fcalc if X (fractional
C coordinates) is a special positions.
C
C local
      integer           ir, iSym
      double precision  Delta(3)
C
C externals
      logical   IsUTr
      external  IsUTr
C
C begin
      SpFac = 0
C
      do iSym = 1, nSymMx
        do ir = 1, 3
          Delta(ir) =   SymMx(iSym, ir, 1) * X(1)
     &                + SymMx(iSym, ir, 2) * X(2)
     &                + SymMx(iSym, ir, 3) * X(3)
     &                + SymMx(iSym, ir, 4) / dfloat(STBF)
     &                - X(ir)
        end do
        if (      IsUTr(Delta(1))
     &      .and. IsUTr(Delta(2))
     &      .and. IsUTr(Delta(3))) SpFac = SpFac + 1
      end do
C
      if (SpFac .eq. 0 .or. mod(nSymMx, max(SpFac, 1)) .ne. 0) then
        call WrnDie(-5, 'TSMAP', 'Fatal Coding Error.')
        call Die
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine ConvTS(nRef, ptxrh, ptxrk, ptxrl, TSel, ptTSel,
     &                  ptP1Fcalc, ptTrFcalc,
     &                  QHerm,
     &                  SymMx, ITSymMx, mSymMx, nSymMx, STBF,
     &                  Rn1, Rn2, Rn3, rRho, FlagMap,
     &                  QSpecPos, QSpecOnly,
     &                  ShowProgress, nIndFM,
     &                  LessMemory,
     &                  xSFNum, xSFNam, xSFType, hpSF,
     &                  XRE,
     &                  MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &                  XRTR,XRSCAL,XCVTEST,
     &                  TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &                  TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &                  HPMULT,HPTYPE,XRMREF,XRCELL,XRVOL)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      integer  nRef
      integer  ptxrh, ptxrk, ptxrl
      integer  TSel(*)
      integer  ptTSel
      integer  ptP1Fcalc, ptTrFcalc
      logical  QHerm
      integer  mSymMx, nSymMx, STBF
      integer  SymMx(mSymMx, 3, 4), ITSymMx(mSymMx, 3, 3)
      integer  Rn1, Rn2, Rn3
      real     rRho(*)
      integer  FlagMap(*)
      logical  QSpecPos, QSpecOnly
      logical  ShowProgress
      integer  nIndFM
      logical  LessMemory
C
      integer        xSFNum
      character*(*)  xSFNam(*), xSFType(*)
      integer        hpSF(*)
C
      DOUBLE PRECISION  XRE
      INTEGER MBINS
      DOUBLE PRECISION XBINLOW, XBINHIGH, BINSHELL(*)
      DOUBLE PRECISION XRTR(3,3), XRSCAL
      LOGICAL XCVTEST
      INTEGER TRPNMX, TRPNX
      CHARACTER*(*) TRPN(4,TRPNMX)
      INTEGER TRPNL(4,TRPNMX), TRPNN
      DOUBLE COMPLEX TRPNDB(4,TRPNMX)
      INTEGER TRPNMLT(TRPNMX)
      CHARACTER*2 TRPNTYP(TRPNMX), TRPNDOM(TRPNMX)
      INTEGER TDEPTH
      INTEGER TRPNLEV(TRPNMX)
      INTEGER HPMULT,HPTYPE,XRMREF
      DOUBLE PRECISION XRCELL(3,3), XRVOL
C
C Main procedure Conventional (Direct) Translation Search.
C
C local
      logical           QBuf
      integer           iMap, jx, jy, jz
      integer           SpFac
      double precision  X(3)
C
      integer           nDone
      double precision  tmEntr, tmExit, tmStart, tmNow, tmWrite, tmWait
C
C pointer
      integer  ptHAsyPt, ptTjRef
      INTEGER  XDERIV
C
C parameters
      double complex    CZero
      double precision  Zero
      parameter(Zero = 0.0D0)
C
C begin
      CZero = dcmplx(Zero, Zero)
C
      if (ShowProgress .or. Timer .gt. 0) call VCPU(tmEntr)
C
      QBuf = .false.
C
      if (ptP1Fcalc .eq. ptTrFcalc) then
        QBuf = .true.
        ptP1Fcalc = AllHp(icplx8(nRef))
        call CopyC8(heap(ptTrFcalc), heap(ptP1Fcalc), nRef)
      end if
C
      ptHAsyPt = AllHp(integ4(2 * nRef))
      ptTjRef  = AllHp(integ4(nRef * nSymMx))
C
      call hklAsyPt(nRef, heap(ptxrh), heap(ptxrk), heap(ptxrl),
     &              TSel, heap(ptHAsyPt), QHerm,
     &              mSymMx, SymMx, nSymMx, STBF,
     &              WrnLev)
C
      if (ShowProgress .or. Timer .gt. 0) then
        call VCPU(tmStart)
        write(6, '(1X, A)')
     &    'TSMAP: Preparing direct translation search...'
      end if
C
      call SetTjRef(nRef, heap(ptxrh), heap(ptxrk), heap(ptxrl),
     &              TSel, QHerm,
     &              mSymMx, SymMx, nSymMx, STBF,
     &              heap(ptTjRef), LessMemory)
C
      if (ShowProgress .or. Timer .gt. 0) then
        call VCPU(tmNow)
        write(6, '(1X, A, F10.4, A)')
     &    'TSMAP: CPU-time: ', tmNow - tmStart,
     &    ' seconds Preparing direct translation search'
C
        nDone = 0
        tmStart = tmNow
        tmWrite = tmNow
        tmWait = dfloat(60)
      end if
C
      iMap = 1
      SpFac = 1
C
      do jz = 0, Rn3-1
        X(3) = jz / dfloat(Rn3)
      do jy = 0, Rn2-1
        X(2) = jy / dfloat(Rn2)
      do jx = 0, Rn1-1
        X(1) = jx / dfloat(Rn1)
C
        if (QSpecPos)
     &    call XSpFac(X, SpFac, SymMx, mSymMx, nSymMx, STBF)
C
        if (.not. QSpecOnly .or. SpFac .ne. 1) then
          if (FlagMap(iMap) .le. 0) then
            call fcAtX(nRef, heap(ptxrh), heap(ptxrk), heap(ptxrl),
     &                 TSel, heap(ptP1Fcalc),
     &                 heap(ptHAsyPt), heap(ptTjRef),
     &                 SymMx, mSymMx, nSymMx, STBF,
     &                 X, SpFac, heap(ptTrFcalc))
C
C Call target function. No derivatives. No monitor calculation.
C Working set only. Using dummy parameters for dtarget and monitor
C expression.
            CALL XTARGETS(.FALSE.,.FALSE.,1,XRE,XDERIV,
     &                    ptTSel, nRef,
     &                    MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &                    XRTR, XRSCAL,XCVTEST,
     &                    TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &                    TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &                    1,1,' ',1,0,CZERO,1,1,' ',' ',1,
     &                    1,1,' ',1,0,CZERO,1,1,' ',' ',1,
     &                    ptxrh, ptxrk, ptxrl, xSFNum, xSFNam,
     &                    xSFType, hpSF, HPMULT,HPTYPE,XRMREF,
     &                    QHerm,
     &                    nSymMx, mSymMx, STBF, SymMx, ITSymMx,0,
     &                    XRCELL,XRVOL)
C
            if (ShowProgress) then
              nDone = nDone + 1
              call VCPU(tmNow)
              if (tmNow - tmWrite .ge. tmWait) then
                write(6, '(1X, A, F10.4, A, I11)')
     &            'TSMAP: CPU-time: ', tmNow - tmStart,
     &            ' seconds: Grid points done: ', nDone
                write(6, '(1X, A, I11, A, F10.0, A)')
     &            'TSMAP: Estimated CPU-time for all ',
     &            nIndFM, ' grid points: ',
     &            (tmNow - tmStart) / nDone * nIndFM, ' seconds'
                tmWrite = tmNow
                tmWait = dfloat(3600)
              end if
            end if
C
            rRho(iMap) = 1 - XRE/XRSCAL
          end if
        end if
        iMap = iMap + 1
      end do
      end do
      end do
C
      do iMap = 1, Rn1 * Rn2 * Rn3
        if (FlagMap(iMap) .gt. 0) then
          rRho(iMap) = rRho(FlagMap(iMap))
        end if
      end do
C
      call FreHp(ptTjRef,  integ4(nRef * nSymMx))
      call FreHp(ptHAsyPt, integ4(2 * nRef))
C
      if (QBuf) then
        call CopyC8(heap(ptP1Fcalc), heap(ptTrFcalc), nRef)
        call FreHp(ptP1Fcalc, icplx8(nRef))
        ptP1Fcalc = ptTrFcalc
      end if
C
      if (ShowProgress .or. Timer .gt. 0) then
        call VCPU(tmExit)
        write(6, '(1X, A, F10.4, A)')
     &    'TSMAP: CPU-time: ', tmExit - tmEntr,
     &    ' seconds Direct Translation Search'
      end if
C
      return
      end
C}
C=======================================================================
C{
      integer function H2iH(H, NH, FullRange)
      IMPLICIT NONE
C I/O
      integer  H, NH
      logical  FullRange
C
C Return array index for h or k or l (1-dimensional)
C
C local
      integer  m
C
C begin
      H2iH = H
C
      if (FullRange) then
        m = (NH - 1) / 2
C
        if (-m .gt. H .or. H .gt. m) then
          H2iH = -1
        else if (H .lt. 0) then
          H2iH = H2iH + NH
        end if
      else
        if (0 .gt. H .or. H .ge. NH) then
          H2iH = -1
        end if
      end if
C
      return
      end
C}
C=======================================================================
C{
      logical function hkl2ihkl(h, k, l, nh, nk, nl, QHerm, ih, ik, il)
      IMPLICIT NONE
C I/O
      integer  h, k, l, nh, nk, nl
      logical  QHerm
      integer  ih, ik, il
C
C Compute array indices for hkl (3-dimensional)
C If QHerm is .true., reflections with l < 0 are mapped onto their
C conjugate complex counterpart with l > 0. In this case, hkl2ihkl
C returns .true. to indicate the mapping.
C
C externals
      integer   H2iH
      external  H2iH
C
C begin
      hkl2ihkl = .false.
C
      if (QHerm) then
        if (l .ge. 0) then
          ih = H2iH( h, nh, .true.)
          ik = H2iH( k, nk, .true.)
          il = H2iH( l, nl, .false.)
        else
          ih = H2iH(-h, nh, .true.)
          ik = H2iH(-k, nk, .true.)
          il = H2iH(-l, nl, .false.)
          hkl2ihkl = .true.
        end if
      else
        ih = H2iH(h, nh, .true.)
        ik = H2iH(k, nk, .true.)
        il = H2iH(l, nl, .true.)
      end if
C
      return
      end
C}
C=======================================================================
C{
      logical function hkl2Wi(h, k, l, Wn, Wi)
      IMPLICIT NONE
C I/O
      INCLUDE 'timer.inc'
      integer  h, k, l, Wn(3), Wi(3)
C
C Return index in TabH for hkl
C   I:  Wm = physical dimensions of work array
C   I:  Wn = natural  dimensions, Wn(i) <= Wm(i)
C
C local
      integer  i, m
C
C externals
      integer   H2iH
      external  H2iH
C
C begin
      hkl2Wi = .false.
C
        Wi(1) = H2iH( h, Wn(1) / 2, .false.)
      if (Wi(1) .ge. 0) then
        Wi(2) = H2iH( k, Wn(2),     .true.)
        Wi(3) = H2iH( l, Wn(3),     .true.)
C
        do i = 1, 3
          if (Wi(i) .lt. 0) then
            write(6, '(1X, A)')
     &        'TSMAP: Fatal Error:',
     &        ' Workspace resolution is not sufficient.'
            call WrnDie(-5, 'TSMAP',
     &           'Set of reflections not consistent with'
     &        // ' local TSMAp symmetry.')
            call Die
          end if
        end do
C
        if (WrnLev .ge. 15) then
          m = 2
          do i = 1, 3
            if (m * Wi(i) .ge. Wn(i)) then
              call WrnDie(-5, 'hkl2Wi', 'Fatal Coding Error.')
              call Die
            end if
            m = 1
          end do
        end if
C
        hkl2Wi = .true.
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine ftilSet(h, k, l,
     &                   mSymMx, SymMx, iSym,
     &                   nh, nk, nl, fcg, QHerm,
     &                   hM, kM, lM, Ht, STBF, ftil, NeedAll)
      IMPLICIT NONE
C I/O
      integer         h, k, l
      integer         mSymMx, SymMx(mSymMx, 3, 4), iSym
      integer         nh, nk, nl
      double complex  fcg(0:nh-1, 0:nk-1, 0:*)
      logical         QHerm
      integer         hM, kM, lM, Ht, STBF
      double complex  ftil
      logical         NeedAll
C
C Compute f~(hs) (Navaza paper, p. 447, right after Eq. (14))
C for given hM, kM, lM, Ht:
C
C   ftil = fcg(hM, kM, lM) * exp(2*pi*i*Ht)
C
C local
      logical         IsConjg
      integer         ih, ik, il
      double complex  rot
C
C parameters
      double precision  Zero
      parameter(Zero = 0.0D0)
C
C externals
      logical   hkl2ihkl
      external  hkl2ihkl
C
C begin
      call HmulMt(h, k, l,
     &            mSymMx, SymMx, iSym,
     &            hM, kM, lM, Ht, STBF)
C
      IsConjg = hkl2ihkl(hM, kM, lM, nh, nk, nl, QHerm, ih, ik, il)
C
      if (ih .ge. nh .or. ik .ge. nk .or. il .ge. nl) then
        call WrnDie(-5, 'ftilSet', 'Fatal Coding Error.')
        call Die
      end if
C
      if (ih .lt. 0 .or. ik .lt. 0 .or. il .lt. 0) then
        if (NeedAll) then
          call WrnDie(-5, 'ftilSet', 'Fatal Coding Error.')
          call Die
        end if
        ftil = cmplx(Zero, Zero)
      else
        ftil = fcg(ih, ik, il)
        if (IsConjg) ftil = dconjg(ftil)
        call e2PIix(dfloat(Ht) / dfloat(STBF), rot)
        ftil = ftil * rot
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine ftsCw1(Wm1, Wm2, ws1,
     &                  SmH, SmHdI2, SmHFF, SmHF2F2,
     &                  Rn1, Rn2, Rn3, rRho)
      IMPLICIT NONE
C I/O
      integer           Wm1, Wm2
      real              ws1(0:Wm1-1, 0:Wm2-1, 0:*)
      integer           SmH
      double precision  SmHdI2, SmHFF, SmHF2F2
      integer           Rn1, Rn2, Rn3
      real              rRho(0:Rn1-1, 0:Rn2-1, 0:*)
C
C Combine result of first (ws3, stored in rRho) and second (ws1)
C transform (Navaza paper, p. 447, Eq. (12)). Result is in rRho.
C   I: ws1      = result from FFT
C   I: Wm       = physical dimensions of ws1
C   I: Sm*      = various sums (see ftsSmH)
C   +: rRho     = real (as opposed to imaginary) part of real space object
C   I: Rn       = natural (and also physical) dimensions of rRho
C
C local
      integer           jx, jy, jz
      double precision  ws1p, ws3p, d
C
C parameters
      double precision  Zero
      parameter(Zero = 0.0D0)
C
C begin
      do jz = 0, Rn3-1
      do jy = 0, Rn2-1
      do jx = 0, Rn1-1
        if (SmH .eq. 0) then
          rRho(jx, jy, jz) = Zero
        else
C Navaza paper, p. 447, Eq. (13)
          ws1p =  ws1(jx, jy, jz) + SmHFF
          ws3p = rRho(jx, jy, jz) + SmHF2F2
          d = ws3p - ws1p * ws1p / dfloat(SmH)
          if (d .gt. Zero) then
            rRho(jx, jy, jz) = sqrt(d) * sqrt(SmHdI2)
          else
            rRho(jx, jy, jz) = Zero
          end if
        end if
      end do
      end do
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine ftsCw2(Wm1, Wm2, ws2,
     &                  SmHdIFF,
     &                  Rn1, Rn2, Rn3, rRho)
      IMPLICIT NONE
C I/O
      integer           Wm1, Wm2
      real              ws2(0:Wm1-1, 0:Wm2-1, 0:*)
      double precision  SmHdIFF
      integer           Rn1, Rn2, Rn3
      real              rRho(0:Rn1-1, 0:Rn2-1, 0:*)
C
C Combine result of third transform (ws2) and previous results stored
C in rRho (Navaza paper, p. 447, Eq. (12)). Result is again in rRho.
C   I: ws1      = result from FFT
C   I: Wm       = physical dimensions of ws1
C   I: SmHdIFF  = see ftsSmH
C   +: rRho     = real      part of real space object
C   I: Rn       = natural (and also physical) dimensions of rRho
C
C local
      integer           jx, jy, jz
      double precision  ws2p
C
C parameters
      double precision  Zero, BigCC
      parameter(Zero = 0.0D0, BigCC = 1.0D6)
C
C begin
      do jz = 0, Rn3-1
      do jy = 0, Rn2-1
      do jx = 0, Rn1-1
        ws2p = ws2(jx, jy, jz) + SmHdIFF
        if (abs(ws2p / BigCC) .lt. rRho(jx, jy, jz)) then
C Navaza paper, p. 447, Eq. (12)
          rRho(jx, jy, jz) = ws2p / rRho(jx, jy, jz)
        else
          rRho(jx, jy, jz) = Zero
        end if
      end do
      end do
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine fcgSet(xrh, xrk, xrl, fcalc, TSel,
     &                  nRef, QHerm,
     &                  nh, nk, nl, fcg)
      IMPLICIT NONE
C I/O
      integer         xrh(*), xrk(*), xrl(*)
      double complex  fcalc(*)
      integer         TSel(*), nRef
      logical         QHerm
      integer         nh, nk, nl
      double complex  fcg(0:nh-1, 0:nk-1, 0:*)
C
C Copy fcalc to 3-dimensional fcalc work array (fcg)
C
C local
      logical  IsConjg
      integer  ih, ik, il
      integer  iRef
C
C parameter
      double complex    CZero
      double precision  Zero
      parameter(Zero = 0.0D0)
C
C externals
      logical   hkl2ihkl
      external  hkl2ihkl
C
C begin
      CZero = dcmplx(Zero, Zero)
C
      call FillC8(fcg, nh * nk * nl, CZero)
C
      do iRef = 1, nRef
        if (TSel(iRef) .ne. 0) then
          IsConjg = hkl2ihkl(xrh(iRef), xrk(iRef), xrl(iRef),
     &                       nh, nk, nl, QHerm,
     &                       ih, ik, il)
C
          if (     IsConjg
     &        .or. ih .lt. 0  .or. ik .lt.  0 .or. il .lt.  0
     &        .or. ih .ge. nh .or. ik .ge. nk .or. il .ge. nl) then
            call WrnDie(-5, 'fcgSet', 'Fatal Coding Error.')
            call Die
          end if
C
          fcg(ih, ik, il) = fcalc(iRef)
C
          if (QHerm .and. xrl(iRef) .eq. 0) then
C If reflection is in hk0 plane, also set the (redundant)
C complex-conjugate Friedel mate. This simplifies the use of
C fcg in the actual summation.
            IsConjg = hkl2ihkl(-xrh(iRef), -xrk(iRef), 0,
     &                         nh, nk, nl, QHerm,
     &                         ih, ik, il)
C
            if (     IsConjg
     &          .or. ih .lt. 0  .or. ik .lt.  0 .or. il .lt.  0
     &          .or. ih .ge. nh .or. ik .ge. nk .or. il .ge. nl) then
              call WrnDie(-5, 'fcgSet', 'Fatal Coding Error.')
              call Die
            end if
C
            fcg(ih, ik, il) = dconjg(fcalc(iRef))
          end if
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine SetfFFTgrid(Lbl, fFFTgrid, hfcg, Rn, Fac,
     &                       FTBase, FTCRBx, FTPrim)
      IMPLICIT NONE
C I/O
      character  Lbl*(*)
      integer    fFFTgrid(3), hfcg(3), Rn(3), Fac
      integer    FTBase, FTCRBx, FTPrim
C
C Compute factor by which work array has to be larger
C   I: Lbl      = label (variable name) for warning message
C   O: fFFTgrid = the computed factors
C   I: hfcg     = max. hkl to be stored in work array
C   I: Rn       = dimensions of real space object for final storage
C   I: Fac      = resolution of work space / res. of real space obj.
C   I: FTBase   = smallest prime allowed for FFT array dimensions
C   I: FTCRBx   = smallest prime allowed for FFT array dimension x
C   I: FTPrim   = largest  prime allowed for FFT array dimensions
C
C local
      integer  H, j
      integer  xFFTgrid, Base
C
C externals
      integer   iLCM
      external  iLCM
C
C begin
      do j = 1, 3
        H = hfcg(j) * 2 * Fac + 1
        fFFTgrid(j) = H / Rn(j)
        if (Rn(j) * fFFTgrid(j) .lt. H)
     &    fFFTgrid(j) = fFFTgrid(j) + 1
        call xprime(fFFTgrid(j), Rn(j), FTBase, FTPrim)
      end do
C
      if (FTCRBx .gt. 1) then
        xFFTgrid = fFFTgrid(1)
        Base = iLCM(FTBase, FTCRBx)
        call xprime(xFFTgrid, Rn(1), Base, FTPrim)
C
        if (xFFTgrid .ne. fFFTgrid(1)) then
          write(6, '(1X, 2A)')
     &      'TSMAP-WRN: Constraint on grid size for',
     &      ' FFT array along x forces'
          write(6, '(1X, 3A, I3, A, I3, A)')
     &      'TSMAP-WRN: ', Lbl, '(x) to be increased from ',
     &      fFFTgrid(1), ' to ', xFFTgrid, '.'
          write(6, '(1X, 2A)')
     &      'TSMAP-WRN: This will cause artificially',
     &      ' increased memory allocation.'
          write(6, '(1X, A, /, 1X, A, I3, A)')
     &      'TSMAP-WRN: A better solution is to set',
     &      'TSMAP-WRN: XRAY FFT XGRIdfactor to a multiple of ',
     &      Base, '.'
C
          fFFTgrid(1) = xFFTgrid
        end if
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine ftsMult(nRef, xrh, xrk, xrl, Mult, TSel, QHerm,
     &                   mSymMx, SymMx, nSymMx,
     &                   Verbose, WrnLev)
      IMPLICIT NONE
C I/O
      integer  nRef, xrh(*), xrk(*), xrl(*), Mult(*), TSel(*)
      logical  QHerm
      integer  mSymMx, SymMx(mSymMx, 3, 4), nSymMx
      logical  Verbose
      integer  WrnLev
C
C Compute Mult(iplicity)
C
C local
      integer  iRef, h, k, l
      integer  iSym, hM, kM, lM, Ht
      integer  R, m, c
      integer  nSel, nUse, nSym
C
C externals
      integer   CmpHKL, RecHem
      external  CmpHKL, RecHem
C
C begin
      nSel = 0
      nUse = 0
      nSym = 0
C
      do iRef = 1, nRef
        Mult(iRef) = 0
C
        if (TSel(iRef) .ne. 0) then
          nSel = nSel + 1
C
          h = xrh(iRef)
          k = xrk(iRef)
          l = xrl(iRef)
C
          if (QHerm .and. RecHem(h, k, l) .eq. -1) then
            call WrnDie(-5, 'ftsMult', 'Fatal Coding Error')
            call Die
          end if
C
          R = 0
          m = 0
C
          do iSym = 1, nSymMx
            call HmulMt(h, k, l, mSymMx, SymMx, iSym, hM, kM, lM, Ht, 1)
C
            if (QHerm .and. RecHem(hM, kM, lM) .eq. -1) then
              hM = -hM
              kM = -kM
              lM = -lM
            end if
C
            c = CmpHKL(h, k, l, hM, kM, lM)
C
            if      (c .eq. 0) then
              m = m + 1
            else if (c .gt. 0) then
              nSym = nSym + 1
              goto 999
            else if (-h .eq. hM .and. -k .eq. kM .and. -l .eq. lM) then
              R = 1
            end if
          end do
C
          if (m .eq. 0) then
            call WrnDie(-5, 'ftsMult', 'Corrupt Symmetry Operations.')
            call Die
          end if
          if (mod(nSymMx, m) .ne. 0) then
            call WrnDie(-5, 'ftsMult', 'Corrupt Symmetry Operations.')
            call Die
          end if
C
          Mult(iRef) = nSymMx / m
C
          if (QHerm .and. R .eq. 0)
     &      Mult(iRef) = Mult(iRef) * 2
C
          nUse = nUse + 1
C
 999      continue
        end if
      end do
C
      if (Verbose .or. WrnLev .ge. 5) then
        write(6,'(1X,A,I7)') 'TSMAP: Total number of reflections', nRef
        write(6,'(1X,A,I7)') 'TSMAP: Selected reflections       ', nSel
        write(6,'(1X,A,I7)') 'TSMAP: Used reflections           ', nUse
        write(6,'(1X,A,I7)') 'TSMAP: Sym. equiv. reflections    ', nSym
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine ftsCheckMult(nRef, xrh, xrk, xrl, TSel, Mult,
     &                        QHerm, mSymMx, SymMx, nSymMx)
      IMPLICIT NONE
C I/O
      integer  nRef, xrh(*), xrk(*), xrl(*), TSel(*), Mult(*)
      logical  QHerm
      integer  mSymMx, SymMx(mSymMx, 3, 4), nSymMx
C
C This subroutine simply checks if the result of ftsMult is consistent.
C Only called if MESSage = DEBUg.
C
C local
      integer  iRef, jRef, h, k, l
      integer  iSym, hM, kM, lM, Ht
C
C externals
      integer   RecHem
      external  RecHem
C
C begin
      do iRef = 1, nRef
        if (TSel(iRef) .ne. 0 .and. Mult(iRef) .eq. 0) then
          h = xrh(iRef)
          k = xrk(iRef)
          l = xrl(iRef)
C
          do iSym = 2, nSymMx
            call HmulMt(h, k, l, mSymMx, SymMx, iSym, hM, kM, lM, Ht, 1)
C
            if (QHerm .and. RecHem(hM, kM, lM) .eq. -1) then
              hM = -hM
              kM = -kM
              lM = -lM
            end if
C
            do jRef = 1, nRef
              if (     Mult(jRef) .gt. 0
     &            .and. xrh(jRef) .eq. hM
     &            .and. xrk(jRef) .eq. kM
     &            .and. xrl(jRef) .eq. lM) goto 999
            end do
          end do
C
          write(6, '(1X, A, 3(1X, I4))')
     &      'ERROR-DIAG: hkl  =', h, k, l
          write(6, '(1X, A, 3(1X, I4))')
     &      'ERROR-DIAG: hklM =', hM, kM, lM
C
          call WrnDie(-5, 'ftsCheckMult', 'Fatal Coding Error.')
          call Die
C
 999      continue
        end if
      end do
C
      write(6, '(1X, A)') 'TSMAP: CheckMult OK.'
C
      return
      end
C}
C=======================================================================
C{
      subroutine ftsSmH(Mult, fobs, fpart, dIobs, nRef, QFpart,
     &                  SmH, SmHdI2, SmHFF, SmHdIFF, SmHF2F2)
      IMPLICIT NONE
C I/O
      integer           Mult(*)
      double precision  fobs(*)
      double complex    fpart(*)
      double precision  dIobs(*)
      integer           nRef
      logical           QFpart
      integer           SmH
      double precision  SmHdI2, SmHFF, SmHdIFF, SmHF2F2
C
C Compute various sums for use in FFT Translations Search.
C
C local
      integer           iRef, mH
      double precision  MeanIoH, ImH
C
C parameters
      double precision  Zero
      parameter(Zero = 0.0D0)
C
C begin
      SmH     = 0
      MeanIoH = Zero
C
      do iRef = 1, nRef
        mH = Mult(iRef)
        if (mH .eq. 0) then
          dIobs(iRef) = dfloat(-1)
        else
          dIobs(iRef) = fobs(iRef)**2
          SmH     = SmH     + mH
          MeanIoH = MeanIoH + mH * dIobs(iRef)
        end if
      end do
C
      if (SmH .ne. 0) MeanIoH = MeanIoH / SmH
C
      SmHdI2  = Zero
      SmHFF   = Zero
      SmHdIFF = Zero
      SmHF2F2 = Zero
C
      do iRef = 1, nRef
        mH = Mult(iRef)
        if (mH .gt. 0) then
          dIobs(iRef) = dIobs(iRef) - MeanIoH
          SmHdI2 = SmHdI2 + mH * dIobs(iRef) * dIobs(iRef)
C
          if (QFpart) then
            ImH = DBLE(fpart(iRef))**2 + dimag(fpart(iRef))**2
            SmHFF   = SmHFF   + mH * ImH
            SmHdIFF = SmHdIFF + mH * ImH * dIobs(iRef)
            SmHF2F2 = SmHF2F2 + mH * ImH * ImH
          end if
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine ftsSetXm(Rn, fFFTgrid, Xn, Xm, PXm, FTAvoi,
     &                    Msg, Verbose, WrnLev)
      IMPLICIT NONE
C I/O
      integer    Rn(3), fFFTgrid(3), Xn(3), Xm(3), PXm, FTAvoi
      character  Msg*(*)
      logical    Verbose
      integer    WrnLev
C
C Set physical dimensions of work array
C
C local
      integer  j
C
C begin
      do j = 1, 3
        Xn(j) = Rn(j) * fFFTgrid(j)
      end do
C
      call mFFTcr(3, Xn, FTAvoi, Xm)
C
      PXm = Xm(1) * Xm(2) * Xm(3)
C
      if ((WrnLev .ge. 5 .or. Verbose) .and. Msg .ne. ' ')
     &  write(6, '(1X, 3A, 3(1X, I4, A, I4, A))')
     &    'TSMAP: Dimensions of ', Msg, ' = ',
     &    Xn(1), '(', Xm(1), ')',
     &    Xn(2), '(', Xm(2), ')',
     &    Xn(3), '(', Xm(3), ')'
C
      return
      end
C}
C=======================================================================
C{
      subroutine ShrFFT(Rn1, Rn2, Rn3,
     &                  Vm1, Vm2, Vm3,
     &                  Vn1, Vn2, Vn3,
     &                  vFFTgrid, VVVisDbl, fred,
     &                  dVVV, fVVR, sVVV, VVR, VRR, RRR)
      IMPLICIT NONE
C I/O
      integer           Rn1, Rn2, Rn3
      integer           Vm1, Vm2, Vm3
      integer           Vn1, Vn2, Vn3
      integer           vFFTgrid(3)
      logical           VVVisDbl
      integer           fred
      double complex    dVVV(0:Vm1/2-1,      0:Vm2-1, 0:Vm3-1)
             complex    fVVR(0:Vm1/2*fred-1, 0:Vm2-1, 0:Rn3-1)
             complex    sVVV(0:Vm1/2-1,      0:Vm2-1, 0:Vm3-1)
             complex     VVR(0:Vm1/2-1,      0:Vm2-1, 0:Rn3-1)
             complex     VRR(0:Vm1/2-1,      0:Rn2-1, 0:Rn3-1)
      real               RRR(0:Rn1  -1,      0:Rn2-1, 0:Rn3-1)
C
C "Shrinking FFT" - 3D FFT with built-in down-sampling.
C Strategy:
C   1D FFT all rows along z, store results in array smaller by vFFTgrid(3)
C   1D FFT all rows along y, store results in array smaller by vFFTgrid(2)
C   1D FFT all rows along x, store results in array smaller by vFFTgrid(1)
C
C VVVisDbl = true:  Use dVVV (double precision) as starting array
C          = false: Use sVVV (single precision) as starting array
C
C All transforms are done in single precision.
C
C local
      integer  ix, iy, iz, i
      logical  FFTERR
      integer  mFFTw(3)
      integer  ptFFTw(3)
C
      integer   MaxSeq
      parameter(MaxSeq = 1024)
      complex  cSeq(0:  MaxSeq-1)
      real     rSeq(0:2*MaxSeq-1)
      equivalence(cSeq(0), rSeq(0))
C
C begin
      if (max(Vn1/2+1, Vn2, Vn3) .gt. MaxSeq) then
        call WrnDie(-5, 'ShrFFT',
     &              'MaxSeq too small --> recompile program')
        call Die
      end if
C
      do i = 1, 3
         mFFTw(i) = 0
        ptFFTw(i) = 0
      end do
C
      if (VVVisDbl) then
        do iy = 0, Vn2 - 1
        do ix = 0, Vn1 / 2
          do iz = 0, Vn3 - 1
            cSeq(iz) = dVVV(ix, iy, iz)
          end do
          call sfft1c(Vn3, cSeq, FFTERR, mFFTw(3), ptFFTw(3))
          if (FFTERR) then
            call WrnDie(-5, 'ShrFFT', 'Fatal: sfft1c failed.')
            call Die
          end if
          do iz = 0, Rn3 - 1
            fVVR(ix, iy, iz) = cSeq(iz * vFFTgrid(3))
          end do
        end do
        end do
      else
        do iy = 0, Vn2 - 1
        do ix = 0, Vn1 / 2
          do iz = 0, Vn3 - 1
            cSeq(iz) = sVVV(ix, iy, iz)
          end do
          call sfft1c(Vn3, cSeq, FFTERR, mFFTw(3), ptFFTw(3))
          if (FFTERR) then
            call WrnDie(-5, 'ShrFFT', 'Fatal: sfft1c failed.')
            call Die
          end if
          do iz = 0, Rn3 - 1
            VVR(ix, iy, iz) = cSeq(iz * vFFTgrid(3))
          end do
        end do
        end do
      end if
C
      do iz = 0, Rn3 - 1
        if (VVVisDbl .and. fred .ne. 1) then
          do iy = 0, Vn2 - 1
          do ix = 0, Vn1 / 2
            VVR(ix, iy, iz) = fVVR(ix, iy, iz)
          end do
          end do
        end if
C
        do ix = 0, Vn1 / 2
          do iy = 0, Vn2 - 1
            cSeq(iy) = VVR(ix, iy, iz)
          end do
          call sfft1c(Vn2, cSeq, FFTERR, mFFTw(2), ptFFTw(2))
          if (FFTERR) then
            call WrnDie(-5, 'ShrFFT', 'Fatal: sfft1c failed.')
            call Die
          end if
          do iy = 0, Rn2 - 1
            VRR(ix, iy, iz) = cSeq(iy * vFFTgrid(2))
          end do
        end do
C
        do iy = 0, Rn2 - 1
          do ix = 0, Vn1 / 2
            cSeq(ix) = VRR(ix, iy, iz)
          end do
          call sfft1cr(Vn1, cSeq, FFTERR, mFFTw(1), ptFFTw(1))
          if (FFTERR) then
            call WrnDie(-5, 'ShrFFT', 'Fatal: sfft1cr failed.')
            call Die
          end if
          do ix = 0, Rn1 - 1
            RRR(ix, iy, iz) = rSeq(ix * vFFTgrid(1))
          end do
        end do
      end do
C
      call sfft1cr(0, 0, FFTERR, mFFTw(1), ptFFTw(1))
      if (FFTERR) then
        call WrnDie(-5, 'ShrFFT', 'Fatal: sfft1cr failed.')
        call Die
      end if
      call sfft1c (0, 0, FFTERR, mFFTw(2), ptFFTw(2))
      if (FFTERR) then
        call WrnDie(-5, 'ShrFFT', 'Fatal: sfft1c failed.')
        call Die
      end if
      call sfft1c (0, 0, FFTERR, mFFTw(3), ptFFTw(3))
      if (FFTERR) then
        call WrnDie(-5, 'ShrFFT', 'Fatal: sfft1c failed.')
        call Die
      end if
C
      return
      end
C}
C======================================================================
C{
      subroutine FastTS(xrh, xrk, xrl, nRef, TSel, QHerm,
     &                  Fobs, Fcalc, Fpart, QFpart,
     &                  mSymMx, SymMx, nSymMx, STBF,
     &                  Rn, rRho,
     &                  FTBase, FTCRBx, FTPrim, FTAvoi, FPrec,
     &                  LessMemory, Verbose)
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'tsmap.inc'
      integer           xrh(*), xrk(*), xrl(*)
      integer           nRef, TSel(*)
      logical           QHerm
      double precision  Fobs(*)
      double complex    Fcalc(*), Fpart(*)
      logical           QFpart
      integer           mSymMx, SymMx(mSymMx, 3, 4), nSymMx, STBF
      integer           Rn(3)
      real              rRho(*)
      integer           FTBase, FTCRBx, FTPrim, FTAvoi
      character*(*)     FPrec
      logical           LessMemory, Verbose
C
C Main procedure FFT Translation Search.
C
C local
      integer  j, k, fred
      integer  uFFTgrid(3), vFFTgrid(3)
      integer  hfcg(3, 2), nfcg(3), Pnfcg
      integer  Um(3), Un(3), PUm
      integer  Vm(3), Vn(3), PVm
C
      integer           SmH
      double precision  SmHdI2, SmHFF, SmHdIFF, SmHF2F2
C
      double precision  tm(5, 4)
      character         tmLbl(4)*11
C
C pointer
      integer  ptMult, ptdIobs
      integer  ptfcg
      integer  pthMs, ptkMs, ptlMs, ptfts
      integer  ptwsp
C
C parameters
      double precision  Zero
      parameter(Zero = 0.0D0)
C
C begin
C
C List of reflections is in P1. To speed up the summations,
C reduce this list to asymmetric unit. Mult > 0 for asym. reflections,
C Mult = 0 for sym. equiv. reflections or unselected reflections.
      ptMult = AllHp(integ4(nRef))
C
      call ftsMult(nRef, xrh, xrk, xrl, heap(ptMult), TSel, QHerm,
     &             mSymMx, SymMx, nSymMx,
     &             Verbose, WrnLev)
      if (WrnLev .ge. 15)
     &  call ftsCheckMult(nRef, xrh, xrk, xrl, TSel, heap(ptMult),
     &                    QHerm, mSymMx, SymMx, nSymMx)
C
C.......................................................................
C
      ptdIobs = AllHp(ireal8(nRef))
C
      call ftsSmH(heap(ptMult), Fobs, Fpart, heap(ptdIobs),
     &            nRef, QFpart,
     &            SmH, SmHdI2, SmHFF, SmHdIFF, SmHF2F2)
C
      if (WrnLev .ge. 15 .or. Verbose) then
        write(6, '(1X, A, I12)')   'TSMAP: SmH =     ', SmH
        write(6, '(1X, A, G12.6)') 'TSMAP: SmHdI2 =  ', SmHdI2
        write(6, '(1X, A, G12.6)') 'TSMAP: SmHFF =   ', SmHFF
        write(6, '(1X, A, G12.6)') 'TSMAP: SmHdIFF = ', SmHdIFF
        write(6, '(1X, A, G12.6)') 'TSMAP: SmHF2F2 = ', SmHF2F2
      end if
C
C.......................................................................
C
C To speed up the summations, copy fcalc to 3D matrix.
C Determine size for 3D matrix to store fcalc
      call DimH3D(nRef, xrh, xrk, xrl, TSel, QHerm, hfcg, nfcg)
C
      if (WrnLev .ge. 15 .or. Verbose) then
        write(6, '(1X,A,3(1X,I5))') 'TSMAP: hfcg = ',
     &    (hfcg(j, 2), j = 1, 3)
        write(6, '(1X,A,3(1X,I5))') 'TSMAP: nfcg = ',
     &    (nfcg(j),    j = 1, 3)
      end if
C Allocate 3D matrix
      Pnfcg = nfcg(1) * nfcg(2) * nfcg(3)
      ptfcg = AllHp(icplx8(Pnfcg))
C
C Copy fcalc to 3D matrix
      call fcgSet(xrh, xrk, xrl, Fcalc, TSel, nRef, QHerm,
     &            nfcg(1), nfcg(2), nfcg(3), heap(ptfcg))
C
C.......................................................................
C
C Compute integer factors [uv]FFTgrid by which the temporary workspaces
C have to be bigger than the real space object for final storage.
C
C Factor for ws1 & ws2 (Navaza paper, p. 447: "This involves reciprocal
C vectors up to twice the data resolution.")
      call SetfFFTgrid('uFFTgrid', uFFTgrid, hfcg(1, 2), Rn, 2,
     &                 FTBase, FTCRBx, FTPrim)
      if (WrnLev .ge. 5 .or. Verbose)
     &  write(6, '(1X,A,3(1X,I3))') 'TSMAP: uFFTgrid set to ',
     &    (uFFTgrid(j), j = 1, 3)
C
C Factor for ws3 (Navaza paper, p. 447: "These terms involve reciprocal
C vectors up to four times the data resolution.")
      call SetfFFTgrid('vFFTgrid', vFFTgrid, hfcg(1, 2), Rn, 4,
     &                 FTBase, FTCRBx, FTPrim)
      if (WrnLev .ge. 5 .or. Verbose)
     &  write(6, '(1X,A,3(1X,I3))') 'TSMAP: vFFTgrid set to ',
     &    (vFFTgrid(j), j = 1, 3)
C
C.......................................................................
C
C Determine physical dimensions for workspaces
      call ftsSetXm(Rn, uFFTgrid, Un, Um, PUm, FTAvoi,
     &              'scratch array U', Verbose, WrnLev)
      call ftsSetXm(Rn, vFFTgrid, Vn, Vm, PVm, FTAvoi,
     &              'scratch array V', Verbose, WrnLev)
C
C.......................................................................
C
C Allocate temp space to store sym. equiv. hkl in the summation
C procedure
      pthMs = AllHp(integ4(nSymMx))
      ptkMs = AllHp(integ4(nSymMx))
      ptlMs = AllHp(integ4(nSymMx))
      ptfts = AllHp(icplx8(nSymMx))
C
C.......................................................................
C
      if (.not. LessMemory .and. FPrec .eq. 'SING') then
C Single precision
C First (huge) transform: ws3
        if (Verbose) write(6, '(1X, A)')
     &    'TSMAP: Allocating scratch array V (single precision)'
        ptwsp = AllHp(ireal4(PVm))
        if (Timer .ge. 0) call VCPU(tm(1, 1))
        call sftsSw3(Vm(1), Vm(2), Vm(3), Vn, heap(ptwsp),
     &               heap(ptfcg), nfcg,
     &               xrh, xrk, xrl,
     &               heap(ptMult), Fpart, nRef, QHerm, QFpart,
     &               SymMx, mSymMx, nSymMx, STBF,
     &               heap(pthMs), heap(ptkMs), heap(ptlMs), heap(ptfts))
        if (Timer .ge. 0) call VCPU(tm(2, 1))
        call ShrFFT(Rn(1), Rn(2), Rn(3),
     &              Vm(1), Vm(2), Vm(3),
     &              Vn(1), Vn(2), Vn(3),
     &              vFFTgrid, .false., 1,
     &              0,
     &              0,
     &              heap(ptwsp),
     &              heap(ptwsp),
     &              heap(ptwsp),
     &              rRho)
        if (Timer .ge. 0) call VCPU(tm(3, 1))
        if (Timer .ge. 0) call VCPU(tm(4, 1))
        call FreHp(ptwsp, ireal4(PVm))
C
C Second (smaller) transform: ws1
        if (Verbose) write(6, '(1X, A)')
     &    'TSMAP: Allocating scratch array U (single precision)'
        ptwsp = AllHp(ireal4(PUm))
        if (Timer .ge. 0) call VCPU(tm(1, 2))
        call sftsSw1(Um(1), Um(2), Um(3), Un, heap(ptwsp),
     &               heap(ptfcg), nfcg,
     &               xrh, xrk, xrl,
     &               heap(ptMult), Fpart, nRef, QHerm, QFpart,
     &               SymMx, mSymMx, nSymMx, STBF,
     &               heap(pthMs), heap(ptkMs), heap(ptlMs), heap(ptfts))
        if (Timer .ge. 0) call VCPU(tm(2, 2))
        call ShrFFT(Rn(1), Rn(2), Rn(3),
     &              Um(1), Um(2), Um(3),
     &              Un(1), Un(2), Un(3),
     &              uFFTgrid, .false., 1,
     &              0,
     &              0,
     &              heap(ptwsp),
     &              heap(ptwsp),
     &              heap(ptwsp),
     &              heap(ptwsp))
        if (Timer .ge. 0) call VCPU(tm(3, 2))
        call ftsCw1(Rn(1), Rn(2), heap(ptwsp),
     &              SmH, SmHdI2, SmHFF, SmHF2F2,
     &              Rn(1), Rn(2), Rn(3), rRho)
        if (Timer .ge. 0) call VCPU(tm(4, 2))
C
C Third transform: ws2
        if (Timer .ge. 0) call VCPU(tm(1, 3))
        call sftsSw2(Um(1), Um(2), Um(3), Un, heap(ptwsp),
     &               heap(ptfcg), nfcg,
     &               xrh, xrk, xrl,
     &               heap(ptMult), heap(ptdIobs), Fpart,
     &               nRef, QHerm, QFpart,
     &               SymMx, mSymMx, nSymMx, STBF,
     &               heap(pthMs), heap(ptkMs), heap(ptlMs), heap(ptfts))
        if (Timer .ge. 0) call VCPU(tm(2, 3))
        call ShrFFT(Rn(1), Rn(2), Rn(3),
     &              Um(1), Um(2), Um(3),
     &              Un(1), Un(2), Un(3),
     &              uFFTgrid, .false., 1,
     &              0,
     &              0,
     &              heap(ptwsp),
     &              heap(ptwsp),
     &              heap(ptwsp),
     &              heap(ptwsp))
        if (Timer .ge. 0) call VCPU(tm(3, 3))
        call ftsCw2(Rn(1), Rn(2), heap(ptwsp),
     &              SmHdIFF,
     &              Rn(1), Rn(2), Rn(3), rRho)
        if (Timer .ge. 0) call VCPU(tm(4, 3))
        call FreHp(ptwsp, ireal4(PUm))
C
      else if (.not. LessMemory) then
C Double precision
        fred = ireal8(256) / ireal4(256)
        if (fred * ireal4(256) .ne. ireal8(256)) then
          call WrnDie(-5, 'FastTS', 'Fatal Coding Error')
          call Die
        end if
        if (Verbose) write(6, '(1X, A)')
     &    'TSMAP: Allocating scratch array V (double precision)'
        ptwsp = AllHp(ireal8(PVm))
        if (Timer .ge. 0) call VCPU(tm(1, 1))
        call dftsSw3(Vm(1), Vm(2), Vm(3), Vn, heap(ptwsp),
     &               heap(ptfcg), nfcg,
     &               xrh, xrk, xrl,
     &               heap(ptMult), Fpart, nRef, QHerm, QFpart,
     &               SymMx, mSymMx, nSymMx, STBF,
     &               heap(pthMs), heap(ptkMs), heap(ptlMs), heap(ptfts))
        if (Timer .ge. 0) call VCPU(tm(2, 1))
        call ShrFFT(Rn(1), Rn(2), Rn(3),
     &              Vm(1), Vm(2), Vm(3),
     &              Vn(1), Vn(2), Vn(3),
     &              vFFTgrid, .true., fred,
     &              heap(ptwsp),
     &              heap(ptwsp),
     &              0,
     &              heap(ptwsp),
     &              heap(ptwsp),
     &              rRho)
        if (Timer .ge. 0) call VCPU(tm(3, 1))
        if (Timer .ge. 0) call VCPU(tm(4, 1))
        call FreHp(ptwsp, ireal8(PVm))
C
C Second (smaller) transform: ws1
        if (Verbose) write(6, '(1X, A)')
     &    'TSMAP: Allocating scratch array U (double precision)'
        ptwsp = AllHp(ireal8(PUm))
        if (Timer .ge. 0) call VCPU(tm(1, 2))
        call dftsSw1(Um(1), Um(2), Um(3), Un, heap(ptwsp),
     &               heap(ptfcg), nfcg,
     &               xrh, xrk, xrl,
     &               heap(ptMult), Fpart, nRef, QHerm, QFpart,
     &               SymMx, mSymMx, nSymMx, STBF,
     &               heap(pthMs), heap(ptkMs), heap(ptlMs), heap(ptfts))
        if (Timer .ge. 0) call VCPU(tm(2, 2))
        call ShrFFT(Rn(1), Rn(2), Rn(3),
     &              Um(1), Um(2), Um(3),
     &              Un(1), Un(2), Un(3),
     &              uFFTgrid, .true., fred,
     &              heap(ptwsp),
     &              heap(ptwsp),
     &              0,
     &              heap(ptwsp),
     &              heap(ptwsp),
     &              heap(ptwsp))
        if (Timer .ge. 0) call VCPU(tm(3, 2))
        call ftsCw1(Rn(1), Rn(2), heap(ptwsp),
     &              SmH, SmHdI2, SmHFF, SmHF2F2,
     &              Rn(1), Rn(2), Rn(3), rRho)
        if (Timer .ge. 0) call VCPU(tm(4, 2))
C
C Third transform: ws2
        if (Timer .ge. 0) call VCPU(tm(1, 3))
        call dftsSw2(Um(1), Um(2), Um(3), Un, heap(ptwsp),
     &               heap(ptfcg), nfcg,
     &               xrh, xrk, xrl,
     &               heap(ptMult), heap(ptdIobs), Fpart,
     &               nRef, QHerm, QFpart,
     &               SymMx, mSymMx, nSymMx, STBF,
     &               heap(pthMs), heap(ptkMs), heap(ptlMs), heap(ptfts))
        if (Timer .ge. 0) call VCPU(tm(2, 3))
        call ShrFFT(Rn(1), Rn(2), Rn(3),
     &              Um(1), Um(2), Um(3),
     &              Un(1), Un(2), Un(3),
     &              uFFTgrid, .true., fred,
     &              heap(ptwsp),
     &              heap(ptwsp),
     &              0,
     &              heap(ptwsp),
     &              heap(ptwsp),
     &              heap(ptwsp))
        if (Timer .ge. 0) call VCPU(tm(3, 3))
        call ftsCw2(Rn(1), Rn(2), heap(ptwsp),
     &              SmHdIFF,
     &              Rn(1), Rn(2), Rn(3), rRho)
        if (Timer .ge. 0) call VCPU(tm(4, 3))
        call FreHp(ptwsp, ireal8(PUm))
C
      else
C Compact double precision
        call tslm(Verbose, Timer, tm(1, 1),
     &            3, Vm, Vn, vFFTgrid, FTAvoi,
     &            heap(ptfcg), nfcg,
     &            nRef, xrh, xrk, xrl,
     &            heap(ptMult), heap(ptdIobs), Fpart, QHerm, QFpart,
     &            SymMx, mSymMx, nSymMx, STBF,
     &            heap(pthMs), heap(ptkMs), heap(ptlMs), heap(ptfts),
     &            Rn, rRho,
     &            PrevVn, MStatV, ptRInV)
        if (Timer .ge. 0) call VCPU(tm(4, 1))
C
        ptwsp = AllHp(ireal4(Rn(1) * Rn(2) * Rn(3)))
C
        call tslm(Verbose, Timer, tm(1, 2),
     &            1, Um, Un, uFFTgrid, FTAvoi,
     &            heap(ptfcg), nfcg,
     &            nRef, xrh, xrk, xrl,
     &            heap(ptMult), heap(ptdIobs), Fpart, QHerm, QFpart,
     &            SymMx, mSymMx, nSymMx, STBF,
     &            heap(pthMs), heap(ptkMs), heap(ptlMs), heap(ptfts),
     &            Rn, heap(ptwsp),
     &            PrevUn, MStatU, ptRInU)
        call ftsCw1(Rn(1), Rn(2), heap(ptwsp),
     &              SmH, SmHdI2, SmHFF, SmHF2F2,
     &              Rn(1), Rn(2), Rn(3), rRho)
        if (Timer .ge. 0) call VCPU(tm(4, 2))
C
        call tslm(Verbose, Timer, tm(1, 3),
     &            2, Um, Un, uFFTgrid, FTAvoi,
     &            heap(ptfcg), nfcg,
     &            nRef, xrh, xrk, xrl,
     &            heap(ptMult), heap(ptdIobs), Fpart, QHerm, QFpart,
     &            SymMx, mSymMx, nSymMx, STBF,
     &            heap(pthMs), heap(ptkMs), heap(ptlMs), heap(ptfts),
     &            Rn, heap(ptwsp),
     &            PrevUn, MStatU, ptRInU)
        call ftsCw2(Rn(1), Rn(2), heap(ptwsp),
     &              SmHdIFF,
     &              Rn(1), Rn(2), Rn(3), rRho)
        if (Timer .ge. 0) call VCPU(tm(4, 3))
C
        call FreHp(ptwsp, ireal4(Rn(1) * Rn(2) * Rn(3)))
      end if
C
C.......................................................................
C
      call FreHp(ptfts, icplx8(nSymMx))
      call FreHp(ptlMs, integ4(nSymMx))
      call FreHp(ptkMs, integ4(nSymMx))
      call FreHp(pthMs, integ4(nSymMx))
C
      call FreHp(ptfcg, icplx8(Pnfcg))
      call FreHp(ptdIobs, ireal8(nRef))
      call FreHp(ptMult, integ4(nRef))
C
C.......................................................................
C
C Print CPU-time matrix
      if (Timer .gt. 0) then
        do k = 2, 5
          tm(k, 4) = Zero
        end do
C
        do j = 1, 3
          tm(5, j) = Zero
          do k = 4, 2, -1
            tm(k, j) = tm(k, j) - tm(k - 1, j)
            tm(k, 4) = tm(k, 4) + tm(k, j)
            tm(5, j) = tm(5, j) + tm(k, j)
          end do
        end do
C
        do k = 2, 4
          tm(5, 4) = tm(5, 4) + tm(k, 4)
        end do
C
        write(6, '(1X, 2A)')
     &    'TSMAP: CPU-time:                  ',
     &     'V          U1         U2    |    V+U1+U2'
        write(6, '(1X, 2A)')
     &    'TSMAP: CPU-time: -----------------',
     &     '----------------------------|-----------'
C
        tmLbl(1) = 'Summation'
        tmLbl(2) = 'FFT'
        tmLbl(3) = 'Combination'
        tmLbl(4) = 'Sum'
C
        do k = 2, 5
          write(6, '(1X, 2A, 3(1X, F10.4), A, F10.4)')
     &      'TSMAP: CPU-time: ', tmLbl(k - 1),
     &      (tm(k, j), j = 1, 3), ' | ', tm(k, 4)
          if (k .eq. 4)
     &      write(6, '(1X, 2A)')
     &        'TSMAP: CPU-time: -----------------',
     &         '----------------------------|-----------'
        end do
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine Comp2Real(n, c, r)
      IMPLICIT NONE
C I/O
      integer           n
      double complex    c(*)
      double precision  r(*)
C
C Copy abs(double complex) to double precision array
C
C local
      integer  i
C
C begin
      do i = 1, n
        r(i) = sqrt(DBLE(c(i))**2 + dimag(c(i))**2)
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine CountSelRef(nRef, TSel, nSel)
      IMPLICIT NONE
C I/O
      integer  nRef, TSel(*), nSel
C
C Count selected reflections
C
C local
      integer  iRef
C
C begin
      nSel = 0
      do iRef = 1, nRef
        if (TSel(iRef) .ne. 0) nSel = nSel + 1
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine tsmapp(Reset,
     &                  iFobs, iP1Fcalc, iTrFcalc, iFpart, iTo,
     &                  Method, QSpecPos,
     &                  xRhoNum, xRhoNam,
     &                  xSFNum, xSFNam, xSFType,
     &                  FPrec, LessMemory, Verbose, BailOut)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      logical        Reset
      integer        iFobs, iP1Fcalc, iTrFcalc, iFpart, iTo
      character*(*)  Method
      logical        QSpecPos
      integer        xRhoNum
      character*(*)  xRhoNam(*)
      integer        xSFNum
      character*(*)  xSFNam(*), xSFType(*)
      character*(*)  FPrec
      logical        LessMemory, Verbose, BailOut
C
C This is the parser for tsmap
C
C begin
      Reset      = .false.
      iFobs      = 0
      iP1Fcalc   = 0
      iTrFcalc   = 0
      iFpart     = 0
      iTo        = 0
      Method     = 'FFT'
      QSpecPos   = .false.
      FPrec      = 'DOUB'
      LessMemory = .false.
      Verbose    = .true.
      BailOut    = .true.
C
      CALL PUSEND('TSMAP>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('TSMAP>')
      CALL MISCOM('TSMAP>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-search-tsmap')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'RESE') THEN
        CALL NEXTLO('RESEt=', Reset)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FOBS') THEN
        call GetiObj(iFobs, xSFNum, xSFNam, ' ',
     &               'FOBSfrom=', 'reciprocal', 'TSMAP')
        if (iFobs .gt. 0) then
          if (      xSFType(iFobs) .ne. 'REAL'
     &        .and. xSFType(iFobs) .ne. 'COMP')
     &      call ChkSFType(iFobs, xSFNam, xSFType, 'COMP', 'TSMAP')
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'P1FC') THEN
        call GetiObj(iP1Fcalc, xSFNum, xSFNam, ' ',
     &               'P1FCalcfrom=', 'reciprocal', 'TSMAP')
        call ChkSFType(iP1Fcalc, xSFNam, xSFType, 'COMP', 'TSMAP')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TRFC') THEN
        call GetiObj(iTrFcalc, xSFNum, xSFNam, ' ',
     &               'TRFCalc=', 'reciprocal', 'TSMAP')
        call ChkSFType(iTRFcalc, xSFNam, xSFType, 'COMP', 'TSMAP')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FPAR') THEN
        call GetiObj(iFpart, xSFNum, xSFNam, ' ',
     &               'FPARtfrom=', 'reciprocal', 'TSMAP')
        call ChkSFType(iFpart, xSFNam, xSFType, 'COMP', 'TSMAP')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TO  ') THEN
        call GetiObj(iTO, xRhoNum, xRhoNam, ' ',
     &               'TO=', 'real', 'TSMAP')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'METH') THEN
        CALL NEXTST('METHod=', WD)
        if (      (WD(1:4) .ne. 'DIRE' .or. WDLEN .lt. 4)
     &      .and. (WD(1:4) .ne. 'FFT ' .or. WDLEN .lt. 3)) then
          CALL DSPERR('TSMAP', 'unknown qualifier')
        else
          METHod = WD(1:4)
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FPRE') THEN
        CALL NEXTST('FPREcision=', WD)
        if      (WD(1:4) .eq. 'SING') then
          FPrec = 'SING'
        else if (WD(1:4) .eq. 'DOUB') then
          FPrec = 'DOUB'
        else
          CALL DSPERR('TSMAP', 'unknown qualifier')
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SPEC') THEN
        CALL NEXTLO('SPECialpositions=',QSpecPos)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'BAIL') THEN
        CALL NEXTLO('BAILout=', BailOut)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'LESS') THEN
        CALL NEXTLO('LESSmemory=', LessMemory)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'VERB') THEN
        CALL NEXTLO('VERBose=', Verbose)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
        WRITE(6,'(2A)') ' -----------tsmap-parameters',
     &   '----------------------------------------------------'
        if (Reset) then
          write(6, '(1X, A)') 'RESEt=TRUE'
        else
          write(6, '(1X, A)') 'RESEt=FALSE'
        end if
        call EchoObj(iFobs,  xSFNam, 'FOBSfrom=')
        call EchoObj(iP1Fcalc, xSFNam, 'P1FCalcfrom=')
        call EchoObj(iTrFcalc, xSFNam, 'TRFCalc=')
        call EchoObj(iFpart, xSFNam, 'FPARtfrom=')
        call EchoObj(iTo, xRhoNam, 'TO=')
        if (Method .eq. 'DIRE') write(6, '(1X, A)') 'METHod=DIREct'
        if (Method .eq. 'FFT ') write(6, '(1X, A)') 'METHod=FFT'
        if (FPrec .eq. 'SING') write(6, '(1X, A)') 'FPREcision=SINGle'
        if (FPrec .eq. 'DOUB') write(6, '(1X, A)') 'FPREcision=DOUBle'
        if (LessMemory) then
          write(6, '(1X, A)') 'LESSmemory=TRUE'
        else
          write(6, '(1X, A)') 'LESSmemory=FALSE'
        end if
        if (QSpecPos) then
          write(6, '(1X, A)') 'SPECialpositions=TRUE'
        else
          write(6, '(1X, A)') 'SPECialpositions=FALSE'
        end if
        if (Verbose) then
          write(6, '(1X, A)') 'VERBose=TRUE'
        else
          write(6, '(1X, A)') 'VERBose=FALSE'
        end if
        WRITE(6,'(2A)') ' ---------------------------',
     &   '----------------------------------------------------'
C=====================================================================
      ELSE
      CALL CHKEND('TSMAP>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      return
      end
C}
C======================================================================
C{
      subroutine tsmap(ptxrh, ptxrk, ptxrl, nRef, ptTSel, QHerm,
     &                 na, nb, nc, nrRho, niRho,
     &                 xRhoNum, xRhoNam, hprrho, hpirho,
     &                 xSFNum, xSFNam, xSFType, hpSF,
     &                 xrNsym, STBF,
     &                 XRE,
     &                 MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &                 XRTR,XRSCAL,XCVTEST,
     &                 TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &                 TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &                 HPMULT,HPTYPE,XRMREF,XRCELL,XRVOL)
      IMPLICIT NONE
C
C Unified Conventional and F2F2 FFT Translation Function main procedure.
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'flagmap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'xfft.inc'
      INCLUDE 'timer.inc'
      integer        ptxrh, ptxrk, ptxrl
      integer        nRef
      integer        ptTSel
      logical        QHerm
      integer        na, nb, nc, nrRho, niRho
      integer        xRhoNum
      character*(*)  xRhoNam(*)
      integer        hprrho(*), hpirho(*)
      integer        xSFNum
      character*(*)  xSFNam(*), xSFType(*)
      integer        hpSF(*)
      integer        xrNsym, STBF
C
      DOUBLE PRECISION  XRE
      INTEGER MBINS
      DOUBLE PRECISION XBINLOW, XBINHIGH
      DOUBLE PRECISION BINSHELL(*)
      DOUBLE PRECISION XRTR(3,3)
      DOUBLE PRECISION XRSCAL
      LOGICAL XCVTEST
      INTEGER TRPNMX, TRPNX
      CHARACTER*(*) TRPN(4,TRPNMX)
      INTEGER TRPNL(4,TRPNMX), TRPNN
      DOUBLE COMPLEX TRPNDB(4,TRPNMX)
      INTEGER TRPNMLT(TRPNMX)
      CHARACTER*2 TRPNTYP(TRPNMX), TRPNDOM(TRPNMX)
      INTEGER TRPNLEV(TRPNMX)
      INTEGER TDEPTH
      INTEGER HPMULT
      INTEGER HPTYPE
      INTEGER XRMREF
      DOUBLE PRECISION XRCELL(3,3), XRVOL
C
C local
      logical    QFpart, QFast, QSpecPos
      logical    Reset, LessMemory, Verbose, BailOut
      character  Method*4, FPrec*4
      integer    iFobs, iP1Fcalc, iTrFcalc, iFpart, iTo
      integer    nSel, PRn
C
      double precision  tmEntr, tmExit
C
C pointer
      integer  ptFobs, ptP1Fcalc, ptTrFcalc, ptFpart
      integer  ptRealFobs
C
C parameters
      double precision  Zero
      parameter(Zero = 0.0D0)
      real      SZero
      parameter(SZero=0.0)
C
C begin
      call tsmapp(Reset,
     &            iFobs, iP1Fcalc, iTrFcalc, iFpart, iTo,
     &            Method, QSpecPos,
     &            xRhoNum, xRhoNam,
     &            xSFNum, xSFNam, xSFType,
     &            FPrec, LessMemory, Verbose, BailOut)
C
      if (Reset) call tslmRst
C
      if (      iFobs    .lt. 1
     &    .and. iP1Fcalc .lt. 1
     &    .and. iTrFcalc .lt. 1
     &    .and. iFpart   .lt. 1
     &    .and. iTo      .lt. 1) return
C
      if (Timer .gt. 0) call VCPU(tmEntr)
C
      if (Method .eq. 'DIRE') then
        QFast = .false.
      else
        QFast = .true.
      end if
C
      if (iFobs .lt. 1) then
        call WrnDie(0, 'TSMAP', 'FOBSfrom not defined.')
        return
      end if
          ptFobs = hpSF(iFobs)
      if (ptFobs .eq. 0) then
        call WrnDie(0, 'TSMAP', 'FOBSfrom object not set.')
        return
      end if
C
      if (iP1Fcalc .lt. 1) then
        call WrnDie(0, 'TSMAP', 'P1FCalcfrom not defined.')
        return
      end if
          ptP1Fcalc = hpSF(iP1Fcalc)
      if (ptP1Fcalc .eq. 0) then
        call WrnDie(0, 'TSMAP', 'P1FCalcfrom object not set.')
        return
      end if
C
      ptTrFcalc = 0
C
      if (iTrFcalc .lt. 1) then
        if (.not. QFast .or. QSpecPos) then
          call WrnDie(0, 'TSMAP', 'TRFCalc not defined.')
          return
        end if
      else
        if (hpSF(iTrFcalc) .eq. 0)
     &    call xSFAl(hpSF(iTrFcalc), XRMREF, xSFType(iTrFcalc))
        ptTrFcalc = hpSF(iTrFcalc)
      end if
C
      if (iFpart .lt. 1) then
        ptFpart = 0
        QFpart = .false.
      else
            ptFpart = hpSF(iFpart)
        if (ptFpart .eq. 0) then
          call WrnDie(0, 'TSMAP', 'FPARtfrom object not set.')
          return
        end if
        QFpart = .true.
      end if
C
      if (iTo .lt. 1) then
        call WrnDie(0, 'TSMAP', 'Real space object not specified.')
        return
      end if
C
      if (xrNsym .ne. 1) then
        call WrnDie(0, 'TSMAP',
     &    'This module can only be used in P1.')
        return
      end if
C
      if (nRef .lt. 1) then
        call WrnDie(0, 'TSMAP',
     &    'Number of reflections is zero.')
        return
      end if
C
      call CountSelRef(nRef, heap(ptTSel), nSel)
      if (nSel .lt. 1) then
        call WrnDie(0, 'TSMAP',
     &    'Number of selected reflections is zero.')
        return
      end if
C
      if (     .not. ValidFlagMap
     &    .or. RnFM(1) .ne. na
     &    .or. RnFM(2) .ne. nb
     &    .or. RnFM(3) .ne. nc) then
        call WrnDie(0, 'TSMAP',
     &    'FlagMap undefined or out of date: call FMAP first.')
        return
      end if
C
      PRn = RnFM(1) * RnFM(2) * RnFM(3)
C
      if (nrRho .ne. PRn) then
        call WrnDie(-5, 'TSMAP', 'Fatal Coding Error.')
        call Die
      end if
C
      if (hprrho(iTo) .eq. 0) then
        call xMapAl(hprrho(iTo), hpirho(iTo), QHerm, nrRho, niRho)
      else if (.not. QHerm) then
        call FillR4(heap(hpirho(iTo)), niRho, SZERO)
      end if
C
      if (QFast) then
        if      (xSFType(iFobs) .eq. 'COMP') then
          ptRealFobs = AllHp(ireal8(nRef))
          call Comp2Real(nRef, heap(ptFobs), heap(ptRealFobs))
        else if (xSFType(iFobs) .eq. 'REAL') then
          ptRealFobs = ptFobs
        else
          call WrnDie(-5, 'TSMAP', 'Fatal Coding Error.')
          call Die
        end if
        if (WrnLev .ge. 10)
     &    write(6, '(1X, A)') 'TSMAP: Running FFT Translation Search'
        call FastTS(heap(ptxrh), heap(ptxrk), heap(ptxrl),
     &              nRef, heap(ptTSel), QHerm,
     &              heap(ptRealFobs), heap(ptP1Fcalc),
     &              heap(ptFpart), QFpart,
     &              mSymMx, SymMx, nSymMx, STBF,
     &              RnFM, heap(hprrho(iTo)),
     &              FTBase, FTCRBx, FTPrim, FTAvoi, FPrec,
     &              LessMemory, Verbose)
        if (ptRealFobs .ne. ptFobs) then
          call FreHP(ptRealFobs, ireal8(nRef))
          ptRealFobs = 0
        end if
C Check if results are consistent with FlagMap
        call MapInd(heap(hprrho(iTo)), heap(ptFlagMap), RnFM,
     &              .true., .true., BailOut, WrnLev, 'TSMAP')
      end if
C
      if (.not. QFpart .and. QSpecPos) then
        QSpecPos = .false.
        if (WrnLev .ge. 10) then
          write(6, '(1X, 2A)')
     &      'TSMAP: FPARtfrom=<undefined> ==>',
     &      ' SPECialpositions=TRUE ignored.'
        end if
      end if
C
      if (.not. QFast .or. QSpecPos) then
        if (Verbose .or. WrnLev .ge. 10) then
          if (QFast .and. QSpecPos) then
            write(6, '(1X, 2A)')
     &        'TSMAP: Running Direct Translation Search',
     &        ' for special positions'
          else
            write(6, '(1X, A, I9, A)')
     &        'TSMAP: Running Direct Translation Search for ',
     &        nIndFM, ' independent grid point(s)'
          end if
        end if
        call ConvTS(nRef, ptxrh, ptxrk, ptxrl, heap(ptTSel), ptTSel,
     &              ptP1Fcalc, ptTrFcalc,
     &              QHerm,
     &              SymMx, ITSymMx, mSymMx, nSymMx, STBF,
     &              RnFM(1), RnFM(2), RnFM(3),
     &              heap(hprrho(iTo)), heap(ptFlagMap),
     &              QSpecPos, QFast,
     &              Verbose, nIndFM,
     &              LessMemory,
     &              xSFNum, xSFNam, xSFType, hpSF,
     &              XRE,
     &              MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &              XRTR,XRSCAL,XCVTEST,
     &              TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &              TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &              HPMULT,HPTYPE,XRMREF,XRCELL,XRVOL)
      end if
C
      if (WrnLev .ge. 5)
     &  call MapStat(PRn, heap(hprrho(iTo)), 'TSMAP: t', .true.)
C
      if (Timer .gt. 0) then
        call VCPU(tmExit)
        write(6, '(1X, A, F10.4, A)')
     &    'TSMAP: CPU-time: ', tmExit - tmEntr, ' seconds total'
      end if
C
      return
      end
C}
