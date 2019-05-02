C=======================================================================
C{
      subroutine FindInRow(nc, Row, iz, ic, found)
      IMPLICIT NONE
C I/O
      integer  nc, Row(*), iz, ic
      logical  found
C
C Find the index ic of the element with the value iz in the integer
C array Row. Row has nc sorted elements.
C The algorithm below works in two stages. For search intervals
C with greater than 20 elements, the search interval is bisected
C (cut in half). If the interval has less than or equal to 20 elements,
C the value iz is searched for by a simple loop over the remaining interval.
C One pass of the simple loop is faster than one pass of the bisecting
C algorithm, therefore small intervals are faster searched that way. The
C break-even point in interval size for the two search strategies is
C machine- and compiler-specific. The value 20 is only an estimate.
C
C local
      integer  is, ie, id
C
C begin
      found = .false.
C
      ic = 1
      is = 1
      ie = nc
 10   continue
        id = ie - is
        if (id .lt. 20) then
          do ic = is, ie
            if (Row(ic) .ge. iz) then
              if (Row(ic) .eq. iz) found = .true.
              return
            end if
          end do
          return
        end if
        ic = is + id / 2
        if (Row(ic) .eq. iz) then
          found = .true.
          return
        end if
        if (Row(ic) .lt. iz) then
          ic = ic + 1
          is = ic
        else
          ie = ic - 1
        end if
      goto 10
C
      end
C}
C=======================================================================
C{
      subroutine Add2Row(RowInfoEl, iz, cf)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      integer           RowInfoEl(4)
      integer           iz
CCC modification ATB 4/27/08
      double complex  cf
C
C This procedure is called from the summation loops of the
C "fast translation search with less memory."
C The coefficient cf is added to the row element with the value iz.
C RowInfoEl is one element of the 2D array RowInfo decribed in the
C subroutine tslmIni (this means a structure of 2 integer and 2 pointers).
C
C local
      integer  mc, nc, ic
      logical  found
C
C begin
          nc = RowInfoEl(1)
          mc = RowInfoEl(2)
      if (mc .eq. 0) then
        RowInfoEl(1) = 1
        RowInfoEl(2) = 1
        RowInfoEl(3) = AllHp(integ4(1))
        RowInfoEl(4) = AllHp(icplx8(1))
        call SetI4Elem(heap(RowInfoEl(3)), 1, iz)
        call SetC8Elem(heap(RowInfoEl(4)), 1, cf)
      else
        ic = 1
        found = .false.
        if (nc .ne. 0)
     &    call FindInRow(nc, heap(RowInfoEl(3)), iz, ic, found)
        if (.not. found) then
          if (mc .eq. nc) then
            RowInfoEl(2) = mc + 1
            RowInfoEl(3) = ReaHp(RowInfoEl(3), integ4(mc), integ4(mc+1))
            RowInfoEl(4) = ReaHp(RowInfoEl(4), icplx8(mc), icplx8(mc+1))
          end if
          RowInfoEl(1) = nc + 1
          if (ic .le. nc) then
            call ShUpI4(heap(RowInfoEl(3)), ic, nc + 1)
            call ShUpC8(heap(RowInfoEl(4)), ic, nc + 1)
          end if
          call SetI4Elem(heap(RowInfoEl(3)), ic, iz)
          call SetC8Elem(heap(RowInfoEl(4)), ic, cf)
        else
          call AddC8Elem(heap(RowInfoEl(4)), ic, cf)
        end if
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine RowInfoStat(Wn1, Wn2, Wn3, RowInfo)
      IMPLICIT NONE
C I/O
      integer  Wn1, Wn2, Wn3
      integer  RowInfo(4, 0:Wn2-1, 0:(Wn1/2)-1)
C
C Compute and print some characteristics of the final state (after finishing
C the summation loop) of the 2D array described in the subroutine tslmIni.
C
C local
      integer  ix, iy, nc
      integer  nNon0Rows, Sum_nCoeff, Min_nCoeff, Max_nCoeff
C
C begin
      nNon0Rows = 0
      Sum_nCoeff = 0
      Min_nCoeff = Wn3
      Max_nCoeff = 0
C
      do ix = 0, (Wn1/2)-1
      do iy = 0,  Wn2-1
            nc = RowInfo(1, iy, ix)
        if (nc .gt. 0) then
          nNon0Rows = nNon0Rows + 1
          Sum_nCoeff = Sum_nCoeff + nc
          Min_nCoeff = min(Min_nCoeff, nc)
          Max_nCoeff = max(Max_nCoeff, nc)
        end if
      end do
      end do
C
      write(6, '(1X, A, I11, A, F8.2, A)')
     &  'TSMAP: Number of non-empty rows = ', nNon0Rows,
     &  ' = ', dfloat(nNon0Rows) / dfloat((Wn1/2) * Wn2)
     &                           * dfloat(100),
     &  ' % of 2D array'
      write(6, '(1X, A, I11, A, F8.2, A)')
     &  'TSMAP: Total number of coefficients = ', Sum_nCoeff,
     &  ' = ', dfloat(Sum_nCoeff) / dfloat((Wn1/2) * Wn2 * Wn3)
     &                            * dfloat(100),
     &  ' % of 3D array'
      write(6, '(1X, A, I11)')
     &  'TSMAP: Minimum number of coefficients per non-empty row = ',
     &  Min_nCoeff
      write(6, '(1X, A, I11)')
     &  'TSMAP: Maximum number of coefficients per non-empty row = ',
     &  Max_nCoeff
      write(6, '(1X, A, G14.6)')
     &  'TSMAP:    Mean number of coefficients per non-empty row = ',
     &  dfloat(Sum_nCoeff) / dfloat(max(1, nNon0Rows))
      write(6, '(1X, A, I11)')
C
      return
      end
C}
C=======================================================================
C{
      subroutine FreIzCf(Wn1, Wn2, RowInfo)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      integer  Wn1, Wn2
      integer  RowInfo(4, 0:Wn2-1, 0:(Wn1/2)-1)
C
C Free memory for indices and coefficients. See subroutine tslmIni.
C
C local
      integer  ix, iy, mc
C
C begin
      do ix = 0, (Wn1/2)-1
      do iy = 0,  Wn2-1
            mc = RowInfo(2, iy, ix)
        if (mc .gt. 0) then
          call FreHp(RowInfo(4, iy, ix), icplx8(mc))
          RowInfo(4, iy, ix) = 0
        end if
      end do
      end do
C
      do ix = 0, (Wn1/2)-1
      do iy = 0,  Wn2-1
            mc = RowInfo(2, iy, ix)
        if (mc .gt. 0) then
          call FreHp(RowInfo(3, iy, ix), integ4(mc))
          RowInfo(2, iy, ix) = 0
          RowInfo(3, iy, ix) = 0
        end if
      end do
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine AllIzCf(Wn1, Wn2, RowInfo)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      integer  Wn1, Wn2
      integer  RowInfo(4, 0:Wn2-1, 0:(Wn1/2)-1)
C
C Allocate memory for indices and coefficients if the row sizes are
C known from previous uses. See subroutine tslmIni.
C
C local
      integer  ix, iy, mc, nc
C
C begin
      do ix = 0, (Wn1/2)-1
      do iy = 0,  Wn2-1
            nc = RowInfo(1, iy, ix)
        if (nc .gt. 0) then
          RowInfo(3, iy, ix) = AllHp(integ4(nc))
          RowInfo(2, iy, ix) = nc
          RowInfo(1, iy, ix) = 0
        end if
      end do
      end do
C
      do ix = 0, (Wn1/2)-1
      do iy = 0,  Wn2-1
            mc = RowInfo(2, iy, ix)
        if (mc .gt. 0) then
          RowInfo(4, iy, ix) = AllHp(icplx8(mc))
        end if
      end do
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine tslmTr(Rn1, Rn2, Rn3,
     &                  Wn1, Wn2, Wn3,
     &                  fFFTgrid,
     &                  RowInfo,
     &                  Cm1,                Co,
     &                  Sm1, Sm2,      PSm, Sl,
     &                  Bm1, Bm2, Bm3, PBm, Bc, Br,
     &                  rRho)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      integer         Rn1, Rn2, Rn3
      integer         Wn1, Wn2, Wn3
      integer         fFFTgrid(3)
      integer         RowInfo(4, 0:Wn2-1, 0:(Wn1/2)-1)
      integer         Cm1
      complex         Co(0:Cm1-1)
      integer         Sm1, Sm2,      PSm
      complex         Sl(0:Sm1-1, 0:Sm2-1)
      integer         Bm1, Bm2, Bm3, PBm
      complex         Bc(0:  Bm1-1, 0:Bm2-1, 0:Bm3-1)
      real            Br(0:2*Bm1-1, 0:Bm2-1, 0:Bm3-1)
      real            rRho(0:Rn1-1, 0:Rn2-1, 0:Rn3-1)
C
C Compute 3D Fast Fourier Transform with coefficients tabulated in
C the 2D array described in the subroutine tslmIni.
C Store results in rRho. Co (Column), Sl (Slice) and Bc and Br (Box)
C are 1D, 2D and 3D scratch arrays, respectively.
C Bc and Br must point to the same memory area.
C fFFTgrid holds the sampling factors: Wn(i) = Rn(i) * fFFTgrid(i)
C
C local
      integer         ix, iy, iz, nc, ic, i
      double complex  cf
      logical         ZeroSlice, FFTERR
      integer         mFFTw(3)
      integer         ptFFTw(3)
C
C parameter
      complex   C4Zero
      real      R4Zero
      parameter(R4Zero = 0.0)
C
C begin
      C4Zero = cmplx(R4Zero, R4Zero)
C
      do i = 1, 3
         mFFTw(i) = 0
        ptFFTw(i) = 0
      end do
C
      call FillC4(Bc, PBm, C4Zero)
C
      do ix = 0, Wn1 / 2 - 1
        call FillC4(Sl, PSm, C4Zero)
        ZeroSlice = .true.
C
        do iy = 0, Wn2 - 1
          call FillC4(Co, Cm1, C4Zero)
C
              nc = RowInfo(1, iy, ix)
          if (nc .gt. 0) then
            ZeroSlice = .false.
            do ic = 1, nc
              call GetI4Elem(heap(RowInfo(3, iy, ix)), ic, iz)
              call GetC8Elem(heap(RowInfo(4, iy, ix)), ic, cf)
              Co(iz) = cf
            end do
C
            call sfft1c(Wn3, Co, FFTERR, mFFTw(3), ptFFTw(3))
            if (FFTERR) then
              call WrnDie(-5, 'tslmTr', 'Fatal: sfft1c failed.')
              call Die
            end if
            do iz = 0, Rn3 - 1
              Sl(iy, iz) = Co(iz * fFFTgrid(3))
            end do
            ZeroSlice = .false.
          end if
        end do
C
        if (.not. ZeroSlice) then
          do iz = 0, Rn3 - 1
            call sfft1c(Wn2, Sl(0, iz), FFTERR, mFFTw(2), ptFFTw(2))
            if (FFTERR) then
              call WrnDie(-5, 'tslmTr', 'Fatal: sfft1c failed.')
              call Die
            end if
            do iy = 0, Rn2 - 1
              Bc(ix, iy, iz) = Sl(iy * fFFTgrid(2), iz)
            end do
          end do
        end if
      end do
C
      do iz = 0, Rn3 - 1
        do iy = 0, Rn2 - 1
          call sfft1cr(Wn1, Bc(0, iy, iz), FFTERR,
     &                 mFFTw(1), ptFFTw(1))
          if (FFTERR) then
            call WrnDie(-5, 'tslmTr', 'Fatal: sfft1cr failed.')
            call Die
          end if
          do ix = 0, Rn1 - 1
            rRho(ix, iy, iz) = Br(ix * fFFTgrid(1), iy, iz)
          end do
        end do
      end do
C
      call sfft1cr(0, 0, FFTERR, mFFTw(1), ptFFTw(1))
      if (FFTERR) then
        call WrnDie(-5, 'tslmTr', 'Fatal: sfft1cr failed.')
        call Die
      end if
      call sfft1c (0, 0, FFTERR, mFFTw(2), ptFFTw(2))
      if (FFTERR) then
        call WrnDie(-5, 'tslmTr', 'Fatal: sfft1c failed.')
        call Die
      end if
      call sfft1c (0, 0, FFTERR, mFFTw(3), ptFFTw(3))
      if (FFTERR) then
        call WrnDie(-5, 'tslmTr', 'Fatal: sfft1c failed.')
        call Die
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine tslmIni
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'tsmap.inc'
C
C Initialize "fast translation search with less memory".
C
C The coefficients for the fast translation search are not stored
C in a huge 3D array, but tabulated in a 2D array "RowInfo".
C In this file, the variable Wn (3 integer) is always used for the
C dimensions of the (non-existing) huge 3D array. The 2D array RowInfo
C has dimensions (Wn(1), Wn(2)), and the indices ix and iy are implicit
C like in the 3D array. The index iz however is tabulated.
C This choice of implicit and tabulated indices is arbitrary.
C
C More source code is required if one wants a dynamic choice, for
C example to minimize the memory requirements [implicit indices are
C determined by selecting the two dimensions with mimimum product
C Wn(i)*Wn(j)] or optimize for speed [tabulated index is determined
C by selecting the dimension with smallest W(i)].
C
C Each element "RowInfoEl" of RowInfo is a structure of two integer and
C two pointers:
C   RowInfoEl(1) = number of tabulated index+coefficient pairs
C   RowInfoEl(2) = space alloacted for indices+coefficients
C   RowInfoEl(3) = pointer to array for indices (iz)
C   RowInfoEl(4) = pointer to array for coefficients (cf)
C
C The RowInfo Array is kept between invocations of the fast translation
C search for each of the two types of coefficient arrays ("U" and "V")
C until reset explicitly. The variable "MStat" indicates:
C   MStat = 0: RowInfo not allocated
C           1: RowInfo allocated
C
C The dimensions of the (non-existing) 3D coefficient arrays are
C stored in PrevUn and PrevVn. If these dimensions are unchanged between
C two invocations of the fast translation search, the RowInfo arrays
C are re-used. This greatly reduces the number of memory allocation
C calls and results in approximately 30% less run-time. Re-using the
C RowInfo arrays also greatly reduces memory fragmentation.
C
C This subroutine simply intitializes the book-keeping variables
C and the pointers to the RowInfo arrays ptRU and ptRV.
C
C begin
      PrevUn(1) = -1
      PrevUn(2) = -1
      PrevUn(3) = -1
      MStatU = 0
      ptRInU = 0
C
      PrevVn(1) = -1
      PrevVn(2) = -1
      PrevVn(3) = -1
      MStatV = 0
      ptRInV = 0
C
      return
      end
C}
C=======================================================================
C{
      subroutine tslmRst
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'tsmap.inc'
C
C Release memory allocated for RowInfo.
C See subroutine tslmIni.
C
C begin
      if (MStatU .gt. 0)
     &  call FreHp(ptRInU, integ4(4 * (PrevUn(1) / 2) * PrevUn(2)))
      if (MStatV .gt. 0)
     &  call FreHp(ptRInV, integ4(4 * (PrevVn(1) / 2) * PrevVn(2)))
C
      call tslmIni
C
      return
      end
C}
C=======================================================================
C{
      subroutine tslm(Verbose, Timer, tm,
     &                WhichW, Wm, Wn, fFFTgrid, FTAvoi,
     &                fcg, nfcg,
     &                nRef, xrh, xrk, xrl,
     &                Mult, dIobs, fpart, QHerm, QFpart,
     &                SymMx, mSymMx, nSymMx, STBF,
     &                hMs, kMs, lMs, fts,
     &                Rn, rRho,
     &                PrevWn, CurrMStat, ptRowInfo)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      logical           Verbose
      integer           Timer
      double precision  tm(*)
      integer           WhichW, Wm(3), Wn(3)
      integer           fFFTgrid(3), FTAvoi
      double complex    fcg(*)
      integer           nfcg(3)
      integer           nRef, xrh(*), xrk(*), xrl(*), Mult(*)
      double precision  dIobs(*)
      double complex    fpart(*)
      logical           QHerm, QFpart
      integer           mSymMx, SymMx(mSymMx, 3, 4), nSymMx, STBF
      integer           hMs(*), kMs(*), lMs(*)
      double complex    fts(*)
      integer           Rn(3)
      real              rRho(*)
      integer           PrevWn(3), CurrMStat
      integer           ptRowInfo
C
C "fast translation search with less memory"
C See subroutine tslmIni.
C
C local
      integer  i
      logical  TrashPrev
      integer  Cn(1), Cm(1)
      integer  Sn(2), Sm(2), PSm
      integer         Bm(3), PBm
      integer  ptC, ptS, ptB
C
C begin
      if (CurrMstat .gt. 0) then
        TrashPrev = .false.
        do i = 1, 3
          if (PrevWn(i) .ne. Wn(i)) TrashPrev = .true.
        end do
        if (TrashPrev) then
          call FreHp(ptRowInfo, integ4(4 * (PrevWn(1) / 2) * PrevWn(2)))
          CurrMStat = 0
          ptRowInfo = 0
        end if
      end if
C
      do i = 1, 3
        PrevWn(i) = Wn(i)
      end do
C
      if (CurrMStat .eq. 0) then
        ptRowInfo = AllHp(integ4(4 * (Wn(1) / 2) * Wn(2)))
        call Fill4(heap(ptRowInfo), 4 * (Wn(1) / 2) * Wn(2), 0)
        CurrMStat = 1
      else
        call AllIzCf(Wn(1), Wn(2), heap(ptRowInfo))
      end if
C
      if (Timer .ge. 0) call VCPU(tm(1))
      if      (WhichW .eq. 3) then
        call tslmSw3(Wn, Wn(1), Wn(2), heap(ptRowInfo),
     &               fcg, nfcg,
     &               xrh, xrk, xrl,
     &               Mult, fpart, nRef, QHerm, QFpart,
     &               SymMx, mSymMx, nSymMx, STBF,
     &               hMs, kMs, lMs, fts)
      else if (WhichW .eq. 1) then
        call tslmSw1(Wn, Wn(1), Wn(2), heap(ptRowInfo),
     &               fcg, nfcg,
     &               xrh, xrk, xrl,
     &               Mult, fpart, nRef, QHerm, QFpart,
     &               SymMx, mSymMx, nSymMx, STBF,
     &               hMs, kMs, lMs, fts)
      else if (WhichW .eq. 2) then
        call tslmSw2(Wn, Wn(1), Wn(2), heap(ptRowInfo),
     &               fcg, nfcg,
     &               xrh, xrk, xrl,
     &               Mult, dIobs, fpart, nRef, QHerm, QFpart,
     &               SymMx, mSymMx, nSymMx, STBF,
     &               hMs, kMs, lMs, fts)
      else
        call WrnDie(-5, 'WhichW', 'Fatal Coding Error')
        call Die
      end if
      if (Timer .ge. 0) call VCPU(tm(2))
C
      Cn(1) = Wn(3)
      call mFFTcr(1, Cn, FTAvoi, Cm)
      ptC = AllHp(icplx4(Cm(1)))
C
      Sn(1) = Wn(2)
      Sn(2) = Rn(3)
      call mFFTcr(2, Sn, FTAvoi, Sm)
      PSm = Sm(1) * Sm(2)
      ptS = AllHp(icplx4(PSm))
C
      Bm(1) = Wm(1) / 2
      Bm(2) = Rn(2)
      Bm(3) = Rn(3)
      PBm = Bm(1) * Bm(2) * Bm(3)
      ptB = AllHp(icplx4(PBm))
C
      call tslmTr(Rn(1), Rn(2), Rn(3),
     &            Wn(1), Wn(2), Wn(3),
     &            fFFTgrid,
     &            heap(ptRowInfo),
     &            Cm(1),                    heap(ptC),
     &            Sm(1), Sm(2),        PSm, heap(ptS),
     &            Bm(1), Bm(2), Bm(3), PBm, heap(ptB), heap(ptB),
     &            rRho)
      if (Timer .ge. 0) call VCPU(tm(3))
C
      call FreHp(ptB, icplx4(PBm))
      call FreHp(ptS, icplx4(PSm))
      call FreHp(ptC, icplx4(Cm(1)))
C
      if (Verbose)
     &  call RowInfoStat(Wn(1), Wn(2), Wn(3), heap(ptRowInfo))
      call FreIzCf(Wn(1), Wn(2), heap(ptRowInfo))
C
      return
      end
C}
