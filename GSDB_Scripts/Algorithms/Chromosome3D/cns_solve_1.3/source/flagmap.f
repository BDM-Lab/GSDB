C{
      subroutine MapInd(DensMap, FlagMap, Rn,
     &                  ChkInDep, MaxInDep, BailOut,
     &                  WrnLev, ErrID)
      implicit none
C I/O
      real       DensMap(*)
      integer    FlagMap(*), Rn(3)
      logical    ChkInDep, MaxInDep
      integer    WrnLev
      character  ErrID*(*)
      logical    BailOut
C
C Check (optional) correlation of dependent and independent
C density values and copy independent to dependent grid points.
C
C local
      integer           PRn, nDep, iMap, jMap
      double precision  minx, maxx, miny, maxy
      double precision  x, y, Sx, Sx2, Sy, Sy2, Sxy
      double precision  LRb, LRm, LRc
C
C parameters
      real  WrnLRc, ErrLRc
      parameter(WrnLRc = 0.99, ErrLRc = 0.98)
C
C begin
      PRn = Rn(3) * Rn(2) * Rn(1)
C
      nDep = 0
C
      if (ChkInDep) then
C Compute the Linear Correlation Coefficient of the dependent vs.
C independent density values. Print a warning if it deviates too
C much from 1.
        minx = 0.
        maxx = 0.
        miny = 0.
        maxy = 0.
        Sx  = 0.
        Sx2 = 0.
        Sy  = 0.
        Sy2 = 0.
        Sxy = 0.
        do iMap = 1, PRn
           jMap = FlagMap(iMap)
          if (jMap .gt. 0) then
            x = DensMap(iMap)
            y = DensMap(jMap)
            if (nDep .eq. 0) then
              minx = x
              maxx = x
              miny = y
              maxy = y
            else
              minx = min(minx, x)
              maxx = max(maxx, x)
              miny = min(miny, y)
              maxy = max(maxy, y)
            end if
            Sx  = Sx  + x
            Sx2 = Sx2 + x * x
            Sy  = Sy  + y
            Sy2 = Sy2 + y * y
            Sxy = Sxy + x * y
            nDep = nDep + 1
          end if
        end do
C
        if (nDep .gt. 0) then
          call Sum2bmc(nDep, minx, maxx, miny, maxy,
     &                 Sx, Sx2, Sy, Sy2, Sxy,
     &                 LRb, LRm, LRc)
C Check if all grid points are exactly equal
          if (      LRm .eq. 0.
     &        .and. LRc .eq. 0.
     &        .and.      maxy - miny
     &              .lt. 1.D-6 * max(abs(miny), abs(maxy)))
     &      LRc = 1.
C
          if (LRc .le. WrnLRc .or. WrnLev .ge. 10)
     &      write(6, '(1X, 3A, F8.5)')
     &        ErrID, ': Correlation of dependent and',
     &        ' independent grid points = ', LRc
C
          if (LRc .le. ErrLRc) then
            write(6, '(1X, 4A)')
     &        '%', ErrID, '-ERR: Poor correlation of dependent and',
     &        ' independent grid points.'
            if (BailOut)
     &        call WrnDie(-5, ErrID, 'Check FMAP parameters.')
          end if
        end if
      end if
C
      if (MaxInDep .and. (.not. ChkInDep .or. nDep .gt. 0)) then
C Set density of independent grid point to maximum of itself
C and all dependent grid points.
        do iMap = 1, PRn
           jMap = FlagMap(iMap)
          if (jMap .gt. 0)
     &      DensMap(jMap) = max(DensMap(jMap), DensMap(iMap))
        end do
      end if
C
C Use FlagMap pointers to make dependent Densities exactly equal.
      do iMap = 1, PRn
         jMap = FlagMap(iMap)
        if (jMap .gt. 0) DensMap(iMap) = DensMap(jMap)
      end do
C
      return
      end
C}
C=======================================================================
C{
      integer function SetNextPt(Rn, GridMis, V, M, f, PivotPt, NextPt)
      implicit none
C I/O
      integer  Rn(3), GridMis(3), V(3), M, f, PivotPt(3), NextPt(3)
C
C local
      integer  iv
C
C external
      integer   f3dix, iModPositive
      external  f3dix, iModPositive
C
C begin
      SetNextPt = 0
C
      do iv = 1, 3
        NextPt(iv) = V(iv) * f * Rn(iv)
        if (mod(NextPt(iv), M) .eq. 0) then
          NextPt(iv)
     &      = iModPositive(PivotPt(iv) + NextPt(iv) / M, Rn(iv))
        else
          GridMis(iv) = 1
          SetNextPt = -1
        end if
      end do
C
      if (SetNextPt .eq. 0) SetNextPt = f3dix(NextPt, Rn)
C
      return
      end
C}
C=======================================================================
C{
      subroutine CoreMarkOSs(FlagMap, Rn, GridMis, PivotPt, pMap,
     &                       ssV, ssM, n_ssVM)
      implicit none
C I/O
      integer  FlagMap(*), Rn(3), GridMis(3), PivotPt(3, 4), pMap
      integer  ssV(3, 3), ssM(3), n_ssVM
C
C local
      integer  mMap, iMap
      integer  iVM, f(3)
C
C Mark origin shift related grid points for a given pivot grid point
C
C external
      integer   SetNextPt
      external  SetNextPt
C
C begin
      if (n_ssVM .eq. 0) return
C
      mMap = Rn(1) * Rn(2) * Rn(3)
C
      f(1) = 0
C
      iVM = 1
C
 10   continue
            iMap = SetNextPt(Rn, GridMis,
     &                       ssV(1, iVM), ssM(iVM), f(iVM),
     &                       PivotPt(1, iVM), PivotPt(1, iVM + 1))
        if (iMap .gt. 0) then
          do while (FlagMap(iMap) .ne. 0)
            iMap = FlagMap(iMap)
          end do
          if (iMap .ne. pMap) FlagMap(iMap) = pMap
        end if
C
        if (iVM + 1 .le. n_ssVM) then
          iVM = iVM + 1
          f(iVM) = 0
          goto 10
        end if
C
 20     continue
          if (f(iVM) + 1 .lt. ssM(iVM)) then
            f(iVM) = f(iVM) + 1
            goto 10
          end if
C
          if (iVM .eq. 1) return
C
          iVM = iVM - 1
        goto 20
C
      end
C}
C=======================================================================
C{
      subroutine MarkOriginShifts(FlagMap, Rn,
     &                            ssV, ssM, n_ssVM,
     &                            WrnLev)
      implicit none
C I/O
      integer  FlagMap(*), Rn(3)
      integer  ssV(3, 3), ssM(3), n_ssVM
      integer  WrnLev
C
C Loop all independent grid points and mark origin shift related grid points.
C
C local
      integer    GridMis(3)
      integer    PivotPt(3, 4)
      integer    iMap, ia, ib, ic, iVM, i
      integer    nIndep, AdjM(3)
      character  xyz*3
C
C external
      integer   iLCM
      external  iLCM
C
C begin
      xyz = 'xyz'
C
      if (WrnLev .ge. 5)
     &  write(6, '(1X, 2A)')
     &    'FMAP: ',
     &    'Removing grid points related by allowed origin shifts'
C
      if (n_ssVM .gt. 0) then
        do iVM = 1, n_ssVM
          if (ssM(iVM) .ne. 0) then
            AdjM(iVM) = ssM(iVM)
          else
            AdjM(iVM) = 1
            do i = 1, 3
              if (ssV(i, iVM) .ne. 0) AdjM(iVM) = iLCM(AdjM(iVM), Rn(i))
            end do
          end if
        end do
C
        if (WrnLev .ge. 5) then
          write(6, '(1X, A)')
     &      'FMAP: s.s.Vector  Grid-Adjusted-Modulus'
          do iVM = 1, n_ssVM
            write(6, '(1X, A, 3(1X,I2), 1X, I3)')
     &        'FMAP: ',
     &        (ssV(i, iVM), i = 1, 3), AdjM(iVM)
          end do
        end if
      end if
C
      do i = 1, 3
        GridMis(i) = 0
      end do
C
      nIndep = 0
C
      iMap = 1
      do ic = 0, Rn(3) - 1
         PivotPt(3, 1) = ic
      do ib = 0, Rn(2) - 1
         PivotPt(2, 1) = ib
      do ia = 0, Rn(1) - 1
         PivotPt(1, 1) = ia
        if (FlagMap(iMap) .eq. 0) then
          call CoreMarkOSs(FlagMap, Rn, GridMis, PivotPt, iMap,
     &                     ssV, AdjM, n_ssVM)
          nIndep = nIndep + 1
        end if
        iMap = iMap + 1
      end do
      end do
      end do
C
      if (WrnLev .ge. 5) then
        write(6, '(1X, 2A, I9)') 'FMAP: ',
     &    '      Total number of grid points = ', Rn(1) * Rn(2) * Rn(3)
        write(6, '(1X, 2A, I9)') 'FMAP: ',
     &    'Remaining independent grid points = ', nIndep
      end if
C
      do i = 1, 3
        if (GridMis(i) .ne. 0) then
          write(6, '(1X, 3A)')
     &      'FMAP-WRN: Some origin shift mates not on grid ',
     &      xyz(i:i), '.'
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine SetColFac(Rn, STBF, LCM, ColFac, GridMis)
      implicit none
C
C I/O
      integer  Rn(3), STBF, LCM, ColFac(4), GridMis(3)
C
C local
      integer  i
C
C external
      integer   iLCM
      external  iLCM
C
C begin
      LCM = STBF
C
      do i = 1, 3
        LCM = iLCM(LCM, Rn(i))
      end do
C
      do i = 1, 3
        ColFac(i) = LCM / Rn(i)
      end do
C
      ColFac(4) = LCM / STBF
C
      do i = 1, 3
        GridMis(i) = 0
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine SymEquiv(SymMx, mSymMx, iSymMx,
     &                    LCM, ColFac, Ri, Re)
      implicit none
C I/O
      integer  mSymMx, SymMx(mSymMx, 3, 4), iSymMx
      integer  LCM, ColFac(4)
      integer  Ri(3), Re(3)
C
C Determine the indices Re for the given indices Ri and symmetry operation
C iSymMx.
C
C local
      integer  i
C
C begin
      do i = 1, 3
        Re(i) =   SymMx(iSymMx, i, 1) * ColFac(1) * Ri(1)
     &          + SymMx(iSymMx, i, 2) * ColFac(2) * Ri(2)
     &          + SymMx(iSymMx, i, 3) * ColFac(3) * Ri(3)
     &          + SymMx(iSymMx, i, 4) * ColFac(4)
        Re(i) = mod(Re(i), LCM)
        if (Re(i) .lt. 0) Re(i) = Re(i) + LCM
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine MarkSymDep(FlagMap, Rn, Ri, iMap,
     &                      SymMx, mSymMx, iSymMx,
     &                      LCM, ColFac, GridMis)
      implicit none
C I/O
      integer  FlagMap(*), Rn(3), Ri(3), iMap
      integer  mSymMx, SymMx(mSymMx, 3, 4), iSymMx
      integer  LCM, ColFac(4), GridMis(3)
C
C Mark sym. equiv. grid points for a given pivot grid point Ri.
C
C local
      integer  eMap, i
      integer  Re(3)
      logical  OK
C
C external
      integer   f3dix
      external  f3dix
C
C begin
      call SymEquiv(SymMx, mSymMx, iSymMx,
     &              LCM, ColFac, Ri, Re)
C
      OK = .true.
C
      do i = 1, 3
        if (mod(Re(i), ColFac(i)) .eq. 0) then
          Re(i) = Re(i) / ColFac(i)
        else
          GridMis(i) = 1
          OK = .false.
        end if
      end do
C
      if (OK) then
        eMap = f3dix(Re, Rn)
        do while (FlagMap(eMap) .ne. 0)
          eMap = FlagMap(eMap)
        end do
        if (eMap .ne. iMap) FlagMap(eMap) = iMap
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine CoreMapFM2RhoMask(i, n, f, l, j)
      implicit none
C I/O
      integer  i, n, f, l, j
C
C Map one index in Flagmap back to RhoMask.
C
C begin
      j = i
      if      (j .gt. l) then
        j = j - n
      else if (j .eq. l) then
        if (f + n .eq. l) j = f
      else if (j .lt. f) then
        j = j + n
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine MapFM2RhoMask(Ri, Rn,
     &                         SMAASY, SNAASY,
     &                         SMBASY, SNBASY,
     &                         SMCASY, SNCASY,
     &                         ja, jb, jc)
      implicit none
C I/O
      integer  Ri(3), Rn(3)
      integer  SMAASY, SNAASY
      integer  SMBASY, SNBASY
      integer  SMCASY, SNCASY
      integer  ja, jb, jc
C
C Map independent points in Flagmap back to RhoMask.
C
C begin
      call CoreMapFM2RhoMask(Ri(1), Rn(1), SMAASY, SNAASY, ja)
      call CoreMapFM2RhoMask(Ri(2), Rn(2), SMBASY, SNBASY, jb)
      call CoreMapFM2RhoMask(Ri(3), Rn(3), SMCASY, SNCASY, jc)
C
      return
      end
C}
C=======================================================================
C{
      subroutine MarkSymEquiv(FlagMap, Rn,
     &                        SymMx, mSymMx, nSymMx, STBF,
     &                        SMAASY, SMBASY, SMCASY,
     &                        SNAASY, SNBASY, SNCASY,
     &                        SRHOMA, SNMASK,
     &                        WrnLev)
      implicit none
C I/O
      integer  FlagMap(*), Rn(3)
      integer  mSymMx, SymMx(mSymMx, 3, 4), nSymMx, STBF
      integer  SMAASY, SMBASY, SMCASY
      integer  SNAASY, SNBASY, SNCASY
      integer  SRHOMA(SMAASY:SNAASY, SMBASY:SNBASY, SMCASY:SNCASY)
      integer  SNMASK
      integer  WrnLev
C
C Loop all independent grid points and mark sym. equiv. grid points.
C
C local
      integer    PRn, iMap, ia, ib, ic, iSymMx, i
      integer    ja, jb, jc
      integer    Ri(3)
      integer    nIndep, LCM, ColFac(4), GridMis(3)
      character  xyz*3
C
C external
      integer   iLCM, f3dix, iModPositive
      external  iLCM, f3dix, iModPositive
C
C begin
      xyz = 'xyz'
C
      if (WrnLev .ge. 5)
     &  write(6, '(1X, 2A)') 'FMAP: ',
     &    'Removing grid points related by symmetry operations'
C
      call SetColFac(Rn, STBF, LCM, ColFac, GridMis)
C
      PRn = Rn(1) * Rn(2) * Rn(3)
C
      do ic = SMCASY, SNCASY
      do ib = SMBASY, SNBASY
      do ia = SMAASY, SNAASY
        if (SRHOMA(ia, ib, ic) .ne. 0) then
          Ri(1) = iModPositive(ia, Rn(1))
          Ri(2) = iModPositive(ib, Rn(2))
          Ri(3) = iModPositive(ic, Rn(3))
          iMap = f3dix(Ri, Rn)
          do iSymMx = 1, nSymMx
            call MarkSymDep(FlagMap, Rn, Ri, iMap,
     &                      SymMx, mSymMx, iSymMx,
     &                      LCM, ColFac, GridMis)
          end do
        end if
      end do
      end do
      end do
C
C Consistency check
      do ic = SMCASY, SNCASY
      do ib = SMBASY, SNBASY
      do ia = SMAASY, SNAASY
        if (SRHOMA(ia, ib, ic) .ne. 0) then
          Ri(1) = iModPositive(ia, Rn(1))
          Ri(2) = iModPositive(ib, Rn(2))
          Ri(3) = iModPositive(ic, Rn(3))
          call MapFM2RhoMask(Ri, Rn,
     &                       SMAASY, SNAASY,
     &                       SMBASY, SNBASY,
     &                       SMCASY, SNCASY,
     &                       ja, jb, jc)
          if (ia .ne. ja .or. ib .ne. jb .or. ic .ne. jc) then
            write(6, '(1X, A)') 'FATAL ERROR: MapFM2RhoMask fails.'
            call WrnDie(-5, 'FMAP', 'Fatal Coding Error.')
            call Die
          end if
          iMap = f3dix(Ri, Rn)
          if (FlagMap(iMap) .ne. 0) then
            call WrnDie(-5, 'FMAP',
     &        'RHOMASK not consistent with local FMAP symmetry.')
            call Die
          end if
        end if
      end do
      end do
      end do
C
      nIndep = 0
      do iMap = 1, PRn
        if (FlagMap(iMap) .eq. 0) nIndep = nIndep + 1
      end do
C
      if (WrnLev .ge. 5) then
        write(6, '(1X, 2A, I9)') 'FMAP: ',
     &    '      Total number of grid points = ', PRn
        write(6, '(1X, 2A, I9)') 'FMAP: ',
     &    'Remaining independent grid points = ', nIndep
      end if
C
      do i = 1, 3
        if (GridMis(i) .ne. 0) then
          write(6, '(1X, 3A)')
     &      'FMAP-WRN: Some symmetry mates not on grid ',
     &      xyz(i:i), '.'
        end if
      end do
C
      if (SNMASK .ne. 0 .and. nIndep .ne. SNMASK) then
        call WrnDie(-5, 'FMAP', 'Fatal Coding Error.')
        call Die
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine MarkAddl(FlagMap, Rn,
     &                    AddlGen, mAddlGen, nAddlGen, STBF,
     &                    WrnLev)
      implicit none
C I/O
      integer  FlagMap(*), Rn(3)
      integer  mAddlGen, AddlGen(mAddlGen, 3, 4), nAddlGen, STBF
      integer  WrnLev
C
C Loop all independent grid points and mark grid points related
C by AddlGen(erators).
C
C local
      integer    PRn, iMap, ia, ib, ic, iAddlGen, i
      integer    Ri(3)
      integer    nIndep, LCM, ColFac(4), GridMis(3)
      character  xyz*3
C
C external
      integer   iLCM
      external  iLCM
C
C begin
      xyz = 'xyz'
C
      if (WrnLev .ge. 5)
     &  write(6, '(1X, 3A)') 'FMAP: ',
     &    'Removing grid points related by',
     &    ' operations due to ADDLgenerators'
C
      call SetColFac(Rn, STBF, LCM, ColFac, GridMis)
C
      PRn = Rn(1) * Rn(2) * Rn(3)
C
      do iAddlGen = 1, nAddlGen
        iMap = 1
        do ic = 0, Rn(3) - 1
           Ri(3) = ic
        do ib = 0, Rn(2) - 1
           Ri(2) = ib
        do ia = 0, Rn(1) - 1
           Ri(1) = ia
          if (FlagMap(iMap) .eq. 0) then
            call MarkSymDep(FlagMap, Rn, Ri, iMap,
     &                      AddlGen, mAddlGen, iAddlGen,
     &                      LCM, ColFac, GridMis)
          end if
          iMap = iMap + 1
        end do
        end do
        end do
      end do
C
      nIndep = 0
      do iMap = 1, PRn
        if (FlagMap(iMap) .eq. 0) nIndep = nIndep + 1
      end do
C
      if (WrnLev .ge. 5) then
        write(6, '(1X, 2A, I9)') 'FMAP: ',
     &    '      Total number of grid points = ', PRn
        write(6, '(1X, 2A, I9)') 'FMAP: ',
     &    'Remaining independent grid points = ', nIndep
      end if
C
      do i = 1, 3
        if (GridMis(i) .ne. 0) then
          write(6, '(1X, 3A)')
     &      'FMAP-WRN: Some symmetry mates not on grid ',
     &      xyz(i:i), '.'
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      integer function OptFlagMap(FlagMap, nFlagMap)
      implicit none
C I/O
      integer  FlagMap(*), nFlagMap
C
C Optimize FlagMap: If FlagMap(i) > 0, always point to final
C independent grid point. In each dependency group, make the
C grid point with the lowest index the independent grid point.
C
C local
      integer  iMap, jMap
      logical  cloop
C
C begin
      OptFlagMap = 0
C
      do iMap = 1, nFlagMap
        if (FlagMap(iMap) .gt. 0) then
          jMap = iMap
                    cloop = .true.
          do while (cloop)
            jMap = FlagMap(jMap)
            if (FlagMap(jMap) .le. 0) cloop = .false.
          end do
          FlagMap(iMap) = jMap
        else
          OptFlagMap = OptFlagMap + 1
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      logical function IsNewAddlGen(mMx, Mx, n)
      implicit none
C I/O
      integer  mMx, Mx(mMx, 3, 4), n
C
C local
      integer  i, ir, ic
      logical  Same
C
C begin
      IsNewAddlGen = .false.
C
      do i = 1, n - 1
        Same = .true.
        do ir = 1, 3
        do ic = 1, 4
          if (Mx(i,ir,ic) .ne. Mx(n,ir,ic)) Same = .false.
        end do
        end do
        if (Same) return
      end do
C
      IsNewAddlGen = .true.
C
      return
      end
C}
C=======================================================================
C{
      subroutine AddInvAddlGen(mMx, Mx, i, j, STBF, OK)
      implicit none
C I/O
      integer  mMx, Mx(mMx, 3, 4), i, j, STBF
      logical  OK
C
C Compute and add inverse of input generators.
C
C local
      integer  det, trace, ir, ic
C
C external
      integer   iModPositive
      external  iModPositive
      logical   IsNewAddlGen
      external  IsNewAddlGen
C
C begin
      Mx(j,1,1) =  Mx(i,2,2) * Mx(i,3,3) - Mx(i,2,3) * Mx(i,3,2)
      Mx(j,1,2) = -Mx(i,1,2) * Mx(i,3,3) + Mx(i,1,3) * Mx(i,3,2)
      Mx(j,1,3) =  Mx(i,1,2) * Mx(i,2,3) - Mx(i,1,3) * Mx(i,2,2)
      Mx(j,2,1) = -Mx(i,2,1) * Mx(i,3,3) + Mx(i,2,3) * Mx(i,3,1)
      Mx(j,2,2) =  Mx(i,1,1) * Mx(i,3,3) - Mx(i,1,3) * Mx(i,3,1)
      Mx(j,2,3) = -Mx(i,1,1) * Mx(i,2,3) + Mx(i,1,3) * Mx(i,2,1)
      Mx(j,3,1) =  Mx(i,2,1) * Mx(i,3,2) - Mx(i,2,2) * Mx(i,3,1)
      Mx(j,3,2) = -Mx(i,1,1) * Mx(i,3,2) + Mx(i,1,2) * Mx(i,3,1)
      Mx(j,3,3) =  Mx(i,1,1) * Mx(i,2,2) - Mx(i,1,2) * Mx(i,2,1)
C
      OK = .false.
C
      det =   Mx(i,1,1) * Mx(j,1,1)
     &      + Mx(i,1,2) * Mx(j,2,1)
     &      + Mx(i,1,3) * Mx(j,3,1)
C
      if (abs(det) .ne. 1) then
        call WrnDie(0, 'FMAP', 'Illegal ADDLgenerators.')
        return
      end if
C
      trace = Mx(i,1,1) + Mx(i,2,2) + Mx(i,3,3)
C
      if (trace .ne. -det .and. trace .ne. -3) then
        call WrnDie(0, 'FMAP',
     &       'Program restriction: Rotation parts of ADDLgenerators'
     &    // ' can only be -1, 2, m.')
        return
      end if
C
      OK = .true.
C
      if (det .eq. -1) then
        do ir = 1, 3
        do ic = 1, 3
          Mx(j,ir,ic) = -Mx(j,ir,ic)
        end do
        end do
      end if
C
      do ir = 1, 3
        Mx(j,ir,4) = 0
        do ic = 1, 3
          Mx(j,ir,4) = Mx(j,ir,4) - Mx(j,ir,ic) * Mx(i,ic,4)
        end do
        Mx(j,ir,4) = iModPositive(Mx(j,ir,4), STBF)
      end do
C
      if (IsNewAddlGen(mMx, Mx, j)) j = j + 1
C
      return
      end
C}
C=======================================================================
C{
      subroutine IniFMCommon(ColdStart, STBF)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'flagmap.inc'
      logical  ColdStart
      integer  STBF
C
C Initialize flagmap.inc common blocks.
C
C local
      integer  i, j
C
C begin
C mrt changed for the sake of run-time checks (lf95 --chk eux)
      if (.not. ColdStart ) then
        if( ptFlagMap .ne. 0) then
          call FreHp(ptFlagMap, integ4(nFlagMap))
        end if
      end if
C
      ptFlagMap = 0
      nFlagMap = 0
      do j = 1, 3
      do i = 1, 2
        nAsyFM(i, j) = 0
      end do
      end do
      call InitSym(SymMx, ITSymMx, mSymMx, nSymMx, STBF)
      n_ssVM = 0
      n1AddlGen = 0
      n2AddlGen = 0
      UseSym  = .false.
      Use_ss  = .false.
      UseAddl = .false.
      ValidFlagMap = .false.
C
      return
      end
C}
C=======================================================================
C{
      subroutine CopyFMtoRSObj(DensMap, FlagMap, PRn)
      implicit none
C I/O
      real       DensMap(*)
      integer    FlagMap(*), PRn
C
C Copy FlagMap to real-space object.
C
C local
      integer  i
C
C begin
      do i = 1, PRn
        if      (FlagMap(i) .eq. 0) then
          DensMap(i) = dfloat( 0)
        else if (FlagMap(i) .lt. 0) then
          DensMap(i) = dfloat(-1)
        else
          DensMap(i) = dfloat( 1)
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine ReadExMap(FileName, n, w)
      implicit none
C I/O
      character  FileName*(*)
      integer    n
      real       w(*)
C
C Read external map, one value per line (useful for debugging)
C
C local
      integer  iUnit, i
      logical  Err
C
C begin
      call AssFil(FileName, iUnit, 'READ', 'FORMATTED', Err)
      if (Err) call Die
C
      do i = 1, n
        read(iUnit, *, err = 91, end = 92) w(i)
      end do
C
      call VClose(iUnit, 'KEEP', Err)
      if (Err) call Die
      return
C
 91   continue
      call WrnDie(-5, 'FMAP', 'Error reading ExMap.')
      call Die
C
 92   continue
      call WrnDie(-5, 'FMAP', 'Unexpected end-of-file reading ExMap.')
      call Die
C
      end
C}
C=======================================================================
C{
      subroutine fmappa(ActionType, iRSObj, ExMap,
     &                  ChkInDep, MaxInDep, BailOut,
     &                  SymMx, ITSymMx, mSymMx, nSymMx, STBF,
     &                  mAddlGen, n1AddlGen, n2AddlGen,
     &                  AddlGen, ITAddlGen,
     &                  n_ssVM, ssV, ssM,
     &                  UseSym, Use_ss, UseAddl, ValidFlagMap,
     &                  xRhoNum, xRhoNam)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      character      ActionType*(*)
      integer        iRSObj
      character      ExMap*(*)
      logical        ChkInDep, MaxInDep, BailOut
      integer        mSymMx, nSymMx, STBF
      integer        SymMx(mSymMx, 3, 4), ITSymMx(mSymMx, 3, 3)
      integer        n_ssVM, ssV(3, 3), ssM(3)
      integer        mAddlGen, n1AddlGen, n2AddlGen
      integer        AddlGen(mAddlGen, 3, 4), ITAddlGen(mAddlGen, 3, 4)
      logical        UseSym, Use_ss, UseAddl, ValidFlagMap
      integer        xRhoNum
      character*(*)  xRhoNam(*)
C
C This is the parser for fmap
C
C local
      integer  i, j, V(4), nV, iAddlGen
      logical  HadSym, Had_ssVM, Had_AddlG, TmpLo, OK
C
      integer   StrMax
      parameter(StrMax = 132)
      character*(StrMax)  Str
C
C external
      integer   iModPositive
      external  iModPositive
      logical   IsNewAddlGen
      external  IsNewAddlGen
C
C begin
      ActionType = ' '
      iRSObj = 0
      ExMap = ' '
      ChkInDep = .true.
      MaxInDep = .true.
      BailOut = .true.
C
      HadSym = .false.
      Had_ssVM = .false.
      Had_AddlG = .false.
C
C parsing
      CALL PUSEND('FMAP>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('FMAP>')
      CALL MISCOM('FMAP>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-fmap')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ACTI') THEN
        CALL NEXTST('ACTIon=', WD)
        if (     WDLEN .lt. 4
     &      .or. (      WD(1:4) .ne. 'BUIL'
     &            .and. WD(1:4) .ne. 'MAPI'
     &            .and. WD(1:4) .ne. 'COPY'
     &            .and. WD(1:4) .ne. 'RESE')) then
          CALL DSPERR('FMAP', 'unknown qualifier')
        else
          ActionType = WD(1:4)
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'RSOB') THEN
        call GetiObj(iRSObj, xRhoNum, xRhoNam, ' ',
     &               'RSOBject=', 'real', 'FMAP')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'EXMA') THEN
        CALL NEXTST('EXMAp=', WD)
        ExMap = ' '
        if (WDLEN .gt. 0) ExMap = WD(1:WDLEN)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SYMM') THEN
        CALL NEXTWD('SYMMetry=')
        IF (WD(1:4).EQ.'=   ') CALL NEXTWD('SYMMetry=')
        if (wd(1:1) .eq. '?') then
          write(6,'(2A)') ' --------------FMAP-symmetry',
     &     '----------------------------------------------------'
          call xSyPri(mSymMx, nSymMx, SymMx, ITSymMx, STBF)
          write(6,'(2A)') ' ---------------------------',
     &     '----------------------------------------------------'
        else
          if (.not. HadSym) then
            call InitSym(SymMx, ITSymMx, mSymMx, nSymMx, STBF)
            HadSym = .true.
            ValidFlagMap = .false.
          end if
          CALL SAVEWD
          CALL NEXTEX('SYMMetry=',WDD,WDMAX,WDDLEN)
          call xrSyPA(WDD, WDDLEN,
     &                mSymMx, nSymMx, SymMx, ITSymMx, STBF)
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ASYM') THEN
        CALL NEXTEX('ASYMm=', Str, StrMax, WDDLEN)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SSVM') THEN
        if (.not. Had_ssVM) then
          n_ssVM = 0
          Had_ssVM = .true.
          ValidFlagMap = .false.
        end if
        CALL NEXTST('SSVM=', WD)
        call IVecPa(WD(1:WDLEN), V, 4, nV)
        if      (nV .ne. 4) then
          call WrnDie(0, 'FMAP', 'Illegal format for SSVM.')
        else if (V(1) .eq. 0 .and. V(2) .eq. 0 .and. V(3) .eq. 0) then
          call WrnDie(0, 'FMAP', 'Illegal s.s.Vector.')
        else if (V(4) .lt. 0 .or. V(4) .eq. 1) then
          call WrnDie(0, 'FMAP', 'Illegal s.s.Modulus.')
        else if (n_ssVM .eq. 3) then
          call WrnDie(0, 'FMAP', 'Too many SSVM.')
        else
          n_ssVM = n_ssVM + 1
          do i = 1, 3
            ssV(i, n_ssVM) = V(i)
          end do
          ssM(n_ssVM) = V(4)
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ADDL') THEN
        CALL NEXTWD('ADDLgenerator=')
        IF (WD(1:4).EQ.'=   ') CALL NEXTWD('ADDLgenerator=')
        if (wd(1:1) .eq. '?') then
          write(6,'(2A)') ' --------------FMAP-ADDLgene',
     &     'rators----------------------------------------------'
          if (wd(2:2) .eq. '?') then
            call xSyPri(mAddlGen, n2AddlGen, AddlGen, ITAddlGen, STBF)
          else
            call xSyPri(mAddlGen, n1AddlGen, AddlGen, ITAddlGen, STBF)
          end if
          write(6,'(2A)') ' ---------------------------',
     &     '----------------------------------------------------'
        else
          if (.not. Had_AddlG) then
            n1AddlGen = 0
            n2AddlGen = 0
            Had_AddlG = .true.
            ValidFlagMap = .false.
          end if
          CALL SAVEWD
          CALL NEXTST('ADDLgenerator=', WD)
          call ToUpper(WD(:WDLEN))
          if (n1AddlGen .ge. mAddlGen / 2) then
            call WrnDie(0, 'FMAP', 'Too many additional generators.')
          else
            call xrSyPA(WD, WDLEN,
     &                  mAddlGen, n1AddlGen, AddlGen, ITAddlGen, STBF)
C
            if (n1AddlGen .gt. 0) then
              do i = 1, 3
                AddlGen(n1AddlGen,i,4)
     &            = iModPositive(AddlGen(n1AddlGen,i,4), STBF)
              end do
            end if
C
            if (.not. IsNewAddlGen(mAddlGen, AddlGen, n1AddlGen))
     &        n1AddlGen = n1AddlGen - 1
C
            n2AddlGen = n1AddlGen + 1
C
                      iAddlGen = 1
            do while (iAddlGen .le. n1AddlGen)
              call AddInvAddlGen(mAddlGen, AddlGen,
     &                           iAddlGen, n2AddlGen, STBF, OK)
              if (.not. OK) then
                n2AddlGen = n1AddlGen
                n1AddlGen = n1AddlGen - 1
                iAddlGen = 0
              end if
              iAddlGen = iAddlGen + 1
            end do
C
            n2AddlGen = n2AddlGen - 1
          end if
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'USES') THEN
        CALL NEXTLO('USESymmetry=', TmpLo)
        if (      TmpLo .and. .not. UseSym) ValidFlagMap = .false.
        if (.not. TmpLo .and.       UseSym) ValidFlagMap = .false.
        UseSym = TmpLo
C=====================================================================
      ELSE IF (WD(1:4).EQ.'USE_') THEN
        CALL NEXTLO('USE_ss=', TmpLo)
        if (      TmpLo .and. .not. USE_ss) ValidFlagMap = .false.
        if (.not. TmpLo .and.       USE_ss) ValidFlagMap = .false.
        USE_ss = TmpLo
C=====================================================================
      ELSE IF (WD(1:4).EQ.'USEA') THEN
        CALL NEXTLO('USEAddl=', TmpLo)
        if (      TmpLo .and. .not. UseAddl) ValidFlagMap = .false.
        if (.not. TmpLo .and.       UseAddl) ValidFlagMap = .false.
        UseAddl = TmpLo
C=====================================================================
      ELSE IF (WD(1:4).EQ.'CHKI') THEN
        CALL NEXTLO('CHKIndep=', ChkInDep)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MAXI') THEN
        CALL NEXTLO('MAXIndep=', MaxInDep)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'BAIL') THEN
        CALL NEXTLO('BAILout=', BailOut)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
        WRITE(6,'(2A)') ' --------------FMAP-parameters',
     &   '--------------------------------------------------'
        if (ActionType .eq. 'BUIL')
     &    write(6, '(1X, A)') 'ActionType=BUILd'
        if (ActionType .eq. 'MAPI')
     &    write(6, '(1X, A)') 'ActionType=MAPIndependent'
        if (ActionType .eq. 'COPY')
     &    write(6, '(1X, A)') 'ActionType=COPY'
        if (ActionType .eq. 'RESE')
     &    write(6, '(1X, A)') 'ActionType=RESEt'
        if (iRSObj .gt. 0)
     &    call EchoObj(iRSObj, xRhoNam, 'RSOBject=')
        write(6, '(1X, A, I3)') 'Number of symmetry operations = ',
     &    nSymMx
        write(6, '(1X, A, I3)') 'Number of s.s.Vectors & Moduli = ',
     &    n_ssVM
        if (n_ssVM .gt. 0) then
          write(6, '(3X, A)') 's.s.Vector  Modulus'
          do i = 1, n_ssVM
            write(6, '(3X, 3(1X,I2), 1X, I3)')
     &        (ssV(j, i), j = 1, 3), ssM(i)
          end do
        end if
        write(6, '(1X, A, I3)')
     &    'Number of ADDLgenerators = ', n1AddlGen
        if (UseSym) then
          write(6, '(1X, A)') 'USESym=TRUE'
        else
          write(6, '(1X, A)') 'USESym=FALSE'
        end if
        if (Use_ss) then
          write(6, '(1X, A)') 'USE_ss=TRUE'
        else
          write(6, '(1X, A)') 'USE_ss=FALSE'
        end if
        if (UseAddl) then
          write(6, '(1X, A)') 'USEAddl=TRUE'
        else
          write(6, '(1X, A)') 'USEAddl=FALSE'
        end if
        WRITE(6,'(2A)') ' -----------------------------',
     &   '--------------------------------------------------'
C=====================================================================
      ELSE
        CALL CHKEND('FMAP>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      return
      end
C}
C=======================================================================
C{
      subroutine fmap(Mode, na, nb, nc, nrRho, niRho, QHerm,
     &                SMAASY, SMBASY, SMCASY,
     &                SNAASY, SNBASY, SNCASY,
     &                SHPRHOMA, SNMASK,
     &                xRhoNum, xRhoNam, hprrho, hpirho,
     &                xrNsym, STBF)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'flagmap.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      integer        Mode
      integer        na, nb, nc, nrRho, niRho
      logical        QHerm
      integer        SMAASY, SMBASY, SMCASY
      integer        SNAASY, SNBASY, SNCASY
      integer        SHPRHOMA, SNMASK
      integer        xRhoNum
      character*(*)  xRhoNam(*)
      integer        hprrho(*), hpirho(*)
      integer        xrNsym, STBF
C
C XRAY FMAP main procedure.
C
C local
      character  ActionType*4
      integer    iRSObj
      character  ExMap*128
      logical    ChkInDep, MaxInDep, BailOut
      integer    PRn
C
      double precision  tmEntr, tmExit, tm0, tm1
      double precision  dbprec
      double complex    dbcomp
C
C parameters
      real      SZero
      parameter(SZero=0.0)
C
C external
      integer   lnblnk, OptFlagMap
      external  lnblnk, OptFlagMap
C
C begin
      if (Mode .eq. 0) then
C Initialize the non-volatile settings and return
        call IniFMCommon(.true., STBF)
        return
      end if
C
      if (Mode .lt. 0) then
C Deallocate FlagMap and return
        call IniFMCommon(.false., STBF)
        return
      end if
C
C Call the parser
      call fmappa(ActionType, iRSObj, ExMap,
     &            ChkInDep, MaxInDep, BailOut,
     &            SymMx, ITSymMx, mSymMx, nSymMx, STBF,
     &            mAddlGen, n1AddlGen, n2AddlGen,
     &            AddlGen, ITAddlGen,
     &            n_ssVM, ssV, ssM,
     &            UseSym, Use_ss, UseAddl, ValidFlagMap,
     &            xRhoNum, xRhoNam)
C
      if (Timer .gt. 0) call VCPU(tmEntr)
C
      if (.not. ValidFlagMap) then
        if (ptFlagMap .ne. 0) call FreHp(ptFlagMap, integ4(nFlagMap))
        ptFlagMap = 0
        nFlagMap = 0
      end if
C
      if (ActionType .eq. ' ') return
C
      if (ActionType .eq. 'RESE') then
        call IniFMCommon(.false., STBF)
        return
      end if
C
      if (      (     ActionType .eq. 'MAPI'
     &           .or. ActionType .eq. 'COPY')
     &    .and. iRSObj .lt. 1) then
        call WrnDie(0, 'FMAP', 'Real space object not specified.')
        return
      end if
C
      if (iRSObj .ge. 1) then
        if (hprrho(iRSObj) .eq. 0) then
          if (ActionType .eq. 'MAPI') then
            write(6, '(3A)') ' %FMAP-ERR: Real space object "',
     &        xrhonam(iRSObj)(1:lnblnk(xrhonam(iRSObj))),
     &        '" not defined.'
            call WrnDie(0, 'FMAP', 'object undefined.')
            return
          end if
          call xMapAl(hprrho(iRSObj), hpirho(iRSObj), QHerm,
     &                nrRho, niRho)
        end if
C
        if (ExMap .ne. ' ')
     &    call ReadExMap(ExMap, PRn, heap(hprrho(iRSObj)))
      end if
C
      if (xrNsym .ne. 1) then
        call WrnDie(0, 'FMAP', 'This module can only be used in P1.')
        return
      end if
C
      if (ValidFlagMap) then
        if (     RnFM(1) .ne. na
     &      .or. RnFM(2) .ne. nb
     &      .or. RnFM(3) .ne. nc
     &      .or. nAsyFM(1, 1) .ne. SMAASY
     &      .or. nAsyFM(2, 1) .ne. SNAASY
     &      .or. nAsyFM(1, 2) .ne. SMBASY
     &      .or. nAsyFM(2, 2) .ne. SNBASY
     &      .or. nAsyFM(1, 3) .ne. SMCASY
     &      .or. nAsyFM(2, 3) .ne. SNCASY)
     &    ValidFlagMap = .false.
      end if
C
      RnFM(1) = na
      RnFM(2) = nb
      RnFM(3) = nc
      PRn = RnFM(1) * RnFM(2) * RnFM(3)
      nAsyFM(1, 1) = SMAASY
      nAsyFM(2, 1) = SNAASY
      nAsyFM(1, 2) = SMBASY
      nAsyFM(2, 2) = SNBASY
      nAsyFM(1, 3) = SMCASY
      nAsyFM(2, 3) = SNCASY
C
      if (nrRho .ne. PRn) then
        call WrnDie(-5, 'FMAP', 'Fatal Coding Error.')
        call Die
      end if
C
      if (ValidFlagMap) then
        if (ActionType .eq. 'BUIL' .and. WrnLev .ge. 5)
     &    write(6, '(1X, A)') 'FMAP: FlagMap still valid.'
      else
        if (WrnLev .ge. 5) then
          if (ActionType .eq. 'MAPI' .or. ActionType .eq. 'COPY') then
            write(6, '(1X, A)')
     &        'FMAP: FlagMap undefined or out of date.'
          end if
          write(6, '(1X, A)') 'FMAP: Building FlagMap'
        end if
C
        if (Timer .gt. 0) call VCPU(tm0)
C
C Allocate space for FlagMap and set all flags to zero
        if (ptFlagMap .ne. 0) call FreHp(ptFlagMap, integ4(nFlagMap))
        nFlagMap = PRn
        ptFlagMap = AllHp(integ4(nFlagMap))
        call Fill4(heap(ptFlagMap), PRn, 0)
C
        if (UseSym) then
C Remove grid points related by symmetry operations
          call MarkSymEquiv(heap(ptFlagMap), RnFM,
     &                      SymMx, mSymMx, nSymMx, STBF,
     &                      SMAASY, SMBASY, SMCASY,
     &                      SNAASY, SNBASY, SNCASY,
     &                      heap(SHPRHOMA), SNMASK,
     &                      WrnLev)
        end if
C
        if (Use_ss) then
C Remove grid points related by allowed origin shifts
          call MarkOriginShifts(heap(ptFlagMap), RnFM,
     &                          ssV, ssM, n_ssVM, WrnLev)
        end if
C
        if (UseAddl) then
C Remove grid points related by additional generators
          call MarkAddl(heap(ptFlagMap), RnFM,
     &                  AddlGen, mAddlGen, n2AddlGen, STBF,
     &                  WrnLev)
        end if
C
        if (UseSym .or. Use_ss .or. UseAddl) then
          nIndFM = OptFlagMap(heap(ptFlagMap), nFlagMap)
        else
          nIndFM = nFlagMap
        end if
C
        ValidFlagMap = .true.
C
        if (Timer .gt. 0) then
          call VCPU(tm1)
          write(6, '(1X, A, F10.4, A)')
     &      'FMAP: CPU-time: ',
     &      tm1 - tm0, ' seconds Building FlagMap'
        end if
      end if
C
      if (WrnLev .ge. 5)
     &  write(6, '(1X, 2A, I9)') 'FMAP: ',
     &    'Independent grid points = ', nIndFM
C
      dbprec = dfloat(nIndFM)
      call declar('RESULT', 'DP', ' ', dbcomp, dbprec)
C
      if (iRSObj .ge. 1 .and. ActionType .eq. 'MAPI') then
        call MapInd(heap(hprrho(iRSObj)), heap(ptFlagMap), RnFM,
     &              ChkInDep, MaxInDep, BailOut, WrnLev, 'FMAP')
      end if
C
      if (iRSObj .ge. 1 .and. ActionType .eq. 'COPY') then
        call CopyFMtoRSObj(heap(hprrho(iRSObj)), heap(ptFlagMap), PRn)
        if (QHerm) call FillR4(heap(hpirho(iRSObj)), niRho, SZERO)
      end if
C
      if (Timer .gt. 0) then
        call VCPU(tmExit)
        write(6, '(1X, A, F10.4, A)')
     &    'FMAP: CPU-time: ', tmExit - tmEntr, ' seconds total'
      end if
C
      return
      end
C}
