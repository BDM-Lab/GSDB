C{
      subroutine PeakSearch(DensMap, FlagMap, Rn,
     &                      Level, MaxPeaks, DensCutOff, ToTPeaks,
     &                      WrnLev)
      implicit none
C I/O
      real     DensMap(*)
      integer  FlagMap(*), Rn(3)
      integer  Level, MaxPeaks
      real     DensCutOff
      integer  TotPeaks
      integer  WrnLev
C
C Peak search in P1 with symmetry and/or origin shift relations
C defined by FlagMap.
C
C   I: DensMap    = map to be search for peaks, e.g. electron density etc.
C   I: FlagMap    < 1 => independent grid point, >= 1 => dependent grid point
C   O:    "       = -1 for independent grid points which are peaks
C                 =  0 for independent grid points which are not peaks
C   I: Rn         = dimensions of DensMap/FlagMap
C   I: Level      < 1 determine the DensCutOff for MaxPeaks highest grid points
C                 = 1 look at 6 nearest neighbors of each grid point
C                 = 2 also look at the 12 2nd nearest neighbors
C                 > 2 also look at the 8 3rd nearest neighbors
C                   (For a simple illustration of the nearest neighbor
C                   relations refer to C. Giacovazzo, Ed. (1992),
C                   Fundamentals of Crystallography, p. 90)
C   I: MaxPeaks   = maximum number of peaks to be stored later
C   O: DensCutOff = Cut-off value such that <= MaxPeaks peaks will be
C                   collected later
C   O: TotPeaks   = total number of peaks
C
C DensMap is scanned for peaks in a 3-deep loop.
C Each time a peak is found, its value is added to a histogram of
C the peak heights.
C At the end of the subroutine, DensCutOff is determined from the
C histogram.
C
C In 1D, a grid point would be a peak if the DensMap value left to it is
C less than or equal to the value at the pivot grid point, and the value
C right to it is also less then or equal to the value at the pivot grid
C point.
C The 3D generalization of this definition for a peak is used here.
C
C local
      integer          nk
      integer  im, jm, km
      integer  i0, j0, k0
      integer  ip, jp, kp
      integer  ibreak, jbreak, kbreak
      integer  nj_nk, ni_nj_nk
C
      integer   nBox, iLastBox
      parameter(nBox = 1024, iLastBox = nBox - 1)
      integer    Box(0:iLastBox)
      integer   iBox, cum
      real      dBox, dDens
C
      integer  iMap, iAsy
      real     Dens, MinDens, MaxDens
C
C parameters
      real      SZero
      parameter(SZero = 0.0)
C
C begin
C Initialize histogram
      do iBox = 0, iLastBox
        Box(iBox) = 0
      end do
C
            nk =                 Rn(1)
         nj_nk =         Rn(2) * Rn(1)
      ni_nj_nk = Rn(3) * Rn(2) * Rn(1)
C
      MinDens = DensMap(1)
      MaxDens = DensMap(1)
C
      do iMap = 1, ni_nj_nk
        if (FlagMap(iMap) .le. 0) then
          FlagMap(iMap) = 0
          if (MinDens .gt. DensMap(iMap)) MinDens = DensMap(iMap)
          if (MaxDens .lt. DensMap(iMap)) MaxDens = DensMap(iMap)
        end if
      end do
C
      dBox = (MaxDens - MinDens) / nBox
C
      TotPeaks = 0
      iMap = 1
C
      im = ni_nj_nk - nj_nk
      i0 = 0
      ip = nj_nk
      ibreak = ni_nj_nk
      do while (ip .lt. ibreak)
        jm = nj_nk - nk
        j0 = 0
        jp = nk
        jbreak = nj_nk
        do while (jp .lt. jbreak)
          km = nk - 1
          k0 = 0
          kp = 1
          kbreak = nk
          do while (kp .lt. kbreak)
            if (FlagMap(iMap) .le. 0) then
              iAsy = iMap
            else
              iAsy = FlagMap(iMap)
            end if
C
            if (FlagMap(iAsy) .eq. 0) then
              Dens = DensMap(iMap)
C
              if (Level .ge. 1) then
C               m00 0m0 00m
C               p00 0p0 00p
C
                if (Dens .lt. DensMap(im + j0 + k0 + 1)) goto 20
                if (Dens .lt. DensMap(ip + j0 + k0 + 1)) goto 20
                if (Dens .lt. DensMap(i0 + jm + k0 + 1)) goto 20
                if (Dens .lt. DensMap(i0 + jp + k0 + 1)) goto 20
                if (Dens .lt. DensMap(i0 + j0 + km + 1)) goto 20
                if (Dens .lt. DensMap(i0 + j0 + kp + 1)) goto 20
C
                if (Level .ge. 2) then
C                 mm0 m0m 0mm mp0 m0p 0mp
C                 pp0 p0p 0pp pm0 p0m 0pm
C
                  if (Dens .lt. DensMap(im + jm + k0 + 1)) goto 20
                  if (Dens .lt. DensMap(ip + jp + k0 + 1)) goto 20
                  if (Dens .lt. DensMap(im + j0 + km + 1)) goto 20
                  if (Dens .lt. DensMap(ip + j0 + kp + 1)) goto 20
                  if (Dens .lt. DensMap(i0 + jm + km + 1)) goto 20
                  if (Dens .lt. DensMap(i0 + jp + kp + 1)) goto 20
                  if (Dens .lt. DensMap(im + jp + k0 + 1)) goto 20
                  if (Dens .lt. DensMap(ip + jm + k0 + 1)) goto 20
                  if (Dens .lt. DensMap(im + j0 + kp + 1)) goto 20
                  if (Dens .lt. DensMap(ip + j0 + km + 1)) goto 20
                  if (Dens .lt. DensMap(i0 + jm + kp + 1)) goto 20
                  if (Dens .lt. DensMap(i0 + jp + km + 1)) goto 20
C
                  if (Level .ge. 3) then
C                   mmm mmp mpm mpp
C                   ppp ppm pmp pmm
C
                    if (Dens .lt. DensMap(im + jm + km + 1)) goto 20
                    if (Dens .lt. DensMap(ip + jp + kp + 1)) goto 20
                    if (Dens .lt. DensMap(im + jm + kp + 1)) goto 20
                    if (Dens .lt. DensMap(ip + jp + km + 1)) goto 20
                    if (Dens .lt. DensMap(im + jp + km + 1)) goto 20
                    if (Dens .lt. DensMap(ip + jm + kp + 1)) goto 20
                    if (Dens .lt. DensMap(im + jp + kp + 1)) goto 20
                    if (Dens .lt. DensMap(ip + jm + km + 1)) goto 20
                  end if
                end if
              end if
C
C A peak was found
              FlagMap(iAsy) = -1
              TotPeaks = TotPeaks + 1
C
C Update the histogram
                  dDens = Dens - MinDens
              if (dDens .eq. SZero .or. dDens .lt. dBox) then
                iBox = 0
              else
                         iBox = int(dDens / dBox)
                if      (iBox .gt. iLastBox) then
                         iBox = iLastBox
                else if (iBox .lt. 0) then
                         iBox = 0
                end if
              end if
              Box(iBox) = Box(iBox) + 1
            end if
C
 20         continue
C
            iMap = iMap + 1
C
            km = k0
            k0 = kp
            kp = kp + 1
            if (kp .eq. nk) then
              kp = 0
              kbreak = 1
            end if
          end do
          jm = j0
          j0 = jp
          jp = jp + nk
          if (jp .eq. nj_nk) then
            jp = 0
            jbreak = nk
          end if
        end do
        im = i0
        i0 = ip
        ip = ip + nj_nk
        if (ip .eq. ni_nj_nk) then
          ip = 0
          ibreak = nj_nk
        end if
      end do
C
C Determine DensCutOff from the histogram
      cum = 0
      iBox = nBox
      do while (iBox .gt. 0)
        cum = cum + Box(iBox - 1)
        if (cum .gt. MaxPeaks) goto 30
        iBox = iBox - 1
      end do
 30   continue
C
      if (iBox .eq. 0) then
        DensCutOff =  MinDens
      else
C                                                     ARBITRARY
        DensCutOff = (MinDens + iBox * dBox) + dBox * 1.e-4
C                                            compensate round-off errors
      end if
C
      if (WrnLev .ge. 10)
     &  write(6, '(1X, A, G14.6)') 'PSEARCH: DensCutOff = ', DensCutOff
C
      return
      end
C}
C=======================================================================
C{
      logical function SiteOrderF(ia, ib, Sites,
     &                            j1, j2, j3, j4, j5, j6, j7)
      implicit none
C I/O
      integer           ia, ib
      double precision  Sites(7, *)
      integer           j1, j2, j3, j4, j5, j6, j7
C
C This functions tells the SORT procedure what is larger and
C what is smaller.
C
C local
      integer  i, ii(3)
C
C begin
      ii(1) = 5
      ii(2) = 6
      ii(3) = 4
C
      do i = 1, 3
        SiteOrderF = .false.
        if (Sites(ii(i), ia) .lt. Sites(ii(i), ib)) return
        SiteOrderF = .true.
        if (Sites(ii(i), ia) .gt. Sites(ii(i), ib)) return
      end do
C
      SiteOrderF = .true.
      return
      end
C}
C=======================================================================
C{
      subroutine SiteExChF(ia, ib, Sites, j1, j2, j3, j4, j5, j6, j7)
      implicit none
C I/O
      integer           ia, ib
      double precision  Sites(7, *)
      integer           j1, j2, j3, j4, j5, j6, j7
C
C This subroutine is called from the SORT procedure each time two list
C entries need to be swapped.
C
C local
      integer           i
      double precision  tmp
C
C begin
      do i = 1, 7
        tmp = Sites(i, ia)
        Sites(i, ia) = Sites(i, ib)
        Sites(i, ib) = tmp
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine GetCloseGridPt(Dn, X, FCTrMx,
     &                          DiClose, DistClose)
      implicit none
C I/O
      integer           Dn(3)
      double precision  X(3)
      double precision  FCTrMx(3, 3)
      integer           DiClose(3)
      double precision  DistClose
C
C Find grid point closest to X
C
C local
      integer           Xi(3), Di(3), ia, ib, ic, i
      double precision  Xg(3), Xd(0:1, 3), DC(3), Dist2
      logical           First
C
C begin
      do i = 1, 3
        Xg(i) = mod(X(i), dfloat(1))
        if (Xg(i) .lt. dfloat(0)) Xg(i) = Xg(i) + dfloat(1)
        Xg(i) = Xg(i) * dfloat(Dn(i))
        Xi(i) = int(Xg(i))
        Xd(0, i) = (Xg(i) - dfloat(Xi(i)    )) / dfloat(Dn(i))
        Xd(1, i) = (Xg(i) - dfloat(Xi(i) + 1)) / dfloat(Dn(i))
      end do
C
      DistClose = dfloat(0)
      First = .true.
C
      do ic = 0, 1
        Di(3) = mod(Xi(3) + ic, Dn(3))
      do ib = 0, 1
        Di(2) = mod(Xi(2) + ib, Dn(2))
      do ia = 0, 1
        Di(1) = mod(Xi(1) + ia, Dn(1))
        call TrCoor(FCTrMX,
     &              Xd(ia, 1),
     &              Xd(ib, 2),
     &              Xd(ic, 3),
     &              DC(1), DC(2), DC(3))
        Dist2 = DC(1)**2 + DC(2)**2 + DC(3)**2
        if (DistClose .gt. Dist2 .or. First) then
          DistClose = Dist2
          do i = 1, 3
            DiClose(i) = Di(i)
          end do
          First = .false.
        end if
      end do
      end do
      end do
C
      DistClose = sqrt(DistClose)
C
      return
      end
C}
C=======================================================================
C{
      subroutine GetDensCloseToSites(Dn, Dn1, Dn2, Dens,
     &                               nSites, Sites, FCTrMX)
      implicit none
C I/O
      integer           Dn(3), Dn1, Dn2
      real              Dens(0:Dn1-1, 0:Dn2-1, 0:*)
      integer           nSites
      double precision   Sites(7, *)
      double precision  FCTrMx(3, 3)
C
C local
      integer           iSite
      integer           DiClose(3)
      double precision  DistClose
C
C externals
      logical   SiteOrderF
      external  SiteOrderF, SiteExChF
C
C begin
      do iSite = 1, nSites
        call GetCloseGridPt(Dn, Sites(1, iSite), FCTrMx,
     &                          DiClose, DistClose)
        Sites(5, iSite) = Dens(DiClose(1), DiClose(2), DiClose(3))
        Sites(6, iSite) = DistClose
        call Get8PLIDens(Dn, Dn(1), Dn(2),
     &                   Dens, Sites(1, iSite), Sites(7, iSite))
      end do
C
      if (nSites .gt. 1) then
        call sort(nSites, SiteExChF, SiteOrderF,
     &            Sites, 0, 0, 0, 0, 0, 0, 0)
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine PSShowSite(Site, QFrac, FCTrMx)
      implicit none
C I/O
      double precision  Site(7)
      logical           QFrac
      double precision  FCTrMx(3, 3)
C
C local
      double precision  Scar(3)
C
C begin
      if (QFrac) then
        write(6, '(1X, A, I5, A, G14.6, 3(A, F8.5), A, 2(1X, G14.6))')
     &    'PSEARCH: ',
     &           nint(Site(4)),
     &    ' Dens = ', Site(5),
     &    ' Pos = (', Site(1),
     &           ',', Site(2),
     &           ',', Site(3), ')',
     &                Site(6),
     &                Site(7)
      else
        call TrVec3(FCTrMX, Site, Scar)
        write(6, '(1X, A, I5, A, G14.6, 3(A, F8.2), A, 2(1X, G14.6))')
     &    'PSEARCH: ',
     &           nint(Site(4)),
     &    ' Dens = ', Site(5),
     &    ' Pos = (', Scar(1),
     &           ',', Scar(2),
     &           ',', Scar(3), ')',
     &                Site(6),
     &                Site(7)
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine PL2Frac(PL, Rn, nAsyFM, Frac)
      implicit none
C I/O
      integer           PL, Rn(3), nAsyFM(2, 3)
      double precision  Frac(3)
C
C Map index stored in PeakList back to RhoMask and convert to
C fractional coordinates.
C
C local
      integer  i
      integer  Ri(3), Rj(3)
C
C begin
      call f3di(PL, Rn, Ri)
      call MapFM2RhoMask(Ri, Rn,
     &                   nAsyFM(1, 1), nAsyFM(2, 1),
     &                   nAsyFM(1, 2), nAsyFM(2, 2),
     &                   nAsyFM(1, 3), nAsyFM(2, 3),
     &                   Rj(1), Rj(2), Rj(3))
      do i = 1, 3
        Frac(i) = dfloat(Rj(i)) / Rn(i)
      end do
C
      return
      end
C}
C=======================================================================
C{
      logical function PsPlOrderF(ia, ib, PeakList, DensMap,
     &                            j1, j2, j3, j4, j5, j6)
      implicit none
C I/O
      integer  ia, ib
      integer  PeakList(*)
      real     DensMap(*)
      integer  j1, j2, j3, j4, j5, j6
C
C This functions tells the SORT procedure what is larger and
C what is smaller.
C
C begin
      PsPlOrderF = .false.
      if (DensMap(PeakList(ia)) .lt. DensMap(PeakList(ib))) return
      PsPlOrderF = .true.
      if (DensMap(PeakList(ia)) .gt. DensMap(PeakList(ib))) return
      PsPlOrderF = .false.
      if (PeakList(ia) .gt. PeakList(ib)) return
      PsPlOrderF = .true.
      return
      end
C}
C=======================================================================
C{
      subroutine PsPlExChF(ia, ib, PeakList, j1, j2, j3, j4, j5, j6, j7)
      implicit none
C I/O
      integer  ia, ib
      integer  PeakList(*)
      integer  j1, j2, j3, j4, j5, j6, j7
C
C This subroutine is called from the SORT procedure each time two list
C entries need to be swapped.
C
C local
      integer  tmp
C
C begin
      tmp = PeakList(ia)
      PeakList(ia) = PeakList(ib)
      PeakList(ib) = tmp
C
      return
      end
C}
C=======================================================================
C{
      subroutine StorePeaks(DensMap, FlagMap, Rn, nAsyFM,
     &                      PeakList, MaxPeaks, DensCutOff,
     &                      nPeaks, nList, QFrac, Symbol, SymLen,
     &                      nSites, Sites,
     &                      FCTrMx, WrnLev)
      implicit none
C I/O
      real              DensMap(*)
      integer           FlagMap(*), Rn(3), nAsyFM(2, 3)
      integer           PeakList(*), MaxPeaks
      real              DensCutOff
      integer           nPeaks, nList
      logical           QFrac
      character         Symbol*(*)
      integer           SymLen
      integer           nSites
      double precision   Sites(7, *)
      double precision  FCTrMx(3, 3)
      integer           WrnLev
C
C Collect, sort and print the peaks.
C Also declare symbols if Symbol is not an empty string.
C
C local
      logical           ContLoop
      integer           PRn, iMap, iPeak, iSite, i
      character         xyz*3
      double precision  Dens, MaxF(3), MaxC(3)
      double complex    dcval
C
      integer   StrMax, StrLen
      parameter(StrMax = 132)
      character*(StrMax)  Str
      character*(StrMax)  buf
C
C externals
      logical   PsPlOrderF
      external  PsPlOrderF, PsPlExChF
C
C begin
      xyz = 'XYZ'
C
      PRn = Rn(1) * Rn(2) * Rn(3)
C
      nPeaks = 0
C
      do iMap = 1, PRn
        if (      FlagMap(iMap) .lt. 0
     &      .and. DensMap(iMap) .ge. DensCutOff) then
          if (nPeaks .ge. MaxPeaks) then
            nPeaks = -1
            return
          end if
                   nPeaks = nPeaks + 1
          PeakList(nPeaks) = iMap
        end if
      end do
C
      if (nPeaks .gt. 1) then
        call sort(nPeaks, PsPlExChF, PsPlOrderF,
     &            PeakList, DensMap, 0, 0, 0, 0, 0, 0)
      end if
C
      if (WrnLev .ge. 10)
     &  write(6, '(1X, A, I8)')
     &    'PSEARCH: Number of peaks in overstore buffer = ',
     &    nPeaks
C
      if (WrnLev .ge. 5) then
        write(6, '(1X, A, I8)')
     &    'PSEARCH: Number of peaks listed              = ',
     &    min(nList, nPeaks)
        if (.not. QFrac) then
          write(6, '(1X, 2A)')
     &      'PSEARCH:                                   ',
     &      '(  Cartesian Coordinates   )'
        else
          write(6, '(1X, 2A)')
     &      'PSEARCH:                                   ',
     &      '(  Fractional Coordinates  )'
        end if
      end if
C
      if (SymLen .gt. 0) then
        buf = Symbol(1:SymLen) // '_NLIST'
        Dens = min(nList, nPeaks)
        call declar(buf, 'DP', ' ', dcval, Dens)
      end if
C
      iSite = 1
C
      do iPeak = 1, min(nList, nPeaks)
        call PL2Frac(PeakList(iPeak), Rn, nAsyFM, MaxF)
        Dens = DensMap(PeakList(iPeak))
C
        call TrVec3(FCTrMX, MaxF, MaxC)
C
        if (WrnLev .ge. 5) then
                    ContLoop = .true.
          do while (ContLoop .and. iSite .le. nSites)
            if (Sites(5, iSite) .le. Dens) then
              ContLoop = .false.
            else
              call PSShowSite(Sites(1, iSite), QFrac, FCTrMx)
              iSite = iSite + 1
            end if
          end do
C
          if (.not. QFrac) then
            write(6, '(1X, A, I5, A, G14.6, 3(A, F8.2), A)')
     &        'PSEARCH: ', iPeak,
     &        ' Dens = ', Dens,
     &        ' Pos = (', MaxC(1), ',', MaxC(2), ',', MaxC(3), ')'
          else
            write(6, '(1X, A, I5, A, G14.6, 3(A, F8.5), A)')
     &        'PSEARCH: ', iPeak,
     &        ' Dens = ', Dens,
     &        ' Pos = (', MaxF(1), ',', MaxF(2), ',', MaxF(3), ')'
          end if
        end if
C
        if (SymLen .gt. 0) then
          call encodi(iPeak, Str, StrMax, StrLen)
C
          do i = 1, 3
            buf =    Symbol(1:SymLen)
     &            // '_'
     &            // xyz(i:i)
     &            // '_'
     &            // Str(1:StrLen)
            call declar(buf, 'DP', ' ', dcval, MaxC(i))
          end do
C
          buf =    Symbol(1:SymLen)
     &          // '_DENS_'
     &          // Str(1:StrLen)
          call declar(buf, 'DP', ' ', dcval, Dens)
        end if
      end do
C
      if (nPeaks .gt. 0) then
        call PL2Frac(PeakList(1), Rn, nAsyFM, MaxF)
        Dens = DensMap(PeakList(1))
      else
        do i = 1, 3
          MaxF(i) = 0.
        end do
        Dens = 0.
      end if
C
      call TrVec3(FCTrMX, MaxF, MaxC)
C
      call declar('RESULT',  'DP', ' ', dcval, Dens)
      call declar('PSMAX',   'DP', ' ', dcval, Dens)
      call declar('PSMAX_X', 'DP', ' ', dcval, MaxC(1))
      call declar('PSMAX_Y', 'DP', ' ', dcval, MaxC(2))
      call declar('PSMAX_Z', 'DP', ' ', dcval, MaxC(3))
C
      if (SymLen .gt. 0 .and. nList .gt. 0 .and. nPeaks .eq. 0)
     &  write(6, '(1X, A)')
     &    'PSEARCH-WRN: No peaks found => NO SYMBOLS DECLARED.'
C
      if (WrnLev .ge. 5) then
        do while (iSite .le. nSites)
          call PSShowSite(Sites(1, iSite), QFrac, FCTrMx)
          iSite = iSite + 1
        end do
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine QueryPStore(nPeaks, PeakList, iPeak,
     &                       Rn, nAsyFM,
     &                       QFrac, Symbol, SymLen,
     &                       FCTrMx, WrnLev)
      implicit none
C I/O
      integer           nPeaks
      integer           PeakList(*)
      integer           iPeak
      integer           Rn(3), nAsyFM(2, 3)
      logical           QFrac
      character         Symbol*(*)
      integer           SymLen
      double precision  FCTrMx(3, 3)
      integer           WrnLev
C
C Copy coordinates of peak number iPeak to
C symbols if Symbol is not an empty string.
C
C local
      integer           i
      character         xyz*3
      double precision  MaxF(3), MaxC(3)
      double precision  dpval
      double complex    dcval
C
      integer   StrMax
      parameter(StrMax = 132)
      character*(StrMax)  buf
C
C begin
      dpval = dfloat(0)
      dcval = dcmplx(dpval, dpval)
C
      xyz = 'XYZ'
C
      if (iPeak .gt. 0 .and. iPeak .le. nPeaks) then
        call PL2Frac(PeakList(iPeak), Rn, nAsyFM, MaxF)
        call TrVec3(FCTrMX, MaxF, MaxC)
C
        if (WrnLev .ge. 5) then
          if (.not. QFrac) then
            write(6, '(1X, A, I5, 3(A, F8.2), A)')
     &        'PSEARCH: ', iPeak,
     &        ' Pos = (', MaxC(1), ',', MaxC(2), ',', MaxC(3), ')'
          else
            write(6, '(1X, A, I5, 3(A, F8.5), A)')
     &        'PSEARCH: ', iPeak,
     &        ' Pos = (', MaxF(1), ',', MaxF(2), ',', MaxF(3), ')'
          end if
        end if
C
        if (SymLen .gt. 0) then
          do i = 1, 3
            buf =    Symbol(1:SymLen)
     &            // '_'
     &            // xyz(i:i)
            call declar(buf, 'DP', ' ', dcval, MaxC(i))
          end do
        end if
C
        call declar('RESULT', 'LO', 'TRUE', dcval, dpval)
      else
        if (WrnLev .ge. 5) then
            buf = '.'
          if (SymLen .gt. 0)
     &      buf = ' => NO SYMBOLS DECLARED.'
          write(6, '(1X, A, I5, 2A)')
     &      'PSEARCH-WRN: peak no.', iPeak,
     &      ' does not exist', buf
        end if
C
        call declar('RESULT', 'LO', 'FALSE', dcval, dpval)
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine IniPSTCommon(ColdStart)
      implicit none
C I/O
      INCLUDE 'funct.inc'
      INCLUDE 'psearch.inc'
      logical  ColdStart
C
C Initialize psearch.inc common blocks.
C
C local
      integer  iPStore, i, j
C
C begin
      do iPStore = 1, MaxPStores
        if (.not. ColdStart ) then
C mrt for the sake of undefined vars checking (lf95 --chk eux)
          if( npPStore(iPStore) .gt. 0) then
            call FreHp(ptPStore(iPStore), integ4(npPStore(iPStore)))
          end if
        end if
        do i = 1, 3
          rnPStore(i, iPStore) = 0
          do j = 1, 2
            naPStore(j, i, iPStore) = 0
          end do
        end do
        npPStore(iPStore) = 0
        ptPStore(iPStore) = 0
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine psrchp(iFrom, Level, nList, QFrac, Symbol, SymLen,
     &                  iPStore, MaxPStores,
     &                  xRhoNum, xRhoNam,
     &                  mSites, nSites, Sites, FCTrMx, CFTrMx)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      integer           iFrom, Level, nList
      logical           QFrac
      character         Symbol*(*)
      integer           SymLen
      integer           iPStore, MaxPStores
      integer           xRhoNum
      character*(*)     xRhoNam(*)
      integer           mSites, nSites
      double precision   Sites(7, *)
      double precision  FCTrMx(3, 3), CFTrMx(3, 3)
C
C local
      integer           iSite, i
      double precision  S(3)
C
C This is the parser for psearch
C
C begin
      iFrom = 0
      Level = 1
      nList = 101
      QFrac = .false.
      Symbol = ' '
      iPStore = 0
      SymLen = 0
      nSites = 0
C
C parsing
      CALL PUSEND('PSEARCH>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('PSEARCH>')
      CALL MISCOM('PSEARCH>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-psearch')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FROM') THEN
        call GetiObj(iFrom, xRhoNum, xRhoNam, ' ',
     &               'FROM=', 'real', 'PSEARCH')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'LEVE') THEN
        CALL NEXTI('LEVEl=', i)
        if (i .lt. 0 .or. i .gt. 3) then
          CALL DSPERR('PSEARCH', 'parameter out of range')
        else
          Level = i
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NLIS') THEN
        CALL NEXTI('NLISt=', i)
        if (i .lt. 0) then
          CALL DSPERR('PSEARCH', 'parameter out of range')
        else
          nList = i
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FRAC') THEN
        CALL NEXTLO('FRACtional=', QFrac)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SYMB') THEN
        CALL NEXTST('Symbol=', WD)
        SymLen = max(0, min(len(Symbol), WDLEN))
        Symbol = ' '
        if (SymLen .gt. 0) then
          Symbol = WD(1:WDLEN)
C If Symbol was placed in quotes, it might have lower case characters.
          call ToUpper(Symbol(:SymLen))
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PSTO') THEN
        CALL NEXTI('PSTOre=', i)
        if (i .lt. 0 .or. i .gt. MaxPStores) then
          CALL DSPERR('PSEARCH', 'parameter out of range')
        else
          iPStore = i
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SITE') THEN
        CALL NEXTF('SITE(X)=', S(1))
        CALL NEXTF('SITE(Y)=', S(2))
        CALL NEXTF('SITE(Z)=', S(3))
        if (nSites .ge. mSites) then
          call WrnDie(0, 'PSEARCH', 'Too many SITE definitions.')
        else
          nSites = nSites + 1
C
          if (QFrac) then
            do i = 1, 3
              Sites(i, nSites) = S(i)
            end do
          else
            call TrVec3(CFTrMX, S, Sites(1, nSites))
          end if
C
          Sites(4, nSites) = dfloat(-nSites)
          Sites(5, nSites) = dfloat(-1)
          Sites(6, nSites) = dfloat(-1)
          Sites(7, nSites) = dfloat(-1)
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'RESE') THEN
        call IniPSTCommon(.false.)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
        WRITE(6,'(2A)') ' -----------psearch-parameters',
     &   '--------------------------------------------------'
        call EchoObj(iFrom, xRhoNam, 'FROM=')
        write(6, '(1X, A, I3)') 'LEVEl=', Level
        write(6, '(1X, A, I6)') 'NLISt=', nList
        if (QFrac) then
          write(6, '(1X, A)') 'FRACtional=TRUE'
        else
          write(6, '(1X, A)') 'FRACtional=FALSE'
        end if
        if (SymLen .lt. 1) then
          write(6, '(1X, A)') 'Symbol=<undefined>'
        else
          write(6, '(1X, 2A)') 'SYMBol=', Symbol(1:SymLen)
        end if
        write(6, '(1X, A, I6)') 'PSTOre=', iPStore
        do iSite = 1, nSites
          if (QFrac) then
            write(6, '(1X, A, 3(1X, F8.5))')
     &        'SITE', (Sites(i, iSite), i = 1, 3)
          else
            call TrVec3(FCTrMX, Sites(1, iSite), S)
            write(6, '(1X, A, 3(1X, F8.2))')
     &        'SITE', (S(i), i = 1, 3)
          end if
        end do
        WRITE(6,'(2A)') ' -----------------------------',
     &   '--------------------------------------------------'
C=====================================================================
      ELSE
        CALL CHKEND('PSEARCH>',DONE)
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
      subroutine psearch(na, nb, nc, nrRho,
     &                   xRhoNum, xRhoNam, hprrho,
     &                   xrNsym, FCTrMx, CFTrMx)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'flagmap.inc'
      INCLUDE 'psearch.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      integer           na, nb, nc, nrRho
      integer           xRhoNum
      character*(*)     xRhoNam(*)
      integer           hprrho(*)
      integer           xrNsym
      double precision  FCTrMx(3, 3), CFTrMx(3, 3)
C
C local
      integer    iFrom, Level, nList, i, j
      logical    QFrac
      character  Symbol*64
      integer    SymLen
      integer    iPStore
C
      integer  PRn, MaxPeaks, TotPeaks, nPeaks
      real     DensCutOff
C
      integer           mSites, nSites
      parameter        (mSites = 100)
      double precision   Sites(7, mSites)
C
      double precision  dpval
      double complex    dcval
      double precision  tmEntr, tmExit
C
C pointer
      integer  PeakList
C
C externals
      integer   lnblnk
      external  lnblnk
C
C begin
      dpval = dfloat(0)
      dcval = dcmplx(dpval, dpval)
C
      if (Timer .gt. 0) call VCPU(tmEntr)
C
C Call the parser
      call  psrchp(iFrom, Level, nList, QFrac, Symbol, SymLen,
     &             iPStore, MaxPStores,
     &             xRhoNum, xRhoNam,
     &             mSites, nSites, Sites, FCTrMx, CFTrMx)
C
      if (iFrom .gt. 0) then
        if (hprrho(iFrom) .eq. 0) then
          write(6, '(3A)') ' %PSEARCH-ERR: Real space object "',
     &      xrhonam(iFrom)(1:lnblnk(xrhonam(iFrom))), '" not defined.'
          call WrnDie(0, 'PSEARCH', 'object undefined.')
          return
        end if
C
        if (xrNsym .ne. 1) then
          call WrnDie(0, 'PSEARCH',
     &      'This module can only be used in P1.')
          return
        end if
C
        if (     .not. ValidFlagMap
     &      .or. RnFM(1) .ne. na
     &      .or. RnFM(2) .ne. nb
     &      .or. RnFM(3) .ne. nc) then
          call WrnDie(0, 'PSEARCH',
     &      'FlagMap undefined or out of date: call FMAP first.')
          return
        end if
C
        PRn = RnFM(1) * RnFM(2) * RnFM(3)
C
        if (nrRho .ne. PRn) then
          call WrnDie(-5, 'PSEARCH', 'Fatal Coding Error.')
          call Die
        end if
C
C Allow for some extra peaks
        if (iPStore .eq. 0) then
          MaxPeaks = nList + 100
        else
          MaxPeaks = PRn
        end if
C
C Peak search in P1
        call PeakSearch(heap(hprrho(iFrom)), heap(ptFlagMap), RnFM,
     &                  Level, MaxPeaks, DensCutOff, TotPeaks, WrnLev)
C
C Allocate space for peak list
        if (iPStore .eq. 0) then
          MaxPeaks = min(MaxPeaks, TotPeaks)
        else
          MaxPeaks = TotPeaks
        end if
          PeakList = 0
        if (MaxPeaks .gt. 0)
     &    PeakList = AllHp(integ4(MaxPeaks))
C
C Get densities closest to SITEs
        call GetDensCloseToSites(RnFM, RnFM(1), RnFM(2),
     &                           heap(hprrho(iFrom)),
     &                           nSites, Sites, FCTrMx)
C
C Collect, sort and print the peaks
        call StorePeaks(heap(hprrho(iFrom)), heap(ptFlagMap), RnFM,
     &                  nAsyFM,
     &                  heap(PeakList), MaxPeaks, DensCutOff,
     &                  nPeaks, nList, QFrac, Symbol, SymLen,
     &                  nSites, Sites,
     &                  FCTrMx, WrnLev)
C
        if (nPeaks .lt. 0) then
          call WrnDie(-5, 'PSEARCH', 'Fatal Coding Error.')
          call Die
        end if
C
        if (iPStore .eq. 0) then
          if (MaxPeaks .gt. 0)
     &      call FreHp(PeakList, integ4(MaxPeaks))
        else
          if (npPStore(iPStore) .gt. 0)
     &      call FreHp(ptPStore(iPStore), integ4(npPStore(iPStore)))
          do i = 1, 3
            rnPStore(i, iPStore) = RnFM(i)
            do j = 1, 2
              naPStore(j, i, iPStore) = nAsyFM(j, i)
            end do
          end do
          npPStore(iPStore) = MaxPeaks
          ptPStore(iPStore) = PeakList
          if (WrnLev .ge. 5) then
            write(6, '(1X, A, I5, A, I2, A)')
     &        'PSEARCH: ', MaxPeaks,
     &        ' peak(s) stored in PSTORE ', iPStore, '.'
          end if
          dpval = dfloat(MaxPeaks)
          call declar('PSTORE', 'DP', ' ', dcval, dpval)
        end if
      else if (iPStore .gt. 0) then
        call QueryPStore(npPStore(iPStore), heap(ptPStore(iPStore)),
     &                   nList,
     &                   rnPStore(1, iPStore), naPStore(1, 1, iPStore),
     &                   QFrac, Symbol, SymLen,
     &                   FCTrMx, WrnLev)
      else
        call WrnDie(0, 'PSEARCH', 'Real space object or PSTORE needed.')
        return
      end if
C
      if (Timer .gt. 0) then
        call VCPU(tmExit)
        write(6, '(1X, A, F10.4, A)')
     &    'PSEARCH: CPU-time: ', tmExit - tmEntr, ' seconds total'
      end if
C
      return
      end
C}
