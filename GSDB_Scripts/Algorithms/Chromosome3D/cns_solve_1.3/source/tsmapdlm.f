C=======================================================================
C{
      subroutine tslmSw3(Wn, Wn1, Wn2, RowInfo,
     &                   fcg, nfcg,
     &                   xrh, xrk, xrl,
     &                   Mult, fpart, nRef, QHerm, QFpart,
     &                   SymMx, mSymMx, nSymMx, STBF,
     &                   hMs, kMs, lMs, fts)
      IMPLICIT NONE
C I/O
      integer           Wn(3), Wn1, Wn2
      integer           RowInfo(4, 0:Wn2-1, 0:(Wn1/2)-1)
      double complex    fcg(*)
      integer           nfcg(3)
      integer           xrh(*), xrk(*), xrl(*), Mult(*)
      double complex    fpart(*)
      integer           nRef
      logical           QHerm, QFpart
      integer           mSymMx, SymMx(mSymMx, 3, 4), nSymMx, STBF
      integer           hMs(*), kMs(*), lMs(*)
      double complex    fts(*)
C
C Compute sum for Fast Translation Search
C   ws3 = Eq. (15) of Navaza paper, p. 447
C
C local
      integer         iRef, mH, Hts
      integer         iS, iS0, iS1, iS2, iS3
      integer         hM0, kM0, lM0
      integer         hM01, kM01, lM01
      integer         hM0p1, kM0p1, lM0p1
      integer         hM01p2, kM01p2, lM01p2
      integer         hM01p23, kM01p23, lM01p23
      logical         vWi, vWj
      integer         Wi(3), Wj(3)
      double complex  ftil0c, ftil1, ftil2c
      double complex  cf
C
C externals
      logical   hkl2Wi
      external  hkl2Wi
C
C begin
      do iRef = 1, nRef
        mH = Mult(iRef)
        if (mH .gt. 0) then
C Tabulate sym. equiv. hkl for given reflection
          do iS = 1, nSymMx
            call ftilSet(xrh(iRef), xrk(iRef), xrl(iRef),
     &                   mSymMx, SymMx, iS,
     &                   nfcg(1), nfcg(2), nfcg(3), fcg, QHerm,
     &                   hMs(iS), kMs(iS), lMs(iS), Hts, STBF,
     &                   fts(iS), .false.)
          end do
C
C Enter 4-deep summation loop
          do iS0 = 1, nSymMx
            hM0 = hMs(iS0)
            kM0 = kMs(iS0)
            lM0 = lMs(iS0)
            ftil0c = dconjg(fts(iS0))
C
            if (QFpart) then
              vWi = hkl2Wi(-hM0, -kM0, -lM0, Wn, Wi)
              vWj = hkl2Wi( hM0,  kM0,  lM0, Wn, Wj)
              cf =   (2 * mH) * ftil0c
     &             * (fpart(iRef)**2)
     &             * dconjg(fpart(iRef))
              if (vWi)
     &          call Add2Row(RowInfo(1,Wi(2),Wi(1)),Wi(3),cf)
              if (vWj)
     &          call Add2Row(RowInfo(1,Wj(2),Wj(1)),Wj(3),dconjg(cf))
            end if
C
            do iS1 = 1, nSymMx
              hM01 = hM0 - hMs(iS1)
              kM01 = kM0 - kMs(iS1)
              lM01 = lM0 - lMs(iS1)
              ftil1 = fts(iS1)
C
              if (QFpart) then
                vWi = hkl2Wi(-hM01, -kM01, -lM01, Wn, Wi)
                if (vWi) then
                  cf =   (4 * mH) * ftil0c * ftil1
     &                 * dconjg(fpart(iRef)) * fpart(iRef)
                  call Add2Row(RowInfo(1,Wi(2),Wi(1)),Wi(3),cf)
                end if
C
                hM0p1 = hM0 + hMs(iS1)
                kM0p1 = kM0 + kMs(iS1)
                lM0p1 = lM0 + lMs(iS1)
                vWi = hkl2Wi(-hM0p1, -kM0p1, -lM0p1, Wn, Wi)
                vWj = hkl2Wi( hM0p1,  kM0p1,  lM0p1, Wn, Wj)
                cf = mH * ftil0c * dconjg(ftil1) * (fpart(iRef)**2)
                if (vWi)
     &            call Add2Row(RowInfo(1,Wi(2),Wi(1)),Wi(3),cf)
                if (vWj)
     &            call Add2Row(RowInfo(1,Wj(2),Wj(1)),Wj(3),dconjg(cf))
              end if
C
              do iS2 = 1, nSymMx
                hM01p2 = hM01 + hMs(iS2)
                kM01p2 = kM01 + kMs(iS2)
                lM01p2 = lM01 + lMs(iS2)
                ftil2c = dconjg(fts(iS2))
C
                if (QFpart) then
                  vWi = hkl2Wi(-hM01p2, -kM01p2, -lM01p2, Wn, Wi)
                  vWj = hkl2Wi( hM01p2,  kM01p2,  lM01p2, Wn, Wj)
                  cf = (2 * mH) * ftil0c * ftil1 * ftil2c * fpart(iRef)
                  if (vWi)
     &             call Add2Row(RowInfo(1,Wi(2),Wi(1)),Wi(3),cf)
                  if (vWj)
     &             call Add2Row(RowInfo(1,Wj(2),Wj(1)),Wj(3),dconjg(cf))
                end if
C
                do iS3 = 1, nSymMx
                  hM01p23 = hM01p2 - hMs(iS3)
                  kM01p23 = kM01p2 - kMs(iS3)
                  lM01p23 = lM01p2 - lMs(iS3)
C
                  vWi = hkl2Wi(-hM01p23,-kM01p23,-lM01p23, Wn, Wi)
                  if (vWi) then
                    cf = mH * ftil0c * ftil1 * ftil2c * fts(iS3)
                    call Add2Row(RowInfo(1,Wi(2),Wi(1)),Wi(3),cf)
                  end if
                end do
              end do
            end do
          end do
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine tslmSw1(Wn, Wn1, Wn2, RowInfo,
     &                   fcg, nfcg,
     &                   xrh, xrk, xrl,
     &                   Mult, fpart, nRef, QHerm, QFpart,
     &                   SymMx, mSymMx, nSymMx, STBF,
     &                   hMs, kMs, lMs, fts)
      IMPLICIT NONE
C I/O
      integer           Wn(3), Wn1, Wn2
      integer           RowInfo(4, 0:Wn2-1, 0:(Wn1/2)-1)
      double complex    fcg(*)
      integer           nfcg(3)
      integer           xrh(*), xrk(*), xrl(*), Mult(*)
      double complex    fpart(*)
      integer           nRef
      logical           QHerm, QFpart
      integer           mSymMx, SymMx(mSymMx, 3, 4), nSymMx, STBF
      integer           hMs(*), kMs(*), lMs(*)
      double complex    fts(*)
C
C Compute sum for Fast Translation Search
C   ws1 = "similar expression" (Navaza paper, p. 447, right after Eq. (14))
C
C local
      integer         iRef, mH, Hts
      integer         iS, iS0, iS1
      integer         hM0, kM0, lM0
      integer         hM01, kM01, lM01
      logical         vWi, vWj
      integer         Wi(3), Wj(3)
      double complex  ftil0c
      double complex  cf
C
C externals
      logical   hkl2Wi
      external  hkl2Wi
C
C begin
      do iRef = 1, nRef
        mH = Mult(iRef)
        if (mH .gt. 0) then
C Tabulate sym. equiv. hkl for given reflection
          do iS = 1, nSymMx
            call ftilSet(xrh(iRef), xrk(iRef), xrl(iRef),
     &                   mSymMx, SymMx, iS,
     &                   nfcg(1), nfcg(2), nfcg(3), fcg, QHerm,
     &                   hMs(iS), kMs(iS), lMs(iS), Hts, STBF,
     &                   fts(iS), .false.)
          end do
C
C Enter 2-deep summation loop
          do iS0 = 1, nSymMx
            hM0 = hMs(iS0)
            kM0 = kMs(iS0)
            lM0 = lMs(iS0)
            ftil0c = dconjg(fts(iS0))
C
            if (QFpart) then
              vWi = hkl2Wi(-hM0, -kM0, -lM0, Wn, Wi)
              vWj = hkl2Wi( hM0,  kM0,  lM0, Wn, Wj)
              cf = mH * ftil0c * fpart(iRef)
              if (vWi)
     &          call Add2Row(RowInfo(1,Wi(2),Wi(1)),Wi(3),cf)
              if (vWj)
     &          call Add2Row(RowInfo(1,Wj(2),Wj(1)),Wj(3),dconjg(cf))
            end if
C
            do iS1 = 1, nSymMx
              hM01 = hM0 - hMs(iS1)
              kM01 = kM0 - kMs(iS1)
              lM01 = lM0 - lMs(iS1)
C
              vWi = hkl2Wi(-hM01, -kM01, -lM01, Wn, Wi)
              if (vWi) then
                cf = mH * ftil0c * fts(iS1)
                call Add2Row(RowInfo(1,Wi(2),Wi(1)),Wi(3),cf)
              end if
            end do
          end do
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine tslmSw2(Wn, Wn1, Wn2, RowInfo,
     &                   fcg, nfcg,
     &                   xrh, xrk, xrl,
     &                   Mult, dIobs, fpart, nRef, QHerm,
     &                   QFpart,
     &                   SymMx, mSymMx, nSymMx, STBF,
     &                   hMs, kMs, lMs, fts)
      IMPLICIT NONE
C I/O
      integer           Wn(3), Wn1, Wn2
      integer           RowInfo(4, 0:Wn2-1, 0:(Wn1/2)-1)
      double complex    fcg(*)
      integer           nfcg(3)
      integer           xrh(*), xrk(*), xrl(*), Mult(*)
      double precision  dIobs(*)
      double complex    fpart(*)
      integer           nRef
      logical           QHerm,  QFpart
      integer           mSymMx, SymMx(mSymMx, 3, 4), nSymMx, STBF
      integer           hMs(*), kMs(*), lMs(*)
      double complex    fts(*)
C
C Compute sum for Fast Translation Search
C   ws2 = Eq. (14) of Navaza paper, p. 447
C
C local
      integer         iRef, mH, Hts
      integer         iS, iS0, iS1
      integer         hM0, kM0, lM0
      integer         hM01, kM01, lM01
      logical         vWi, vWj
      integer         Wi(3), Wj(3)
      double complex  ftil0c
      double complex  cf
C
C externals
      logical   hkl2Wi
      external  hkl2Wi
C
C begin
      do iRef = 1, nRef
        mH = Mult(iRef)
        if (mH .gt. 0) then
C Tabulate sym. equiv. hkl for given reflection
          do iS = 1, nSymMx
            call ftilSet(xrh(iRef), xrk(iRef), xrl(iRef),
     &                   mSymMx, SymMx, iS,
     &                   nfcg(1), nfcg(2), nfcg(3), fcg, QHerm,
     &                   hMs(iS), kMs(iS), lMs(iS), Hts, STBF,
     &                   fts(iS), .false.)
          end do
C
C Enter 2-deep summation loop
          do iS0 = 1, nSymMx
            hM0 = hMs(iS0)
            kM0 = kMs(iS0)
            lM0 = lMs(iS0)
            ftil0c = dconjg(fts(iS0))
C
            if (QFpart) then
              vWi = hkl2Wi(-hM0, -kM0, -lM0, Wn, Wi)
              vWj = hkl2Wi( hM0,  kM0,  lM0, Wn, Wj)
              cf = mH * ftil0c * fpart(iRef) * dIobs(iRef)
              if (vWi)
     &          call Add2Row(RowInfo(1,Wi(2),Wi(1)),Wi(3),cf)
              if (vWj)
     &          call Add2Row(RowInfo(1,Wj(2),Wj(1)),Wj(3),dconjg(cf))
            end if
C
            do iS1 = 1, nSymMx
              hM01 = hM0 - hMs(iS1)
              kM01 = kM0 - kMs(iS1)
              lM01 = lM0 - lMs(iS1)
C
              vWi = hkl2Wi(-hM01, -kM01, -lM01, Wn, Wi)
              if (vWi) then
                cf = mH * ftil0c * fts(iS1) * dIobs(iRef)
                call Add2Row(RowInfo(1,Wi(2),Wi(1)),Wi(3),cf)
              end if
            end do
          end do
        end if
      end do
C
      return
      end
C}
