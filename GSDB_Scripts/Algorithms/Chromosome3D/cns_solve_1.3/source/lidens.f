C=======================================================================
C{
      subroutine LICorner(Xfrac, Gn, Gi, Xw)
      implicit none
C I/O
      double precision  Xfrac(3)
      integer           Gn(3), Gi(3)
      double precision  Xw(3, 0:1)
C
C local
      integer           i
      double precision  Xg
C
C begin
      do i = 1, 3
        Xg = mod(Xfrac(i), dfloat(1))
        Xg = mod(Xg + dfloat(2), dfloat(1))
        Xg = Xg * dfloat(Gn(i))
        Gi(i) = int(Xg)
        Xw(i, 1) = Xg - dfloat(Gi(i))
        Xw(i, 0) = dfloat(1) - Xw(i, 1)
        if (Gi(i) .lt. 0 .or. Gi(i) .ge. Gn(i)) then
          call WrnDie(-5, 'LICorner', 'Fatal Coding Error.')
          call Die
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine Get8PLIDens(Dn, Dn1, Dn2, Dens, X, LIDens)
      implicit none
C I/O
      integer           Dn(3), Dn1, Dn2
      real              Dens(0:Dn1-1, 0:Dn2-1, 0:*)
      double precision  X(3), LIDens
C
C 8-point linear interpolation to determine density at position X
C
C local
      integer           Di(3), Dg(3), ia, ib, ic
      double precision  Xw(3, 0:1)
C
C begin
      call LICorner(X, Dn, Di, Xw)
C
      LIDens = dfloat(0)
C
      do ic = 0, 1
        Dg(3) = mod(Di(3) + ic, Dn(3))
      do ib = 0, 1
        Dg(2) = mod(Di(2) + ib, Dn(2))
      do ia = 0, 1
        Dg(1) = mod(Di(1) + ia, Dn(1))
        LIDens = LIDens +   Dens(Dg(1), Dg(2), Dg(3))
     &                    * Xw(1, ia) * Xw(2, ib) * Xw(3, ic)
      end do
      end do
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine PrintLIDens(Sfrac, result, QFrac, FCTrMX)
      implicit none
C I/O
      double precision  Sfrac(3), result
      logical           QFrac
      double precision  FCTrMX(3, 3)
C
C local
      integer           i
      double precision  Sorth(3)
C
C begin
      if (QFrac) then
        write(6, '(1X, A, 3(1X, F8.5, A, G12.6))')
     &    'LIDens: Density at site',
     &    (Sfrac(i), i = 1, 3), ' = ', result
      else
        call TrVec3(FCTrMX, Sfrac, Sorth)
        write(6, '(1X, A, 3(1X, F8.2), A, G12.6)')
     &    'LIDens: Density at site',
     &    (Sorth(i), i = 1, 3), ' = ', result
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine LIDensPa(iDens,
     &                    xRhoNum, xRhoNam,
     &                    nAtSele, AtSele, CoordX, CoordY, CoordZ,
     &                    To,
     &                    QFrac, QSite, Site,
     &                    FCTrMx, CFTrMx)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      integer           iDens
      integer           xRhoNum
      character         xRhoNam(*)*(*)
      integer           nAtSele, AtSele(*)
      double precision  CoordX(*), CoordY(*), CoordZ(*)
      character         To*(*)
      logical           QFrac, QSite
      double precision  Site(3)
      double precision  FCTrMx(3, 3), CFTrMx(3, 3)
C
C This is the parser for LIDens
C
C local
      integer           i
      double precision  Sorth(3)
C
C begin
      iDens = 0
      nAtSele = 0
      To = ' '
      QFrac = .false.
      QSite = .false.
      do i = 1, 3
        Site(i) = dfloat(0)
      end do
C
C parsing
      CALL PUSEND('LIDens>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('LIDens>')
      CALL MISCOM('LIDens>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-lidens')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DENS') THEN
        call GetiObj(iDens, xRhoNum, xRhoNam, ' ',
     &               'density=', 'real', 'LIDens')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ATOM') THEN
        call SelctA(AtSele, nAtSele, CoordX, CoordY, CoordZ, .true.)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TO  ') THEN
        CALL NEXTASS('to=')
        To = ' '
        To = WD(1:WDLEN)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FRAC') THEN
        CALL NEXTLO('fractional=', QFrac)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SITE') THEN
        CALL NEXTF('site(x)=', Site(1))
        CALL NEXTF('site(y)=', Site(2))
        CALL NEXTF('site(z)=', Site(3))
        if (.not. QFrac) then
          do i = 1, 3
            Sorth(i) = Site(i)
          end do
          call TrVec3(CFTrMX, Sorth, Site)
        end if
        QSite = .true.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
        WRITE(6,'(2A)') ' -----------lidens-parameters-',
     &   '--------------------------------------------------'
        call EchoObj(iDens,  xRhoNam, 'density=')
        write(6, '(1X, A, I6)') 'Number of atoms selected: ', nAtSele
        call EchoObj(1,  To, 'to=')
        if (QFrac) then
          write(6, '(1X, A)') 'fractional=true'
        else
          write(6, '(1X, A)') 'fractional=false'
        end if
        if (QSite) then
          if (QFrac) then
            write(6, '(1X, A, 3(1X, F8.5))')
     &        'LIDens: Density at site',
     &        ( Site(i), i = 1, 3)
          else
            call TrVec3(FCTrMX, Site, Sorth)
            write(6, '(1X, A, 3(1X, F8.2))')
     &        'LIDens: Density at site',
     &        (Sorth(i), i = 1, 3)
          end if
        end if
        WRITE(6,'(2A)') ' -----------------------------',
     &   '--------------------------------------------------'
C=====================================================================
      ELSE
        CALL CHKEND('LIDens>',DONE)
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
      subroutine LIDens(na, nb, nc, nrRho,
     &                  xRhoNum, xRhoNam, hprrho,
     &                  xrNsym, FCTrMx, CFTrMx,
     &                  nAtom, AtSele, CoordX, CoordY, CoordZ)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'numbers.inc'
      integer           na, nb, nc, nrRho
      integer           xRhoNum
      character*(*)     xRhoNam(*)
      integer           hprrho(*)
      integer           xrNsym
      double precision  FCTrMx(3, 3), CFTrMx(3, 3)
      integer           nAtom, AtSele(*)
      double precision  CoordX(*), CoordY(*), CoordZ(*)
C
C XRAY LIDens main procedure.
C
C local
      integer           iDens, Dn(3), nAtSele, iAtom
      character         To*(WORD_SIZE)
      logical           QFrac, QSite
      double precision  Site(3), Sfrac(3), result
      double complex    dbcomp
C
C pointers
      integer  ptResults
C
C externals
      integer   lnblnk
      external  lnblnk
C
C begin
C Call the parser
      call LIDensPa(iDens,
     &              xRhoNum, xRhoNam,
     &              nAtSele, AtSele, CoordX, CoordY, CoordZ,
     &              To,
     &              QFrac, QSite, Site,
     &              FCTrMx, CFTrMx)
C
      if (iDens .lt. 1) then
        call WrnDie(0, 'LIDens', 'Density object not specified.')
        return
      end if
C
      if (hprrho(iDens) .eq. 0) then
        write(6, '(3A)') ' %LIDens-ERR: Real space object "',
     &    xrhonam(iDens)(1:lnblnk(xrhonam(iDens))), '" not defined.'
        call WrnDie(0, 'LIDens', 'object undefined.')
        return
      end if
C
      if (nAtSele .gt. 0 .and. To .eq. ' ') then
        call WrnDie(-5, 'LIDens', 'atom object missing.')
        return
      end if
C
      if (xrNsym .ne. 1) then
        call WrnDie(0, 'LIDens',
     &    'This module can only be used in P1.')
        return
      end if
C
      if (nrRho .ne. na * nb * nc) then
        call WrnDie(-5, 'LIDens', 'Fatal Coding Error.')
        call Die
      end if
C
      Dn(1) = na
      Dn(2) = nb
      Dn(3) = nc
C
      if (nAtSele .gt. 0) then
        ptResults = AllHp(ireal8(nAtom))
        call FillR8(heap(ptResults), nAtom, Zero)
        do iAtom = 1, nAtom
          if (AtSele(iAtom) .eq. 1) then
            call TrCoor(CFTrMX, CoordX(iAtom),
     &                          CoordY(iAtom),
     &                          CoordZ(iAtom),
     &                          Sfrac(1), Sfrac(2), Sfrac(3))
            call Get8PLIDens(Dn, Dn(1), Dn(2), heap(hprrho(iDens)),
     &                       Sfrac, result)
            call XCOPYR(heap(ptResults), iAtom, result, 1)
            if (WrnLev .ge. 10)
     &        call PrintLIDens(Sfrac, result, QFrac, FCTrMX)
          end if
        end do
        call CpAtPr(To, heap(ptResults))
        call FreHp(ptResults, ireal8(nAtom))
      end if
C
      if (QSite) then
        call Get8PLIDens(Dn, Dn(1), Dn(2), heap(hprrho(iDens)),
     &                   Site, result)
        call declar('RESULT',  'DP', ' ', dbcomp, result)
        if (WrnLev .ge. 5)
     &    call PrintLIDens(Site, result, QFrac, FCTrMX)
      end if
C
      return
      end
C}
