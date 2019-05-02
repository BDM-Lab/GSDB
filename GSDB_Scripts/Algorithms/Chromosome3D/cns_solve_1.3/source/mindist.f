C=======================================================================
C{
      subroutine NormFrac(nCfrac, Cfrac)
      implicit none
C I/O
      integer          nCfrac
      double precision  Cfrac(*)
C
C local
      integer  i
C
C begin
      do i = 1, nCfrac
        Cfrac(i) = mod(Cfrac(i), dfloat(1)) + 2.
        Cfrac(i) = mod(Cfrac(i), dfloat(1))
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine ShortFrac(nCfrac, Cfrac)
      implicit none
C I/O
      integer          nCfrac
      double precision  Cfrac(*)
C
C local
      integer  i
C
C begin
      do i = 1, nCfrac
        Cfrac(i) = mod(Cfrac(i), dfloat(1)) + 2.
        Cfrac(i) = mod(Cfrac(i), dfloat(1))
        if (Cfrac(i) .gt. 0.5) Cfrac(i) = Cfrac(i) - 1.
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine GetMaxDist2(FCTrMx, MaxDist2)
      implicit none
C I/O
      double precision  FCTrMx(3, 3)
      double precision  MaxDist2
C
C local
      integer           i, j, k, l
      double precision  Dfrac(3), Dorth(3), D2
C
C begin
      MaxDist2 = -1.
C
      do i =  0, 1
        Dfrac(1) = dfloat(i)
      do j = -1, 1
        Dfrac(2) = dfloat(j)
      do k = -1, 1
        Dfrac(3) = dfloat(k)
        if (i .ne. 0 .or. j .ne. 0 .or. k .ne. 0) then
          call TrVec3(FCTrMX, Dfrac, Dorth)
          D2 = 0.
          do l = 1, 3
            D2 = D2 + Dorth(l)**2
          end do
          if (MaxDist2 .lt. 0. .or. MaxDist2 .gt. D2)
     &      MaxDist2 = D2
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
      subroutine SMxMulCfrac(mSymMx, SymMx, iSymMx, STBF,
     &                       Cfrac, SymCfrac)
      implicit none
C I/O
      integer           mSymMx, SymMx(mSymMx, 3, 4), iSymMx, STBF
      double precision  Cfrac(3), SymCfrac(3)
C
C local
      integer  i
C
C begin
      do i = 1, 3
        SymCfrac(i) =   SymMx(iSymMx, i, 1) * Cfrac(1)
     &                + SymMx(iSymMx, i, 2) * Cfrac(2)
     &                + SymMx(iSymMx, i, 3) * Cfrac(3)
     &                + SymMx(iSymMx, i, 4) / dfloat(STBF)
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine SetEquiv(CoordX, CoordY, CoordZ,
     &                    nNCSym, mNCSym, NCSop,
     &                    nSymMx, mSymMx, STBF, SymMx,
     &                    CFTrMx,
     &                    Equiv)
      implicit none
C I/O
      double precision  CoordX, CoordY, CoordZ
      integer           nNCSym, mNCSym
      double precision  NCSop(mNCSym, 3, 4)
      integer           nSymMx, mSymMx, STBF
      integer           SymMx(mSymMx, 3, 4)
      double precision  CFTrMx(3, 3)
      double precision  Equiv(3, *)
C
C local
      integer           iEquiv, iNCS, iSymMx, i
      double precision  Corth(3), Cfrac(3)
C
C begin
      iEquiv = 0
C
      do iNCS = 1, nNCSym
        do i = 1, 3
          Corth(i) =   NCSop(iNCS, i, 1) * CoordX
     &               + NCSop(iNCS, i, 2) * CoordY
     &               + NCSop(iNCS, i, 3) * CoordZ
     &               + NCSop(iNCS, i, 4)
        end do
        call TrCoor(CFTrMX, Corth(1), Corth(2), Corth(3),
     &                      Cfrac(1), Cfrac(2), Cfrac(3))
        do iSymMx = 1, nSymMx
          iEquiv = iEquiv + 1
          call SMxMulCfrac(mSymMx, SymMx, iSymMx, STBF,
     &                     Cfrac, Equiv(1, iEquiv))
          call NormFrac(3, Equiv(1, iEquiv))
        end do
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine GetShortDist2(iAtom, Equivi,
     &                         jAtom, Equivj,
     &                         nNCSym, nSymMx,
     &                         FCTrMx, MaxDist2,
     &                         NBQSPC2, QSpecPos,
     &                         Dist2)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      integer           iAtom, jAtom
      double precision  Equivi(3, *), Equivj(3, *)
      integer           nNCSym, nSymMx
      double precision  FCTrMx(3, 3)
      double precision  MaxDist2
      double precision  NBQSPC2
      logical           QSpecPos
      double precision  Dist2
C
C local
      integer           jEq, i
      double precision  Dfrac(3), Dorth(3), D2
C
C begin
      Dist2 = MaxDist2
C
      jEq = 0
      if (iAtom .eq. jAtom) jEq = 1
      do while (jEq .lt. nNCSym * nSymMx)
        jEq = jEq + 1
        do i = 1, 3
          Dfrac(i) = Equivi(i, 1) - Equivj(i, jEq)
        end do
        call ShortFrac(3, Dfrac)
        call TrVec3(FCTrMX, Dfrac, Dorth)
        D2 = 0.
        do i = 1, 3
          D2 = D2 + Dorth(i)**2
        end do
        if ((     .not. QSpecPos
     &       .or. D2 .gt. NBQSPC2
     &       .or. jEq .gt. nSymMx)
     &      .and. Dist2 .gt. D2) Dist2 = D2
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine CheckDistances(nAtom, AtSele,
     &                          nNCSym, nSymMx,
     &                          ptEquiv,
     &                          FCTrMx, MaxDist2, MinDist2,
     &                          NBQSPC2, QSpecPos,
     &                          Result, WrnLev)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      integer           nAtom
      integer           AtSele(nAtom, 2)
      integer           nNCSym, nSymMx
      integer           ptEquiv(nAtom)
      double precision  FCTrMx(3, 3)
      double precision  MaxDist2, MinDist2
      double precision  NBQSPC2
      logical           QSpecPos
      logical           Result
      integer           WrnLev
C
C local
      integer           iAtom, jAtom
      double precision  Dist2
C
C begin
      Result = .false.
C
C Loop over all pairs (i,j) with i >= j and at at least one of i or j
C in set2.
      do iAtom = 1, nAtom
        if (AtSele(iAtom, 1) .ne. 0 .or. AtSele(iAtom, 2) .ne. 0) then
          do jAtom = iAtom, nAtom
            if (            AtSele(jAtom, 2) .ne. 0
     &          .or. (      AtSele(iAtom, 2) .ne. 0
     &                .and. AtSele(jAtom, 1) .ne. 0)) then
              call GetShortDist2(iAtom, heap(ptEquiv(iAtom)),
     &                           jAtom, heap(ptEquiv(jAtom)),
     &                           nNCSym, nSymMx,
     &                           FCTrMx, MaxDist2,
     &                           NBQSPC2, QSpecPos,
     &                           Dist2)
              if (WrnLev .ge. 10)
     &          write(6, '(1X, A, 2(1X, I8), 1X, F8.2)')
     &            'MinDistance', iAtom, jAtom, sqrt(Dist2)
              if (Dist2 .lt. MinDist2) return
            end if
          end do
        end if
      end do
C
      Result = .true.
C
      return
      end
C}
C=======================================================================
C{
      subroutine MinDistPa(nAtom, nAtSele, AtSele,
     &                     CoordX, CoordY, CoordZ,
     &                     MinDist2, QSpecPos, NBQSPC2)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      integer           nAtom
      integer           nAtSele(2), AtSele(nAtom, 2)
      double precision  CoordX(*), CoordY(*), CoordZ(*)
      double precision  MinDist2
      logical           QSpecPos
      double precision  NBQSPC2
C
C This is the parser for MinDist
C
C begin
      nAtSele(1) = 0
      nAtSele(2) = 0
      MinDist2 = 1.
      QSpecPos = .false.
C
C parsing
      CALL PUSEND('MinDistance>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('MinDistance>')
      CALL MISCOM('MinDistance>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
CCC modification ATB 4/27/08
      CALL CNSHELP('cns-mindist')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SET1') THEN
        call SelctA(AtSele(1, 1), nAtSele(1),
     &    CoordX, CoordY, CoordZ, .true.)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SET2') THEN
        call SelctA(AtSele(1, 2), nAtSele(2),
     &    CoordX, CoordY, CoordZ, .true.)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MIND') THEN
        MinDist2 = 1.
        CALL NEXTF('MinDistance=', MinDist2)
        MinDist2 = MinDist2 * MinDist2
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SPEC') THEN
        CALL NEXTLO('SpecialPositions=', QSpecPos)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
        WRITE(6,'(2A)') ' ----------------------------mindistance-',
     &   'parameters----------------------------'
        write(6, '(1X, A, I6)') 'Number of atoms in set 1: ',
     &    nAtSele(1)
        write(6, '(1X, A, I6)') 'Number of atoms in set 2: ',
     &    nAtSele(2)
        write(6, '(1X, A, F8.2)') 'MinDistance = ', sqrt(MinDist2)
        call EchoLogical('SpecialPositions', QSpecPos)
        if (QSpecPos) then
          write(6, '(1X, A, F8.2)')
     &      'Special position tolerance = ', sqrt(NBQSPC2)
        end if
        WRITE(6,'(2A)') ' ----------------------------------------',
     &   '--------------------------------------'
C=====================================================================
      ELSE
        CALL CHKEND('MinDistance>',DONE)
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
      subroutine MinDist(nAtom, CoordX, CoordY, CoordZ,
     &                   nNCSym, mNCSym, NCSop,
     &                   nSymMx, mSymMx, STBF, SymMx,
     &                   CFTrMx, FCTrMx, NBQSPC2,
     &                   AtSele, ptEquiv)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      integer           nAtom
      double precision  CoordX(*), CoordY(*), CoordZ(*)
      integer           nNCSym, mNCSym
      double precision  NCSop(mNCSym, 3, 4)
      integer           nSymMx, mSymMx, STBF
      integer           SymMx(mSymMx, 3, 4)
      double precision  CFTrMx(3, 3), FCTrMx(3, 3)
      double precision  NBQSPC2
      integer           AtSele(nAtom, 2)
      integer           ptEquiv(nAtom)
C
C MinDist main procedure.
C
C local
      integer           nAtSele(2)
      double precision  MaxDist2, MinDist2
      logical           QSpecPos
      integer           iSet, iAtom, mEquiv, nUnknownCoord
      logical           Result
      character         buf*80
      double precision  dpval
      double complex    dcval
C
C begin
C Call the parser
      dpval = dfloat(0)
      dcval = dcmplx(dpval, dpval)
C
      call MinDistPa(nAtom, nAtSele, AtSele,
     &               CoordX, CoordY, CoordZ,
     &               MinDist2, QSpecPos, NBQSPC2)
C
      do iSet = 1, 2
        if (nAtSele(iSet) .lt. 1) then
          write(buf, '(A, I1, A)')
     &      'No atoms selected for set ', iSet, '.'
          call WrnDie(0, 'MinDistance', buf)
          return
        end if
      end do
C
      do iAtom = 1, nAtom
        ptEquiv(iAtom) = 0
      end do
C
      mEquiv = 3 * nSymMx * nNCSym
      nUnknownCoord = 0
C
      do iAtom = 1, nAtom
        ptEquiv(iAtom) = 0
        if (AtSele(iAtom, 1) .ne. 0 .or. AtSele(iAtom, 2) .ne. 0) then
          if (     CoordX(iAtom) .eq. 9999.
     &        .or. CoordY(iAtom) .eq. 9999.
     &        .or. CoordZ(iAtom) .eq. 9999.) then
            nUnknownCoord = nUnknownCoord + 1
          else
            ptEquiv(iAtom) = AllHp(ireal8(mEquiv))
            call SetEquiv(CoordX(iAtom), CoordY(iAtom), CoordZ(iAtom),
     &                    nNCSym, mNCSym, NCSop,
     &                    nSymMx, mSymMx, STBF, SymMx,
     &                    CFTrMx,
     &                    heap(ptEquiv(iAtom)))
          end if
        end if
      end do
C
      if (nUnknownCoord .gt. 0) then
        write(buf, '(I8, A)')
     &    nUnknownCoord, ' selected atoms have unknown coordinates.'
        call WrnDie(0, 'MinDistance', buf)
      else
        call GetMaxDist2(FCTrMx, MaxDist2)
        call CheckDistances(nAtom, AtSele,
     &                      nNCSym, nSymMx,
     &                      ptEquiv,
     &                      FCTrMx, MaxDist2, MinDist2,
     &                      NBQSPC2, QSpecPos,
     &                      Result, WrnLev)
        if (Result) then
          call declar('RESULT', 'LO', 'TRUE',  dcval, dpval)
        else
          call declar('RESULT', 'LO', 'FALSE', dcval, dpval)
        end if
      end if
C
      do iAtom = 1, nAtom
        if (ptEquiv(iAtom) .ne. 0)
     &    call FreHp(ptEquiv(iAtom), ireal8(mEquiv))
      end do
C
      return
      end
C}
C=======================================================================
C{
      SUBROUTINE MINDWR
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
C MinDist wrapper routine.
C
C local
      INTEGER           PTATOMSELE, PTPTEQUIV
      DOUBLE PRECISION  NBQSPC2
C
C begin
      IF (LNCSRE) THEN
        CALL WRNDIE(0, 'MinDistance', 'only strict NCS is allowed')
        RETURN
      END IF
C
      IF (XNNSYM .GT. 1 .AND. .NOT. LNCSST) THEN
        CALL WRNDIE(-5, 'MinDistance', 'Fatal Coding Error')
        CALL DIE
      END IF
C
      PTATOMSELE = ALLHP(INTEG4(NATOM * 2))
      PTPTEQUIV  = ALLHP(INTEG4(NATOM))
C
      NBQSPC2 = NBQSPC**2
      CALL MINDIST(NATOM, X, Y, Z,
     &             XNNSYM, MNCSYM, NCSOP,
     &             XRNSYM, XRMSYM, XRSYTH, XRSYMM,
     &             XRTR, XRINTR, NBQSPC2,
     &             HEAP(PTATOMSELE), HEAP(PTPTEQUIV))
C
      CALL FREHP(PTATOMSELE, INTEG4(NATOM * 2))
      CALL FREHP(PTPTEQUIV,  INTEG4(NATOM))
C
      RETURN
      END
C}
