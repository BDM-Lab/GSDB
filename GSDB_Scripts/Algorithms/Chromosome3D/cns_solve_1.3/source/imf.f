C{
      logical function CmpRotMx(mSymMx, SymMx, iS, jS)
      implicit none
C I/O
      integer  mSymMx, SymMx(mSymMx, 3, 4), iS, jS
C
C local
      integer  ir, ic
C
C begin
      CmpRotMx = .false.
C
      do ir = 1, 3
      do ic = 1, 3
        if (SymMx(iS, ir, ic) .ne. SymMx(jS, ir, ic)) return
      end do
      end do
C
      CmpRotMx = .true.
C
      return
      end
C}
C=======================================================================
C{
      subroutine SetSymFlags(mSymMx, nSymMx, SymMx, SymFlags)
      implicit none
C I/O
      integer  mSymMx, nSymMx, SymMx(mSymMx, 3, 4), SymFlags(*)
C
C local
      integer  iS, jS, nCent
C
C externals
      logical   CmpRotMx
      external  CmpRotMx
C
C begin
      do iS = 1, nSymMx
        SymFlags(iS) = 0
      end do
C
      nCent = 0
C
      do iS = 1, nSymMx - 1
        if (SymFlags(iS) .eq. 0) then
          do jS = iS + 1, nSymMx
            if (SymFlags(jS) .eq. 0) then
              if (CmpRotMx(mSymMx, SymMx, iS, jS)) then
                SymFlags(jS) = iS
                nCent = nCent + 1
              end if
            end if
          end do
        end if
      end do
C
      if (mod(nSymMx, nSymMx - nCent) .ne. 0) then
        call WrnDie(-5, 'SetSymFlags', 'Fatal Coding Error.')
        call Die
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine CalcIMFmap(Rn, Rn1, Rn2, FlagMap,
     &                      PatMap, IMFmap,
     &                      nAtoms, Cx, Cy, Cz, Wght, CFTrMx,
     &                      mSymMx, nSymMx, SymMx, STBF, SymFlags,
     &                      UnitWeights,
     &                      ChkInDep)
      implicit none
C I/O
      integer           Rn(3), Rn1, Rn2, FlagMap(*)
      real              PatMap(0:Rn1-1, 0:Rn2-1, 0:*), IMFmap(*)
      integer           nAtoms
      double precision  Cx(*), Cy(*), Cz(*), Wght(*)
      double precision  CFTrMx(3, 3)
      integer           mSymMx, nSymMx, SymMx(mSymMx, 3, 4), STBF
      integer           SymFlags(*)
      logical           UnitWeights
      logical           ChkInDep
C
C Calculate Image Seeking Minimum Function
C Algorithm:
C   (IMF map has to be initialized with max value of Patterson map)
C   Loop over known atoms
C     Loop over SymOps -> generate symmetry mate of atoms
C       Shift origin of Patterson to that position
C         Loop over grid points of IMF map
C            Get density of shifted Patterson map at each grid point
C            by doing a 8-point linear interpolation.
C            Update IMF map if that density is less than current
C            value in IMF map.
C
C 8-point linear interpolation: determine the grid box in the
C Patterson map where the grid point of the IMF map is located.
C Compute the weighted average of the 8 vertices of the grid
C box. The naive approach is to do this 8-point interpolation
C including determination of the grid box for each grid point
C of the IMF map. However, for trigonal and hexagonal space
C groups the resulting IMF map will in general not have the
C space group symmetry.
C The algorithm below does two things different than the naive
C approach: The obvious optimization is that the grid box is
C determined only once per position. For each grid point of the
C IMF map, the box can simply be shifted, since the relative
C position of the IMF grid point in the Patterson grid box is
C always the same for a given atom position (and so are the
C weights for the averaging which are stored in Xw below).
C The second measure is to make sure that the resulting IMF map
C has exactly space group symmetry in all cases. Instead of
C applying the symmetry operations to the atom positions and
C then determining the Patterson grid box, the symmetry
C operations are applied to the whole grid box (8 individual
C grid points). The transformed grid box can be shifted around
C as usual. The weights stored in Xw can be used for all
C symmetry mates without modification. (Caution: the SymOps
C have to be applied to -position of each of the 8 vertices.
C See the two "Also make sign..." comments below.)
C
C local
      integer           LCM, ColFac(4), GridMis(3)
      integer           iAtom, iSymMx, iMap, ia, ib, ic, i8, i
      integer           Ri(3), Va(3, 0:7), Vs(3, 0:7), Vg(3, 0:7)
      double precision  AsyF(3), Xw(3, 0:1), W(0:7), LIDens
      real              wLIDens
C
C begin
      call SetColFac(Rn, STBF, LCM, ColFac, GridMis)
C
      do iAtom = 1, nAtoms
        call TrCoor(CFTrMX, Cx(iAtom), Cy(iAtom), Cz(iAtom),
     &                      AsyF(1), AsyF(2), AsyF(3))
        do i = 1, 3
          AsyF(i) = -AsyF(i)
        end do
C
        call LICorner(AsyF, Rn, Ri, Xw)
C
        i8 = 0
        do ic = 0, 1
        do ib = 0, 1
        do ia = 0, 1
C         Also make sign of Va right for application of SymOps
          Va(1, i8) = -(Ri(1) + ia)
          Va(2, i8) = -(Ri(2) + ib)
          Va(3, i8) = -(Ri(3) + ic)
          W(i8) = Xw(1, ia) * Xw(2, ib) * Xw(3, ic)
          i8 = i8 + 1
        end do
        end do
        end do
C
        do iSymMx = 1, nSymMx
          if (SymFlags(iSymMx) .eq. 0) then
            do i8 = 0, 7
              call SymEquiv(SymMx, mSymMx, iSymMx,
     &                      LCM, ColFac, Va(1, i8), Vs(1, i8))
              do i = 1, 3
                if (Vs(i, i8) .ne. 0) then
                  if (mod(Vs(i, i8), ColFac(i)) .eq. 0) then
C                   Also make sign of Vs right
                    Vs(i, i8) = Rn(i) - Vs(i, i8) / ColFac(i)
                  else
                    call WrnDie(-5, 'IMF',
     &                'Symmetry mate of grid point is not on grid.')
                    call Die
                  end if
                end if
              end do
            end do
C
            iMap = 1
            do ic = 0, Rn(3) - 1
              do i8 = 0, 7
                Vg(3, i8) = mod(Vs(3, i8) + ic, Rn(3))
              end do
            do ib = 0, Rn(2) - 1
              do i8 = 0, 7
                Vg(2, i8) = mod(Vs(2, i8) + ib, Rn(2))
              end do
            do ia = 0, Rn(1) - 1
              if (FlagMap(iMap) .le. 0 .or. ChkInDep) then
                LIDens = dfloat(0)
                do i8 = 0, 7
                  Vg(1, i8) = mod(Vs(1, i8) + ia, Rn(1))
                  LIDens = LIDens
     &              + PatMap(Vg(1, i8), Vg(2, i8), Vg(3, i8)) * W(i8)
                end do
                if (UnitWeights) then
                  wLIDens = LIDens
                else
                  wLIDens = LIDens / Wght(iAtom)
                end if
                IMFmap(iMap) = min(IMFmap(iMap), wLIDens)
              end if
              iMap = iMap + 1
            end do
            end do
            end do
          end if
        end do
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine ShowIMFxyzw(nAtoms, Cx, Cy, Cz, Wght, CFTrMx)
      implicit none
C I/O
      integer           nAtoms
      double precision  Cx(*), Cy(*), Cz(*), Wght(*)
      double precision  CFTrMx(3, 3)
C
C local
      integer           iAtom, i
      double precision  FX(3)
C
C begin
      write(6, '(1X, A, i5)')
     &  'IMF: Number of positions = ', nAtoms
      write(6, '(1X, 2A)')
     &  'IMF:    Cartesian Coordinates',
     &  '        Fractional Coordinates    Weight'
C
      do iAtom = 1, nAtoms
        call TrCoor(CFTrMX, Cx(iAtom), Cy(iAtom), Cz(iAtom),
     &                      FX(1), FX(2), FX(3))
        write(6, '(1X, A, 3(1X, F8.2), 1X, 3(1X, F8.5), 1X, G14.6)')
     &    'IMF: ',
     &    Cx(iAtom), Cy(iAtom), Cz(iAtom),
     &    (FX(i), i = 1, 3), Wght(iAtom)
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine SetIMFxyzw(CoordX, CoordY, CoordZ, Occ,
     &                      xrSm, xrSa, xrSc, xrFP,
     &                      nSelAtoms, xrAtom, xrIndx,
     &                      Cx, Cy, Cz, Wght)
      implicit none
C I/O
      double precision  CoordX(*), CoordY(*), CoordZ(*), Occ(*)
      integer           xrSm
      double precision  xrSa(xrSm, 4), xrSc(xrSm), xrFP(xrSm)
      integer           nSelAtoms, xrAtom(*), xrIndx(*)
      double precision  Cx(*), Cy(*), Cz(*)
      double precision  Wght(*)
C
C local
      integer  iSelAtom
C
C begin
      do iSelAtom = 1, nSelAtoms
        Cx(iSelAtom) = CoordX(xrAtom(iSelAtom))
        Cy(iSelAtom) = CoordY(xrAtom(iSelAtom))
        Cz(iSelAtom) = CoordZ(xrAtom(iSelAtom))
C
        Wght(iSelAtom) =   xrSc(xrIndx(iSelAtom))
     &                   + xrSa(xrIndx(iSelAtom), 1)
     &                   + xrSa(xrIndx(iSelAtom), 2)
     &                   + xrSa(xrIndx(iSelAtom), 3)
     &                   + xrSa(xrIndx(iSelAtom), 4)
     &                   + xrFP(xrIndx(iSelAtom))
        Wght(iSelAtom) = Wght(iSelAtom) * Occ(xrAtom(iSelAtom))
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine IMFpa(iPatMap, iIMFMap,
     &                 xRhoNum, xRhoNam,
     &                 nAtSele, AtSele, CoordX, CoordY, CoordZ,
     &                 UnitWeights,
     &                 ChkInDep, BailOut)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      integer           iPatMap, iIMFMap
      integer           xRhoNum
      character*(*)     xRhoNam(*)
      integer           nAtSele, AtSele(*)
      double precision  CoordX(*), CoordY(*), CoordZ(*)
      logical           UnitWeights
      logical           ChkInDep, BailOut
C
C This is the parser for IMF
C
C begin
      iPatMap = 0
      iIMFMap = 0
      nAtSele = 0
      UnitWeights = .false.
      ChkInDep = .false.
      BailOut = .true.
C
C parsing
      CALL PUSEND('IMF>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('IMF>')
      CALL MISCOM('IMF>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-imf')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PATM') THEN
        call GetiObj(iPatMap, xRhoNum, xRhoNam, ' ',
     &               'PATMap=', 'real', 'IMF')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'IMFM') THEN
        call GetiObj(iIMFMap, xRhoNum, xRhoNam, ' ',
     &               'IMFMap=', 'real', 'IMF')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ATOM') THEN
        call SelctA(AtSele, nAtSele, CoordX, CoordY, CoordZ, .true.)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'UNIT') THEN
        CALL NEXTLO('UnitWeights=', UnitWeights)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'CHKI') THEN
        CALL NEXTLO('CHKIndep=', ChkInDep)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'BAIL') THEN
        CALL NEXTLO('BAILout=', BailOut)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
        WRITE(6,'(2A)') ' --------------imf-parameters-',
     &   '--------------------------------------------------'
        call EchoObj(iPatMap, xRhoNam, 'PATMap=')
        call EchoObj(iIMFMap, xRhoNam, 'IMFMap=')
        write(6, '(1X, A, I6)') 'Number of atoms selected: ', nAtSele
        call EchoLogical('UnitWeights', UnitWeights)
        WRITE(6,'(2A)') ' -----------------------------',
     &   '--------------------------------------------------'
C=====================================================================
      ELSE
        CALL CHKEND('IMF>',DONE)
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
      subroutine IMF(QHerm, na, nb, nc, nrRho,
     &               xRhoNum, xRhoNam, hprrho,
     &               nAtom, xrNAtF, xrAtoF, xrIndF,
     &               CoordX, CoordY, CoordZ,
     &               xrSm, xrSa, xrSc, xrFP, xrFDP,
     &               xrNsym, XRMSYM, STBF, XRSYMM, XRITSY,
     &               CFTrMx, FCTrMx,
     &               nNCSsym,
     &               AtSele)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'flagmap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      logical           QHerm
      integer           na, nb, nc, nrRho
      integer           xRhoNum
      character*(*)     xRhoNam(*)
      integer           hprrho(*)
      integer           nAtom
      integer           xrNAtF, xrAtoF(*), xrIndF(*)
      double precision  CoordX(*), CoordY(*), CoordZ(*)
      integer           xrSm
      double precision  xrSa(xrSm, 4), xrSc(xrSm)
      double precision  xrFP(*), xrFDP(*)
      integer           xrNsym, XRMSYM, STBF
      integer           XRSYMM(XRMSYM, 3, 4), XRITSY(XRMSYM, 3, 3)
      double precision  CFTrMx(3, 3), FCTrMx(3, 3)
      integer           nNCSsym
      integer           AtSele(*)
C
C XRAY IMF main procedure.
C
C local
      integer  iPatMap, iIMFMap
      integer  PRn
      integer  nAtSele, nNorm, nAnom
      logical  UnitWeights
      logical  ChkInDep, BailOut
      integer  SymFlags(mSymMx)
C
C pointer
      integer  ptAnomFlag, ptFqs, ptAtom, ptIndx
      integer  ptCx, ptCy, ptCz, ptWght
C
      double precision  tmEntr, tmExit
C
C parameters
      real      SZero
      parameter(SZero=0.0)
C
C externals
      integer   lnblnk
      external  lnblnk
C
C begin
C Call the parser
      call IMFpa(iPatMap, iIMFMap,
     &           xRhoNum, xRhoNam,
     &           nAtSele, AtSele,
     &           CoordX, CoordY, CoordZ,
     &           UnitWeights,
     &           ChkInDep, BailOut)
C
      if (Timer .gt. 0) call VCPU(tmEntr)
C
      if (iPatMap .lt. 1 .and. iIMFMap .lt. 1) return
C
      if (iPatMap .lt. 1) then
        call WrnDie(0, 'IMF', 'PATMap not specified.')
        return
      end if
C
      if (nAtSele .lt. 1) then
        call WrnDie(0, 'IMF', 'zero atoms selected.')
        return
      end if
C
      if (hprrho(iPatMap) .eq. 0) then
        write(6, '(3A)') ' %IMF-ERR: Map "',
     &    xrhonam(iPatMap)(1:lnblnk(xrhonam(iPatMap))),
     &    '" not defined.'
        call WrnDie(0, 'IMF', 'object undefined.')
        return
      end if
C
      if (iIMFMap .lt. 1) then
        call WrnDie(0, 'IMF', 'IMFMap not specified.')
        return
      end if
C
      if (hprrho(iIMFMap) .eq. 0) then
        write(6, '(3A)') ' %IMF-ERR: Map "',
     &    xrhonam(iIMFMap)(1:lnblnk(xrhonam(iIMFMap))),
     &    '" not defined.'
        call WrnDie(0, 'IMF', 'object undefined.')
        return
      end if
C
      if (xrNsym .ne. 1) then
        call WrnDie(0, 'IMF',
     &    'This module can only be used in P1.')
        return
      end if
C
      if (     .not. ValidFlagMap
     &    .or. RnFM(1) .ne. na
     &    .or. RnFM(2) .ne. nb
     &    .or. RnFM(3) .ne. nc) then
        call WrnDie(0, 'IMF',
     &    'FlagMap undefined or out of date: call FMAP first.')
        return
      end if
C
      PRn = RnFM(1) * RnFM(2) * RnFM(3)
C
      if (nrRho .ne. PRn) then
        call WrnDie(-5, 'IMF', 'Fatal Coding Error.')
        call Die
      end if
C
      if (nNCSsym .gt. 1) then
        write(6, '(2A)')
     &    ' %IMF-WRN: Non-crystallographic symmetry',
     &    ' not applied to selected atoms!'
      end if
C
      ptAnomFlag = AllHp(integ4(nAtom))
      ptAtom     = AllHp(integ4(nAtom))
      ptIndx     = AllHp(integ4(nAtom))
      ptFqs      = AllHp(ireal8(nAtom))
C
      call Fill4(heap(ptAnomFlag), nAtom, 0)
C
      call xrAssoc(AtSele,
     &             heap(ptAnomFlag), xrAtoF, xrIndF,
     &             heap(ptAtom), heap(ptIndx), heap(ptFqs),
     &             nNorm, nAnom,
     &             3, 'IMF', .true., xrNAtF, xrFDP, QHerm,
     &             xrNsym, XRMSYM, STBF, XRSYMM, XRITSY,
     &             CFTrMx, FCTrMx)
C
      if (nNorm .ne. nAtSele .or. nNorm .ne. nAnom) then
        call WrnDie(-5, 'IMF', 'Fatal Coding Error.')
        call Die
      end if
C
      call FreHp(ptAnomFlag, integ4(nAtom))
C
      ptCx   = AllHp(ireal8(nAtSele))
      ptCy   = AllHp(ireal8(nAtSele))
      ptCz   = AllHp(ireal8(nAtSele))
      ptWght = AllHp(ireal8(nAtSele))
C
      call SetIMFxyzw(CoordX, CoordY, CoordZ, heap(ptFqs),
     &                xrSm, xrSa, xrSc, xrFP,
     &                nAtSele, heap(ptAtom), heap(ptIndx),
     &                heap(ptCx), heap(ptCy), heap(ptCz), heap(ptWght))
C
      call FreHp(ptAtom,   integ4(nAtom))
      call FreHp(ptIndx,   integ4(nAtom))
      call FreHp(ptFqs,    ireal8(nAtom))
C
      if (WrnLev .ge. 5)
     &  call ShowIMFxyzw(nAtSele,
     &                   heap(ptCx), heap(ptCy), heap(ptCz),
     &                   heap(ptWght), CFTrMx)
C
      call SetSymFlags(mSymMx, nSymMx, SymMx, SymFlags)
C
      call CalcIMFmap(RnFM, RnFM(1), RnFM(2), heap(ptFlagMap),
     &                heap(hprrho(iPatMap)), heap(hprrho(iIMFMap)),
     &                nAtSele, heap(ptCx), heap(ptCy), heap(ptCz),
     &                heap(ptWght), CFTrMx,
     &                mSymMx, nSymMx, SymMx, STBF, SymFlags,
     &                UnitWeights,
     &                ChkInDep)
C
      call FreHp(ptCx,       ireal8(nAtSele))
      call FreHp(ptCy,       ireal8(nAtSele))
      call FreHp(ptCz,       ireal8(nAtSele))
      call FreHp(ptWght,     ireal8(nAtSele))
C
      call MapInd(heap(hprrho(iIMFMap)), heap(ptFlagMap), RnFM,
     &            ChkInDep, .false., BailOut, WrnLev, 'IMF')
C
      if (WrnLev .ge. 5) then
        call MapStat(PRn, heap(hprrho(iPatMap)), 'IMF: PATmap ', .true.)
        call MapStat(PRn, heap(hprrho(iIMFMap)), 'IMF: IMFmap ', .true.)
      end if
C
      if (Timer .gt. 0) then
        call VCPU(tmExit)
        write(6, '(1X, A, F10.4, A)')
     &    'IMF: CPU-time: ', tmExit - tmEntr, ' seconds total'
      end if
C
      return
      end
C}
