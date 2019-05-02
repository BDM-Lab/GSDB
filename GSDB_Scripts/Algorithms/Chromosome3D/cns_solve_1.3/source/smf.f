C{
      subroutine HarkEquiv(HarkOp,
     &                     LCM, ColFac, Ri, Re)
      implicit none
C I/O
      integer  HarkOp(3, 4)
      integer  LCM, ColFac(4)
      integer  Ri(3), Re(3)
C
C Determine the indices Re for given indices Ri and Harker operator
C
C local
      integer  i
C
C begin
      do i = 1, 3
        Re(i) =   HarkOp(i, 1) * ColFac(1) * Ri(1)
     &          + HarkOp(i, 2) * ColFac(2) * Ri(2)
     &          + HarkOp(i, 3) * ColFac(3) * Ri(3)
     &          + HarkOp(i, 4) * ColFac(4)
        Re(i) = mod(Re(i), LCM)
        if (Re(i) .lt. 0) Re(i) = Re(i) + LCM
      end do
C
      return
      end
C}
C=======================================================================
C{
      real function SMFatX(Rn, Ri, PatMap,
     &                     SgOrderP,
     &                     nHarkOps, HarkOps,
     &                     UnitWeights,
     &                     LCM, ColFac, GridMis)
      implicit none
C I/O
      integer  Rn(3), Ri(3)
      real     PatMap(*)
      integer  SgOrderP
      integer  nHarkOps, HarkOps(13, *)
      logical  UnitWeights
      integer  LCM, ColFac(4), GridMis(3)
C
C local
      integer  iHarkOps, Re(3), eMap, m, i
      logical  OK
C
C Determine SMF value at a given position X, where X(i)=Ri(i)/Rn(i)
C
C external
      integer   f3dix
      external  f3dix
C
C begin
      SMFatX = PatMap(1)
      if (SgOrderP .ne. 0) SMFatX = SMFatX / float(SgOrderP)
C
      m = 1
C
      do iHarkOps = 1, nHarkOps
        call HarkEquiv(HarkOps(2, iHarkOps),
     &                 LCM, ColFac, Ri, Re)
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
          if (.not. UnitWeights) m = HarkOps(1, iHarkOps)
          SMFatX = min(SMFatX, PatMap(eMap) / float(m))
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine CalcSMFmap(Rn, FlagMap,
     &                      PatMap, SMFmap,
     &                      SgOrderP,
     &                      nHarkOps, HarkOps,
     &                      UnitWeights,
     &                      ChkInDep)
      implicit none
C I/O
      integer  Rn(3), FlagMap(*)
      real     PatMap(*), SMFmap(*)
      integer  SgOrderP
      integer  nHarkOps, HarkOps(13, *)
      logical  UnitWeights
      logical  ChkInDep
C
C Calculation of the Symmetry Minimum Function
C
C local
      integer    PRn, iMap, ia, ib, ic, i
      integer    Ri(3)
      integer    LCM, ColFac(4), GridMis(3)
      real       PatMax, PatMin
      logical    OK
      character  xyz*3
C
C external
      real      SMFatX
      external  SMFatX
C
C begin
      xyz = 'xyz'
C
      PRn = Rn(3) * Rn(2) * Rn(1)
C
      PatMin = PatMap(1)
      PatMax = PatMap(1)
C
      do iMap = 2, PRn
        PatMin = min(PatMin, PatMap(iMap))
        PatMax = max(PatMax, PatMap(iMap))
      end do
C
      call SetColFac(Rn, 12, LCM, ColFac, GridMis)
C
      iMap = 1
      do ic = 0, Rn(3) - 1
         Ri(3) = ic
      do ib = 0, Rn(2) - 1
         Ri(2) = ib
      do ia = 0, Rn(1) - 1
         Ri(1) = ia
        if (FlagMap(iMap) .le. 0 .or. ChkInDep) then
          SMFmap(iMap) = SMFatX(Rn, Ri, PatMap,
     &                          SgOrderP, nHarkOps, HarkOps,
     &                          UnitWeights,
     &                          LCM, ColFac, GridMis)
        end if
        iMap = iMap + 1
      end do
      end do
      end do
C
      OK = .true.
C
      do i = 1, 3
        if (GridMis(i) .ne. 0) then
          write(6, '(1X, 3A)')
     &      '%SMF-ERR: Some Harker vectors not on grid ',
     &      xyz(i:i), '.'
          OK = .false.
        end if
      end do
C
      if (.not. OK)
     &  call WrnDie(0, 'SMF', 'Inappropriate grid.')
C
      return
      end
C}
C=======================================================================
C{
      subroutine SMFpa(iPatMap, iSMFMap,
     &                 SgOrderP,
     &                 mHarkOps, nHarkOps, HarkOps,
     &                 xRhoNum, xRhoNam,
     &                 UnitWeights,
     &                 ChkInDep, BailOut)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      integer        iPatMap, iSMFMap
      integer        SgOrderP
      integer        mHarkOps, nHarkOps, HarkOps(13, *)
      integer        xRhoNum
      character*(*)  xRhoNam(*)
      logical        UnitWeights
      logical        ChkInDep, BailOut
C
C This is the parser for SMF
C
C local
      integer  iHarkOps, i, j, nOps, Mx(3, 4), strused
      logical  OK
C
C externals
      integer   DECODI
      external  DECODI
C
C begin
      iPatMap = 0
      iSMFMap = 0
      SgOrderP = 1
      nHarkOps = 0
      UnitWeights = .false.
      ChkInDep = .false.
      BailOut = .true.
C
C parsing
      CALL PUSEND('SMF>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('SMF>')
      CALL MISCOM('SMF>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-smf')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PATM') THEN
        call GetiObj(iPatMap, xRhoNum, xRhoNam, ' ',
     &               'PATMap=', 'real', 'SMF')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SMFM') THEN
        call GetiObj(iSMFMap, xRhoNum, xRhoNam, ' ',
     &               'SMFMap=', 'real', 'SMF')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SGOR') THEN
        CALL NEXTI('SGORderP=', SgOrderP)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'HARK') THEN
        CALL NEXTST('HARKerOperator=', WD)
        i = index(WD(1:WDLEN), '(')
        j = index(WD(1:WDLEN), ')')
        if (i .lt. 2 .or. i .gt. j) then
          call WrnDie(0, 'SMF', 'Illegal format for HARKerOperator.')
        else
          WD(i:i) = ' '
          WD(j:j) = ' '
          nOps = DECODI(WD(1:i), i, OK)
          if (.not. OK) then
            call WrnDie(0, 'SMF', 'Illegal HARKerOperator.')
          else
            call xyzparse(WD(i:WDLEN), WDLEN-i+1, Mx, 1, 12, strused)
            if (strused .lt. 0) then
              call WrnDie(0, 'SMF', 'Illegal HARKerOperator.')
            else if (nHarkOps .eq. mHarkOps) then
              call WrnDie(0, 'SMF', 'Too many Harker operators.')
            else
              nHarkOps = nHarkOps + 1
              HarkOps(1, nHarkOps) = nOps
              call CopyI4(Mx, HarkOps(2, nHarkOps), 12)
            end if
          end if
        end if
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
        WRITE(6,'(2A)') ' --------------smf-parameters-',
     &   '--------------------------------------------------'
        call EchoObj(iPatMap, xRhoNam, 'PATMap=')
        call EchoObj(iSMFMap, xRhoNam, 'SMFMap=')
        write(6, '(1X, A, I2)') 'SGORderP=', SgOrderP
        write(6, '(1X, A, I2)') 'Number of unique Harker operators = ',
     &    nHarkOps
        do iHarkOps = 1, nHarkOps
          write(6, '(1X, A, 13I3, A)')
     &      'HARKerOperator=(',
     &      (HarkOps(i, iHarkOps), i = 1, 13),
     &      ')'
        end do
        call EchoLogical('UnitWeights', UnitWeights)
        WRITE(6,'(2A)') ' -----------------------------',
     &   '--------------------------------------------------'
C=====================================================================
      ELSE
        CALL CHKEND('SMF>',DONE)
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
      subroutine SMF(na, nb, nc, nrRho, niRho, QHerm,
     &               xRhoNum, xRhoNam, hprrho, hpirho,
     &               xrNsym)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'flagmap.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      integer        na, nb, nc, nrRho, niRho
      logical        QHerm
      integer        xRhoNum
      character*(*)  xRhoNam(*)
      integer        hprrho(*), hpirho(*)
      integer        xrNsym
C
C XRAY SMF main procedure.
C
C local
      integer    iPatMap, iSMFMap
      integer    SgOrderP
      integer    mHarkOps, nHarkOps
      parameter (mHarkOps = 34)
      integer     HarkOps(13, mHarkOps)
      integer    PRn
      logical    UnitWeights
      logical    ChkInDep, BailOut
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
      call SMFpa(iPatMap, iSMFMap,
     &           SgOrderP,
     &           mHarkOps, nHarkOps, HarkOps,
     &           xRhoNum, xRhoNam,
     &           UnitWeights,
     &           ChkInDep, BailOut)
C
      if (Timer .gt. 0) call VCPU(tmEntr)
C
      if (iPatMap .lt. 1 .and. iSMFMap .lt. 1) return
C
      if (iPatMap .lt. 1) then
        call WrnDie(0, 'SMF', 'PATMap not specified.')
        return
      end if
C
      if (hprrho(iPatMap) .eq. 0) then
        write(6, '(3A)') ' %SMF-ERR: Map "',
     &    xrhonam(iPatMap)(1:lnblnk(xrhonam(iPatMap))),
     &    '" not defined.'
        call WrnDie(0, 'SMF', 'object undefined.')
        return
      end if
C
      if (iSMFMap .lt. 1) then
        call WrnDie(0, 'SMF', 'SMFMap not specified.')
        return
      end if
C
      if (xrNsym .ne. 1) then
        call WrnDie(0, 'SMF',
     &    'This module can only be used in P1.')
        return
      end if
C
      if (     .not. ValidFlagMap
     &    .or. RnFM(1) .ne. na
     &    .or. RnFM(2) .ne. nb
     &    .or. RnFM(3) .ne. nc) then
        call WrnDie(0, 'SMF',
     &    'FlagMap undefined or out of date: call FMAP first.')
        return
      end if
C
      PRn = RnFM(1) * RnFM(2) * RnFM(3)
C
      if (nrRho .ne. PRn) then
        call WrnDie(-5, 'SMF', 'Fatal Coding Error.')
        call Die
      end if
C
      if (hprrho(iSMFMap) .eq. 0) then
        call xMapAl(hprrho(iSMFMap), hpirho(iSMFMap), QHerm,
     &              nrRho, niRho)
      else if (.not. QHerm) then
        call FillR4(heap(hpirho(iSMFMap)), niRho, SZERO)
      end if
C
      call CalcSMFmap(RnFM, heap(ptFlagMap),
     &                heap(hprrho(iPatMap)),
     &                heap(hprrho(iSMFMap)),
     &                SgOrderP,
     &                nHarkOps, HarkOps,
     &                UnitWeights,
     &                ChkInDep)
C
      call MapInd(heap(hprrho(iSMFMap)), heap(ptFlagMap), RnFM,
     &            ChkInDep, .false., BailOut, WrnLev, 'SMF')
C
      if (WrnLev .ge. 5) then
        call MapStat(PRn, heap(hprrho(iPatMap)), 'SMF: PATmap ', .true.)
        call MapStat(PRn, heap(hprrho(iSMFMap)), 'SMF: SMFmap ', .true.)
      end if
C
      if (Timer .gt. 0) then
        call VCPU(tmExit)
        write(6, '(1X, A, F10.4, A)')
     &    'SMF: CPU-time: ', tmExit - tmEntr, ' seconds total'
      end if
C
      return
      end
C}
