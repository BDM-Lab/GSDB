      SUBROUTINE COPYR4(A,COPY,N)
C
C  Copies array A into COPY, applies to REAL arrays
      REAL A(*), COPY(*)
      INTEGER N, I
      IF (N.GT.0) THEN
      DO I=1,N
      COPY(I)=A(I)
      END DO
      END IF
      RETURN
      END
C=======================================================================
C{
      subroutine ClrYard(iYard)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mapyard.inc'
      integer  iYard
C
C local
      integer  PRn, i
C
C begin
      if (LbMapYard(iYard) .ne. ' ') then
        LbMapYard(iYard) = ' '
C
        PRn = 1
C
        do i = 1, 3
          PRn = PRn * RnMapYard(i, iYard)
                      RnMapYard(i, iYard) = 0
        end do
C
        call FreHp(PtMapYard(iYard), ireal4(PRn))
        PtMapYard(iYard) = 0
      end if
C
      return
      end
C}
C=======================================================================
C{
      subroutine IniMYCommon(ColdStart)
      implicit none
C I/O
      INCLUDE 'mapyard.inc'
      logical  ColdStart
C
C Initialize mapyard.inc common blocks.
C
C local
      integer  iYard, i
C
C begin
      do iYard = 1, mMapYard
        if (ColdStart) then
          LbMapYard(iYard) = ' '
          do i = 1, 3
            RnMapYard(i, iYard) = 0
          end do
          PtMapYard(iYard) = 0
        else
          call ClrYard(iYard)
        end if
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine mapypa(na, nb, nc, nrRho, niRho, QHerm,
     &                  xRhoNum, xRhoNam, hprrho, hpirho)
      implicit none
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'mapyard.inc'
      integer        na, nb, nc, nrRho, niRho
      logical        QHerm
      integer        xRhoNum
      character*(*)  xRhoNam(*)
      integer        hprrho(*), hpirho(*)
C
C This is the parser for mapyard
C
C local
      integer  iRSObj, iYard, i
      logical  Break
C
C parameters
      real      SZero
      parameter(SZero=0.0)
C
C external
      integer   lnblnk, FndiObj
      external  lnblnk, FndiObj
C
C begin
      iRSObj = 0
C
C parsing
      CALL PUSEND('MAPYard>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('MAPYard>')
      CALL MISCOM('MAPYard>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-mapyard')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'STOR') THEN
        call GetiObj(iRSObj, xRhoNum, xRhoNam, ' ',
     &               'STORe(real-space-object)=', 'real', 'MAPYard')
        if (iRSObj .gt. 0) then
          if (hprrho(iRSObj) .eq. 0) then
            write(6, '(3A)') ' %MAPYARD-ERR: Real space object "',
     &        xrhonam(iRSObj)(1:lnblnk(xrhonam(iRSObj))),
     &        '" not defined.'
            call WrnDie(0, 'MAPYARD', 'object undefined.')
          else
            CALL NEXTWD('STORe(yard-label)=')
            iYard = FndiObj(mMapYard, LbMapYard, WD(1:WDLEN))
            if (iYard .gt. 0) then
              call ClrYard(iYard)
              if (WrnLev .ge. 5)
     &          write(6, '(1X, A)') 'MAPYARD: Old map was removed.'
            else
              Break = .false.
              iYard = 1
              do while (.not. Break .and. iYard .le. mMapYard)
                if (PtMapYard(iYard) .eq. 0) then
                  Break = .true.
                else
                  iYard = iYard + 1
                end if
              end do
              if (.not. Break) then
                call WrnDie(0, 'MAPYard',
     &            'MapYard is full, cannot store new map.')
                iYard = 0
              end if
            end if
            if (iYard .gt. 0) then
              LbMapYard(iYard) = WD(1:WDLEN)
              RnMapYard(1, iYard) = na
              RnMapYard(2, iYard) = nb
              RnMapYard(3, iYard) = nc
              PtMapYard(iYard) = AllHp(ireal4(nrRho))
              call CopyR4(heap(hprrho(iRSObj)),
     &                    heap(PtMapYard(iYard)), nrRho)
            end if
          end if
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'REST') THEN
        CALL NEXTWD('RESTore(yard-label)=')
        iYard = FndiObj(mMapYard, LbMapYard, WD(1:WDLEN))
        if (iYard .lt. 1) then
          call WrnDie(0, 'MAPYard', 'Unknown map label.')
        else
          call GetiObj(iRSObj, xRhoNum, xRhoNam, ' ',
     &                 'RESTore(real-space-object)=',
     &                 'real', 'MAPYard')
          if (iRSObj .gt. 0) then
            if (     RnMapYard(1, iYard) .ne. na
     &          .or. RnMapYard(2, iYard) .ne. nb
     &          .or. RnMapYard(3, iYard) .ne. nc) then
              call WrnDie(0, 'MAPYARD',
     &          'map dimensions do not match.')
            else
              if (hprrho(iRSObj) .eq. 0) then
                call xMapAl(hprrho(iRSObj), hpirho(iRSObj), QHerm,
     &                      nrRho, niRho)
              else if (.not. QHerm) then
                call FillR4(heap(hpirho(iRSObj)), niRho, SZERO)
              end if
              call CopyR4(heap(PtMapYard(iYard)),
     &                    heap(hprrho(iRSObj)), nrRho)
            end if
          end if
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'CLEA') THEN
        CALL NEXTWD('CLEAr=')
        iYard = FndiObj(mMapYard, LbMapYard, WD(1:WDLEN))
        if (iYard .gt. 0) then
          call ClrYard(iYard)
          if (WrnLev .ge. 10)
     &      write(6, '(1X, A)') 'MAPYARD: Map was removed.'
        else if (WD(1:4) .eq. 'ALL') then
          do iYard = 1, mMapYard
            call ClrYard(iYard)
          end do
          if (WrnLev .ge. 10)
     &      write(6, '(1X, A)') 'MAPYARD: All maps were removed.'
        else
          call WrnDie(0, 'MAPYard', 'Unknown map label.')
        end if
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
        WRITE(6,'(2A)') ' -----------mapyard-inventory-',
     &   '--------------------------------------------------'
        WRITE(6,'(1X, A)') 'Label              na   nb   nc'
        do iYard = 1, mMapYard
          if (PtMapYard(iYard) .ne. 0) then
            write(6, '(1X, A, 3(1X, I4))')
     &        LbMapYard(iYard),
     &        (RnMapYard(i, iYard), i = 1, 3)
          end if
        end do
        WRITE(6,'(2A)') ' -----------------------------',
     &   '--------------------------------------------------'
C=====================================================================
      ELSE
        CALL CHKEND('MAPYard>',DONE)
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
      subroutine mapyard(Mode, na, nb, nc, nrRho, niRho, QHerm,
     &                   xRhoNum, xRhoNam, hprrho, hpirho,
     &                   xrNsym)
      implicit none
C I/O
      integer        Mode
      integer        na, nb, nc, nrRho, niRho
      logical        QHerm
      integer        xRhoNum
      character*(*)  xRhoNam(*)
      integer        hprrho(*), hpirho(*)
      integer        xrNsym
C
C XRAY MAPYard main procedure.
C
C begin
      if (Mode .eq. 0) then
C Initialize common block and return
        call IniMYCommon(.false.)
        return
      end if
C
      if (Mode .lt. 0) then
C Deallocate maps and return
        call IniMYCommon(.false.)
        return
      end if
C
      if (nrRho .ne. na * nb * nc) then
        call WrnDie(-5, 'MAPYARD', 'Fatal Coding Error.')
        call Die
      end if
C
      if (xrNsym .ne. 1) then
        call WrnDie(0, 'MAPYARD',
     &    'This module can only be used in P1.')
        return
      end if
C
C Call the parser
      call mapypa(na, nb, nc, nrRho, niRho, QHerm,
     &            xRhoNum, xRhoNam, hprrho, hpirho)
C
      return
      end
C}
