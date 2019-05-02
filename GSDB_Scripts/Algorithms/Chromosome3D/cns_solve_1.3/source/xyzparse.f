C=======================================================================
C{
      subroutine Discretize(fVal, iVal, Fac, Err)
      IMPLICIT NONE
C I/O
      double precision  fVal
      integer           iVal
      integer           Fac
      logical           Err
C
C Convert the floating point value fVal to the
C rational value iVal / Fac.
C If the nearest rational value deviates too much from fVal,
C Err is set to true.
C
C local
      double precision  fV
C
C begin
      Err = .true.
C
      if (Fac .eq. 0) return
C
          fV = fVal * dfloat(Fac)
      if (fV .lt. dfloat(0)) then
        iVal = int(fV - .5)
      else
        iVal = int(fV + .5)
      end if
C
          fV = fV - iVal
      if (fV .lt. dfloat(0)) fV = -fV
      if (fV .gt. .0001) return
C
      Err = .false.
C
      return
      end
C}
C=======================================================================
C{
      logical function IsDigit(str)
      IMPLICIT NONE
C I/O
      character  str*(1)
C
C Return true if the character in str is a digit, false otherwise.
C
C begin
      IsDigit = .true.
      if (str .eq. '0') return
      if (str .eq. '1') return
      if (str .eq. '2') return
      if (str .eq. '3') return
      if (str .eq. '4') return
      if (str .eq. '5') return
      if (str .eq. '6') return
      if (str .eq. '7') return
      if (str .eq. '8') return
      if (str .eq. '9') return
      IsDigit = .false.
      return
      end
C}
C=======================================================================
C{
      subroutine StrToVal(str, strlen, val, strused)
      IMPLICIT NONE
C I/O
      character         str*(*)
      integer           strlen
      double precision  val
      integer           strused
C
C Decode the floating point number at the beginning of str.
C strused is negative if an error occurs. In this case,
C str(-strused:-strused) is the position where this subroutine
C generated the error.
C If strused is positive, val is the decoded floating point number,
C and str(strused:strused) is the last character used by this
C subroutine. For example, if str is '1.2E10*X', val will be 1.2E10
C and strused will be 6.
C
C local
      integer    iM, iE, i, j
      logical    MDot, MDgt, EDgt
      logical    contloop
      character  buf*11, fmt*16
C
C externals
      logical   IsDigit
      external  IsDigit
C
C begin
      strlen = min(strlen, len(str))
C
C skip white space
      i = 1
      contloop = .true.
      do while (contloop .and. i .le. strlen)
        if (str(i:i) .eq. ' ') then
          i = i + 1
        else
          contloop = .false.
        end if
      end do
C
      MDot = .false.
      MDgt = .false.
      EDgt = .false.
C
C Scan mantissa
      iM = i
C
      if (i .le. strlen) then
        if (str(i:i) .eq. '+' .or. str(i:i) .eq. '-') i = i + 1
      end if
      contloop = .true.
      do while (contloop .and. i .le. strlen)
        if (str(i:i) .eq. '.') then
          if (MDot) then
            contloop = .false.
          else
            MDot = .true.
          end if
        else if (IsDigit(str(i:i))) then
          MDgt = .true.
        else
          contloop = .false.
        end if
        if (contloop) i = i + 1
      end do
C
      iE = i
C
      if (i .le. strlen) then
        if (     str(i:i) .eq. 'E' .or. str(i:i) .eq. 'e'
     &      .or. str(i:i) .eq. 'D' .or. str(i:i) .eq. 'd') then
          i = i + 1
C Scan exponent
          if (i .le. strlen) then
            if (str(i:i) .eq. '+' .or. str(i:i) .eq. '-') i = i + 1
          end if
          contloop = .true.
          do while (contloop .and. i .le. strlen)
            if (IsDigit(str(i:i))) then
              EDgt = .true.
              i = i + 1
            else
              contloop = .false.
            end if
          end do
        end if
      end if
C
      val = dfloat(0)
      strused = i - 1
C
      if (MDgt) then
        if (.not. EDgt) i = iE
        write(buf, '(I11)') i - iM
        j = 1
        contloop = .true.
        do while (contloop)
          if (buf(j:j) .eq. ' ') then
            j = j + 1
          else
            contloop = .false.
          end if
        end do
        write(fmt, '(3A)') '(F', buf(j:), '.0)'
        read(str(iM:i-1), fmt, end = 99, err=99) val
      end if
C
      return
C
 99   continue
      call WrnDie(-5, 'StrToVal', 'Fatal Coding Error')
      call Die
C
      return
      end
C}
C=======================================================================
C{
      subroutine xyzparse(str, strlen,
     &                    RTMx, FacR, FacT, strused)
      IMPLICIT NONE
C I/O
      character  str*(*)
      integer    strlen
      integer    RTMx(0:2, 0:3), FacR, FacT
      integer    strused
C
C Versatile parser for strings encoding symmetry operations, Harker
C operators, transformation matrices or other operators with a 3 by 3
C rotational part and translational part with 3 components.
C Examples for str:
C   -y,x-y,z+1/3
C   2*x+1/4,y-z+1/4,y+z+3/4
C   2/3*x-1/3*y-1/3*z, 1/3*x+1/3*y-2/3*z, 1/3*x+1/3*y+1/3*z
C   0.25x, 1/4y, 4z
C The rotational part of the operator is stored in the integer matrix
C RTMx(0:2, 0:2), the translational part is stored in RTMx(0:2, 3).
C Floating point values are converted to rational numbers using the
C Discretize subroutine (see there). The values for the rotational part
C are RTMx(0:2, 0:2) / FacR. The values for the translational part are
C RTMx(0:2, 3) / FacT.
C strused is negative if an error occurs. In this case,
C str(-strused:-strused) is the position where this subroutine
C generated the error.
C A positive strused signals success. (In this case strused will always
C be equal to strlen.)
C
C local
      integer           istr, i, j
      integer           Row, Column, Sign, Mult
      integer           Mx(0:2, 0:3)
      double precision  ValR(0:2), ValT, Value, V
      logical           P_Add, P_Mult, P_Value, P_XYZ, P_Comma
      logical           contloop, Err
      character  c*(1)
C
C externals
      logical   IsDigit
      external  IsDigit
C
C begin
      strlen = min(strlen, len(str))
C
      do j = 0, 3
      do i = 0, 2
        RTMx(i, j) = 0
          Mx(i, j) = 0
      end do
      end do
C
      Row    = 0
      Column = -1
      Sign   = 1
      Mult   = 0
      do i = 0, 2
        ValR(i) = dfloat(0)
      end do
      ValT   = dfloat(0)
      Value  = dfloat(0)
      P_Add   = .true.
      P_Mult  = .false.
      P_Value = .true.
      P_XYZ   = .true.
      P_Comma = .false.
C
      istr = 1
      contloop = .true.
      do while (contloop)
        c = ' '
        if (istr .gt. strlen) then
          contloop = .false.
        else
          c = str(istr:istr)
        end if
        if (.not. contloop .or. .not. (c .eq. ' ' .or. c .eq. '_')) then
          if (c .eq. '+' .or. c .eq. '-') then
            if (.not. P_Add) goto 99
            if (c .eq. '+') then
              Sign =  1
            else
              Sign = -1
            end if
            if (Column .ge. 0) then
              ValR(Column) = ValR(Column) + Value
            else
              ValT         = ValT         + Value
            end if
            Value = dfloat(0)
            Column = -1
            Mult = 0
            P_Add   = .false.
            P_Mult  = .false.
            P_Value = .true.
            P_XYZ   = .true.
            P_Comma = .false.
          else if (c .eq. '*') then
            if (.not. P_Mult) goto 99
            Mult = 1
            P_Add   = .false.
            P_Mult  = .false.
            P_Value = .true.
            P_XYZ   = .true.
            P_Comma = .false.
          else if (c .eq. '/' .or. c .eq. ':') then
            if (.not. P_Mult) goto 99
            Mult = -1
            P_Add   = .false.
            P_Mult  = .false.
            P_Value = .true.
            P_XYZ   = .false.
            P_Comma = .false.
          else if (c .eq. '.' .or. IsDigit(c)) then
            if (.not. P_Value) goto 99
            call StrToVal(str(istr:strlen), strlen - istr + 1,
     &                    V, strused)
            istr = istr + strused - 1
            if (Sign .eq. -1) then
              V = -V
              Sign = 1
            end if
            if      (Mult .eq.  1) then
              Value = Value * V
            else if (Mult .eq. -1) then
              if      (V     .ne. dfloat(0)) then
                Value = Value / V
              else if (Value .ne. dfloat(0)) then
                goto 99
              end if
            else
              Value = V
            end if
            P_Add   = .true.
            P_Mult  = .true.
            P_Value = .false.
            P_XYZ   = .true.
            P_Comma = .true.
          else if (     c .eq. 'X' .or. c .eq. 'x'
     &             .or. c .eq. 'Y' .or. c .eq. 'y'
     &             .or. c .eq. 'Z' .or. c .eq. 'z') then
            if (.not. P_XYZ) goto 99
            if      (c .eq. 'X' .or. c .eq. 'x') then
              Column = 0
            else if (c .eq. 'Y' .or. c .eq. 'y') then
              Column = 1
            else
              Column = 2
            end if
            if (Value .eq. dfloat(0)) then
              Value = Sign
              Sign = 1
            end if
            P_Add   = .true.
            P_Mult  = .true.
            P_Value = .false.
            P_XYZ   = .false.
            P_Comma = .true.
          else if (.not. contloop .or. c .eq. ',' .or. c .eq. ';') then
            if (Row .eq. 2 .and. contloop) goto 99
            if (.not. P_Comma) goto 99
            if (Column .ge. 0) then
              ValR(Column) = ValR(Column) + Value
            else
              ValT         = ValT         + Value
            end if
            do i = 0, 2
              call Discretize(ValR(i), Mx(Row, i), FacR, Err)
              if (Err) goto 99
            end do
            call Discretize(ValT, Mx(Row, 3), FacT, Err)
            if (Err) goto 99
            Row = Row + 1
            Column = -1
            Sign   = 1
            Mult   = 0
            do i = 0, 2
              ValR(i) = dfloat(0)
            end do
            ValT   = dfloat(0)
            Value  = dfloat(0)
            P_Add   = .true.
            P_Mult  = .false.
            P_Value = .true.
            P_XYZ   = .true.
            P_Comma = .false.
          else
            goto 99
          end if
        end if
        if (contloop) istr = istr + 1
      end do
C
      if (Row .ne. 3) goto 99
C
      do j = 0, 3
      do i = 0, 2
        RTMx(i, j) = Mx(i, j)
      end do
      end do
C
      strused =  istr
      return
C
 99   continue
      strused = -istr
C
      return
      end
C}
