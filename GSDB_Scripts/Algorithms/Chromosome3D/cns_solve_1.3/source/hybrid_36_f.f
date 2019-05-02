C Fortran 77 port of the hy36encode() and hy36decode() functions in the
C hybrid_36.py Python prototype/reference implementation.
C See the Python script for more information.
C
C This file has no external dependencies.
C
C This file is unrestricted Open Source (cctbx.sf.net).
C Please send corrections and enhancements to cctbx@cci.lbl.gov .
C
C See also: http://cci.lbl.gov/hybrid_36/
C
C Ralf W. Grosse-Kunstleve, Feb 2007.

C This subroutine is an implementation detail.
      subroutine encode_pure(
     &  digit_set,
     &  width,
     &  value,
     &  result)
      implicit none
C Input
      character digit_set*(*)
      integer width
      integer value
C Output
      character result*(*)
C Local
      character buf*16
      integer digits_size
      integer i, j, k, rest, val
C
      digits_size = len(digit_set)
      val = value
      i = 0
      j = 0
      if (val .lt. 0) then
        j = 1
        val = -val
      endif
    1 continue
        rest = val / digits_size
        k = val - rest * digits_size + 1
        i = i + 1
        buf(i:i) = digit_set(k:k)
        if (rest .ne. 0) then
          val = rest
          goto 1
        endif
      if (j .ne. 0) then
        i = i + 1
        buf(i:i) = '-'
      endif
      if (i .ge. width) then
        k = 1
      else
        k = width - i
        result(1:k) = ' '
        k = k + 1
      endif
    2 continue
        result(k:k) = buf(i:i)
        i = i - 1
        k = k + 1
        if (i .gt. 0) goto 2
      return
      end

C This subroutine is an implementation detail.
      subroutine decode_pure(
     &  digits_values,
     &  digits_size,
     &  s,
     &  result,
     &  errmsg,
     &  errmsg_len)
      implicit none
C Input
      integer digits_values(0:*)
      integer digits_size
      character s*(*)
C Output
      integer result
      character errmsg*(*)
      integer errmsg_len
C Local
      logical first_call
      save first_call
      data first_call /.true./
      integer iblank
      save iblank
      integer iminus
      save iminus
      integer si, dv
      logical have_minus
      logical have_non_blank
      integer value
      integer i
C
      if (first_call) then
        first_call = .false.
        iblank = ichar(' ')
        iminus = ichar('-')
      endif
      have_minus = .false.
      have_non_blank = .false.
      value = 0
      do 1 i=1,len(s)
        si = ichar(s(i:i))
        if (si .lt. 0 .or. si .gt. 127) goto 2
        if (si .eq. iblank) then
          if (.not. have_non_blank) goto 1
          value = value * digits_size
        else if (si .eq. iminus) then
          if (have_non_blank) goto 2
          have_non_blank = .true.
          have_minus = .true.
          goto 1
        else
          have_non_blank = .true.
          dv = digits_values(si)
          if (dv .lt. 0 .or. dv .ge. digits_size) goto 2
          value = value * digits_size
          value = value + dv
        endif
    1 continue
      if (have_minus) value = -value
      result = value
      errmsg_len = 0
      return
    2 result = 0
      errmsg = 'invalid number literal.'
      errmsg_len = 23
      return
      end

C
C hybrid-36 encoder: converts integer value to string result
C
C   width: must be 4 (e.g. for residue sequence numbers)
C               or 5 (e.g. for atom serial numbers)
C
C   value: integer value to be converted
C
C   result: string holding the result (len(result) must be >= width)
C           result(width+1:len(result)) is NOT modified in order
C           to maximize runtime performance!
C
C   errmsg: string holding error message, if any
C           len(errmsg) should be >= 80
C           errmsg is not modified if errmsg_len below is 0
C           DO NOT use if (errmsg .eq. ' ') to check for errors!
C           It is not supported because it would be too slow.
C
C   errmsg_len: length of error message, or 0 on success
C               use if (errmsg_len .ne. 0) to check for errors
C
      subroutine hy36encode(
     &  width,
     &  value,
     &  result,
     &  errmsg,
     &  errmsg_len)
      implicit none
C Input
      integer width
      integer value
C Output
      character result*(*)
      character errmsg*(*)
      integer errmsg_len
C Local
      character digits_upper*36
      character digits_lower*36
      data digits_upper /'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data digits_lower /'0123456789abcdefghijklmnopqrstuvwxyz'/
      save digits_upper
      save digits_lower
      integer i
C
      i = value
      if (width .eq. 4) then
        if (i .ge. -999) then
          if (i .lt. 10000) then
            call encode_pure(digits_upper(1:10), 4, i, result)
            errmsg_len = 0
            return
          endif
          i = i - 10000
          if (i .lt. 1213056) then
            i = i + 466560
            call encode_pure(digits_upper, 0, i, result)
            errmsg_len = 0
            return
          endif
          i = i - 1213056
          if (i .lt. 1213056) then
            i = i + 466560
            call encode_pure(digits_lower, 0, i, result)
            errmsg_len = 0
            return
          endif
        endif
      else if (width .eq. 5) then
        if (i .ge. -9999) then
          if (i .lt. 100000) then
            call encode_pure(digits_upper(1:10), 5, i, result)
            errmsg_len = 0
            return
          endif
          i = i - 100000
          if (i .lt. 43670016) then
            i = i + 16796160
            call encode_pure(digits_upper, 0, i, result)
            errmsg_len = 0
            return
          endif
          i = i - 43670016
          if (i .lt. 43670016) then
            i = i + 16796160
            call encode_pure(digits_lower, 0, i, result)
            errmsg_len = 0
            return
          endif
        endif
      else
        errmsg = 'unsupported width.'
        errmsg_len = 18
        goto 1
      endif
      errmsg = 'value out of range.'
      errmsg_len = 19
    1 do 2, i=1,min(max(1,width),len(result))
        result(i:i) = '*'
    2 continue
      return
      end

C
C hybrid-36 decoder: converts string s to integer result
C
C   width: must be 4 (e.g. for residue sequence numbers)
C               or 5 (e.g. for atom serial numbers)
C
C   s: string to be converted
C        len(s) must be equal to width, or an error message is
C        returned otherwise
C
C   result: integer holding the conversion result
C
C   errmsg, errmsg_len: see hy36encode documentation above
C
      subroutine hy36decode(
     &  width,
     &  s,
     &  result,
     &  errmsg,
     &  errmsg_len)
      implicit none
C Input
      integer width
      character s*(*)
C Output
      integer result
      character errmsg*(*)
      integer errmsg_len
C Local
      character digits_upper*36
      character digits_lower*36
      data digits_upper /'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data digits_lower /'0123456789abcdefghijklmnopqrstuvwxyz'/
      save digits_upper
      save digits_lower
      logical first_call
      integer digits_values_upper(0:127)
      integer digits_values_lower(0:127)
      save first_call
      save digits_values_upper
      save digits_values_lower
      data first_call /.true./
      data digits_values_upper /128*-1/
      data digits_values_lower /128*-1/
      character ie_range*57
      save ie_range
      data ie_range /
     &  'internal error hy36decode: integer value out of range.'/
      integer i, di
C
      if (first_call) then
        first_call = .false.
        do 1, i=1,36
          di = ichar(digits_upper(i:i))
          if (di .lt. 0 .or. di .gt. 127) then
            result = 0
            errmsg = ie_range
            errmsg_len = len(ie_range)
            return
          endif
          digits_values_upper(di) = i-1
    1   continue
        do 2, i=1,36
          di = ichar(digits_lower(i:i))
          if (di .lt. 0 .or. di .gt. 127) then
            result = 0
            errmsg = ie_range
            errmsg_len = len(ie_range)
            return
          endif
          digits_values_lower(di) = i-1
    2   continue
      endif
      if (len(s) .eq. width) then
        di = ichar(s(1:1))
        if (di .ge. 0 .and. di .le. 127) then
          if (digits_values_upper(di) .ge. 10) then
            call decode_pure(
     &        digits_values_upper, 36, s, result, errmsg, errmsg_len)
            if (errmsg_len .ne. 0) return
            if (width .eq. 4) then
C                             - 10*36**(width-1) + 10**width
              result = result - 456560
              errmsg_len = 0
            else if (width .eq. 5) then
              result = result - 16696160
              errmsg_len = 0
            else
              goto 3
            endif
            return
          else if (digits_values_lower(di) .ge. 10) then
            call decode_pure(
     &        digits_values_lower, 36, s, result, errmsg, errmsg_len)
            if (errmsg_len .ne. 0) return
            if (width .eq. 4) then
C                             + 16*36**(width-1) + 10**width
              result = result + 756496
              errmsg_len = 0
            else if (width .eq. 5) then
              result = result + 26973856
              errmsg_len = 0
            else
              goto 3
            endif
            return
          else
            call decode_pure(
     &        digits_values_upper, 10, s, result, errmsg, errmsg_len)
            if (.not. (width .eq. 4 .or. width .eq. 5)) goto 3
            return
          endif
        endif
      endif
    3 result = 0
      errmsg = 'unsupported width.'
      errmsg_len = 18
      return
      end
C