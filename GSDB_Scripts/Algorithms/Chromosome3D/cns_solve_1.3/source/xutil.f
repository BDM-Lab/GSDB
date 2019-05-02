C{
      SUBROUTINE E2PIIX(X, RES)
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      DOUBLE PRECISION  X
      DOUBLE COMPLEX    RES
C
C Return exp(2*PI*i*x)
C
C local
      DOUBLE PRECISION  TPX
C parameters
      DOUBLE PRECISION  TWO
      PARAMETER(TWO = 2.0D0)
C begin
      TPX = TWO * PI * X
      RES = DCMPLX(DCOS(TPX), DSIN(TPX))
C
      RETURN
      END
C}
C=======================================================================
C{
      INTEGER FUNCTION LNBLNK(STRING)
      IMPLICIT NONE
C I/O
      CHARACTER STRING*(*)
C
C Return length of string without trailing blanks
C
C begin
      DO LNBLNK = LEN(STRING), 1, -1
        IF (STRING(LNBLNK:LNBLNK) .NE. ' ') RETURN
      END DO
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE TOUPPER(STRING)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER STRING*(*)
C
C Map all characters of string to upper case
C
C local
      INTEGER  I
C begin
      DO I = 1, LEN(STRING)
        STRING(I:I) = CHAR(ASCIIM(ICHAR(STRING(I:I))))
      END DO
C
      RETURN
      END
C}
C=======================================================================
C{
      INTEGER FUNCTION RECHEM(H, K, L)
      IMPLICIT NONE
C I/O
      INTEGER  H, K, L
C
C Return +1 if hkl is in "positive" reciprocal hemisphere (l > 0 ...)
C        -1 if hkl is in "negative" reciprocal hemisphere (l < 0)
C         0 if hkl is 000
C
C begin
      RECHEM =  1
C
      IF (L .GT. 0) RETURN
      IF (L .EQ. 0) THEN
        IF (K .GT. 0) RETURN
        IF (K .EQ. 0) THEN
          IF (H .GT. 0) RETURN
          IF (H .EQ. 0) THEN
            RECHEM = 0
            RETURN
          END IF
        END IF
      END IF
C
      RECHEM = -1
C
      END
C}
C=======================================================================
C{
      INTEGER FUNCTION CMPHKL(H1, K1, L1, H2, K2, L2)
      IMPLICIT NONE
C I/O
      INTEGER  H1, K1, L1, H2, K2, L2
C
C Define asymmetric unit in reciprocal space by deciding which
C hkl is "larger".
C   Return value =  1 => hkl1 is larger
C                =  0 => hkl1 = hkl2
C                = -1 => hkl2 is larger
C local
      LOGICAL  TAKEIT
C
C begin
C Go through the cascade of criteria for the asymmetric unit
      IF ((K1 .LT. 0 .OR. K2 .LT. 0) .AND. K1 .NE. K2) THEN
        TAKEIT = K2 .GT. K1
      ELSE IF (H1 .NE. H2) THEN
        TAKEIT = H2 .GT. H1
      ELSE IF (L1 .NE. L2) THEN
        TAKEIT = L2 .LT. L1
      ELSE
        TAKEIT = K2 .GT. K1
      END IF
C
      IF (TAKEIT) THEN
        CMPHKL =  1
      ELSE IF (H1 .EQ. H2 .AND. K1 .EQ. K2 .AND. L1 .EQ. L2) THEN
        CMPHKL =  0
      ELSE
        CMPHKL = -1
      END IF
C
      END
C}
C=======================================================================
C{
      INTEGER FUNCTION FNDIOBJ(NOBJS, OBJLBL, WTDLBL)
      IMPLICIT NONE
C I/O
      INTEGER        NOBJS
      CHARACTER*(*)  OBJLBL(*)
      CHARACTER*(*)  WTDLBL
C
C Find index of real space object with given label WTDLBL.
C
C local
      INTEGER  IOBJ
C
C begin
      FNDIOBJ = 0
C
      DO IOBJ = 1, NOBJS
        IF (OBJLBL(IOBJ) .EQ. WTDLBL) THEN
          FNDIOBJ = IOBJ
          RETURN
        END IF
      END DO
C
      RETURN
      END
C}
C=======================================================================
C{
      INTEGER FUNCTION F3DIX(I, N)
      IMPLICIT NONE
C I/O
      INTEGER  I(3), N(3)
C
C Return 1D index given 3D indices I(j) and array dimensions N(j)
C
C begin
      F3DIX = (I(3) * N(2) + I(2)) * N(1) + I(1) + 1
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE F3DI(IX, N, I)
      IMPLICIT NONE
C I/O
      INTEGER  IX, N(3), I(3)
C
C Reverse of F3DIX.
C
C begin
      I(1) = IX - 1
      I(3) = I(1)        / (N(1) * N(2))
      I(1) = I(1) - I(3) * (N(1) * N(2))
      I(2) = I(1)        / N(1)
      I(1) = I(1) - I(2) * N(1)
C
      RETURN
      END
C}
C=======================================================================
C{
      INTEGER FUNCTION IMODPOSITIVE(IX, IY)
      IMPLICIT NONE
C I/O
      INTEGER  IX, IY
C
C Example: iModPositive(7, 3) = mod(7, 3) = 1
C          iModPositive(-7, 3) = 3 + mod(-7, 3) = 3 - 1 = 2
C
C local
      INTEGER  IMP
C
C begin
      IMP = IX
C
      IF (IY .GT. 0) THEN
        IMP = MOD(IMP, IY)
        IF (IMP .LT. 0) IMP = IMP + IY
      END IF
C
      IMODPOSITIVE = IMP
C
      RETURN
      END
C}
C=======================================================================
C{
      INTEGER FUNCTION ILCM(A, B)
      IMPLICIT NONE
C I/O
      INTEGER  A, B
C
C Return Least Common Multiple of A and B.
C
C local
      INTEGER  AA, RI, RJ, RK
C
C begin
          AA = A
      IF (AA .EQ. 0) AA = 1
C
          RI = AA
          RJ = B
      IF (RJ .NE. 0) THEN
 10     CONTINUE
          RK = MOD(RI, RJ)
          IF (RK .EQ. 0) THEN
            RI = RJ
            GOTO 20
          END IF
          RI = MOD(RJ, RK)
          IF (RI .EQ. 0) THEN
            RI = RK
            GOTO 20
          END IF
          RJ = MOD(RK, RI)
          IF (RJ .EQ. 0) GOTO 20
        GOTO 10
C
 20     RI = AA / RI * B
      END IF
C
      IF (RI .LT. 0) RI = -RI
C
      ILCM = RI
C
      RETURN
      END
C}
C=======================================================================
C{
      LOGICAL FUNCTION ISUTR(DFRAC)
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION  DFRAC
C
C "Is Unit Translation":
C Return .true. if mod(DFrac, 1) is close to 0, .false. otherwise.
C
C local
      DOUBLE PRECISION  D
C parameters
      DOUBLE PRECISION  ONE, ONEHALF, MAXDELTA
      PARAMETER(ONE = 1.0D0, ONEHALF = 0.5D0, MAXDELTA = 1.0D-6)
C
C begin
      D = MOD(DFRAC, ONE)
C
      IF      (D .LT. -ONEHALF) THEN
        D = D + ONE
      ELSE IF (D .GT.  ONEHALF) THEN
        D = D - ONE
      END IF
C
      IF (ABS(D) .GT. MAXDELTA) THEN
        ISUTR = .FALSE.
      ELSE
        ISUTR = .TRUE.
      END IF
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE INITSYM(SYMMX, ITSYMMX, MSYMMX, NSYMMX, STBF)
      IMPLICIT NONE
C I/O
      INTEGER  MSYMMX, NSYMMX, STBF
      INTEGER  SYMMX(MSYMMX, 3, 4), ITSYMMX(MSYMMX, 3, 3)
C
C Initialize list of symmetry operations.
C
C begin
C Identity symmetry operation
      SYMMX(1, 1, 1) = 1
      SYMMX(1, 1, 2) = 0
      SYMMX(1, 1, 3) = 0
      SYMMX(1, 1, 4) = 0
      SYMMX(1, 2, 1) = 0
      SYMMX(1, 2, 2) = 1
      SYMMX(1, 2, 3) = 0
      SYMMX(1, 2, 4) = 0
      SYMMX(1, 3, 1) = 0
      SYMMX(1, 3, 2) = 0
      SYMMX(1, 3, 3) = 1
      SYMMX(1, 3, 4) = 0
      NSYMMX = 1
C
C Call xrSyPA to set ITSymMx
      CALL XRSYPA('(X,Y,Z)', 7, MSYMMX, NSYMMX, SYMMX, ITSYMMX, STBF)
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE TRCOOR(TRMX, X, Y, Z, XT, YT, ZT)
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION  TRMX(3, 3)
      DOUBLE PRECISION  X, Y, Z
      DOUBLE PRECISION  XT, YT, ZT
C
C (XYZ)T = TRMX*(XYZ)
C
C E.g., convert fractional coordinates to cartesian coordinates
C or vice versa.
C
C begin
      XT =   TRMX(1, 1) * X
     &     + TRMX(1, 2) * Y
     &     + TRMX(1, 3) * Z
      YT =   TRMX(2, 1) * X
     &     + TRMX(2, 2) * Y
     &     + TRMX(2, 3) * Z
      ZT =   TRMX(3, 1) * X
     &     + TRMX(3, 2) * Y
     &     + TRMX(3, 3) * Z
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE TRVEC3(TRMX, V, TV)
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION  TRMX(3, 3), V(3), TV(3)
C
C TV = TRMX * V
C
C E.g., convert fractional coordinates to cartesian coordinates
C or vice versa.
C
C local
      INTEGER  I
C
C begin
      DO I = 1, 3
        TV(I) =   TRMX(I, 1) * V(1)
     &          + TRMX(I, 2) * V(2)
     &          + TRMX(I, 3) * V(3)
      END DO
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE SUM2BMC(NXY, MINX, MAXX, MINY, MAXY,
     &                   SX, SX2, SY, SY2, SXY,
     &                   B, M, C)
      IMPLICIT NONE
C I/O
      INTEGER           NXY
      DOUBLE PRECISION  MINX, MAXX, MINY, MAXY
      DOUBLE PRECISION  SX, SX2, SY, SY2, SXY
      DOUBLE PRECISION  B, M, C
C
C Compute Linear Regression Coefficients b & m and Correlation
C Coefficient c from pre-computed min/maxxy and the sums Sx ... Syx
C
C local
      DOUBLE PRECISION  N, D
      DOUBLE PRECISION  DX, DY
C
C parameter
      DOUBLE PRECISION  SDYN
      PARAMETER(SDYN = 1.D-6)
C
C begin
      B = 0.
      M = 0.
      C = 0.
      IF (NXY .LT. 1) RETURN
C
      IF (MINX .EQ. MAXX) RETURN
      IF (MINY .EQ. MAXY) THEN
        B = MINY
        RETURN
      END IF
C
      N = DFLOAT(NXY)
      DX = MAX(ABS(MINX - SX / N), ABS(MAXX - SX / N))
      DY = MAX(ABS(MINY - SY / N), ABS(MAXY - SY / N))
      IF (DX .EQ. 0.) RETURN
      IF (DY .EQ. 0.) THEN
        B = SY / N
        RETURN
      END IF
C
      IF (DX .LT. DY * SDYN) RETURN
      IF (DY .LT. DX * SDYN) THEN
        B = SY / N
        RETURN
      END IF
C
      D = N * SX2 - SX * SX
      IF (D .NE. 0.) THEN
        B = (SX2 * SY - SX * SXY) / D
        M = (N * SXY - SX * SY) / D
      END IF
C
      D =   (SX2 - SX * SX / N)
     &    * (SY2 - SY * SY / N)
      IF (D .GT. 0.)
     &  C = (SXY - SX * SY / N) / SQRT(D)
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE MAPSTAT(N, W, LBL, FULL)
      IMPLICIT NONE
C I/O
      INTEGER    N
      REAL       W(*)
      CHARACTER  LBL*(*)
      LOGICAL    FULL
C
C Compute and print map statistics.
C
C local
      INTEGER           I
      REAL              WMIN, WMAX
      DOUBLE PRECISION  WAVE, WAVE2, SWGHT, WSIGMA
C
C parameters
      DOUBLE PRECISION  ZERO
      PARAMETER(ZERO = 0.0D0)
C
C begin
      WMIN   = 0.
      WMAX   = 0.
      WAVE   = ZERO
      WSIGMA = ZERO
C
      IF (N .GT. 0) THEN
        WMIN  = W(1)
        WMAX  = W(1)
C
        IF (FULL) THEN
          WAVE  = W(1)
          WAVE2 = W(1)**2
C
          DO I = 2, N
            WMIN  = MIN(WMIN, W(I))
            WMAX  = MAX(WMAX, W(I))
            WAVE  = WAVE  + W(I)
            WAVE2 = WAVE2 + W(I)**2
          END DO
C
          SWGHT = N
          WAVE = WAVE / SWGHT
          WSIGMA = WAVE2 / SWGHT - WAVE**2
          IF (WSIGMA .LT. ZERO) WSIGMA = ZERO
          WSIGMA = SQRT(WSIGMA)
        ELSE
          DO I = 2, N
            WMIN  = MIN(WMIN, W(I))
            WMAX  = MAX(WMAX, W(I))
          END DO
        END IF
      END IF
C
        WRITE(6, '(1X, 2A, G14.6)') LBL, 'min =   ', WMIN
        WRITE(6, '(1X, 2A, G14.6)') LBL, 'max =   ', WMAX
      IF (FULL) THEN
        WRITE(6, '(1X, 2A, G14.6)') LBL, 'ave =   ', WAVE
        WRITE(6, '(1X, 2A, G14.6)') LBL, 'sigma = ', WSIGMA
      END IF
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE ECHOLOGICAL(NAME, VALUE)
      IMPLICIT NONE
C I/O
      CHARACTER  NAME*(*)
      LOGICAL    VALUE
C
C Show name and value of logical variable.
C
C begin
      IF (VALUE) THEN
        WRITE(6,'(1X,2A)') NAME, '=true'
      ELSE
        WRITE(6,'(1X,2A)') NAME, '=false'
      END IF
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE ECHOOBJ(IOBJ, OBJNAM, WDID)
      IMPLICIT NONE
C I/O
      INTEGER        IOBJ
      CHARACTER*(*)  OBJNAM(*), WDID
C
C Show name of real space or reciprocal object. E.g., "FROM=FCALc"
C or "FROM=<undefined>".
C
C externals
      INTEGER   LNBLNK
      EXTERNAL  LNBLNK
C
C begin
      IF (IOBJ .GT. 0) THEN
        WRITE(6,'(1X,2A)') WDID, OBJNAM(IOBJ)(1:LNBLNK(OBJNAM(IOBJ)))
      ELSE
        WRITE(6,'(1X,2A)') WDID, '<undefined>'
      END IF
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE GETIOBJ(IOBJ, OBJNUM, OBJNAM, ADDLNAM,
     &                   WDID, DOMID, ERRID)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER        IOBJ, OBJNUM
      CHARACTER*(*)  OBJNAM(*), ADDLNAM, WDID, DOMID, ERRID
C
C Get NextWD from input and use it to identify a real space or
C reciprocal space object. Return index of object if found, 0
C otherwise.
C
C local
      INTEGER  I
C
C externals
      INTEGER   FNDIOBJ
      EXTERNAL  FNDIOBJ
C
C begin
      CALL NEXTASS(WDID)
      IF (WD(1:4).EQ.'?   ') THEN
        CALL ECHOOBJ(IOBJ, OBJNAM, WDID)
      ELSE
            IOBJ = FNDIOBJ(OBJNUM, OBJNAM, WD(1:WDLEN))
        IF (IOBJ .LE. 0) THEN
          IF (ADDLNAM .EQ. 'HKL') THEN
            DO I = 1, LEN(ADDLNAM)
              IF (WD(1:WDLEN) .EQ. ADDLNAM(I:I)) THEN
                IOBJ = -I
                RETURN
              END IF
            END DO
          END IF
        END IF
C
        IF (IOBJ .EQ. 0) THEN
          WRITE(6, '(7A)')
     &      ' %', ERRID, '-ERR: ',
     &      DOMID, ' space object "',
     &      WD(1:WDLEN), '" not declared.'
          CALL WRNDIE(0, ERRID, 'object undeclared.')
        END IF
      END IF
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE CHKSFTYPE(IOBJ, OBJNAM, SFTYPE, WTDTYPE, ERRID)
      IMPLICIT NONE
C I/O
      INTEGER        IOBJ
      CHARACTER*(*)  OBJNAM(*), SFTYPE(*), WTDTYPE, ERRID
C
C Check if reciprocal space object is of the required
C "Wanted Type" WtdType.
C
C externals
      INTEGER   LNBLNK
      EXTERNAL  LNBLNK
C
C begin
      IF (IOBJ .GT. 0) THEN
        IF (SFTYPE(IOBJ) .NE. WTDTYPE) THEN
          WRITE(6, '(5A)')
     &      ' %', ERRID, '-ERR: reciprocal space object "',
     &      OBJNAM(IOBJ)(1:LNBLNK(OBJNAM(IOBJ))),
     &      '" is of the wrong type.'
          CALL WRNDIE(5, ERRID, 'TYPE MISMATCH.')
          IOBJ = 0
        END IF
      END IF
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE IVECPA(STR, V, MV, NV)
      IMPLICIT NONE
C I/O
      CHARACTER*(*)  STR
      INTEGER        V(*), MV, NV
C
C Parse strings like "(1 2 3 4)", return integers in V, number of
C integers in nV. mV is the maximum number of integers expected.
C nV=0 indicates an error condition. nV=1,2,...mV indicates success.
C
C local
      INTEGER           IV, I, S, E
      LOGICAL           ALLOWCOMMA, EOE, OK
      DOUBLE PRECISION  F
C
C externals
      DOUBLE PRECISION  DECODF
      EXTERNAL          DECODF
C
C begin
      NV = 0
C
      DO IV = 1, MV
        V(IV) = 0
      END DO
C
      I = 1
      OK = .FALSE.
C
      DO WHILE (I .LE. LEN(STR) .AND. .NOT. OK)
        IF      (STR(I:I) .EQ. '(') THEN
          OK = .TRUE.
        ELSE IF (STR(I:I) .NE. ' ') THEN
          RETURN
        END IF
        I = I + 1
      END DO
C
      IF (.NOT. OK) RETURN
C
      IV = 0
      ALLOWCOMMA = .FALSE.
      EOE = .FALSE.
C
      DO WHILE (I .LE. LEN(STR) .AND. OK)
        IF (EOE) THEN
          IF (STR(I:I) .NE. ' ') OK = .FALSE.
          I = I + 1
        ELSE IF (STR(I:I) .EQ. ')') THEN
          EOE = .TRUE.
          I = I + 1
        ELSE IF (STR(I:I) .EQ. ',') THEN
          IF (.NOT. ALLOWCOMMA) OK = .FALSE.
          ALLOWCOMMA = .FALSE.
          I = I + 1
        ELSE IF (STR(I:I) .EQ. ' ') THEN
          I = I + 1
        ELSE IF (IV .GE. MV) THEN
          OK = .FALSE.
        ELSE
          S = I
          E = 0
          DO WHILE (I .LE. LEN(STR) .AND. E .EQ. 0)
            IF (     STR(I:I) .EQ. ')'
     &          .OR. STR(I:I) .EQ. ','
     &          .OR. STR(I:I) .EQ. ' ') THEN
              E = I - 1
            ELSE
              I = I + 1
            END IF
          END DO
          IF (E .LT. S) THEN
            OK = .FALSE.
          ELSE
            F = DECODF(STR(S:E), E-S+1, OK)
            IF (OK) THEN
                IV = IV + 1
              V(IV) = NINT(F)
              IF (ABS(F - DFLOAT(V(IV))) .GT. 1.0D-6) OK = .FALSE.
            END IF
          END IF
          ALLOWCOMMA = .TRUE.
        END IF
      END DO
C
      IF (EOE .AND. OK) NV = IV
C
      RETURN
      END
C}
C======================================================================
      SUBROUTINE FILLC4(A,N,VALUE)
C
C Fills A with value, for C4 arrays
C
      IMPLICIT NONE
C
C input/ output
             COMPLEX A(*), VALUE
      INTEGER N
C
C local
      INTEGER I
C
C begin
C
      IF (N.GT.0) THEN
        DO I=1,N
          A(I)=VALUE
        END DO
      END IF
C
      RETURN
      END
C=======================================================================
C{
      subroutine SetI4Elem(Array, iArray, Val)
      IMPLICIT NONE
C I/O
      integer  Array(*), Val
      integer  iArray
C
C Assign a value to the iArray'th element of the integer Array
C Intended use: call SetI4Elem(heap(pointer), i, h)
C
C begin
      Array(iArray) = Val
C
      return
      end
C}
C=======================================================================
C{
      subroutine SetC8Elem(Array, iArray, Val)
      IMPLICIT NONE
C I/O
      double complex  Array(*), Val
      integer         iArray
C
C Assign a value to the iArray'th element of the double complex Array
C Intended use: call SetC8Elem(heap(pointer), i, f)
C
C begin
      Array(iArray) = Val
C
      return
      end
C}
C=======================================================================
C{
      subroutine AddC8Elem(Array, iArray, Val)
      IMPLICIT NONE
C I/O
      double complex  Array(*), Val
      integer         iArray
C
C Add a value to the iArray'th element of the double complex Array
C Intended use: call AddC8Elem(heap(pointer), i, f)
C
C begin
      Array(iArray) = Array(iArray) + Val
C
      return
      end
C}
C=======================================================================
C{
      subroutine GetI4Elem(Array, iArray, Val)
      IMPLICIT NONE
C I/O
      integer  Array(*), Val
      integer  iArray
C
C Retrieve the value to the iArray'th element of the integer Array
C Intended use: call GetI4Elem(heap(pointer), i, h)
C
C begin
      Val = Array(iArray)
C
      return
      end
C}
C=======================================================================
C{
      subroutine GetC8Elem(Array, iArray, Val)
      IMPLICIT NONE
C I/O
      double complex  Array(*), Val
      integer         iArray
C
C Retrieve the value to the iArray'th element of the double complex Array
C Intended use: call GetC8Elem(heap(pointer), i, f)
C
C begin
      Val = Array(iArray)
C
      return
      end
C}
C=======================================================================
C{
      subroutine ShUpI4(Array, ic, nc)
      IMPLICIT NONE
C I/O
      integer  Array(*)
      integer  ic, nc
C
C "Shift up" elements in a sorted integer Array to make space for a
C new element.
C nc is the number of elements in the Array including the new element.
C ic is the position for the new element.
C
C local
      integer  jr
C
C begin
      do jr = nc, ic + 1, -1
        Array(jr) = Array(jr - 1)
      end do
C
      return
      end
C}
C=======================================================================
C{
      subroutine ShUpC8(Array, ic, nc)
      IMPLICIT NONE
C I/O
      double complex  Array(*)
      integer         ic, nc
C
C "Shift up" elements in a (sorted) double complex Array to make space for a
C new element.
C nc is the number of elements in the Array including the new element.
C ic is the position for the new element.
C
C local
      integer  jr
C
C begin
      do jr = nc, ic + 1, -1
        Array(jr) = Array(jr - 1)
      end do
C
      return
      end
C}
C=======================================================================
