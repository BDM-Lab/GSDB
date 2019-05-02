      SUBROUTINE CHKNUM(WD,WDLEN,FOUND)
C
C checks whether WD is a number (real or integer)
C A number is defined as having a digit (1,2,3,4,5,6,7,8,9,0) as
C the first character or having (.1,.2,.3,.4,.5,.6,.7,.8,.9,.0) as
C the first two characters.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      CHARACTER*(*) WD
      INTEGER WDLEN
      LOGICAL FOUND
C local
      CHARACTER*1 TEMP
C begin
      FOUND=.FALSE.
      IF (WD(1:1).EQ.'.'.AND.WDLEN.GE.2) THEN
      TEMP=WD(2:2)
      ELSE
      TEMP=WD(1:1)
      END IF
      IF (INDEX('0123456789',TEMP).GT.0) FOUND=.TRUE.
      RETURN
      END
C
      SUBROUTINE CHKNUMN(WD,WDLEN,FOUND)
C
C checks whether WD is a number (real or integer) negative sign dealt
C with correctly
C A number is defined as having a digit (1,2,3,4,5,6,7,8,9,0) as
C the first character or having (.1,.2,.3,.4,.5,.6,.7,.8,.9,.0) as
C the first two characters.
C
C Author: Axel T. Brunger and Paul Adams
C ======================================
C
      IMPLICIT NONE
C I/O
      CHARACTER*(*) WD
      INTEGER WDLEN
      LOGICAL FOUND
C local
      CHARACTER*1 TEMP
C begin
      FOUND=.FALSE.
      IF (WD(1:1).EQ.'.'.AND.WDLEN.GE.2) THEN
      TEMP=WD(2:2)
      ELSE IF (WD(1:2).EQ.'-.'.AND.WDLEN.GE.3) THEN
      TEMP=WD(3:3)
      ELSE IF (WD(1:1).EQ.'-'.AND.WDLEN.GE.2) THEN
      TEMP=WD(2:2)
      ELSE
      TEMP=WD(1:1)
      END IF
      IF (INDEX('0123456789',TEMP).GT.0) FOUND=.TRUE.
      RETURN
      END
C
      SUBROUTINE CHKWRD(WD,WDLEN,OK)
C
C checks whether WD begins with a letter and that it
C does not contain any special characters
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      CHARACTER*(*) WD
      INTEGER WDLEN
      LOGICAL OK
C local
      INTEGER I
      CHARACTER*53 LETTER
      CHARACTER*10 NUMBER
C begin
      LETTER='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_'
      NUMBER='1234567890'
      OK=.TRUE.
      IF (INDEX(LETTER,WD(1:1)).EQ.0) OK=.FALSE.
      DO I=2,WDLEN
      IF (INDEX(LETTER,WD(I:I)).EQ.0.AND.INDEX(NUMBER,WD(I:I)).EQ.0)
     &   OK=.FALSE.
      END DO
      RETURN
      END
C
      SUBROUTINE COPYST(ST2,ST2MAX,ST2LEN,ST1,ST1LEN)
C
C copies string ST1 to string ST2, takes care of string lengths and
C string ST2 is set to blanks before copying.
C
C Author: Axel Brunger
C ====================
C
      IMPLICIT NONE
C input/ouput
      CHARACTER*(*) ST1
      INTEGER ST1LEN
      CHARACTER*(*) ST2
      INTEGER ST2MAX, ST2LEN
C begin
      ST2LEN=0
      ST2=' '
      IF (ST1LEN.GT.ST2MAX) THEN
      WRITE(6,'(A)')
     &  ' %COPYST-ERR: ST2MAX too small. Check input file.'
      WRITE(6,'(3A)')
     &  '           Offending string:"',ST1(1:ST1LEN),'"'
      WRITE(6,'(A,I5)')
     &  '                with length=',ST1LEN
      WRITE(6,'(A,I5)')
     &  '           Max allowed length of string=',ST2MAX
      ELSE
      ST2LEN=ST1LEN
      IF (ST1LEN.GT.0) ST2(1:ST1LEN)=ST1(1:ST1LEN)
      END IF
      RETURN
      END
C
      SUBROUTINE ADDST(ST,STMAX,STLEN,ADST,ADLEN)
C
C concatenates string ADST onto ST
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      CHARACTER*(*) ST
      INTEGER STMAX, STLEN
      CHARACTER*(*) ADST
      INTEGER ADLEN
C begin
      IF (ADLEN+STLEN.GT.STMAX) THEN
      WRITE(6,'(A)') ' %ADDST-ERR: STMAX too small. Check code'
      ELSE
      IF (ADLEN.GT.0) ST(STLEN+1:ADLEN+STLEN)=ADST(1:ADLEN)
      STLEN=STLEN+ADLEN
      END IF
      RETURN
      END
C
      SUBROUTINE TRIMM(ST,STLEN)
C routine removes trailing blanks in string ST
C
C Axel T. Brunger
C ===============
C
      IMPLICIT NONE
C input/output
      CHARACTER*(*) ST
      INTEGER STLEN
C begin
      IF (STLEN.GT.0) THEN
      DO WHILE (STLEN.GE.1.AND.ST(STLEN:STLEN).EQ.' ')
      STLEN=STLEN-1
      END DO
      END IF
      RETURN
      END
C
      SUBROUTINE TRIML(ST,STLEN)
C routine removes leading blanks in string ST and fills the
C remainder on the right site with blanks
C
C Axel T. Brunger
C ===============
C
      IMPLICIT NONE
C input/output
      CHARACTER*(*) ST
      INTEGER STLEN
C local
      INTEGER IS, I, STOLD
      CHARACTER*1 TEMP
C begin
      IF (STLEN.GT.0) THEN
      STOLD=STLEN
      IS=1
      DO WHILE (IS.LE.STLEN.AND.ST(IS:IS).EQ.' ')
      IS=IS+1
      END DO
      STLEN=STLEN-IS+1
      IF (STLEN.LT.STOLD) THEN
      IS=IS-1
C
C The following code avoids a self-assignment ( ST(..) = ST(..) )
C and therefore makes the CFT 1.13 compiler on the CRAY happy.
C As soon as this is fixed (in the compiler) this inefficient loop
C should be replaced !
      DO I=1,STLEN
      TEMP=ST(I+IS:I+IS)
      ST(I:I)=TEMP
      END DO
      ST(STLEN+1:STOLD)=' '
      END IF
      END IF
      RETURN
      END
C
      SUBROUTINE SPLITI(STRING,NUM,ICODE,STLEN,OK)
C
C     Split the STRING containing a residue number with optional iCode
C     into integer NUM and letter ICODE, where the integer part may be
C     4-character Hybrid-36 encoded. The last character is accepted as
C     an insertion code if alphabetic and the string is not exactly 4
C     characters with a leading alphabetic character. This allows
C     Hybrid-36 to be used with an iCode, and allows input of normal
C     base-10 integers larger than 9999 in selection expressions.
C
C     The string is assumed to contain only alphanumeric, '-' or ' '
C     characters.  If decoding fails then NUM is set to 0, iCode is set
C     to ' ' and no error message is issued.
C
C MODIFICATION: re-written from general routine to split integer+alpha
C characters to specific routine for resSeq+iCode with Hybrid-36 resSeq.
C Argument changess: renamed ALPHA to ICODE. ICODE is always one
C character, and STLEN is input only.
C
C Joseph M. Krahn
C ===============
C
      IMPLICIT NONE
C input/output
      CHARACTER*(*) STRING
      CHARACTER*1 ICODE
      INTEGER NUM
      LOGICAL OK
C local
      INTEGER LEN, STLEN, DECODI
      CHARACTER*1 ERRMSG
      INTEGER ERRMSG_LEN
C begin
      OK=.TRUE.
      ICODE=' '
      NUM=0
      LEN=STLEN
      CALL TRIMM(STRING,LEN)
      IF (LEN.EQ.0) RETURN
      IF (LGE(STRING(LEN:LEN),'A') .AND. .NOT.
     &    (LGE(STRING(1:1),'A').AND.LEN.EQ.4)) THEN
      ICODE=STRING(LEN:LEN)
      LEN=LEN-1
      ENDIF
      IF (LGE(STRING(1:1),'A').AND.LEN.EQ.4) THEN
      CALL HY36DECODE(4,STRING(1:4),NUM,ERRMSG,ERRMSG_LEN)
      IF (ERRMSG_LEN.NE.0) OK=.FALSE.
      ELSE IF (LEN.GT.0) THEN
      NUM=DECODI(STRING,LEN,OK)
      ELSE
      OK=.FALSE.
      END IF
      RETURN
      END
C
      SUBROUTINE ENCODI(I,ST,STMAX,STLEN)
C
C Converts integer I into string ST.
C Restricted to maximum of 20 digits.
C
C By Axel T. Brunger
C ==================
C
      IMPLICIT NONE
C input/output
      INTEGER I
      CHARACTER*(*) ST
      INTEGER STMAX, STLEN
C local
      CHARACTER*20 TEMP
      INTEGER P
C begin
C
      ST=' '
C
C encode integer I using (right-justified) I20 format.
      WRITE(TEMP,'(I20)') I
C
C remove leading blanks
      P=1
      DO WHILE (P.LE.20.AND.TEMP(P:P).EQ.' ')
      P=P+1
      END DO
      STLEN=20-P+1
      IF (STMAX.LT.STLEN) THEN
      WRITE(6,'(A)') ' %ENCODI-ERR: check STMAX'
      ELSE
      ST(1:STLEN)=TEMP(P:20)
      END IF
      RETURN
      END
C
      SUBROUTINE ENCODF(R,ST,STMAX,STLEN)
C
C This will encode the real number r into the string st and will
C attempt to shorten the encoding as much as possible:
C removal of leading blanks, removal of trailing
C blanks, removal of trailing zeros in the mantissa.
C The conversion format is stored in symbol.inc.
C
C Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'symbol.inc'
      DOUBLE PRECISION R
      CHARACTER*(*) ST
      INTEGER STMAX, STLEN
C local
      INTEGER P, PEXP, STNEW
      LOGICAL CLOOP
      CHARACTER*40 TEMP
C begin
C
      ST=' '
C
C encode real R using (right-justified) CCFORM format.
      IF (R.EQ.0.0D0) THEN
      TEMP='0.'
      ELSE
      WRITE(TEMP(1:CCLENG),CCFORM) R
      END IF
C
C remove leading and trailing blanks
      STLEN=CCLENG
      CALL TRIML(TEMP,STLEN)
      CALL TRIMM(TEMP,STLEN)
C
C remove trailing zeros in the mantissa
      P=STLEN+1
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      P=P-1
      IF (P.EQ.1.OR.TEMP(P:P).EQ.'E') CLOOP=.FALSE.
      END DO
      IF (TEMP(P:P).EQ.'E') THEN
      PEXP=P
      ELSE
      PEXP=STLEN+1
      END IF
      P=PEXP
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      P=P-1
      IF (P.EQ.1.OR.TEMP(P:P).NE.'0') CLOOP=.FALSE.
      END DO
      IF (TEMP(P:P).EQ.'.') P=P-1
      STNEW=STLEN-(PEXP-P-1)
C
      IF (STMAX.LT.STNEW) THEN
      WRITE(6,'(A)') ' %ENCODF-ERR: check STMAX'
      ELSE
      ST(1:P)=TEMP(1:P)
      IF (PEXP.LE.STLEN) ST(P+1:STNEW)=TEMP(PEXP:STLEN)
      STLEN=STNEW
      END IF
C
      RETURN
      END
C
      INTEGER FUNCTION DECODI(ST,STLEN,OK)
C
C Converts the string into an integer, if conversion not
C successful OK=.FALSE.
C Restricted to an input field with a maximum of 20 characters.
C
C Axel T. Brunger
C ===============
C
      IMPLICIT NONE
C input/output
      CHARACTER*(*) ST
      INTEGER STLEN
      LOGICAL OK
C local
      INTEGER I
      CHARACTER*20 TEMP
C
C begin
      IF (STLEN.LE.0) THEN
      OK=.FALSE.
      DECODI=0
      ELSE
      OK=.TRUE.
      TEMP=' '
C
C right-justify the integer in the 20 character string TEMP
      TEMP(20-STLEN+1:20)=ST(1:STLEN)
C
      READ(TEMP,'(I20)',ERR=9999) I
      DECODI=I
C
      GOTO 8888
9999  OK=.FALSE.
      DECODI=0
8888  CONTINUE
C
      END IF
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION DECODF(ST,STLEN,OK)
C
C Decodes the string into a real number and returns its value. It
C uses a F20.0 format.
C
C Axel T. Brunger
C ===============
C
      IMPLICIT NONE
C input/output
      CHARACTER*(*) ST
      INTEGER STLEN
      LOGICAL OK
C local
      DOUBLE PRECISION R
      CHARACTER*20 TEMP
C begin
C
      OK=.TRUE.
      TEMP=' '
      IF (ST(1:STLEN).EQ.'D'.OR.ST(1:STLEN).EQ.'E') GOTO 9999
      IF (ST(1:STLEN).EQ.'d'.OR.ST(1:STLEN).EQ.'e') GOTO 9999
C
C right-justify the integer in the 20 character string TEMP
      TEMP(20-STLEN+1:20)=ST(1:STLEN)
C
      READ(TEMP,'(F20.0)',ERR=9999) R
      DECODF=R
C
      GOTO 8888
9999  OK=.FALSE.
      R=0.0D0
8888  CONTINUE
      RETURN
      END
C
      LOGICAL FUNCTION LTSTEQ(ST1,ST1LEN,ST2,ST2LEN,LEQFG)
C
C Determines whether string st1 is lexicographically less (or equal)
C than string st2, i.e. st1 <= st2. If one string is equal to the
C beginning of the other string it is considered to be less than it.
C If flag LEQFG is false LTSTEQ is flase if strings are equal.
C
C Axel T. Brunger
C ===============
C
      IMPLICIT NONE
C input/output
      CHARACTER*(*) ST1
      INTEGER   ST1LEN
      CHARACTER*(*) ST2
      INTEGER   ST2LEN
      LOGICAL LEQFG
C begin
      IF (LEQFG) THEN
      LTSTEQ=LLE(ST1(1:ST1LEN),ST2(1:ST2LEN))
      ELSE
      LTSTEQ=LLT(ST1(1:ST1LEN),ST2(1:ST2LEN))
      END IF
C
      RETURN
      END
      SUBROUTINE EQSTWC(ST,STLEN,WC,WCLEN,START,DIM,Q)
C
C This logical function matchs the word in ST against a wildcard
C specification.
C
C The wildcard characters are:
C
C  * - matches any string of arbitary length
C  % - matches any single character
C  # - matches any string of digits of arbitrary length
C  + - matches any single digit
C
C Trailing blanks are permitted on either the wildcard or string.
C
C The sequence <backslash><char> forces <char> to be interpreted
C as a literal character. If <backslash> is at the end of a string,
C it is ignored as if it specified a trailing blank.
C
C This algorithm requires recursion to operate. Whenever the wild
C card involves matching any number of characters, a loop must enter
C where the rest of the pattern is matched against succesive
C portions of the string. A stack is used to maintain local
C variables and return points.
C
C
      IMPLICIT NONE
C input/output
      INTEGER DIM
      INTEGER STLEN
      CHARACTER*(*) ST(DIM)
      INTEGER WCLEN
      CHARACTER*(*) WC
      INTEGER START
      LOGICAL Q(*)
C local
      INTEGER RET, IST, IWC, STL, WCL, OLDLST, K, KK
      LOGICAL COND
      INTEGER STKSIZ
      PARAMETER (STKSIZ=99)
      INTEGER STACK(STKSIZ), LSTUSD
      CHARACTER*1 CHAR
      INTEGER :: I, WCLEN_
      CHARACTER*512 :: WC_
      CHARACTER*1 WC_STRING, WC_CHAR, WC_DIGITS, WC_DIGIT, BACKSL
C parameter
      WC_STRING=ACHAR(17)
      WC_CHAR=ACHAR(18)
      WC_DIGITS=ACHAR(19)
      WC_DIGIT=ACHAR(20)
C begin
C Convert wildcard characters to special non-printable character
C constants, and convert <backslash><char> to <char>.
      CALL SETBSL(BACKSL)
      WCLEN_=0
      I=0
      DO WHILE(I<WCLEN)
        I=I+1
        IF (WC(I:I).EQ.BACKSL) THEN
          IF (I.EQ.WCLEN) EXIT
          I=I+1
          CHAR = WC(I:I)
        ELSE IF (WC(I:I).EQ.'*') THEN
          CHAR = WC_STRING
        ELSE IF (WC(I:I).EQ.'%') THEN
          CHAR = WC_CHAR
        ELSE IF (WC(I:I).EQ.'#') THEN
          CHAR = WC_DIGITS
        ELSE IF (WC(I:I).EQ.'+') THEN
          CHAR = WC_DIGIT
        ELSE
          CHAR = WC(I:I)
        END IF
        WCLEN_=WCLEN_+1
        WC_(WCLEN_:WCLEN_) = CHAR
      END DO
! From here on, wildcard patter is in WC_,WCLEN_ instead of WC,WCLEN,
! and WC parameter constants instead of printable wildcard chars.
      LSTUSD=0
C
C if dim greater than 5 make a quick check whether WC_ is a
C wildcard at all
      COND=.TRUE.
      IF (DIM.GT.5) THEN
         COND=(INDEX(WC_(1:WCLEN_),WC_STRING).GT.0.OR.
     @         INDEX(WC_(1:WCLEN_),WC_CHAR).GT.0.OR.
     @         INDEX(WC_(1:WCLEN_),WC_DIGITS).GT.0.OR.
     @         INDEX(WC_(1:WCLEN_),WC_DIGIT).GT.0 )
      END IF
C
C     IF (.NOT.COND) THEN
      IF (COND) GOTO 8765
         DO 17788 K=START,START+DIM-1
            Q(K)=WC_(1:WCLEN_).EQ.ST(K-START+1)(1:STLEN)
17788 CONTINUE
         GOTO 8764
C
C loop over all elements of character array to make wildcard checking
8765     CONTINUE
         WCL = WCLEN_
         CALL TRIMM(WC_,WCL)
         DO 17789 K=START,START+DIM-1
            KK=K-START+1
            IST = 1
            IWC = 1
            STL = STLEN
17790 IF (STL.GE.1.AND.ST(KK)(STL:STL).EQ.' ') THEN
               STL=STL-1
      GOTO 17790
      ENDIF
            OLDLST = LSTUSD
   10       CONTINUE
            IF ( .NOT. (IST .GT. STL .AND. IWC .GT. WCL) ) GOTO 9875
               Q(K) = .TRUE.
               GOTO 9876
9875        IF ( .NOT. (IWC .GT. WCL) ) GOTO 9874
               Q(K) = .FALSE.
               GOTO 9876
9874        IF ( .NOT. (WC_(IWC:IWC) .EQ. WC_STRING ) ) GOTO 9873
               IWC = IWC + 1
               RET=20
               GOTO 55555
   20          IF (Q(K)) GOTO 44444
17791 CONTINUE
                  IST = IST + 1
                  RET=30
                  GOTO 55555
   30             IF (Q(K)) GOTO 44444
      IF(.NOT. (IST .GT. STL)) GOTO 17791
               Q(K) = .FALSE.
               GOTO 9876
9873        IF ( .NOT. (WC_(IWC:IWC) .EQ. WC_DIGITS) ) GOTO 9872
               IWC = IWC + 1
               RET=40
               GOTO 55555
   40          IF (Q(K)) GOTO 44444
17792 CONTINUE
                  CHAR=ST(KK)(IST:IST)
                     IF (INDEX('0123456789',CHAR).EQ.0) THEN
                     Q(K) = .FALSE.
                     GOTO 44444
                  END IF
                  IST = IST + 1
                  RET=50
                  GOTO 55555
   50             IF (Q(K)) GOTO 44444
      IF(.NOT. (IST .GT. STL)) GOTO 17792
               Q(K) = .FALSE.
               GOTO 9876
9872        IF ( .NOT. (IST .GT. STL) ) GOTO 9871
               Q(K) = .FALSE.
               GOTO 9876
9871        IF ( .NOT. (WC_(IWC:IWC) .EQ. WC_CHAR ) ) GOTO 9870
               IWC = IWC + 1
               IST = IST + 1
               GOTO 10
9870        IF ( .NOT. ( WC_(IWC:IWC) .EQ. WC_DIGIT) ) GOTO 9869
C
               CHAR=ST(KK)(IST:IST)
               IF (INDEX('0123456789',CHAR).EQ.0) THEN
                  Q(K) = .FALSE.
                  GOTO 44444
               END IF
               IST = IST + 1
               IWC = IWC + 1
               GOTO 10
9869        IF(.NOT.(ST(KK)(IST:IST).EQ.WC_(IWC:IWC)))GOTO 9867
               IWC = IWC + 1
               IST = IST + 1
               GOTO 10
9867        CONTINUE
               Q(K) = .FALSE.
9876        CONTINUE
            GOTO 44444
C
C ==================================================================
C BEGIN
55555       CONTINUE
            LSTUSD = LSTUSD + 3
            IF (LSTUSD .GT. STKSIZ) THEN
      CALL WRNDIE(-5,'EQSTWC',
     &  ' STKSIZ (routine EQSTWC) exceeded.  --> recompile.')
               RETURN
            END IF
            STACK(LSTUSD-2) = IST
            STACK(LSTUSD-1) = IWC
            STACK(LSTUSD) = RET
            GOTO 10
C END
C ==================================================================
C
C ==================================================================
C BEGIN
44444       CONTINUE
            IF (LSTUSD .EQ. OLDLST) GOTO 9999
            IF (LSTUSD .LT. OLDLST) THEN
               WRITE (6,'(A)')
     @ ' %EQSTWC-ERR: Stack Underflow. Check code.'
               GOTO 9999
            END IF
            IST = STACK(LSTUSD-2)
            IWC = STACK(LSTUSD-1)
            RET = STACK(LSTUSD)
            LSTUSD = LSTUSD - 3
            IF (RET.EQ.20) GOTO 20
            IF (RET.EQ.30) GOTO 30
            IF (RET.EQ.40) GOTO 40
            IF (RET.EQ.50) GOTO 50
            WRITE(6,'(A)') ' %EQSTWC-ERR: Unknown return address'
C END
C ====================================================================
C
9999  CONTINUE
17789 CONTINUE
8764  CONTINUE
      RETURN
      END
C
      SUBROUTINE ENCODC(C,ST,STMAX,STLEN)
C
      IMPLICIT NONE
C input/output
      INCLUDE 'symbol.inc'
      DOUBLE COMPLEX C
      CHARACTER*(*) ST
      INTEGER STMAX, STLEN
C local
      INTEGER TMPMAX, TMPLEN
      PARAMETER( TMPMAX=40 )
      CHARACTER*(TMPMAX) TEMP
C begin
C
      ST=' '
C
      ST(1:1) = '('
      STLEN = 1
      CALL ENCODF(DBLE(C),TEMP,TMPMAX,TMPLEN)
      ST(2:1+TMPLEN) = TEMP(1:TMPLEN)
      STLEN = STLEN + TMPLEN + 1
      ST(STLEN:STLEN) = ','
      CALL ENCODF(DIMAG(C),TEMP,TMPMAX,TMPLEN)
      ST(STLEN+1:STLEN+TMPLEN) = TEMP(1:TMPLEN)
      STLEN = STLEN + TMPLEN + 1
      ST(STLEN:STLEN) = ')'
C
      RETURN
      END
C
      SUBROUTINE DECODC(CCC,ST,STLEN,OK)
C
C Decodes the string into a double complex number and returns its value.
C It uses two calls to decodf.
C
      IMPLICIT NONE
C input/output
      DOUBLE COMPLEX CCC
      CHARACTER*(*) ST
      INTEGER STLEN
      LOGICAL OK
C local
      INTEGER I, DPLEN1, DPLEN2, START, STOP
      DOUBLE PRECISION DPVAL1, DPVAL2, DECODF
      LOGICAL DONE
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
      OK=.TRUE.
C
      IF ( ST(1:1) .EQ. '(' ) THEN
      IF ( ST(STLEN:STLEN) .NE. ')' ) THEN
      OK = .FALSE.
      ELSE
      START = 2
      STOP = STLEN - 1
      ENDIF
      ELSE
      START = 1
      STOP = STLEN
      ENDIF
C
      DONE = .FALSE.
      DPLEN1 = 0
      I = START
      DO WHILE ( I .LE. STOP .AND. .NOT. DONE )
      IF ( ST(I:I) .EQ. ',' ) THEN
      DPVAL1 = DECODF( ST(START:START+DPLEN1-1), DPLEN1, OK )
      IF ( OK ) THEN
      DPLEN2 = STOP - I
      IF ( DPLEN2 .GT. 0 ) THEN
      DPVAL2 = DECODF( ST(I+1:I+DPLEN2), DPLEN2, OK )
      ELSE
      OK = .FALSE.
      ENDIF
      ENDIF
      DONE = .TRUE.
      ELSE
      DPLEN1 = DPLEN1 + 1
      ENDIF
      I = I + 1
      ENDDO
C
      IF ( OK ) THEN
      CCC= DCMPLX( DPVAL1, DPVAL2 )
      ELSE
      CCC = DCMPLX(ZERO,ZERO)
      ENDIF
C
      RETURN
      END
C ====================================================================
      FUNCTION RINDEX(STR,SUBSTR)
C Same as INDEX, but last occurance first.
C Equivalent to Fortran90 INDEX(STR,SUBSTR,BACK=.TRUE.), or
C the RINDEX extension in some Fortran compilers.
C
C Author: Joseph M. Krahn
C =======================
      IMPLICIT  NONE
C I/O
      INTEGER RINDEX
      CHARACTER*(*) STR,SUBSTR
C local
      INTEGER I
C begin
      I = MAX(LEN(STR)-LEN(SUBSTR),0)
      DO WHILE(I.GT.0)
        IF (STR(I:I+LEN(SUBSTR)-1) .EQ. SUBSTR) EXIT
        I=I-1
      END DO
      RINDEX = I
      RETURN
      END
C ====================================================================
      SUBROUTINE STRFMT(STR,STRMAX,STRLEN,OK,TYPE,DCVAL,DPVAL,
     &                  STVAL,STVAL_LEN,FORM,LFORM,CFORM)
C Author: Joseph M. Krahn
C based on WDSUB2 by Mark McCallum and Axel T. Brunger
C ==========================================
C Requires: LEN(CFORM) >= LFORM+6
      IMPLICIT  NONE
C I/O
      INCLUDE 'funct.inc'
      CHARACTER*(*) STR
      INTEGER STRLEN, STRMAX
      LOGICAL OK
      CHARACTER*2 TYPE
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
      CHARACTER*(*) STVAL, FORM, CFORM
      INTEGER STVAL_LEN,LFORM
      CHARACTER*1 MARK
C local
      INTEGER I
      CHARACTER*1 ERRMSG
C parameter
      MARK=ACHAR(30)
C begin
      OK=.TRUE.
C=============================================================
      IF (TYPE .EQ. 'DP') THEN
      IF (FORM.EQ.' ' .OR. FORM.EQ.'*') THEN
C use free-field output
      CALL ENCODF(DPVAL,STR,STRMAX,STRLEN)
C Use hybrid-36 encoding for integer output when the format is
C exactly 'I4' or 'I5', and the value would otherwise overflow.
      ELSE IF (FORM(1:LFORM).EQ.'I4'.AND.INT(DPVAL).GT.9999) THEN
      CALL HY36ENCODE(4,INT(DPVAL),STR,ERRMSG,I)
      IF (I.NE.0) STR='****'
      STRLEN=4
      ELSE IF (FORM(1:LFORM).EQ.'I5'.AND.INT(DPVAL).GT.99999) THEN
      CALL HY36ENCODE(5,INT(DPVAL),STR,ERRMSG,I)
      IF (I.NE.0) STR='*****'
      STRLEN=5
      ELSE
C use formatted output
      CFORM='('//FORM(1:LFORM)//',A1)'
      IF (INDEX(FORM(1:LFORM),'I').NE.0) THEN
      WRITE(STR(1:STRMAX),CFORM,ERR=99) INT(DPVAL),MARK
      ELSE
      WRITE(STR(1:STRMAX),CFORM,ERR=99) DPVAL,MARK
      END IF
      STRLEN=RINDEX(STR(1:STRMAX),MARK)-1
      IF (STRLEN.LE.0) GOTO 99
      END IF
C=============================================================
      ELSEIF (TYPE .EQ. 'DC') THEN
      IF (FORM.EQ.' ' .OR. FORM.EQ.'*') THEN
C use free-field output
      CALL ENCODC(DCVAL,STR,STRMAX,STRLEN)
      ELSE
C use formatted output
      CFORM='(2'//FORM(1:LFORM)//',A1)'
      WRITE(STR(1:STRMAX),CFORM,ERR=99) DBLE(DCVAL),DIMAG(DCVAL),MARK
      STRLEN=INDEX(STR(1:STRMAX),MARK)-1
      IF (STRLEN.LE.0) GOTO 99
      END IF
C=============================================================
      ELSEIF ( TYPE .EQ. 'ST' ) THEN
      IF (FORM.EQ.' ' .OR. FORM.EQ.'*') THEN
C use free-field output
      CALL COPYST(STR,STRMAX,STRLEN,STVAL,STVAL_LEN)
      ELSE
C use formatted output
      CFORM='('//FORM(1:LFORM)//',A1)'
      WRITE(STR(1:STRMAX),CFORM,ERR=99) STVAL(1:STVAL_LEN),MARK
      STRLEN=INDEX(STR(1:STRMAX),MARK)-1
      IF (STRLEN.LE.0) GOTO 99
      END IF
C=============================================================
      ELSEIF (TYPE .EQ. 'LO') THEN
C logical
      IF (FORM.EQ.' ' .OR. FORM.EQ.'*') THEN
C use free-field output
      CALL COPYST(STR,STRMAX,STRLEN,STVAL,STVAL_LEN)
      ELSE
C use formatted output
      CFORM='('//FORM(1:LFORM)//',A1)'
      WRITE(STR(1:STRMAX),CFORM,ERR=99) STVAL(1:1).EQ."T",MARK
      STRLEN=RINDEX(STR(1:STRMAX),MARK)-1
      IF (STRLEN.LE.0) GOTO 99
      END IF
C=============================================================
      ELSE
      CALL DSPERR('STRFMT','Corrupt variable tables')
      GOTO 100
      ENDIF
      RETURN
99    CONTINUE
      CALL DSPERR('STRFMT','Incorrect FORMAT statement')
100   CONTINUE
      OK=.FALSE.
      STRLEN=1
      STR(1:STRMAX)='?'
      RETURN
      END
