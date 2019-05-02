      SUBROUTINE EXRPN(PROM,RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,
     &                RPNTYP,RPNDOM,RPNLEV,TYPE,DOMAIN,DEPTH,FUNC,OPER,
     &                TYPEDEF,QHERM,ERR,
     &                XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
C
C Converts an expression into Reverse Polish Notation.
C
C PROM is the prompt for parsing,
C RPN is the command stack,
C RPNMX is the maximum number of commands,
C RPNN is the actual number of commands,
C RPNL is the actual length of the command string,
C RPNDB are double precision/double complex constants,
C RPNMLT is the number of arguments for functions,
C RPNTYP is the data type of the command,
C RPNDOM is the domain type of the command,
C RPNLEV is the level of the operation,
C TYPE is the data type of the expression,
C DOMAIN is the domain type of the expression,
C DEPTH is the required depth of the variable stack,
C FUNC is a subroutine that declares all functions,
C OPER is a subroutine that declares all operands,
C TYPEDEF is a subroutine that declares data types and domain types,
C QHERM, is a logical flag passed to TYPEDEF,
C ERR is an error flag.
C XSFNUM, XSFNAM, XSFTYPE, XRHONUM, XRHONAM, XRHOTYP is information
C that are passed to subroutines OPER and TYPEDEF (info for operands)
C
C useful test cases:
C ==================
C ((2+3)*(4+5/(3.+5^2))^(2*3)/(5+6.^3/3.1))*4/5^3^2 =0.000729924
C ((2+3)(4+5/(3.+5^2))^(2 3)/(5+6.^3/3.1))4/5^3^2 =0.000729924
C (2^3)^4 =4096
C 2^3^4  =2.417851639229260E+024
C max(2,max(10,50),max(2,3)) =50.0
C 0<=h<=1/2 and 0<=k<=1/2 and 0<=l<=1/4 and k<h and l<min(k,h/3)
C -1+8 = 7
C +5+4 = 9
C true and ( false or true and 1=3+4) or false and sin(3.)=5.
C true and not 1=2 or not not false
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) PROM
      INTEGER RPNMX, RPNN, RPNX
      CHARACTER*(*) RPN(4,*)
      INTEGER RPNL(4,*)
      DOUBLE COMPLEX RPNDB(4,*)
      INTEGER RPNMLT(*)
      CHARACTER*2 RPNTYP(*), RPNDOM(*)
      INTEGER RPNLEV(*)
      CHARACTER*2 TYPE, DOMAIN
      INTEGER DEPTH
      EXTERNAL FUNC, OPER, TYPEDEF
      LOGICAL QHERM, ERR
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*), XRHOTYP(*)
C local
      LOGICAL OK, QDEBUG
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
      CHARACTER*2 WDTYP
      INTEGER I
C return addresses
      INTEGER LEVEL, RETMAX
      PARAMETER (RETMAX=100)
      INTEGER RETADR(RETMAX)
C counter for number of function arguments
      INTEGER MLTARG(RETMAX)
C function name store
      CHARACTER*(WDMAX) FT(4,RETMAX)
      INTEGER FTL(4,RETMAX)
C function parameter store
      DOUBLE COMPLEX FTV(4,RETMAX)
C function #args store
      INTEGER FTA(RETMAX)
C function depth
      INTEGER FDEPTH(RETMAX)
C type stack
      INTEGER VMAX, VLEVEL
      PARAMETER (VMAX=100)
      CHARACTER*2 ATYPE(VMAX), ADOMAIN(VMAX), FTYPE, FDOMAIN
C parameters
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
C====6====1=========2=========3=========4=========5=========6=========72
C
      IF (RPNX.GT.WDMAX) THEN
      CALL WRNDIE(-5,'EXPRESSION','RPNX is larger than WDMAX')
      END IF
C
      CALL NEXTDO(PROM)
      IF (WD(1:WDLEN).EQ.'DEBUG') THEN
      QDEBUG=.TRUE.
      ELSE
      QDEBUG=.FALSE.
      CALL SAVEWD
      END IF
C
C initialize RPN command stack
      DO RPNN=1,RPNMX
      RPN(1,RPNN)=' '
      RPN(2,RPNN)=' '
      RPN(3,RPNN)=' '
      RPN(4,RPNN)=' '
      RPNL(1,RPNN)=1
      RPNL(2,RPNN)=1
      RPNL(3,RPNN)=1
      RPNL(4,RPNN)=1
      RPNDB(1,RPNN)=DCMPLX(ZERO,ZERO)
      RPNDB(2,RPNN)=DCMPLX(ZERO,ZERO)
      RPNDB(3,RPNN)=DCMPLX(ZERO,ZERO)
      RPNDB(4,RPNN)=DCMPLX(ZERO,ZERO)
      RPNMLT(RPNN)=0
      RPNTYP(RPNN)='??'
      RPNDOM(RPNN)='??'
      RPNLEV(RPNN)=0
      END DO
C set RPN counter
      RPNN=0
C
C
C set initial expression level
      LEVEL=1
C
C set initial type stack level
      VLEVEL=0
      DEPTH=0
C
C ** invoke procedure priority5 **
      RETADR(LEVEL)=5001
      GOTO 5000
5001  CONTINUE
C ** return label                **
C
      IF (VLEVEL.NE.1) THEN
      CALL DSPERR(PROM,'Expression incomplete')
      ERR=.TRUE.
      END IF
      TYPE=ATYPE(1)
      DOMAIN=ADOMAIN(1)
C
      IF (RPNN.EQ.RPNMX) THEN
      CALL WRNDIE(-5,'EXRPN',
     &            'exceeded RPNMX. Expression too complicated.')
      ERR=.TRUE.
      END IF
C
      IF (QDEBUG) THEN
      WRITE(6,'(2A)') ' expression type=',TYPE
      WRITE(6,'(A,I10)') ' expression depth=',DEPTH
      WRITE(6,'(2A)')
     & ' level | command     |   type |  #args | ',
     & '             constants'
      DO I=1,RPNN
      WRITE(6,'(A,I3,1X,13A,I2,4(A,F5.1,F5.1),A)') ' ',
     &         RPNLEV(I),'[',RPN(1,I)(1:4),';',
     &         RPN(2,I)(1:4),';',RPN(3,I)(1:4),';',RPN(4,I)(1:4),'] ',
     &         RPNTYP(I),' ',RPNDOM(I),' ',RPNMLT(I),
     &         ' [',DBLE(RPNDB(1,I)),DIMAG(RPNDB(1,I)),
     &         ' ;',DBLE(RPNDB(2,I)),DIMAG(RPNDB(2,I)),
     &         ' ;',DBLE(RPNDB(3,I)),DIMAG(RPNDB(3,I)),
     &         ' ;',DBLE(RPNDB(4,I)),DIMAG(RPNDB(4,I)),']'
      END DO
      END IF
C
C ** return      **
      GOTO 9999
C ** return      **
C======================================================================
C
C
C======================================================================
C======================================================================
C priority5: logical "or".
C priority6: logical "and".
C priority7: relational operators.
C priority1: add and subtract.
C priority2: multiply and divide.
C priority3: exponentiation.
C priority4: sign, logical "not", constants, symbols,
C            functions, operands.
C======================================================================
C======================================================================
C---- BEGIN PROCEDURE priority5
5000  CONTINUE
C
      IF (LEVEL.GE.RETMAX) THEN
      CALL WRNDIE(-5,PROM,'exceeded RETMAX -- level too deep.')
      ERR=.TRUE.
      ELSE
      LEVEL=LEVEL+1
      END IF
C
C ** invoke procedure priority6    **
      RETADR(LEVEL)=6001
      GOTO 6000
6001  CONTINUE
C ** return label             **
C
50001 CONTINUE
C
      CALL NEXTDO(PROM)
C
C========================================
      IF (WD(1:WDLEN).NE.'OR') GOTO 50002
C========================================
C
C ** invoke procedure priority6    **
      RETADR(LEVEL)=6003
      GOTO 6000
6003  CONTINUE
C ** return label             **
C
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'OR',2)
      CALL TYPEDEF(ERR,'OR',2,VLEVEL,2,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &             QHERM,
     &             XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      CALL EXOPER(PROM,VLEVEL,VMAX,DEPTH,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &            RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
C
      GOTO 50001
C
50002 CONTINUE
      CALL SAVEWD
      LEVEL=LEVEL-1
C
C     return to address:
      IF (RETADR(LEVEL).EQ.5001) GOTO 5001
      IF (RETADR(LEVEL).EQ.5002) GOTO 5002
      IF (RETADR(LEVEL).EQ.5003) GOTO 5003
      IF (RETADR(LEVEL).EQ.5004) GOTO 5004
      WRITE(6,'(A)') ' DO-ERR: Return address unknown. Internal error'
C---- END PROCEDURE priority5 -----------------------------------------
C======================================================================
C======================================================================
C======================================================================
C
C
C
C
C======================================================================
C======================================================================
C======================================================================
C---- BEGIN PROCEDURE priority6 ---------------------------------------
6000  CONTINUE
C
      IF (LEVEL.GE.RETMAX) THEN
      CALL WRNDIE(-5,PROM,'exceeded RETMAX -- level too deep.')
      ERR=.TRUE.
      ELSE
      LEVEL=LEVEL+1
      END IF
C
C ** invoke procedure priority7 **
      RETADR(LEVEL)=7001
      GOTO 7000
7001  CONTINUE
C ** return label              **
C
60001 CONTINUE
C
      CALL NEXTDO(PROM)
C
C========================================
      IF (WD(1:WDLEN).NE.'AND') GOTO 60002
C========================================
C
C ** invoke procedure priority7    **
      RETADR(LEVEL)=7003
      GOTO 7000
7003  CONTINUE
C ** return label             **
C
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'AND',3)
      CALL TYPEDEF(ERR,'AND',3,VLEVEL,2,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &             QHERM,
     &             XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      CALL EXOPER(PROM,VLEVEL,VMAX,DEPTH,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &            RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
C
      GOTO 60001
C
C
60002 CONTINUE
      CALL SAVEWD
      LEVEL=LEVEL-1
C
C     return to address:
      IF (RETADR(LEVEL).EQ.6001) GOTO 6001
      IF (RETADR(LEVEL).EQ.6003) GOTO 6003
      WRITE(6,'(A)') ' DO-ERR: Return address unknown. Internal error'
C---- END PROCEDURE priority6 -----------------------------------------
C======================================================================
C======================================================================
C======================================================================
C
C
C
C
C
C
C======================================================================
C======================================================================
C======================================================================
C---- BEGIN PROCEDURE priority7 ---------------------------------------
7000  CONTINUE
C
      IF (LEVEL.GE.RETMAX) THEN
      CALL WRNDIE(-5,PROM,'exceeded RETMAX -- level too deep.')
      ERR=.TRUE.
      ELSE
      LEVEL=LEVEL+1
      END IF
C
C implicit ATTRibute command
      CALL NEXTDO(PROM)
      IF (WD(1:WDLEN).NE.'ATTR') THEN
      CALL SAVEWD
      END IF
C
C ** invoke procedure priority1 **
      RETADR(LEVEL)=1001
      GOTO 1000
1001  CONTINUE
C ** return label              **
C
C========================================
C relational operator
C
      CALL NEXTDO(PROM)
      IF (WD(1:WDLEN).NE.'='.AND.WD(1:WDLEN).NE.'#'
     &    .AND.WD(1:WDLEN).NE.'>'.AND.WD(1:WDLEN).NE.'<') GOTO 70001
C
C found a relational operator, store it in FT
      FT(1,LEVEL)(1:1)=WD(1:1)
      FTL(1,LEVEL)=1
C
C check if it is >= or <=
      CALL NEXTDO(PROM)
      IF (WD(1:WDLEN).EQ.'=') THEN
      FT(1,LEVEL)(2:2)='='
      FTL(1,LEVEL)=2
      ELSE
      CALL SAVEWD
      END IF
C
C ** invoke procedure priority1 **
      RETADR(LEVEL)=1002
      GOTO 1000
1002  CONTINUE
C ** return label              **
C
C store the relational operator in the command stack
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),FT(1,LEVEL),
     &            FTL(1,LEVEL))
C
      CALL TYPEDEF(ERR,FT(1,LEVEL),FTL(1,LEVEL),VLEVEL,2,ATYPE,FTYPE,
     &             ADOMAIN,FDOMAIN,QHERM,
     &             XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      CALL EXOPER(PROM,VLEVEL,VMAX,DEPTH,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &            RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
C
C Implicit "and": check for commands of the type x < y < z which need
C to be converted into x < y and y < z
      CALL NEXTDO(PROM)
      IF (WD(1:WDLEN).NE.'='.AND.WD(1:WDLEN).NE.'#'
     &    .AND.WD(1:WDLEN).NE.'>'.AND.WD(1:WDLEN).NE.'<') THEN
      CALL SAVEWD
      GOTO 77777
      END IF
C
C put RECALL in the command stack
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'RECALL',6)
      VLEVEL=VLEVEL+1
      RPNTYP(RPNN)=ATYPE(VLEVEL)
      RPNDOM(RPNN)=ADOMAIN(VLEVEL)
      RPNLEV(RPNN)=VLEVEL
C
C found a relational operator, store it in FT
      FT(1,LEVEL)(1:1)=WD(1:1)
      FTL(1,LEVEL)=1
C
C check if it is >= or <=
      CALL NEXTDO(PROM)
      IF (WD(1:WDLEN).EQ.'=') THEN
      FT(1,LEVEL)(2:2)='='
      FTL(1,LEVEL)=2
      ELSE
      CALL SAVEWD
      END IF
C
C ** invoke procedure priority1 **
      RETADR(LEVEL)=1007
      GOTO 1000
1007  CONTINUE
C ** return label              **
C
C store the relational operator in the command stack
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),FT(1,LEVEL),
     &            FTL(1,LEVEL))
      CALL TYPEDEF(ERR,FT(1,LEVEL),FTL(1,LEVEL),VLEVEL,2,ATYPE,FTYPE,
     &             ADOMAIN,FDOMAIN,QHERM,
     &             XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      CALL EXOPER(PROM,VLEVEL,VMAX,DEPTH,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &            RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
C
C add the implicit AND operator
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'AND',3)
      CALL TYPEDEF(ERR,'AND',3,VLEVEL,2,ATYPE,FTYPE,
     &             ADOMAIN,FDOMAIN,QHERM,
     &             XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      CALL EXOPER(PROM,VLEVEL,VMAX,DEPTH,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &            RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
      GOTO 77777
C
C========================================
70001 CONTINUE
      CALL SAVEWD
C
77777 CONTINUE
C
      LEVEL=LEVEL-1
C
C     return to address:
      IF (RETADR(LEVEL).EQ.7001) GOTO 7001
      IF (RETADR(LEVEL).EQ.7002) GOTO 7002
      IF (RETADR(LEVEL).EQ.7003) GOTO 7003
      WRITE(6,'(A)') ' DO-ERR: Return address unknown. Internal error'
C---- END PROCEDURE priority7 -----------------------------------------
C======================================================================
C======================================================================
C
C
C======================================================================
C======================================================================
C---- BEGIN PROCEDURE priority1 ---------------------------------------
1000  CONTINUE
C
      IF (LEVEL.GE.RETMAX) THEN
      CALL WRNDIE(-5,PROM,'exceeded RETMAX -- level too deep.')
      ERR=.TRUE.
      ELSE
      LEVEL=LEVEL+1
      END IF
C
      CALL NEXTDO(PROM)
C========================================
C sign
      IF (WD(1:WDLEN).NE.'-') GOTO 10001
C========================================
C
C ** invoke procedure priority2 **
      RETADR(LEVEL)=2002
      GOTO 2000
2002  CONTINUE
C ** return label              **
C
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'CHS',3)
      CALL TYPEDEF(ERR,'CHS',3,VLEVEL,0,ATYPE,FTYPE,
     &             ADOMAIN,FDOMAIN,QHERM,
     &             XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      RPNTYP(RPNN)=FTYPE
      RPNDOM(RPNN)=FDOMAIN
      RPNLEV(RPNN)=VLEVEL
C
      GOTO 10003
C
C========================================
C sign
10001 IF (WD(1:WDLEN).NE.'+') GOTO 10002
C========================================
C
C ** invoke procedure priority2 **
      RETADR(LEVEL)=2003
      GOTO 2000
2003  CONTINUE
C ** return label              **
C
      CALL TYPEDEF(ERR,'CHS',3,VLEVEL,0,ATYPE,FTYPE,
     &             ADOMAIN,FDOMAIN,QHERM,
     &             XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      RPNTYP(RPNN)=FTYPE
      RPNDOM(RPNN)=FDOMAIN
      RPNLEV(RPNN)=VLEVEL
C
      GOTO 10003
C
10002 CONTINUE
      CALL SAVEWD
C
C ** invoke procedure priority2    **
      RETADR(LEVEL)=2001
      GOTO 2000
2001  CONTINUE
C ** return label             **
C
10003 CONTINUE
C
      CALL NEXTDO(PROM)
C
C========================================
      IF (WD(1:WDLEN).NE.'+') GOTO 10004
C========================================
C
C ** invoke procedure priority2    **
      RETADR(LEVEL)=2005
      GOTO 2000
2005  CONTINUE
C ** return label             **
C
      RPNN=MIN(RPNMX,RPNN+1)
C
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'+',1)
      CALL TYPEDEF(ERR,'+',1,VLEVEL,2,ATYPE,FTYPE,
     &             ADOMAIN,FDOMAIN,QHERM,
     &             XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      CALL EXOPER(PROM,VLEVEL,VMAX,DEPTH,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &            RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
C
      GOTO 10003
C
C========================================
10004 IF (WD(1:WDLEN).NE.'-') GOTO 10005
C========================================
C
C ** invoke procedure priority2    **
      RETADR(LEVEL)=2006
      GOTO 2000
2006  CONTINUE
C ** return label             **
C
      RPNN=MIN(RPNMX,RPNN+1)
C
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'-',1)
      CALL TYPEDEF(ERR,'-',1,VLEVEL,2,ATYPE,FTYPE,
     &             ADOMAIN,FDOMAIN,QHERM,
     &             XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      CALL EXOPER(PROM,VLEVEL,VMAX,DEPTH,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &            RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
C
      GOTO 10003
C
10005 CONTINUE
      CALL SAVEWD
      LEVEL=LEVEL-1
C
C     return to address:
      IF (RETADR(LEVEL).EQ.1001) GOTO 1001
      IF (RETADR(LEVEL).EQ.1002) GOTO 1002
      IF (RETADR(LEVEL).EQ.1007) GOTO 1007
      WRITE(6,'(A)') ' DO-ERR: Return address unknown. Internal error'
C---- END PROCEDURE priority1 -----------------------------------------
C======================================================================
C======================================================================
C======================================================================
C
C
C
C
C======================================================================
C======================================================================
C======================================================================
C---- BEGIN PROCEDURE priority2 ---------------------------------------
2000  CONTINUE
C
      IF (LEVEL.GE.RETMAX) THEN
      CALL WRNDIE(-5,PROM,'exceeded RETMAX -- level too deep.')
      ERR=.TRUE.
      ELSE
      LEVEL=LEVEL+1
      END IF
C
C ** invoke procedure priority3 **
      RETADR(LEVEL)=3001
      GOTO 3000
3001  CONTINUE
C ** return label              **
C
20001 CONTINUE
C
      CALL NEXTDO(PROM)
C
C========================================
      IF (WD(1:WDLEN).NE.'*') GOTO 20002
C========================================
C
C ** invoke procedure priority3    **
      RETADR(LEVEL)=3002
      GOTO 3000
3002  CONTINUE
C ** return label             **
C
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'*',1)
      CALL TYPEDEF(ERR,'*',1,VLEVEL,2,ATYPE,FTYPE,
     &             ADOMAIN,FDOMAIN,QHERM,
     &              XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      CALL EXOPER(PROM,VLEVEL,VMAX,DEPTH,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &            RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
C
      GOTO 20001
C
C=================================
C implicit '*'
20002 CONTINUE
C check the next item.
C constant
      CALL CHKNUM(WD,WDLEN,OK)
      OK=OK.OR.WD(1:WDLEN).EQ.'INFINITY'
C numerical functions
      IF (.NOT.OK) THEN
      CALL FUNC(OK,FT(1,LEVEL),FTL(1,LEVEL),FTV(1,LEVEL),
     &          FTA(LEVEL),.FALSE.,FDEPTH(LEVEL))
      END IF
C numerical operands
      IF (.NOT.OK) CALL OPER(OK,FT(1,LEVEL),FTL(1,LEVEL),.FALSE.,
     &          XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
C symbols and opening parenthesis
      IF (.NOT.OK.AND.WD(1:1).NE.'$'.AND.WD(1:WDLEN).NE.'(') GOTO 20003
C=================================
C
      CALL SAVEWD
C
C ** invoke procedure priority3    **
      RETADR(LEVEL)=3003
      GOTO 3000
3003  CONTINUE
C ** return label             **
C
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'*',1)
      CALL TYPEDEF(ERR,'*',1,VLEVEL,2,ATYPE,FTYPE,
     &             ADOMAIN,FDOMAIN,QHERM,
     &             XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      CALL EXOPER(PROM,VLEVEL,VMAX,DEPTH,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &            RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
C
      GOTO 20001
C========================================
20003 IF (WD(1:WDLEN).NE.'/') GOTO 20004
C========================================
C
C ** invoke procedure priority2    **
      RETADR(LEVEL)=3005
      GOTO 3000
3005  CONTINUE
C ** return label             **
C
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'/',1)
      CALL TYPEDEF(ERR,'/',1,VLEVEL,2,ATYPE,FTYPE,
     &             ADOMAIN,FDOMAIN,QHERM,
     &             XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      CALL EXOPER(PROM,VLEVEL,VMAX,DEPTH,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &            RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
C
      GOTO 20001
C
20004 CONTINUE
      CALL SAVEWD
      LEVEL=LEVEL-1
C
C     return to address:
      IF (RETADR(LEVEL).EQ.2001) GOTO 2001
      IF (RETADR(LEVEL).EQ.2002) GOTO 2002
      IF (RETADR(LEVEL).EQ.2003) GOTO 2003
      IF (RETADR(LEVEL).EQ.2005) GOTO 2005
      IF (RETADR(LEVEL).EQ.2006) GOTO 2006
      WRITE(6,'(A)') ' DO-ERR: Return address unknown. Internal error'
C---- END PROCEDURE priority2 -----------------------------------------
C======================================================================
C======================================================================
C======================================================================
C
C
C
C======================================================================
C======================================================================
C======================================================================
C---- BEGIN PROCEDURE priority3 ---------------------------------------
3000  CONTINUE
C
      IF (LEVEL.GE.RETMAX) THEN
      CALL WRNDIE(-5,PROM,'exceeded RETMAX -- level too deep.')
      ERR=.TRUE.
      ELSE
      LEVEL=LEVEL+1
      END IF
C
C ** invoke procedure priority4 **
      RETADR(LEVEL)=4001
      GOTO 4000
4001  CONTINUE
C ** return label  **
C
      CALL NEXTDO(PROM)
C
C========================================
      IF (WD(1:WDLEN).NE.'^') GOTO 30002
C========================================
C
C ** invoke procedure priority3    **
      RETADR(LEVEL)=3004
      GOTO 3000
3004  CONTINUE
C ** return label             **
C
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'^',1)
      CALL TYPEDEF(ERR,'^',1,VLEVEL,2,ATYPE,FTYPE,
     &             ADOMAIN,FDOMAIN,QHERM,
     &             XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      CALL EXOPER(PROM,VLEVEL,VMAX,DEPTH,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &            RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
C
      GOTO 30001
C
30002 CONTINUE
      CALL SAVEWD
C
30001 CONTINUE
C
      LEVEL=LEVEL-1
C
C     return to address:
      IF (RETADR(LEVEL).EQ.3001) GOTO 3001
      IF (RETADR(LEVEL).EQ.3002) GOTO 3002
      IF (RETADR(LEVEL).EQ.3003) GOTO 3003
      IF (RETADR(LEVEL).EQ.3004) GOTO 3004
      IF (RETADR(LEVEL).EQ.3005) GOTO 3005
      WRITE(6,'(A)') ' DO-ERR: Return address unknown. Internal error'
C---- END PROCEDURE priority3 -----------------------------------------
C======================================================================
C======================================================================
C======================================================================
C
C
C
C
C======================================================================
C======================================================================
C======================================================================
C---- BEGIN PROCEDURE priority4 ---------------------------------------
4000  CONTINUE
C
      IF (LEVEL.GE.RETMAX) THEN
      CALL WRNDIE(-5,PROM,'exceeded RETMAX -- level too deep.')
      ERR=.TRUE.
      ELSE
      LEVEL=LEVEL+1
      END IF
C
      CALL NEXTDO(PROM)
C
C
C========================================
C expression in parenthesis
      IF (WD(1:WDLEN).NE.'(') GOTO 40004
C========================================
C
C ** invoke procedure priority5 **
      RETADR(LEVEL)=5002
      GOTO 5000
5002  CONTINUE
C ** return label                **
C
      CALL NEXTDO(PROM)
C
      IF (.NOT.WD(1:WDLEN).EQ.')') THEN
      CALL DSPERR(PROM,'Expecting ")".')
      ERR=.TRUE.
      END IF
      GOTO 44444
C========================================
C constant.
40004 CONTINUE
C
C infinity
      IF (WD(1:WDLEN).EQ.'INFINITY') THEN
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'CONS',4)
      RPNDB(1,RPNN)=DCMPLX(R4BIG,ZERO)
      CALL EXPUSH(PROM,VLEVEL,VMAX,DEPTH,ATYPE,'DP',ADOMAIN,
     &            '  ',RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
C
C I
      ELSEIF (WD(1:WDLEN).EQ.'I') THEN
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'CONS',4)
      RPNDB(1,RPNN)=DCMPLX(ZERO,ONE)
      CALL EXPUSH(PROM,VLEVEL,VMAX,DEPTH,ATYPE,'DC',ADOMAIN,
     &            '  ',RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
      ELSE
C
C number (integer or floating point)
      CALL CHKNUM(WD,WDLEN,OK)
      IF (.NOT.OK) THEN
      GOTO 40005
      ELSE
      DPVAL=DECODF(WD,WDLEN,OK)
      IF (.NOT.OK) THEN
      CALL DSPERR(PROM,'Error converting constant.')
      ERR=.TRUE.
      ELSE
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'CONS',4)
      RPNDB(1,RPNN)=DCMPLX(DPVAL,ZERO)
      CALL EXPUSH(PROM,VLEVEL,VMAX,DEPTH,ATYPE,'DP',ADOMAIN,
     &            '  ',RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
      END IF
      END IF
      END IF
      GOTO 44444
C =======================================
C symbols
40005 IF (WD(1:1).NE.'$') GOTO 40006
      CALL WDSUB(WD,WDMAX,WDLEN,OK,WDTYP,DPVAL,DCVAL)
      IF (.NOT.OK) THEN
      CALL DSPERR('WDSUB','symbol not found')
      ERR=.TRUE.
      ELSE
      IF (WDTYP.EQ.'DP') THEN
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'CONS',4)
      RPNDB(1,RPNN)=DCMPLX(DPVAL,ZERO)
      CALL EXPUSH(PROM,VLEVEL,VMAX,DEPTH,ATYPE,'DP',ADOMAIN,
     &             '  ',RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
      ELSEIF (WDTYP.EQ.'DC') THEN
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'CONS',4)
      RPNDB(1,RPNN)=DCVAL
      CALL EXPUSH(PROM,VLEVEL,VMAX,DEPTH,ATYPE,'DC',ADOMAIN,
     &            '  ',RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
      ELSEIF (WDTYP.EQ.'ST') THEN
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'"',1)
      CALL ADDST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),WD,WDLEN)
      CALL ADDST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'"',1)
      CALL EXPUSH(PROM,VLEVEL,VMAX,DEPTH,ATYPE,'ST',ADOMAIN,
     &             '  ',RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
      ELSEIF (WDTYP.EQ.'LO') THEN
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),WD,WDLEN)
      CALL EXPUSH(PROM,VLEVEL,VMAX,DEPTH,ATYPE,'LO',ADOMAIN,
     &            '  ',RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
      ELSE
      CALL DSPERR(PROM,'Incorrect data type.')
      ERR=.TRUE.
      END IF
      END IF
      GOTO 44444
C========================================
C functions
40006 CONTINUE
      CALL FUNC(OK,FT(1,LEVEL),FTL(1,LEVEL),FTV(1,LEVEL),
     &          FTA(LEVEL),.TRUE.,FDEPTH(LEVEL))
      IF (.NOT.OK) GOTO 40007
C found a match, have to parse the arguments
C
      MLTARG(LEVEL)=0
C
C bypass argument list if the number of arguments is zero
      IF (FTA(LEVEL).EQ.0) GOTO 99995
C
      CALL NEXTDO(PROM)
      IF (.NOT.WD(1:WDLEN).EQ.'(') THEN
      CALL DSPERR(PROM,'Missing "(".')
      ERR=.TRUE.
      END IF
C
C
      MLTARG(LEVEL)=MLTARG(LEVEL)+1
C ** invoke procedure priority5 **
      RETADR(LEVEL)=5003
      GOTO 5000
5003  CONTINUE
C ** return label                **
C
99992 CALL NEXTDO(PROM)
      IF (WD(1:WDLEN).NE.',') GOTO 99991
C
      MLTARG(LEVEL)=MLTARG(LEVEL)+1
C ** invoke procedure priority5 **
      RETADR(LEVEL)=5004
      GOTO 5000
5004  CONTINUE
C ** return label                **
C
      GOTO 99992
C
C
99991 IF (WD(1:WDLEN).NE.')') THEN
      CALL DSPERR(PROM,'Missing "(".')
      ERR=.TRUE.
      GOTO 44444
      END IF
C
C check if number of arguments match
      IF (FTA(LEVEL).NE.-1) THEN
      IF (FTA(LEVEL).NE.MLTARG(LEVEL)) THEN
      CALL DSPERR(PROM,
     &  'incorrect number of arguments for function '//FT(1,LEVEL))
      ERR=.TRUE.
      GOTO 44444
      END IF
      END IF
C
99995 RPNN=MIN(RPNMX,RPNN+1)
C
C copy function name and auxiliary arguments
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),FT(1,LEVEL),
     &            FTL(1,LEVEL))
      CALL COPYST(RPN(2,RPNN),RPNX,RPNL(2,RPNN),FT(2,LEVEL),
     &            FTL(2,LEVEL))
      CALL COPYST(RPN(3,RPNN),RPNX,RPNL(3,RPNN),FT(3,LEVEL),
     &            FTL(3,LEVEL))
      CALL COPYST(RPN(4,RPNN),RPNX,RPNL(4,RPNN),FT(4,LEVEL),
     &            FTL(4,LEVEL))
C
C copy function parameters
      RPNDB(1,RPNN)=FTV(1,LEVEL)
      RPNDB(2,RPNN)=FTV(2,LEVEL)
      RPNDB(3,RPNN)=FTV(3,LEVEL)
      RPNDB(4,RPNN)=FTV(4,LEVEL)
C
      RPNMLT(RPNN)=MLTARG(LEVEL)
C
      CALL TYPEDEF(ERR,FT(1,LEVEL),FTL(1,LEVEL),VLEVEL,MLTARG(LEVEL),
     &                   ATYPE,FTYPE,ADOMAIN,FDOMAIN,QHERM,
     &                   XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
C
C increase the depth by the depth of the function
      DEPTH=DEPTH+FDEPTH(LEVEL)
C
      CALL EXFNCT(PROM,VLEVEL,VMAX,DEPTH,MLTARG(LEVEL),
     &                  ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &                  RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
C
      GOTO 44444
C========================================
C operands
40007 CALL OPER(OK,FT(1,LEVEL),FTL(1,LEVEL),.TRUE.,
     &          XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      IF (.NOT.OK) GOTO 40008
      RPNN=MIN(RPNMX,RPNN+1)
C
C copy all optional fields
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),
     &            FT(1,LEVEL),FTL(1,LEVEL))
      CALL COPYST(RPN(2,RPNN),RPNX,RPNL(2,RPNN),
     &            FT(2,LEVEL),FTL(2,LEVEL))
      CALL COPYST(RPN(3,RPNN),RPNX,RPNL(3,RPNN),
     &            FT(3,LEVEL),FTL(3,LEVEL))
      CALL COPYST(RPN(4,RPNN),RPNX,RPNL(4,RPNN),
     &            FT(4,LEVEL),FTL(4,LEVEL))
C
      CALL TYPEDEF(ERR,FT(1,LEVEL),FTL(1,LEVEL),VLEVEL,0,
     &                   ATYPE,FTYPE,ADOMAIN,FDOMAIN,QHERM,
     &                   XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      CALL EXPUSH(PROM,VLEVEL,VMAX,DEPTH,ATYPE,FTYPE,ADOMAIN,FDOMAIN,
     &                  RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
      GOTO 44444
C========================================
C quoted strings
40008 IF (WD(1:1).NE.'"') GOTO 40009
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),WD,WDLEN)
      CALL EXPUSH(PROM,VLEVEL,VMAX,DEPTH,ATYPE,'ST',ADOMAIN,'  ',
     &                  RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
      GOTO 44444
C
C========================================
C not
40009 IF (WD(1:WDLEN).NE.'NOT') GOTO 40010
C========================================
C
C ** invoke procedure priority7 **
      RETADR(LEVEL)=7002
      GOTO 7000
7002  CONTINUE
C ** return label              **
C
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'NOT',3)
      CALL TYPEDEF(ERR,'NOT',3,VLEVEL,0,ATYPE,FTYPE,
     &          ADOMAIN,FDOMAIN,QHERM,
     &          XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,XRHOTYP)
      RPNTYP(RPNN)=FTYPE
      RPNDOM(RPNN)=FDOMAIN
      RPNLEV(RPNN)=VLEVEL
      GOTO 44444
C
C========================================
C constant ("TRUE" or "ALL").
40010 CONTINUE
      IF (WD(1:WDLEN).NE.'TRUE'.AND.WD(1:WDLEN).NE.'ALL') GOTO 40011
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'TRUE',4)
      CALL EXPUSH(PROM,VLEVEL,VMAX,DEPTH,ATYPE,'LO',ADOMAIN,'  ',
     &                  RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
      GOTO 44444
C
C========================================
C constant ("FALSE" or "NONE").
40011 CONTINUE
      IF (WD(1:WDLEN).NE.'FALSE'.AND.WD(1:WDLEN).NE.'NONE') GOTO 40012
      RPNN=MIN(RPNMX,RPNN+1)
      CALL COPYST(RPN(1,RPNN),RPNX,RPNL(1,RPNN),'FALSE',5)
      CALL EXPUSH(PROM,VLEVEL,VMAX,DEPTH,ATYPE,'LO',ADOMAIN,'  ',
     &                  RPNTYP(RPNN),RPNDOM(RPNN),RPNLEV(RPNN),ERR)
      GOTO 44444
C========================================
C termination
40012 CONTINUE
      ERR=.TRUE.
      CALL DSPERR(PROM,
     &  'unrecognized statement or variable/type mismatch')
      GOTO 44444
C========================================
C
44444 CONTINUE
C
      LEVEL=LEVEL-1
C
C     return to address:
      IF (RETADR(LEVEL).EQ.4001) GOTO 4001
      WRITE(6,'(A)') ' DO-ERR: Return address unknown. Internal error'
C---- END PROCEDURE priority4 -----------------------------------------
C======================================================================
C======================================================================
C======================================================================
C
C
C
9999  CONTINUE
      RETURN
      END
C======================================================================
      SUBROUTINE EXOPER(PROM,VLEVEL,VMAX,DEPTH,
     &           ATYPE,FTYPE,ADOMAIN,FDOMAIN,RPNTYP,RPNDOM,RPNLEV,ERR)
C
C Performs an operation involving two arguments on the
C type (ATYPE) and domain (ADOMAIN) stack.
C
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      CHARACTER*(*) PROM
      INTEGER VLEVEL, VMAX, DEPTH
      CHARACTER*2 ATYPE(*), FTYPE, ADOMAIN(*), FDOMAIN, RPNTYP
      CHARACTER*2 RPNDOM
      INTEGER RPNLEV
      LOGICAL ERR
C begin
      IF (VLEVEL.LE.1) THEN
      CALL DSPERR(PROM,'Expression incomplete')
      ERR=.TRUE.
      ELSE
      VLEVEL=VLEVEL-1
      ATYPE(VLEVEL)=FTYPE
      ADOMAIN(VLEVEL)=FDOMAIN
      RPNTYP=FTYPE
      RPNDOM=FDOMAIN
      RPNLEV=VLEVEL
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE EXPUSH(PROM,VLEVEL,VMAX,DEPTH,
     &           ATYPE,FTYPE,ADOMAIN,FDOMAIN,RPNTYP,RPNDOM,RPNLEV,ERR)
C
C Performs an operation involving an operand without arguments
C on the type (ATYPE) and domain (ADOMAIN) stack.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      CHARACTER*(*) PROM
      INTEGER VLEVEL, VMAX, DEPTH
      CHARACTER*2 ATYPE(*), FTYPE, ADOMAIN(*), FDOMAIN, RPNTYP
      CHARACTER*2 RPNDOM
      INTEGER RPNLEV
      LOGICAL ERR
C begin
      IF (VLEVEL.GE.VMAX) THEN
      CALL DSPERR(PROM,'Expression is too complicated')
      ERR=.TRUE.
      ELSE
      VLEVEL=VLEVEL+1
      DEPTH=MAX(DEPTH,VLEVEL)
      ATYPE(VLEVEL)=FTYPE
      ADOMAIN(VLEVEL)=FDOMAIN
      RPNTYP=FTYPE
      RPNDOM=FDOMAIN
      RPNLEV=VLEVEL
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE EXFNCT(PROM,VLEVEL,VMAX,DEPTH,NARGS,
     &           ATYPE,FTYPE,ADOMAIN,FDOMAIN,RPNTYP,RPNDOM,RPNLEV,ERR)
C
C Performs a function operation involving NARGS arguments on the
C type (ATYPE) and domain (ADOMAIN) stack.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      CHARACTER*(*) PROM
      INTEGER VLEVEL, VMAX, DEPTH, NARGS
      CHARACTER*2 ATYPE(*), FTYPE, ADOMAIN(*), FDOMAIN, RPNTYP
      CHARACTER*2 RPNDOM
      INTEGER RPNLEV
      LOGICAL ERR
C begin
      VLEVEL=VLEVEL-NARGS
      IF (VLEVEL.GE.VMAX) THEN
      CALL DSPERR(PROM,'Expression is too complicated')
      ERR=.TRUE.
      ELSE
      VLEVEL=VLEVEL+1
      DEPTH=MAX(DEPTH,VLEVEL)
      ATYPE(VLEVEL)=FTYPE
      ADOMAIN(VLEVEL)=FDOMAIN
      RPNTYP=FTYPE
      RPNDOM=FDOMAIN
      RPNLEV=VLEVEL
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE NEXTDO(PROMPT)
C
C Same as NEXTWD expect that the special set of single-character
C words are used during parsing. No symbol substitutions.  Quoted
C strings are returned in quotes.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) PROMPT
C local
      CHARACTER*2 TYPE
      INTEGER I, ICHR
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
      LOGICAL OK, VARFLG(NUMVARFLG), QUPPER
      CHARACTER*(WDMAX) TMPWD
      INTEGER TMPLEN, OFFSET
C begin
C
C no symbol substitutions (except when referenced
C as $$ and of type string, see below)
      QSUBS=.FALSE.
C
C expression mode
      QEXPRS=.TRUE.
      QUPPER=.FALSE.
C
      CALL NEXTWD(PROMPT)
C
C check if symbol is referenced as $$ and is of type string
C
      IF (WD(1:2).EQ.'$$') THEN
      CALL COPYST(WDD,WDDMAX,WDDLEN,WD,WDLEN)
      CALL WDSUB(WDD,WDDMAX,WDDLEN,OK,TYPE,DPVAL,DCVAL)
      IF (TYPE.EQ.'ST'.AND.OK) QUPPER = .TRUE.
C
C check if symbol has STRIP% directive and is of type string
      ELSE IF(WD(1:1).EQ.'$') THEN
      IF(WDLEN.GT.1) THEN
      TMPLEN = WDLEN - 1
      TMPWD(1:TMPLEN) = WD(2:WDLEN)
      CALL DEFGETFLG(TMPWD,TMPLEN,.FALSE.,OFFSET,VARFLG)
      IF(VARFLG(2)) THEN
C has STRIP% directive
      CALL COPYST(WDD,WDDMAX,WDDLEN,WD,WDLEN)
      CALL WDSUB(WDD,WDDMAX,WDDLEN,OK,TYPE,DPVAL,DCVAL)
      IF (TYPE.EQ.'ST'.AND.OK) QUPPER = .TRUE.
      END IF
      END IF
      END IF
C
      IF(QUPPER) THEN
      QQUOT=.FALSE.
C
C convert into upper case
      DO I=1,WDDLEN
      ICHR=ASCIIM(ICHAR(WDD(I:I)))
      WDD(I:I)=CHAR(ICHR)
      END DO
      CALL COPYST(WD,WDMAX,WDLEN,WDD,WDDLEN)
      END IF
C
C
C don't add quotes if this was already done (could happen
C when using SAVEWD).
      IF (QQUOT.AND.(WD(1:1).NE.'"'.OR.WD(WDLEN:WDLEN).NE.'"')) THEN
      WD='"'//WD(1:WDLEN)//'"'
      WDLEN=WDLEN+2
      END IF
      QEXPRS=.FALSE.
      QSUBS=.TRUE.
      RETURN
      END
