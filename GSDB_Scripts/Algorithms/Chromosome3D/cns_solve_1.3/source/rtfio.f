      SUBROUTINE RTFRDR
C
C residue topology file reader
C
C For syntax see main parsing loop "HELP"
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/ouput
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'rtf.inc'
      INCLUDE 'funct.inc'
C local
      LOGICAL OK, AUTANG, AUTDIH, COND, UNDFMA
      INTEGER I, NN, MULT
      INTEGER NRX, NAT, NBO, NTH, NPH, NIM, ND, NA, NGR
      INTEGER NRXPT, NATPT, NBOPT, NTHPT, NPHPT, NIMPT
      INTEGER RTYPE
      CHARACTER*5 W1, W2, W3, W4, NAME5
      CHARACTER*4 OPTION, NAME, BRA, KET
      INTEGER MARK
      DOUBLE PRECISION ANUM, ZERO
      PARAMETER (MARK=-9999, ANUM=9999.0D0, BRA='(   ',KET=')   ')
      PARAMETER (ZERO=0.0D0)
C
C begin
C
C defaults (some other values are initialized in main program)
      AUTANG=.FALSE.
      AUTDIH=.FALSE.
C
C parsing
      CALL PUSEND('RTFRDR>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('RTFRDR>')
      CALL MISCOM('RTFRDR>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-topology')
C
C parse RESEt option
      ELSE IF (WD(1:4).EQ.'RESE') THEN
      CALL RTFINI
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      CALL RTWRIT('OUTPUT')
      ELSE IF (WD(1:4).EQ.'MASS') THEN
C
C this section parses the mass / chemical-type connections (optional)
      CALL NEXTA4('ATOM-TYPE=',NAME)
      I=SRCHC4(NAME,ATCT,1,IATC,MARK)
      IF (I.GT.0) THEN
      CALL NEXTF('ATOM-MASS=',MASDEF(I))
      ELSE
      IF (IATC.GE.MAXAT2) THEN
      CALL WRNDIE(-5,'RTFRDR',
     & 'exceeded MAXAT2 (RTF) parameter --> recompile program')
      ELSE
      IATC=IATC+1
      ATCT(IATC)=NAME
      CALL NEXTF('ATOM-MASS=',MASDEF(IATC))
      END IF
      END IF
C
      ELSE IF (WD(1:4).EQ.'RESI'.OR.WD(1:4).EQ.'PRES') THEN
      UNDFMA=.FALSE.
C
C parsing of (p-) residue specifications
      IF (WD(1:4).EQ.'RESI') THEN
      RTYPE=0
      ELSE
      RTYPE=1
      END IF
      CALL NEXTA4('NAME=',NAME)
      DO I=1,NRTRS
      IF (AA(I).EQ.NAME.AND.RTRTYP(I).EQ.RTYPE) THEN
      WRITE(6,'(2A)') ' %RTFRDR-ERR: duplicate (P-)RESIdue name ',NAME
      END IF
      END DO
C
C set up residue pointers
      IF (NRTRS.GT.0) THEN
      NATPT=NIC(1,NRTRS)
      NBOPT=NIC(2,NRTRS)
      NTHPT=NIC(3,NRTRS)
      NPHPT=NIC(4,NRTRS)
      NIMPT=NIC(5,NRTRS)
      NRXPT=NIC(6,NRTRS)
      ELSE
      NATPT=0
      NBOPT=0
      NTHPT=0
      NPHPT=0
      NIMPT=0
      NRXPT=0
      END IF
C
      IF (NRTRS.GE.MXRTRS) THEN
      CALL WRNDIE(-5,'RTFRDR',
     & 'exceeded MXRTRS (RTF) parameter --> recompile program')
      ELSE
      NRTRS=NRTRS+1
      END IF
      AA(NRTRS)=NAME
      RTRTYP(NRTRS)=RTYPE
C
C defaults
      NAT=0
      NBO=0
      NTH=0
      NPH=0
      NIM=0
      NRX=0
      NA=0
      ND=0
      NGR=0
C
C parsing
      CALL PUSEND('RESIDUE>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('RESIDUE>')
      COND=.TRUE.
      IF (WD(1:4).EQ.'OMIT') THEN
      OPTION='OMIT'
      ELSE IF (WD(1:4).EQ.'ADD '.AND.RTYPE.EQ.1) THEN
      OPTION='ADD '
      ELSE IF (WD(1:4).EQ.'DELE'.AND.RTYPE.EQ.1) THEN
      OPTION='DELE'
      ELSE IF (WD(1:4).EQ.'MODI'.AND.RTYPE.EQ.1) THEN
      OPTION='MODI'
      ELSE
      COND=.FALSE.
      IF (RTYPE.EQ.0) THEN
      OPTION='    '
      ELSE
      OPTION='ADD '
      END IF
      END IF
      IF (COND) CALL NEXTWD('RESIdue>')
      IF (WD(1:4).EQ.'ATOM') THEN
C
C parse atom specifications
      CALL NEXTST('ATOM-iupac-name=',NAME5)
      IF (SRCHC4(NAME5,FTP,NATPT+1,NATPT+NAT,MARK).GT.0) THEN
      WRITE(6,'(2A)') ' %RTFIO-ERR: Duplicate atom name ',NAME5
      END IF
      IF (NAT+NATPT.GE.MXRTA) THEN
      CALL WRNDIE(-5,'RTFRDR',
     & 'exceeded MXRTA (RTF) parameter --> recompile program')
      ELSE
      NAT=NAT+1
      END IF
      FTP(NAT+NATPT)=NAME5
      DELAT(NAT+NATPT)=OPTION
C
C fix to make sure that atoms are at least in one group for RESIdues
C to account for the case when the group statement has been left out
C the RESIdue completely
      IF (RTYPE.EQ.0.AND.NGR.EQ.0) NGR=1
C
      GRPR(NAT+NATPT)=NGR
C
C defaults
      MXN(NAT+NATPT)=NRX+NRXPT
      MAC(NAT+NATPT)='    '
      IF (OPTION.EQ.'MODI') THEN
      CHG(NAT+NATPT)=ANUM
      ELSE
      CHG(NAT+NATPT)=ZERO
      END IF
      ARMASS(NAT+NATPT)=ANUM
C
C parsing
      CALL PUSEND('ATOM>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('ATOM>')
      IF (WD(1:4).EQ.'CHAR') THEN
      CALL NEXTF('CHARge=',CHG(NAT+NATPT))
      ELSE IF (WD(1:4).EQ.'MASS') THEN
      CALL NEXTF('MASS=',ARMASS(NAT+NATPT))
      ELSE IF (WD(1:4).EQ.'TYPE') THEN
      CALL NEXTA4('ATOM-TYPE=',MAC(NAT+NATPT))
      ELSE IF (WD(1:4).EQ.'EXCL') THEN
      CALL NEXTST('EXCLusion=',NAME5)
      IF (NAME5.NE.BRA) THEN
      WRITE(6,'(A)') ' %EXCL-ERR: "(" expected'
      ELSE
      CALL NEXTST('EXCLusions:',NAME5)
      DO WHILE (NAME5.NE.KET)
      IF (NRX+NRXPT.GE.MXRTX) THEN
      CALL WRNDIE(-5,'RTFRDR',
     & 'exceeded MXRTX (RTF) parameter --> recompile program')
      ELSE
      NRX=NRX+1
      END IF
      MNB(NRX+NRXPT)=NAME5
      CALL NEXTST('EXCLusions:',NAME5)
      END DO
      MXN(NAT+NATPT) = NRX+NRXPT
      END IF
C
      ELSE
      CALL CHKEND('ATOM>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
C if required process mass default
      IF (ABS(ARMASS(NAT+NATPT)).GE.ANUM-1.0D0) THEN
      IF (OPTION.EQ.'    '.OR.MAC(NAT+NATPT).NE.'    ') THEN
      I=SRCHC4(MAC(NAT+NATPT),ATCT,1,IATC,MARK)
      IF (I.NE.MARK) THEN
      ARMASS(NAT+NATPT)=MASDEF(I)
      ELSE
      UNDFMA=.TRUE.
      END IF
      END IF
      END IF
C
      ELSE IF (WD(1:4).EQ.'BOND') THEN
C
C parse bond specifications
      CALL NEXTST('BOND-I=',W1)
      CALL NEXTST('BOND-J=',W2)
      IF (W1 .EQ. W2) THEN
      WRITE(6,'(5A)') ' %RTFRDR-ERR: self bond ',W1,' ',W2,' ignored'
      ELSE
      IF (NBO+NBOPT.GE.MXRTB) THEN
      CALL WRNDIE(-5,'RTFRDR',
     & 'exceeded MXRTB (RTF) parameter --> recompile program')
      ELSE
      NBO=NBO+1
      END IF
      DELBD(NBO+NBOPT)=OPTION
      MIB(NBO+NBOPT) = W1
      MJB(NBO+NBOPT) = W2
      END IF
C
      ELSE IF (WD(1:4).EQ.'ANGL') THEN
C
C parse angle specifications
      CALL NEXTST('ANGLE-I=',W1)
      CALL NEXTST('ANGLE-J=',W2)
      CALL NEXTST('ANGLE-K=',W3)
      OK=(W1.NE.W2.AND.W2.NE.W3.AND.W3.NE.W1)
      IF (.NOT.OK) THEN
      WRITE(6,'(7A)') ' %RTFRDR-ERR: self angle ',W1,' ',W2,' ',W3,
     &               ' ignored'
      ELSE
      IF (NTH+NTHPT.GE.MXRTT) THEN
      CALL WRNDIE(-5,'RTFRDR',
     & 'exceeded MXRTT (RTF) parameter --> recompile program')
      ELSE
      NTH=NTH+1
      END IF
      DELAN(NTH+NTHPT)=OPTION
      MIT(NTH+NTHPT)=W1
      MJT(NTH+NTHPT)=W2
      MKT(NTH+NTHPT)=W3
      END IF
C
      ELSE IF (WD(1:4).EQ.'DIHE') THEN
C
C parse dihedral specifications
      CALL NEXTST('DIHEdral-I=',W1)
      CALL NEXTST('DIHEdral-J=',W2)
      CALL NEXTST('DIHEdral-K=',W3)
      CALL NEXTST('DIHEdral-L=',W4)
C get multiplicity
      CALL NEXTWD('RESIDUE>')
      IF (WD(1:4).EQ.'MULT') THEN
      CALL NEXTI('MULTiplity=',MULT)
      ELSE
      CALL SAVEWD
      MULT=1
      END IF
C
      OK=W1.NE.W2.AND.W1.NE.W3.AND.W1.NE.W4.AND.W2.NE.W3.AND.
     &   W2.NE.W4.AND.W3.NE.W4
      IF (.NOT.OK) THEN
      WRITE(6,'(9A)') ' %RTFRDR-ERR: self dihedral ',W1,' ',
     &                W2,' ',W3,' ',W4,' ignored'
      ELSE
C
      DO NN=1,MULT
      IF (NPH+NPHPT.GE.MXRTP) THEN
      CALL WRNDIE(-5,'RTFRDR',
     & 'exceeded MXRTP (RTF) parameter --> recompile program')
      ELSE
      NPH=NPH+1
      END IF
      DELPT(NPH+NPHPT)=OPTION
      MIP(NPH+NPHPT)=W1
      MJP(NPH+NPHPT)=W2
      MKP(NPH+NPHPT)=W3
      MLP(NPH+NPHPT)=W4
      END DO
C
      END IF
C
      ELSE IF (WD(1:4).EQ.'IMPR') THEN
C
C parse improper specifications
      CALL NEXTST('IMPRoper-I=',W1)
      CALL NEXTST('IMPRoper-J=',W2)
      CALL NEXTST('IMPRoper-K=',W3)
      CALL NEXTST('IMPRoper-L=',W4)
C get multiplicity
      CALL NEXTWD('RESIDUE>')
      IF (WD(1:4).EQ.'MULT') THEN
      CALL NEXTI('MULTiplity=',MULT)
      ELSE
      CALL SAVEWD
      MULT=1
      END IF
C
      OK=W1.NE.W2.AND.W1.NE.W3.AND.W1.NE.W4.AND.
     &   W2.NE.W3.AND.W2.NE.W4.AND.W3.NE.W4
      IF (.NOT.OK) THEN
      WRITE(6,'(9A)') ' %RTFRDR-ERR: self improper',W1,' ',
     &               W2,' ',W3,' ',W4,' ignored'
      ELSE
C
      DO NN=1,MULT
      IF (NIM+NIMPT.GE.MXRTI) THEN
      CALL WRNDIE(-5,'RTFRDR',
     & 'exceeded MXRTI (RTF) parameter --> recompile program')
      ELSE
      NIM=NIM+1
      END IF
      DELMT(NIM+NIMPT)=OPTION
      MIM(NIM+NIMPT)=W1
      MJM(NIM+NIMPT)=W2
      MKM(NIM+NIMPT)=W3
      MLM(NIM+NIMPT)=W4
      END DO
C
      END IF
C
      ELSE IF (WD(1:4).EQ.'DONO') THEN
C
C parse donor specifications - for backwards compatibility only
      CALL NEXTST('DONOr-hydrogen=',W1)
      CALL NEXTST('DONOr-atom=',W2)
C
      ELSE IF (WD(1:4).EQ.'ACCE') THEN
C
C parse acceptor specifications - for backwards compatibility only
      CALL NEXTST('ACCEptor-atom=',W1)
      CALL NEXTST('ANTEcedent-atom=',W2)
C
C parse group boundary
      ELSE IF (WD(1:4).EQ.'GROU') THEN
      NGR=NGR+1
C
      ELSE
      CALL CHKEND('RESIDUE>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
C autogenerate angles and dihedrals if requested
      IF (AUTANG) THEN
      CALL AUTOAN(NBO,NBOPT,DELBD,MIB,MJB,MXRTT,NTH,NTHPT,DELAN,
     &                  MIT,MJT,MKT)
      END IF
      IF (AUTDIH) THEN
      CALL AUTODI(NBO,NBOPT,DELBD,MIB,MJB,MXRTP,NPH,NPHPT,DELPT,
     &                  MIP,MJP,MKP,MLP)
      END IF
C
C update pointers
      NIC(1,NRTRS)=NATPT+NAT
      NIC(2,NRTRS)=NBOPT+NBO
      NIC(3,NRTRS)=NTHPT+NTH
      NIC(4,NRTRS)=NPHPT+NPH
      NIC(5,NRTRS)=NIMPT+NIM
      NIC(6,NRTRS)=NRXPT+NRX
      NIC(7,NRTRS)=0
      NIC(8,NRTRS)=0
      NIC(9,NRTRS)=NGR
C
      IF (UNDFMA) THEN
      WRITE(6,'(2A)') ' %RTFIO-warning: undefined masses in residue ',
     & NAME
      CALL WRNDIE(+5,'RTFIO',
     & 'there are undefined masses.')
      END IF
C
C
      ELSE IF (WD(1:4).EQ.'AUTO') THEN
C
C parsing of autogeneration options
      CALL PUSEND('AUTOGENERATE>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('AUTOGENERATE>')
      IF (WD(1:4).EQ.'ANGL') THEN
      CALL NEXTLO('ANGLe=',AUTANG)
      ELSE IF (WD(1:4).EQ.'DIHE') THEN
      CALL NEXTLO('DIHEdral=',AUTDIH)
      ELSE
      CALL CHKEND('AUTOGENERATE>',DONE)
      END IF
      END DO
      DONE=.FALSE.
      ELSE
      CALL CHKEND('RTFRDR>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
      RETURN
      END
C
      SUBROUTINE RTWRIT(FILE)
C
C Writes topology file information
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/ouput
      INCLUDE 'rtf.inc'
      CHARACTER*(*) FILE
C local
      INTEGER I, J, K, I0
      INTEGER NATPT, NBOPT, NTHPT, NPHPT, NIMPT
      INTEGER KXSTRT, CURRGR
      INTEGER MARK, UNIT
      LOGICAL ERROR
      DOUBLE PRECISION ANUM
      CHARACTER*40 OP, TY
      CHARACTER*18 CH, CMA
      INTEGER OPLEN, TYLEN, CHLEN, CMALEN
      PARAMETER (I0=0,MARK=-9999,ANUM=9999.0D0)
C begin
      CALL ASSFIL(FILE,UNIT,'WRITE','FORMATTED',ERROR)
      IF (.NOT.ERROR) THEN
C
C print topology file
      IF (FILE.NE.'OUTPUT') CALL WRTITL(UNIT,0)
      WRITE(UNIT,'(A)') ' '
      DO K=1,IATC
      WRITE(UNIT,'(3A,F8.4)') ' MASS ',ATCT(K),' ',MASDEF(K)
      END DO
C
      DO I=1,NRTRS
      WRITE(UNIT,'(A)')
     1' !-----------------------------------------------------------'
      WRITE(UNIT,'(A)') ' '
      IF (RTRTYP(I).GT.0) THEN
      WRITE(UNIT,'(2A)') ' PRESidue ',AA(I)
      ELSE
      WRITE(UNIT,'(2A)') ' RESIdue ',AA(I)
      END IF
      WRITE(UNIT,'(A)') ' '
C
      IF (I.GT.1) THEN
      NATPT=NIC(1,I-1)
      NBOPT=NIC(2,I-1)
      NTHPT=NIC(3,I-1)
      NPHPT=NIC(4,I-1)
      NIMPT=NIC(5,I-1)
      ELSE
      NATPT=0
      NBOPT=0
      NTHPT=0
      NPHPT=0
      NIMPT=0
      END IF
C
      CURRGR=0
      DO K=NATPT+1,NIC(1,I)
      IF (K.GT.1) THEN
      KXSTRT=MXN(K-1)
      ELSE
      KXSTRT=0
      END IF
      IF (GRPR(K).NE.CURRGR) THEN
      CURRGR=GRPR(K)
      WRITE(UNIT,'(A)') ' GROUp'
      END IF
      IF (DELAT(K).NE.'    ') THEN
      OP='      '//DELAT(K)
      OPLEN=10
      ELSE
      OP='      '
      OPLEN=6
      END IF
      IF (MAC(K).EQ.'    ') THEN
      TY=' '
      TYLEN=1
      ELSE
      TY=' TYPE='//MAC(K)
      TYLEN=10
      END IF
      IF (ABS(CHG(K)).GE.ANUM-1.0D0) THEN
      CH=' '
      CHLEN=1
      ELSE
      WRITE(CH,'(A,F10.4)') ' CHARge=',CHG(K)
      CHLEN=18
      END IF
      IF (ABS(ARMASS(K)).GE.ANUM-1.0D0) THEN
      CMA=' '
      CMALEN=1
      ELSE
      WRITE(CMA,'(A,F10.4)') ' MASS=',ARMASS(K)
      CMALEN=18
      END IF
C
      IF (MXN(K)-KXSTRT.GE.1) THEN
      WRITE(UNIT,'(6A/,A)')  OP(1:OPLEN),' ATOM  ',FTP(K),TY(1:TYLEN),
     &                  CH(1:CHLEN),CMA(1:CMALEN),'        EXCLude=('
      WRITE(UNIT,'(10X,26A)') (MNB(J),' ',J=KXSTRT+1,MXN(K))
      WRITE(UNIT,'(A/,A)') '        )','      END !ATOM'
      ELSE
      WRITE(UNIT,'(7A)')  OP(1:OPLEN),' ATOM  ',FTP(K),TY(1:TYLEN),
     &                    CH(1:CHLEN),CMA(1:CMALEN),' END'
      END IF
      KXSTRT=MXN(K)+1
      END DO
      IF (CURRGR.NE.0) WRITE(UNIT,'(A)') ' !END GROUp'
C
      IF (NIC(2,I)-NBOPT.GT.0) WRITE(UNIT,'(A)') ' '
      DO K=NBOPT+1,NIC(2,I)
      WRITE(UNIT,'(6A)') ' ',DELBD(K),' BOND  ',MIB(K),' ',MJB(K)
      END DO
C
      IF (NIC(3,I)-NTHPT.GT.0) WRITE(UNIT,'(A)') ' '
      DO K=NTHPT+1,NIC(3,I)
      WRITE(UNIT,'(8A)') ' ',DELAN(K),' ANGLe  ',MIT(K),' ',
     &                    MJT(K),' ',MKT(K)
      END DO
C
      IF (NIC(4,I)-NPHPT.GT.0) WRITE(UNIT,'(A)') ' '
      DO K=NPHPT+1,NIC(4,I)
      WRITE(UNIT,'(10A)') ' ',DELPT(K),' DIHEdral  ',MIP(K),' ',
     &                    MJP(K),' ',MKP(K),' ',MLP(K)
      END DO
C
      IF (NIC(5,I)-NIMPT.GT.0) WRITE(UNIT,'(A)') ' '
      DO K=NIMPT+1,NIC(5,I)
      WRITE(UNIT,'(10A)') ' ',DELMT(K),' IMPRoper  ',MIM(K),' ',
     &                    MJM(K),' ',MKM(K),' ',MLM(K)
      END DO
C
      WRITE(UNIT,'(A)') ' '
      WRITE(UNIT,'(A)') ' '
      WRITE(UNIT,'(A,A)') ' END ! RESIdue ',AA(I)
C
      END DO
C
      CALL VCLOSE(UNIT,'KEEP',ERROR)
      END IF
C
      RETURN
      END
C
      SUBROUTINE AUTOAN(NBO,NBOPT,DELBD,MIB,MJB,MXRTT,NTH,NTHPT,DELAN,
     &                  MIT,MJT,MKT)
C
C Automatically generates an angle list for a residue. All possible
C angles are generated. The generation
C is entirely based on the bond list of the residue. A specific angle
C may be omitted by including the entry "OMIT ANGLe a b c" in the
C residue specifications.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER NBO, NBOPT
      CHARACTER*4 DELBD(*)
      CHARACTER*5 MIB(*), MJB(*)
      INTEGER MXRTT, NTH, NTHPT
      CHARACTER*4 DELAN(MXRTT)
      CHARACTER*5 MIT(MXRTT), MJT(MXRTT), MKT(MXRTT)
C local
      INTEGER I, IOLD, II, J, NA
      CHARACTER*5 WI, WJ, WK, STEMP
      LOGICAL COND
C begin
      NA=NTH+NTHPT
C
C make triangular search through (bond-list * bond-list)
      DO I=NBOPT+1,NBOPT+NBO
      IF (DELBD(I).NE.'DELE') THEN
      DO J=I+1,NBOPT+NBO
      IF (DELBD(I).NE.'DELE') THEN
      COND=.TRUE.
      IF (MIB(I).EQ.MIB(J)) THEN
      WI=MJB(I)
      WJ=MIB(I)
      WK=MJB(J)
      ELSE IF (MIB(I).EQ.MJB(J)) THEN
      WI=MJB(I)
      WJ=MIB(I)
      WK=MIB(J)
      ELSE IF (MJB(I).EQ.MIB(J)) THEN
      WI=MIB(I)
      WJ=MJB(I)
      WK=MJB(J)
      ELSE IF (MJB(I).EQ.MJB(J)) THEN
      WI=MIB(I)
      WJ=MJB(I)
      WK=MIB(J)
      ELSE
      COND=.FALSE.
      END IF
      IF (COND) THEN
C
C now check whether the angle WI, WJ, WK is already present,
C which means it should be omitted during the autogeneration
      IOLD=NTHPT+1
      DO WHILE (IOLD.LE.NTHPT+NTH.AND.COND)
C note the permutation!
      COND=
     &       (MIT(IOLD).NE.WI
     &    .OR.MJT(IOLD).NE.WJ
     &    .OR.MKT(IOLD).NE.WK)
     &  .AND.
     &        (MIT(IOLD).NE.WK
     &     .OR.MJT(IOLD).NE.WJ
     &     .OR.MKT(IOLD).NE.WI)
      IOLD=IOLD+1
      END DO
      IF (COND) THEN
      IF (NA.GE.MXRTT) THEN
      CALL WRNDIE(-5,'AUTOAN',
     & 'exceeded MXRTT (RTF) parameter --> recompile program')
      ELSE
      NA=NA+1
      MIT(NA)=WI
      MJT(NA)=WJ
      MKT(NA)=WK
      DELAN(NA)='    '
      END IF
      END IF
      END IF
      END IF
      END DO
      END IF
      END DO
C
C finally delete all entries with "OMIT"
      II=NTHPT
      DO I=NTHPT+1,NA
      IF (DELAN(I).NE.'OMIT') THEN
      II=II+1
      STEMP=MIT(I)
      MIT(II)=STEMP
      STEMP=MJT(I)
      MJT(II)=STEMP
      STEMP=MKT(I)
      MKT(II)=STEMP
      STEMP=DELAN(I)
      DELAN(II)=STEMP
      END IF
      END DO
      NTH=II-NTHPT
C
      RETURN
      END
C
      SUBROUTINE AUTODI(NBO,NBOPT,DELBD,MIB,MJB,MXRTP,NPH,NPHPT,DELPT,
     &                  MIP,MJP,MKP,MLP)
C
C Automatically generates an dihedral list for a residue. All
C possible dihedrals are generated.
C The generation is entirely based on the bond list of the residue.
C A specific dihedral may be omitted by including the entry
C "OMIT DIHEdral a b c d" in the residue specifications. Already
C present dihedrals are not repeated.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER NBO, NBOPT
      CHARACTER*4 DELBD(*)
      CHARACTER*5 MIB(*), MJB(*)
      INTEGER MXRTP, NPH, NPHPT
      CHARACTER*4 DELPT(MXRTP)
      CHARACTER*5 MIP(MXRTP), MJP(MXRTP), MKP(MXRTP), MLP(MXRTP)
C local
      INTEGER I, II, J, K, IOLD,  NA
      LOGICAL COND
      CHARACTER*5 WI, WJ, WK, WL, STEMP
C begin
      NA=NPH+NPHPT
C
      DO I=NBOPT+1,NBOPT+NBO
      IF (DELBD(I).NE.'DELE') THEN
      DO J=NBOPT+1,NBOPT+NBO
      IF (DELBD(J).NE.'DELE'.AND.J.NE.I) THEN
      COND=.TRUE.
      IF (MIB(I).EQ.MIB(J)) THEN
      WK=MJB(I)
      WJ=MIB(I)
      WI=MJB(J)
      ELSE IF (MIB(I).EQ.MJB(J)) THEN
      WK=MJB(I)
      WJ=MIB(I)
      WI=MIB(J)
      ELSE IF (MJB(I).EQ.MIB(J)) THEN
      WK=MIB(I)
      WJ=MJB(I)
      WI=MJB(J)
      ELSE IF (MJB(I).EQ.MJB(J)) THEN
      WK=MIB(I)
      WJ=MJB(I)
      WI=MIB(J)
      ELSE
      COND=.FALSE.
      END IF
      IF (COND) THEN
      DO K=J+1,NBOPT+NBO
      IF (DELBD(K).NE.'DELE'.AND.K.NE.I) THEN
      COND=.TRUE.
      IF (MIB(K).EQ.WK) THEN
      WL=MJB(K)
      ELSE IF (MJB(K).EQ.WK) THEN
      WL=MIB(K)
      ELSE
      COND=.FALSE.
      END IF
      IF (COND) THEN
C
C now check whether the dihedral WI, WJ, WK, WL is already present,
C which means it should be omitted during the autogeneration
      IOLD=NPHPT+1
      DO WHILE (IOLD.LE.NPHPT+NPH.AND.COND)
C note the permutation!
      COND=(MIP(IOLD).NE.WI
     &  .OR.MJP(IOLD).NE.WJ
     &  .OR.MKP(IOLD).NE.WK
     &  .OR.MLP(IOLD).NE.WL)
     &     .AND.
     &      (MIP(IOLD).NE.WL
     &   .OR.MJP(IOLD).NE.WK
     &   .OR.MKP(IOLD).NE.WJ
     &   .OR.MLP(IOLD).NE.WI)
      IOLD=IOLD+1
      END DO
      IF (COND) THEN
      IF (NA.GE.MXRTP) THEN
      CALL WRNDIE(-5,'AUTOAN',
     & 'exceeded MXRTP (RTF) parameter --> recompile program')
      ELSE
      NA=NA+1
      MIP(NA)=WI
      MJP(NA)=WJ
      MKP(NA)=WK
      MLP(NA)=WL
      DELPT(NA)='    '
      END IF
      END IF
      END IF
      END IF
      END DO
      END IF
      END IF
      END DO
      END IF
      END DO
C
C finally, delete all entries with with "OMIT"
      II=NPHPT
      DO I=NPHPT+1,NA
      IF (DELPT(I).NE.'OMIT') THEN
      II=II+1
      STEMP=MIP(I)
      MIP(II)=STEMP
      STEMP=MJP(I)
      MJP(II)=STEMP
      STEMP=MKP(I)
      MKP(II)=STEMP
      STEMP=MLP(I)
      MLP(II)=STEMP
      STEMP=DELPT(I)
      DELPT(II)=STEMP
      END IF
      END DO
      NPH=II-NPHPT
C
      RETURN
      END
C
      SUBROUTINE RTFINI
C
C Routine  initializes the RTF topology data structure
C
C Author: Axel T. Brunger
      IMPLICIT NONE
C I/O
      INCLUDE 'rtf.inc'
C begin
      NRTRS=0
      IATC=0
      RETURN
      END
