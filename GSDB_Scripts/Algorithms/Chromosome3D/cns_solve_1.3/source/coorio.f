      SUBROUTINE CWRITE
C
C Write selected coordinates on specified file in PDB format.
C for syntax see "HELP"
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'deriv.inc'
C local
      CHARACTER*4 SFROM, FORM
      INTEGER FLAGS, NS
C defaults
      OFILE='OUTPUT'
      FLAGS=ALLHP(INTEG4(NATOM))
      CALL FILL4(HEAP(FLAGS),NATOM,1)
      NS=NATOM
      SFROM='MAIN'
      FORM='CNS'
C
C parsing
      CALL PUSEND('WRITE-COOR>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('WRITE-COOR>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-write-coordinate')
C
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
      ELSE IF (WD(1:4).EQ.'FROM') THEN
      CALL NEXTA4('FROM=',SFROM)
      ELSE IF (WD(1:4).EQ.'FORM') THEN
      CALL NEXTA4('FORMat=',FORM)
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(HEAP(FLAGS),NS,X,Y,Z,.TRUE.)
      ELSE
      CALL CHKEND('WRITE-COOR>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      IF (SFROM.EQ.'MAIN') THEN
      CALL CWRIT2(OFILE,X,Y,Z,WMAIN,QMAIN,HEAP(FLAGS),NS,FORM)
      ELSE IF (SFROM.EQ.'COMP') THEN
      CALL CWRIT2(OFILE,XCOMP,YCOMP,ZCOMP,WCOMP,QCOMP,HEAP(FLAGS),
     &            NS,FORM)
      ELSE IF (SFROM.EQ.'REFE') THEN
      CALL CWRIT2(OFILE,REFX,REFY,REFZ,KCNSTR,KCNSTR,HEAP(FLAGS),
     &            NS,FORM)
      ELSE IF (SFROM.EQ.'DERI') THEN
      CALL CWRIT2(OFILE,DX,DY,DZ,FBETA,FBETA,HEAP(FLAGS),NS,FORM)
      ELSE
      WRITE(6,'(A)') ' %WRITE-COOR-ERR: unknown FROM parameter:',SFROM
      END IF
C
      CALL FREHP(FLAGS,INTEG4(NATOM))
C
      RETURN
      END
C
      SUBROUTINE CWRIT2(FILE,X,Y,Z,WMAIN,QMAIN,ISLCT,NSLCT,FORM)
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) FILE
      DOUBLE PRECISION X(*), Y(*), Z(*), WMAIN(*), QMAIN(*)
      INTEGER ISLCT(*), NSLCT
      CHARACTER*4 FORM
C local
      INTEGER I, J, UNIT, L, SEGLEN
      CHARACTER*4 NAME, SEGTMP, ELCH
      CHARACTER*5 A
      CHARACTER*6 KEYWORD
      LOGICAL QERROR, FOUND
      CHARACTER*5 I36
      CHARACTER*100 ERRMSG
      INTEGER ERRMSG_LEN
      CHARACTER*2 ELEM
C data
      INTEGER ELEM2_SIZE
      PARAMETER(ELEM2_SIZE=97)
      CHARACTER*2 ELEM2_LIST(ELEM2_SIZE)
      DATA ELEM2_LIST /
     &  'HE', 'LI', 'BE', 'NE', 'NA', 'MG', 'AL', 'SI', 'CL', 'AR', 
     &  'CA', 'SC', 'TI', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 
     &  'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'ZR', 'NB', 
     &  'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 
     &  'TE', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 
     &  'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 
     &  'TA', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 
     &  'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'NP', 'PU', 
     &  'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 
     &  'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'RG' /
C
C begin
      IF (NSLCT.LE.0) THEN
      WRITE(6,'(A)') ' CWRITE: zero atoms specified. No file created'
      ELSE
      CALL ASSFIL(FILE,UNIT,'WRITE','FORMATTED',QERROR)
      IF (.NOT.QERROR) THEN
      IF(NSLCT.LT.NATOM) WRITE(6,'(A)') ' CWRITE: using atom subset.'
C
C Brookhaven PDB format
C
CTOM   1223  O   GLY   153     -11.704  -9.200    .489  1.00  0.80           O
C
C     ccccc ''''Iyyy O,,,,L   ........>>>>>>>>////////ppppppiiiiii iii      ee
C234567890123456789012345678901234567890123456789012345678901234567890123456789
C                         ^ insertion character
C                    ^ chain identifier
C               ^ altid
C           ^ additional character for some atom names (mostly h's)
C
C record name: 1-6
C
C atom serial number: 7-11
C
C atom name: 13-16, position depends on element: 13-14 is right-justified element name
C
C alternate location identifier: 17  (should be A or B)
C
C residue name: 18-20, left justified
C
C chain identifier: 22
C
C residue number: 23-26, right justified
C
C code for insertion of residues: 27
C
C element symbol name, right justified: 77-78
C 
C
      IF (FORM(1:3).EQ.'CNS'.OR.FORM(1:4).EQ.'PDBO') THEN
      CALL WRTITL(UNIT,1)
      END IF
      DO I=1,NATOM
      IF (ISLCT(I).EQ.1) THEN
C
C right-adjust the resid
      L=4
      CALL TRIMM(RESID(I),L)
      A='     '
C
C special treatment for RESSEQ and ICODE, hybrid-36, ATB 02/15/10
C test if right-most character in RESID is a letter -
C   if so, it is an insertion character but only if
C   the left-most character is blank or a number.
C   If the left-most character is a letter, then the
C   RESID is a hybrid-36 field and is not shifted. 
      IF (LGE(RESID(I)(L:L),'A').AND.
     &    (RESID(I)(1:1).EQ.' '.OR..NOT.
     &     LGE(RESID(I)(1:1),'A') ) ) THEN
C insertion character found: shift it by one unit to the right 
      A(5-L+1:5)=RESID(I)(1:L)
      ELSE
      A(4-L+1:4)=RESID(I)(1:L)
      END IF
C
C shift atom names when they exceed 3 characters
      IF (TYPE(I)(4:4).EQ.' ') THEN
      NAME=' '//TYPE(I)(1:3)
      ELSE
      NAME=TYPE(I)
      END IF
C 
C convert ATOM serial number into hybrid_36 format ATB 2/15/09
      CALL HY36ENCODE(5,I,I36,ERRMSG,ERRMSG_LEN)
      IF (ERRMSG_LEN.NE.0) THEN
      WRITE(6,'(A,A)') ' COORIO-err', ERRMSG(1:ERRMSG_LEN) 
      CALL WRNDIE(-5,'COORIO','Fatal error in hybrid36 conversion.')
      END IF
C
      IF (FORM.EQ.'PDBO'.OR.FORM.EQ.'PDBA'.OR.FORM.EQ.'PDBH') THEN
      IF (FORM.EQ.'PDBO'.OR.FORM.EQ.'PDBA') THEN
      KEYWORD='ATOM  '
      ELSE
      KEYWORD='HETATM'
      END IF
C
C detect 2-letter element types by comparing NAME and CHEMical type.
      ELEM(1:1)=CHAR(ASCIIM(ICHAR(IAC(I)(1:1))))
      ELEM(2:2)=CHAR(ASCIIM(ICHAR(IAC(I)(2:2))))
C First, set ELCH for a single-letter element
      ELCH=' '//ELEM(1:1)
C The first two characters of the name and CHEMical type must match,
C case-insensitive.
      IF (ASCIIM(ICHAR(TYPE(I)(1:1))).EQ.ICHAR(ELEM(1:1)).AND.
     &    ASCIIM(ICHAR(TYPE(I)(2:2))).EQ.ICHAR(ELEM(2:2))) THEN
C they must also be a valid 2-letter element symbol.
      J=0
      FOUND=.FALSE.
      DO WHILE(J.LT.ELEM2_SIZE .AND. .NOT.FOUND)
      J=J+1
      FOUND = ELEM2_LIST(J).EQ.ELEM
      END DO
      IF (FOUND) THEN
C must also satisfy one of:
C 1) type and name are both exactly 2 letters
C 2) 2nd type character is lower case
C 3) type characters 3,4 indicate a formal charge [+-][0-9]
      IF (TYPE(I)(3:4).EQ.'  ' .OR.
     &    (IAC(I)(2:2).GE.'a'.AND.IAC(I)(3:3).LE.'z')) THEN
      ELCH = ELEM
      ELSE IF ((IAC(I)(3:3).EQ.'+'.OR.IAC(I)(3:3).EQ.'-') .AND.
     &         (IAC(I)(4:4).GE.'0'.AND.IAC(I)(4:4).LE.'9')) THEN
      ELCH = ELEM//IAC(I)(3:4)
      ELSE
      FOUND=.FALSE.
      END IF
      IF (FOUND .AND. NAME(1:1).EQ.' ') NAME = NAME(2:4)
      END IF
      END IF
C use segid as chainid if it is one character
      SEGTMP=SEGID(I)
      SEGLEN=4
      CALL TRIMM(SEGTMP,SEGLEN)
      CALL TRIML(SEGTMP,SEGLEN)
      IF (SEGLEN.EQ.1) THEN
      WRITE(UNIT,
     & '(A,A,1X,A4,1X,A4,A1,A5,3X,3F8.3,2F6.2,6X,2A4)') KEYWORD,
     & I36,NAME,RES(I),SEGTMP(1:1),A,
     & X(I),Y(I),Z(I),QMAIN(I),WMAIN(I),SEGID(I),ELCH
      ELSE
      WRITE(UNIT,
     & '(A,A,1X,A4,1X,A4,A1,A5,3X,3F8.3,2F6.2,6X,2A4)') KEYWORD,
     & I36,NAME,RES(I),' ',A,
     & X(I),Y(I),Z(I),QMAIN(I),WMAIN(I),SEGID(I),ELCH
      END IF
      ELSE
C old-style CNS record
      WRITE(UNIT,
     & '(A,A,1X,A4,1X,A4,1X,A5,3X,3F8.3,2F6.2,6X,A4)') 'ATOM  ',
     & I36,NAME,RES(I),A,
     & X(I),Y(I),Z(I),QMAIN(I),WMAIN(I),SEGID(I)
      END IF
C
      END IF
      END DO
      IF (FORM(1:3).EQ.'CNS'.OR.FORM(1:4).EQ.'PDBO') THEN
      WRITE(UNIT,'(A)') 'END'
      CALL VCLOSE(UNIT,'KEEP',QERROR)
      END IF
C
      END IF
      END IF
C
      RETURN
      END
C
      SUBROUTINE CREAD(X,Y,Z,WMAIN,QMAIN,FLAGS,QCHAIN)
C
C reads coordinates.
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      DOUBLE PRECISION X(*), Y(*), Z(*), WMAIN(*), QMAIN(*)
      INTEGER FLAGS(*)
      LOGICAL QCHAIN
C local
      LOGICAL COND, ECHOLD, QATOM
      INTEGER NMISS, NMULT, I, IATOM, IOFFS
      DOUBLE PRECISION XIN, YIN, ZIN, WIN, QIN
      CHARACTER*4 SID, RID, REN, IUP
      CHARACTER*5 A
      CHARACTER*1 CHAIN
      INTEGER NUMBER, IJ, CURRSTRM
C parameter
      INTEGER MARK
      DOUBLE PRECISION ONE
      PARAMETER (MARK=-9999, ONE=1.0D0)
C begin
C
C initialize various counters
      NMISS=0
      NMULT=0
      IOFFS=0
      ECHOLD=QECHO
      USED=.TRUE.
      NUMBER=-1
C
      CURRSTRM=NSTRM
C
      DO WHILE (USED)
      QATOM=.FALSE.
C
      IF (WD(1:4).EQ.'ATOM'.OR.WD(1:4).EQ.'HETA') THEN
C ==================================================
C this section reads the Brookhaven atom coordinates
C ==================================================
      READ(COMLYN,'(12X,A4,1X,A4,A1,A5,3X,3F8.3,F6.2,F6.2,6X,A4)',
     &     ERR=8888) IUP,REN,CHAIN,A,XIN,YIN,ZIN,QIN,WIN,SID
C
C special handling of PDB RESSEQ and ICODE, ATB 02/15/10
C
C if the ICODE character is non-blank, then incorporate
C into the RESID as the right-most character, but only if
C   the RESSEQ has only three characters. 
      IF (A(5:5).EQ.' ') then
         RID=A(1:4)
      ELSE
         IF (A(1:1).NE.' ') then
           WRITE(6,'(A,A,A)')     
     & ' %CREAD-ERR: residue ID and insertion character ',
     & A,' exceed 4 characters.'           
           CALL WRNDIE(-5,'COORIO',
     &      'Unsupported PDB fields.')
         END IF
         RID=A(2:5) 
      END IF
C
C make RID left-justified
      IJ=4
      CALL TRIML(RID,IJ)
      CALL TRIMM(RID,IJ)
C make RESID left-justified and remove trailing blanks
CCC      IJ=5
CCC      CALL TRIML(A,IJ)
CCC      CALL TRIMM(A,IJ)
CCC      IF (IJ.GT.4) THEN
CCC      WRITE(6,'(A,A,A)')
CCC     & ' %CREAD-ERR: residue ID and insertion character ',
CCC     & A,' exceed 4 characters.'
CCC      END IF
CCC      RID=A(1:4)
C
C make IUP left-justified
      IJ=4
      CALL TRIML(IUP,IJ)
C
C make REN left-justified
      IJ=4
      CALL TRIML(REN,IJ)
C
      IF (QCHAIN) THEN
         IF (CHAIN.NE.' ') THEN
            SID=CHAIN
         END IF
      END IF
C
      QATOM=.TRUE.
C
      ELSE IF (WD(1:3).EQ.'TER') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'BREAK') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'CRYST') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:4).EQ.'ORIG') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'SCALE') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:4).EQ.'REMA') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'HEADER') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'ANISOU') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'COMPND') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'TITLE') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'SOURCE') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'KEYWDS') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'EXPDTA') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'DBREF') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'AUTHOR') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'REVDAT') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:4).EQ.'JRNL') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'SEQADV') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'SEQRES') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'MODRES') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'FTNOTE') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'HET   ') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'HETNAM') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'HETSYN') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'CAVEAT') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'FORMUL') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'HELIX') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'SHEET') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:4).EQ.'TURN') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'SSBOND') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'SPLIT') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:4).EQ.'SITE') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:4).EQ.'LINK') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'MTRIX') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'CISPEP') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'CONECT') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'MASTER') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'SPRSDE') THEN
CC      QECHO=.TRUE.
      CURSOR=COMLEN
C
      ELSE
      USED=.FALSE.
C
C check for missing end statement, ATB 2/21/10
      IF (CURRSTRM.NE.NSTRM.AND.WD(1:WDLEN).NE.'END') THEN
C
      WRITE(6,'(A)') 
     & ' COORIO-err: End of coordinate file without END statement.'
      WRITE(6,'(A)') 
     & '             Please insert END statement and re-run.'
      CALL WRNDIE(-5,'CREAD','Fatal error in reading coordinates.')
      END IF
C
      END IF
C
      IF (QATOM) THEN
C
C ======================================
C try to find atom in molecular topology
C ======================================
      IOFFS=IOFFS+1
C
C make RID right-justified
      IJ=4
      CALL TRIML(RID,IJ)
C
C get atom number (note the shortcut for coordinate files
C with correct sequence !!)
      COND=SID.EQ.SEGID(IOFFS).AND.RID.EQ.RESID(IOFFS).AND.
     &     REN.EQ.RES(IOFFS).AND.IUP.EQ.TYPE(IOFFS)
      IF (COND) THEN
      IATOM=IOFFS
      ELSE
      IATOM=GETATN(SID,RID,REN,IUP,MARK)
C
CCCC modification ATB 4/27/08
      IF (IATOM.NE.MARK) IOFFS=IATOM
      END IF
      IF (IATOM.EQ.MARK) THEN
      WRITE(6,'(9A)') ' %READC-ERR: atom ',SID,' ',RID,' ',
     &      REN,' ',IUP,' not found in molecular structure'
      ELSE
      IF (FLAGS(IATOM).EQ.1) THEN
      IF (INITIA(IATOM,X,Y,Z)) NMULT=NMULT+1
      X(IATOM)=XIN
      Y(IATOM)=YIN
      Z(IATOM)=ZIN
      QMAIN(IATOM)=QIN
      WMAIN(IATOM)=WIN
      END IF
      END IF
C
      CURSOR=COMLEN
      END IF
C
      IF (USED) THEN
      CALL NEXTWD('COOR>')
      QECHO=.FALSE.
      END IF
C
      END DO
C
C check multiple coordinate assignments
      IF (NMULT.GT.0.AND.WRNLEV.GE.5) THEN
      WRITE(6,'(A,I6,A)') ' %READC-WRN: multiple coordinates for ',
     &     NMULT,' atoms'
      END IF
C
C check still undefined atoms
      IF (IOFFS.GT.0) THEN
      NMISS=0
      DO I=1,NATOM
      IF (.NOT.INITIA(I,X,Y,Z).AND.FLAGS(I).EQ.1) NMISS=NMISS+1
      END DO
      IF (NMISS.GT.0.AND.WRNLEV.GE.5) THEN
      WRITE(6,'(A,I6,A)') ' %READC-WRN: still ',NMISS,
     &              ' missing coordinates (in selected subset)'
      END IF
      END IF
C
C-ERROR-LABEL
      GOTO 7777
8888  ERROR=.TRUE.
      WRITE(6,'(A)') ' %COOR-ERR: ERROR during reading coordinates'
7777  CONTINUE
C
      QECHO=ECHOLD
      RETURN
      END
C
