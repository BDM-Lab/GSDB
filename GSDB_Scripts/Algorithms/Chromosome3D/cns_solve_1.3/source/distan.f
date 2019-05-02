C======================================================================
      SUBROUTINE DISTAN
C
C Nonbonded distance analysis routine.  This routine
C analyses the nonbonded distances as provided by the
C nonbonded list routines.  Selected parts of the nonbonded
C list may be printed by specifying a upper and lower cutoff
C and atom selections.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'pick.inc'
C begin
      IISEL=ALLHP(INTEG4(NATOM))
      IJSEL=ALLHP(INTEG4(NATOM))
      CALL DISTA2(HEAP(IISEL),HEAP(IJSEL))
      CALL FREHP(IJSEL,INTEG4(NATOM))
      CALL FREHP(IISEL,INTEG4(NATOM))
      RETURN
      END
C======================================================================
      SUBROUTINE DISTA2(ISEL,JSEL)
C
C See routine DISTAN above
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'pick.inc'
      INTEGER ISEL(*), JSEL(*)
C local
      DOUBLE PRECISION EVDW,ELEC,EPVDW,EPELEC,VEVDW,VELEC,VEPVDW,VEPELE
      INTEGER I, N
      LOGICAL ERR
C parameter
      DOUBLE PRECISION ZERO, LARGE
      PARAMETER (ZERO=0.0D0, LARGE=9999.D0)
C begin
C
C defaults
C
      CALL FILL4(ISEL,NATOM,1)
      NISEL=NATOM
      CALL FILL4(JSEL,NATOM,1)
      NJSEL=NATOM
      RCUTON=ZERO
      RCUTOF=LARGE
      OFILE='OUTPUT'
      SDISP='PRIN'
      RESLOW=-1
      RESDELTA=-1
      NDISTANCE=-1
C
C parsing
      CALL PUSEND('DISTANCE>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('DISTANCE>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-distance')
C
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' --------nonbonded-distance-analysis---',
     & '--------------------'
      WRITE(6,'(A,I6,A,I6)')
     & ' No. of selected atoms: "FROM" set=',NISEL,' "TO" set',NJSEL
      WRITE(6,'(A,F8.3,A,F8.3)')
     & ' CUTOn=',RCUTON,' CUTOff=',RCUTOF
      IF (RESDELTA.GT.0) THEN
      WRITE(6,'(A,I10)')
     & ' RESDelta=',RESDELTA
      END IF
      IF (RESLOW.GT.0) THEN
      WRITE(6,'(A,I10)')
     & ' RESLow=',RESLOW
      END IF
      IF (NDISTANCE.GT.0) THEN
      WRITE(6,'(A,I10)')
     & ' NDIStance=',NDISTANCE
      END IF
      WRITE(6,'(2A)')
     & ' DISPosition=',SDISP
      WRITE(6,'(2A)') ' --------------------------------------',
     & '--------------------'
      ELSE IF (WD(1:4).EQ.'FROM') THEN
      CALL SELCTA(ISEL,NISEL,X,Y,Z,.TRUE.)
      ELSE IF (WD(1:2).EQ.'TO') THEN
      CALL SELCTA(JSEL,NJSEL,X,Y,Z,.TRUE.)
      ELSE IF (WD(1:5).EQ.'CUTON') THEN
      CALL NEXTF('CUTON=',RCUTON)
      ELSE IF (WD(1:5).EQ.'CUTOF') THEN
      CALL NEXTF('CUTOFf=',RCUTOF)
      ELSE IF (WD(1:4).EQ.'RESD') THEN
      CALL NEXTI('RESDelta=',RESDELTA)
      ELSE IF (WD(1:4).EQ.'RESL') THEN
      CALL NEXTI('RESLow=',RESLOW)
      ELSE IF (WD(1:4).EQ.'NDIS') THEN
      CALL NEXTI('NDIStance=',NDISTANCE)
      ELSE IF (WD(1:4).EQ.'DISP') THEN
      CALL NEXTA4('disposition=',SDISP)
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
      ELSE
      CALL CHKEND('DISTANCE>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      CALL ASSFIL(OFILE,RUNIT,'WRITE','FORMATTED',ERROR)
C
      IF (NISEL.GT.0.AND.NJSEL.GT.0) THEN
C
C write some info
      IF (SDISP.EQ.'PRIN') THEN
      WRITE(6,'(A)')' DISTAN: nonbonded distances printed'
      ELSEIF (SDISP.EQ.'DEN') THEN
      WRITE(6,'(A)')' DISTAN: DEN distance assignments generated'
      ELSE IF (SDISP.EQ.'RMSD') THEN
      WRITE(6,'(A)')
     &  ' DISTAN: minimum nonbonded distances copied to RMSD array'
      ELSE IF (SDISP.EQ.'MATR') THEN
      WRITE(6,'(A)') ' DISTAN: distances written into matrix '
      END IF
C
C initialize RMSD array if required
      IF (SDISP.EQ.'RMSD') THEN
      DO I=1,NATOM
      IF (ISEL(I).EQ.1) THEN
      RMSD(I)=LARGE
      ELSE
      RMSD(I)=ZERO
      END IF
      END DO
      END IF
C
C allocate space for SMAT and initialize if required
      IF (SDISP.EQ.'MATR') THEN
      ISMAT=ALLHP(IREAL8(NISEL*NJSEL))
      IMAPI=ALLHP(INTEG4(NATOM))
      IMAPJ=ALLHP(INTEG4(NATOM))
      CALL DISTA5('INIT',NISEL,NJSEL,ISEL,JSEL,HEAP(ISMAT),RUNIT,NATOM,
     &            HEAP(IMAPI),HEAP(IMAPJ))
      END IF
C
C loop over all active Pairs of Interacting Groups (PIGs)
      DO N=1,NPIG
C
C check for unknown coordinates
      CALL ATMCHK(HEAP(IINTER(N)),ERR)
      IF (.NOT.ERR) THEN
C
C call energy function with analysis flag set
      CALL ENBOND(EVDW,ELEC,EPVDW,EPELEC,VEVDW,VELEC,VEPVDW,VEPELE,
     &            N,NPIG,'ANAL')
      END IF
      END DO
C
C modify RMSD array if required
      IF (SDISP.EQ.'RMSD') THEN
      DO I=1,NATOM
      IF (ISEL(I).EQ.1.AND.RMSD(I).GE.LARGE-RSMALL) THEN
      RMSD(I)=ZERO
      END IF
      END DO
      END IF
C
C print SMAT and de-allocate space if required
      IF (SDISP.EQ.'MATR') THEN
      CALL DISTA5('PRIN',NISEL,NJSEL,ISEL,JSEL,HEAP(ISMAT),RUNIT,NATOM,
     &            HEAP(IMAPI),HEAP(IMAPJ))
      CALL FREHP(IMAPI,INTEG4(NATOM))
      CALL FREHP(IMAPJ,INTEG4(NATOM))
      CALL FREHP(ISMAT,IREAL8(NISEL*NJSEL))
      END IF
C
      END IF
C
      RETURN
      END
C====================================================================
      SUBROUTINE ENBANL(EVDW,ELEC,EVDWV,ELECV,NATOM,JNB,
     &          LIST,CG,
     &          MAXCN,CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,
     &          EVWGHT,VVWGHT,QIMAG,QFORC,XJ,YJ,ZJ,JND,CGL,
     &          CNBAL,CNBBL,CNBIL,XJL,YJL,ZJL,JSYM)
C
C nonbonded distance analysis routine
C
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'funct.inc'
      DOUBLE PRECISION EVDW, ELEC
      DOUBLE PRECISION EVDWV, ELECV
      DOUBLE PRECISION EVWGHT, VVWGHT
      INTEGER NATOM, JNB(*), LIST(4)
      DOUBLE PRECISION CG(*)
      INTEGER MAXCN
      DOUBLE PRECISION CNBA(MAXCN,MAXCN), CNBB(MAXCN,MAXCN)
      DOUBLE PRECISION CNBVR(MAXCN,MAXCN)
      INTEGER LOOKUP(*)
      DOUBLE PRECISION ESCALE, VSCALE
      LOGICAL QIMAG, QFORC
      DOUBLE PRECISION XJ(*), YJ(*), ZJ(*)
      INTEGER JND(*)
      DOUBLE PRECISION CGL(*), CNBAL(*), CNBBL(*), CNBIL(*), XJL(*)
      DOUBLE PRECISION YJL(*), ZJL(*)
      INTEGER JSYM
C pointer
      INTEGER PINTER, PERM, RINTER, RESNUM
C local
      INTEGER NINTER, I
C 
C determine total number of interactions
      NINTER=0
      DO I=1,NATOM
      NINTER=NINTER+JNB(I)
      END DO
      PINTER=ALLHP(INTEG4(NINTER))
      PERM=ALLHP(INTEG4(NINTER))
      RINTER=ALLHP(IREAL8(NINTER))
      RESNUM=ALLHP(INTEG4(NATOM))
C
      CALL ENBAN2(SDISP,JNB,LIST,HEAP(IISEL),HEAP(IJSEL),
     &            RCUTON,RCUTOF,
     &            X,Y,Z,RMSD,RUNIT,
     &            NISEL,NJSEL,HEAP(IMAPI),HEAP(IMAPJ),HEAP(ISMAT),
     &            QIMAG,LNCSST,
     &            JSYM,XJ,YJ,ZJ,XJL,YJL,ZJL,JND,RESDELTA,RESLOW,
     &            NDISTANCE,NINTER,HEAP(PINTER),HEAP(PERM),
     &            HEAP(RINTER),HEAP(RESNUM))
C
      CALL FREHP(RESNUM,INTEG4(NATOM))
      CALL FREHP(RINTER,IREAL8(NINTER))
      CALL FREHP(PERM,INTEG4(NINTER))
      CALL FREHP(PINTER,INTEG4(NINTER))
C
      RETURN
      END
C======================================================================
      SUBROUTINE ENBAN2(SDISP,JNB,LIST,ISEL,JSEL,RCUTON,RCUTOF,
     &                  X,Y,Z,RMSD,RUNIT,
     &                  NISEL,NJSEL,MAPI,MAPJ,SMAT,QIMAG,LNCSST,
     &                  ISYM,XJ,YJ,ZJ,XJL,YJL,ZJL,JND,RESDELTA,
     &                  RESLOW,
     &                  NDISTANCE,NINTER,PINTER,PERM,RINTER,RESNUM)
C
C See routine ENBANL above
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      CHARACTER*4 SDISP
      INTEGER JNB(*), LIST(4), ISEL(*), JSEL(*)
      DOUBLE PRECISION RCUTON, RCUTOF, X(*),Y(*),Z(*), RMSD(*)
      INTEGER RUNIT
      INTEGER NISEL, NJSEL, MAPI(*), MAPJ(*)
      DOUBLE PRECISION SMAT(NISEL,NJSEL)
      LOGICAL QIMAG, LNCSST
      INTEGER ISYM
      DOUBLE PRECISION XJ(*), YJ(*), ZJ(*), XJL(*), YJL(*), ZJL(*)
      INTEGER JND(*), RESDELTA, RESLOW
      INTEGER NDISTANCE, NINTER, PINTER(*), PERM(*)
      DOUBLE PRECISION RINTER(*)
      INTEGER RESNUM(*)
C local
      DOUBLE PRECISION RANDOM
      INTEGER I, N, J, JJ, MINTER, LINTER
      CHARACTER*12 SOPS
      CHARACTER*1 ICODE
      LOGICAL OK, COND
C parameter
C begin
C
C initialize the string for the printed output line (SOPS=12chars)
      IF (QIMAG) THEN
      WRITE(SOPS,'(A,I3,A)') '"(XSYM#',ISYM,') '
      ELSE
      IF (LNCSST.AND.ISYM.GT.1) THEN
      WRITE(SOPS,'(A,I3,A)') '"(NCS#',ISYM,')  '
      ELSE
      WRITE(SOPS,'(A)') '"           '
      END IF
      END IF
C
C generate residue integers from RESId character string
      DO I=1,NATOM
      CALL SPLITI(RESID(I),RESNUM(I),ICODE,4,OK)
      IF (.NOT.OK) THEN
      CALL WRNDIE(-5,'DISTAN',
     & ' Non-standard resids (must be integer or hybrid36).')
      GOTO 9999
      END IF
      END DO
C
C
C if NDISTANCE specified then generate an array that sspecifies
C a random selection of exactly NDISTANCE interactions
      IF (NDISTANCE.GT.0) THEN
C
C do dry run
C main loop over atom I
      LINTER=0
      DO I=1,NATOM
C
      N=JNB(I)
      IF (N.GT.0) THEN
C
C define local index list JND
      CALL GETNB(N,JND,LIST)
C
C gather "j" coordinates
      DO J=1,N
      XJL(J)=XJ(JND(J))
      YJL(J)=YJ(JND(J))
      ZJL(J)=ZJ(JND(J))
      END DO
C
C compute differences to X(I), Y(I), Z(I)
      IF (QIMAG) THEN
      CALL XPIMAG(X(I),Y(I),Z(I),1,N,XJL,YJL,ZJL,XJL,YJL,ZJL,.FALSE.,0)
      ELSE
      DO J=1,N
      XJL(J)=X(I)-XJL(J)
      YJL(J)=Y(I)-YJL(J)
      ZJL(J)=Z(I)-ZJL(J)
      END DO
      END IF
C
C compute distance
      DO J=1,N
      XJL(J)=SQRT(XJL(J)**2+YJL(J)**2+ZJL(J)**2)
      END DO
C
C dispose distances
      DO J=1,N
      JJ=JND(J)
C
C RCUTOF criterion
      IF (XJL(J).GE.RCUTON.AND.XJL(J).LE.RCUTOF) THEN
C
C residue-delta criterion 
      IF (RESDELTA.LT.0.OR.
     & (ABS(RESNUM(I)-RESNUM(JJ)).LE.RESDELTA
     & .AND.
     &  SEGID(I).EQ.SEGID(JJ))) THEN
C
C residue-low criterion 
      IF (RESLOW.LT.0.OR.
     & (ABS(RESNUM(I)-RESNUM(JJ)).GE.RESLOW
     & .AND.
     &  SEGID(I).EQ.SEGID(JJ))) THEN
C 
C selection criterion
      IF ((JSEL(I).EQ.1.AND.ISEL(JJ).EQ.1).OR.
     &    (ISEL(I).EQ.1.AND.JSEL(JJ).EQ.1)) THEN
      LINTER=LINTER+1
      END IF
      END IF
      END IF
      END IF
      END DO
      END IF
      END DO
C
      NDISTANCE=MIN(NDISTANCE,LINTER)
C
C end do dry run
C
C fill array with random numbers between 0 and 1
      DO I=1,LINTER
      CALL GGUBFS(RANDOM)
      RINTER(I)=RANDOM
      END DO
C
C sort array
      CALL SORTP(LINTER,PERM,ORDERR,RINTER,1,0,0,0,0,0,0)
      DO I=1,LINTER
      PINTER(I)=0
      END DO
C
C pick NDISTANCE first indices and set flags to 1
      DO I=1,NDISTANCE
      PINTER(PERM(I))=1
      END DO
C
C reset nonbonded list pointer
      CALL RESNB(LIST)
      END IF
CCC
C
C main loop over atom I
      MINTER=0
      DO I=1,NATOM
C
      N=JNB(I)
      IF (N.GT.0) THEN
C
C define local index list JND
      CALL GETNB(N,JND,LIST)
C
C gather "j" coordinates
      DO J=1,N
      XJL(J)=XJ(JND(J))
      YJL(J)=YJ(JND(J))
      ZJL(J)=ZJ(JND(J))
      END DO
C
C compute differences to X(I), Y(I), Z(I)
      IF (QIMAG) THEN
      CALL XPIMAG(X(I),Y(I),Z(I),1,N,XJL,YJL,ZJL,XJL,YJL,ZJL,.FALSE.,0)
      ELSE
      DO J=1,N
      XJL(J)=X(I)-XJL(J)
      YJL(J)=Y(I)-YJL(J)
      ZJL(J)=Z(I)-ZJL(J)
      END DO
      END IF
C
C compute distance
      DO J=1,N
      XJL(J)=SQRT(XJL(J)**2+YJL(J)**2+ZJL(J)**2)
      END DO
C
C dispose distances
      DO J=1,N
      JJ=JND(J)
C
C RCUTOF criterion
      IF (XJL(J).GE.RCUTON.AND.XJL(J).LE.RCUTOF) THEN
C
C residue-delta criterion 
      IF (RESDELTA.LT.0.OR.
     & (ABS(RESNUM(I)-RESNUM(JJ)).LE.RESDELTA
     & .AND.
     &  SEGID(I).EQ.SEGID(JJ))) THEN
C
C residue-low criterion 
      IF (RESLOW.LT.0.OR.
     & (ABS(RESNUM(I)-RESNUM(JJ)).GE.RESLOW
     & .AND.
     &  SEGID(I).EQ.SEGID(JJ))) THEN
C
C selection criterion
      IF ((JSEL(I).EQ.1.AND.ISEL(JJ).EQ.1).OR.
     &    (ISEL(I).EQ.1.AND.JSEL(JJ).EQ.1)) THEN
C
      COND=.TRUE.
C ndistance criterion
      IF (NDISTANCE.GT.0) THEN
      MINTER=MINTER+1
      IF (PINTER(MINTER).EQ.0) COND=.FALSE.
      END IF
C
      IF (COND) THEN
C
      IF (SDISP.EQ.'RMSD') THEN
      IF (RMSD(JJ).GT.XJL(J)) RMSD(JJ)=XJL(J)
      IF (RMSD(I).GT.XJL(J)) RMSD(I)=XJL(J)
C
      ELSE IF (SDISP.EQ.'PRIN') THEN
      WRITE(RUNIT,'(17A,F7.4,A)')
     &  ' atoms "', SEGID(I),'-',RESID(I),
     &  '-',RES(I),'-',TYPE(I),
     & '" and "',SEGID(JJ),'-',RESID(JJ),
     & '-',RES(JJ),'-',TYPE(JJ),
     & SOPS,XJL(J),' A apart'
C    
      ELSE IF (SDISP.EQ.'DEN') THEN
      WRITE(RUNIT,'(13A,F7.4,A)')
     &  ' ASSIgn ( KNOWn and ATOM "', SEGID(I),
     &          '" "',RESID(I),'" "',
     &  TYPE(I),'" ) ( KNOWn and ATOM "',SEGID(JJ),
     &          '" "',RESID(JJ),'" "',
     &  TYPE(JJ),'" ) ',XJL(J),' 0. 0.'
C
      ELSE IF (SDISP.EQ.'MATR') THEN
      IF (JSEL(I).EQ.1.AND.ISEL(JJ).EQ.1) THEN
      SMAT(MAPI(JJ),MAPJ(I))=XJL(J)
      ELSE IF (ISEL(I).EQ.1.AND.JSEL(JJ).EQ.1) THEN      
      SMAT(MAPI(I),MAPJ(JJ))=XJL(J)
      END IF
C
      END IF
      END IF
      END IF
      END IF
      END IF
      END IF
      END DO
C
      END IF
C
C end of main loop
      END DO
9999  CONTINUE
      RETURN
      END
C======================================================================
      SUBROUTINE DISTA5(MODE,NISEL,NJSEL,ISEL,JSEL,SMAT,RUNIT,NATOM,
     &                  MAPI,MAPJ)
C
C initializes and prints SMAT distance matrix
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      CHARACTER*4 MODE
      INTEGER NISEL, NJSEL, ISEL(*), JSEL(*)
      DOUBLE PRECISION SMAT(NISEL,NJSEL)
      INTEGER RUNIT, NATOM, MAPI(*), MAPJ(*)
C local
      CHARACTER*4 ISID,IRID,IREN,IAT,JSID,JRID,JREN,JAT
      INTEGER I, J, IM
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      IF (MODE.EQ.'INIT') THEN
      DO I=1,NISEL
      DO J=1,NJSEL
      SMAT(I,J)=ZERO
      END DO
      END DO
C
C define map functions
      IM=0
      DO I=1,NATOM
      IF (ISEL(I).EQ.1) THEN
      IM=IM+1
      MAPI(I)=IM
      ELSE
      MAPI(I)=0
      END IF
      END DO
C
      IM=0
      DO I=1,NATOM
      IF (JSEL(I).EQ.1) THEN
      IM=IM+1
      MAPJ(I)=IM
      ELSE
      MAPJ(I)=0
      END IF
      END DO
C
      ELSE
      CALL WRTITL(RUNIT,0)
      WRITE(RUNIT,'(A)') ' '
      WRITE(RUNIT,'(A)') ' DISTANCE MATRIX'
C
      WRITE(RUNIT,'(A,I10)') ' NSET1=',NISEL
      DO I=1,NATOM
      IF (ISEL(I).EQ.1) THEN
      CALL ATOMID(I,ISID,IRID,IREN,IAT)
      WRITE(RUNIT,'(1X,I6,9A)') I,'=(',ISID,' ',IRID,' ',IREN,
     &      ' ',IAT,')'
      END IF
      END DO
      WRITE(RUNIT,'(A)') ' '
C
      WRITE(RUNIT,'(A,I10)') ' NSET2=',NJSEL
      DO J=1,NATOM
      IF (JSEL(J).EQ.1) THEN
      CALL ATOMID(J,JSID,JRID,JREN,JAT)
      WRITE(RUNIT,'(1X,I6,9A)') J,'=(',JSID,' ',JRID,' ',JREN,
     &      ' ',JAT,')'
      END IF
      END DO
      WRITE(RUNIT,'(A)') ' '
C
C
      DO I=1,NISEL
      DO J=1,NJSEL
      WRITE(RUNIT,'(1X,G15.8)') SMAT(I,J)
      END DO
      WRITE(RUNIT,'(A)')
      END DO
C
      END IF
      RETURN
      END
