      SUBROUTINE ARISET
C
C Sets up parameters for ARIA assignment, calibration and
C violation analysis.
C Nilges, M., Macias, M., O'Donoghue, S., & Oschkinat, M. (1997).
C Automated NOESY Interpretation with Ambiguous Distance Restraints
C The Refined Solution Structure of the Pleckstrin Homology Domain
C from beta-Spectrin.
C J. Mol. Biol. 269, 408-422.
C
C Author: Michael Nilges
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'noe.inc'
C local
      INTEGER ISLCT, JSLCT
C begin
C set up some arrays of size NATOM
      IF (HPNHGL.EQ.0) THEN
      HPNHGL=ALLHP(INTEG4(NATOM))
      CALL FILL4(HEAP(HPNHGL),NATOM,1)
      CALL MAKIND(HEAP(HPNHGL),NATOM,ISLCT)
      END IF
      IF (HPNCAL.EQ.0) THEN
      HPNCAL=ALLHP(INTEG4(NATOM))
      CALL FILL4(HEAP(HPNCAL),NATOM,0)
      END IF
      IF (HPNPSE.EQ.0) THEN
      HPNPSE=ALLHP(IREAL8(NATOM))
      CALL FILLR8(HEAP(HPNPSE),NATOM,0.0D0)
      END IF
      IF (HPNDIA.EQ.0) THEN
      HPNDIA=ALLHP(IREAL8(NATOM))
      CALL FILLR8(HEAP(HPNDIA),NATOM,1.0D0)
      END IF
C
      ISLCT=ALLHP(INTEG4(NATOM))
      JSLCT=ALLHP(INTEG4(NATOM))
      CALL ARISE2(HEAP(ISLCT),HEAP(JSLCT),HEAP(HPNHGL),
     &            HEAP(HPNCAL),HEAP(HPNPSE),HEAP(HPNDIA))
      CALL FREHP(JSLCT,INTEG4(NATOM))
      CALL FREHP(ISLCT,INTEG4(NATOM))
      RETURN
      END
C====================================================================
      SUBROUTINE ARISE2(ISLCT,JSLCT,NOEHGL,CALPTR,NOEPSE,DIAGON)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER ISLCT(*), JSLCT(*), NOEHGL(*), CALPTR(*)
      DOUBLE PRECISION NOEPSE(*), DIAGON(*)
C local
      INTEGER I, NISLCT, NMISS
      CHARACTER*4 CLASS, ACTION
      DOUBLE PRECISION VIOUPL, VIOLOL, TEMPX, TEMPY, TEMPZ
      DOUBLE PRECISION OLDNOE, NEWNOE
      LOGICAL MATCH
C parameter
      DOUBLE PRECISION ZERO, ONE, SIX
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, SIX=6.0D0)
C function
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
C all defaults are set in NOERES!
      CALL PUSEND('ARIA>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('ARIA>')
      CALL MISCOM('ARIA>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
         CALL CNSHELP('cns-aria')
C
C====================================================================
      ELSE IF (WD(1:2).EQ.'DO') THEN
      CALL ARIMAN(HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),
     &            HEAP(HPNWGH),HEAP(HPNVOL),HEAP(HPNCV))
C====================================================================
      ELSE IF (WD(1:4).EQ.'CFLI') THEN
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      IF (NISLCT.GT.0) THEN
      NMISS=0
      DO I=1,NATOM
      IF (ISLCT(I).NE.1) THEN
      ELSE IF (.NOT.INITIA(I,X,Y,Z)) THEN
      ISLCT(I)=0
      NMISS=NMISS+1
      NISLCT=NISLCT-1
      END IF
      END DO
      END IF
      IF (NISLCT.EQ.2) THEN
      CALL MAKIND(ISLCT,NATOM,NISLCT)
      TEMPX=X(ISLCT(1))
      TEMPY=Y(ISLCT(1))
      TEMPZ=Z(ISLCT(1))
      X(ISLCT(1))=X(ISLCT(2))
      Y(ISLCT(1))=Y(ISLCT(2))
      Z(ISLCT(1))=Z(ISLCT(2))
      X(ISLCT(2))=TEMPX
      Y(ISLCT(2))=TEMPY
      Z(ISLCT(2))=TEMPZ
      TEMPX=XV(ISLCT(1))
      TEMPY=YV(ISLCT(1))
      TEMPZ=ZV(ISLCT(1))
      XV(ISLCT(1))=XV(ISLCT(2))
      YV(ISLCT(1))=YV(ISLCT(2))
      ZV(ISLCT(1))=ZV(ISLCT(2))
      XV(ISLCT(2))=TEMPX
      YV(ISLCT(2))=TEMPY
      ZV(ISLCT(2))=TEMPZ
      END IF
C
      ELSE IF (WD(1:4).EQ.'VFLI') THEN
C
      IF (NISLCT.GT.0) THEN
      NMISS=0
      DO I=1,NATOM
      IF (ISLCT(I).NE.1) THEN
      ELSE IF (.NOT.INITIA(I,X,Y,Z)) THEN
      ISLCT(I)=0
      NMISS=NMISS+1
      NISLCT=NISLCT-1
      END IF
      END DO
      END IF
      IF (NISLCT.EQ.2) THEN
      CALL MAKIND(ISLCT,NATOM,NISLCT)
      TEMPX=XV(ISLCT(1))
      TEMPY=YV(ISLCT(1))
      TEMPZ=ZV(ISLCT(1))
      XV(ISLCT(1))=XV(ISLCT(2))
      YV(ISLCT(1))=YV(ISLCT(2))
      ZV(ISLCT(1))=ZV(ISLCT(2))
      XV(ISLCT(2))=TEMPX
      YV(ISLCT(2))=TEMPY
      ZV(ISLCT(2))=TEMPZ
      END IF
C
C====================================================================
      ELSE IF (WD(1:4).EQ.'FLIP') THEN
C
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      IF (NISLCT.GT.0) THEN
      NMISS=0
      DO I=1,NATOM
      IF (ISLCT(I).NE.1) THEN
      ELSE IF (.NOT.INITIA(I,X,Y,Z)) THEN
      ISLCT(I)=0
      NMISS=NMISS+1
      NISLCT=NISLCT-1
      ELSE IF (.NOT.FRSTEL(I,NOEHGL)) THEN
      ISLCT(I)=0
      NMISS=NMISS+1
      NISLCT=NISLCT-1
      END IF
      END DO
      END IF
C
      IF (NISLCT.EQ.2) THEN
      CALL MAKIND(ISLCT,NATOM,NISLCT)
      CALL ENOE(OLDNOE,'NOAN',0)
      WRITE(6,'(A,F15.5)') ' ARIFLI: ', OLDNOE
      CALL ARIFLI(ISLCT,JSLCT,HEAP(HPNILS),HEAP(HPNJLS),NOEHGL)
      CALL ENOE(NEWNOE,'NOAN',0)
      WRITE(6,'(A,F15.5)') ' ARIFLI: ', NEWNOE
      IF (OLDNOE.GT.NEWNOE) THEN
      CALL ARIFLI(ISLCT,JSLCT,HEAP(HPNILS),HEAP(HPNJLS),NOEHGL)
      END IF
      END IF
C
C===================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(A,I6,A,I6,A,/,A,F8.3,A,I6)')
     & ' NOE: total number of restraints:',NOENUM,
     & ' partitioned into ',NOECCN,' classes',
     & ' NOE: ceiling=',NOECEI,' current allocation=',NOEMAX
      IF (NOEICV.GT.0) THEN
      WRITE(6,'(A,/,A,I6)')
     & ' NOE: data are partitioned into working set and test set.',
     & ' NOE: test set number=',NOEICV
      END IF
C===================================================================
      ELSE IF (WD(1:4).EQ.'COUN') THEN
C
C get the class name:
      CALL PUSEND('COUNtviolations>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('COUNtviolations>')
C
      IF (WD(1:4).EQ.'HELP') THEN
         CALL CNSHELP('cns-aria-countviolations')
C
      ELSE IF (WD(1:4).EQ.'RESE') THEN
      VIOSTR = 0
      VIOEXC = ZERO
      VIOTHR = ZERO
      CALL NOEEXC(0,NOENUM,VIOEXC,VIOSTR,'INIT',ZERO,ZERO,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
C
      ELSE IF (WD(1:4).EQ.'THRE') THEN
      CALL NEXTF('THREshold=',VIOTHR)
      VIOSTR=VIOSTR+1
      CALL ARIVIO(NOENUM,VIOTHR,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
C
      ELSE IF (WD(1:4).EQ.'EXCL') THEN
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTF('EXCLude=',VIOEXC)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      CALL NOEEXC(I,NOENUM,VIOEXC,VIOSTR,'EXCL',ZERO,ZERO,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      END IF
      END DO
C
      ELSE IF (WD(1:3).EQ.'SET') THEN
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTF('ratio=',VIOEXC)
C get lower limit. check for TOKEN.
      CALL NEXTA4('lower_limit=',ACTION)
      IF (ACTION.EQ.'TOKE') THEN
      VIOLOL = -ONE
      ELSE
      CALL SAVEWD
      CALL NEXTF('lower_limit=',VIOLOL)
      END IF
C get upper limit. check for TOKEN.
      CALL NEXTA4('upper_limit=',ACTION)
      IF (ACTION.EQ.'TOKE') THEN
      VIOUPL = -ONE
      ELSE
      CALL SAVEWD
      CALL NEXTF('upper_limit=',VIOUPL)
      END IF
C
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      CALL NOEEXC(I,NOENUM,VIOEXC,VIOSTR,'SET ',VIOLOL,VIOUPL,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      END IF
      END DO
C
      ELSE IF (WD(1:4).EQ.'MULT') THEN
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTF('ratio=',VIOEXC)
C get lower limit. check for TOKEN.
      CALL NEXTA4('lower_limit=',ACTION)
      IF (ACTION.EQ.'TOKE') THEN
      VIOLOL = -ONE
      ELSE
      CALL SAVEWD
      CALL NEXTF('lower_limit=',VIOLOL)
      END IF
C get upper limit. check for TOKEN.
      CALL NEXTA4('upper_limit=',ACTION)
      IF (ACTION.EQ.'TOKE') THEN
      VIOUPL = -ONE
      ELSE
      CALL SAVEWD
      CALL NEXTF('upper_limit=',VIOUPL)
      END IF
C
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      CALL NOEEXC(I,NOENUM,VIOEXC,VIOSTR,'MULT',VIOLOL,VIOUPL,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      END IF
      END DO
C
      ELSE IF (WD(1:4).EQ.'WEIG') THEN
      CALL NEXTA4('class-name=',CLASS)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      CALL NOEEXC(I,NOENUM,VIOEXC,VIOSTR,'WEIG',VIOLOL,VIOUPL,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      END IF
      END DO
C
      ELSE IF (WD(1:4).EQ.'LIST') THEN
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTF('EXCLude=',VIOEXC)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      CALL NVIOLS(NOENUM,VIOEXC,VIOSTR,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      END IF
      END DO
C
      ELSE
      CALL CHKEND('COUNtviolations>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C====================================================================
      ELSE IF (WD(1:4).EQ.'CALI') THEN
      CALL CALIBR(ISLCT,CALPTR,DIAGON)
C====================================================================
      ELSE IF (WD(1:4).EQ.'ANAL') THEN
      CALL ANALRS(NOEHGL,ISLCT,JSLCT)
C====================================================================
      ELSE
      CALL CHKEND('ARIA>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
      RETURN
      END
C
C====================================================================
C
      SUBROUTINE ARIFLI(ISLCT,JSLCT,NOEILS,NOEJLS,NOEHGL)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'noe.inc'
      INTEGER NOEILS(*),NOEJLS(*)
      INTEGER NOEHGL(*),ISLCT(*),JSLCT(*)
C local
      INTEGER I,J,IAT,JAT,INOE
      LOGICAL QNEXT
C
C begin
C
      I=ISLCT(1)
      J=ISLCT(2)
      IAT=0
      QNEXT=.TRUE.
      DO WHILE (QNEXT)
        IAT=IAT+1
        ISLCT(IAT)=I
        CALL NEXTHY(I, NOEHGL,QNEXT)
      END DO
      JAT=0
      QNEXT=.TRUE.
      DO WHILE (QNEXT)
        JAT=JAT+1
        JSLCT(JAT)=J
        CALL NEXTHY(J, NOEHGL,QNEXT)
      END DO
      IF (IAT.EQ.JAT) THEN
      DO INOE=1,NOEMAX
        DO I=1,IAT
          IF (NOEILS(INOE).EQ.ISLCT(I)) THEN
            NOEILS(INOE)=JSLCT(I)
          ELSEIF (NOEILS(INOE).EQ.JSLCT(I)) THEN
            NOEILS(INOE)=ISLCT(I)
          END IF
          IF (NOEJLS(INOE).EQ.ISLCT(I)) THEN
            NOEJLS(INOE)=JSLCT(I)
          ELSEIF (NOEJLS(INOE).EQ.JSLCT(I)) THEN
            NOEJLS(INOE)=ISLCT(I)
          END IF
        END DO
      END DO
      END IF
C
      RETURN
      END
C
C====================================================================
C
      SUBROUTINE NEXTHY(IAT, HGLNK, QNEXT)
C
C returns next element in a group in IAT (cyclically).
C QNEXT true if IAT is an element on input
C
      IMPLICIT NONE
C
C global
C
C input/ output
      INTEGER IAT, HGLNK(*)
      LOGICAL QNEXT
C
C local
      INTEGER JAT
C
C begin
C
      JAT=HGLNK(IAT)
      QNEXT=(HGLNK(JAT).LT.JAT)
      IAT=JAT
C
      END
C
C======================================================================
C
      SUBROUTINE EQVGRP(NATOM,HGLNK)
C
C group hydrogen atoms into groups of unresolved hydrogens.
C
C Author: Michael Nilges, HHMI and Yale University
C
      IMPLICIT NONE
C global
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'coord.inc'
C input/output
      INTEGER NATOM,HGLNK(*)
C local
      INTEGER IISLCT,INHYBO,IPATBO,NSLCT,NHYBMX,IHGR
C
C begin
C
      NHYBMX=4
C
C initialise group link array
C
C allocate space for the bond lists
      IISLCT=ALLHP(INTEG4(NATOM))
      INHYBO=ALLHP(INTEG4(NATOM))
      IPATBO=ALLHP(INTEG4(NHYBMX*NATOM))
      IHGR=ALLHP(INTEG4(NATOM))
C
C parsing
      CALL PUSEND('EQUIvalent>')
      DO WHILE (.NOT. DONE)
      CALL NEXTWD('EQUIvalent>')
      CALL MISCOM('EQUIvalent>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
         CALL CNSHELP('cns-aria-equivalent')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'METH') THEN
      CALL ALLGRP(3,.FALSE.,HGLNK,HEAP(INHYBO),4,
     @       HEAP(IPATBO),HEAP(IISLCT),HEAP(IHGR))
C=====================================================================
      ELSE IF (WD(1:4).EQ.'1-2 ') THEN
      CALL ALLGRP(2,.FALSE.,HGLNK,HEAP(INHYBO),4,
     @       HEAP(IPATBO),HEAP(IISLCT),HEAP(IHGR))
C=====================================================================
      ELSE IF ((WD(1:4).EQ.'NONE').OR.(WD(1:4).EQ.'INIT')) THEN
      CALL HGRINI(NATOM,HGLNK)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(HEAP(IISLCT),NSLCT,X,Y,Z,.TRUE.)
      CALL NEWGRP(HEAP(IISLCT),NSLCT,.TRUE.,HGLNK,HEAP(IHGR))
C=====================================================================
      ELSE
      CALL CHKEND('EQUIvalent>',DONE)
C
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      CALL FREHP(IHGR,INTEG4(NATOM))
      CALL FREHP(IPATBO,INTEG4(NHYBMX*NATOM))
      CALL FREHP(INHYBO,INTEG4(NATOM))
      CALL FREHP(IISLCT,INTEG4(NATOM))
C
      RETURN
      END
C
C=======================================================================
C
      SUBROUTINE HGRINI(NATOM,HGLNK)
C
C initialises the hydrogen group array.
C Checks if atoms are known and interact.
C All subsequent hydrogen tests are based on the hydrogen group array.
C
C Author: Michael Nilges, HHMI and Yale University
C
      IMPLICIT NONE
C global
      INCLUDE 'cns.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
C input/output
      INTEGER NATOM,HGLNK(*)
C local
      INTEGER I,NH
C
C begin
C
      NH=0
      DO I=1,NATOM
      HGLNK(I)=0
      END DO
      DO I=1,NATOM
      IF (HYDROG(I)) THEN
      HGLNK(I)=I
      NH=NH+1
      END IF
      END DO
C
      RETURN
      END
C
C=======================================================================
C
      SUBROUTINE ALLGRP(HMIN,OVRWRT,HGLNK,NHYBON,NHYBMX,
     &                  PATBON,ISLCT,LHGR)
C
C defines hydrogen groups by connectivity and number.
C if OVRWRT any existing overlapping groups are overwritten.
C
C Author: Michael Nilges, HHMI and Yale University
C
      IMPLICIT NONE
C global
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
C input/output
      INTEGER HMIN,HGLNK(*),NHYBON(*)
      INTEGER NHYBMX,PATBON(NHYBMX,*),ISLCT(*),LHGR(*)
      LOGICAL OVRWRT
C local
      INTEGER NHYDRO,I,J,IBT,JBT,I1,J1
C function
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
C set up list of bonds to hydrogens
      DO I=1,NATOM
      NHYBON(I)=0
      END DO
C
      DO I=1,NBOND
      IBT=IB(I)
      JBT=JB(I)
      IF (IBT.GT.0 .AND. JBT.GT.0) THEN
      IF (HGLNK(JBT).GT.0) THEN
      NHYBON(IBT)=NHYBON(IBT)+1
      PATBON(NHYBON(IBT),IBT)=JBT
      END IF
      IF (HGLNK(IBT).GT.0) THEN
      NHYBON(JBT)=NHYBON(JBT)+1
      PATBON(NHYBON(JBT),JBT)=IBT
      END IF
      END IF
      END DO
C
C insert groups
      DO I=1,NATOM
C
      IF (HGLNK(I).EQ.0) THEN
      NHYDRO=NHYBON(I)
      DO I1=1,NATOM
      ISLCT(I1)=0
      END DO
C
      IF (NHYDRO.GE.HMIN) THEN
C found a group
      DO J=1,NHYDRO
      J1=PATBON(J,I)
      ISLCT(J1)=1
      END DO
      CALL NEWGRP(ISLCT,NHYDRO,OVRWRT,HGLNK,LHGR)
C
      ELSEIF (NHYDRO.GT.0) THEN
C one group for each h-atom
      DO J=1,NHYDRO
      J1=PATBON(J,I)
      ISLCT(J1)=1
      CALL NEWGRP(ISLCT,1,OVRWRT,HGLNK,LHGR)
      ISLCT(J1)=0
      END DO
C
      END IF
      END IF
      END DO
C
      RETURN
      END
C
C=======================================================================
C
      SUBROUTINE NEWGRP(ISLCT,NSLCT,OVRWRT,HGLNK,LGRP)
C
C groups are defined by array ISELCT.
C non-hydrogen atoms are removed.
C the groups are connected in a linked list:
C HGLNK(I1ST) = ILAST; => HGLNK(I1ST) > I1ST
C HGLNK(I2ND) = I1ST;  => HGLNK(I2ND) < I2ND
C HGLNK(I3RD) = I2ND;  ...
C ...
C HGLNK(I) = I  <=>  I is the only atom in group
C HGLNK(I) = 0  <=>  I is not a hydrogen
C if OVRWRT any existing overlapping groups are overwritten.
C
C Author: Michael Nilges, HHMI and Yale University
C
C global variables
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
C input/output
      INTEGER ISLCT(*),NSLCT,HGLNK(*),LGRP(*)
      LOGICAL OVRWRT
C local variables
      INTEGER NHYDRO,I,J,J1,JTEMP
      LOGICAL CLEAR
C function
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
C remove non-hydrogen and unknown or non-interacting atoms in selection
      DO I=1,NATOM
      IF (HGLNK(I).GT.0) THEN
      ISLCT(I)=ISLCT(I)
      ELSE
      ISLCT(I)=0
      END IF
      END DO
C gather selected group indices
      NHYDRO=0
      DO I=1,NATOM
      IF ((ISLCT(I)).GT.0) THEN
      NHYDRO=NHYDRO+1
      LGRP(NHYDRO)=I
      END IF
      END DO
C
      CLEAR=.TRUE.
C check if any groups exist and delete if OVRWRT
C
      DO J=1,NHYDRO
      J1=HGLNK(LGRP(J))
      IF (HGLNK(J1).NE.J1) THEN
C group
      IF (OVRWRT) THEN
C remove group
      DO WHILE ((HGLNK(J1).NE.J1))
      JTEMP=J1
      HGLNK(J1)=J1
      J1=HGLNK(JTEMP)
      END DO
      ELSE
      CLEAR=.FALSE.
      END IF
      END IF
      END DO
C
      IF (CLEAR) THEN
C add a new group
      IF (NHYDRO.GT.0) THEN
      DO J=1,NHYDRO
      IF (J.GT.1) THEN
      J1=J-1
      ELSE
      J1=NHYDRO
      END IF
      HGLNK(LGRP(J))=LGRP(J1)
      END DO
      END IF
      END IF
C
      RETURN
      END
C
C====================================================================
C
      LOGICAL FUNCTION FRSTEL(IAT, HGLNK)
C
C true if IAT is first element on input.
C
      IMPLICIT NONE
C
C input/ output
      INTEGER IAT, HGLNK(*)
C
C begin
C
      FRSTEL=((HGLNK(IAT).GE.IAT).AND.(HGLNK(IAT).GT.0))
C
      END
C
C====================================================================
C
