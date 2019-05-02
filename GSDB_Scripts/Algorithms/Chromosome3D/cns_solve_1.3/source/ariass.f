C======================================================================
C========== assignment of ambiguous NOEs ==============================
C======================================================================
      SUBROUTINE ANALRS(NOEHGL,ISLCT,JSLCT)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'numbers.inc'
      INTEGER NOEHGL(*),ISLCT(*),JSLCT(*)
C local
      INTEGER I, NISLCT, NJSLCT, ICL1, ICL2
      CHARACTER*4 CLASS, CLASS2
      CHARACTER*4 DISPO, ACTION, ACCMOD
      DOUBLE PRECISION DELPP1, DELPP2
      LOGICAL MATCH, QSORT
      INTEGER HPNDDD, HPNDII, HPNDJJ, HPNDD2, HPNDI2, HPNDJ2
      INTEGER HPNFRW, HPNFR2
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      CALL FILL4(ISLCT,NATOM,1)
      CALL FILL4(JSLCT,NATOM,1)
      DISPO = 'SATI'
      DELPP2 = ZERO
      DELPP1 = ZERO
      ACTION = 'PARA'
      ACCMOD = 'PARA'
      QSORT=.FALSE.
C
      IF (HPNMAT.EQ.0) THEN
      CALL ALLDSP (NOEHGL,
     &       HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &       HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &       HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &       HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &       HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      HPNMAT=ALLHP(IREAL8(MATDIM))
      HPNMA2=ALLHP(INTEG4(MATDIM))
      CALL FILLR8(HEAP(HPNMAT),MATDIM,ZERO)
      CALL FILL4(HEAP(HPNMA2),MATDIM,1)
      END IF
C
      CALL PUSEND('ANALyse_restraints>')
C
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('ANALyse_restraints>')
C
C[
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-aria-analyse')
C
      ELSE IF ((WD(1:4).EQ.'EQUI').OR.(WD(1:4).EQ.'UNRE')) THEN
      CALL EQVGRP(NATOM,NOEHGL)
C
      ELSE IF (WD(1:4).EQ.'LEVE') THEN
      CALL NEXTF('assignment-LEVEl: ',AMBLEV)
C
      ELSE IF (WD(1:4).EQ.'CUTO') THEN
      CALL NEXTF('CUTOff: ',AMBCUT)
C
      ELSE IF (WD(1:4).EQ.'DEL1') THEN
      CALL NEXTF('DEL1 (F1): ',DELPP1)
C
      ELSE IF (WD(1:4).EQ.'DEL2') THEN
      CALL NEXTF('DEL2 (F2): ',DELPP2)
C
      ELSE IF (WD(1:4).EQ.'DISP') THEN
      CALL NEXTA4('DISPosition (SATI/VIOL/ALL): ',DISPO)
C
      ELSE IF (WD(1:4).EQ.'MODE') THEN
      CALL NEXTA4('MODE: ',AMBMOD)
C
      ELSE IF (WD(1:4).EQ.'MINN') THEN
      CALL NEXTI('Minimum number of possibilities: ',AMBMIN)
      ELSE IF (WD(1:4).EQ.'MAXN') THEN
      CALL NEXTI('Maximum number of possibilities: ',AMBMAX)
      ELSE IF ((WD(1:4).EQ.'NUMB').OR.(WD(1:4).EQ.'NPOS')) THEN
      CALL PUSEND('NumberPOSsible>')
C
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('NumberPOSsible>')
C[
      IF (WD(1:4).EQ.'MINI') THEN
      CALL NEXTI('Minimum number of possibilities: ',AMBMIN)
      ELSE IF (WD(1:4).EQ.'MAXI') THEN
      CALL NEXTI('Maximum number of possibilities: ',AMBMAX)
      ELSE IF (WD(1:4).EQ.'NUMB') THEN
      CALL NEXTI('Number of possibilities: ',AMBMAX)
      AMBMIN=AMBMAX
      ELSE
      CALL CHKEND('NumberPOSsible>',DONE)
      END IF
C
      END DO
C
      DONE=.FALSE.
C
      ELSE IF ((WD(1:4).EQ.'RANG')) THEN
      CALL PUSEND('residueRANGe>')
C
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('residueRANGe>')
C[
      IF (WD(1:4).EQ.'MINI') THEN
      CALL NEXTI('MINImum residue range: ',RNGMIN)
      ELSE IF (WD(1:4).EQ.'MAXI') THEN
      CALL NEXTI('MAXImum residue range: ',RNGMAX)
      ELSE IF (WD(1:4).EQ.'NUMB') THEN
      CALL NEXTI('Residue range: ',RNGMAX)
      RNGMIN=RNGMAX
      ELSE
      CALL CHKEND('residueRANGe>',DONE)
      END IF
C
      END DO
C
      DONE=.FALSE.
C
      ELSE IF (WD(1:4).EQ.'FROM') THEN
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      ELSE IF (WD(1:4).EQ.'TO  ') THEN
      CALL SELCTA(JSLCT,NJSLCT,X,Y,Z,.TRUE.)
C
      ELSE IF (WD(1:4).EQ.'RESE') THEN
C
      IF (HPNMA2.NE.0) CALL FREHP(HPNMA2,INTEG4(MATDIM))
      IF (HPNMAT.NE.0) CALL FREHP(HPNMAT,IREAL8(MATDIM))
      HPNMAT=0
      HPNMA2=0
      AMBSTR=0
      AMBLEV=ONE
      AMBMAX=1
      AMBMIN=1
      AMBMOD='ALL '
      CALL FILL4(ISLCT,NATOM,1)
      CALL FILL4(JSLCT,NATOM,1)
C
      ELSE IF (WD(1:4).EQ.'ACCU') THEN
      AMBSTR=AMBSTR+1
      AVEMOD='R-6 '
      ACCMOD='ACCU'
      ACTION='ACCU'
C
      ELSE IF (WD(1:4).EQ.'AVER') THEN
      AMBSTR=AMBSTR+1
      AVEMOD='AVER'
      ACCMOD='AVER'
      ACTION='ACCU'
C
      ELSE IF (WD(1:4).EQ.'MINI') THEN
      AMBSTR=1
      AVEMOD='MINI'
      ACCMOD='MINI'
      ACTION='ACCU'
C
      ELSE IF (WD(1:4).EQ.'SING') THEN
      AMBSTR=1
      AVEMOD='AVER'
      ACCMOD='SING'
      ACTION='ACCU'
C
      ELSE IF (WD(1:4).EQ.'FREQ') THEN
      AMBSTR=1
      AVEMOD='AVER'
      ACCMOD='FREQ'
      ACTION='ACCU'
C
      ELSE IF (WD(1:4).EQ.'FRQU') THEN
      CALL NEXTA4('class-name=',CLASS)
C
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
C[
      IF (MATCH) THEN
      HPNDDD=ALLHP(IREAL8(MAXSEL))
      HPNDII=ALLHP(INTEG4(MAXSEL))
      HPNDJJ=ALLHP(INTEG4(MAXSEL))
      HPNFRW=ALLHP(IREAL8(NATOM))
      CALL FRQUPD(NOEHGL,HEAP(HPNMAT),HEAP(HPNDDD),HEAP(HPNDII),
     &            HEAP(HPNDJJ),HEAP(HPNFRW),
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      CALL FREHP(HPNDJJ,INTEG4(MAXSEL))
      CALL FREHP(HPNDII,INTEG4(MAXSEL))
      CALL FREHP(HPNFRW,IREAL8(NATOM))
      CALL FREHP(HPNDDD,IREAL8(MAXSEL))
      END IF
C
      END DO
C
C
      ELSE IF (WD(1:4).EQ.'VSTA') THEN
      CALL NEXTA4('class-name=',CLASS)
C
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
C[
      IF (MATCH) THEN
      HPNDDD=ALLHP(IREAL8(MAXSEL))
      HPNDII=ALLHP(INTEG4(MAXSEL))
      HPNDJJ=ALLHP(INTEG4(MAXSEL))
      HPNFRW=ALLHP(IREAL8(NATOM))
      CALL VISTAT(NOEHGL,HEAP(HPNMAT),HEAP(HPNDDD),HEAP(HPNDII),
     &            HEAP(HPNDJJ),HEAP(HPNFRW),
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      CALL FREHP(HPNFRW,IREAL8(NATOM))
      CALL FREHP(HPNDJJ,INTEG4(MAXSEL))
      CALL FREHP(HPNDII,INTEG4(MAXSEL))
      CALL FREHP(HPNDDD,IREAL8(MAXSEL))
      END IF
C
      END DO
C
      ELSE IF (WD(1:4).EQ.'RSTA') THEN
      CALL NEXTA4('class-name=',CLASS)
C
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
C[
      IF (MATCH) THEN
      HPNDDD=ALLHP(IREAL8(MAXSEL))
      HPNDII=ALLHP(INTEG4(MAXSEL))
      HPNDJJ=ALLHP(INTEG4(MAXSEL))
      HPNFRW=ALLHP(IREAL8(NATOM))
      CALL RESTAT(NOEHGL,HEAP(HPNMAT),HEAP(HPNDDD),HEAP(HPNDII),
     &            HEAP(HPNDJJ),HEAP(HPNFRW),
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      CALL FREHP(HPNFRW,IREAL8(NATOM))
      CALL FREHP(HPNDJJ,INTEG4(MAXSEL))
      CALL FREHP(HPNDII,INTEG4(MAXSEL))
      CALL FREHP(HPNDDD,IREAL8(MAXSEL))
      END IF
C
      END DO
C
      ELSE IF (WD(1:4).EQ.'CHEC') THEN
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      HPNDDD=ALLHP(IREAL8(MAXSEL))
      HPNDII=ALLHP(INTEG4(MAXSEL))
      HPNDJJ=ALLHP(INTEG4(MAXSEL))
      HPNDD2=ALLHP(IREAL8(MAXSEL))
      HPNDI2=ALLHP(INTEG4(MAXSEL))
      HPNDJ2=ALLHP(INTEG4(MAXSEL))
C
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
C[
      IF (MATCH) THEN
      CALL NOECHK(NOEHGL,HEAP(HPNMAT),
     &            HEAP(HPNDDD),HEAP(HPNDII),HEAP(HPNDJJ),
     &            HEAP(HPNDD2),HEAP(HPNDI2),HEAP(HPNDJ2),
     &            AVEMOD,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      END IF
C
      END DO
C
      CALL FREHP(HPNDJ2,INTEG4(MAXSEL))
      CALL FREHP(HPNDI2,INTEG4(MAXSEL))
      CALL FREHP(HPNDD2,IREAL8(MAXSEL))
      CALL FREHP(HPNDJJ,INTEG4(MAXSEL))
      CALL FREHP(HPNDII,INTEG4(MAXSEL))
      CALL FREHP(HPNDDD,IREAL8(MAXSEL))
C
C
      ELSE IF (WD(1:4).EQ.'PURG') THEN
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTA4('class-name=',CLASS2)
      ICL1=0
      ICL2=0
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) ICL1=I
      CALL EQSTWC(NOECNM(I),4,CLASS2,4,1,1,MATCH)
      IF (MATCH) ICL2=I
      END DO
      HPNDDD=ALLHP(IREAL8(MAXSEL))
      HPNDII=ALLHP(INTEG4(MAXSEL))
      HPNDJJ=ALLHP(INTEG4(MAXSEL))
      HPNDD2=ALLHP(IREAL8(MAXSEL))
      HPNDI2=ALLHP(INTEG4(MAXSEL))
      HPNDJ2=ALLHP(INTEG4(MAXSEL))
C
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
C[
      IF (MATCH) THEN
      CALL NOEPRG
     &          (ICL1,ICL2,NOEHGL,HEAP(HPNMA2),HEAP(HPNMAT),
     &           HEAP(HPNDDD),HEAP(HPNDII),HEAP(HPNDJJ),
     &           HEAP(HPNDD2),HEAP(HPNDI2),HEAP(HPNDJ2),
     &           AVEMOD,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      END IF
C
      END DO
C
      CALL FREHP(HPNDJ2,INTEG4(MAXSEL))
      CALL FREHP(HPNDI2,INTEG4(MAXSEL))
      CALL FREHP(HPNDD2,IREAL8(MAXSEL))
      CALL FREHP(HPNDJJ,INTEG4(MAXSEL))
      CALL FREHP(HPNDII,INTEG4(MAXSEL))
      CALL FREHP(HPNDDD,IREAL8(MAXSEL))

      ELSE IF (WD(1:4).EQ.'TRAN') THEN
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTA4('class-name=',CLASS2)
      ICL1=0
      ICL2=0
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) ICL1=I
      CALL EQSTWC(NOECNM(I),4,CLASS2,4,1,1,MATCH)
      IF (MATCH) ICL2=I
      END DO
      HPNDDD=ALLHP(IREAL8(MAXSEL))
      HPNDII=ALLHP(INTEG4(MAXSEL))
      HPNDJJ=ALLHP(INTEG4(MAXSEL))
      HPNDD2=ALLHP(IREAL8(MAXSEL))
      HPNDI2=ALLHP(INTEG4(MAXSEL))
      HPNDJ2=ALLHP(INTEG4(MAXSEL))
C
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
C[
      IF (MATCH) THEN
      CALL NOETRF
     &          (ICL1,ICL2,NOEHGL,HEAP(HPNMA2),HEAP(HPNMAT),
     &           HEAP(HPNDDD),HEAP(HPNDII),HEAP(HPNDJJ),
     &           HEAP(HPNDD2),HEAP(HPNDI2),HEAP(HPNDJ2),
     &           AVEMOD,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      END IF
C
      END DO
C
      CALL FREHP(HPNDJ2,INTEG4(MAXSEL))
      CALL FREHP(HPNDI2,INTEG4(MAXSEL))
      CALL FREHP(HPNDD2,IREAL8(MAXSEL))
      CALL FREHP(HPNDJJ,INTEG4(MAXSEL))
      CALL FREHP(HPNDII,INTEG4(MAXSEL))
      CALL FREHP(HPNDDD,IREAL8(MAXSEL))
C
      ELSE IF (WD(1:4).EQ.'SILE') THEN
      CALL NEXTA4('class-name=',CLASS)
      AMBFRM='NONE'
      QSORT=.FALSE.
      ACTION='ASSI'
C
      ELSE IF (WD(1:4).EQ.'LIST') THEN
      CALL NEXTA4('class-name=',CLASS)
      AMBFRM='LIST'
      QSORT=.FALSE.
      ACTION='ASSI'
C
      ELSE IF (WD(1:4).EQ.'SORT') THEN
      CALL NEXTA4('class-name=',CLASS)
      AMBFRM='LIST'
      QSORT=.TRUE.
      ACTION='ASSI'
C
      ELSE IF (WD(1:4).EQ.'REST') THEN
      CALL NEXTA4('class-name=',CLASS)
      AMBFRM='REST'
      QSORT=.TRUE.
      ACTION='ASSI'
C
      ELSE IF (WD(1:4).EQ.'OR-R') THEN
      CALL NEXTA4('class-name=',CLASS)
      AMBFRM='OR-R'
      QSORT=.TRUE.
      ACTION='ASSI'
      ELSE
      CALL CHKEND('ANALyse_restraints>',DONE)
      END IF
C
      END DO
C
      DONE=.FALSE.
C
C
      IF (ACTION.EQ.'ASSI') THEN
C
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
C[
      IF (MATCH) THEN
      HPNDDD=ALLHP(IREAL8(MAXSEL))
      HPNDII=ALLHP(INTEG4(MAXSEL))
      HPNDJJ=ALLHP(INTEG4(MAXSEL))
      HPNFR2=ALLHP(IREAL8(MAXSEL))
      CALL AMBLIS(ISLCT,JSLCT,NOEHGL,HEAP(HPNMA2),HEAP(HPNMAT),
     &       HEAP(HPNDDD),HEAP(HPNFR2),HEAP(HPNDII),HEAP(HPNDJJ),
     &       QSORT,AMBFRM,DISPO,I,DELPP1,DELPP2,
     &       HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &       HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &       HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &       HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &       HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      CALL FREHP(HPNFR2,IREAL8(MAXSEL))
      CALL FREHP(HPNDJJ,INTEG4(MAXSEL))
      CALL FREHP(HPNDII,INTEG4(MAXSEL))
      CALL FREHP(HPNDDD,IREAL8(MAXSEL))
      END IF
C
      END DO
C
      ELSEIF (ACTION.EQ.'ACCU') THEN
      WRITE (6,'(A,I8,A,I8)') ' using space of ', MATDIM,
     &  ' at ', HPNMAT
C
      CALL AMBACC(NOEHGL,HEAP(HPNMAT),ACCMOD,
     &       HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &       HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &       HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &       HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &       HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
C
      END IF
C
C
      RETURN
      END
C
C========================================================================
C
      SUBROUTINE AMBACC
     &          (NOEHGL,NOEMAT,MODE,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      DOUBLE PRECISION NOEMAT(*)
      INTEGER NOEHGL(*)
      CHARACTER*4 MODE
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)
C local
      DOUBLE PRECISION XI, YI, ZI, XJ, YJ, ZJ
      DOUBLE PRECISION XIJ, YIJ, ZIJ, SIJ
      DOUBLE PRECISION SIX, SEVEN,ZERO, ONE
      DOUBLE PRECISION TWO, THREE, FOUR
      DOUBLE PRECISION SAVEXP
      INTEGER I, J, N, NSTART, NSTOP, II, JJ, NI, NJ
      LOGICAL QNEXT, QNEXT1, QNEXT2
C for printing 'assignments' of ambiguous NOEs
      INTEGER IAT, JAT, NMAT, NORR
      DOUBLE PRECISION RNL, TOTPRB
C parameter
      PARAMETER (SIX=6.0D0, SEVEN=7.0D0, ONE=1.0D0, TWO=2.0D0)
      PARAMETER (THREE=3.0D0, FOUR=4.0D0, ZERO=0.0D0)
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
      NMAT=0
C
      DO N=NSTART,NSTOP
      SAVEXP=NAVEXP(NOECND(N))
      TOTPRB=ZERO
C
C
      DO NORR=NOEORR(N)+1, NOEORR(N+1)
C
C loop over all pairs of groups of atoms belonging to restraint N
C
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
C
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
C[
      IF (MODE.EQ.'FREQ') THEN
      RNL=ONE
      ELSEIF (NOEAVE(NOECND(N)).EQ.NOECEN) THEN
      NI=0
      XI=ZERO
      YI=ZERO
      ZI=ZERO
      NJ=0
      XJ=ZERO
      YJ=ZERO
      ZJ=ZERO
      I=IAT
C
      QNEXT=.TRUE.
      DO WHILE (QNEXT)
      XI=XI+X(I)
      YI=YI+Y(I)
      ZI=ZI+Z(I)
      NI=NI+1
      CALL NEXTHY(I,NOEHGL,QNEXT)
      END DO
C
      XI=XI/NI
      YI=YI/NI
      ZI=ZI/NI
      J=JAT
C
      QNEXT=.TRUE.
      DO WHILE (QNEXT)
      XJ=XJ+X(J)
      YJ=YJ+Y(J)
      ZJ=ZJ+Z(J)
      NJ=NJ+1
      CALL NEXTHY(J, NOEHGL,QNEXT)
      END DO
C
      XJ=XJ/NJ
      YJ=YJ/NJ
      ZJ=ZJ/NJ
      RNL=(XI-XJ)*(XI-XJ)+(YI-YJ)*(YI-YJ)+(ZI-ZJ)*(ZI-ZJ)
      RNL=RNL**(-SAVEXP/TWO)
      ELSE
      RNL=ZERO
      I=IAT
C
      QNEXT1=.TRUE.
      DO WHILE (QNEXT1)
C
      IF (.NOT.(INITIA(I,X,Y,Z)))
     & CALL WRNDIE(-5,'ATMCHK','Unknown coordinates')
      J=JAT
C
      QNEXT2=.TRUE.
      DO WHILE (QNEXT2)
C
      IF (.NOT.(INITIA(J,X,Y,Z)))
     & CALL WRNDIE(-5,'ATMCHK','Unknown coordinates')
      XIJ=X(I)-X(J)
      YIJ=Y(I)-Y(J)
      ZIJ=Z(I)-Z(J)
      SIJ=XIJ**2+YIJ**2+ZIJ**2
      RNL=RNL+ONE/(SIJ**(SAVEXP/TWO))
      CALL NEXTHY(J, NOEHGL,QNEXT2)
      END DO
C
      CALL NEXTHY(I, NOEHGL,QNEXT1)
      END DO
C
      END IF
C
      NMAT=NMAT+1
      TOTPRB=TOTPRB+RNL
C[
      IF (NMAT.GT.MATDIM) THEN
      WRITE(6,'(A,2I8)') ' AMBACC: incompatible array dimensions: ',
     &      NMAT,MATDIM
      ELSE
C[
      IF ((MODE.EQ.'FREQ')) THEN
      NOEMAT(NMAT)=RNL
      ELSEIF ((MODE.EQ.'SING').OR.(AMBSTR.EQ.1)) THEN
      NOEMAT(NMAT)=RNL
      ELSEIF (MODE.EQ.'MINI') THEN
      NOEMAT(NMAT)=MAX(NOEMAT(NMAT),RNL)
      ELSEIF (MODE.EQ.'ACCU') THEN
      NOEMAT(NMAT)=(NOEMAT(NMAT)*(AMBSTR-1)+RNL)/(AMBSTR)
      ELSEIF (MODE.EQ.'AVER') THEN
      NOEMAT(NMAT)=
     &  (((NOEMAT(NMAT)**(-ONE/SAVEXP))*(AMBSTR-1)
     &     +(RNL**(-ONE/SAVEXP)))/AMBSTR)**(-SAVEXP)
      ELSE
      WRITE(6,'(A,A)') ' AMBACC: unknown mode: ', MODE
      END IF
C
      END IF
C
      END IF
C
      END DO
C
      END IF
C
      END DO
C
      END DO
C
      END DO
      WRITE(6,'(A,I5,A,A)') ' AMBACC: read in ', NMAT,
     &  ' distances in mode ', MODE
C
      RETURN
      END
C====================================================================
      SUBROUTINE AMBLIS
     &     (ISLCT,JSLCT,NOEHGL,NOEFLG,NOEMAT,
     &      ASPROB,FRQWGT,SELI,SELJ,SORT,FORM,DISPO,ICL,
     &      DELPP1,DELPP2,
     &      NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &      NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &      NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'mtf.inc'
      INTEGER ISLCT(*),JSLCT(*),NOEHGL(*)
      INTEGER NOEFLG(*)
      DOUBLE PRECISION NOEMAT(*)
      DOUBLE PRECISION ASPROB(*)
      DOUBLE PRECISION FRQWGT(*)
      INTEGER SELI(*),SELJ(*)
      LOGICAL SORT
      CHARACTER*4 FORM,DISPO
      INTEGER ICL
      DOUBLE PRECISION DELPP1,DELPP2
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)
C local
      DOUBLE PRECISION SAVEXP
      INTEGER I, J, K, N, L, NSTART, NSTOP, II, JJ
      INTEGER LSEL, NCV
C for printing 'assignments' of ambiguous NOEs
      INTEGER IAT, JAT, NMAT, IAMB, TEMPI, TEMPJ, NORR
      INTEGER NRESPR,NAMB
      DOUBLE PRECISION RNL, TOTPRB, ACCLEV
      DOUBLE PRECISION TEMPD, TEMPF, FRQDI1, FRQDI2
      CHARACTER*3 TSTRNG
      CHARACTER*4 PSNAMI,PSNAMJ
      LOGICAL FIRSTI, FIRSTJ, MATCH, QISLCT, QJSLCT, QWRITE, FWRITE
      LOGICAL PSOK
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
      NMAT=0
C
      FWRITE=.TRUE.
      DO N=NSTART,NSTOP
        IF (NOECV(N).EQ.NOEICV) THEN
          NCV=NOECV(N)
          NOECV(N)=0
        END IF
        SAVEXP=NAVEXP(NOECND(N))
C
C loop over all pairs of atoms belonging to restraint N
        L=0
        LSEL=0
        DO NORR=NOEORR(N)+1, NOEORR(N+1)
          DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
            IAT=NOEILS(II)
            IF (FRSTEL(IAT,NOEHGL)) THEN
              DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
                JAT=NOEJLS(JJ)
                IF (FRSTEL(JAT,NOEHGL)) THEN
                  NMAT=NMAT+1
                  IF (NMAT.GT.MATDIM) THEN
                      CALL WRNDIE(-5,'AMBLIS',
     &                  'programming error: MATDIM')
                  END IF
                  L=L+1
                  IF (NOEFLG(NMAT).GT.0) THEN
                    LSEL=LSEL+1
                    IF (LSEL.GT.MAXSEL) THEN
                      CALL WRNDIE(-5,'AMBLIS',
     &                  'programming error: MAXSEL')
                    END IF
                    ASPROB(LSEL)=NOEMAT(NMAT)/NMONO(NOECND(N))
                    FRQWGT(LSEL)=ONE
                    SELI(LSEL)=IAT
                    SELJ(LSEL)=JAT
                  END IF
                END IF
              END DO
            END IF
          END DO
        END DO
C
        IF (LSEL.EQ.0) THEN
          WRITE (6, '(A, 2I8)')
     &  ' %AMBLIS-ERR: problem with restraint ',
     &    N, NOEPID(N)
        ELSE
C
C multiply the distances with frequency weighting
          IF ((DELPP1*DELPP2 .NE. ZERO)
     &      .AND.(NOEPP1(N).GT.(-999.0))
     &      .AND.(NOEPP2(N).GT.(-999.0))) THEN
            DO I=1,LSEL
              IF ((RMSD(SELI(I)).LE.(-999.0))
     &        .OR.(RMSD(SELJ(I)).LE.(-999.0))) THEN
                FRQWGT(I)=ONE
              ELSE
                FRQWGT(I)=EXP(-(RMSD(SELI(I))-NOEPP1(N))**2/DELPP1**2
     &                      -(RMSD(SELJ(I))-NOEPP2(N))**2/DELPP2**2)
              END IF
            END DO
          ELSE
            DO I=1,LSEL
              FRQWGT(I)=ONE
            END DO
          END IF
C
C calculate total probability to normalize
          TOTPRB=ZERO
          DO I=1,LSEL
            TOTPRB=TOTPRB+ASPROB(I)*FRQWGT(I)
          END DO
C
C sort the ASPROB, frequency and atom pointer arrays
          IF (SORT) THEN
            DO I=1,LSEL-1
              DO J=I+1,LSEL
                IF ((ASPROB(J)*FRQWGT(J)).GT.(ASPROB(I)*FRQWGT(I))) THEN
                  TEMPD=ASPROB(I)
                  ASPROB(I)=ASPROB(J)
                  ASPROB(J)=TEMPD
                  TEMPF=FRQWGT(I)
                  FRQWGT(I)=FRQWGT(J)
                  FRQWGT(J)=TEMPF
                  TEMPI=SELI(I)
                  SELI(I)=SELI(J)
                  SELI(J)=TEMPI
                  TEMPJ=SELJ(I)
                  SELJ(I)=SELJ(J)
                  SELJ(J)=TEMPJ
                END IF
                IF (FORM.EQ.'DUMP') THEN
                  WRITE(PUNIT,
     &              '(I5,I5,8(1X,A),F6.2,2F9.3,2E10.3)')
     &            I,J,SEGID(SELI(I)),RESID(SELI(I)),
     &            RES(SELI(I)),TYPE(SELI(I)),
     &            SEGID(SELJ(I)),RESID(SELJ(I)),RES(SELJ(I)),
     &            TYPE(SELJ(I)),RNL,
     &            RMSD(SELI(I)), RMSD(SELJ(I)),
     &            ASPROB(I)*FRQWGT(I)/TOTPRB,FRQWGT(I)
                END IF
              END DO
            END DO
          END IF
C
C calculate weighted ideal frequency
          FRQDI1=ZERO
          FRQDI2=ZERO
          IF ((RMSD(SELI(1)).GT.(-999.0))
     &        .AND.(RMSD(SELJ(1)).GT.(-999.0))
     &        .AND.(NOEPP1(N).GT.(-999.0))
     &        .AND.(NOEPP2(N).GT.(-999.0))) THEN
            DO I=1,LSEL
              FRQDI1=FRQDI1+ASPROB(I)*FRQWGT(I)/TOTPRB*RMSD(SELI(I))
              FRQDI2=FRQDI2+ASPROB(I)*FRQWGT(I)/TOTPRB*RMSD(SELJ(I))
            END DO
            NOERAV(N)=EXP(-(FRQDI1-NOEPP1(N))**2/DELPP1**2
     &                   -(FRQDI2-NOEPP2(N))**2/DELPP2**2)
          ELSE
            FRQDI1=-ONE
            FRQDI2=-ONE
            NOERAV(N)=ZERO
          END IF
C
C establish how many restraints one needs to exceed acclev or ambcut
C ambcut is a simple distance cutoff
C amblev is the r-6 sum over the accumulated average or minimum distances
          IAMB=0
          ACCLEV=ZERO
          RNL=ZERO
          DO WHILE ((ACCLEV.LT.AMBLEV)
     &                .AND.(IAMB.LT.LSEL)
     &                .AND.(RNL.LE.AMBCUT))
            IAMB=IAMB+1
            ACCLEV=ACCLEV+ASPROB(IAMB)*FRQWGT(IAMB)/TOTPRB
            RNL=MAX(ASPROB(IAMB),RSMALL)**(-ONE/SAVEXP)
          END DO
          IF (RNL.GT.AMBCUT) IAMB=IAMB-1
          IF ((IAMB.LE.0).AND.(NOECV(N).GT.0)) THEN
          WRITE (6, '(A, 2I8)')
     &      ' %AMBLIS-ERR: no distance < ambcut, one kept, at ',
     &      N, NOEPID(N)
            IAMB=1
          END IF
C
C count how many different residue pairs are involved
          IF ((AMBMOD.EQ.'ALL ').OR.(IAMB.LE.1)) THEN
            NAMB=IAMB
          ELSE
            NRESPR=1
            DO K=1,IAMB-1
              MATCH=.FALSE.
              DO LSEL=K+1,IAMB
                IF ((   (SEGID(SELI(K)).EQ.SEGID(SELI(LSEL)))
     &             .AND.(RESID(SELI(K)).EQ.RESID(SELI(LSEL)))
     &             .AND.(SEGID(SELJ(K)).EQ.SEGID(SELJ(LSEL)))
     &             .AND.(RESID(SELJ(K)).EQ.RESID(SELJ(LSEL))))
     &           .OR.(  (SEGID(SELI(K)).EQ.SEGID(SELJ(LSEL)))
     &             .AND.(RESID(SELI(K)).EQ.RESID(SELJ(LSEL)))
     &             .AND.(SEGID(SELJ(K)).EQ.SEGID(SELI(LSEL)))
     &             .AND.(RESID(SELJ(K)).EQ.RESID(SELI(LSEL))))) THEN
                  MATCH=.TRUE.
                END IF
              END DO
              IF (.NOT.MATCH) NRESPR=NRESPR+1
            END DO
            NAMB=NRESPR
          END IF
C
C check is the restraint is selected
          QISLCT=.FALSE.
          QJSLCT=.FALSE.
          DO I=1, IAMB
            IF ((ISLCT(SELI(I)).GT.0)) QISLCT=.TRUE.
            IF ((JSLCT(SELJ(I)).GT.0)) QJSLCT=.TRUE.
          END DO
          IF (.NOT.(QISLCT.AND.QJSLCT)) THEN
            QISLCT=.FALSE.
            QJSLCT=.FALSE.
            DO I=1, IAMB
              IF ((ISLCT(SELJ(I)).GT.0)) QISLCT=.TRUE.
              IF ((JSLCT(SELI(I)).GT.0)) QJSLCT=.TRUE.
            END DO
          END IF
C
          IF (NOECV(N).EQ.0) THEN
            TSTRNG=' - '
          ELSE
            TSTRNG=' + '
          END IF
C
C write a distance restraint with ambiguous atom selections
C here the inter mode does not make much sense (yet?)
          IF ((FORM.EQ.'REST')) THEN
            QWRITE =
     &      (((DISPO.EQ.'ALL ').OR.
     &        ((NOECV(N).NE.0).AND.(DISPO.EQ.'SATI')).OR.
     &        ((NOECV(N).EQ.0).AND.(DISPO.EQ.'VIOL')))
     &         .AND.QISLCT.AND.QJSLCT
     &         .AND.(NAMB.LE.AMBMAX).AND.(NAMB.GE.AMBMIN)
     &         .AND.(ICL.EQ.NOECND(N)))

c     3    .AND.(NRANGE.LE.RNGMAX).AND.(NRANGE.GE.RNGMIN))
C
C sort the seli and selj arrays until iamb
            IF (QWRITE) THEN
C
              DO I=1,IAMB-1
                DO J=I+1,IAMB
                  IF (SELI(J).GT.SELI(I)) THEN
                    TEMPI=SELI(I)
                    SELI(I)=SELI(J)
                    SELI(J)=TEMPI
                  END IF
                  IF (SELJ(J).GT.SELJ(I)) THEN
                    TEMPJ=SELJ(I)
                    SELJ(I)=SELJ(J)
                    SELJ(J)=TEMPJ
                  END IF
                END DO
              END DO
C
              TEMPI=-9999
              TEMPJ=-9999
              FIRSTI=.TRUE.
              FIRSTJ=.TRUE.
C write I atoms
              WRITE(PUNIT,'(1X,A)') 'ASSI ('
              DO I=1,IAMB
                IF (SELI(I).NE.TEMPI) THEN
                  TEMPI=SELI(I)
                  IF (.NOT.FIRSTI) THEN
                    WRITE(PUNIT,'(1X,A)') ' OR '
                  END IF
                  CALL RXWRA2(PUNIT,(SELI(I)),NOEHGL,'PSEU')
                  FIRSTI=.FALSE.
                END IF
              END DO
C
C write J atoms
              WRITE(PUNIT,'(1X,A)') ')('
              DO I=1,IAMB
                IF (SELJ(I).NE.TEMPJ) THEN
                  TEMPJ=SELJ(I)
                  IF (.NOT.FIRSTJ) THEN
                    WRITE(PUNIT,'(1X,A)') ' OR '
                  END IF
                  CALL RXWRA2(PUNIT,(SELJ(I)),NOEHGL,'PSEU')
                  FIRSTJ=.FALSE.
                END IF
              END DO
              WRITE(PUNIT,'(1X,A)') ')'
C
C write restraint
              WRITE(PUNIT,'(1X,3F10.3,A,I5,2(A,E12.5),2(A,F10.3))')
     &        NOEDIS(N),NOELOW(N),NOEHIG(N),
     &        ' peak ',NOEPID(N),' weight ',NOEWGH(N),
     &        ' volume ',NOEVOL(N),
     &        ' ppm1 ', NOEPP1(N), ' ppm2 ', NOEPP2(N)
            END IF
C
C write an ambiguous distance restraint in "OR" format
          ELSE IF ((FORM.EQ.'OR-R')) THEN
            QWRITE =
     &      (((DISPO.EQ.'ALL ').OR.
     &        ((NOECV(N).NE.0).AND.(DISPO.EQ.'SATI')).OR.
     &        ((NOECV(N).EQ.0).AND.(DISPO.EQ.'VIOL')))
     &         .AND.QISLCT.AND.QJSLCT
     &         .AND.(NAMB.LE.AMBMAX).AND.(NAMB.GE.AMBMIN)
     &         .AND.(ICL.EQ.NOECND(N)))
C
            IF (QWRITE) THEN
              FIRSTI=.TRUE.
              WRITE(PUNIT,'(1X,A, I5, A)') 'ASSI {', NOEPID(N), '}'
              DO I=1,IAMB
                IF (FIRSTI) THEN
                  CALL RXWRA2(PUNIT,(SELI(I)),NOEHGL,'PSEU')
                  CALL RXWRA2(PUNIT,(SELJ(I)),NOEHGL,'PSEU')
                  WRITE(PUNIT,'(1X,3F10.3,A,I5,2(A,E12.5),2(A,F10.3))')
     &            NOEDIS(N),NOELOW(N),NOEHIG(N),
     &            ' peak ',NOEPID(N),' weight ',NOEWGH(N),
     &            ' volume ',NOEVOL(N),
     &            ' ppm1 ', NOEPP1(N), ' ppm2 ', NOEPP2(N)
                    FIRSTI=.FALSE.
                ELSE
                  WRITE(PUNIT,'(1X,A, I5, A)') 'OR {', NOEPID(N), '}'
                  CALL RXWRA2(PUNIT,(SELI(I)),NOEHGL,'PSEU')
                  CALL RXWRA2(PUNIT,(SELJ(I)),NOEHGL,'PSEU')
                END IF
              END DO
            END IF
C
          ELSEIF (FORM.EQ.'LIST') THEN
C
            QWRITE = (QISLCT.AND.QJSLCT
     &     .AND.(NAMB.LE.AMBMAX).AND.(NAMB.GE.AMBMIN)
     &     .AND.(ICL.EQ.NOECND(N)))
C
C TODO write first time
            IF (FWRITE) THEN
              WRITE(PUNIT,'(A,A,A,A,A,A,A,A,A,A,A,A,A,A)')
     & ' X ', 'peakno', '     assignment 1   ', '     assignment 2   ',
     & 'd_ave ', '   volume ', 'lower ', 'upper ', '   Hppm1 ',
     & '   Hppm2 ', '   Pppm1 ', '   Pppm2 ', '   Contr', ' Nps'
              FWRITE=.FALSE.
            END IF
            IF (QWRITE) THEN
              DO I=1,IAMB
                CALL PSUNAM((SELI(I)),NOEHGL,PSNAMI,PSOK)
                IF (.NOT.PSOK) PSNAMI=TYPE(SELI(I))
                CALL PSUNAM((SELJ(I)),NOEHGL,PSNAMJ,PSOK)
                IF (.NOT.PSOK) PSNAMJ=TYPE(SELJ(I))
                RNL=MAX(ASPROB(I),RSMALL)**(-ONE/SAVEXP)
                WRITE(PUNIT,
     &          '(A,I5,8(1X,A),F6.2,E10.3,2F6.2,4F9.3,2E10.3,I4)')
     &          TSTRNG,NOEPID(N),SEGID(SELI(I)),RESID(SELI(I)),
     &          RES(SELI(I)),PSNAMI,
     &          SEGID(SELJ(I)),RESID(SELJ(I)),RES(SELJ(I)),
     &          PSNAMJ,RNL,NOEVOL(N),
     &          NOEDIS(N)-NOELOW(N),NOEDIS(N)+NOEHIG(N),
     &          RMSD(SELI(I)), RMSD(SELJ(I)),NOEPP1(N),NOEPP2(N),
     &          ASPROB(I)*FRQWGT(I)/TOTPRB,FRQWGT(I), NAMB
              END DO
            END IF
          ELSE
          END IF
          IF (NOECV(N).EQ.0) NOECV(N)=NCV
        END IF
      END DO
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE FRQUPD
     &          (NOEHGL,NOEMAT,ASPROB,SELI,SELJ,FRQWHT,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'mtf.inc'
      DOUBLE PRECISION NOEMAT(*), ASPROB(*)
      INTEGER NOEHGL(*)
      INTEGER SELI(*),SELJ(*)
      DOUBLE PRECISION FRQWHT(*)
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)
      INTEGER I,II,IAT,JJ,JAT,L,N,NORR,NSTART,NSTOP,NMAT
      DOUBLE PRECISION RN,RNL,TOTPRB,SAVEXP
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
      NMAT=0
C
C now we update the RMSD array
C at a later stage, this will be a PPM array
C
      DO I=1,NATOM
      RMSD(I)=ZERO
      FRQWHT(I)=ZERO
      END DO
C
      DO N=NSTART,NSTOP
      SAVEXP=NAVEXP(NOECND(N))
C
C loop over all pairs of atoms belonging to restraint N
      L=0
C
      DO NORR=NOEORR(N)+1, NOEORR(N+1)
C
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
C
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      NMAT=NMAT+1
      L=L+1
      ASPROB(L)=NOEMAT(NMAT)/NMONO(NOECND(N))
      SELI(L)=IAT
      SELJ(L)=JAT
      END IF
C
      END DO
C
      END IF
C
      END DO
C
      END DO
C[
      IF (L.EQ.0) THEN
      WRITE (6, '(A, 2I8)') ' %AMBLIS-ERR: problem with restraint ',
     &  N, NOEPID(N)
C
      ELSE
C calculate the r-6 summed distance^-6
      RN=ZERO
C
      DO I=1,L
      RN=RN+ASPROB(I)
      END DO
C
      TOTPRB=MAX(RN**(-ONE/SAVEXP),RSMALL)
C
C calculate new set of ideal frequencies
C[
      IF ((NOECV(N).GT.0)
     &    .AND.(NOEPP1(N).GT.(-9999.0))
     &    .AND.(NOEPP2(N).GT.(-9999.0))) THEN
C
      DO I=1,L
      RNL=MAX(ASPROB(I)**(-ONE/SAVEXP),RSMALL)
      RMSD(SELI(I))=RMSD(SELI(I))+(TOTPRB/RNL)**SAVEXP*NOEPP1(N)
      FRQWHT(SELI(I))=FRQWHT(SELI(I))+(TOTPRB/RNL)**SAVEXP
      RMSD(SELJ(I))=RMSD(SELJ(I))+(TOTPRB/RNL)**SAVEXP*NOEPP2(N)
      FRQWHT(SELJ(I))=FRQWHT(SELJ(I))+(TOTPRB/RNL)**SAVEXP
      END DO
C
      END IF
C
      END IF
C
      END DO
C
      DO I=1,NATOM
C[
      IF (FRQWHT(I).GT.RSMALL) THEN
      RMSD(I)=RMSD(I)/FRQWHT(I)
      ELSE
      RMSD(I)=-9999.0D0
      END IF
C
      END DO
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE RESTAT
     &          (NOEHGL,NOEMAT,ASPROB,SELI,SELJ,FRQWHT,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'mtf.inc'
      DOUBLE PRECISION NOEMAT(*), ASPROB(*)
      INTEGER NOEHGL(*)
      INTEGER SELI(*),SELJ(*)
      DOUBLE PRECISION FRQWHT(*)
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)
      INTEGER I,II,IAT,JJ,JAT,L,N,NORR,NSTART,NSTOP,NMAT
      DOUBLE PRECISION RN,RNL,TOTPRB,SAVEXP
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
      NMAT=0
C
C now we update the RMSD array
C at a later stage, this will be a PPM array
C
      DO I=1,NATOM
      RMSD(I)=ZERO
      FRQWHT(I)=ZERO
      END DO
C
      DO N=NSTART,NSTOP
      SAVEXP=NAVEXP(NOECND(N))
C
C loop over all pairs of atoms belonging to restraint N
      L=0
C
      DO NORR=NOEORR(N)+1, NOEORR(N+1)
C
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
C
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      NMAT=NMAT+1
      L=L+1
      ASPROB(L)=NOEMAT(NMAT)/NMONO(NOECND(N))
      SELI(L)=IAT
      SELJ(L)=JAT
      END IF
C
      END DO
C
      END IF
C
      END DO
C
      END DO
C
C[
      IF (L.EQ.0) THEN
      WRITE (6, '(A, 2I8)') ' %AMBLIS-ERR: problem with restraint ',
     &  N, NOEPID(N)
C
      ELSE
C calculate the r-6 summed distance^-6
      RN=ZERO
C
      DO I=1,L
      RN=RN+ASPROB(I)
      END DO
C
      TOTPRB=MAX(RN**(-ONE/SAVEXP),RSMALL)
C
C calculate new set of ideal frequencies
C[
      IF ((NOECV(N).NE.0)) THEN
C
      DO I=1,L
      RNL=MAX(ASPROB(I)**(-ONE/SAVEXP),RSMALL)
      RMSD(SELI(I))=RMSD(SELI(I))+(TOTPRB/RNL)**SAVEXP
      FRQWHT(SELI(I))=FRQWHT(SELI(I))+(TOTPRB/RNL)**SAVEXP
      RMSD(SELJ(I))=RMSD(SELJ(I))+(TOTPRB/RNL)**SAVEXP
      FRQWHT(SELJ(I))=FRQWHT(SELJ(I))+(TOTPRB/RNL)**SAVEXP
      END DO
C
      END IF
C
      END IF
C
      END DO
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE VISTAT
     &          (NOEHGL,NOEMAT,ASPROB,SELI,SELJ,FRQWHT,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'mtf.inc'
      DOUBLE PRECISION NOEMAT(*), ASPROB(*)
      INTEGER NOEHGL(*)
      INTEGER SELI(*),SELJ(*)
      DOUBLE PRECISION FRQWHT(*)
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)
      INTEGER I,II,IAT,JJ,JAT,L,N,NORR,NSTART,NSTOP,NMAT
      DOUBLE PRECISION RN,RNL,TOTPRB,SAVEXP
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
      NMAT=0
C
C now we update the RMSD array
C at a later stage, this will be a PPM array
C
      DO I=1,NATOM
        RMSD(I)=ZERO
        FRQWHT(I)=ZERO
      END DO
C
      DO N=NSTART,NSTOP
      SAVEXP=NAVEXP(NOECND(N))
C
C loop over all pairs of atoms belonging to restraint N
      L=0
C
      DO NORR=NOEORR(N)+1, NOEORR(N+1)
C
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
C
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      NMAT=NMAT+1
      L=L+1
      ASPROB(L)=NOEMAT(NMAT)/NMONO(NOECND(N))
      SELI(L)=IAT
      SELJ(L)=JAT
      END IF
C
      END DO
C
      END IF
C
      END DO
C
      END DO
C
C[
      IF (L.EQ.0) THEN
      WRITE (6, '(A, 2I8)') ' %AMBLIS-ERR: problem with restraint ',
     &  N, NOEPID(N)
C
      ELSE
C calculate the r-6 summed distance^-6
      RN=ZERO
C
      DO I=1,L
      RN=RN+ASPROB(I)
      END DO
C
      TOTPRB=MAX(RN**(-ONE/SAVEXP),RSMALL)
C
C calculate new set of ideal frequencies
C[
      IF ((NOECV(N).EQ.0)) THEN
C
      DO I=1,L
        RNL=MAX(ASPROB(I)**(-ONE/SAVEXP),RSMALL)
        RMSD(SELI(I))=RMSD(SELI(I))+(TOTPRB/RNL)**SAVEXP
        FRQWHT(SELI(I))=FRQWHT(SELI(I))+(TOTPRB/RNL)**SAVEXP
        RMSD(SELJ(I))=RMSD(SELJ(I))+(TOTPRB/RNL)**SAVEXP
        FRQWHT(SELJ(I))=FRQWHT(SELJ(I))+(TOTPRB/RNL)**SAVEXP
      END DO
C
      END IF
C
      END IF
C
      END DO
C
      RETURN
      END
C
C====================================================================
C
      SUBROUTINE ALLDSP(NOEHGL,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
C
      INTEGER NOEHGL(*)
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)
C local
      DOUBLE PRECISION SIX, SEVEN,ZERO, ONE
      DOUBLE PRECISION TWO, THREE, FOUR
      INTEGER N, NSTART, NSTOP, II, JJ
C for printing 'assignments' of ambiguous NOEs
      INTEGER IAT, JAT, NMAT,NSEL, NORR
C parameter
      PARAMETER (SIX=6.0D0, SEVEN=7.0D0, ONE=1.0D0, TWO=2.0D0)
      PARAMETER (THREE=3.0D0, FOUR=4.0D0, ZERO=0.0D0)
C
C function
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
      NMAT=0
      MAXSEL=0
C
      DO N=NSTART,NSTOP
C
C loop over all pairs of atoms belonging to restraint N
        NSEL=0
        DO NORR=NOEORR(N)+1, NOEORR(N+1)
          DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
            IAT=NOEILS(II)
            IF (FRSTEL(IAT,NOEHGL)) THEN
              DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
                JAT=NOEJLS(JJ)
                IF (FRSTEL(JAT,NOEHGL)) THEN
                  NMAT=NMAT+1
                  NSEL=NSEL+1
                  MAXSEL=MAX(MAXSEL,NSEL)
                END IF
              END DO
            END IF
          END DO
        END DO
      END DO
C
      WRITE(6,'(A,I5,A)') ' ALLDSP: making space for ',
     & NMAT,' distances'
      MATDIM=NMAT
      WRITE(6,'(A,I5,A)') ' ALLDSP: maximum number of selections ',
     & MAXSEL
C
      RETURN
      END
C
C======================================================================
C========== checking of NOE list ======================================
C======================================================================
C
      SUBROUTINE NOECHK
     &          (NOEHGL,NOEMAT,
     &           ASPROB,SELI,SELJ,SELDI2,SELI2,SELJ2,
     &           MODE,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
c new noechec
c ideas:
c comparison of involved atom pairs
c use basically amblis code
c first accumulate ensemble into noemat array
c then compare all pairs of restraints (this may take a while)
c sort the restraints wrt distance
c cutoff at specified cutoff
c then compare the two sets of atom pairs for equality
C
      IMPLICIT NONE
C I/O
      INCLUDE 'noe.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'numbers.inc'
      DOUBLE PRECISION NOEMAT(*)
      DOUBLE PRECISION ASPROB(*), SELDI2(*)
      INTEGER NOEHGL(*)
      INTEGER SELI(*),SELJ(*),SELI2(*),SELJ2(*)
      CHARACTER*4 MODE
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)

C local
      DOUBLE PRECISION SAVEXP, RN
      INTEGER I, J, N, L, NSTART, NSTOP, II, JJ
      INTEGER IAT, JAT, NMAT, IAMB, TEMPI, TEMPJ, NORR
      INTEGER NMAT2, N2, I1, I2, IAMB2, NEQUAL
      DOUBLE PRECISION RNL, TOTPRB, ACCLEV, TEMPD
      LOGICAL QDIFF
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
      NMAT=0
C
      DO N=NSTART,NSTOP-1
        SAVEXP=NAVEXP(NOECND(N))
C
C loop over all pairs of atoms belonging to restraint N
        L=0
        DO NORR=NOEORR(N)+1, NOEORR(N+1)
          DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
            IAT=NOEILS(II)
            IF (FRSTEL(IAT,NOEHGL)) THEN
              DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
                JAT=NOEJLS(JJ)
                IF (FRSTEL(JAT,NOEHGL)) THEN
                  NMAT=NMAT+1
                  L=L+1
                  IF (NMAT.GT.MATDIM) THEN
                    CALL WRNDIE(-5,'NOECHK',
     &                  'programming error: MATDIM')
                  END IF
                  IF (L.GT.MAXSEL) THEN
                    CALL WRNDIE(-5,'NOECHK',
     &                  'programming error: MAXSEL')
                  END IF
                  ASPROB(L)=NOEMAT(NMAT)/NMONO(NOECND(N))
                  SELI(L)=IAT
                  SELJ(L)=JAT
                END IF
              END DO
            END IF
          END DO
        END DO
C
        DO I=1,L-1
          DO J=I+1,L
            IF (ASPROB(J).GT.ASPROB(I)) THEN
              TEMPD=ASPROB(I)
              TEMPI=SELI(I)
              TEMPJ=SELJ(I)
              ASPROB(I)=ASPROB(J)
              SELI(I)=SELI(J)
              SELJ(I)=SELJ(J)
              ASPROB(J)=TEMPD
              SELI(J)=TEMPI
              SELJ(J)=TEMPJ
            END IF
          END DO
        END DO
C
C calculate the r-6 summed distance^-6
        RN=ZERO
        DO I=1,L
          RN=RN+ASPROB(I)
        END DO
        TOTPRB=MAX(RN**(-ONE/SAVEXP),RSMALL)
C
C establish how many one needs to exceed acclev or ambcut
C ambcut is a simple distance cutoff
C amblev is the r-6 sum over the accumulated average or minimum distances
        IAMB=0
        ACCLEV=ZERO
        RNL=ZERO
        DO WHILE ((ACCLEV.LT.AMBLEV)
     &          .AND.(IAMB.LT.L)
     &          .AND.(RNL.LE.AMBCUT))
          IAMB=IAMB+1
          RNL=MAX(ASPROB(IAMB)**(-ONE/SAVEXP),RSMALL)
          ACCLEV=ACCLEV+(TOTPRB/RNL)**SAVEXP
        END DO
        IF (RNL.GT.AMBCUT) IAMB=IAMB-1
        IF (IAMB.EQ.0) THEN
          WRITE (6, '(A, 2I8)') ' no distance < ambcut, one kept, at ',
     &    N, NOEPID(N)
          IAMB=1
        END IF
C
C now we loop over the whole matrix again
C only if the restraint is active
C
        NMAT2=NMAT
          DO N2=N+1,NSTOP
            IF ((NOECV(N2).GT.0).AND.(NOECV(N).GT.0)) THEN
              SAVEXP=NAVEXP(NOECND(N2))
C
C loop over all pairs of atoms belonging to restraint N
              L=0
              DO NORR=NOEORR(N2)+1, NOEORR(N2+1)
                DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
                  IAT=NOEILS(II)
                  IF (FRSTEL(IAT,NOEHGL)) THEN
                    DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
                      JAT=NOEJLS(JJ)
                      IF (FRSTEL(JAT,NOEHGL)) THEN
                        NMAT2=NMAT2+1
                        L=L+1
                        IF (NMAT2.GT.MATDIM) THEN
                          CALL WRNDIE(-5,'NOECHK',
     &                    'programming error: MATDIM')
                        END IF
                        IF (L.GT.MAXSEL) THEN
                          CALL WRNDIE(-5,'NOECHK',
     &                    'programming error: MAXSEL')
                        END IF
                        SELDI2(L)=NOEMAT(NMAT2)/NMONO(NOECND(N2))
                        SELI2(L)=IAT
                        SELJ2(L)=JAT
                      END IF
                    END DO
                  END IF
                END DO
              END DO
C
              DO I=1,L-1
                DO J=I+1,L
                  IF (SELDI2(J).GT.SELDI2(I)) THEN
                    TEMPD=SELDI2(I)
                    TEMPI=SELI2(I)
                    TEMPJ=SELJ2(I)
                    SELDI2(I)=SELDI2(J)
                    SELI2(I)=SELI2(J)
                    SELJ2(I)=SELJ2(J)
                    SELDI2(J)=TEMPD
                    SELI2(J)=TEMPI
                    SELJ2(J)=TEMPJ
                  END IF
                END DO
              END DO
C
C calculate the r-6 summed distance^-6
              RN=ZERO
              DO I=1,L
                RN=RN+SELDI2(I)
              END DO
              TOTPRB=MAX(RN**(-ONE/SAVEXP),RSMALL)
C
C establish how many one needs to exceed acclev or ambcut
C ambcut is a simple distance cutoff
C amblev is the r-6 sum over the accumulated average or minimum distances
              IAMB2=0
              ACCLEV=ZERO
              RNL=ZERO
              DO WHILE ((ACCLEV.LT.AMBLEV)
     &          .AND.(IAMB2.LT.L)
     &          .AND.(RNL.LE.AMBCUT))
                IAMB2=IAMB2+1
                RNL=MAX(SELDI2(IAMB2)**(-ONE/SAVEXP),RSMALL)
                ACCLEV=ACCLEV+(TOTPRB/RNL)**SAVEXP
              END DO
              IF (RNL.GT.AMBCUT) IAMB2=IAMB2-1
              IF (IAMB2.EQ.0) THEN
                IAMB2=1
              END IF
C
C now compare the two sets of SEL arrays
C if they differ in at least one position, we accept both
              QDIFF=.FALSE.
              IF (IAMB.NE.IAMB2) THEN
                QDIFF=.TRUE.
              ELSE
              NEQUAL=0
              DO I1=1,IAMB
                DO I2=1,IAMB
                  IF (((SELI(I1).EQ.SELI2(I2))
     &                .AND.(SELJ(I1).EQ.SELJ2(I2)))
     &                .OR.((SELI(I1).EQ.SELJ2(I2))
     &                .AND.(SELJ(I1).EQ.SELI2(I2))))
     &            THEN
                    NEQUAL=NEQUAL+1
                  END IF
                END DO
              END DO
              IF (NEQUAL.LT.IAMB) QDIFF=.TRUE.
            END IF
C
C if qdiff, we accept both restraints
C otherwise we have to decide which one to keep
C ideas:
C we keep the more restraining one of the two
            IF (.NOT.QDIFF) THEN
              IF ((NOEHIG(N)+NOELOW(N))
     &             .GT.(NOEHIG(N2)+NOELOW(N2))) THEN
                NOECV(N)=0
                WRITE (6,'(A,2I8, 2X, 2I8)')
     &  ' Duplicate restraint removed (1): ',
     &  N, NOEPID(N), N2, NOEPID(N2)
              ELSE
                NOECV(N2)=0
                WRITE (6,'(A,2I8, 2X, 2I8)')
     &  ' Duplicate restraint removed (2): ',
     &  N, NOEPID(N), N2, NOEPID(N2)
              END IF
            END IF
          END IF
        END DO
      END DO
C
      RETURN
      END
C
C==========================================================================
C
      SUBROUTINE NOEPRG
     &          (ICL1,ICL2,NOEHGL,NOEFLG,NOEMAT,
     &           ASPROB,SELI,SELJ,SELDI2,SELI2,SELJ2,
     &           MODE,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
C removes any occurances of an unambiguous assignment in class icl1
C from the list of possibilities in class icl2
C
      IMPLICIT NONE
C I/O
      INCLUDE 'noe.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'numbers.inc'
      INTEGER ICL1,ICL2
      INTEGER NOEFLG(*)
      DOUBLE PRECISION NOEMAT(*)
      DOUBLE PRECISION ASPROB(*), SELDI2(*)
      INTEGER NOEHGL(*)
      INTEGER SELI(*),SELJ(*),SELI2(*),SELJ2(*)
      CHARACTER*4 MODE
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)

C local
      DOUBLE PRECISION SAVEXP, RN
      INTEGER I, J, N, L, NSTART, NSTOP, II, JJ
      INTEGER IAT, JAT, NMAT, IAMB, TEMPI, TEMPJ, NORR
      INTEGER NMAT2, N2, I2
      DOUBLE PRECISION RNL, TOTPRB, ACCLEV, TEMPD
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
      NMAT=0
C
      DO N=NSTART,NSTOP-1
        SAVEXP=NAVEXP(NOECND(N))
C
C loop over all pairs of atoms belonging to restraint N
        L=0
        DO NORR=NOEORR(N)+1, NOEORR(N+1)
          DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
            IAT=NOEILS(II)
            IF (FRSTEL(IAT,NOEHGL)) THEN
              DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
                JAT=NOEJLS(JJ)
                  IF (FRSTEL(JAT,NOEHGL)) THEN
        NMAT=NMAT+1
        L=L+1
        ASPROB(L)=NOEMAT(NMAT)/NMONO(NOECND(N))
        SELI(L)=IAT
        SELJ(L)=JAT
                  END IF
              END DO
            END IF
          END DO
        END DO
C
C sort the sel arrayS
        DO I=1,L-1
          DO J=I+1,L
            IF (ASPROB(J).GT.ASPROB(I)) THEN
              TEMPD=ASPROB(I)
              TEMPI=SELI(I)
              TEMPJ=SELJ(I)
              ASPROB(I)=ASPROB(J)
              SELI(I)=SELI(J)
              SELJ(I)=SELJ(J)
              ASPROB(J)=TEMPD
              SELI(J)=TEMPI
              SELJ(J)=TEMPJ
            END IF
          END DO
        END DO
C
C calculate the r-6 summed distance^-6
        RN=ZERO
        DO I=1,L
          RN=RN+ASPROB(I)
        END DO
        TOTPRB=MAX(RN**(-ONE/SAVEXP),RSMALL)
C
C establish how many one needs to exceed acclev or ambcut
C ambcut is a simple distance cutoff
C amblev is the r-6 sum over the accumulated average or minimum distances
        IAMB=0
        ACCLEV=ZERO
        RNL=ZERO
        DO WHILE ((ACCLEV.LT.AMBLEV)
     &          .AND.(IAMB.LT.L)
     &          .AND.(RNL.LE.AMBCUT))
          IAMB=IAMB+1
          RNL=MAX(ASPROB(IAMB)**(-ONE/SAVEXP),RSMALL)
          ACCLEV=ACCLEV+(TOTPRB/RNL)**SAVEXP
        END DO
        IF (RNL.GT.AMBCUT) IAMB=IAMB-1
        IF (IAMB.EQ.0) THEN
          WRITE (6, '(A, 2I8)') ' no distance < ambcut, one kept, at ',
     &      N, NOEPID(N)
          IAMB=1
        END IF
C
C now we loop over the whole matrix again
C only if the restraint is active
C
        NMAT2=NMAT
        DO N2=N+1,NSTOP
          SAVEXP=NAVEXP(NOECND(N2))
C
C loop over all pairs of atoms belonging to restraint N2
          L=0
          DO NORR=NOEORR(N2)+1, NOEORR(N2+1)
            DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
              IAT=NOEILS(II)
              IF (FRSTEL(IAT,NOEHGL)) THEN
                DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
                  JAT=NOEJLS(JJ)
                  IF (FRSTEL(JAT,NOEHGL)) THEN
        NMAT2=NMAT2+1
        L=L+1
        SELDI2(L)=NOEMAT(NMAT2)/NMONO(NOECND(N2))
        SELI2(L)=IAT
        SELJ2(L)=JAT
                  END IF
                END DO
              END IF
            END DO
          END DO
C
C we remove a possibility under the following conditions
C  1- the fist restraint is from class icl1, the second from icl2
C  2- both restraints are active (noecnd > 0)
C  3- restraint 1 is unambiguous (iamb = 1)
C  4- the unambiguous restraint 1 is a possibility in a restraint from icl2
C  5- restraint 2 is ambiguous (L>1)
          IF ((NOECND(N).EQ.ICL1).AND.(NOECND(N2).EQ.ICL2)
     &     .AND.(NOECV(N).GT.0).AND.(NOECV(N2).GT.0)
     &     .AND.(IAMB.EQ.1).AND.(L.GT.1)) THEN
            DO I2=1,L
              IF (((SELI(1).EQ.SELI2(I2))
     &            .AND.(SELJ(1).EQ.SELJ2(I2)))
     &         .OR.((SELI(1).EQ.SELJ2(I2))
     &            .AND.(SELJ(1).EQ.SELI2(I2))))
     &        THEN
                NOEFLG(NMAT2-L+I2)=0
                WRITE (6, '(A, 2I8,A, 2I8)')
     &          ' NOEPRG: ', N, NOEPID(N), ' FROM ', N2, NOEPID(N2)
              END IF
            END DO
          END IF
        END DO
      END DO
C
      RETURN
      END
C
C==========================================================================
C
      SUBROUTINE NOETRF
     &          (ICL1,ICL2,NOEHGL,NOEFLG,NOEMAT,
     &           ASPROB,SELI,SELJ,SELDI2,SELI2,SELJ2,
     &           MODE,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
C transfers an ambiguity from a better resolved spectrum (class icl1) to
C a spectrum of less resolution (class icl2)
C this is not working yet
C
      IMPLICIT NONE
C I/O
      INCLUDE 'noe.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'numbers.inc'
      INTEGER ICL1,ICL2
      INTEGER NOEFLG(*)
      DOUBLE PRECISION NOEMAT(*)
      DOUBLE PRECISION ASPROB(*), SELDI2(*)
      INTEGER NOEHGL(*)
      INTEGER SELI(*),SELJ(*),SELI2(*),SELJ2(*)
      CHARACTER*4 MODE
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)

C local
      DOUBLE PRECISION SAVEXP, RN
      INTEGER I, J, N, L, NSTART, NSTOP, II, JJ
      INTEGER IAT, JAT, NMAT, IAMB, TEMPI, TEMPJ, NORR
      INTEGER NMAT2, N2, I2
      DOUBLE PRECISION RNL, TOTPRB, ACCLEV, TEMPD
      LOGICAL QDIFF
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
      NMAT=0
C
      DO N=NSTART,NSTOP-1
      SAVEXP=NAVEXP(NOECND(N))
C
C loop over all pairs of atoms belonging to restraint N
      L=0
C
      DO NORR=NOEORR(N)+1, NOEORR(N+1)
C
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
C
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      NMAT=NMAT+1
      L=L+1
      ASPROB(L)=NOEMAT(NMAT)/NMONO(NOECND(N))
      SELI(L)=IAT
      SELJ(L)=JAT
      END IF
C
      END DO
C
      END IF
C
      END DO
C
      END DO
C
C sort the sel arrayS
      DO I=1,L-1
C
      DO J=I+1,L
C[
      IF (ASPROB(J).GT.ASPROB(I)) THEN
      TEMPD=ASPROB(I)
      TEMPI=SELI(I)
      TEMPJ=SELJ(I)
      ASPROB(I)=ASPROB(J)
      SELI(I)=SELI(J)
      SELJ(I)=SELJ(J)
      ASPROB(J)=TEMPD
      SELI(J)=TEMPI
      SELJ(J)=TEMPJ
      END IF
C
      END DO
C
      END DO
C
C calculate the r-6 summed distance^-6
      RN=ZERO
C
      DO I=1,L
      RN=RN+ASPROB(I)
      END DO
C
      TOTPRB=MAX(RN**(-ONE/SAVEXP),RSMALL)
C
C establish how many one needs to exceed acclev or ambcut
C ambcut is a simple distance cutoff
C amblev is the r-6 sum over the accumulated average or minimum distances
      IAMB=0
      ACCLEV=ZERO
      RNL=ZERO
C
      DO WHILE ((ACCLEV.LT.AMBLEV)
     &          .AND.(IAMB.LT.L)
     &          .AND.(RNL.LE.AMBCUT))
      IAMB=IAMB+1
      RNL=MAX(ASPROB(IAMB)**(-ONE/SAVEXP),RSMALL)
      ACCLEV=ACCLEV+(TOTPRB/RNL)**SAVEXP
      END DO
C
      IF (RNL.GT.AMBCUT) IAMB=IAMB-1
C[
      IF (IAMB.EQ.0) THEN
      WRITE (6, '(A, 2I8)') ' no distance < ambcut, one kept, at ',
     &  N, NOEPID(N)
      IAMB=1
      END IF
C
C now we loop over the whole matrix again
C only if the restraint is active
C
      NMAT2=NMAT
C
      DO N2=N+1,NSTOP
C[
      IF ((NOECND(N).EQ.ICL1).AND.(NOECND(N2).EQ.ICL2)
     &     .AND.(NOECV(N2).GT.0).AND.(NOECV(N).GT.0)
     &     .AND.(IAMB.EQ.1)) THEN
      SAVEXP=NAVEXP(NOECND(N2))
C
C loop over all pairs of atoms belonging to restraint N
      L=0
C
      DO NORR=NOEORR(N2)+1, NOEORR(N2+1)
C
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
C
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      NMAT2=NMAT2+1
      L=L+1
      SELDI2(L)=NOEMAT(NMAT2)/NMONO(NOECND(N2))
      SELI2(L)=IAT
      SELJ2(L)=JAT
      END IF
C
      END DO
C
      END IF
C
      END DO
C
      END DO
C
      QDIFF=.TRUE.
      DO I2=1,L
C[
      IF (((SELI(1).EQ.SELI2(I2)).AND.(SELJ(1).EQ.SELJ2(I2)))
     &   .OR.((SELI(1).EQ.SELJ2(I2)).AND.(SELJ(1).EQ.SELI2(I2))))
     &  THEN
      QDIFF=.FALSE.
      END IF
C
      END DO
C[
      IF (.NOT.QDIFF) THEN
C
      DO I2=1,L
      NOEFLG(NMAT2-L+I2)=0
C[
      IF (.NOT.(((SELI(1).EQ.SELI2(I2)).AND.(SELJ(1).EQ.SELJ2(I2)))
     &   .OR.((SELI(1).EQ.SELJ2(I2)).AND.(SELJ(1).EQ.SELI2(I2)))))
     &  THEN
      NOEFLG(NMAT2-L+I2)=1
      END IF
C
      END DO
C
      END IF
C
      END IF
C
      END DO
C
      END DO
C
      RETURN
      END
C
C====6====1=========2=========3=========4=========5=========6=========72
C
      SUBROUTINE RXWRA2(UNIT,IAT,HGLNK,MODE)
C
C writes an unresolved group in NOE ASSIgn format.
C Author: M. Nilges
C
      IMPLICIT NONE
C
C global
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
C
C input/ output
      INTEGER UNIT, IAT, HGLNK(*)
      CHARACTER*4 MODE
C
C local
      INTEGER I1, I2
      CHARACTER*6 PARAN
      CHARACTER*4 PSNAME
      LOGICAL PSOK, QNEXT
C
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      IF (MODE.EQ.'PSEU') THEN
        CALL PSUNAM(IAT,HGLNK,PSNAME,PSOK)
        IF (PSOK) THEN
          WRITE(UNIT,'(8A)')
     & '   (  ','segid "', SEGID(IAT),
     & '" and resid ', RESID(IAT),' and name ', PSNAME,')'
        END IF
      END IF
C
      PARAN='   (( '
      I1=IAT
      I2=I1
      IF ((.NOT.PSOK).OR.(MODE.EQ.'ALL ')) THEN
        CALL NEXTHY(I1, HGLNK,QNEXT)
        DO WHILE (QNEXT)
          CALL NEXTHY(I1, HGLNK,QNEXT)
          WRITE(UNIT,'(8A)')
     & PARAN,'segid "', SEGID(I2),
     & '" and resid ', RESID(I2),' and name ', TYPE(I2),')'
          WRITE(UNIT,'(A)') '     OR '
          PARAN='    ( '
          I2=I1
        END DO
        WRITE(UNIT,'(8A)')
     & PARAN,'segid "', SEGID(I2),
     & '" and resid ', RESID(I2),' and name ', TYPE(I2),'))'
      END IF
C
      RETURN
      END
C
C====6====1=========2=========3=========4=========5=========6=========72
C
      SUBROUTINE PSUNAM(IAT,HGLNK,PSNAME,PSOK)
C
C writes an unresolved group in NOE ASSIgn format.
C Author: M. Nilges
C
      IMPLICIT NONE
C
C global
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
C
C input/ output
      INTEGER IAT, HGLNK(*)
      CHARACTER*4 PSNAME
C
C local
      INTEGER I1, I2, C, POSC, POSB1
      CHARACTER*4 PSUNA1,PSUNA2
      LOGICAL PSOK, QNEXT
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      PSOK=.TRUE.
      I1=IAT
      I2=I1
      PSUNA2=TYPE(I2)
      PSNAME=PSUNA2
      POSC=0
C find first position that is different
        CALL NEXTHY(I1,HGLNK,QNEXT)
        DO WHILE (QNEXT)
          PSUNA1=TYPE(I1)
          DO C=1,LEN(PSUNA1)
            IF ((POSC.EQ.0).AND.(PSUNA2(1:C).NE.PSUNA1(1:C))) THEN
              POSC=C
            END IF
          END DO
          CALL NEXTHY(I1,HGLNK,QNEXT)
        END DO
C do we have a difference at all? we should!
        IF (POSC.EQ.0) PSOK=.FALSE.
C look if the position is a blank
        I1=IAT
        IF (PSOK) THEN
          IF (PSUNA2(POSC:POSC).EQ.' ') PSOK=.FALSE.
          CALL NEXTHY(I1,HGLNK,QNEXT)
          DO WHILE (QNEXT)
            PSUNA1=TYPE(I1)
            IF (PSUNA1(POSC:POSC).EQ.' ') PSOK=.FALSE.
            CALL NEXTHY(I1,HGLNK,QNEXT)
          END DO
        END IF
C make sure it is the last position
        IF (POSC.LT.LEN(PSUNA2)) THEN
          POSB1=POSC+1
          IF (PSUNA2(POSB1:POSB1).NE.' ') PSOK=.FALSE.
        END IF
C change atom name to pseudo
        IF (PSOK) THEN
          PSNAME(POSC:POSC)='%'
        END IF
      RETURN
      END
C
C=========================================================================
C
