C======================================================================
      SUBROUTINE PIG
C
C Pair of Interacting Groups setup parser.
C
C For syntax see main parsing loop "HELP"
C
C Authors: Axel T. Brunger and Thomas Simonson
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C pointer
      INTEGER IFLAGS, JFLAGS
C begin
      IFLAGS=ALLHP(INTEG4(NATOM))
      JFLAGS=ALLHP(INTEG4(NATOM))
      CALL PIG2(HEAP(IFLAGS),HEAP(JFLAGS))
      CALL FREHP(JFLAGS,INTEG4(NATOM))
      CALL FREHP(IFLAGS,INTEG4(NATOM))
      RETURN
      END
C======================================================================
      SUBROUTINE PIG2(IFLAGS,JFLAGS)
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'update.inc'
      INCLUDE 'heap.inc'
      INTEGER IFLAGS(*), JFLAGS(*)
C local
      CHARACTER*4 WDKEEP
      INTEGER ISELCT, JSELCT, I, N
      LOGICAL FOUND, MATCH, FIRST, GETWGT
      DOUBLE PRECISION WGHT
C parameter
      DOUBLE PRECISION ZERO, RAD, ONE
      PARAMETER (RAD=PI/180.0D0, ZERO=0.0D0, ONE=1.0D0)
C begin
C
C flag for setting up the first Pair of Interacting Groups (PIG)
      FIRST=.TRUE.
C parsing
      CALL PUSEND('IGROup>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('IGROup>')
      CALL MISCOM('IGROup>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-igroup')
C
      ELSE IF (WD(1:1).EQ.'?') THEN
C
      WRITE (6,'(A,I3,A)') ' INTE: ',NPIG,
     &' active pairs of interacting groups'
      DO I=1,NPIG
      DO N=1,NENERT
      REN2(N)=PIGWGHT(I,N)
      RENV(N)=PIGAVWT(I,N)
      ENDDO
      WRITE(6,'(A,I6,A)')
     & ' ---------- Interacting Groups: pair ',I,
     & ' ------------------------------------'
      WRITE (6,'(2A)')' | WEIGhts :         ',
     & '                                                          |'
      CALL PRINT2(.FALSE.,REN2)
      WRITE(6,'(2A)') ' --------------------------------------------',
     & '-----------------------------------'
      ENDDO
C
C===
      ELSE IF (WD(1:4).EQ.'INTE') THEN
C
C for each new call to the PIG INTE command, initialize the PIGs;
C for subsequent calls reserve INTERE from heap.
      IF (FIRST) THEN
      CALL PIGRSET('RESE')
      FIRST=.FALSE.
      QPIGRST=.TRUE.
      ELSE
      NPIG=NPIG+1
      IF (NPIG.GT.NPIGMAX) CALL WRNDIE(-5,'IGROup',
     &     'NPIGMAX (cnst.inc) exceeded --> recompile program')
      XINTER=NATOM
      IINTER(NPIG)=ALLHP(INTEG4(XINTER))
      ENDIF
C
      DO I=1,NATOM
      IFLAGS(I)=0
      END DO
      CALL SELCTA(IFLAGS,ISELCT,X,Y,Z,.TRUE.)
      CALL COPYI4(IFLAGS,JFLAGS,NATOM)
      CALL SELCTA(JFLAGS,JSELCT,X,Y,Z,.TRUE.)
C
      IF (NATOM.GT.0) THEN
C
C set up INTERE
      CALL PIGINT(IFLAGS,JFLAGS,HEAP(IINTER(NPIG)))
C reset all nonbonded list flags, i.e., generate exclusion and 1-4
C list, generate atom pair list.
      UPNBLS(NPIG)=.TRUE.
      UPNBSS(NPIG)=.TRUE.
      UPNBEX(NPIG)=.TRUE.
C default weights
      DO I=1,NENERT
      PIGWGHT(NPIG,I)=1.
      ENDDO
C
      END IF
C
C look for Hamiltonian weighting scheme for current PIG
      CALL NEXTWD('IGROup>')
      IF (WD(1:4).NE.'WEIG') THEN
      CALL SAVEWD
      ELSE
C parse weights
      CALL PUSEND('WEIGHT>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('WEIGHT>')
      CALL MISCOM('WEIGHT>',USED)
      IF (.NOT.USED) THEN
      WDKEEP=WD(1:4)
      FOUND=.FALSE.
      GETWGT=.TRUE.
      DO I=1,NENERT
      CALL EQSTWC(ANER(I),4,WDKEEP,4,1,1,MATCH)
      IF (MATCH) THEN
      FOUND=.TRUE.
      IF (GETWGT) THEN
      CALL NEXTF('WEIGHT>',WGHT)
      GETWGT=.FALSE.
      ENDIF
      PIGWGHT(NPIG,I)=WGHT
      ENDIF
      ENDDO
      IF (.NOT.FOUND) CALL CHKEND('WEIGHT>',DONE)
      END IF
      ENDDO
      DONE=.FALSE.
C
      ENDIF
C
C default weights for obsolete perturbation term
      DO I=1,NENERT
      PIGAVWT(NPIG,I)=ZERO
      ENDDO
      ELSE
      CALL CHKEND('IGROup>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      RETURN
      END
C====================================================================
      SUBROUTINE PIGINT(IFLAGS,JFLAGS,INTERE)
C
C Set up INTERE
C The coding for array INTERE is done in such a way that a pair (i,j)
C will be computed if and only if ABS(INTERE(i)+INTERE(j)).le.1
C which is equivalent to the statement ((i in A and j in B) or (i in B
C and j in A)).
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INTEGER IFLAGS(*), JFLAGS(*), INTERE(*)
C local
      INTEGER I
C begin
      DO I=1,NATOM
      IF (IFLAGS(I).EQ.0.AND.JFLAGS(I).EQ.0) THEN
      INTERE(I)=+9999
      ELSE IF (IFLAGS(I).EQ.1.AND.JFLAGS(I).EQ.1) THEN
      INTERE(I)=0
      ELSE IF (IFLAGS(I).EQ.1.AND.JFLAGS(I).EQ.0) THEN
      INTERE(I)=+1
      ELSE IF (IFLAGS(I).EQ.0.AND.JFLAGS(I).EQ.1) THEN
      INTERE(I)=-1
      END IF
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE PIGRSET(ACTION)
C
C INITialize : set up 1st default Pair of Interacting Groups (PIG)
C RESEt : delete all but first Pair of Interacting Groups
C FREE : delete all Pairs of Interacting Groups
C
C Author: Thomas Simonson
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'heap.inc'
      CHARACTER*4 ACTION
C local
      INTEGER I,N
C parameter
      DOUBLE PRECISION ONE, ZERO
      PARAMETER (ONE=1.0D0, ZERO=0.0D0)
C begin
      IF (ACTION.EQ.'INIT') THEN
      NPIG=1
      IINTER(1)=ALLHP(INTEG4(MAXA2))
      CALL FILL4(HEAP(IINTER(1)),MAXA2,0)
      DO I=1,NENERT
      PIGWGHT(1,I)=ONE
      PIGAVWT(1,I)=ZERO
      ENDDO
      DO N=2,NPIGMAX
      IINTER(N)=0
      ENDDO
      ENDIF
      QPIGRST=.FALSE.
C
      IF (ACTION.EQ.'RESE'.OR.ACTION.EQ.'FREE') THEN
      DO N=2,NPIG
      IF (IINTER(N).NE.0) CALL FREHP(IINTER(N),INTEG4(XINTER))
      IINTER(N)=0
      DO I=1,NENERT
      PIGWGHT(N,I)=ZERO
      PIGAVWT(N,I)=ZERO
      ENDDO
      ENDDO
      NPIG=1
      CALL FILL4(HEAP(IINTER(1)),MAXA2,0)
      DO I=1,NENERT
      PIGWGHT(1,I)=ONE
      PIGAVWT(1,I)=ZERO
      ENDDO
      ENDIF
      QPIGRST=.FALSE.
C
      IF (ACTION.EQ.'FREE') THEN
      CALL FREHP(IINTER(1),INTEG4(MAXA2))
      ENDIF
C
      RETURN
      END
