      SUBROUTINE SKELETON(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &  MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &  XRITSY,QHERM,XRCELL,XRINTR,XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &  HPRHOMA,NRHO,IRHO,NMASK)
C   TOPLEVEL SKELETONIZATION ROUTINE
C  Based on on algorithm in
C  PRISM:  Topologically Constrained Phase Refinement for
C  Macromolecular Crystallography,  by David Baker, Christopher Bystroff,
C  Robert J. Fletterick and David A. Agard, Acta Cryst. D49 429-439, 1993.
C  Authors:  Carl Elkin and Axel T. Brunger
C  Modifications: remove duplicate outputs and declare symbols - Jiang
C modifcation, ATB 6/18/08: changed STACK dimension to 0:SEARCHDETH
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'comand.inc'
C I/O
      INTEGER NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP
      INTEGER MAASY,MBASY,MCASY,NAASY,NBASY,NCASY
      INTEGER XRNSYM,XRMSYM,XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4),XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9),XRINTR(3,3)
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*),HPRHOMA,IRHO,NRHO,NMASK
C local
      DOUBLE PRECISION PEAKDEN,B,MINDEN,MAXDEN,EPDEN,MINOUT,SMALLED
      INTEGER MIN_GRAPH,IFROM,I,MAXGRAPHSTATS,MAXNODEFRAC
      INTEGER NNSELE,OUNIT,SEARCHDEPTH
      LOGICAL ERR, COND
      CHARACTER*(WORD_SIZE) CFILE
C pointer
      INTEGER HPRMAP, MPACK
C begin
C skeleton default values
      ERR=.FALSE.
      IFROM=1
      MINOUT=0
      MINDEN = 1.2
      MAXDEN = 6.0
      EPDEN = 2.0
      PEAKDEN = 2.0
      B = 5.0
      MAXGRAPHSTATS=500
      MAXNODEFRAC=2
      MIN_GRAPH = 15
      SEARCHDEPTH=10000
      SMALLED = 0.0005
C
C default output files
      OFILE='OUTPUT'
      CFILE=' '
C
C make default selection:
      IF (HPRHOMA.EQ.0) THEN
         CALL WRNDIE(-5,'SKELETON','No maps are defined.')
         MPACK=0
         ERR=.TRUE.
      ELSE
         MPACK=ALLHP(INTEG4(NRHO))
         CALL XMPPCK0(HEAP(MPACK),NNSELE,
     &  MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,HEAP(HPRHOMA))
      END IF
C
C parsing
      CALL PUSEND('SKELeton>')
      DO WHILE (.NOT.DONE)
         CALL NEXTWD('SKELeton>')
         CALL MISCOM('SKELeton>',USED)
         IF (.NOT.USED) THEN
            IF (WD(1:4).EQ.'HELP') THEN
C
               CALL CNSHELP('cns-xray-skeletonize')
C
C=====================================================================
            ELSE IF (WD(1:4).EQ.'FROM') THEN
               CALL NEXTWD('FROM=')
               IF (WD(1:4).EQ.'=   ') CALL NEXTWD('FROM=')
               IF (WD(1:4).EQ.'?   ') THEN
                  WRITE(6,'(1X,2A)') 'FROM=',XRHONAM(IFROM)
               ELSE
                  COND=.FALSE.
                  DO I=1,XRHONUM
                     IF (WD(1:WDLEN).EQ.XRHONAM(I)) THEN
                        COND=.TRUE.
                        IFROM=I
                     END IF
                  END DO
                  IF (.NOT.COND) THEN
                     WRITE(6,'(3A)') ' %XSKEL-ERR: real space object ',
     &                    WD(1:WDLEN),' not declared.'
                     CALL WRNDIE(-5,'XSKEL','object undeclared.')
                     ERR=.TRUE.
                  END IF
               END IF
C=====================================================================
            ELSE IF (WD(1:4).EQ.'B') THEN
               CALL NEXTF('B=',B)
               B=MAX(RSMALL,B)
C=====================================================================
            ELSE IF (WD(1:4).EQ.'EPDE') THEN
               CALL NEXTF('EPDEn=',EPDEN)
C=====================================================================
            ELSE IF (WD(1:4).EQ.'PEAK') THEN
               CALL NEXTF('PEAKden=',PEAKDEN)
C=====================================================================
            ELSE IF (WD(1:4).EQ.'MINO') THEN
               CALL NEXTF('MINout=',MINOUT)
C=====================================================================
            ELSE IF (WD(1:4).EQ.'SMED') THEN
               CALL NEXTF('SMallElectronDensity=',SMALLED)
C=====================================================================
            ELSE IF (WD(1:4).EQ.'MXNF') THEN
               CALL NEXTI('MaXNodeFrac=',MAXNODEFRAC)
C=====================================================================
            ELSE IF (WD(1:4).EQ.'MAGS') THEN
               CALL NEXTI('MAxGraphStats=',MAXGRAPHSTATS)
               IF (MIN_GRAPH.GE.MAXGRAPHSTATS) THEN
                  WRITE (6,'(2A)') 'WARNING:  value of mingraph must',
     &                 ' be less than the'
                  WRITE (6,'(2A,I4,A,I4)') 'value of maxgraphstats.',
     &                 '  Now, MINGraph=',
     &                 MIN_GRAPH,'  MAGS=',MAXGRAPHSTATS
               END IF
C=====================================================================
            ELSE IF (WD(1:4).EQ.'MXSD') THEN
               CALL NEXTI('MaXSearchDepth=',SEARCHDEPTH)
C=====================================================================
            ELSE IF (WD(1:4).EQ.'MING') THEN
               CALL NEXTI('MINGraph=',MIN_GRAPH)
               IF (MIN_GRAPH.GE.MAXGRAPHSTATS) THEN
                  WRITE (6,'(2A)') 'WARNING:  value of mingraph must',
     &                 ' be less than the'
                  WRITE (6,'(2A,I4,A,I4)') 'value of maxgraphstats.',
     &                 '  Now, MINGraph=',
     &                 MIN_GRAPH,'  MAGS=',MAXGRAPHSTATS
               END IF
C=====================================================================
            ELSE IF (WD(1:4).EQ.'MIND') THEN
               CALL NEXTF('MINDen=',MINDEN)
C=====================================================================
            ELSE IF (WD(1:4).EQ.'MAXD') THEN
               CALL NEXTF('MAXDen=',MAXDEN)
C=====================================================================
            ELSE IF (WD(1:4).EQ.'SELE'.OR.WD(1:1).EQ.'(') THEN
               CALL XMPSELE(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,
     &              NCASY,QHERM,XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &              HPRHOMA,NRHO,NMASK,
     &              XRNSYM,HEAP(MPACK),NNSELE)
C=====================================================================
            ELSE IF (WD(1:4).EQ.'?   ') THEN
               WRITE(6,'(2A)') ' |-----------------------',
     &              ' SKELeton ------------------------------------|'
               WRITE(6,'(A,F4.2,A,F5.2,A,I4,A,I4,A)')' | B =',
     &              B, '   EPDEn=',EPDEN,'   MINGraph=',MIN_GRAPH,
     &              '   MAxGraphStats = ',MAXGRAPHSTATS,'        |'
               WRITE(6,'(A,F5.2,A,F5.2,A,F5.2,A,I6,A)') ' | MAXDen= ',
     &              MAXDEN,'  MINDen= ',MINDEN,'  PEAKden= ',PEAKDEN,
     &              ' MaXSearchDepth= ',SEARCHDEPTH,' |'
               WRITE(6,'(3A,I4,A,F7.5,A)') ' | FROM = ',XRHONAM(IFROM),
     &              '  MaXNodeFrac=',
     &              MAXNODEFRAC,'   SMallED=',SMALLED,
     &              '     |'
               WRITE(6,'(2A)') ' |-----------------------',
     &              '----------------------------------------------|'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PDBF') THEN
      CALL NEXTFI('PDBFile=',CFILE)
C=====================================================================
            ELSE
               CALL CHKEND('SKELeton>',DONE)
            END IF
         END IF
      END DO
      DONE=.FALSE.
C
      CALL ASSFIL(OFILE,OUNIT,'WRITE','FORMATTED',ERROR)
C
      IF (MPACK.NE.0) THEN
         IF (.NOT.ERR.AND.HPRRHO(IFROM).EQ.0) THEN
            WRITE(6,'(3A)') ' %XSKELeton-ERR: real space object ',
     &           XRHONAM(IFROM),' undefined.'
            CALL WRNDIE(-5,'XSKELeton','object undefined.')
            ERR=.TRUE.
         ELSEIF (.NOT.ERR) THEN
            HPRMAP=HPRRHO(IFROM)
            WRITE(OUNIT,'(2A)') ' |-----------------------',
     &           ' SKELeton ------------------------------------|'
            WRITE(OUNIT,'(A,F4.2,A,F5.2,A,I4,A,I4,A)')' | B =',
     &           B, '   EPDEn=',EPDEN,'   MINGraph=',MIN_GRAPH,
     &           '   MAxGraphStats = ',MAXGRAPHSTATS,'        |'
            WRITE(OUNIT,'(A,F5.2,A,F5.2,A,F5.2,A,I6,A)') ' | MAXDen= ',
     &           MAXDEN,'  MINDen= ',MINDEN,'  PEAKden= ',PEAKDEN,
     &           ' MaXSearchDepth= ',SEARCHDEPTH,' |'
            WRITE(OUNIT,'(3A,I4,A,F7.5,A)') ' | FROM = ',XRHONAM(IFROM),
     &           '  MaXNodeFrac=',
     &           MAXNODEFRAC,'   SMallED=',SMALLED,
     &           '     |'
            WRITE(OUNIT,'(2A)') ' |-----------------------',
     &           '----------------------------------------------|'
            CALL XSKELETONIZE (NA,NB,NC,NAASY,NBASY,NCASY,
     &           MAASY,MBASY,MCASY,
     &           HPRMAP,MINDEN,MAXDEN,EPDEN,MIN_GRAPH,
     &           MAXGRAPHSTATS,PEAKDEN,B,XRCELL,XRNSYM,XRMSYM,
     &           XRSYTH,XRSYMM,HPRHOMA,NRHO,XRINTR,HEAP(MPACK),
     &           OUNIT,MINOUT,MAXNODEFRAC,SEARCHDEPTH,SMALLED,CFILE)
         END IF
      END IF
C
      IF (MPACK.NE.0) THEN
      CALL FREHP(MPACK,INTEG4(NRHO))
      END IF
C
      RETURN
      END
C=======================================================================
      SUBROUTINE XSKELETONIZE (NA,NB,NC,NAASY,NBASY,NCASY,MAASY,MBASY,
     &     MCASY,HPRMAP,MINDEN,MAXDEN,EPDEN,MIN_GRAPH,MAXGRAPHSTATS,
     &     PEAKDEN,B,XRCELL,XRNSYM,XRMSYM,XRSYTH,XRSYMM,HPRHOMA,NRHO,
     &     XRINTR,MPACK,OUNIT,MINOUT,MAXNODEFRAC,SEARCHDEPTH,
     &     SMALLED,CFILE)
C     SKELETONIZES MAP
C  Authors: Carl Elkin and Axel T. Brunger
C     density is 3*3 grid of density values
C     minden is the minimum value for a node
C     nodes with density > maxden will not be thinned
C     endpoints with density greater than epden will not be removed
C     min_graph is the minimum graph size which will not be pruned
C     peakden is the denisty on the skeleton itself
C     b is the temperature factor for the gaussian
C     minout is the minimum density value to be outputted
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'xfft.inc'
      INCLUDE 'comand.inc'
C I/O
      INTEGER NA,NB,NC,HPRHOMA,NRHO,OUNIT,HPRMAP,SEARCHDEPTH
      INTEGER NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,MIN_GRAPH
      INTEGER MAXGRAPHSTATS,DENSCUSH
      INTEGER XRNSYM,XRMSYM,XRSYTH,XRSYMM(XRMSYM,3,4)
      INTEGER MPACK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION PEAKDEN,B,MINDEN,MAXDEN,EPDEN,MINOUT,XRINTR(3,3)
      DOUBLE PRECISION SMALLED
      CHARACTER*(*) CFILE
C local
      INTEGER NUMNODES,CUR,PERM,NEIBNUM,MARK,REMOVED,NEIB,NODE,THINCOUNT
      INTEGER NREMOVED,NREMAIN,MAXNODEFRAC,STACK,BOXDEN
      INTEGER XNPTMP,YNPTMP,ZNPTMP,DENTMP,XNODEPOS,YNODEPOS,ZNODEPOS
      INTEGER AMAX,BMAX,CMAX,CUSHSIZE,CUNIT
      DOUBLE PRECISION XRCELL(9),RCUT
C begin
C   allocate HEAP space for variables
      XNPTMP=ALLHP(INTEG4(NRHO/MAXNODEFRAC))
      YNPTMP=ALLHP(INTEG4(NRHO/MAXNODEFRAC))
      ZNPTMP=ALLHP(INTEG4(NRHO/MAXNODEFRAC))
C   Convert from temperature factor to B factor
      B=MAX(RSMALL,10./B)
C   Calculate Rcut-- the radius at which the Gaussian ceases to be signifigant
      RCUT=SQRT (ABS(-LOG (SMALLED)/B))
C   Calculate size of cushion
C   First, compute norm of A*, B*, C* in grid
      FTAS=XRCELL(7)*NA
      FTBS=XRCELL(8)*NB
      FTCS=XRCELL(9)*NC
C   Then, determine AMAX, BMAX, CMAX which are the size of the cussion that
C   surrounds the box in the three spacial dimensions.
      AMAX=MAX(MIN(NA/2,NINT(RCUT*FTAS)),1)
      BMAX=MAX(MIN(NB/2,NINT(RCUT*FTBS)),1)
      CMAX=MAX(MIN(NC/2,NINT(RCUT*FTCS)),1)
C Identify all nodes, and sort them by density
      CALL XSKIDNODE(NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,MINDEN,
     &     HEAP(HPRMAP),NUMNODES,HEAP(XNPTMP),HEAP(YNPTMP),
     &     HEAP(ZNPTMP),MPACK,NRHO,MAXNODEFRAC,OUNIT)
      XNODEPOS=ALLHP(INTEG4(NUMNODES))
      YNODEPOS=ALLHP(INTEG4(NUMNODES))
      ZNODEPOS=ALLHP(INTEG4(NUMNODES))
      DENTMP=ALLHP(IREAL8(NUMNODES))
      PERM=ALLHP(INTEG4(NUMNODES))
      CALL XSKDENSORT (NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,NUMNODES,
     &     HEAP(HPRMAP),HEAP(DENTMP),HEAP(XNODEPOS),HEAP(YNODEPOS),
     &     HEAP(ZNODEPOS),HEAP(XNPTMP),HEAP(YNPTMP),HEAP(ZNPTMP),
     &     HEAP(PERM),NRHO,MAXNODEFRAC)
C Free some heap space and allocate more
      CALL FREHP(XNPTMP,INTEG4(NRHO/MAXNODEFRAC))
      CALL FREHP(YNPTMP,INTEG4(NRHO/MAXNODEFRAC))
      CALL FREHP(ZNPTMP,INTEG4(NRHO/MAXNODEFRAC))
      CALL FREHP(DENTMP,IREAL8(NUMNODES))
      CALL FREHP(PERM,INTEG4(NUMNODES))
      NEIBNUM=ALLHP(INTEG4(NUMNODES))
      MARK=ALLHP(ILOGIC(NUMNODES))
      REMOVED=ALLHP(ILOGIC(NUMNODES))
      NEIB=ALLHP(INTEG4(NUMNODES*6))
      NODE=ALLHP(INTEG4(NRHO))
      STACK=ALLHP(INTEG4(SEARCHDEPTH+1))
      NREMOVED=ALLHP(INTEG4(MAXGRAPHSTATS))
      NREMAIN=ALLHP(INTEG4(MAXGRAPHSTATS))
C put a list of neighbors of nodes into 2-d array neib
      CALL XSKIDNEIGHBOR (NUMNODES,HEAP(XNODEPOS),HEAP(YNODEPOS),
     &     HEAP(ZNODEPOS),NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,
     &     NCASY,HEAP(NODE),XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &     HEAP(HPRHOMA),MPACK,HEAP(NEIB),HEAP(NEIBNUM))
      CALL FREHP(NODE,INTEG4(NRHO))
C thin graph; i.e., remove all nodes not necesary to preserve connectivity
      CALL XSKTHIN (NUMNODES,HEAP(NEIBNUM),HEAP(MARK),
     &     HEAP(REMOVED),HEAP(NEIB),EPDEN,HEAP(XNODEPOS),
     &     HEAP(YNODEPOS),HEAP(ZNODEPOS),HEAP(HPRMAP),NAASY,NBASY,
     &     NCASY,MAASY,MBASY,MCASY,MAXDEN,THINCOUNT,HEAP(STACK),
     &     SEARCHDEPTH,OUNIT)
C prune all graphs which are too small
      CALL XSKPRUNE (NUMNODES,HEAP(NEIBNUM),HEAP(MARK),
     &     HEAP(REMOVED),HEAP(NEIB),CUR,MIN_GRAPH,MAXGRAPHSTATS,
     &     THINCOUNT,HEAP(NREMOVED),HEAP(NREMAIN),HEAP(STACK),
     &     SEARCHDEPTH,OUNIT)
      CALL FREHP (NEIBNUM,INTEG4(NUMNODES))
      CALL FREHP (MARK,ILOGIC(NUMNODES))
      CUSHSIZE=(NAASY-MAASY+(2*AMAX)+1)*(NBASY-MBASY+(2*BMAX)+1)*
     &     (NCASY-MCASY+(2*CMAX)+1)
      DENSCUSH = ALLHP(IREAL4(CUSHSIZE))
      BOXDEN = ALLHP(IREAL8(((2*AMAX)+1)*((2*BMAX)+1)*((2*CMAX)+1)))
C add cylindrical gaussians to the skeleton
      CALL XSKADDGAUSSIAN (NA,NB,NC,NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,
     &     NUMNODES,HEAP(REMOVED),HEAP(XNODEPOS),HEAP(YNODEPOS),
     &     HEAP(ZNODEPOS),HEAP(HPRMAP),HEAP(DENSCUSH),PEAKDEN,B,XRNSYM,
     &     XRMSYM,XRSYTH,XRSYMM,HEAP(HPRHOMA),MPACK,XRINTR,
     &     AMAX,BMAX,CMAX,HEAP(BOXDEN))
C output some stats
      CALL XSKSTATS (NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,HEAP(REMOVED),
     &     NUMNODES,HEAP(XNODEPOS),HEAP(YNODEPOS),HEAP(ZNODEPOS),
     &     HEAP(HPRMAP),MINDEN,OUNIT,MINOUT)
C write skeleton as coordinates if requested
      IF (CFILE.NE.' ') THEN
        CALL ASSFIL(CFILE,CUNIT,'WRITE','FORMATTED',ERROR)
        CALL XSKCOOR (NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,
     &                HEAP(REMOVED),NUMNODES,HEAP(XNODEPOS),
     &                HEAP(YNODEPOS),HEAP(ZNODEPOS),HEAP(HPRMAP),
     &                CUNIT,XRINTR,NA,NB,NC)
      END IF
C free the heap
      CALL FREHP (REMOVED,ILOGIC(NUMNODES))
      CALL FREHP (XNODEPOS,INTEG4(NUMNODES))
      CALL FREHP (YNODEPOS,INTEG4(NUMNODES))
      CALL FREHP (ZNODEPOS,INTEG4(NUMNODES))
      CALL FREHP (DENSCUSH,IREAL4(CUSHSIZE))
      CALL FREHP (BOXDEN,IREAL8(((2*AMAX)+1)*((2*BMAX)+1)*((2*CMAX)+1)))
      CALL FREHP (NEIB,INTEG4(NUMNODES*6))
      CALL FREHP (NREMOVED,INTEG4(MAXGRAPHSTATS))
      CALL FREHP (NREMAIN,INTEG4(MAXGRAPHSTATS))
      CALL FREHP (STACK,INTEG4(SEARCHDEPTH+1))
      RETURN
      END
C=======================================================================
      SUBROUTINE XSKIDNODE(NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,MINDEN,
     &          DENSITY,NUMNODES,XNPTMP,YNPTMP,ZNPTMP,MPACK,NRHO,
     &          MAXNODEFRAC,OUNIT)
C     counts the total number of nodes, (grid points with density > minden)
C     and stores their x,y,z, coordinates in xnptmp,ynptmp,znptmp
C Author:  Carl Elkin
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C I/O
      INTEGER MAXNODEFRAC,OUNIT
      INTEGER MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,NUMNODES,NRHO
      INTEGER XNPTMP (NRHO/MAXNODEFRAC)
      INTEGER YNPTMP (NRHO/MAXNODEFRAC)
      INTEGER ZNPTMP (NRHO/MAXNODEFRAC)
      INTEGER MPACK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION MINDEN
      REAL DENSITY (MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBPREC
C local
      INTEGER I,J,K
C begin
      NUMNODES = 0
      DO I = MAASY,NAASY
         DO J = MBASY,NBASY
            DO K = MCASY,NCASY
               IF ((DENSITY(I,J,K).GT.MINDEN)
     &              .AND.(MPACK(I,J,K).EQ.1)) THEN
                  NUMNODES=NUMNODES+1
                  XNPTMP(NUMNODES)=I
                  YNPTMP(NUMNODES)=J
                  ZNPTMP(NUMNODES)=K
                  IF (NUMNODES*MAXNODEFRAC.GT.NRHO) THEN
                    WRITE (6,'(A)')'Too many nodes!  Adjust parameters.'
                    CALL WRNDIE (-5,'XSKEleton',
     &                   'Too many nodes.')
                  END IF
               END IF
            END DO
         END DO
      END DO
      WRITE(OUNIT,'(I6,A)') NUMNODES,' nodes were identified.'
      DBPREC=NUMNODES
      CALL DECLAR( 'NUMNODES','DP',' ',DBCOMP,DBPREC)
      RETURN
      END
C=======================================================================
      SUBROUTINE XSKDENSORT (NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,
     &     NUMNODES,DENSITY,DENTMP,XNODEPOS,YNODEPOS,ZNODEPOS,
     &     XNPTMP,YNPTMP,ZNPTMP,PERM,NRHO,MAXNODEFRAC)
C     sorts the nodes by density, arranged from least to greatest
C     outputs stored coordinates in xnodepos,ynodepos,znodepos
C Author:  Carl Elkin
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C I/O
      INTEGER NUMNODES,NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,NRHO
      INTEGER MAXNODEFRAC
      INTEGER PERM(NUMNODES),XNODEPOS(NUMNODES),YNODEPOS(NUMNODES)
      INTEGER ZNODEPOS(NUMNODES)
      INTEGER XNPTMP(NRHO/MAXNODEFRAC),YNPTMP(NRHO/MAXNODEFRAC)
      INTEGER ZNPTMP(NRHO/MAXNODEFRAC)
      REAL DENSITY (MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION DENTMP(NUMNODES)
C local
      EXTERNAL ORDRRR
      INTEGER I
C begin
      DO I = 1,NUMNODES
         DENTMP(I) = DENSITY(XNPTMP(I),YNPTMP(I),ZNPTMP(I))
      END DO
C sort T-FUNCtion
      CALL SORTP(NUMNODES,PERM,ORDRRR,DENTMP,0,0,0,0,0,0,0)
      DO I = 1,NUMNODES
         XNODEPOS(I)=XNPTMP(PERM(I))
         YNODEPOS(I)=YNPTMP(PERM(I))
         ZNODEPOS(I)=ZNPTMP(PERM(I))
      END DO
      RETURN
      END
C=======================================================================
      SUBROUTINE XSKIDNEIGHBOR (NUMNODES,XNODEPOS,YNODEPOS,ZNODEPOS,
     &     NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,NODE,XRNSYM,
     &     XRMSYM,XRSYTH,XRSYMM,RHOMASK,MPACK,NEIB,NEIBNUM)
C     identifies as neighbors nodes at X+-1,Y,Z;X,Y+-1,Z;or X,Y,Z+-1
C     defines neib [nodenum,i] as the ith neighbor of node @nodenum
C     defines neibnum [i] as the number of neighbors of node i
C Author:  Carl Elkin
      IMPLICIT NONE
C
      INTEGER NUMNODES
      INTEGER XNODEPOS(NUMNODES),YNODEPOS(NUMNODES),ZNODEPOS(NUMNODES)
      INTEGER NEIB(NUMNODES,6),NEIBNUM(NUMNODES)
      INTEGER NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY
      INTEGER NODE (MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER MPACK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNSYM,XRMSYM,XRSYTH,XRSYMM(XRMSYM,3,4)
C local
      LOGICAL COND
      INTEGER SYM
      INTEGER COUNT,I,J,K,XNTMP(6),YNTMP(6),ZNTMP(6)
      INTEGER IIS,JJS,KKS,II2,JJ2,KK2
C begin
      DO I = MAASY,NAASY
         DO J = MBASY,NBASY
            DO K = MCASY,NCASY
               NODE(I,J,K)=0
            END DO
         END DO
      END DO
      DO I = 1,NUMNODES
         NODE(XNODEPOS(I),YNODEPOS(I),ZNODEPOS(I))=I
      END DO
      DO I = 1,6
         XNTMP(I)=0
         YNTMP(I)=0
         ZNTMP(I)=0
      END DO
      XNTMP(1)= 1
      XNTMP(2)=-1
      YNTMP(3)= 1
      YNTMP(4)=-1
      ZNTMP(5)= 1
      ZNTMP(6)=-1
      DO I = 1,NUMNODES
         COUNT=0
         DO J=1,6
            II2=XNODEPOS(I)+XNTMP(J)
            JJ2=YNODEPOS(I)+YNTMP(J)
            KK2=ZNODEPOS(I)+ZNTMP(J)
C Apply all symmetry operators to grid point until the symmetry-related grid
C point falls into the asymmetric unit
            COND=(II2.GE.MAASY.AND.II2.LE.NAASY.AND.
     &           JJ2.GE.MBASY.AND.JJ2.LE.NBASY.AND.
     &           KK2.GE.MCASY.AND.KK2.LE.NCASY)
            IF (COND) THEN
               COND=RHOMASK(II2,JJ2,KK2).GT.0
            END IF
            IIS=II2
            JJS=JJ2
            KKS=KK2
            IF (.NOT.(COND)) THEN
               SYM=1
               DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
                  IIS=MOD(XRSYMM(SYM,1,1)*II2+
     &                 XRSYMM(SYM,1,2)*JJ2+XRSYMM(SYM,1,3)*KK2
     &                    +(XRSYMM(SYM,1,4)*NA)/XRSYTH+10000*NA
     &                 -MAASY,NA)+MAASY
                  JJS=MOD(XRSYMM(SYM,2,1)*II2+XRSYMM(SYM,2,2)
     &                 *JJ2+XRSYMM(SYM,2,3)*KK2
     &                 +(XRSYMM(SYM,2,4)*NB)/XRSYTH+10000*NB
     &                 -MBASY,NB)+MBASY
                  KKS=MOD(XRSYMM(SYM,3,1)*II2+XRSYMM(SYM,3,2)*
     &                 JJ2+XRSYMM(SYM,3,3)*KK2
     &                 +(XRSYMM(SYM,3,4)*NC)/XRSYTH+10000*NC
     &                 -MCASY,NC)+MCASY
                  COND=((((IIS.GE.MAASY).AND.(IIS.LE.NAASY)).AND.
     &                 ((JJS.GE.MBASY).AND.(JJS.LE.NBASY))).AND.
     &                 ((KKS.GE.MCASY).AND.(KKS.LE.NCASY)))
                  IF (COND) THEN
                     COND=RHOMASK(IIS,JJS,KKS).GT.0
                  END IF
                  SYM=SYM+1
               END DO
            END IF
            IF (.NOT.COND) CALL WRNDIE
     &           (-5,'XSKELeton','Internal error no. 2')
            IF (NODE(IIS,JJS,KKS).GT.0) THEN
               COUNT=COUNT+1
               NEIB (I,COUNT) = NODE(IIS,JJS,KKS)
            END IF
         END DO
         NEIBNUM(I) = COUNT
      END DO
      RETURN
      END
C=======================================================================
      SUBROUTINE XSKTHIN (NUMNODES,NEIBNUM,MARK,REMOVED,NEIB,EPDEN,
     &     XNODEPOS,YNODEPOS,ZNODEPOS,DENSITY,NAASY,NBASY,NCASY,
     &     MAASY,MBASY,MCASY,MAXDEN,THINCOUNT,STACK,SEARCHDEPTH,OUNIT)
C     Thins graphs by removing all nodes with density < maxden whose removal
C     does not break chains
C Author:  Carl Elkin
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'comand.inc'
C I/O
      INTEGER NUMNODES,SEARCHDEPTH,OUNIT
      INTEGER NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,NEIBNUM (NUMNODES)
      INTEGER STACK(0:SEARCHDEPTH)
      DOUBLE PRECISION EPDEN,MAXDEN
      INTEGER NEIB(NUMNODES,6),THINCOUNT
      LOGICAL MARK(NUMNODES),REMOVED(NUMNODES)
      INTEGER XNODEPOS(NUMNODES),YNODEPOS(NUMNODES),ZNODEPOS(NUMNODES)
      REAL DENSITY (MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
C local
      INTEGER NTHIN,UNREMOVED,I,NMARKED,CUR
      LOGICAL CUTPT
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBPREC
C cutpt is true if node is a cutpoint, i.e., if removing it will break a chain
C     INITIALIZE
      DO I = 1,NUMNODES
         REMOVED(I)=.FALSE.
      END DO
      THINCOUNT=0
C   For each node, we remove the node, then go to a neighbor and mark all nodes
C   in the same chain.
      DO NTHIN = 1,NUMNODES
         IF (DENSITY (XNODEPOS(NTHIN),YNODEPOS(NTHIN),
     &        ZNODEPOS(NTHIN)).LT.MAXDEN) THEN
            UNREMOVED=0
            DO I = 1, NEIBNUM(NTHIN)
               IF ((.NOT.REMOVED(NEIB(NTHIN,I)))
     &              .AND. (UNREMOVED.LE.1)) THEN
                  CUR = NEIB(NTHIN,I)
                  UNREMOVED = UNREMOVED + 1
               END IF
            END DO
            REMOVED (NTHIN) = .TRUE.
            IF (UNREMOVED.GT.1) THEN
C             IF ITS NOT AN ENDPOINT
               CALL XSKSEARCH (NUMNODES,NEIBNUM,MARK,REMOVED,NEIB,
     &              CUR,NMARKED,STACK,SEARCHDEPTH)
C   If one of the other neighbors is unmarked, the original node is a cut
C   point and is not removed
               CUTPT = .FALSE.
               DO I = 1, NEIBNUM(NTHIN)
                  IF (.NOT.REMOVED(NEIB(NTHIN,I))) THEN
                     CUTPT = (CUTPT.OR.(.NOT.(MARK(NEIB(NTHIN,I)))))
                  END IF
               END DO
               REMOVED (NTHIN) = .NOT.CUTPT
            ELSE
C   If the node is an endpoint, it is only removed if density < epden
               REMOVED (NTHIN)=(.NOT.(DENSITY(XNODEPOS(NTHIN),
     &              YNODEPOS(NTHIN),ZNODEPOS(NTHIN)).GT. EPDEN))
            END IF
            IF (REMOVED (NTHIN)) THINCOUNT=THINCOUNT+1
         END IF
      END DO
      WRITE (OUNIT,'(I6,A)')
     &  THINCOUNT,' nodes were removed during thinning.'
      WRITE (OUNIT,'(I6,A)')(NUMNODES-THINCOUNT),
     &     ' nodes remain after thinning.'
      DBPREC=NUMNODES-THINCOUNT
      CALL DECLAR( 'THINNODES','DP',' ',DBCOMP,DBPREC)
      RETURN
      END
C=======================================================================
      SUBROUTINE XSKPRUNE (NUMNODES,NEIBNUM,MARK,REMOVED,NEIB,CURSAVE,
     &     MIN_GRAPH,MAXGRAPHSTATS,THINCOUNT,NREMOVED,NREMAIN,STACK,
     &     SEARCHDEPTH,OUNIT)
C     Removes all graphs with fewer than min_graph points
C Author:  Carl Elkin
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'comand.inc'
C I/O
      INTEGER NUMNODES,MAXGRAPHSTATS,SEARCHDEPTH,OUNIT
      INTEGER NEIBNUM (NUMNODES),CURSAVE,MIN_GRAPH,THINCOUNT,BIGCHAIN
      INTEGER NEIB(NUMNODES,6),NREMOVED(MAXGRAPHSTATS)
      INTEGER STACK(0:SEARCHDEPTH)
      INTEGER NREMAIN(MAXGRAPHSTATS)
      LOGICAL MARK(NUMNODES),REMOVED(NUMNODES)
C local
      INTEGER I,NMARKED,NPRUN,PRUNECOUNT
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBPREC
      BIGCHAIN = 0
      DO I = 1,MAXGRAPHSTATS
         NREMOVED(I)=0
         NREMAIN(I)=0
      END DO
C For each node, we count the number of nodes in the chain, and remove
C the whole chain if the number of nodes is less than min_graph
      DO NPRUN = 1,NUMNODES
         IF (.NOT.REMOVED(NPRUN)) THEN
            CALL XSKSEARCH (NUMNODES,NEIBNUM,MARK,REMOVED,NEIB,
     &           NPRUN,NMARKED,STACK,SEARCHDEPTH)
            IF (WRNLEV.GE.10) THEN
               WRITE (6,'(A,I5,A,I5)') 'Chain about node#',NPRUN,
     &              ' was of length ',NMARKED
            END IF
            IF (NMARKED.LT.MIN_GRAPH) THEN
               CALL XSKDESTROY (NUMNODES,NEIBNUM,REMOVED,NEIB,NPRUN,
     &              STACK,SEARCHDEPTH)
               NREMOVED(NMARKED)=NREMOVED(NMARKED)+1
            ELSE
               IF (NMARKED.LE.MAXGRAPHSTATS)
     &              NREMAIN(NMARKED)=NREMAIN(NMARKED)+1
            END IF
            IF (NMARKED.GT.BIGCHAIN) BIGCHAIN=NMARKED
         END IF
      END DO
C output some stats
      PRUNECOUNT = NREMOVED(1)
      IF (NREMOVED(1).GT.1) THEN
         WRITE (OUNIT,'(I6,A)')
     &     NREMOVED(1),' isolated nodes were pruned.'
      ELSE
         IF (NREMOVED(1).EQ.1) THEN
            WRITE (OUNIT,'(A)') '     1 isolated node was pruned.'
         END IF
      END IF
      DO I = 2,(MIN(MAXGRAPHSTATS,MIN_GRAPH-1))
         PRUNECOUNT=PRUNECOUNT+(NREMOVED(I)*I)
         IF (NREMOVED(I).GT.1) THEN
            WRITE (OUNIT,'(I6,A,I3,A)') NREMOVED(I),' chains with ',I
     &           ,' nodes each were pruned.'
         ELSE
            IF (NREMOVED(I).EQ.1) THEN
               WRITE (OUNIT,'(I6,A,I3,A)') NREMOVED(I),' chain  with ',I
     &              ,' nodes was pruned.'
            END IF
         END IF
      END DO
      DO I = (MIN(MAXGRAPHSTATS,(MIN_GRAPH))),MAXGRAPHSTATS
         IF ((NREMAIN(I)/I).GT.1) THEN
            WRITE (OUNIT,'(I6,A,I3,A)') NREMAIN(I)/I,
     &           ' chains with ',I,' nodes each remain.'
         ELSE
            IF ((NREMAIN(I)/I).EQ.1) THEN
               WRITE (OUNIT,'(I6,A,I3,A)') NREMAIN(I)/I,
     &           ' chain  with ',I,' nodes remains.'
            END IF
         END IF
      END DO
      WRITE (OUNIT,'(A,I5,A)') 'A total of ',PRUNECOUNT,
     &     ' nodes were removed during pruning.'
      IF (NUMNODES.EQ.(THINCOUNT+PRUNECOUNT)) THEN
         WRITE (6,'(A)')
     &     'WARNING:  all nodes were removed, and none remain.'
         WRITE (6,'(A)') 'Try adjusting parameters.'
         CALL WRNDIE (-5,'XSKEleton','Density map will be featureless.')
      END IF
      WRITE (OUNIT,'(A,I6,A,F5.1,A,I6,A)')
     &     'The largest chain has',BIGCHAIN,' nodes,',
     &     (REAL(BIGCHAIN*100)/
     &     REAL(NUMNODES-(THINCOUNT+PRUNECOUNT))),'% of the',
     &     (NUMNODES-(THINCOUNT+PRUNECOUNT)),
     &     ' remaining nodes'
      DBPREC=NUMNODES-(THINCOUNT+PRUNECOUNT)
      CALL DECLAR( 'PRUNENODES','DP',' ',DBCOMP,DBPREC)
      DBPREC=BIGCHAIN
      CALL DECLAR( 'BIGCHAIN','DP',' ',DBCOMP,DBPREC)
      RETURN
      END
C=========================================================================
      SUBROUTINE XSKSEARCH (NUMNODES,NEIBNUM,MARK,REMOVED,NEIB,
     &     CURSAVE,NMARKED,STACK,SEARCHDEPTH)
C     Marks every node in the same graph and counts nodes marked
C Author:  Carl Elkin
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
C I/O
      INTEGER NUMNODES,SEARCHDEPTH,CURSAVE,NMARKED
      INTEGER NEIBNUM(NUMNODES),NEIB(NUMNODES,6),STACK(0:SEARCHDEPTH)
      LOGICAL MARK(NUMNODES),REMOVED(NUMNODES)
C local
      INTEGER PLACE,UNMARK,LOWUNADJ,I,CUR
C begin
      PLACE = 0
      CUR = CURSAVE
      DO I = 1,NUMNODES
         MARK (I) = .FALSE.
      END DO
      NMARKED = 0
      IF (NEIBNUM(CUR).EQ.0)  THEN
         MARK (CUR) = .TRUE.
         NMARKED = NMARKED + 1
      ELSE
         DO WHILE (PLACE .GT. -1)
            IF (.NOT.(MARK(CUR))) THEN
               NMARKED = NMARKED + 1
               MARK(CUR) = .TRUE.
            END IF
C   Count unmarked neighbors, and find the lowest unmarked adjacent one
            UNMARK = 0
            DO I = 1, NEIBNUM(CUR)
               IF (((.NOT.MARK(NEIB(CUR,I))).AND.
     &              (.NOT.REMOVED(NEIB(CUR,I))))
     &              .AND. (UNMARK.LE.1)) THEN
C  WE ONLY CARE IF THERE ARE ZERO, ONE OR MORE UNMARKED ONES
                  IF (UNMARK.EQ.0) THEN
                     LOWUNADJ = I
                  END IF
                  UNMARK = UNMARK + 1
               END IF
            END DO
C       IF MORE THEN ONE ADJACENT POINT
            IF (NEIBNUM (CUR) .GT. 1) THEN
C  If there is still at least one node we haven't seen before...
               IF (UNMARK .GT. 0) THEN
C
C  If at a branch, put the node not taken onto the top of the stack
                  IF ((UNMARK.GT.1).AND.((STACK (PLACE).NE.CUR).OR.
     &                 (PLACE.EQ.0))) THEN
                     PLACE = PLACE + 1
                     IF (PLACE.GE.SEARCHDEPTH) THEN
                        CALL WRNDIE (-5,'XSKEleton',
     &                       'MaXSearchDepth too small.')
                     END IF
                     STACK (PLACE) = CUR
                  END IF
C  Move on to the next node
                  CUR = NEIB(CUR,LOWUNADJ)
               ELSE
C                   If all other nodes have been seen before
                  CUR = STACK (PLACE)
                  PLACE = PLACE - 1
               END IF
            ELSE
C                If it is an endpoint
               IF (UNMARK.EQ.0) THEN
                  IF (PLACE .GE. 0) THEN
                     CUR = STACK (PLACE)
                  END IF
                  PLACE = PLACE - 1
               ELSE
                  CUR = NEIB(CUR,LOWUNADJ)
               END IF
            END IF
         END DO
      END IF
      RETURN
      END
C=======================================================================
      SUBROUTINE XSKDESTROY (NUMNODES,NEIBNUM,REMOVED,NEIB,CURSAVE,
     &     STACK,SEARCHDEPTH)
C Removes all the nodes in the graph of node cursave
C Author:  Carl Elkin
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
C I/O
      INTEGER NUMNODES,SEARCHDEPTH,CURSAVE
      INTEGER NEIBNUM(NUMNODES),NEIB(NUMNODES,6),STACK(0:SEARCHDEPTH)
      LOGICAL REMOVED(NUMNODES)
C local
      INTEGER PLACE,UNREMOVED,LOWUNADJ,I,CUR
C begin
      PLACE = 0
      CUR = CURSAVE
      IF (NEIBNUM(CUR).EQ.0)  THEN
         REMOVED (CUR) = .TRUE.
      ELSE
         DO WHILE (PLACE .GT. -1)
            IF (.NOT.(REMOVED(CUR))) THEN
               REMOVED(CUR) = .TRUE.
            END IF
C   COUNT UNREMOVED NEIGHBORS, AND FIND THE LOWEST UNREMOVED ADJACENT ONE
            UNREMOVED = 0
            DO I = 1, NEIBNUM(CUR)
               IF((.NOT.REMOVED(NEIB(CUR,I))).AND.(UNREMOVED.LE.1))THEN
C  WE ONLY CARE IF THERE ARE ZERO, ONE OR MORE UNMARKED ONES
                  IF (UNREMOVED.EQ.0) THEN
                     LOWUNADJ = I
                  END IF
                  UNREMOVED = UNREMOVED + 1
               END IF
            END DO
C       IF MORE THEN ONE ADJACENT POINT
            IF (NEIBNUM (CUR) .GT. 1) THEN
C  If there is at least one node we haven't removed
               IF (UNREMOVED .GT. 0) THEN
C  If at a branch, put the node not removed on top of the stack
                  IF ((UNREMOVED.GT.1).AND.((STACK(PLACE).NE.CUR)
     &                 .OR.(PLACE.EQ.0))) THEN
                     PLACE = PLACE + 1
                     STACK (PLACE) = CUR
                  END IF
C  Move on to the next node
                  CUR = NEIB(CUR,LOWUNADJ)
               ELSE
C                   IF ALL OTHER NODES HAVE BEEN REMOVED
                  CUR = STACK (PLACE)
                  PLACE = PLACE - 1
               END IF
            ELSE
C                IF IT IS AN ENDPOINT
               IF (UNREMOVED.EQ.0) THEN
                  IF (PLACE .GE. 0) THEN
                     CUR = STACK (PLACE)
                  END IF
                  PLACE = PLACE - 1
               ELSE
                  CUR = NEIB(CUR,LOWUNADJ)
               END IF
            END IF
         END DO
      END IF
      RETURN
      END
C=======================================================================
      SUBROUTINE XSKADDGAUSSIAN (NA,NB,NC,NAASY,NBASY,NCASY,MAASY,MBASY,
     &     MCASY,NUMNODES,REMOVED,XNODEPOS,YNODEPOS,ZNODEPOS,DENSITY,
     &     DENSCUSH,PEAKDEN,B,XRNSYM,XRMSYM,XRSYTH,XRSYMM,RHOMASK,
     &     MPACK,XRINTR,AMAX,BMAX,CMAX,BOXDEN)
C  Adds cylindrical gaussians around skeleton
C  We have a cushion which extends four units beyond the grid
C Author:  Carl Elkin
      IMPLICIT NONE
C
      INTEGER NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,NUMNODES,NA,NB,NC
      INTEGER XNODEPOS(NUMNODES),YNODEPOS(NUMNODES),ZNODEPOS(NUMNODES)
      INTEGER XRNSYM,XRMSYM,XRSYTH,XRSYMM(XRMSYM,3,4)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER MPACK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER AMAX,BMAX,CMAX
      LOGICAL REMOVED(NUMNODES)
      REAL DENSITY (MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL DENSCUSH (MAASY-AMAX:NAASY+AMAX,
     &     MBASY-BMAX:NBASY+BMAX,MCASY-CMAX:NCASY+CMAX)
      DOUBLE PRECISION PEAKDEN,B
      DOUBLE PRECISION XRINTR(3,3)
      DOUBLE PRECISION BOXDEN(-AMAX:AMAX,-BMAX:BMAX,-CMAX:CMAX)
C local
      DOUBLE PRECISION DISTSQ, TEMP
      LOGICAL COND
      INTEGER I,II,JJ,KK,II2,JJ2,KK2,IIS,JJS,KKS,SYM
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C calculate density (ii,jj,kk) away from central point
      DO II = -AMAX,AMAX
         DO JJ = -BMAX,BMAX
            DO KK = -CMAX,CMAX
               CALL XSKSHELLD(II,JJ,KK,XRINTR,DISTSQ,NA,NB,NC)
               BOXDEN(II,JJ,KK) = PEAKDEN * EXP(-B*DISTSQ)
            END DO
         END DO
      END DO
C initialize denscush and density
      DO II = MAASY-AMAX,NAASY+AMAX
         DO JJ = MBASY-BMAX,NBASY+BMAX
            DO KK = MCASY-CMAX,NCASY+CMAX
               DENSCUSH(II,JJ,KK)=ZERO
               IF ((II.GE.MAASY).AND.(II.LE.NAASY)
     &              .AND.(JJ.GE.MBASY).AND.(JJ.LE.NBASY)
     &              .AND.(KK.GE.MCASY).AND.(KK.LE.NCASY))
     &              DENSITY (II,JJ,KK)=ZERO
            END DO
         END DO
      END DO
C  for each node, add on the gaussian to surrounding points
      DO I = 1,NUMNODES
         IF (.NOT.REMOVED(I)) THEN
            DO II = -AMAX,AMAX
               II2 = XNODEPOS(I) + II
               DO JJ = -BMAX,BMAX
                  JJ2 = YNODEPOS(I) + JJ
                  DO KK = -CMAX,CMAX
                     KK2 = ZNODEPOS(I) + KK
                     TEMP=DENSCUSH(II2,JJ2,KK2)
                     DENSCUSH(II2,JJ2,KK2) =
     &                    MAX(TEMP,BOXDEN(II,JJ,KK))
                  END DO
               END DO
            END DO
         END IF
      END DO
C Now, fold the cushion back onto the old density
C Apply all symmetry operators to grid point until the symmetry-related grid
C  point falls into the asymmetric unit
      DO II2 = MAASY-AMAX,NAASY+AMAX
         DO JJ2 = MBASY-BMAX,NBASY+BMAX
            DO KK2 = MCASY-CMAX,NCASY+CMAX
               IIS = II2
               JJS = JJ2
               KKS = KK2
               COND=(IIS.GE.MAASY.AND.IIS.LE.NAASY.AND.
     &              JJS.GE.MBASY.AND.JJS.LE.NBASY.AND.
     &              KKS.GE.MCASY.AND.KKS.LE.NCASY)
               IF (COND) COND=RHOMASK(IIS,JJS,KKS).GT.0
               IF (.NOT.(COND)) THEN
                  SYM=1
                  DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
                     IIS=MOD(XRSYMM(SYM,1,1)*II2+
     &                    XRSYMM(SYM,1,2)*JJ2+XRSYMM(SYM,1,3)*KK2
     &                    +(XRSYMM(SYM,1,4)*NA)/XRSYTH+10000*NA
     &                    -MAASY,NA)+MAASY
                     JJS=MOD(XRSYMM(SYM,2,1)*II2+XRSYMM(SYM,2,2)
     &                    *JJ2+XRSYMM(SYM,2,3)*KK2
     &                    +(XRSYMM(SYM,2,4)*NB)/XRSYTH+10000*NB
     &                    -MBASY,NB)+MBASY
                     KKS=MOD(XRSYMM(SYM,3,1)*II2+XRSYMM(SYM,3,2)*
     &                    JJ2+XRSYMM(SYM,3,3)*KK2
     &                    +(XRSYMM(SYM,3,4)*NC)/XRSYTH+10000*NC
     &                    -MCASY,NC)+MCASY
                     COND=(IIS.GE.MAASY.AND.IIS.LE.NAASY.AND.
     &                    JJS.GE.MBASY.AND.JJS.LE.NBASY.AND.
     &                    KKS.GE.MCASY.AND.KKS.LE.NCASY)
                     IF (COND) COND=RHOMASK(IIS,JJS,KKS).GT.0
                     SYM=SYM+1
                  END DO
               END IF
               IF (.NOT.COND) CALL WRNDIE
     &                 (-5,'XSKELeton','Internal error no. 2')
               IF (MPACK(IIS,JJS,KKS).EQ.1) THEN
                  DENSITY(IIS,JJS,KKS) = MAX(DENSITY(IIS,JJS,KKS),
     &                 DENSCUSH(II2,JJ2,KK2))
               END IF
            END DO
         END DO
      END DO
      RETURN
      END
C=======================================================================
      SUBROUTINE XSKSHELLD(I,J,K,XRINTR,DISTSQ,NA,NB,NC)
C     COMPUTES DISTANCE TO I,J,K
C Author:  Carl Elkin
      IMPLICIT NONE
C I/O
      INTEGER I,J,K,NA,NB,NC
      DOUBLE PRECISION DISTSQ,XRINTR(3,3)
C local
      DOUBLE PRECISION X,Y,Z,I2,J2,K2
C begin
C convert into orthogonal coordinates
      I2=REAL(I)/REAL(NA)
      J2=REAL(J)/REAL(NB)
      K2=REAL(K)/REAL(NC)
      X=XRINTR(1,1)*I2 +XRINTR(1,2)*J2 +XRINTR(1,3)*K2
      Y=XRINTR(2,1)*I2 +XRINTR(2,2)*J2 +XRINTR(2,3)*K2
      Z=XRINTR(3,1)*I2 +XRINTR(3,2)*J2 +XRINTR(3,3)*K2
C calculate distance
      DISTSQ = X*X+Y*Y+Z*Z
      RETURN
      END
C=======================================================================
      SUBROUTINE XSKSTATS (NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,
     &     REMOVED,NUMNODES,XNODEPOS,YNODEPOS,ZNODEPOS,DENSITY,MINDEN,
     &     OUNIT,MINOUT)
C Author:  Carl Elkin
      IMPLICIT NONE
C
      INCLUDE 'timer.inc'
C I/O
      DOUBLE PRECISION MINDEN,MINOUT
      INTEGER NUMNODES,NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,OUNIT
      REAL DENSITY (MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XNODEPOS(NUMNODES),YNODEPOS(NUMNODES),ZNODEPOS(NUMNODES)
      LOGICAL REMOVED(NUMNODES)
C local
      INTEGER I
C begin
      IF (WRNLEV.GE.10) THEN
         WRITE (6,'(A)') 'NODE #, XPOS,YPOS,ZPOS, DENSITY, REMOVED'
         DO I = 1,NUMNODES
            WRITE (6,'(I5,3I3,F5.2,L2)') I,XNODEPOS(I),YNODEPOS(I),
     &           ZNODEPOS(I),DENSITY(XNODEPOS(I),YNODEPOS(I),
     &           ZNODEPOS(I)),REMOVED(I)
         END DO
      END IF
      RETURN
      END
C=======================================================================
      SUBROUTINE XSKCOOR (NAASY,NBASY,NCASY,MAASY,MBASY,MCASY,
     &     REMOVED,NUMNODES,XNODEPOS,YNODEPOS,ZNODEPOS,DENSITY,
     &     OUNIT,XRINTR,NA,NB,NC)
C Author:  Paul Adams
      IMPLICIT NONE
C
      INCLUDE 'timer.inc'
C
      INTEGER NUMNODES,NAASY,NBASY,NCASY,MAASY,MBASY,MCASY
      INTEGER NA, NB, NC, OUNIT
      REAL DENSITY (MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XNODEPOS(NUMNODES),YNODEPOS(NUMNODES),ZNODEPOS(NUMNODES)
      LOGICAL REMOVED(NUMNODES)
      DOUBLE PRECISION XRINTR(3,3)
C local
      INTEGER I
      DOUBLE PRECISION I2, J2, K2, X, Y, Z, B, Q
C begin
      DO I = 1,NUMNODES
         IF (.NOT.REMOVED(I)) THEN
C convert into orthogonal coordinates
            I2=DFLOAT(XNODEPOS(I))/DFLOAT(NA)
            J2=DFLOAT(YNODEPOS(I))/DFLOAT(NB)
            K2=DFLOAT(ZNODEPOS(I))/DFLOAT(NC)
            X=XRINTR(1,1)*I2 +XRINTR(1,2)*J2 +XRINTR(1,3)*K2
            Y=XRINTR(2,1)*I2 +XRINTR(2,2)*J2 +XRINTR(2,3)*K2
            Z=XRINTR(3,1)*I2 +XRINTR(3,2)*J2 +XRINTR(3,3)*K2
            B=DENSITY(XNODEPOS(I),YNODEPOS(I),ZNODEPOS(I))
            Q=1.0D0
C PDB ATOM record
            WRITE(OUNIT,
     &          '(A,I5,1X,A4,1X,A4,1X,I5,3X,3F8.3,2F6.2,6X,A4)')
     &          'ATOM  ',I,'SKEL','SKEL',I,X,Y,Z,Q,B,'    '
         END IF
      END DO
C
      RETURN
      END
C=======================================================================
