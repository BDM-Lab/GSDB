      SUBROUTINE XASYMM(RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,RPNTYP,
     &                  RPNDOM,RPNLEV,DEPTH)
C
C Routine parses an ASYMm statement that defines
C the asymmetric unit for maps.  The expression
C is converted into a RPN command stack.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C
      INTEGER RPNMX, RPNX
      CHARACTER*(*) RPN(4,RPNMX)
      INTEGER RPNL(4,RPNMX), RPNN
      DOUBLE COMPLEX RPNDB(4,RPNMX)
      INTEGER RPNMLT(RPNMX), RPNLEV(RPNMX)
      CHARACTER*2 RPNTYP(RPNMX), RPNDOM(RPNMX)
      INTEGER DEPTH
C local
      INTEGER I
      CHARACTER*11 PROM
      LOGICAL ERR
C
C definition of the functions
      EXTERNAL XDOFUNC
C
C definition of the structure factor operands
      EXTERNAL XDOOPER
C
C definition of function, operand, operation types and domains
      EXTERNAL XDOTYPE
C
      CHARACTER*2 TYPE, DOMAIN
C begin
      PROM='ASYMmetric>'
C
      ERR=.FALSE.
C
C opening parenthesis
      CALL NEXTDO(PROM)
      IF (WD(1:1).EQ.'=') CALL NEXTDO(PROM)
      IF (WD(1:1).EQ.'?') THEN
      WRITE(6,'(A)')' Definition of asymmetric unit (in RPN notation).'
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
      ELSE
C
      IF (WD(1:1).NE.'(') THEN
      CALL DSPERR(PROM,'"(" expected')
      ERR=.TRUE.
      END IF
C
      IF (.NOT.ERR) THEN
C
C parse the selection (it is ok to use Hermitian symmetry here
C since we're not going to interpret structure factor or map
C operands)
      CALL EXRPN(PROM,RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,RPNTYP,
     &     RPNDOM,RPNLEV,TYPE,DOMAIN,DEPTH,XDOFUNC,XDOOPER,XDOTYPE,
     &    .TRUE.,ERR,0,' ',' ',0,' ',' ')
      CALL NEXTDO(PROM)
      IF (WD(1:1).NE.')'.AND..NOT.ERR) THEN
      CALL DSPERR(PROM,
     &  'item cannot be interpreted or ")" expected.')
      ERR=.TRUE.
      END IF
      IF (TYPE.NE.'LO'.AND..NOT.ERR) THEN
      CALL WRNDIE(-5,PROM,
     & 'Data type mismatch.  Selection must be a logical expression.')
      ERR=.TRUE.
      END IF
      END IF
C
      IF (ERR) THEN
      RPNN=0
      DEPTH=1
      END IF
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   HPRHOMA,NRHO,NMASK,XRSYGP,XRSYIV)
C
C Defines an asymmetric unit for maps.
C
C The extent of the map in P1 is defined by the
C range A=1,...,NA, B=1,...,NB, C=1,...,NC.
C The asymmetric unit is defined by the range
C A=MAASY,..,NAASY, B=MBASY,..,NBASY, C=MCASY,..,NCASY for all elements
C for which RHOMASK unequal 0.  Special positions
C are indicated by RHOMASK > 1.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'comand.inc'
      DOUBLE PRECISION MAPR
      INTEGER NA, NB, NC, NAP, NBP, NCP, NAPP, NBPP, NCPP
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9)
      INTEGER ADEPTH, ARPNMX, ARPNN, ARPNX
      CHARACTER*(*) ARPN(4,*)
      INTEGER ARPNL(4,*)
      DOUBLE COMPLEX ARPNDB(4,*)
      INTEGER ARPNMLT(*)
      INTEGER HPRHOMA, NRHO, NMASK
      INTEGER XRSYGP(XRMSYM,XRMSYM),XRSYIV(XRMSYM)
C local
      LOGICAL ERR
C variable stack pointers
      INTEGER VLEVEL, VSTACK, TSTACK
C other pointers
      INTEGER TEMPA, TEMPB, TEMPC, TAS, TBS, TCS
C min/max functions
      INTEGER NASYMM, MASYMM, KASYMM
      PARAMETER (MASYMM=6, KASYMM=5)
      DOUBLE PRECISION ASYMM(MASYMM,KASYMM,4)
      INTEGER IASYMM(MASYMM)
C timer
      DOUBLE PRECISION SECS, SECS2
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
C begin
      IF (MAPR.EQ.ZERO) THEN
      CALL WRNDIE(-5,'XASUDEF','MAPResolution not specified.')
      ELSE
      ERR=.FALSE.
C
C
C define the grid for the maps of the crystal at
C the specified resolution (MAPR)
      CALL XMAPGRD(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,
     &             XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &             XRCELL)
C
C check symmetry operators
      CALL XSYMCHK(XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,ERR,
     &             XRSYGP,XRSYIV)
C
      IF (.NOT.ERR) THEN
C
      IF (ARPNN.EQ.0) THEN
      MAASY=0
      NAASY=NA-1
      MBASY=0
      NBASY=NB-1
      MCASY=0
      NCASY=NC-1
      NASYMM=0
C
      IF (XRNSYM.GT.1) THEN
      WRITE(6,'(A)')
     & ' %XASUDEF-warning: no asymmetric unit for maps specified.'
      WRITE(6,'(A)')
     & ' %XASUDEF-warning: excessive memory and CPU usage may occur.'
      END IF
C
      WRITE(6,'(A)') ' Maps will be stored in P1:'
      WRITE(6,'(6(A,I6))') '   A=',MAASY,',...,',NAASY,
     &      '  B=',MBASY,',...,',NBASY,'  C=',MCASY,',...,',NCASY
C
      ELSE
C
C first pass: get brick size and other expressions.
      VSTACK=ALLHP(IREAL8(ADEPTH*4))
      TSTACK=ALLHP(INTEG4(ADEPTH))
      CALL XASUEVA(ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &             MAASY,NAASY,MBASY,NBASY,MCASY,NCASY,
     &             NA,NB,NC,VLEVEL,ADEPTH,HEAP(VSTACK),HEAP(TSTACK),
     &             NASYMM,MASYMM,KASYMM,ASYMM,IASYMM,ERR)
      CALL FREHP(TSTACK,INTEG4(ADEPTH))
      CALL FREHP(VSTACK,IREAL8(ADEPTH*4))
C
      WRITE(6,'(A)') ' Minimum brick that covers asymmetric unit:'
      WRITE(6,'(6(A,I6))') '   A=',MAASY,',...,',NAASY,
     &      '  B=',MBASY,',...,',NBASY,'  C=',MCASY,',...,',NCASY
      END IF
      END IF
C
C
      IF (.NOT.ERR) THEN
C
C allocate space RHOMASK: number of points in brick
      NRHO=(NAASY-MAASY+1)*(NBASY-MBASY+1)*(NCASY-MCASY+1)
      HPRHOMA=ALLHP(INTEG4(NRHO))
C
C second pass: use the min/max functions
      IF (TIMER.GT.0) CALL VCPU(SECS)
C initialize the RHOMASK array to 1
      CALL FILL4(HEAP(HPRHOMA),NRHO,1)
      CALL XASUMMM(MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &             HEAP(HPRHOMA),NASYMM,MASYMM,KASYMM,ASYMM,IASYMM)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4)')
     & ' XASUDEF: CPU-time:  fill RHOMASK=',SECS2-SECS
      END IF
C
C The asymmetric unit as specified in the International
C tables may contain redundant points at the boundary
C planes.  Furthermore, we have to mark special positions.
      IF (TIMER.GT.0) CALL VCPU(SECS)
      TEMPA=ALLHP(INTEG4(NAASY-MAASY+1))
      TEMPB=ALLHP(INTEG4(NAASY-MAASY+1))
      TEMPC=ALLHP(INTEG4(NAASY-MAASY+1))
      TAS=ALLHP(INTEG4(4*NA+NAASY-MAASY+1))
      TBS=ALLHP(INTEG4(4*NB+NBASY-MBASY+1))
      TCS=ALLHP(INTEG4(4*NC+NCASY-MCASY+1))
      CALL XASUCHK(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &       HEAP(HPRHOMA),NMASK,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &       HEAP(TEMPA),HEAP(TEMPB),HEAP(TEMPC),
     &       HEAP(TAS),HEAP(TBS),HEAP(TCS))
      CALL FREHP(TCS,INTEG4(4*NC+NCASY-MCASY+1))
      CALL FREHP(TBS,INTEG4(4*NB+NBASY-MBASY+1))
      CALL FREHP(TAS,INTEG4(4*NA+NAASY-MAASY+1))
      CALL FREHP(TEMPC,INTEG4(NAASY-MAASY+1))
      CALL FREHP(TEMPB,INTEG4(NAASY-MAASY+1))
      CALL FREHP(TEMPA,INTEG4(NAASY-MAASY+1))
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4)')
     & ' XASUDEF: CPU-time:  boundary, special pos. check=',SECS2-SECS
      END IF
      END IF
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMAPGRD(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL)
C Defines grid and sublattice at the specified map resolution.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'timer.inc'
      INCLUDE 'xfft.inc'
      DOUBLE PRECISION MAPR
      INTEGER NA, NB, NC, NAP, NBP, NCP, NAPP, NBPP, NCPP
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9)
C local
      DOUBLE PRECISION DBPREC
      DOUBLE COMPLEX DBCOMP
      LOGICAL DONE
C begin
C
C Determine the extent of a density map in P1 at the
C specified resolution (MAPR) and grid size parameter (FFT GRID...).
C The sublattice size is also determined.
C
      DONE=.FALSE.
C
      DO WHILE (.NOT.DONE)
      IF (QMEMAUTO.AND.MEMORY.LT.DEFMEMA) THEN
      MEMORY=DEFMEMA
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I10)')
     & ' %XMAPASU-AUTOmem: increasing memory allocation to ',MEMORY
      END IF
      END IF
C
      CALL XFGRID2(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,
     &             XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &             XRCELL)
C
      IF (QMEMAUTO) THEN
      IF (NAPP*NBPP*NCPP.GT.6) THEN
      MEMORY=INT(MEMORY*NAPP*NBPP*NCPP/6.0D0)
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I10)')
     & ' %XMAPASU-AUTOmem: increasing memory allocation to ',MEMORY
      END IF
      ELSE
      DONE=.TRUE.
      END IF
      ELSE
      DONE=.TRUE.
      END IF
      END DO
C
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(6(A,I4),A)')
     & ' XMAPASU: using grid [',NA,',',NB,',',NC,'] and sublattice [',
     & NAP,',',NBP,',',NCP,']'
      END IF
C
      DBPREC=NA
      CALL DECLAR( 'NA','DP',' ',DBCOMP,DBPREC)
      DBPREC=NB
      CALL DECLAR( 'NB','DP',' ',DBCOMP,DBPREC)
      DBPREC=NC
      CALL DECLAR( 'NC','DP',' ',DBCOMP,DBPREC)
C
      IF (NAPP*NBPP*NCPP.GT.12.AND.WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     &  ' *****************************************************'
      CALL WRNDIE(+1,'XFFT','Warning: Extreme factorization of FFT.')
      WRITE(6,'(A,I6,/,A)')
     &  ' %XMAPASU: Factorization of FFT=',NAPP*NBPP*NCPP,
     &  ' %XMAPASU: This may cause excessive CPU time usage.'
      WRITE(6,'(A,I12,A)')
     &  ' %XMAPASU: Recommendation: MEMOry=',
     &    INT(MEMORY*NAPP*NBPP*NCPP/12.D0),' or larger.'
      WRITE(6,'(A)')
     &  ' *****************************************************'
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XFGRID2(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL)
C
C determine integers NAPP, NAP, NBPP, NBP, NCPP, NCP such that
C 1) NA > a/(grid*MAPR)
C    NB > b/(grid*MAPR)
C    NC > c/(grid*MAPR)
C 2) NA, NB, NC as small as possible
C 3) NAP*NAPP=NA, NBP*NBPP=NB, NCP*NCPP=NC,
C 4) NAP*NBP*NCP<MEMORY,
C 5) NAP, NBP, NCP as large as possible,
C 6) space group symmetry operations conserve the grid,
C 7) NCP is an even integer for hermitian symmetry ,
C 8) NAP, NBP, NCP have only prime factors allowed by FFT
C    (see parameters set in machine dependent routine FFTPRP).
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'xfft.inc'
      INCLUDE 'consta.inc'
      DOUBLE PRECISION MAPR
      INTEGER NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9)
C
      INTEGER   ILCM
      EXTERNAL  ILCM
C local
      LOGICAL DONE
      INTEGER BASEA, BASEB, BASEC, SYM
      INTEGER NCC, NBB, NAA
      LOGICAL EQBC, EQAB, EQAC
C begin
C
C Check which cell dimensions are equal
      EQBC = .FALSE.
      EQAB = .FALSE.
      EQAC = .FALSE.
      IF (ABS(XRCELL(2)-XRCELL(3)).LT.R4SMAL) EQBC = .TRUE.
      IF (ABS(XRCELL(1)-XRCELL(2)).LT.R4SMAL) EQAB = .TRUE.
      IF (ABS(XRCELL(1)-XRCELL(3)).LT.R4SMAL) EQAC = .TRUE.
C
C determine the factors required by space-group
C symmetry.
      BASEA=1
      BASEB=1
      BASEC=1
      DO SYM=2,XRNSYM
      IF (XRSYMM(SYM,1,4).NE.0) THEN
      BASEA=MAX(BASEA,XRSYTH/XRSYMM(SYM,1,4))
      END IF
      IF (XRSYMM(SYM,2,4).NE.0) THEN
      BASEB=MAX(BASEB,XRSYTH/XRSYMM(SYM,2,4))
      END IF
      IF (XRSYMM(SYM,3,4).NE.0) THEN
      BASEC=MAX(BASEC,XRSYTH/XRSYMM(SYM,3,4))
      END IF
      END DO
      IF (XGRIDF .GT. 1) BASEA = ILCM(BASEA, XGRIDF)
      IF (YGRIDF .GT. 1) BASEB = ILCM(BASEB, YGRIDF)
      IF (ZGRIDF .GT. 1) BASEC = ILCM(BASEC, ZGRIDF)
C
C Take care of dependent cell dimensions
      IF (EQBC) THEN
        BASEB = ILCM(BASEB, BASEC)
        BASEC = BASEB
      END IF
      IF (EQAB) THEN
        BASEA = ILCM(BASEA, BASEB)
        BASEB = BASEA
      END IF
      IF (EQAC) THEN
        BASEA = ILCM(BASEA, BASEC)
        BASEC = BASEA
      END IF
C
C Check if required base factors and FTPRIM are consistent
      CALL VFYAGF(BASEA, BASEB, BASEC, BASEA, BASEB, BASEC, FTPRIM)
C
C make initial guess
      NA=MAX(2,INT(XRCELL(1)/(MAPR*GRID)+RSMALL))
      NB=MAX(2,INT(XRCELL(2)/(MAPR*GRID)+RSMALL))
      NC=MAX(2,INT(XRCELL(3)/(MAPR*GRID)+RSMALL))
      NAP=NA
      NBP=NB
      NCP=NC
      NAPP=1
      NBPP=1
      NCPP=1
      CALL XPRIME2(NAP,BASEA,NAPP,FTPRIM)
      CALL XPRIME2(NBP,BASEB,NBPP,FTPRIM)
      IF (QHERM) THEN
      CALL XPRIME2(NCP,BASEC,2*NCPP,FTPRIM)
      ELSE
      CALL XPRIME2(NCP,BASEC,NCPP,FTPRIM)
      END IF
C
C make number of grid points equal if cell dimensions are equal
      IF (EQBC) NBP=NCP
      IF (EQAB) NAP=NBP
      IF (EQAC) NAP=NCP
C
C loop until everything fits into memory
      DONE=.FALSE.
      DO WHILE (NAP*NBP*NCP.GT.MEMORY.AND..NOT.DONE)
      IF (NCP.GT.4) THEN
      NCPP=NCPP+1
      NCC=NC
      IF (QHERM) THEN
      CALL XPRIME2(NCC,BASEC,2*NCPP,FTPRIM)
      ELSE
      CALL XPRIME2(NCC,BASEC,NCPP,FTPRIM)
      END IF
      NCP=NCC/NCPP
C
      ELSE IF (NBP.GT.4) THEN
      NBPP=NBPP+1
      NBB=NB
      CALL XPRIME2(NBB,BASEB,NBPP,FTPRIM)
      NBP=NBB/NBPP
C
      ELSE IF (NAP.GT.4) THEN
      NAPP=NAPP+1
      NAA=NA
      CALL XPRIME2(NAA,BASEA,NAPP,FTPRIM)
      NAP=NAA/NAPP
C
      ELSE
      DONE=.TRUE.
      END IF
C
C make number of grid points equal if cell dimensions are equal
      IF (EQBC) NBP=NCP
      IF (EQAB) NAP=NBP
      IF (EQAC) NAP=NCP
      IF (EQBC) NBPP=NCPP
      IF (EQAB) NAPP=NBPP
      IF (EQAC) NAPP=NCPP
      END DO
C
C compute actual NA, NB, NC
      NA=NAP*NAPP
      NB=NBP*NBPP
      NC=NCP*NCPP
C
      RETURN
      END
C======================================================================
      SUBROUTINE XPRIME2(N,BASE1,BASE2,PRIME)
C
C Routine increments integer N until BASE1 is a factor of N
C and BASE2 is a factor of N and N/BASE2 does not contain any
C primes greater than PRIME.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INTEGER N, BASE1, BASE2, PRIME
C local
      INTEGER NN, P
      LOGICAL CLOOP
C begin
      NN=0
      N=N-1
      DO WHILE (NN.NE.1)
C
C increment N until base1 and base2 are factors of N
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      N=N+1
      IF (MOD(N,BASE1).EQ.0.AND.MOD(N,BASE2).EQ.0) CLOOP=.FALSE.
      END DO
C
C divide N/BASE2 by integers less equal PRIME
      NN=N/BASE2
      P=MIN(NN,PRIME)
      DO WHILE (P.GT.1)
      DO WHILE (MOD(NN,P).EQ.0)
      NN=NN/P
      END DO
      P=P-1
      END DO
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XASUMMM(MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   RHOMASK,NASYMM,MASYMM,KASYMM,ASYMM,IASYMM)
C
C uses the min/max functions in order to define RHOMASK
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C I/O
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER NASYMM, MASYMM, KASYMM
      DOUBLE PRECISION ASYMM(MASYMM,KASYMM,4)
      INTEGER IASYMM(MASYMM)
C local
      INTEGER A, B, C, I
C parameters
      INTEGER SCONS, SX, SY, SZ
      PARAMETER (SCONS=1, SX=2, SY=3, SZ=4)
C begin
      DO I=1,NASYMM
C
      IF (IASYMM(I).EQ.2) THEN
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (
     &ASYMM(I,1,SCONS)+ASYMM(I,1,SX)*A+ASYMM(I,1,SY)*B+ASYMM(I,1,SZ)*C
     &    .GT.
     &ASYMM(I,2,SCONS)+ASYMM(I,2,SX)*A+ASYMM(I,2,SY)*B+ASYMM(I,2,SZ)*C
     &  ) THEN
      RHOMASK(A,B,C)=0
      END IF
      END DO
      END DO
      END DO
C
      ELSE IF (IASYMM(I).EQ.3) THEN
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (
     &ASYMM(I,1,SCONS)+ASYMM(I,1,SX)*A+ASYMM(I,1,SY)*B+ASYMM(I,1,SZ)*C
     &    .GT. MIN(
     &ASYMM(I,2,SCONS)+ASYMM(I,2,SX)*A+ASYMM(I,2,SY)*B+ASYMM(I,2,SZ)*C,
     &ASYMM(I,3,SCONS)+ASYMM(I,3,SX)*A+ASYMM(I,3,SY)*B+ASYMM(I,3,SZ)*C)
     &  ) THEN
      RHOMASK(A,B,C)=0
      END IF
      END DO
      END DO
      END DO
C
      ELSE IF (IASYMM(I).EQ.4) THEN
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (
     &ASYMM(I,1,SCONS)+ASYMM(I,1,SX)*A+ASYMM(I,1,SY)*B+ASYMM(I,1,SZ)*C
     &    .GT. MIN(
     &ASYMM(I,2,SCONS)+ASYMM(I,2,SX)*A+ASYMM(I,2,SY)*B+ASYMM(I,2,SZ)*C,
     &ASYMM(I,3,SCONS)+ASYMM(I,3,SX)*A+ASYMM(I,3,SY)*B+ASYMM(I,3,SZ)*C,
     &ASYMM(I,4,SCONS)+ASYMM(I,4,SX)*A+ASYMM(I,4,SY)*B+ASYMM(I,4,SZ)*C)
     &  ) THEN
      RHOMASK(A,B,C)=0
      END IF
      END DO
      END DO
      END DO
C
      ELSE IF (IASYMM(I).EQ.5) THEN
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (
     &ASYMM(I,1,SCONS)+ASYMM(I,1,SX)*A+ASYMM(I,1,SY)*B+ASYMM(I,1,SZ)*C
     &    .GT. MIN(
     &ASYMM(I,2,SCONS)+ASYMM(I,2,SX)*A+ASYMM(I,2,SY)*B+ASYMM(I,2,SZ)*C,
     &ASYMM(I,3,SCONS)+ASYMM(I,3,SX)*A+ASYMM(I,3,SY)*B+ASYMM(I,3,SZ)*C,
     &ASYMM(I,4,SCONS)+ASYMM(I,4,SX)*A+ASYMM(I,4,SY)*B+ASYMM(I,4,SZ)*C,
     &ASYMM(I,5,SCONS)+ASYMM(I,5,SX)*A+ASYMM(I,5,SY)*B+ASYMM(I,5,SZ)*C)
     &  ) THEN
      RHOMASK(A,B,C)=0
      END IF
      END DO
      END DO
      END DO
      END IF
C
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XASUCHK(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                RHOMASK,NMASK,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &                TEMPA,TEMPB,TEMPC,TAS,TBS,TCS)
C
C Removes symmetry-redundant grid points and marks special
C positions.  Checks consistency between ASU and symmetry
C operators.
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C I/O
      INCLUDE 'timer.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER NMASK, XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER TEMPA(MAASY:NAASY), TEMPB(MAASY:NAASY)
      INTEGER TEMPC(MAASY:NAASY)
      INTEGER TAS(-2*NA+MAASY:2*NA+NAASY)
      INTEGER TBS(-2*NB+MBASY:2*NB+NBASY)
      INTEGER TCS(-2*NC+MCASY:2*NC+NCASY)
C local
      INTEGER SYM, A, B, C, AS, BS, CS, NSPEC, CHCKSUM, NDEL
      INTEGER AA1, AA2, AA3, BB1, BB2, BB3
C begin
C
C map all points into a primary cell
      DO C=MCASY+NC,NCASY
      CS=MOD(C+10000*NC-MCASY,NC) +MCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).EQ.1) THEN
      RHOMASK(A,B,CS)=1
      RHOMASK(A,B,C)=0
      END IF
      END DO
      END DO
      END DO
C
      DO A=MAASY+NA,NAASY
      AS=MOD(A+10000*NA-MAASY,NA) +MAASY
      DO B=MBASY,NBASY
      DO C=MCASY,NCASY
      IF (RHOMASK(A,B,C).EQ.1) THEN
      RHOMASK(AS,B,C)=1
      RHOMASK(A,B,C)=0
      END IF
      END DO
      END DO
      END DO
C
      DO B=MBASY+NB,NBASY
      BS=MOD(B+10000*NB-MBASY,NB) +MBASY
      DO A=MAASY,NAASY
      DO C=MCASY,NCASY
      IF (RHOMASK(A,B,C).EQ.1) THEN
      RHOMASK(A,BS,C)=1
      RHOMASK(A,B,C)=0
      END IF
      END DO
      END DO
      END DO
C
      NDEL=0
C
C loop over all symmetry operators except identity (first
C symmetry operator is identity).
      DO SYM=2,XRNSYM
C
C check consistency of the grid with the space group symmetry
      IF  (MOD(XRSYMM(SYM,1,4)*NA,XRSYTH).NE.0
     & .OR.MOD(XRSYMM(SYM,2,4)*NB,XRSYTH).NE.0
     & .OR.MOD(XRSYMM(SYM,3,4)*NC,XRSYTH).NE.0) THEN
      CALL WRNDIE(-5,'XMDOAS3',
     &      'Map grid inconsistent with space group.')
      END IF
C
C
C precompute stuff for the innermost loop
      DO A=MAASY,NAASY
      TEMPA(A)=XRSYMM(SYM,1,1)*A
      TEMPB(A)=XRSYMM(SYM,2,1)*A
      TEMPC(A)=XRSYMM(SYM,3,1)*A
      END DO
C
      DO A=-2*NA+MAASY,2*NA+NAASY
      TAS(A)=MOD(A+10000*NA-MAASY,NA) +MAASY
      END DO
      DO B=-2*NB+MBASY,2*NB+NBASY
      TBS(B)=MOD(B+10000*NB-MBASY,NB) +MBASY
      END DO
      DO C=-2*NC+MCASY,2*NC+NCASY
      TCS(C)=MOD(C+10000*NC-MCASY,NC) +MCASY
      END DO
C
C loop over the asymmetric unit
      DO C=MCASY,NCASY
      AA1=XRSYMM(SYM,1,3)*C+(XRSYMM(SYM,1,4)*NA)/XRSYTH
      AA2=XRSYMM(SYM,2,3)*C+(XRSYMM(SYM,2,4)*NB)/XRSYTH
      AA3=XRSYMM(SYM,3,3)*C+(XRSYMM(SYM,3,4)*NC)/XRSYTH
      DO B=MBASY,NBASY
      BB1=AA1+XRSYMM(SYM,1,2)*B
      BB2=AA2+XRSYMM(SYM,2,2)*B
      BB3=AA3+XRSYMM(SYM,3,2)*B
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).NE.0) THEN
C
C compute symmetry-related indices and project into ASU brick
      AS=TAS(BB1 + TEMPA(A))
      BS=TBS(BB2 + TEMPB(A))
      CS=TCS(BB3 + TEMPC(A))
C
      IF (AS.GT.NAASY.OR.BS.GT.NBASY.OR.CS.GT.NCASY) THEN
C outside ASU brick range
      CONTINUE
      ELSEIF (A.EQ.AS.AND.B.EQ.BS.AND.C.EQ.CS) THEN
C this is a special position unless the symmetry operator is identity
      IF (SYM.NE.1) RHOMASK(A,B,C)=RHOMASK(A,B,C)+1
      ELSEIF (RHOMASK(AS,BS,CS).NE.0) THEN
C this grid point is already present
      NDEL=NDEL+1
      RHOMASK(AS,BS,CS)=0
C
      END IF
C
      END IF
      END DO
      END DO
      END DO
C
      END DO
C
C get the number of non-zero points, special positions in RHOMASK,
C check consistency of ASU and space group, size for smaller
C brick
      NMASK=0
      NSPEC=0
      CHCKSUM=0
      DO A=MAASY,NAASY
      DO B=MBASY,NBASY
      DO C=MCASY,NCASY
      IF (RHOMASK(A,B,C).EQ.1) THEN
      NMASK=NMASK+1
      CHCKSUM=CHCKSUM+XRNSYM
      ELSEIF (RHOMASK(A,B,C).GT.1) THEN
      NMASK=NMASK+1
      NSPEC=NSPEC+1
      CHCKSUM=CHCKSUM+XRNSYM/RHOMASK(A,B,C)
      END IF
      END DO
      END DO
      END DO
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A,I12)') ' Number of deleted points in ASU=',
     &  NDEL
      WRITE(6,'(A,I12)') ' Number of non-zero lattice points in ASU=',
     &  NMASK
      WRITE(6,'(A,I12)') ' Number of special positions in ASU=',
     &  NSPEC
      END IF
C
      IF (CHCKSUM.NE.NA*NB*NC) THEN
      WRITE(6,'(A,I12)') ' XMDOAS3: expected number of points',NA*NB*NC
      WRITE(6,'(A,I12)') ' XMDOAS3: actual number of points',CHCKSUM
      CALL WRNDIE(-5,'XMDOAS3',
     & 'Asymmetric map unit is incompatible with symmetry operators.')
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XASUEVA(RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,
     &                   MAASY,NAASY,MBASY,NBASY,MCASY,NCASY,
     &                   NA,NB,NC,VLEVEL,VMAX,VSTACK,TSTACK,
     &                   NASYMM,MASYMM,KASYMM,ASYMM,IASYMM,ERR)
C
C Evaluates an asymmetric unit expression from the RPN
C command stack.
C
C Expressions that are allowed cover all asymmetric unit
C specifications in the International Tables:
C
C [   max(<expr>,...,<expr>) | <expr    <rel> ]
C       <expr>
C              [ <rel>     <expr> |  min(<expr>,...,<expr>) ]
C
C <rel> :== < | <=
C
C <expr> :== -<expr> | <co> | <arg> | <expr>*<expr> | <expr>+<expr> |
C            <expr>-<expr> | <expr>/<expr> | <expr> <expr>
C <co> := constant
C <arg> :== x,y,z
C
C Restrictions:
C 1. all <expr> must be linear in x, y, z
C 2. max. of four arguments for min, max functions.
C
C
C The brick determined by these a expresssions is determined
C (maasy, naasy; mbasy, nbasy; mcasy, ncasy).
C
C A data structure is defined that represents the min, max
C functions.
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C input/output
      INCLUDE 'consta.inc'
      INTEGER RPNMX, RPNN, RPNX
      CHARACTER*(*) RPN(4,*)
      INTEGER RPNL(4,*)
      DOUBLE COMPLEX RPNDB(4,*)
      INTEGER RPNMLT(*)
      INTEGER MAASY, NAASY, MBASY, NBASY, MCASY, NCASY
      INTEGER NA, NB, NC
      INTEGER VLEVEL, VMAX
      DOUBLE PRECISION VSTACK(VMAX,4)
      INTEGER TSTACK(*)
      INTEGER NASYMM, MASYMM, KASYMM
      DOUBLE PRECISION ASYMM(MASYMM,KASYMM,4)
      INTEGER IASYMM(MASYMM)
      LOGICAL ERR
C local
      INTEGER RPNI, TOL, I
C parameters
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      INTEGER SCONS, SX, SY, SZ, SMIX, SMIN, SMAX
      PARAMETER (SCONS=1, SX=2, SY=3, SZ=4, SMIX=5, SMIN=10, SMAX=20)
C begin
C
C initialize variable stack level
      VLEVEL=0
C
C initialize min max function data structure
      NASYMM=0
C
C reset error
      ERR=.FALSE.
C
C initialize brick size
      MAASY=-NA
      NAASY=+NA
      MBASY=-NB
      NBASY=+NB
      MCASY=-NC
      NCASY=+NC
C
      DO RPNI=1,RPNN
C contants
      IF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'CONS') THEN
      VLEVEL=VLEVEL+1
      VSTACK(VLEVEL,SCONS)=DBLE(RPNDB(1,RPNI))
      VSTACK(VLEVEL,SX)=ZERO
      VSTACK(VLEVEL,SY)=ZERO
      VSTACK(VLEVEL,SZ)=ZERO
      TSTACK(VLEVEL)=SCONS
C change sign
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'CHS'.AND.
     &                     TSTACK(VLEVEL).LE.SMIX) THEN
      VSTACK(VLEVEL,SCONS)=-VSTACK(VLEVEL,SCONS)
      VSTACK(VLEVEL,SX)=-VSTACK(VLEVEL,SX)
      VSTACK(VLEVEL,SY)=-VSTACK(VLEVEL,SY)
      VSTACK(VLEVEL,SZ)=-VSTACK(VLEVEL,SZ)
C plus
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'+'.AND.
     &    TSTACK(VLEVEL).LE.SMIX.AND.TSTACK(VLEVEL-1).LE.SMIX) THEN
      VLEVEL=VLEVEL-1
      VSTACK(VLEVEL,SCONS)=VSTACK(VLEVEL,SCONS)+VSTACK(VLEVEL+1,SCONS)
      VSTACK(VLEVEL,SX)=VSTACK(VLEVEL,SX)+VSTACK(VLEVEL+1,SX)
      VSTACK(VLEVEL,SY)=VSTACK(VLEVEL,SY)+VSTACK(VLEVEL+1,SY)
      VSTACK(VLEVEL,SZ)=VSTACK(VLEVEL,SZ)+VSTACK(VLEVEL+1,SZ)
      IF (TSTACK(VLEVEL).NE.TSTACK(VLEVEL+1)) THEN
      TSTACK(VLEVEL)=SMIX
      END IF
C minus
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'-'.AND.
     &    TSTACK(VLEVEL).LE.SMIX.AND.TSTACK(VLEVEL-1).LE.SMIX) THEN
      VLEVEL=VLEVEL-1
      VSTACK(VLEVEL,SCONS)=VSTACK(VLEVEL,SCONS)-VSTACK(VLEVEL+1,SCONS)
      VSTACK(VLEVEL,SX)=VSTACK(VLEVEL,SX)-VSTACK(VLEVEL+1,SX)
      VSTACK(VLEVEL,SY)=VSTACK(VLEVEL,SY)-VSTACK(VLEVEL+1,SY)
      VSTACK(VLEVEL,SZ)=VSTACK(VLEVEL,SZ)-VSTACK(VLEVEL+1,SZ)
      IF (TSTACK(VLEVEL).NE.TSTACK(VLEVEL+1)) THEN
      TSTACK(VLEVEL)=SMIX
      END IF
C multiply: arg * constant
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'*'
     &        .AND.TSTACK(VLEVEL).EQ.SCONS) THEN
      VLEVEL=VLEVEL-1
      VSTACK(VLEVEL,SCONS)=VSTACK(VLEVEL,SCONS)*VSTACK(VLEVEL+1,SCONS)
      VSTACK(VLEVEL,SX)=VSTACK(VLEVEL,SX)*VSTACK(VLEVEL+1,SCONS)
      VSTACK(VLEVEL,SY)=VSTACK(VLEVEL,SY)*VSTACK(VLEVEL+1,SCONS)
      VSTACK(VLEVEL,SZ)=VSTACK(VLEVEL,SZ)*VSTACK(VLEVEL+1,SCONS)
C multiply: constant * arg
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'*'
     &        .AND.TSTACK(VLEVEL-1).EQ.SCONS) THEN
      VLEVEL=VLEVEL-1
      VSTACK(VLEVEL,SX)=VSTACK(VLEVEL+1,SX)*VSTACK(VLEVEL,SCONS)
      VSTACK(VLEVEL,SY)=VSTACK(VLEVEL+1,SY)*VSTACK(VLEVEL,SCONS)
      VSTACK(VLEVEL,SZ)=VSTACK(VLEVEL+1,SZ)*VSTACK(VLEVEL,SCONS)
      VSTACK(VLEVEL,SCONS)=VSTACK(VLEVEL+1,SCONS)*VSTACK(VLEVEL,SCONS)
      TSTACK(VLEVEL)=TSTACK(VLEVEL+1)
C divide
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'/'
     &       .AND.TSTACK(VLEVEL).EQ.SCONS) THEN
      VLEVEL=VLEVEL-1
      IF (ABS(VSTACK(VLEVEL+1,SCONS)).LT.RSMALL) THEN
      CALL WRNDIE(-5,'XASEVAL',
     &  'Divide by zero in asymmetric unit expression.')
      ERR=.TRUE.
      ELSE
      VSTACK(VLEVEL,SCONS)=VSTACK(VLEVEL,SCONS)/VSTACK(VLEVEL+1,SCONS)
      VSTACK(VLEVEL,SX)=VSTACK(VLEVEL,SX)/VSTACK(VLEVEL+1,SCONS)
      VSTACK(VLEVEL,SY)=VSTACK(VLEVEL,SY)/VSTACK(VLEVEL+1,SCONS)
      VSTACK(VLEVEL,SZ)=VSTACK(VLEVEL,SZ)/VSTACK(VLEVEL+1,SCONS)
      END IF
C multiple-argument functions
C
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'MAX') THEN
      IF (NASYMM.GE.MASYMM) THEN
      CALL WRNDIE(-5,'XASEVAL',
     &   'MASYMM exceeded.  Too many expressions for asymmetric unit.')
      ERR=.TRUE.
      ELSEIF (KASYMM.LE.RPNMLT(RPNI)) THEN
      CALL WRNDIE(-5,'XASEVAL',
     &   'MASYMM exceeded.  Too many arguments for MAX function.')
      ERR=.TRUE.
      ELSE
      NASYMM=NASYMM+1
      VLEVEL=VLEVEL-RPNMLT(RPNI)+1
      TSTACK(VLEVEL)=SMAX
      DO I=1,RPNMLT(RPNI)
C
C negate the expression, e.g., max(y)<x <=> -x<min(-y)
      ASYMM(NASYMM,I+1,SCONS)=-VSTACK(VLEVEL+I-1,SCONS)
      ASYMM(NASYMM,I+1,SX)=-VSTACK(VLEVEL+I-1,SX)
      ASYMM(NASYMM,I+1,SY)=-VSTACK(VLEVEL+I-1,SY)
      ASYMM(NASYMM,I+1,SZ)=-VSTACK(VLEVEL+I-1,SZ)
      END DO
      IASYMM(NASYMM)=RPNMLT(RPNI)+1
      END IF
C
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'MIN') THEN
      IF (NASYMM.GE.MASYMM) THEN
      CALL WRNDIE(-5,'XASEVAL',
     &   'MASYMM exceeded.  Too many expressions for asymmetric unit.')
      ERR=.TRUE.
      ELSEIF (KASYMM.LE.RPNMLT(RPNI)) THEN
      CALL WRNDIE(-5,'XASEVAL',
     &   'MASYMM exceeded.  Too many arguments for MIN function.')
      ERR=.TRUE.
      ELSE
      NASYMM=NASYMM+1
      VLEVEL=VLEVEL-RPNMLT(RPNI)+1
      TSTACK(VLEVEL)=SMIN
      DO I=1,RPNMLT(RPNI)
      ASYMM(NASYMM,I+1,SCONS)=VSTACK(VLEVEL+I-1,SCONS)
      ASYMM(NASYMM,I+1,SX)=VSTACK(VLEVEL+I-1,SX)
      ASYMM(NASYMM,I+1,SY)=VSTACK(VLEVEL+I-1,SY)
      ASYMM(NASYMM,I+1,SZ)=VSTACK(VLEVEL+I-1,SZ)
      END DO
      IASYMM(NASYMM)=RPNMLT(RPNI)+1
      END IF
C
C logical operations
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'AND') THEN
      VLEVEL=VLEVEL-1
C relational operators
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'<'.OR.
     &        RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'<=') THEN
      IF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'<') THEN
      TOL=1
      ELSE
      TOL=0
      END IF
C
C check if this is a relation of the type x<constant or constant<x
      IF ((TSTACK(VLEVEL).LE.SZ.AND.TSTACK(VLEVEL-1).EQ.SCONS).OR.
     &    (TSTACK(VLEVEL-1).LE.SZ.AND.TSTACK(VLEVEL).EQ.SCONS)) THEN
      IF (TSTACK(VLEVEL).EQ.SCONS) THEN
      IF (TSTACK(VLEVEL-1).EQ.SX) THEN
C x<constant
      NAASY=MIN(NAASY,INT(VSTACK(VLEVEL,SCONS)
     &              /VSTACK(VLEVEL-1,SX)+10000+R4SMAL)-TOL-10000)
      ELSEIF (TSTACK(VLEVEL-1).EQ.SY) THEN
C y<constant
      NBASY=MIN(NBASY,INT(VSTACK(VLEVEL,SCONS)
     &              /VSTACK(VLEVEL-1,SY)+10000+R4SMAL)-TOL-10000)
      ELSEIF (TSTACK(VLEVEL-1).EQ.SZ) THEN
C z<constant
      NCASY=MIN(NCASY,INT(VSTACK(VLEVEL,SCONS)
     &              /VSTACK(VLEVEL-1,SZ)+10000+R4SMAL)-TOL-10000)
      ELSE
      ERR=.TRUE.
      END IF
      ELSEIF (TSTACK(VLEVEL-1).EQ.SCONS) THEN
      IF (TSTACK(VLEVEL).EQ.SX) THEN
C constant<x
      MAASY=MAX(MAASY,INT(VSTACK(VLEVEL-1,SCONS)
     &              /VSTACK(VLEVEL,SX)+10000+R4SMAL)+TOL-10000)
      ELSEIF (TSTACK(VLEVEL).EQ.SY) THEN
C constant<y
      MBASY=MAX(MBASY,INT(VSTACK(VLEVEL-1,SCONS)
     &              /VSTACK(VLEVEL,SY)+10000+R4SMAL)+TOL-10000)
      ELSEIF (TSTACK(VLEVEL).EQ.SZ) THEN
C constant<z
      MCASY=MAX(MCASY,INT(VSTACK(VLEVEL-1,SCONS)
     &              /VSTACK(VLEVEL,SZ)+10000+R4SMAL)+TOL-10000)
      ELSE
      ERR=.TRUE.
      END IF
      ELSE
      ERR=.TRUE.
      END IF
C
      ELSEIF((TSTACK(VLEVEL).LE.SMIX.AND.TSTACK(VLEVEL-1).LE.SMIX))THEN
C
C relation of the type x<y
      IF (NASYMM.GE.MASYMM) THEN
      CALL WRNDIE(-5,'XASEVAL',
     &   'MASYMM exceeded.  Too many expressions for asymmetric unit.')
      ERR=.TRUE.
      ELSE
      NASYMM=NASYMM+1
      IF (TOL.EQ.1) THEN
C for "<"
      ASYMM(NASYMM,1,SCONS)=VSTACK(VLEVEL-1,SCONS)+RSMALL
      ELSE
C for "<="
      ASYMM(NASYMM,1,SCONS)=VSTACK(VLEVEL-1,SCONS)-RSMALL
      END IF
      ASYMM(NASYMM,1,SX)=VSTACK(VLEVEL-1,SX)
      ASYMM(NASYMM,1,SY)=VSTACK(VLEVEL-1,SY)
      ASYMM(NASYMM,1,SZ)=VSTACK(VLEVEL-1,SZ)
      ASYMM(NASYMM,2,SCONS)=VSTACK(VLEVEL,SCONS)
      ASYMM(NASYMM,2,SX)=VSTACK(VLEVEL,SX)
      ASYMM(NASYMM,2,SY)=VSTACK(VLEVEL,SY)
      ASYMM(NASYMM,2,SZ)=VSTACK(VLEVEL,SZ)
      IASYMM(NASYMM)=2
      END IF
C
      ELSEIF((TSTACK(VLEVEL).EQ.SMIN.AND.TSTACK(VLEVEL-1).LE.SMIX))THEN
      IF (NASYMM.GE.MASYMM) THEN
      CALL WRNDIE(-5,'XASEVAL',
     &   'MASYMM exceeded.  Too many expressions for asymmetric unit.')
      ERR=.TRUE.
      ELSE
C
C we don't have to increment NASYMM since this was already done when
C the MIN function was executed.
      IF (TOL.EQ.1) THEN
C for "<"
      ASYMM(NASYMM,1,SCONS)=VSTACK(VLEVEL-1,SCONS)+RSMALL
      ELSE
C for "<="
      ASYMM(NASYMM,1,SCONS)=VSTACK(VLEVEL-1,SCONS)-RSMALL
      END IF
      ASYMM(NASYMM,1,SX)=VSTACK(VLEVEL-1,SX)
      ASYMM(NASYMM,1,SY)=VSTACK(VLEVEL-1,SY)
      ASYMM(NASYMM,1,SZ)=VSTACK(VLEVEL-1,SZ)
      END IF
C
      ELSEIF((TSTACK(VLEVEL-1).EQ.SMAX.AND.TSTACK(VLEVEL).LE.SMIX))THEN
      IF (NASYMM.GE.MASYMM) THEN
      CALL WRNDIE(-5,'XASEVAL',
     &   'MASYMM exceeded.  Too many expressions for asymmetric unit.')
      ERR=.TRUE.
      ELSE
C
C we don't have to increment NASYMM since this was already done when
C the MIN function was executed.
      IF (TOL.EQ.1) THEN
C for "<"
      ASYMM(NASYMM,1,SCONS)=-VSTACK(VLEVEL,SCONS)+RSMALL
      ELSE
C for "<="
      ASYMM(NASYMM,1,SCONS)=-VSTACK(VLEVEL,SCONS)-RSMALL
      END IF
      ASYMM(NASYMM,1,SX)=-VSTACK(VLEVEL,SX)
      ASYMM(NASYMM,1,SY)=-VSTACK(VLEVEL,SY)
      ASYMM(NASYMM,1,SZ)=-VSTACK(VLEVEL,SZ)
      END IF
C
      ELSE
      ERR=.TRUE.
      END IF
      VLEVEL=VLEVEL-1
C other map-specific functions
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'X') THEN
      VLEVEL=VLEVEL+1
      VSTACK(VLEVEL,SCONS)=ZERO
      VSTACK(VLEVEL,SX)=ONE/NA
      VSTACK(VLEVEL,SY)=ZERO
      VSTACK(VLEVEL,SZ)=ZERO
      TSTACK(VLEVEL)=SX
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'Y') THEN
      VLEVEL=VLEVEL+1
      VSTACK(VLEVEL,SCONS)=ZERO
      VSTACK(VLEVEL,SX)=ZERO
      VSTACK(VLEVEL,SY)=ONE/NB
      VSTACK(VLEVEL,SZ)=ZERO
      TSTACK(VLEVEL)=SY
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'Z') THEN
      VLEVEL=VLEVEL+1
      VSTACK(VLEVEL,SCONS)=ZERO
      VSTACK(VLEVEL,SX)=ZERO
      VSTACK(VLEVEL,SY)=ZERO
      VSTACK(VLEVEL,SZ)=ONE/NC
      TSTACK(VLEVEL)=SZ
C recall
      ELSEIF (RPN(1,RPNI)(1:RPNL(1,RPNI)).EQ.'RECALL') THEN
      VLEVEL=VLEVEL+1
      ELSE
      ERR=.TRUE.
      END IF
      END DO
C
C error checking
      IF (MAASY.GE.NAASY.OR.MBASY.GE.NBASY.OR.MCASY.GE.NCASY) THEN
      CALL WRNDIE(-5,'XASEVAL',
     &  'Expression defines a zero asymmetric unit.')
      ERR=.TRUE.
      END IF
C
      IF (ERR) THEN
      CALL WRNDIE(-5,'XASEVAL',
     &   'Illegal expression for asymmetric unit.')
      END IF
C
      END
