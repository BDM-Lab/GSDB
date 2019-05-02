C
C=====================================================================
C
      SUBROUTINE NOEEXC
     &          (ICL,NOENUM,VIOEXC,NSTRUC,MODE,VIOLOL,VIOUPL,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
C plan: get frequency differences into the scheme
C 1) calculate all exp(frequency differences) and store
C    (e.g. NOEWGT = NOEWGT * ...)
C    this might be done in ariass somewhere.
C    need to use standard distance-weighted frequency of assignment
C 2) normalize such that the smallest difference gets weight 1
C 3) exclude if probability is > 0.5
C    or cumulative figure of merit < 0.5
C    P(Nviol,DelFreq) = Nviol/Nstruct  / (F(DelFreq) * InitWeight)
C
C excludes datapoints for which NOEVIO exceeds threshold
C VIOEXC*NSTRUC by setting NOECV to 0
C
C Author: Michael Nilges
C
      IMPLICIT NONE
C I/O
C
      INCLUDE 'numbers.inc'
      INCLUDE 'consta.inc'
C
      INTEGER NOENUM,ICL,NSTRUC
      DOUBLE PRECISION VIOEXC,VIOLOL,VIOUPL
      CHARACTER*4 MODE
C
C global NOE arrays on HEAP:
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
C
C local
      INTEGER N, NEXC
      DOUBLE PRECISION AWEIGH,AFIGME,FIGMER
      INTEGER NAVE
C begin
C
      IF (NOENUM.GT.0) THEN
C
      IF (MODE.EQ.'INIT') THEN
        DO N=1,NOENUM
          NOECV(N)=1
          NOEVIO(N)=0
          NOERRV(N)=ZERO
          NOERAV(N)=ZERO
        END DO
C
C else calculate a figure of merit for the peak from the initial weight,
C the frequency figure of merit, and the violation statistic.
      ELSE
        AWEIGH=ZERO
        AFIGME=ZERO
        NAVE=0
        DO N=1,NOENUM
          IF (NOECND(N).EQ.ICL) THEN
            IF (NOEWGH(N).LE.ONE) THEN
              NAVE=NAVE+1
              AWEIGH=AWEIGH+NOEWGH(N)
            END IF
          END IF
        END DO
        AWEIGH=AWEIGH/NAVE
        WRITE(6,'(A,F15.3)')
     &  ' COUNtviolations: average initial weight: ', AWEIGH
        NAVE=0
        DO N=1,NOENUM
          IF (NOECND(N).EQ.ICL) THEN
            IF (NOERAV(N).GT.ZERO) THEN
              NAVE=NAVE+1
              AFIGME=AFIGME+NOERAV(N)
            END IF
          END IF
        END DO
        AFIGME=AFIGME/MAX(NAVE,1)
        WRITE(6,'(A,F15.3)')
     &  ' COUNtviolations: average figure of merit: ', AFIGME
        NEXC=0
        DO N=1,NOENUM
          IF (NOECND(N).EQ.ICL) THEN
            FIGMER=ONE-FLOAT(NOEVIO(N))/MAX(1,NSTRUC)
            IF (NOERAV(N).GT.ZERO) FIGMER=FIGMER*NOERAV(N)/AFIGME
            NOECV(N) = 1
            IF (NOEWGH(N).LE.ONE) THEN
              FIGMER=FIGMER*NOEWGH(N)/AWEIGH
              IF (FIGMER.LE.(ONE-VIOEXC)) THEN
                NEXC=NEXC+1
                IF (MODE.EQ.'EXCL') THEN
                  NOECV(N) = 0
                ELSEIF (MODE.EQ.'SET ') THEN
C
C if TOKEN was specified on input, VIOxxL are -1
                  IF (VIOUPL.GE.ZERO) THEN
                    NOEHIG(N) = MAX(ZERO,(VIOUPL-NOEDIS(N)))
                  END IF
                  IF (VIOLOL.GE.ZERO) THEN
                    NOELOW(N) = MAX(ZERO,(NOEDIS(N)-VIOLOL))
                  END IF
                  WRITE(6,'(A,F15.3,F15.3)')
     &  ' COUNtviolations: new error estimates: ', NOELOW(N), NOEHIG(N)
                ELSEIF (MODE.EQ.'MULT') THEN
C
C if TOKEN was specified on input, VIOxxL are -1
                  IF (VIOUPL.GE.ZERO) THEN
                    NOEHIG(N) = NOEHIG(N)*VIOUPL
                  END IF
                  IF (VIOLOL.GE.ZERO) THEN
                    NOELOW(N) = NOELOW(N)*VIOLOL
                  END IF
                  WRITE(6,'(A,F15.3,F15.3)')
     &  ' COUNtviolations: new error estimates: ', NOELOW(N), NOEHIG(N)
C
                ELSEIF (MODE.EQ.'WEIG') THEN
                  NOEWGH(N)=EXP(-FIGMER**2)
                ELSE
                  WRITE(6,'(A)')
     &  ' %COUNTV-ERR: unknown mode'
                END IF
              END IF
            END IF
          END IF
        END DO
      END IF
      END IF
C
      IF (MODE.EQ.'INIT') THEN
      WRITE (6,'(A)') ' COUNtviolations: database initialized'
      ELSE IF (MODE.EQ.'EXCL') THEN
      WRITE (6,'(A,I5,A)')
     &  ' COUNtviolations: ', NEXC, ' restraints excluded '
      ELSE
      WRITE (6,'(A,I5,A)')
     &  ' COUNtviolations: ', NEXC, ' limits reset '
      END IF
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE ARIVIO
     &          (NOENUM,VIOTHR,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
C
C subroutine counts noe violations and calculates average violation
C
C Author: Michael Nilges
C
      IMPLICIT NONE
C I/O
      INCLUDE 'pick.inc'
      DOUBLE PRECISION VIOTHR
      INTEGER NOENUM
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
      DOUBLE PRECISION EN
      INTEGER N, NVIO
C
C begin
C
      NVIO=0
      IF (NOENUM.GT.0) THEN
      DO N=1,NOENUM
      CALL ENOE(EN,'ANAL',N)
      NOERRV(N)=NOERRV(N)+PCDATA(PCDEVI)*PCDATA(PCDEVI)
      IF (ABS(PCDATA(PCDEVI)).GT.VIOTHR) THEN
      NOEVIO(N)=NOEVIO(N)+1
      NVIO=NVIO+1
      END IF
      END DO
      END IF
      WRITE(6,'(A,I5,A,f15.5)')
     &     ' COUNTV:',NVIO,' violations above ', VIOTHR
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE NVIOLS
     &          (NOENUM,VIOEXC,NSTRUC,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
C lists datapoints for which NOEVIO exceeds threshold
C VIOEXC*NSTRUC
C Author: Michael Nilges
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'numbers.inc'
      INTEGER NOENUM,NSTRUC
      DOUBLE PRECISION VIOEXC
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
C
C local
      INTEGER I,J,K,N, NORR
C
C begin
C
      IF (NOENUM.GT.0) THEN
C{
      DO N=1,NOENUM
C{
      IF (NOEVIO(N).GT.VIOEXC*NSTRUC) THEN
C{
      WRITE(PUNIT,'(A,I5,A)')
     & ' ========== restraint ',NOEPID(N),' =========='

      DO NORR=NOEORR(N)+1, NOEORR(N+1)
      WRITE(PUNIT,'(A)') ' set-i-atoms'
      DO I=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      K=NOEILS(I)
      WRITE(PUNIT,'(9X,4(1X,A))') SEGID(K),RESID(K),RES(K),TYPE(K)
      END DO
      WRITE(PUNIT,'(A)') ' set-j-atoms'
      DO J=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      K=NOEJLS(J)
      WRITE(PUNIT,'(9X,4(1X,A))') SEGID(K),RESID(K),RES(K),TYPE(K)
      END DO
      END DO
C
      WRITE(PUNIT,'(I5,A,F5.2,A,F5.2,A,F5.2,A,F5.2,A)')
     &  NOEVIO(N), ' Violations above threshold, rms violation=',
     &  SQRT(NOERRV(N)/MAX(NOEVIO(N),1)),
     &  ' Target distance=',NOEDIS(N),
     &  ' (-',NOELOW(N),'/+',NOEHIG(N),')'
      END IF
C}
      END DO
C}
      END IF
C}
      RETURN
      END
C
C=====================================================================
C
