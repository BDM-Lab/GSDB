      SUBROUTINE PSFRDR
C
C This routine reads the protein structure file module (formatted file)
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'ener.inc'
C local
      INTEGER I, II, U, NGR2, DUMMY
      INTEGER NATO2, NBON2, NTHET2, NPH2, NIMPH2, NDO2, NAC2, NN2
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      ERROR=.FALSE.
      EOF=.FALSE.
C
      U=ISTRM(NSTRM)
C
C read fixed-field formatted file
      IF (.NOT.ERROR) THEN
      CALL RDTITB(U,'FORMATTED')
      READ(U,'(/I8)',ERR=7,END=6) NATO2
      IF (NATOM+NATO2.GT.MAXA) THEN
      CALL WRNDIE(-5,'PSFRDR',
     & 'exceeded MAXA parameter --> recompile program')
      END IF
      READ(U,'(I8,1X,A,1X,A,1X,A,1X,A,1X,A,1X,2G14.6,I8)',ERR=7,END=6)
     &               (II,SEGID(I),RESID(I),RES(I),TYPE(I),IAC(I),
     &               CG(I),AMASS(I),DUMMY,I=NATOM+1,NATOM+NATO2)
      DO I=NATOM+1,NATOM+NATO2
      QLOOKU(I)=.TRUE.
      END DO
      READ(U,'(/I8)',ERR=7,END=6) NBON2
      IF (NBOND+NBON2.GT.MAXB) THEN
      CALL WRNDIE(-5,'PSFRDR',
     & 'exceeded MAXB parameter --> recompile program')
      END IF
      READ(U,'(8I8)',ERR=7,END=6) (IB(I),JB(I),I=NBOND+1,NBOND+NBON2)
      DO I=NBOND+1,NBOND+NBON2
      QACBB(I)=.TRUE.
      QACBC(I)=.TRUE.
      END DO
      READ(U,'(/I8)',ERR=7,END=6) NTHET2
      IF (NTHETA+NTHET2.GT.MAXT) THEN
      CALL WRNDIE(-5,'PSFRDR',
     & 'exceeded MAXT parameter --> recompile program')
      END IF
      READ(U,'(9I8)',ERR=7,END=6) (IT(I),JT(I),KT(I),I=NTHETA+1,
     &                            NTHETA+NTHET2)
      DO I=NTHETA+1,NTHETA+NTHET2
      QACTB(I)=.TRUE.
      QACTC(I)=.TRUE.
C
C initialize Urey-Bradley terms
      QACTUB(I)=.TRUE.
      QACTUC(I)=.TRUE.
      ACTUB(I)=ZERO
      ACTUC(I)=ZERO
C
      END DO
      READ(U,'(/I8)',ERR=7,END=6) NPH2
      IF (NPHI+NPH2.GT.MAXP) THEN
      CALL WRNDIE(-5,'PSFRDR',
     & 'exceeded MAXP parameter --> recompile program')
      END IF
      READ(U,'(8I8)',ERR=7,END=6) (IP(I),JP(I),KP(I),LP(I),I=NPHI+1,
     &                            NPHI+NPH2)
      DO I=NPHI+1,NPHI+NPH2
      QACPB(I)=.TRUE.
      QACPC(I)=.TRUE.
      QACPD(I)=.TRUE.
      END DO
      READ(U,'(/I8)',ERR=7,END=6) NIMPH2
      IF (NIMPHI+NIMPH2.GT.MAXIMP) THEN
      CALL WRNDIE(-5,'PSFRDR',
     & 'exceeded NIMPHI parameter --> recompile program')
      END IF
      READ(U,'(8I8)',ERR=7,END=6)(IM(I),JM(I),KM(I),LM(I),I=NIMPHI+1,
     &                            NIMPHI+NIMPH2)
      DO I=NIMPHI+1,NIMPHI+NIMPH2
      QACIB(I)=.TRUE.
      QACIC(I)=.TRUE.
      QACID(I)=.TRUE.
      END DO
      READ(U,'(/I8)',ERR=7,END=6) NDO2
      READ(U,'(8I8)',ERR=7,END=6) (DUMMY,DUMMY,I=1,NDO2)
      READ(U,'(/I8)',ERR=7,END=6) NAC2
      READ(U,'(8I8)',ERR=7,END=6) (DUMMY,DUMMY,I=1,NAC2)
      READ(U,'(/I8)',ERR=7,END=6) NN2
      IF (NNB+NN2.GT.MAXNB) THEN
      CALL WRNDIE(-5,'PSFRDR',
     & 'exceeded MAXNB parameter --> recompile program')
      END IF
      READ(U,'(8I8)',ERR=7,END=6) (INB(I),I=NNB+1,NNB+NN2)
      READ(U,'(8I8)',ERR=7,END=6) (IBLO(I),I=NATOM+1,NATOM+NATO2)
      READ(U,'(/I8)',ERR=7,END=6) NGR2
      IF (NGRP+NGR2.GT.MAXGRP) THEN
      CALL WRNDIE(-5,'PSFRDR',
     & 'exceeded MAXGRP parameter --> recompile program')
      END IF
      READ(U,'(9I8)',ERR=7,END=6)
     &       (IGPBS(I),DUMMY,DUMMY,I=NGRP+1,NGRP+NGR2)
      IGPBS(NGRP+NGR2+1)=NATOM+NATO2
      END IF
C-error-labels----
      GOTO 77
7     ERROR=.TRUE.
      GOTO 77
6     EOF=.TRUE.
77    CONTINUE
C-----------------
C
      IF (ERROR) WRITE(6,'(A)') ' %PSFRDR-ERR: error during read'
      IF (EOF) WRITE(6,'(A)') ' %PSFRDR-ERR: EOF during read'
C
      IF (.NOT.ERROR.AND..NOT.EOF) THEN
C
      IF (NATOM.GT.0) THEN
C shift the atom indices
      DO I=NBOND+1,NBOND+NBON2
      IB(I)=IB(I)+NATOM
      JB(I)=JB(I)+NATOM
      END DO
      DO I=NTHETA+1,NTHETA+NTHET2
      IT(I)=IT(I)+NATOM
      JT(I)=JT(I)+NATOM
      KT(I)=KT(I)+NATOM
      END DO
      DO I=NPHI+1,NPHI+NPH2
      IP(I)=IP(I)+NATOM
      JP(I)=JP(I)+NATOM
      KP(I)=KP(I)+NATOM
      LP(I)=LP(I)+NATOM
      END DO
      DO I=NIMPHI+1,NIMPHI+NIMPH2
      IM(I)=IM(I)+NATOM
      JM(I)=JM(I)+NATOM
      KM(I)=KM(I)+NATOM
      LM(I)=LM(I)+NATOM
      END DO
      DO I=NNB+1,NNB+NN2
      INB(I)=INB(I)+NATOM
      END DO
      DO I=NATOM+1,NATOM+NATO2
      IBLO(I)=IBLO(I)+IBLO(NATOM)
      END DO
      DO I=NGRP+1,NGRP+NGR2
      IGPBS(I)=IGPBS(I)+NATOM
      END DO
      END IF
C
C update the max-counters
      NATOM=NATOM+NATO2
      NBOND=NBOND+NBON2
      NTHETA=NTHETA+NTHET2
      NPHI=NPHI+NPH2
      NIMPHI=NIMPHI+NIMPH2
      NNB=NNB+NN2
      NGRP=NGRP+NGR2
      END IF
      RETURN
      END
C
      SUBROUTINE PSFWRT(FILE,DUMMY)
C
C Write PSF file.
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'ener.inc'
C
      CHARACTER*(*) FILE
      INTEGER DUMMY
C local
      INTEGER I, U
      LOGICAL ERROR
C begin
C
      CALL ASSFIL(FILE,U,'WRITE','FORMATTED',ERROR)
      IF (.NOT.ERROR) THEN
      WRITE(U,'(A)',ERR=7) 'PSF '
      CALL WRTITL(U,2)
      WRITE(U,'(/I8,A)',ERR=7) NATOM, ' !NATOM'
      WRITE(U,'(I8,1X,A,1X,A,1X,A,1X,A,1X,A,1X,2G14.6,I8)',ERR=7)
     &                    (I,SEGID(I),RESID(I),RES(I),TYPE(I),IAC(I),
     &                     CG(I),AMASS(I),0,I=1,NATOM)
      WRITE(U,'(/I8,A)',ERR=7) NBOND, ' !NBOND: bonds'
      WRITE(U,'(8I8)',ERR=7) (IB(I),JB(I),I=1,NBOND)
      WRITE(U,'(/I8,A)',ERR=7) NTHETA,' !NTHETA: angles'
      WRITE(U,'(9I8)',ERR=7) (IT(I),JT(I),KT(I),I=1,NTHETA)
      WRITE(U,'(/I8,A)',ERR=7) NPHI,' !NPHI: dihedrals'
      WRITE(U,'(8I8)',ERR=7) (IP(I),JP(I),KP(I),LP(I),I=1,NPHI)
      WRITE(U,'(/I8,A)',ERR=7) NIMPHI,' !NIMPHI: impropers'
      WRITE(U,'(8I8)',ERR=7)(IM(I),JM(I),KM(I),LM(I),I=1,NIMPHI)
      WRITE(U,'(/I8,A)',ERR=7) DUMMY,' !NDON: donors'
      WRITE(U,'(8I8)',ERR=7) (DUMMY,DUMMY,I=1,DUMMY)
      WRITE(U,'(/I8,A)',ERR=7) DUMMY,' !NACC: acceptors'
      WRITE(U,'(8I8)',ERR=7) (DUMMY,DUMMY,I=1,DUMMY)
      WRITE(U,'(/I8,A)',ERR=7) NNB,' !NNB'
      WRITE(U,'(8I8)',ERR=7) (INB(I),I=1,NNB)
      WRITE(U,'(8I8)',ERR=7) (IBLO(I),I=1,NATOM)
      WRITE(U,'(/2I8,A)',ERR=7) NGRP,0,' !NGRP'
      WRITE(U,'(9I8)',ERR=7) (IGPBS(I),0,0,I=1,NGRP)
C
      CALL VCLOSE(U,'KEEP',ERROR)
      END IF
C-error-labels----
      GOTO 77
7     ERROR=.TRUE.
77    CONTINUE
C-----------------
      IF (ERROR) WRITE(6,'(A)') ' %PSFWRT-ERR: ERROR DURING WRITE'
      RETURN
      END
