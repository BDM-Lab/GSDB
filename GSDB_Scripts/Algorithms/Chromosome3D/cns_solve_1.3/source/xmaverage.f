      SUBROUTINE XMAVER(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRTR,XRINTR,
     &     XRHONUM,XRHONAM,HPRRHO,HPIRHO,HPRHOMA,NRHO,IRHO,NMASK)
C
C Map Averageing.
C
C Author: Paul Adams
C ========================================
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'xmask.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'comand.inc'
      INTEGER NA, NB, NC
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRINTR(3,3),XRTR(3,3)
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER HPRHOMA, IRHO, NRHO, NMASK
C local
      DOUBLE PRECISION MAVMAT(3,3),MAVVEC(3),CENVEC(3),DET
      DOUBLE PRECISION MAVTMP(MXMAVG,3,4)
      DOUBLE PRECISION MAVOP(MXMAVG,3,4)
      INTEGER IFROM, IMASK, IMAP, MFROM(MXMAVG), NMAVGR, I, J
      INTEGER AMIN,AMAX,BMIN,BMAX,CMIN,CMAX
      LOGICAL ERR, COND
      LOGICAL QSELE,QMAP,QVECT,QMATR,QCENT
C pointer
      INTEGER HPRMAP, HPIMAP, HPMASK(MXMAVG)
      INTEGER ACCUMU, NINPOL, GRIDRT, MASKTMP, MAPTMP
C parameter
      DOUBLE PRECISION SMALL, ONE
      PARAMETER (SMALL=1.0D-3,ONE=1.0D0)
C begin
C default values
      ERR=.FALSE.
      IFROM=1
      NMAVGR=0
      QSELE=.FALSE.
      QMAP=.FALSE.
      QVECT=.FALSE.
      QMATR=.FALSE.
      QCENT=.FALSE.
      DO I=1,MXMAVG
      HPMASK(I)=0
      MFROM(I)=0
      END DO
C
C make default selection:
      IF (HPRHOMA.EQ.0) THEN
      CALL WRNDIE(-5,'MAPAVER','No maps are defined.')
      END IF
C
C parsing
      CALL PUSEND('MAPaver>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('MAPaver>')
      CALL MISCOM('MAPaver>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-average')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FROM') THEN
      CALL NEXTWD('FROM=')
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD('FROM=')
      IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(1X,2A)') 'FROM=',XRHONAM(IFROM)
      ELSE
C
      COND=.FALSE.
      DO I=1,XRHONUM
      IF (WD(1:WDLEN).EQ.XRHONAM(I)) THEN
      COND=.TRUE.
      IFROM=I
      END IF
      END DO
C
      IF (.NOT.COND) THEN
      WRITE(6,'(3A)') ' %MAPAVER-ERR: real space object ',
     &     WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'MAPAVER','object undeclared.')
      ERR=.TRUE.
      END IF
      END IF
C =====================================================================
      ELSE IF (WD(1:4).EQ.'MANI') THEN
C
      CALL PUSEND('MAPaver-manipulate>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('MAPaver-manipulate>')
      CALL MISCOM('MAPaver-manipulate>',USED)
      IF (.NOT.USED) THEN
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-average-manipulate')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MASK') THEN
      CALL NEXTWD('MASK=')
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD('MASK=')
C
      COND=.FALSE.
      DO I=1,XRHONUM
      IF (WD(1:WDLEN).EQ.XRHONAM(I)) THEN
      COND=.TRUE.
      IMASK=I
      END IF
      END DO
C
      IF (.NOT.COND) THEN
      WRITE(6,'(3A)') ' %MAPAVER-ERR: real space object ',
     &     WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'MAPAVER','object undeclared.')
      ERR=.TRUE.
      END IF
C
      IF (.NOT.ERR.AND.HPRRHO(IMASK).EQ.0) THEN
      WRITE(6,'(3A)') ' %MAPAVER-ERR: real space object ',
     &     XRHONAM(IMASK),' undefined.'
      CALL WRNDIE(-5,'MAPAVER','object undefined.')
      ELSEIF (.NOT.ERR) THEN
      QSELE=.TRUE.
      ENDIF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MAP') THEN
      CALL NEXTWD('MAP=')
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD('MAP=')
C
      COND=.FALSE.
      DO I=1,XRHONUM
      IF (WD(1:WDLEN).EQ.XRHONAM(I)) THEN
      COND=.TRUE.
      IMAP=I
      END IF
      END DO
C
      IF (.NOT.COND) THEN
      WRITE(6,'(3A)') ' %MAPAVER-ERR: real space object ',
     &     WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'MAPAVER','object undeclared.')
      ERR=.TRUE.
      END IF
C
      IF (.NOT.ERR.AND.HPRRHO(IMAP).EQ.0) THEN
      WRITE(6,'(3A)') ' %MAPAVER-ERR: real space object ',
     &     XRHONAM(IMAP),' undefined.'
      CALL WRNDIE(-5,'MAPAVER','object undefined.')
      ELSEIF (.NOT.ERR) THEN
      QMAP=.TRUE.
      ENDIF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TRAN') THEN
      CALL NEXTVF('TRANslation= ',MAVVEC)
      QVECT=.TRUE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'CENT') THEN
      CALL NEXTVF('CENTer= ',CENVEC)
      QCENT=.TRUE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MATR'.OR.WD(1:4).EQ.'EULE'.OR.
     &     WD(1:4).EQ.'LATT'.OR.WD(1:4).EQ.'SPHE'.OR.
     &     WD(1:4).EQ.'AXIS') THEN
      CALL MATPAR(WD(1:4),MAVMAT)
C check determinant of matrix
      DET=(MAVMAT(1,1)*MAVMAT(2,2)-MAVMAT(1,2)*MAVMAT(2,1))*MAVMAT(3,3)
     &   +(MAVMAT(2,1)*MAVMAT(3,2)-MAVMAT(2,2)*MAVMAT(3,1))*MAVMAT(1,3)
     &   +(MAVMAT(3,1)*MAVMAT(1,2)-MAVMAT(3,2)*MAVMAT(1,1))*MAVMAT(2,3)
      IF (ABS(DET-ONE).GT.SMALL) THEN
      WRITE(6,'(A,F8.5)') ' %XMAVRT-ERR: invalid determinant =',DET
      CALL WRNDIE(-5,'MAPAVER',
     &     'unsuitable rotation matrix given.')
      END IF
      QMATR=.TRUE.
C======================================================================
      ELSE
      CALL CHKEND('MAPaver-manipulate>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (QSELE) THEN
C
C if undefined set center of rotation to zero
      IF (.NOT.QCENT) THEN
      DO I=1,3
      CENVEC(I)=0.0D0
      END DO
      END IF
C
C if undefined set rotation to identity operator
      IF (.NOT.QMATR) THEN
      DO I=1,3
      DO J=1,3
      IF (I.EQ.J) THEN
      MAVMAT(I,J)=1.0D0
      ELSE
      MAVMAT(I,J)=0.0D0
      END IF
      END DO
      END DO
      END IF
C
C if undefined set translations to zero
      IF (.NOT.QVECT) THEN
      DO I=1,3
      MAVVEC(I)=0.0D0
      END DO
      END IF
C
C transfer matrix and translation to one array
      DO I=1,3
         DO J=1,3
            MAVTMP(1,I,J)=MAVMAT(I,J)
         END DO
         MAVTMP(1,I,4)=MAVVEC(I)
      END DO
C
      GRIDRT=ALLHP(IREAL8(XRNSYM*3*4))
C
C find the limits of the mask after transformation
      CALL MSKDIM(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,XRTR,XRINTR,
     &     HEAP(HPRHOMA),HEAP(HPRRHO(IMASK)),HEAP(GRIDRT),
     &     MAVTMP,CENVEC,ONE,AMIN,AMAX,BMIN,BMAX,CMIN,CMAX)
C
C allocate space for transformed mask
      MASKTMP=ALLHP(IREAL4((AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
      IF (QMAP) THEN
      MAPTMP=ALLHP(IREAL4((AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
C
C transform the mask
      CALL MSKTRN(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,XRTR,XRINTR,
     &     HEAP(HPRHOMA),HEAP(HPRRHO(IMASK)),HEAP(HPRRHO(IMAP)),
     &     HEAP(GRIDRT),MAVTMP,CENVEC,AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,
     &     HEAP(MASKTMP),QMAP,HEAP(MAPTMP))
C
      CALL FREHP(MAPTMP,
     &           IREAL4((AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
      ELSE
C
C transform the mask
      CALL MSKTRN(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,XRTR,XRINTR,
     &     HEAP(HPRHOMA),HEAP(HPRRHO(IMASK)),0,
     &     HEAP(GRIDRT),MAVTMP,CENVEC,AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,
     &     HEAP(MASKTMP),QMAP,0)
C
      END IF
      CALL FREHP(MASKTMP,
     &           IREAL4((AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
      CALL FREHP(GRIDRT,IREAL8(XRNSYM*3*4))
C
      ELSE
      CALL WRNDIE(-1,'MAPAVER','mask selection undefined.')
      END IF
C
      QSELE=.FALSE.
      QMAP=.FALSE.
      QMATR=.FALSE.
      QVECT=.FALSE.
      QCENT=.FALSE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'GROU') THEN
C
C increment group counter
      NMAVGR=NMAVGR+1
C
C check to see if we have too many groups
      IF (NMAVGR.GT.MXMAVG) THEN
      WRITE(6,'(2A)') ' %MAPAVER-ERR: maximum number of groups ',
     &     'exceeded - increase MXMAVG.'
      CALL WRNDIE(-5,'MAPAVER','too many groups.')
      END IF
C
      HPMASK(NMAVGR)=0
C
      CALL PUSEND('MAPaver-group>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('MAPaver-group>')
      CALL MISCOM('MAPaver-group>',USED)
      IF (.NOT.USED) THEN
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-average-group')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MASK') THEN
      CALL NEXTWD('MASK=')
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD('MASK=')
C
      COND=.FALSE.
      DO I=1,XRHONUM
      IF (WD(1:WDLEN).EQ.XRHONAM(I)) THEN
      COND=.TRUE.
      IMASK=I
      END IF
      END DO
C
      IF (.NOT.COND) THEN
      WRITE(6,'(3A)') ' %MAPAVER-ERR: real space object ',
     &     WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'MAPAVER','object undeclared.')
      ERR=.TRUE.
      END IF
C
      IF (.NOT.ERR.AND.HPRRHO(IMASK).EQ.0) THEN
      WRITE(6,'(3A)') ' %MAPAVER-ERR: real space object ',
     &     XRHONAM(IMASK),' undefined.'
      CALL WRNDIE(-5,'MAPAVER','object undefined.')
      ELSEIF (.NOT.ERR) THEN
      HPMASK(NMAVGR)=HPRRHO(IMASK)
      MFROM(NMAVGR)=IMASK
      QSELE=.TRUE.
      ENDIF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TRAN') THEN
      CALL NEXTVF('TRANslation= ',MAVVEC)
C add this vector to the list
      DO I=1,3
      MAVOP(NMAVGR,I,4)=MAVVEC(I)
      END DO
      QVECT=.TRUE.
C==================================================================
      ELSE IF (WD(1:4).EQ.'MATR'.OR.WD(1:4).EQ.'EULE'.OR.
     &     WD(1:4).EQ.'LATT'.OR.WD(1:4).EQ.'SPHE'.OR.
     &     WD(1:4).EQ.'AXIS') THEN
      CALL MATPAR(WD(1:4),MAVMAT)
C check determinant of matrix
      DET=(MAVMAT(1,1)*MAVMAT(2,2)-MAVMAT(1,2)*MAVMAT(2,1))*MAVMAT(3,3)
     &   +(MAVMAT(2,1)*MAVMAT(3,2)-MAVMAT(2,2)*MAVMAT(3,1))*MAVMAT(1,3)
     &   +(MAVMAT(3,1)*MAVMAT(1,2)-MAVMAT(3,2)*MAVMAT(1,1))*MAVMAT(2,3)
      IF (ABS(DET-ONE).GT.SMALL) THEN
      WRITE(6,'(A,F8.5)') ' %XMAVRT-ERR: invalid determinant =',DET
      CALL WRNDIE(-5,'MAPAVER',
     &     'unsuitable rotation matrix given.')
      END IF
C add this matrix to the list
      DO I=1,3
      DO J=1,3
      MAVOP(NMAVGR,I,J)=MAVMAT(I,J)
      END DO
      END DO
      QMATR=.TRUE.
C======================================================================
      ELSE
      CALL CHKEND('MAPaver-group>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
C check that all things are defined for the group
C if not give errors and decrease the group counter
      IF (.NOT.(QVECT.AND.QMATR)) THEN
      HPMASK(NMAVGR)=0
      NMAVGR=NMAVGR-1
      IF (.NOT.(QVECT)) THEN
      CALL WRNDIE(-1,'MAPAVER','group translation undefined.')
      ELSE IF (.NOT.(QMATR)) THEN
      CALL WRNDIE(-1,'MAPAVER','group matrix undefined.')
      END IF
      ELSE
      QSELE=.FALSE.
      QVECT=.FALSE.
      QMATR=.FALSE.
      END IF
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' |-----------------------',
     &     ' Mapaver -----------------------|'
      WRITE(6,'(A,A10,A)') ' | From = ',XRHONAM(IFROM),
     &     '                                     |'
      DO I=1,NMAVGR
      WRITE(6,'(2A)') ' |-----------------------',
     &     '--------------------------------|'
      WRITE(6,'(A,I4,A)')
     & ' | Group ',I,'                                            |'
      IF (MFROM(I).NE.0) THEN
      WRITE(6,'(A,A10,A)') ' | Mask = ',XRHONAM(MFROM(I)),
     &     '                                     |'
      END IF
      WRITE(6,'(A,3F10.4,A)')
     & ' | Matrix = ( ',(MAVOP(I,1,J),J=1,3),' )           |'
      WRITE(6,'(A,3F10.4,A)')
     & ' |            ',(MAVOP(I,2,J),J=1,3),' )           |'
      WRITE(6,'(A,3F10.4,A)')
     & ' |            ',(MAVOP(I,3,J),J=1,3),' )           |'
      WRITE(6,'(A,3F10.4,A)')
     & ' | Vector = ( ',(MAVOP(I,J,4),J=1,3),' )           |'
      END DO
      WRITE(6,'(2A)') ' |-----------------------',
     &     '--------------------------------|'
C=====================================================================
      ELSE
      CALL CHKEND('MAPaver>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (NMAVGR.GT.0) THEN
C
      IF (.NOT.ERR.AND.HPRRHO(IFROM).EQ.0) THEN
      WRITE(6,'(3A)') ' %MAPAVER-ERR: real space object ',
     &     XRHONAM(IFROM),' undefined.'
      CALL WRNDIE(-5,'MAPAVER','object undefined.')
      ELSEIF (.NOT.ERR) THEN
      HPRMAP=HPRRHO(IFROM)
      HPIMAP=HPIRHO(IFROM)
C
      IF (HPRMAP.EQ.0.OR.HPRHOMA.EQ.0) THEN
      WRITE(6,'(3A)') ' %MAPAVER-ERR: real space object ',
     &     XRHONAM(IFROM),' empty.'
      CALL WRNDIE(-5,'MAPAVER','object empty.')
      END IF
C
C allocate space for temporary arrays
      ACCUMU=ALLHP(IREAL4(NRHO))
      NINPOL=ALLHP(INTEG4(NRHO))
      GRIDRT=ALLHP(IREAL8(XRNSYM*3*4))
C
C do the averageing
      CALL XMAVER2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,XRTR,XRINTR,
     &     HEAP(HPRMAP),HEAP(HPRHOMA),NRHO,HPMASK,
     &     HEAP(ACCUMU),HEAP(NINPOL),HEAP(GRIDRT),MAVOP,NMAVGR)
C
C free-up the space
      CALL FREHP(GRIDRT,IREAL8(XRNSYM*3*4))
      CALL FREHP(NINPOL,INTEG4(NRHO))
      CALL FREHP(ACCUMU,IREAL4(NRHO))
      END IF
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMAVER2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,XRTR,XRINTR,
     &     RRHO,RHOMASK,NRHO,HPMASK,ACCUMU,NINPOL,GRIDRT,MAVOP,NMAVGR)
C
C Average density in map using defined masks
C
C Author: Paul Adams
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'xmask.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
C I/O
      INTEGER NA, NB, NC
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER XRNSYM, XRSYTH, XRMSYM
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER NMAVGR
      DOUBLE PRECISION XRTR(3,3),XRINTR(3,3)
      DOUBLE PRECISION GRIDRT(XRNSYM,3,4), MAVOP(MXMAVG,3,4)
      INTEGER NRHO
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL ACCUMU(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER NINPOL(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER HPMASK(*)
C local
      INTEGER I,J,K,A,B,C
      INTEGER AMIN,AMAX,BMIN,BMAX,CMIN,CMAX
      DOUBLE PRECISION CORREL,ACORREL
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION MATTMP(MXMAVG,3,4), CENVEC(3)
      INTEGER MAXWORD, ADDLEN, WORDLEN
      PARAMETER (MAXWORD=80)
      CHARACTER*(MAXWORD) WORD, ADDWORD
C pointer
      INTEGER TEMP, MASKTMP, TRNMAT, PRIMSK
C parameter
      INTEGER IZERO,ONE
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0,IZERO=0,ONE=1)
C begin
C
C load the electron density into the accumulation map
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      ACCUMU(A,B,C)=RRHO(A,B,C)
      NINPOL(A,B,C)=ONE
      END DO
      END DO
      END DO
      ACORREL=ZERO
C
C add the correlation for primary mask into tables
C (for completeness only - is always one)
      WORD='AV_CORR_OP_'
      WORDLEN=11
      CALL ENCODI(1,ADDWORD,MAXWORD,ADDLEN)
      CALL ADDST(WORD,MAXWORD,WORDLEN,ADDWORD,ADDLEN)
      CALL DECLAR( WORD(1:WORDLEN),'DP',' ', DBCOMP, 1.0D0 )
C
C allocate the temporary arrays needed
      TEMP=ALLHP(IREAL4(NRHO))
C
C loop over groups minus identity
      DO I=2,NMAVGR
C
C precompute the transformation matrices
      CALL PREMAT(NA,NB,NC,XRNSYM,XRINTR,MAVOP,I,
     &            XRTR,GRIDRT,'NONI',XRNSYM,XRMSYM,
     &            XRSYTH,XRSYMM,XRITSY)
C
C transform, interpolate density into primary mask space
      CALL MAVINT(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,RRHO,RHOMASK,NRHO,
     &     HEAP(HPMASK(1)),GRIDRT,
     &     ACCUMU,CORREL,NINPOL,HEAP(TEMP),'FORW')
C
C put the correlation for this map operation into the word table
      WORD='AV_CORR_OP_'
      WORDLEN=11
      CALL ENCODI(I,ADDWORD,MAXWORD,ADDLEN)
      CALL ADDST(WORD,MAXWORD,WORDLEN,ADDWORD,ADDLEN)
      CALL DECLAR( WORD(1:WORDLEN),'DP',' ', DBCOMP, CORREL )
C
C accumulate average correlations - exclude first operator
C as this should always be equal to one
      IF (I.GT.ONE) THEN
         ACORREL=ACORREL+CORREL
      END IF
C
      END DO
C
      ACORREL=ACORREL/(NMAVGR-1)
      WORD='AV_CORR'
      WORDLEN=7
      CALL DECLAR( WORD(1:WORDLEN),'DP',' ', DBCOMP, ACORREL )
C
C average the electron density in accumulation map (in primary mask)
      CALL MAVDEN(MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     RHOMASK,HEAP(HPMASK(1)),ACCUMU,NINPOL)
C
C loop over all groups (1st is identity operator)
      DO I=1,NMAVGR
C
C precompute the inverse transformation matrix
      CALL PREMAT(NA,NB,NC,XRNSYM,XRINTR,MAVOP,I,
     &            XRTR,GRIDRT,'INVE',XRNSYM,XRMSYM,
     &            XRSYTH,XRSYMM,XRITSY)
C
C generate mask if not defined already
      IF (HPMASK(I).EQ.0) THEN
C
C generate space for temporary copy of primary mask
      PRIMSK=ALLHP(IREAL4(NRHO))
C
C copy primary mask to temporary space
      CALL MSKCPY(HEAP(HPMASK(1)),HEAP(PRIMSK),
     &            MAASY,MBASY,MCASY,NAASY,NBASY,NCASY)
C
C operator (group 1 to group I)
      DO K=1,3
      DO J=1,4
      MATTMP(1,K,J)=MAVOP(I,K,J)
      END DO
      CENVEC(K)=ZERO
      END DO
C
      TRNMAT=ALLHP(IREAL8(XRNSYM*3*4))
C
C find the limits of the mask after transformation
      CALL MSKDIM(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,XRTR,XRINTR,
     &     RHOMASK,HEAP(PRIMSK),HEAP(TRNMAT),
     &     MATTMP,CENVEC,1.0D0,AMIN,AMAX,BMIN,BMAX,CMIN,CMAX)
C
C allocate space for transformed mask
      MASKTMP=ALLHP(IREAL4((AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
C
C transform the mask
      CALL MSKTRN(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,XRTR,XRINTR,
     &     RHOMASK,HEAP(PRIMSK),0,
     &     HEAP(TRNMAT),MATTMP,CENVEC,AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,
     &     HEAP(MASKTMP),.FALSE.,0)
C
      CALL FREHP(MASKTMP,
     &           IREAL4((AMAX-AMIN+1)*(BMAX-BMIN+1)*(CMAX-CMIN+1)))
      CALL FREHP(TRNMAT,IREAL8(XRNSYM*3*4))
C
C transform, interpolate density into secondary mask space
      CALL MAVINT(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,ACCUMU,RHOMASK,NRHO,
     &     HEAP(PRIMSK),GRIDRT,
     &     RRHO,CORREL,NINPOL,HEAP(TEMP),'BACK')
C
C free up space for temporary copy of primary mask
      CALL FREHP(PRIMSK,IREAL4(NRHO))
C
      ELSE
C
C transform, interpolate density into secondary mask space
      CALL MAVINT(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,ACCUMU,RHOMASK,NRHO,
     &     HEAP(HPMASK(I)),GRIDRT,
     &     RRHO,CORREL,NINPOL,HEAP(TEMP),'BACK')
      END IF
C
      END DO
C
C free space for temporary arrays
      CALL FREHP(TEMP,IREAL4(NRHO))
C
      RETURN
      END
C======================================================================
      SUBROUTINE MAVINT(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,RRHO,RHOMASK,NRHO,MMASK1,
     &     RTMAT,ACCUMU,CORREL,NINPOL,TEMP,PASS)
C
C transform the points of mask1
C interpolate from the valid points in mask2 into mask1
C
C Author: Paul Adams
C
      IMPLICIT NONE
      INCLUDE 'xmask.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'funct.inc'
C I/O
      CHARACTER*4 PASS
      INTEGER NA, NB, NC
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER XRNSYM, XRSYTH, XRMSYM
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER NRHO
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL MMASK1(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION RTMAT(XRNSYM,3,4),CORREL
      REAL ACCUMU(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL TEMP(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER NINPOL(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
C local
      INTEGER A,B,C,AN1,AN2,BN1,BN2,CN1,CN2
      INTEGER AA,BB,CC,AT,BT,CT
      INTEGER MSO,SYM,AS,BS,CS,NPOINT
      DOUBLE PRECISION AP,BP,CP,DX,DY,DZ
      DOUBLE PRECISION INTERP,SUMINT,SUMORI,TOPSUM,BSUM1,BSUM2
      DOUBLE PRECISION LRRHO(2,2,2),RSYTH
      INTEGER LTX,LTY,LTZ
      DOUBLE PRECISION ATRANS,BTRANS,CTRANS
      LOGICAL COND
C parameter
      INTEGER IZERO,ONE
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0,IZERO=0,ONE=1)
C
C initialise
      RSYTH=XRSYTH
      CORREL=ZERO
      SUMINT=ZERO
      SUMORI=ZERO
      NPOINT=IZERO
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      TEMP(A,B,C)=ZERO
      END DO
      END DO
      END DO
C
C loop over points in primary mask
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
C
      IF (RHOMASK(A,B,C).GT.IZERO.AND.MMASK1(A,B,C).LE.ZERO) THEN
C
C get the symmetry operation which was applied to map this grid point
C into the ASU
      CALL MSKA2M(MMASK1(A,B,C),MSO,LTX,LTY,LTZ,XRMSYM)
C
      ATRANS=A-(LTX*NA)
      BTRANS=B-(LTY*NB)
      CTRANS=C-(LTZ*NC)
C
C apply the transformation to off-grid space
      AP=RTMAT(MSO,1,1)*ATRANS+RTMAT(MSO,1,2)*BTRANS
     &  +RTMAT(MSO,1,3)*CTRANS+RTMAT(MSO,1,4)
      BP=RTMAT(MSO,2,1)*ATRANS+RTMAT(MSO,2,2)*BTRANS
     &  +RTMAT(MSO,2,3)*CTRANS+RTMAT(MSO,2,4)
      CP=RTMAT(MSO,3,1)*ATRANS+RTMAT(MSO,3,2)*BTRANS
     &  +RTMAT(MSO,3,3)*CTRANS+RTMAT(MSO,3,4)
C
C the 8 nearest on-grid points
      AN1=INT(AP)
      AN2=AN1+SIGN(1.0D0,AP)
      BN1=INT(BP)
      BN2=BN1+SIGN(1.0D0,BP)
      CN1=INT(CP)
      CN2=CN1+SIGN(1.0D0,CP)
C
C calculate the distance from the point to the grid
      DX=AP-MIN(AN1,AN2)
      DY=BP-MIN(BN1,BN2)
      DZ=CP-MIN(CN1,CN2)
C
C loop over nearest points on-grid in mask 2
      CT=ONE
      DO CC=MIN(CN1,CN2),MAX(CN1,CN2)
      BT=ONE
      DO BB=MIN(BN1,BN2),MAX(BN1,BN2)
      AT=ONE
      DO AA=MIN(AN1,AN2),MAX(AN1,AN2)
C
C is the point in the ASU box
      COND=(AA.GE.MAASY.AND.AA.LE.NAASY.AND.
     &      BB.GE.MBASY.AND.BB.LE.NBASY.AND.
     &      CC.GE.MCASY.AND.CC.LE.NCASY)
C
C if the point is in the ASU box test to see if its within RHOMASK
      IF (COND) COND=RHOMASK(AA,BB,CC).GT.0
C
C if the above are false then we have to remap back within RHOMASK
C by trying all symmetry operators
      IF (.NOT.COND) THEN
      SYM=ONE
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
      AS=MOD(XRSYMM(SYM,1,1)*AA+XRSYMM(SYM,1,2)*BB+XRSYMM(SYM,1,3)*CC
     &     +(XRSYMM(SYM,1,4)*NA)/XRSYTH+10000*NA-MAASY,NA) +MAASY
      BS=MOD(XRSYMM(SYM,2,1)*AA+XRSYMM(SYM,2,2)*BB+XRSYMM(SYM,2,3)*CC
     &     +(XRSYMM(SYM,2,4)*NB)/XRSYTH+10000*NB-MBASY,NB) +MBASY
      CS=MOD(XRSYMM(SYM,3,1)*AA+XRSYMM(SYM,3,2)*BB+XRSYMM(SYM,3,3)*CC
     &     +(XRSYMM(SYM,3,4)*NC)/XRSYTH+10000*NC-MCASY,NC) +MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND) COND=RHOMASK(AS,BS,CS).GT.0
      SYM=SYM+ONE
      END DO
C
C the point was in the ASU to start with so use the initial values
      ELSE
      AS=AA
      BS=BB
      CS=CC
      END IF
C
C error if we can't map it back
      IF (.NOT.COND) THEN
      CALL WRNDIE(-5,'XMAVER',
     &     'Internal error: grid point outside ASU')
      END IF
C
C store the density locally
      LRRHO(AT,BT,CT)=RRHO(AS,BS,CS)
C
      AT=AT+ONE
      END DO
      BT=BT+ONE
      END DO
      CT=CT+ONE
      END DO
C
C do the 8-point linear interpolation
      INTERP=LRRHO(1,1,1)*(1-DX)*(1-DY)*(1-DZ)+
     &       LRRHO(2,1,1)*DX*(1-DY)*(1-DZ)+
     &       LRRHO(1,2,1)*(1-DX)*DY*(1-DZ)+
     &       LRRHO(2,2,1)*DX*DY*(1-DZ)+
     &       LRRHO(1,1,2)*(1-DX)*(1-DY)*DZ+
     &       LRRHO(2,1,2)*DX*(1-DY)*DZ+
     &       LRRHO(1,2,2)*(1-DX)*DY*DZ+
     &       LRRHO(2,2,2)*DX*DY*DZ
C
      NPOINT=NPOINT+ONE
      IF (PASS.EQ.'FORW') THEN
C forward pass
C update the accumulation, and interpolation counter arrays
      ACCUMU(A,B,C)=ACCUMU(A,B,C)+INTERP
      TEMP(A,B,C)=INTERP
      SUMINT=SUMINT+INTERP
      SUMORI=SUMORI+RRHO(A,B,C)
      NINPOL(A,B,C)=NINPOL(A,B,C)+ONE
C backward pass
      ELSEIF (PASS.EQ.'BACK') THEN
      ACCUMU(A,B,C)=INTERP
      ELSE
      CALL WRNDIE(-5,'XMAVER',
     &     'Programming error: unknown MAVINT mode')
      END IF
C
      END IF
      END DO
      END DO
      END DO
C
      IF (PASS.EQ.'FORW') THEN
C forward pass
C correlation function
      SUMINT=SUMINT/NPOINT
      SUMORI=SUMORI/NPOINT
      TOPSUM=ZERO
      BSUM1=ZERO
      BSUM2=ZERO
      NPOINT=IZERO
C loop over points in primary mask
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).GT.IZERO.AND.MMASK1(A,B,C).LE.ZERO) THEN
         TOPSUM=TOPSUM+((RRHO(A,B,C)*TEMP(A,B,C))-(SUMINT*SUMORI))
         BSUM1=BSUM1+((RRHO(A,B,C)*RRHO(A,B,C))-(SUMORI*SUMORI))
         BSUM2=BSUM2+((TEMP(A,B,C)*TEMP(A,B,C))-(SUMINT*SUMINT))
         NPOINT=NPOINT+ONE
      END IF
      END DO
      END DO
      END DO
C
      CORREL=(TOPSUM/NPOINT)/SQRT((BSUM1/NPOINT)*(BSUM2/NPOINT))
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE MAVDEN(MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     RHOMASK,MMASK,ACCUMU,NINPOL)
C
C average the interpolated density in the mask
C
C Author: Paul Adams
C
      IMPLICIT NONE
C I/O
      INTEGER MAASY,MBASY,MCASY,NAASY,NBASY,NCASY
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL MMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL ACCUMU(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER NINPOL(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
C local
      INTEGER A,B,C
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
C loop over points in mask
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).GT.0.AND.MMASK(A,B,C).LE.ZERO) THEN
      IF (NINPOL(A,B,C).GT.0) THEN
      ACCUMU(A,B,C)=ACCUMU(A,B,C)/NINPOL(A,B,C)
      ELSE
      ACCUMU(A,B,C)=ZERO
      END IF
      END IF
      END DO
      END DO
      END DO
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE MSKTRN(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,XRTR,XRINTR,
     &     RHOMASK,MASK,MAP,RTMAT,MATRIX,CENTER,
     &     AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,TEMP,QMAP,MTEMP)
C
C given a mask and rotation matrix R and translation T
C return the transformed mask and map if defined
C
C Author: Paul Adams
C
      IMPLICIT NONE
      INCLUDE 'xmask.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'funct.inc'
C I/O
      INTEGER NA, NB, NC
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER AMIN,AMAX,BMIN,BMAX,CMIN,CMAX
      INTEGER XRNSYM, XRSYTH, XRMSYM
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL MASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL MAP(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL TEMP(AMIN:AMAX,BMIN:BMAX,CMIN:CMAX)
      REAL MTEMP(AMIN:AMAX,BMIN:BMAX,CMIN:CMAX)
      DOUBLE PRECISION RTMAT(XRNSYM,3,4),XRTR(3,3),XRINTR(3,3)
      DOUBLE PRECISION MATRIX(MXMAVG,3,4),CENTER(3)
      LOGICAL QMAP
C local
      INTEGER A,B,C,AN1,AN2,BN1,BN2,CN1,CN2,CT,CC,BT,BB,AT,AA
      INTEGER SYM,AS,BS,CS,ASS,BSS,CSS,ASSS,BSSS,CSSS
      DOUBLE PRECISION CENA,CENB,CENC,A1,B1,C1
      REAL KEY
      DOUBLE PRECISION AP,BP,CP,DX,DY,DZ,OVPCNT
      DOUBLE PRECISION LMASK(2,2,2),LMAP(2,2,2),INTERP,MINTERP
      INTEGER ATRANS,BTRANS,CTRANS,OVPNT,TOPNT
      LOGICAL COND,QOVER
C parameter
      DOUBLE PRECISION ZERO, QUARTER, HALF
      PARAMETER (ZERO=0.0D0,QUARTER=0.25D0,HALF=0.5D0)
C
C initialise
      CENA=CENTER(1)*XRTR(1,1)*NA
      CENB=CENTER(2)*XRTR(2,2)*NB
      CENC=CENTER(3)*XRTR(3,3)*NC
C
C initialise temporary mask
      DO C=CMIN,CMAX
      DO B=BMIN,BMAX
      DO A=AMIN,AMAX
         TEMP(A,B,C)=0
      END DO
      END DO
      END DO
C
C initialise temporary map if needed
      IF (QMAP) THEN
         DO C=CMIN,CMAX
         DO B=BMIN,BMAX
         DO A=AMIN,AMAX
            MTEMP(A,B,C)=ZERO
         END DO
         END DO
         END DO
      END IF
C
C set up the realspace transformation matrices
      CALL PREMAT(NA,NB,NC,XRNSYM,XRINTR,MATRIX,1,
     &            XRTR,RTMAT,'INVE',XRNSYM,XRMSYM,
     &            XRSYTH,XRSYMM,XRITSY)
C
C loop over points in primary mask
      DO C=CMIN,CMAX
      DO B=BMIN,BMAX
      DO A=AMIN,AMAX
C
C apply center of rotation correction
      A1=A-CENA
      B1=B-CENB
      C1=C-CENC
C
C apply the transformation to off-grid space
      AP=RTMAT(1,1,1)*A1+RTMAT(1,1,2)*B1
     &  +RTMAT(1,1,3)*C1+RTMAT(1,1,4)
      BP=RTMAT(1,2,1)*A1+RTMAT(1,2,2)*B1
     &  +RTMAT(1,2,3)*C1+RTMAT(1,2,4)
      CP=RTMAT(1,3,1)*A1+RTMAT(1,3,2)*B1
     &  +RTMAT(1,3,3)*C1+RTMAT(1,3,4)
C
C apply center of rotation correction
      AP=AP+CENA
      BP=BP+CENB
      CP=CP+CENC
C
C the limits of the nearest grid points
      AN1=INT(AP)
      AN2=AN1+INT(SIGN(1.0D0,AP))
      BN1=INT(BP)
      BN2=BN1+INT(SIGN(1.0D0,BP))
      CN1=INT(CP)
      CN2=CN1+INT(SIGN(1.0D0,CP))
C
C calculate the distance from the point to the grid
      DX=AP-MIN(AN1,AN2)
      DY=BP-MIN(BN1,BN2)
      DZ=CP-MIN(CN1,CN2)
C
C loop over nearest points on-grid in mask 2
      CT=1
      DO CC=MIN(CN1,CN2),MAX(CN1,CN2)
      BT=1
      DO BB=MIN(BN1,BN2),MAX(BN1,BN2)
      AT=1
      DO AA=MIN(AN1,AN2),MAX(AN1,AN2)
C
C is the point in the ASU box
      COND=(AA.GE.MAASY.AND.AA.LE.NAASY.AND.
     &      BB.GE.MBASY.AND.BB.LE.NBASY.AND.
     &      CC.GE.MCASY.AND.CC.LE.NCASY)
C
C if the point is in the ASU box test to see if its within RHOMASK
      IF (COND) COND=RHOMASK(AA,BB,CC).GT.0
C
C if the above are false then we have to remap back within RHOMASK
C by trying all symmetry operators
      IF (.NOT.COND) THEN
      SYM=1
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
      AS=MOD(XRSYMM(SYM,1,1)*AA+XRSYMM(SYM,1,2)*BB+XRSYMM(SYM,1,3)*CC
     &     +(XRSYMM(SYM,1,4)*NA)/XRSYTH+10000*NA-MAASY,NA) +MAASY
      BS=MOD(XRSYMM(SYM,2,1)*AA+XRSYMM(SYM,2,2)*BB+XRSYMM(SYM,2,3)*CC
     &     +(XRSYMM(SYM,2,4)*NB)/XRSYTH+10000*NB-MBASY,NB) +MBASY
      CS=MOD(XRSYMM(SYM,3,1)*AA+XRSYMM(SYM,3,2)*BB+XRSYMM(SYM,3,3)*CC
     &     +(XRSYMM(SYM,3,4)*NC)/XRSYTH+10000*NC-MCASY,NC) +MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND) COND=RHOMASK(AS,BS,CS).GT.0
      SYM=SYM+1
      END DO
C
C the point was in the ASU to start with so use the initial values
      ELSE
      AS=AA
      BS=BB
      CS=CC
      END IF
C
C error if we can't map it back
      IF (.NOT.COND) THEN
      CALL WRNDIE(-5,'MSKTRN',
     &     'Internal error: grid point outside ASU')
      END IF
C
C store the density and mask flag locally
      IF (MASK(AS,BS,CS).LE.0) THEN
C
C check to see if this ASU lattice point actually matches the
C real space one we are concerned with - ie remove any xtal
C symmetry duplication
         CALL MSKA2M(MASK(AS,BS,CS),SYM,ATRANS,BTRANS,CTRANS,
     &               XRMSYM)
         ASS=AS-(ATRANS*NA)-(XRSYMM(SYM,1,4)*NA)/XRSYTH
         BSS=BS-(BTRANS*NB)-(XRSYMM(SYM,2,4)*NB)/XRSYTH
         CSS=CS-(CTRANS*NC)-(XRSYMM(SYM,3,4)*NC)/XRSYTH
         ASSS=XRITSY(SYM,1,1)*ASS+XRITSY(SYM,2,1)*BSS
     &       +XRITSY(SYM,3,1)*CSS
         BSSS=XRITSY(SYM,1,2)*ASS+XRITSY(SYM,2,2)*BSS
     &       +XRITSY(SYM,3,2)*CSS
         CSSS=XRITSY(SYM,1,3)*ASS+XRITSY(SYM,2,3)*BSS
     &       +XRITSY(SYM,3,3)*CSS
C
         IF (ASSS.EQ.AA.AND.BSSS.EQ.BB.AND.CSSS.EQ.CC) THEN
            LMASK(AT,BT,CT)=1
         ELSE
            LMASK(AT,BT,CT)=0
         END IF
      ELSE
         LMASK(AT,BT,CT)=0
      END IF
      IF (QMAP) THEN
         LMAP(AT,BT,CT)=MAP(AS,BS,CS)
      END IF
C
      AT=AT+1
      END DO
      BT=BT+1
      END DO
      CT=CT+1
      END DO
C
C do the 8-point linear interpolation
      INTERP=LMASK(1,1,1)*(1-DX)*(1-DY)*(1-DZ)+
     &       LMASK(2,1,1)*DX*(1-DY)*(1-DZ)+
     &       LMASK(1,2,1)*(1-DX)*DY*(1-DZ)+
     &       LMASK(2,2,1)*DX*DY*(1-DZ)+
     &       LMASK(1,1,2)*(1-DX)*(1-DY)*DZ+
     &       LMASK(2,1,2)*DX*(1-DY)*DZ+
     &       LMASK(1,2,2)*(1-DX)*DY*DZ+
     &       LMASK(2,2,2)*DX*DY*DZ
      IF (QMAP) THEN
         MINTERP=LMAP(1,1,1)*(1-DX)*(1-DY)*(1-DZ)+
     &           LMAP(2,1,1)*DX*(1-DY)*(1-DZ)+
     &           LMAP(1,2,1)*(1-DX)*DY*(1-DZ)+
     &           LMAP(2,2,1)*DX*DY*(1-DZ)+
     &           LMAP(1,1,2)*(1-DX)*(1-DY)*DZ+
     &           LMAP(2,1,2)*DX*(1-DY)*DZ+
     &           LMAP(1,2,2)*(1-DX)*DY*DZ+
     &           LMAP(2,2,2)*DX*DY*DZ
      END IF
C
C is this point inside or outside the mask
      IF (INTERP.GT.QUARTER) THEN
         TEMP(A,B,C)=-1
         IF (QMAP) MTEMP(A,B,C)=MINTERP
      ELSE
         TEMP(A,B,C)=1
      END IF
C
      END DO
      END DO
      END DO
C
C after filling the temporary mask array, remap the points back into
C ASU, overwriting the original mask as we go - should probably make
C some check to see if the transformed mask no longer obeys xtal conditions
C ie. does the mask fit in the ASU with no overlaps.
C
C set the mask to all outside mask
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
         MASK(A,B,C)=1
      END DO
      END DO
      END DO
C
C set overlap flag to false
      QOVER=.FALSE.
      OVPNT=0
      TOPNT=0
C
C loop over points in primary mask
      DO C=CMIN,CMAX
      DO B=BMIN,BMAX
      DO A=AMIN,AMAX
C
C map A, B, C into asymmetric unit
      SYM=1
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
C
      AA=XRSYMM(SYM,1,1)*A+XRSYMM(SYM,1,2)*B+XRSYMM(SYM,1,3)*C
     &       +(XRSYMM(SYM,1,4)*NA)/XRSYTH
      BB=XRSYMM(SYM,2,1)*A+XRSYMM(SYM,2,2)*B+XRSYMM(SYM,2,3)*C
     &       +(XRSYMM(SYM,2,4)*NB)/XRSYTH
      CC=XRSYMM(SYM,3,1)*A+XRSYMM(SYM,3,2)*B+XRSYMM(SYM,3,3)*C
     &       +(XRSYMM(SYM,3,4)*NC)/XRSYTH
C
      AS=MOD(AA+10000*NA-MAASY,NA) +MAASY
      BS=MOD(BB+10000*NB-MBASY,NB) +MBASY
      CS=MOD(CC+10000*NC-MCASY,NC) +MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND)  COND=RHOMASK(AS,BS,CS).GT.0
      SYM=SYM+1
      END DO
C
      IF (COND) THEN
      SYM=SYM-1
C
C update symmetry and lattice translation info
      ATRANS=(AS-AA)/NA
      BTRANS=(BS-BB)/NB
      CTRANS=(CS-CC)/NC
C
C pack the symmetry operator and lattice translation info
      IF (TEMP(A,B,C).LT.0) THEN
        CALL MSKM2A(SYM,ATRANS,BTRANS,CTRANS,KEY,XRMSYM)
        IF (MASK(AS,BS,CS).LE.0) THEN
           QOVER=.TRUE.
           OVPNT=OVPNT+1
        ELSE
           MASK(AS,BS,CS)=KEY
           IF (QMAP) THEN
              MAP(AS,BS,CS)=MTEMP(A,B,C)
           END IF
        END IF
      ELSE
         IF (MASK(AS,BS,CS).GT.0) THEN
            MASK(AS,BS,CS)=TEMP(A,B,C)
         END IF
      END IF
C
      ELSE
      CALL WRNDIE(-5,'MSKTRN','could not remap lattice point to ASU')
      END IF
C
      TOPNT=TOPNT+1
C
      END DO
      END DO
      END DO
C
C if overlap flag set then print a warning
      IF (QOVER) THEN
         OVPCNT=(OVPNT*100.0D0)/TOPNT
         IF (OVPCNT.GT.0.1D0) THEN
            WRITE(6,'(A,F6.2,A)')
     &      ' MSKTRN-WRN: the transformed mask overlaps itself by ',
     &      OVPCNT,'%'
         END IF
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE MSKCPY(MASK1,MASK2,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY)
C
      IMPLICIT NONE
C IO
      INTEGER MAASY,MBASY,MCASY,NAASY,NBASY,NCASY
      REAL MASK1(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL MASK2(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
C local
      INTEGER A,B,C
C begin
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      MASK2(A,B,C)=MASK1(A,B,C)
      END DO
      END DO
      END DO
C
      RETURN
      END
C======================================================================
