      SUBROUTINE XRASSOC(XRFLAG,ANOMFLAG,XRATOF,XRINDF,XRATOM,
     &                   XRINDX,XRFQS,XRNATO,NANOM,XNAMEAS,XASSOC,
     &                   QASELE,XRNATF,XRFDP,QHERM,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &                   XRTR,XRINTR)
C
C Generation of atom flag list for ASSOciate statement given by
C XRFLAG.  Separation into normal and anomalous scatterers.
C
C The array XRFQS defines the occupancy scale factor
C for atoms at special positions.  Routine
C also checks the number of anomalous scattering atoms.  Normal
C scatterers will appear in the XRATOM list from 1 to NANOM.
C Anomalous scatterers will appear from NANOM+1 to XRNATO.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'funct.inc'
      INTEGER XRFLAG(*), ANOMFLAG(*), XRATOF(*), XRINDF(*)
      INTEGER XRATOM(*), XRINDX(*)
      DOUBLE PRECISION XRFQS(*)
      INTEGER XRNATO, NANOM, XNAMEAS
      CHARACTER*(*) XASSOC
      LOGICAL QASELE
      INTEGER XRNATF
      DOUBLE PRECISION XRFDP(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3)
C local
      INTEGER IAT, ISYM, NSPECL, LEN
      DOUBLE PRECISION XX, YY, ZZ, RTH, XXX, YYY, ZZZ
      DOUBLE PRECISION XXS, YYS, ZZS
      LOGICAL QONE, ERR
      DOUBLE PRECISION NBQSPC2
C parameters
      DOUBLE PRECISION ONE, HALF, P01
      PARAMETER (ONE=1.0D0, HALF=0.5D0, P01=0.1D0)
C begin
      NBQSPC2=NBQSPC**2
C
C make default anomalous scatterer selection if required
C (default selection is ( attr scatter_fdp > 0.1 ) )
C
      IF (.NOT.QASELE) THEN
C
C use array XRATOM temporarily to store the indices of all atoms
C with form factor entries
      DO IAT=1,NATOM
      XRATOM(IAT)=0
      END DO
      DO IAT=1,XRNATF
      XRATOM(XRATOF(IAT))=IAT
      END DO
C
C go through all atoms
      DO IAT=1,NATOM
C
C form factor available and atom selected?
      IF (XRATOM(IAT).GT.0.AND.XRFLAG(IAT).EQ.1) THEN
C
C if f'' is greater than 0.1 then declare the atom anomalous
      IF (XRFDP(XRINDF(XRATOM(IAT))).GT.P01) THEN
      ANOMFLAG(IAT)=1
      ELSE
      ANOMFLAG(IAT)=0
      END IF
      END IF
      END DO
      END IF
C
C make list of normal scatterers included in FCALC calculation
      XRNATO=0
      DO IAT=1,XRNATF
      IF (XRFLAG(XRATOF(IAT)).EQ.1.AND.
     &    (ANOMFLAG(XRATOF(IAT)).EQ.0.OR.QHERM)) THEN
C
C temporarily mark the XRFLAG array
      XRFLAG(XRATOF(IAT))=-XRFLAG(XRATOF(IAT))
C
      XRNATO=XRNATO+1
      XRATOM(XRNATO)=XRATOF(IAT)
      XRINDX(XRNATO)=XRINDF(IAT)
      END IF
      END DO
C
      NANOM=XRNATO
C
C now make list of anomalous scatterers: anomalous scatterers
C appear at the end of the list between NANOM and XRNATO
      DO IAT=1,XRNATF
      IF (.NOT.QHERM.AND.
     & (XRFLAG(XRATOF(IAT)).EQ.1.AND.ANOMFLAG(XRATOF(IAT)).EQ.1))
     & THEN
C
C temporarily mark the XRFLAG array
      XRFLAG(XRATOF(IAT))=-XRFLAG(XRATOF(IAT))
C
      XRNATO=XRNATO+1
      XRATOM(XRNATO)=XRATOF(IAT)
      XRINDX(XRNATO)=XRINDF(IAT)
      END IF
      END DO
C
C check if any selected atoms have missing form factors
      ERR=.FALSE.
      DO IAT=1,NATOM
      IF (XRFLAG(IAT).EQ.1) THEN
      ERR=.TRUE.
      WRITE(6,'(10A)')
     & ' %XRASSOC-ERR: missing SCATter definition for ( ',
     & SEGID(IAT),' ',RES(IAT),' ',RESID(IAT),' ',TYPE(IAT),
     & ' ) chemical=', IAC(IAT)
      END IF
      END DO
      IF (ERR) THEN
      CALL WRNDIE(-1,'XRASSOC',
     & 'missing SCATter definition for SELEcted atoms.')
      END IF
C
C reset marks in XRFLAG array
      DO IAT=1,NATOM
      IF (XRFLAG(IAT).LT.0) XRFLAG(IAT)=-XRFLAG(IAT)
      END DO
C
C initialize special position occupancy scale factor
      DO IAT=1,XRNATO
      XRFQS(XRATOM(IAT))=ONE
      END DO
C
C check for unknown coordinates
      ERR=.FALSE.
      DO IAT=1,XRNATO
      IF (.NOT.INITIA(XRATOM(IAT),X,Y,Z)) THEN
      ERR=.TRUE.
      WRITE(6,'(9A)')
     &  ' %XRASSOC-ERR: unknown coordinates for atom "',
     &    SEGID(XRATOM(IAT)),'-',RESID(XRATOM(IAT)),'-',
     &    RES(XRATOM(IAT)),'-',TYPE(XRATOM(IAT)),'"'
      END IF
      END DO
      IF (ERR) THEN
      CALL WRNDIE(-5,'XRASSOC','Unknown coordinates')
      XRNATO=0
      END IF
C
C check occupancies
      QONE=.TRUE.
      DO IAT=1,XRNATO
      QONE=QONE.AND.ABS(QMAIN(XRATOM(IAT))-ONE).LT.RSMALL
      END DO
C
C check for special positions
C ===========================
C For an atom at a special position P we count the number of
C symmetry operators that are identity for P.  Suppose that
C there are m such operators.  It follows that for each
C symmetry mate S of P with S not equal P there are m symmetry
C operators OP such that S = OP [P].  Thus, we can use the general
C structure factor expression for a special position provided we
C divide the structure factors by m.  This is the purpose of the
C XRFQS array.
      NSPECL=0
C
C loop over all symmetry operators and check for special positions
      RTH=XRSYTH
      DO ISYM=2,XRNSYM
      DO IAT=1,XRNATO
C fractionalize coordinates
      XX=XRTR(1,1)*X(XRATOM(IAT))
     &  +XRTR(1,2)*Y(XRATOM(IAT))+XRTR(1,3)*Z(XRATOM(IAT))
      YY=XRTR(2,1)*X(XRATOM(IAT))
     &  +XRTR(2,2)*Y(XRATOM(IAT))+XRTR(2,3)*Z(XRATOM(IAT))
      ZZ=XRTR(3,1)*X(XRATOM(IAT))
     &  +XRTR(3,2)*Y(XRATOM(IAT))+XRTR(3,3)*Z(XRATOM(IAT))
C compute symmetry mate
      XXS=XRSYMM(ISYM,1,1)*XX + XRSYMM(ISYM,1,2)*YY
     & + XRSYMM(ISYM,1,3)*ZZ + XRSYMM(ISYM,1,4)/RTH
      YYS=XRSYMM(ISYM,2,1)*XX + XRSYMM(ISYM,2,2)*YY
     & + XRSYMM(ISYM,2,3)*ZZ + XRSYMM(ISYM,2,4)/RTH
      ZZS=XRSYMM(ISYM,3,1)*XX + XRSYMM(ISYM,3,2)*YY
     & + XRSYMM(ISYM,3,3)*ZZ + XRSYMM(ISYM,3,4)/RTH
C
C compute difference vector
      XX=XX-XXS
      YY=YY-YYS
      ZZ=ZZ-ZZS
C
C compute minimum image difference in fractional coordinates
      XX=INT(ABS(XX)+HALF)*SIGN(ONE,-XX) +XX
      YY=INT(ABS(YY)+HALF)*SIGN(ONE,-YY) +YY
      ZZ=INT(ABS(ZZ)+HALF)*SIGN(ONE,-ZZ) +ZZ
C
C convert into orthogonal coordinates
      XXX=XRINTR(1,1)*XX+XRINTR(1,2)*YY+XRINTR(1,3)*ZZ
      YYY=XRINTR(2,1)*XX+XRINTR(2,2)*YY+XRINTR(2,3)*ZZ
      ZZZ=XRINTR(3,1)*XX+XRINTR(3,2)*YY+XRINTR(3,3)*ZZ
C
C check distance to see if it is a special position
      IF (XXX*XXX+YYY*YYY+ZZZ*ZZZ.LT.NBQSPC2) THEN
      XRFQS(XRATOM(IAT))=XRFQS(XRATOM(IAT))+ONE
      END IF
      END DO
      END DO
C
      DO IAT=1,XRNATO
      XRFQS(XRATOM(IAT))=ONE/XRFQS(XRATOM(IAT))
      IF (ABS(XRFQS(XRATOM(IAT))-ONE).GT.RSMALL) NSPECL=NSPECL+1
      END DO
C
C write some info
      IF (XRNATO.EQ.0) THEN
      WRITE(6,'(A)')
     & ' %XRASSOC-ERR: number of selected scattering atoms is zero'
      ELSE
      IF (WRNLEV.GE.5) THEN
C
      LEN=XNAMEAS
      CALL TRIMM(XASSOC,LEN)
      IF (QONE) THEN
      WRITE(6,'(3A,I7,A,I4,A,I3,A)')
     & ' ',XASSOC(1:LEN),': #scatt.=',XRNATO,
     &         ' #anomalous=',XRNATO-NANOM,
     &         ' #special pos.=',NSPECL,
     &         ' occupancies=1'
      ELSE
      WRITE(6,'(3A,I7,A,I4,A,I3,A)')
     & ' ',XASSOC(1:LEN),': #scatt.=',XRNATO,
     &         ' #anomalous=',XRNATO-NANOM,
     &         ' #special pos.=',NSPECL,
     &         ' occupancies <> 1'
      END IF
      END IF
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE XRENER(E)
C
C Routine computes the energy and derivatives corresponding
C to the specified crystallographic target function.
C Depending on XRLTOL, the routine computes the exact
C derivatives or it uses a linear approximation in x, y, z.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      DOUBLE PRECISION E
C local
      INTEGER XDIFF, YDIFF, ZDIFF
C parameter
      INTEGER ONE
      PARAMETER (ONE=1.0D0)
C begin
      XDIFF=ALLHP(IREAL8(NATOM))
      YDIFF=ALLHP(IREAL8(NATOM))
      ZDIFF=ALLHP(IREAL8(NATOM))
      CALL XRENE2(E,HEAP(HPFX),HEAP(HPFY),HEAP(HPFZ),
     &     HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),HEAP(XDIFF),
     &     HEAP(YDIFF),HEAP(ZDIFF),
     &     HEAP(HPATOF))
      CALL FREHP(ZDIFF,IREAL8(NATOM))
      CALL FREHP(YDIFF,IREAL8(NATOM))
      CALL FREHP(XDIFF,IREAL8(NATOM))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XRENE2(E,XRFX,XRFY,XRFZ,XRDX,XRDY,XRDZ,
     &           XDIFF,YDIFF,ZDIFF,XRATOF)
C
C See routine XRENER above
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      DOUBLE PRECISION E
      DOUBLE PRECISION XRFX(*), XRFY(*), XRFZ(*), XRDX(*)
      DOUBLE PRECISION XRDY(*), XRDZ(*), XDIFF(*), YDIFF(*)
      DOUBLE PRECISION ZDIFF(*)
      INTEGER XRATOF(*)
C local
      DOUBLE PRECISION DELABS
      INTEGER IAT
C parameters
      DOUBLE PRECISION ZERO, TWO
      PARAMETER (ZERO=0.0D0, TWO=2.0D0)
C begin
      E=ZERO
C
C check to make sure that atom selections are still defined
      IF (XRUPAT) THEN
      WRITE(6,'(A)')
     & ' %XRAY-ERR: atom numbers have changed. Information about',
     & '            atomic form factors and ASSOciate ',
     & '            atom selections has been lost.'
      CALL WRNDIE(-1,'XENER>','ASSOciate atom selections lost. ')
C
C check to make sure that there are some reflections
      ELSEIF (XRNREF.EQ.0) THEN
      WRITE(6,'(A)') ' XENER-err: number of reflections is zero.'
      CALL WRNDIE(-5,'XENER','number of reflections is zero. ')
C
C everything's OK
      ELSE
C
C The following section checks whether XCALCS has to be called, either
C because XRQCHK is true or because the difference between
C current and old coordinates is greater than XRTOL.
C The difference between current and old coordinates is returned
C in XDIFF, YDIFF, ZDIFF
C
C check whether XCALC has to be called
      IF (XRQCHK) THEN
      XRQCHK=.FALSE.
      CALL XCALCS(.TRUE.,.FALSE.,XRE,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRINTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           DRPNMX,DRPNX,DRPN,DRPNL,DRPNN,DRPNDB,
     &           DRPNMLT,DRPNLEV,DRPNTYP,DRPNDOM,DDEPTH,
     &           MRPNMX,MRPNX,MRPN,MRPNL,MRPNN,MRPNDB,
     &           MRPNMLT,MRPNLEV,MRPNTYP,MRPNDOM,MDEPTH,
     &           SRPNMX,SRPNX,SRPN,SRPNL,SRPNN,SRPNDB,
     &           SRPNMLT,SRPNLEV,SRPNTYP,SRPNDOM,SDEPTH,
     &           CRPNMX,CRPNX,CRPN,CRPNL,CRPNN,CRPNDB,
     &           CRPNMLT,CRPNLEV,CRPNTYP,CRPNDOM,CDEPTH,
     &           XRNREF,XRMREF,HPH,HPK,HPL,HPMULT,HPTYPE,
     &           HPTSEL,QHERM,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &           XRSYGP, XRSYIV,
     &           XRCELL,XRVOL,
     &           HPANOMFLAG,HPATOF,HPINDF,XNAMEAS,XASSOC,
     &           QASELE,XRNATF,QFFT,
     &           HPDX,HPDY,HPDZ,HPDT,HPDQ,
     &           XSFMX,XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           XSFGNAM,XSFGTYP,XSFGORD,
     &           XRSN,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,QLOOK,
     &           MAPR,IHPFLAG,HPFLAG,MAXHPFLAG,
     &           XRRED,XRREUP)
      DO IAT=1,XRNATF
      XDIFF(XRATOF(IAT))=ZERO
      YDIFF(XRATOF(IAT))=ZERO
      ZDIFF(XRATOF(IAT))=ZERO
      XRFX(XRATOF(IAT))=X(XRATOF(IAT))
      XRFY(XRATOF(IAT))=Y(XRATOF(IAT))
      XRFZ(XRATOF(IAT))=Z(XRATOF(IAT))
      END DO
      ELSE
C
C check difference between current and old coordinates
C if difference greater than tolerance call XCALCS
      DELABS=ZERO
      DO IAT=1,XRNATF
      XDIFF(XRATOF(IAT))=X(XRATOF(IAT))-XRFX(XRATOF(IAT))
      YDIFF(XRATOF(IAT))=Y(XRATOF(IAT))-XRFY(XRATOF(IAT))
      ZDIFF(XRATOF(IAT))=Z(XRATOF(IAT))-XRFZ(XRATOF(IAT))
      DELABS=
     &    MAX(DELABS,ABS(XDIFF(XRATOF(IAT))),
     &               ABS(YDIFF(XRATOF(IAT))),
     &               ABS(ZDIFF(XRATOF(IAT))))
      END DO
C
      IF (DELABS.GE.XRLTOL) THEN
      XRQCHK=.FALSE.
      CALL XCALCS(.TRUE.,.FALSE.,XRE,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRINTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           DRPNMX,DRPNX,DRPN,DRPNL,DRPNN,DRPNDB,
     &           DRPNMLT,DRPNLEV,DRPNTYP,DRPNDOM,DDEPTH,
     &           MRPNMX,MRPNX,MRPN,MRPNL,MRPNN,MRPNDB,
     &           MRPNMLT,MRPNLEV,MRPNTYP,MRPNDOM,MDEPTH,
     &           SRPNMX,SRPNX,SRPN,SRPNL,SRPNN,SRPNDB,
     &           SRPNMLT,SRPNLEV,SRPNTYP,SRPNDOM,SDEPTH,
     &           CRPNMX,CRPNX,CRPN,CRPNL,CRPNN,CRPNDB,
     &           CRPNMLT,CRPNLEV,CRPNTYP,CRPNDOM,CDEPTH,
     &           XRNREF,XRMREF,HPH,HPK,HPL,HPMULT,HPTYPE,
     &           HPTSEL,QHERM,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &           XRSYGP, XRSYIV,
     &           XRCELL,XRVOL,
     &           HPANOMFLAG,HPATOF,HPINDF,XNAMEAS,XASSOC,
     &           QASELE,XRNATF,QFFT,
     &           HPDX,HPDY,HPDZ,HPDT,HPDQ,
     &           XSFMX,XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           XSFGNAM,XSFGTYP,XSFGORD,
     &           XRSN,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,QLOOK,
     &           MAPR,IHPFLAG,HPFLAG,MAXHPFLAG,
     &           XRRED,XRREUP)
      DO IAT=1,XRNATF
      XDIFF(XRATOF(IAT))=ZERO
      YDIFF(XRATOF(IAT))=ZERO
      ZDIFF(XRATOF(IAT))=ZERO
      XRFX(XRATOF(IAT))=X(XRATOF(IAT))
      XRFY(XRATOF(IAT))=Y(XRATOF(IAT))
      XRFZ(XRATOF(IAT))=Z(XRATOF(IAT))
      END DO
      END IF
      END IF
C
C compute the energy
      E=XRE
C
      DO IAT=1,XRNATF
C
C accumulate energy
      E=E +XRDX(XRATOF(IAT))*XDIFF(XRATOF(IAT))
     &    +XRDY(XRATOF(IAT))*YDIFF(XRATOF(IAT))
     &    +XRDZ(XRATOF(IAT))*ZDIFF(XRATOF(IAT))
C
C compute and accumulate derivatives
      DX(XRATOF(IAT))=DX(XRATOF(IAT)) + XRDX(XRATOF(IAT))
      DY(XRATOF(IAT))=DY(XRATOF(IAT)) + XRDY(XRATOF(IAT))
      DZ(XRATOF(IAT))=DZ(XRATOF(IAT)) + XRDZ(XRATOF(IAT))
C
      END DO
C
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE XCALCS(QDXYZ,QDB,XRE,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRINTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           DRPNMX,DRPNX,DRPN,DRPNL,DRPNN,DRPNDB,
     &           DRPNMLT,DRPNLEV,DRPNTYP,DRPNDOM,DDEPTH,
     &           MRPNMX,MRPNX,MRPN,MRPNL,MRPNN,MRPNDB,
     &           MRPNMLT,MRPNLEV,MRPNTYP,MRPNDOM,MDEPTH,
     &           SRPNMX,SRPNX,SRPN,SRPNL,SRPNN,SRPNDB,
     &           SRPNMLT,SRPNLEV,SRPNTYP,SRPNDOM,SDEPTH,
     &           CRPNMX,CRPNX,CRPN,CRPNL,CRPNN,CRPNDB,
     &           CRPNMLT,CRPNLEV,CRPNTYP,CRPNDOM,CDEPTH,
     &           XRNREF,XRMREF,HPH,HPK,HPL,HPMULT,HPTYPE,
     &           HPTSEL,QHERM,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &           XRSYGP, XRSYIV,
     &           XRCELL,XRVOL,
     &           HPANOMFLAG,HPATOF,HPINDF,XNAMEAS,XASSOC,
     &           QASELE,XRNATF,QFFT,
     &           HPDX,HPDY,HPDZ,HPDT,HPDQ,
     &           XSFMX,XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           XSFGNAM,XSFGTYP,XSFGORD,
     &           XRSN,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,QLOOK,
     &           MAPR,IHPFLAG,HPFLAG,MAXHPFLAG,
     &           XRRED,XRREUP)
C
C Routine computes associated structure factor objects
C from specified atom selections, target function,
C and derivatives (only if QDXYZ or QDB are set to true).
C
C If QDXYZ is true, the derivatives are returned in
C XRDX, XRDY, XRDY.
C If QDB is true, the derivatives with respect to b, q, f'
C and f'' are returned in XRDT, XRDQ, XRDX, XRDY, respectively.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      LOGICAL QDXYZ, QDB
      DOUBLE PRECISION XRE
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3)
      DOUBLE PRECISION XRSCAL
      LOGICAL XCVTEST
C
      INTEGER TRPNMX, TRPNX
      CHARACTER*(*) TRPN(4,TRPNMX)
      INTEGER TRPNL(4,TRPNMX), TRPNN
      DOUBLE COMPLEX TRPNDB(4,TRPNMX)
      INTEGER TRPNMLT(TRPNMX), TRPNLEV(TRPNMX)
      CHARACTER*2 TRPNTYP(TRPNMX), TRPNDOM(TRPNMX)
      INTEGER TDEPTH
C
      INTEGER DRPNMX, DRPNX
      CHARACTER*(*) DRPN(4,DRPNMX,*)
      INTEGER DRPNL(4,DRPNMX,*), DRPNN(*)
      DOUBLE COMPLEX DRPNDB(4,DRPNMX,*)
      INTEGER DRPNMLT(DRPNMX,*), DRPNLEV(DRPNMX,*)
      CHARACTER*2 DRPNTYP(DRPNMX,*), DRPNDOM(DRPNMX,*)
      INTEGER DDEPTH(*)
C
      INTEGER MRPNMX, MRPNX
      CHARACTER*(*) MRPN(4,MRPNMX)
      INTEGER MRPNL(4,MRPNMX), MRPNN
      DOUBLE COMPLEX MRPNDB(4,MRPNMX)
      INTEGER MRPNMLT(MRPNMX), MRPNLEV(MRPNMX)
      CHARACTER*2 MRPNTYP(MRPNMX), MRPNDOM(MRPNMX)
      INTEGER MDEPTH
C
      INTEGER SRPNMX, SRPNX
      CHARACTER*(*) SRPN(4,SRPNMX)
      INTEGER SRPNL(4,SRPNMX), SRPNN
      DOUBLE COMPLEX SRPNDB(4,SRPNMX)
      INTEGER SRPNMLT(SRPNMX), SRPNLEV(SRPNMX)
      CHARACTER*2 SRPNTYP(SRPNMX), SRPNDOM(SRPNMX)
      INTEGER SDEPTH
C
      INTEGER CRPNMX, CRPNX
      CHARACTER*(*) CRPN(4,CRPNMX)
      INTEGER CRPNL(4,CRPNMX), CRPNN
      DOUBLE COMPLEX CRPNDB(4,CRPNMX)
      INTEGER CRPNMLT(CRPNMX), CRPNLEV(CRPNMX)
      CHARACTER*2 CRPNTYP(CRPNMX), CRPNDOM(CRPNMX)
      INTEGER CDEPTH
C
      INTEGER XRNREF, XRMREF, HPH, HPK, HPL, HPMULT, HPTYPE
      INTEGER HPTSEL
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
      DOUBLE PRECISION XRCELL(3,3), XRVOL
      INTEGER HPANOMFLAG, HPATOF, HPINDF, XNAMEAS
      CHARACTER*(*) XASSOC(*)
      LOGICAL QASELE
      INTEGER XRNATF
      LOGICAL QFFT
      INTEGER HPDX, HPDY, HPDZ, HPDT, HPDQ
      INTEGER XSFMX, XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER XRSN, XRSM
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(*), XRFDP(*)
      LOGICAL QLOOK
      DOUBLE PRECISION MAPR
      INTEGER IHPFLAG, HPFLAG(*), MAXHPFLAG
      LOGICAL XRRED, XRREUP
C local
C pointer
      INTEGER XDERIV, XRNATO, NANOM, IFCALC, HPATOM, HPINDX
      INTEGER PFCALC
C begin
C
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
C get target selection
      CALL XTSELSET(HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,
     &           HPTYPE,XRREUP,HPTSEL,
     &           XRMREF,XRNREF,XCVTEST,
     &           XRTR,
     &           SRPNMX,SRPNX,SRPN,SRPNL,SRPNN,SRPNDB,
     &           SRPNMLT,SRPNLEV,SRPNTYP,SRPNDOM,SDEPTH,
     &           CRPNMX,CRPNX,CRPN,CRPNL,CRPNN,CRPNDB,
     &           CRPNMLT,CRPNLEV,CRPNTYP,CRPNDOM,CDEPTH,
     &           QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRCELL,XRVOL)
C
C
      IF (IHPFLAG.EQ.0) THEN
      CALL WRNDIE(-5,'XCALCS',
     &     'ASSOciate object(s) undefined.')
      ELSE
      XDERIV=ALLHP(ICPLX8(XRNREF*IHPFLAG))
      XRNATO=ALLHP(INTEG4(MAXHPFLAG))
      NANOM=ALLHP(INTEG4(MAXHPFLAG))
      IFCALC=ALLHP(INTEG4(MAXHPFLAG))
      HPATOM=ALLHP(INTEG4(MAXHPFLAG))
      HPINDX=ALLHP(INTEG4(MAXHPFLAG))
      PFCALC=ALLHP(INTEG4(MAXHPFLAG))
      CALL XCALCS2(QDXYZ,QDB,HEAP(XDERIV),XRE,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRINTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           DRPNMX,DRPNX,DRPN,DRPNL,DRPNN,DRPNDB,
     &           DRPNMLT,DRPNLEV,DRPNTYP,DRPNDOM,DDEPTH,
     &           MRPNMX,MRPNX,MRPN,MRPNL,MRPNN,MRPNDB,
     &           MRPNMLT,MRPNLEV,MRPNTYP,MRPNDOM,MDEPTH,
     &           XRNREF,XRMREF,HPH,HPK,HPL,HPMULT,HPTYPE,
     &           HPTSEL,QHERM,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &           XRCELL,XRVOL,
     &           HPANOMFLAG,HPATOF,HPINDF,XNAMEAS,XASSOC,
     &           QASELE,XRNATF,QFFT,
     &           HPDX,HPDY,HPDZ,HPDT,HPDQ,
     &           XSFMX,XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           XSFGNAM,XSFGTYP,XSFGORD,
     &           XRSN,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,QLOOK,
     &           MAPR,IHPFLAG,HPFLAG,
     &           HEAP(XRNATO),HEAP(NANOM),HEAP(IFCALC),
     &           HEAP(HPATOM),HEAP(HPINDX),
     &           HEAP(PFCALC))
      CALL FREHP(XDERIV,ICPLX8(XRNREF*IHPFLAG))
      CALL FREHP(XRNATO,INTEG4(MAXHPFLAG))
      CALL FREHP(NANOM,INTEG4(MAXHPFLAG))
      CALL FREHP(IFCALC,INTEG4(MAXHPFLAG))
      CALL FREHP(HPATOM,INTEG4(MAXHPFLAG))
      CALL FREHP(HPINDX,INTEG4(MAXHPFLAG))
      CALL FREHP(PFCALC,INTEG4(MAXHPFLAG))
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE XCALCS2(QDXYZ,QDB,XDERIV,XRE,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRINTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           DRPNMX,DRPNX,DRPN,DRPNL,DRPNN,DRPNDB,
     &           DRPNMLT,DRPNLEV,DRPNTYP,DRPNDOM,DDEPTH,
     &           MRPNMX,MRPNX,MRPN,MRPNL,MRPNN,MRPNDB,
     &           MRPNMLT,MRPNLEV,MRPNTYP,MRPNDOM,MDEPTH,
     &           XRNREF,XRMREF,HPH,HPK,HPL,HPMULT,HPTYPE,
     &           HPTSEL,QHERM,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &           XRCELL,XRVOL,
     &           HPANOMFLAG,HPATOF,HPINDF,XNAMEAS,XASSOC,
     &           QASELE,XRNATF,QFFT,
     &           HPDX,HPDY,HPDZ,HPDT,HPDQ,
     &           XSFMX,XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           XSFGNAM,XSFGTYP,XSFGORD,
     &           XRSN,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,QLOOK,
     &           MAPR,IHPFLAG,HPFLAG,
     &           XRNATO, NANOM, IFCALC,HPATOM,HPINDX,
     &           PFCALC)
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      LOGICAL QDXYZ, QDB
      INTEGER XRNREF
      DOUBLE COMPLEX XDERIV(XRNREF,*)
      DOUBLE PRECISION XRE
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3)
      DOUBLE PRECISION XRSCAL
      LOGICAL XCVTEST
C
      INTEGER TRPNMX, TRPNX
      CHARACTER*(*) TRPN(4,TRPNMX)
      INTEGER TRPNL(4,TRPNMX), TRPNN
      DOUBLE COMPLEX TRPNDB(4,TRPNMX)
      INTEGER TRPNMLT(TRPNMX), TRPNLEV(TRPNMX)
      CHARACTER*2 TRPNTYP(TRPNMX), TRPNDOM(TRPNMX)
      INTEGER TDEPTH
C
      INTEGER DRPNMX, DRPNX
      CHARACTER*(*) DRPN(4,DRPNMX,*)
      INTEGER DRPNL(4,DRPNMX,*), DRPNN(*)
      DOUBLE COMPLEX DRPNDB(4,DRPNMX,*)
      INTEGER DRPNMLT(DRPNMX,*), DRPNLEV(DRPNMX,*)
      CHARACTER*2 DRPNTYP(DRPNMX,*), DRPNDOM(DRPNMX,*)
      INTEGER DDEPTH(*)
C
      INTEGER MRPNMX, MRPNX
      CHARACTER*(*) MRPN(4,MRPNMX)
      INTEGER MRPNL(4,MRPNMX), MRPNN
      DOUBLE COMPLEX MRPNDB(4,MRPNMX)
      INTEGER MRPNMLT(MRPNMX), MRPNLEV(MRPNMX)
      CHARACTER*2 MRPNTYP(MRPNMX), MRPNDOM(MRPNMX)
      INTEGER MDEPTH
C
      INTEGER XRMREF, HPH, HPK, HPL, HPMULT, HPTYPE
      INTEGER HPTSEL
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRCELL(3,3), XRVOL
      INTEGER HPANOMFLAG, HPATOF, HPINDF, XNAMEAS
      CHARACTER*(*) XASSOC(*)
      LOGICAL QASELE
      INTEGER XRNATF
      LOGICAL QFFT
      INTEGER HPDX, HPDY, HPDZ, HPDT, HPDQ
      INTEGER XSFMX, XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER XRSN, XRSM
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(*), XRFDP(*)
      LOGICAL QLOOK
      INTEGER IHPFLAG, HPFLAG(*)
      DOUBLE PRECISION MAPR
      INTEGER XRNATO(*), NANOM(*), IFCALC(*)
      INTEGER HPATOM(*), HPINDX(*)
      INTEGER PFCALC(*)
C local
      INTEGER MAXHPFLAG
      PARAMETER (MAXHPFLAG=10)
      LOGICAL QDIRCT
      INTEGER II
      INTEGER FNDIOBJ
      EXTERNAL FNDIOBJ
      DOUBLE COMPLEX DBCOMP
C pointers
      INTEGER HPFQS
      INTEGER XRNCSX, XRNCSY, XRNCSZ
      INTEGER XNCS, YNCS, ZNCS, DXNCS, DYNCS, DZNCS, EXPTBL
C parameters
      INTEGER NTBL
      PARAMETER (NTBL=4000)
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
      IF (QDXYZ.AND.QDB) THEN
      CALL WRNDIE(-5,'XCALCS',
     &     'programming error: QDB and QDXYZ true')
      END IF
C
C
C allocate space for temporary atom flag arrays
      DO II=1,IHPFLAG
      HPATOM(II)=ALLHP(INTEG4(NATOM))
      HPINDX(II)=ALLHP(INTEG4(NATOM))
      END DO
      HPFQS=ALLHP(IREAL8(NATOM))
C
      DO II=1,IHPFLAG
C get atom selections from associate statements
      CALL XRASSOC(HEAP(HPFLAG(II)),HEAP(HPANOMFLAG),
     &            HEAP(HPATOF),HEAP(HPINDF),HEAP(HPATOM(II)),
     &            HEAP(HPINDX(II)),HEAP(HPFQS),XRNATO(II),NANOM(II),
     &            XNAMEAS,XASSOC(II),QASELE,XRNATF,XRFDP,QHERM,
     &            XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &            XRTR,XRINTR)
      END DO
C
C decide whether direct summation method is required
      QDIRCT=.NOT.QFFT
      DO II=1,IHPFLAG
      QDIRCT=QDIRCT.OR.(XRNATO(II).GT.NANOM(II))
      END DO
C
C allocate space for the FFT routines if required
      IF (QFFT) THEN
      XNCS=ALLHP(IREAL8(NATOM))
      YNCS=ALLHP(IREAL8(NATOM))
      ZNCS=ALLHP(IREAL8(NATOM))
      DXNCS=ALLHP(IREAL8(NATOM))
      DYNCS=ALLHP(IREAL8(NATOM))
      DZNCS=ALLHP(IREAL8(NATOM))
      EXPTBL=ALLHP(IREAL8(NTBL+1))
      END IF
C
C allocate space for the direct summation method if required
      IF (QDIRCT) THEN
      XRNCSX=ALLHP(IREAL8(NATOM))
      XRNCSY=ALLHP(IREAL8(NATOM))
      XRNCSZ=ALLHP(IREAL8(NATOM))
      END IF
C
C initialize derivative arrays
      IF (QDXYZ.OR.QDB) THEN
      CALL FILLR8(HEAP(HPDX),NATOM,ZERO)
      CALL FILLR8(HEAP(HPDY),NATOM,ZERO)
      CALL FILLR8(HEAP(HPDZ),NATOM,ZERO)
      CALL FILLR8(HEAP(HPDT),NATOM,ZERO)
      CALL FILLR8(HEAP(HPDQ),NATOM,ZERO)
      END IF
C
C loop through associate statements
      DO II=1,IHPFLAG
C
C get indices of reciprocal space objects corresponding to each
C associate statement.  If non-existent, declare a structure factor object.
      PFCALC(II)=0
      IFCALC(II)=FNDIOBJ(XSFNUM, XSFNAM, XASSOC(II))
      IF (IFCALC(II).LE.0) THEN
C declare object
      IF (XSFNUM.GE.XSFMX) THEN
      CALL WRNDIE(-5,'XCALCS',
     & 'exceeded XSFMX parameter --> recompile program')
      ELSE
      XSFNUM=XSFNUM+1
      XSFNAM(XSFNUM)=XASSOC(II)
      XSFTYPE(XSFNUM)='COMPLEX'
      XSFGNAM(XSFNUM)=0
      XSFGTYP(XSFNUM)=' '
      XSFGORD(XSFNUM)=0
      HPSF(XSFNUM)=0
      IFCALC(II)=XSFNUM
      END IF
      END IF
      CALL CHKSFTYPE(IFCALC(II), XSFNAM, XSFTYPE, 'COMP', 'XCALCS')
C
      IF (IFCALC(II).GT.0.AND.XRNREF.GT.0) THEN
C
C allocate space if necessary
      IF (HPSF(IFCALC(II)).EQ.0) THEN
      CALL XSFAL(HPSF(IFCALC(II)), XRMREF, XSFTYPE(IFCALC(II)))
      END IF
C
C define PFCALC pointer
      PFCALC(II)=HPSF(IFCALC(II))
C
C initialize FCALC
      CALL FILLC8(HEAP(PFCALC(II)),XRNREF,DCMPLX(ZERO,ZERO))
      END IF
      END DO
C
C compute ASSOciate objects from selected atom sets
C
C loop through associate statements
      DO II=1,IHPFLAG
C
      IF (XRNATO(II).GT.0.AND.PFCALC(II).NE.0) THEN
C
C normal scatterers: call FFT method or direct summation
      IF (QFFT.AND.XRNATO(II).GT.0) THEN
      CALL XFFTUP(.FALSE.,.FALSE.,1,XRNATO(II),
     &     HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(PFCALC(II)),
     &     HEAP(HPATOM(II)),HEAP(HPINDX(II)),HEAP(HPFQS),
     &     HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),HEAP(HPDT),
     &     HEAP(HPDQ),HEAP(XNCS),HEAP(YNCS),HEAP(ZNCS),
     &     HEAP(DXNCS),HEAP(DYNCS),HEAP(DZNCS),
     &     NTBL,HEAP(EXPTBL),QLOOK,XRSM,XRSA,XRSB,XRSC,XRFP,
     &     HEAP(HPTSEL),XRNREF,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &     QHERM,XRCELL,MAPR,XRTR,XRINTR,XRVOL)
C
C compute anomalous contribution: use direct summation
      IF (NANOM(II)+1.LE.XRNATO(II)) THEN
      CALL XRDRCT(HEAP(HPTSEL),.FALSE.,.FALSE.,NANOM(II)+1,XRNATO(II),
     &     HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(PFCALC(II)),
     &     HEAP(HPATOM(II)),HEAP(HPINDX(II)),HEAP(HPFQS),
     &     HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),HEAP(HPDT),HEAP(HPDQ),
     &     HEAP(XRNCSX),HEAP(XRNCSY),HEAP(XRNCSZ),
     &     .TRUE.,.TRUE.,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRNREF)
      END IF
      ELSEIF (.NOT.QFFT) THEN
      IF (NANOM(II).GT.0) THEN
      CALL XRDRCT(HEAP(HPTSEL),.FALSE.,.FALSE.,1,NANOM(II),
     &     HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(PFCALC(II)),
     &     HEAP(HPATOM(II)),HEAP(HPINDX(II)),HEAP(HPFQS),
     &     HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),HEAP(HPDT),HEAP(HPDQ),
     &     HEAP(XRNCSX),HEAP(XRNCSY),HEAP(XRNCSZ),
     &     .FALSE.,.FALSE.,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRNREF)
      END IF
C
C compute anomalous contribution: use direct summation
      IF (NANOM(II)+1.LE.XRNATO(II)) THEN
      CALL XRDRCT(HEAP(HPTSEL),.FALSE.,.FALSE.,NANOM(II)+1,XRNATO(II),
     &     HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(PFCALC(II)),
     &     HEAP(HPATOM(II)),HEAP(HPINDX(II)),HEAP(HPFQS),
     &     HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),HEAP(HPDT),HEAP(HPDQ),
     &     HEAP(XRNCSX),HEAP(XRNCSY),HEAP(XRNCSZ),
     &     .TRUE.,.FALSE.,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRNREF)
      END IF
      END IF
      END IF
      END DO
C
C use the reciprocal space targets and derivatives
C (dtarget/d<associate_object) if required.
C
C for cross-validation print the target value for the test set
C (we must not compute the derivatives for the test set!)
      IF (XCVTEST) THEN
C
C test set;  don't compute derivatives
      CALL XTARGETS(.FALSE.,.TRUE.,-1,XRE,XDERIV,
     &           HPTSEL,XRNREF,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           1,1,' ',1,0,DBCOMP,1,1,' ',' ',1,
     &           MRPNMX,MRPNX,MRPN,MRPNL,MRPNN,MRPNDB,
     &           MRPNMLT,MRPNLEV,MRPNTYP,MRPNDOM,MDEPTH,
     &           HPH,HPK,HPL,XSFNUM,XSFNAM,
     &           XSFTYPE,HPSF,HPMULT,HPTYPE,XRMREF,
     &           QHERM,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,0,
     &           XRCELL,XRVOL)
      END IF
C
C working set;  compute derivatives. They'll come back
C in XDERIV(*,II)
      CALL XTARGETS(QDXYZ.OR.QDB,.TRUE.,1,XRE,XDERIV,
     &           HPTSEL,XRNREF,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           DRPNMX,DRPNX,DRPN,DRPNL,DRPNN,DRPNDB,
     &           DRPNMLT,DRPNLEV,DRPNTYP,DRPNDOM,DDEPTH,
     &           MRPNMX,MRPNX,MRPN,MRPNL,MRPNN,MRPNDB,
     &           MRPNMLT,MRPNLEV,MRPNTYP,MRPNDOM,MDEPTH,
     &           HPH,HPK,HPL,XSFNUM,XSFNAM,
     &           XSFTYPE,HPSF,HPMULT,HPTYPE,XRMREF,
     &           QHERM,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,IHPFLAG,
     &           XRCELL,XRVOL)
C
      IF (QDXYZ.OR.QDB) THEN
C
C compute derivatives with respect to atomic parameters
C
C loop through associate statements
      DO II=1,IHPFLAG
      IF (XRNATO(II).GT.0.AND.PFCALC(II).NE.0) THEN
C All scatterers: call FFT method or direct summation
      IF (QFFT) THEN
C
C note that XFFTUP modifies the XDERIV(*,II) array, so we
C have to call this routine last!!!
C
C compute anomalous contribution: use direct summation
      IF (NANOM(II)+1.LE.XRNATO(II)) THEN
      CALL XRDRCT(HEAP(HPTSEL),QDXYZ,QDB,NANOM(II)+1,XRNATO(II),
     &     HEAP(HPH),HEAP(HPK),HEAP(HPL),XDERIV(1,II),
     &     HEAP(HPATOM(II)),HEAP(HPINDX(II)),HEAP(HPFQS),
     &     HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),HEAP(HPDT),HEAP(HPDQ),
     &     HEAP(XRNCSX),HEAP(XRNCSY),HEAP(XRNCSZ),
     &     .TRUE.,.TRUE.,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRNREF)
      END IF
      CALL XFFTUP(QDXYZ,QDB,1,XRNATO(II),
     &     HEAP(HPH),HEAP(HPK),HEAP(HPL),XDERIV(1,II),
     &     HEAP(HPATOM(II)),HEAP(HPINDX(II)),HEAP(HPFQS),
     &     HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),HEAP(HPDT),
     &     HEAP(HPDQ),HEAP(XNCS),HEAP(YNCS),HEAP(ZNCS),
     &     HEAP(DXNCS),HEAP(DYNCS),HEAP(DZNCS),
     &     NTBL,HEAP(EXPTBL),QLOOK,XRSM,XRSA,XRSB,XRSC,XRFP,
     &     HEAP(HPTSEL),XRNREF,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &     QHERM,XRCELL,MAPR,XRTR,XRINTR,XRVOL)
      ELSEIF (.NOT.QFFT) THEN
C
C compute anomalous contribution: use direct summation
      IF (NANOM(II)+1.LE.XRNATO(II)) THEN
      CALL XRDRCT(HEAP(HPTSEL),QDXYZ,QDB,NANOM(II)+1,XRNATO(II),
     &     HEAP(HPH),HEAP(HPK),HEAP(HPL),XDERIV(1,II),
     &     HEAP(HPATOM(II)),HEAP(HPINDX(II)),HEAP(HPFQS),
     &     HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),HEAP(HPDT),HEAP(HPDQ),
     &     HEAP(XRNCSX),HEAP(XRNCSY),HEAP(XRNCSZ),
     &     .TRUE.,.FALSE.,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRNREF)
      END IF
      IF (NANOM(II).GT.0) THEN
      CALL XRDRCT(HEAP(HPTSEL),QDXYZ,QDB,1,NANOM(II),
     &     HEAP(HPH),HEAP(HPK),HEAP(HPL),XDERIV(1,II),
     &     HEAP(HPATOM(II)),HEAP(HPINDX(II)),HEAP(HPFQS),
     &     HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),HEAP(HPDT),HEAP(HPDQ),
     &     HEAP(XRNCSX),HEAP(XRNCSY),HEAP(XRNCSZ),
     &     .FALSE.,.FALSE.,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRNREF)
      END IF
      END IF
      END IF
      END DO
C
C
      END IF
C
C Free up space
      IF (QFFT) THEN
      CALL FREHP(EXPTBL,IREAL8(NTBL+1))
      CALL FREHP(DZNCS,IREAL8(NATOM))
      CALL FREHP(DXNCS,IREAL8(NATOM))
      CALL FREHP(DYNCS,IREAL8(NATOM))
      CALL FREHP(ZNCS,IREAL8(NATOM))
      CALL FREHP(XNCS,IREAL8(NATOM))
      CALL FREHP(YNCS,IREAL8(NATOM))
      END IF
C
C free up heap space for direct summation method
      IF (QDIRCT) THEN
      CALL FREHP(XRNCSZ,IREAL8(NATOM))
      CALL FREHP(XRNCSY,IREAL8(NATOM))
      CALL FREHP(XRNCSX,IREAL8(NATOM))
      END IF
C
C
      DO II=1,IHPFLAG
      CALL FREHP(HPATOM(II),INTEG4(NATOM))
      CALL FREHP(HPINDX(II),INTEG4(NATOM))
      END DO
      CALL FREHP(HPFQS,IREAL8(NATOM))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XRDRCT(TSEL,QDXYZ,QDB,START,STOP,
     &            XRH,XRK,XRL,FCALC,XRATOM,
     &            XRINDX,XRFQS,XRDX,XRDY,
     &            XRDZ,XRDT,XRDQ,XRNCSX,XRNCSY,XRNCSZ,
     &            QANOM,QIMAG,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,
     &            XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRNREF)
C
C Computes the structure factor equation by direct summation.
C
C Mode:
C    QDXYZ and QDB false: compute fcalcs.
C    QDXYZ or QDB true: compute derivatives with respect to
C                       x,y,z (QDXYZ) or to b,q,f',f'' (QDB).
C                       The QDXYZ and QDB options are mutually
C                       exclusive.  On input, the routine assumes that
C                       the derivatives of the target function
C                       are stored in FCALC.
C
C All atoms between START and STOP are used in the summation.
C
C QIMAG flags the calculation of only the imaginary (f'') part
C of the form factors.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'mtf.inc'
      INTEGER TSEL(*)
      LOGICAL QDXYZ, QDB
      INTEGER START, STOP, XRH(*), XRK(*), XRL(*)
      DOUBLE COMPLEX FCALC(*)
      INTEGER XRATOM(*), XRINDX(*)
      DOUBLE PRECISION XRFQS(*)
      DOUBLE PRECISION XRDX(NATOM), XRDY(NATOM), XRDZ(NATOM)
      DOUBLE PRECISION XRDT(NATOM), XRDQ(NATOM)
      DOUBLE PRECISION XRNCSX(*), XRNCSY(*), XRNCSZ(*)
      LOGICAL QANOM, QIMAG
      INTEGER XRSM
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(XRSM), XRFDP(XRSM)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRNREF
C local
      INTEGER REFLCT, ISYM, INSYM, IAT, B, C, U, V
      DOUBLE PRECISION RTH, SYMM(3,4), FA, FB, H2PI, K2PI, L2PI
      DOUBLE PRECISION ARG, SSQ, SSQ4, DX, DY, DZ
      DOUBLE PRECISION DXA, DYA, DZA, DXB, DYB, DZB, DTA, DTB, DQA, DQB
      DOUBLE COMPLEX CTEMP
      DOUBLE PRECISION XRFFL
      DOUBLE COMPLEX FLOCF, FLOCFQ, FLOCFP, FLOCFDP, CARG
      DOUBLE PRECISION EXPB, XRFDPX
      DOUBLE PRECISION SECS, SECS2
C parameters
      DOUBLE PRECISION ZERO, ONE, TWOPI, PI1000, FOUR
      PARAMETER (ZERO=0.0D0, TWOPI=2.0D0*PI, PI1000=10000.0D0*PI)
      PARAMETER (FOUR=4.0D0, ONE=1.0D0)
C begin
C
      IF (TIMER.GT.0) CALL VCPU(SECS)
C loop over non-crystallographic symmetry operators
C
      DO INSYM=1,XNNSYM
C
C loop over crystallographic symmetry operators
      DO ISYM=1,XRNSYM
C
C Concatenate non-crystallographic symmetry operator, orthogonal
C to fractional operator and the crystallographic symmetry operator.
C The result will be in fractional coordintes.
C   SYMM=[XRSYMM(ISYM,*,*)*XRTR] *NCSOP(INSYM,*,*)
      RTH=XRSYTH
      DO U=1,3
      DO V=1,4
      SYMM(U,V)=ZERO
      DO B=1,3
      DO C=1,3
      SYMM(U,V)=SYMM(U,V)
     &         +XRSYMM(ISYM,U,B)*XRTR(B,C)*NCSOP(INSYM,C,V)
      END DO
      END DO
      END DO
      SYMM(U,4)=SYMM(U,4)
     &         +XRSYMM(ISYM,U,4)/RTH
      END DO
C
C apply the concatenated operators to the atomic coordinates.
C store the result in XRNCSX, XRNCSY, XRNCSZ
      DO IAT=START,STOP
      XRNCSX(IAT)=SYMM(1,1)*X(XRATOM(IAT)) + SYMM(1,2)*Y(XRATOM(IAT))
     &          + SYMM(1,3)*Z(XRATOM(IAT)) + SYMM(1,4)
      XRNCSY(IAT)=SYMM(2,1)*X(XRATOM(IAT)) + SYMM(2,2)*Y(XRATOM(IAT))
     &          + SYMM(2,3)*Z(XRATOM(IAT)) + SYMM(2,4)
      XRNCSZ(IAT)=SYMM(3,1)*X(XRATOM(IAT)) + SYMM(3,2)*Y(XRATOM(IAT))
     &          + SYMM(3,3)*Z(XRATOM(IAT)) + SYMM(3,4)
      END DO
C
C
C branch into appropriate part of the routine
      IF (.NOT.QDXYZ.AND..NOT.QDB) THEN
C
C ************************
C structure factor section
C ========================
C
C loop over all reflections
!$omp parallel do default(shared) 
!$omp& schedule(dynamic)
!$omp& private(reflct,h2pi,k2pi,l2pi,ssq,ssq4)
!$omp& private(iat,xrffl,xrfdpx,flocf,arg)
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).NE.0) THEN
C
C compute H*2*PI, K*2*PI, L*2*PI for this reflection
      H2PI=TWOPI*XRH(REFLCT)
      K2PI=TWOPI*XRK(REFLCT)
      L2PI=TWOPI*XRL(REFLCT)
C
C compute s**2 for this reflection
      CALL XRSSQ(XRH(REFLCT),XRK(REFLCT),XRL(REFLCT),SSQ,XRTR)
C
C compute s**2/4 for this reflection
      SSQ4=SSQ/FOUR
C
      DO IAT=START,STOP
C
C compute atomic from factors
C not needed if only imaginary component
      IF (.NOT.QIMAG) THEN
      XRFFL=XRSC(XRINDX(IAT))+XRFP(XRINDX(IAT))
     &     +XRSA(XRINDX(IAT),1) *EXP( -XRSB(XRINDX(IAT),1)*SSQ4 )
     &     +XRSA(XRINDX(IAT),2) *EXP( -XRSB(XRINDX(IAT),2)*SSQ4 )
     &     +XRSA(XRINDX(IAT),3) *EXP( -XRSB(XRINDX(IAT),3)*SSQ4 )
     &     +XRSA(XRINDX(IAT),4) *EXP( -XRSB(XRINDX(IAT),4)*SSQ4 )
      ELSE
      XRFFL=ZERO
      END IF
C
C in anomalous mode we include the f'' contribution to the SF
      IF (QANOM) THEN
      XRFDPX=XRFDP(XRINDX(IAT))
      ELSE
      XRFDPX=ZERO
      END IF
C
      FLOCF=XRFQS(XRATOM(IAT))*QMAIN(XRATOM(IAT))
     &           *EXP(-WMAIN(XRATOM(IAT))*SSQ4)
     &          *DCMPLX(XRFFL,XRFDPX)
C
      ARG=H2PI*XRNCSX(IAT)+K2PI*XRNCSY(IAT)+L2PI*XRNCSZ(IAT)
      FCALC(REFLCT)=FCALC(REFLCT)
     &                 +FLOCF*DCMPLX(COS(ARG),SIN(ARG))
      END DO
      END IF
      END DO
C
      ELSEIF (QDXYZ.OR.QDB) THEN
C
C ******************
C derivative section
C ==================
C
C
C loop over all reflections
!$omp parallel do default(shared) 
!$omp& schedule(dynamic)
!$omp& private(reflct,h2pi,k2pi,l2pi,ssq,ssq4)
!$omp& private(iat,xrffl,expb,xrfdpx,arg,carg,flocf)
!$omp$ private(ctemp,fa,fb)
!$omp& private(dx,dy,dz,dxa,dya,dza,dxb,dyb,dzb)
!$omp& private(dta,dtb,dqa,dqb,flocfq,flocfp,flocfdp)
!$omp& reduction(+:xrdx,xrdy,xrdz,xrdq,xrdt)
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).NE.0) THEN
C
C compute H*2*PI, K*2*PI, L*2*PI for this reflection
      H2PI=TWOPI*XRH(REFLCT)
      K2PI=TWOPI*XRK(REFLCT)
      L2PI=TWOPI*XRL(REFLCT)
C
C compute s**2 for this reflection
      CALL XRSSQ(XRH(REFLCT),XRK(REFLCT),XRL(REFLCT),SSQ,XRTR)
C
C compute s**2/4 for this reflection
      SSQ4=SSQ/FOUR
C
C compute derivative of symmetry operation for this reflection
C (apply transpose of symmetry operator index)
      DX= SYMM(1,1)*H2PI + SYMM(2,1)*K2PI + SYMM(3,1)*L2PI
      DY= SYMM(1,2)*H2PI + SYMM(2,2)*K2PI + SYMM(3,2)*L2PI
      DZ= SYMM(1,3)*H2PI + SYMM(2,3)*K2PI + SYMM(3,3)*L2PI
C
C multiply by derivatives of target w/respect to real and imag FCALC
      DXA=DX*DBLE(FCALC(REFLCT))
      DYA=DY*DBLE(FCALC(REFLCT))
      DZA=DZ*DBLE(FCALC(REFLCT))
      DXB=DX*DIMAG(FCALC(REFLCT))
      DYB=DY*DIMAG(FCALC(REFLCT))
      DZB=DZ*DIMAG(FCALC(REFLCT))
      DTA=-DBLE(FCALC(REFLCT))*SSQ4
      DTB=-DIMAG(FCALC(REFLCT))*SSQ4
      DQA=DBLE(FCALC(REFLCT))
      DQB=DIMAG(FCALC(REFLCT))
C
C accumulate the derivatives with respect to atomic coordinates
      DO IAT=START,STOP
C
C compute atomic from factors
C not needed if only imaginary component
      IF (.NOT.QIMAG) THEN
      XRFFL=XRSC(XRINDX(IAT))+XRFP(XRINDX(IAT))
     &     +XRSA(XRINDX(IAT),1) *EXP( -XRSB(XRINDX(IAT),1)*SSQ4 )
     &     +XRSA(XRINDX(IAT),2) *EXP( -XRSB(XRINDX(IAT),2)*SSQ4 )
     &     +XRSA(XRINDX(IAT),3) *EXP( -XRSB(XRINDX(IAT),3)*SSQ4 )
     &     +XRSA(XRINDX(IAT),4) *EXP( -XRSB(XRINDX(IAT),4)*SSQ4 )
      ELSE
      XRFFL=ZERO
      END IF
C
      EXPB=EXP(-WMAIN(XRATOM(IAT))*SSQ4)
C
C in anomalous mode we include the f'' contribution to the SF
      IF (QANOM) THEN
      XRFDPX=XRFDP(XRINDX(IAT))
      ELSE
      XRFDPX=ZERO
      END IF
C
      ARG=H2PI*XRNCSX(IAT)+K2PI*XRNCSY(IAT)+L2PI*XRNCSZ(IAT)
      CARG=DCMPLX(COS(ARG),SIN(ARG))
C
      IF (QDXYZ) THEN
C
      FLOCF=XRFQS(XRATOM(IAT))
     &      *QMAIN(XRATOM(IAT))*EXPB*DCMPLX(XRFFL,XRFDPX)
      CTEMP=FLOCF*CARG
      FA=DBLE(CTEMP)
      FB=DIMAG(CTEMP)
C
C compute derivative with respect to X, Y, Z
C x
      XRDX(XRATOM(IAT))=XRDX(XRATOM(IAT))-FB*DXA+FA*DXB
C
C y
      XRDY(XRATOM(IAT))=XRDY(XRATOM(IAT))-FB*DYA+FA*DYB
C
C z
      XRDZ(XRATOM(IAT))=XRDZ(XRATOM(IAT))-FB*DZA+FA*DZB
C
      ELSEIF (QDB) THEN
C
C compute derivatives with respect to B, Q, F', and F''
C
C occupancy
      FLOCFQ=XRFQS(XRATOM(IAT))*EXPB*DCMPLX(XRFFL,XRFDPX)
      CTEMP=FLOCFQ*CARG
      FA=DBLE(CTEMP)
      FB=DIMAG(CTEMP)
      XRDQ(XRATOM(IAT))=XRDQ(XRATOM(IAT))+FA*DQA+FB*DQB
C
C isotropic B-factor
      FLOCF=FLOCFQ*QMAIN(XRATOM(IAT))
      CTEMP=FLOCF*CARG
      FA=DBLE(CTEMP)
      FB=DIMAG(CTEMP)
      XRDT(XRATOM(IAT))=XRDT(XRATOM(IAT))+ FA*DTA+FB*DTB
C
C f'
      FLOCFP=XRFQS(XRATOM(IAT))
     &      *QMAIN(XRATOM(IAT))*EXPB*DCMPLX(ONE,ZERO)
      CTEMP=FLOCFP*CARG
      FA=DBLE(CTEMP)
      FB=DIMAG(CTEMP)
      XRDX(XRATOM(IAT))=XRDX(XRATOM(IAT))+FA*DQA+FB*DQB
C
C f''
      IF (QANOM) THEN
      FLOCFDP=XRFQS(XRATOM(IAT))
     &      *QMAIN(XRATOM(IAT))*EXPB*DCMPLX(ZERO,ONE)
      CTEMP=FLOCFDP*CARG
      FA=DBLE(CTEMP)
      FB=DIMAG(CTEMP)
      XRDY(XRATOM(IAT))=XRDY(XRATOM(IAT))+FA*DQA+FB*DQB
      END IF
C
      END IF
C
      END DO
      END IF
      END DO
C
      END IF
C
      END DO
      END DO
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4)')
     & ' XPRED: CPU-time:  direct structure factor calculation=',
     &   SECS2-SECS
      END IF
C
      RETURN
      END
C====================================================================
      SUBROUTINE XRSSQ(H,K,L,SSQ,XRTR)
C
C                            *      *      *
C Compute the square of S=H*A  + K*B  + L*C .
C We make use of XRINTR=(A,B,C) ==> XRTR=transpose(A*,B*,C*).
C
C Author: Axel T. Brunger
      IMPLICIT NONE
C I/O
      INTEGER H, K, L
      DOUBLE PRECISION SSQ, XRTR(3,3)
C begin
      SSQ=  (XRTR(1,1)*H + XRTR(2,1)*K + XRTR(3,1)*L)**2
     &    + (XRTR(1,2)*H + XRTR(2,2)*K + XRTR(3,2)*L)**2
     &    + (XRTR(1,3)*H + XRTR(2,3)*K + XRTR(3,3)*L)**2
      RETURN
      END
C====================================================================
      SUBROUTINE XPHASE(F,AMPLTD,PHASE)
C
C computes the phase in radians and amplitude from F=A+iB
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      DOUBLE COMPLEX F
      DOUBLE PRECISION AMPLTD, PHASE
C local
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, A, B
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
C compute the phase
      A=DBLE(F)
      B=DIMAG(F)
      AMPLTD=SQRT(A**2+B**2)
      IF (ABS(A).LT.RSMALL.AND.ABS(B).LT.RSMALL) THEN
      PHASE=ZERO
      ELSE IF (A.EQ.ZERO.AND.B.GT.ZERO) THEN
      PHASE=PI/TWO
      ELSE IF (A.EQ.ZERO.AND.B.LE.ZERO) THEN
      PHASE=THREE*PI/TWO
      ELSE
      PHASE=ATAN2(B,A)+TWO*PI
      END IF
      PHASE=MOD(PHASE+RSMALL,TWO*PI)
      RETURN
      END
C=====================================================================
      SUBROUTINE XAB(F,AMPLTD,PHASE)
C
C does the inverse of XPHASE, i.e. computes F=A+iB from phase in
C radians and amplitude.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      DOUBLE COMPLEX F
      DOUBLE PRECISION AMPLTD, PHASE
C local
C begin
      F=AMPLTD*DCMPLX(DCOS(PHASE),DSIN(PHASE))
      RETURN
      END
