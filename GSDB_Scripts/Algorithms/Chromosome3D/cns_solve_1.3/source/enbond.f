      SUBROUTINE ENBOND(VDW,ELEC,PVDW,PELEC,VDWV,ELECV,PVDWV,PELECV,N,
     &                  NNPIG,ANALYS)
C
C this is the front-end for the non-bonding routines for
C intra-molecular and symmetry-related interactions.
C
C Authors: Axel T. Brunger and Thomas Simonson
C ============================================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'param.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'fastnb.inc'
      DOUBLE PRECISION VDW, ELEC, PVDW, PELEC
      DOUBLE PRECISION VDWV, ELECV, PVDWV, PELECV
      INTEGER NNPIG,N
      CHARACTER*4 ANALYS
C local
      LOGICAL QUPDAT
      DOUBLE PRECISION DEVI, ESCALE, VSCALE
      DOUBLE PRECISION EVWGHT
      INTEGER XF, YF, ZF
C parameter
      DOUBLE PRECISION ZERO, ONE, NNNN
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, NNNN=9999.D0)
C begin
C
C make sure lookup tables for parameters are updated
      IF (UPNBLK) THEN
      UPNBLK=.FALSE.
      CALL NBUPDA
      IF (TIMER.GT.1) CALL DSPCPU(' ENBOND: after parameter lookup')
      END IF
C
C enforce update and bypass "previous position" test if any
C of the update flags are set.
      QUPDAT=UPNBEX(N)
     &       .OR. (UPNBLS(N).AND.(QENER(SSVDW).OR.QENER(SSELEC)))
     &       .OR. (UPNBSS(N).AND.(QENER(SSPVDW).OR.QENER(SSPELE)))
      IF (QUPDAT) UPNBEX(N)=.TRUE.
C
C do automatic update based on "previous position" test.
      IF (.NOT.QUPDAT.AND.NBXTOL.LT.NNNN) THEN
C
C check whether atoms have moved too far since last time list was
C obtained
      CALL XPDEVI(NATOM,X,Y,Z,HEAP(INBX),HEAP(INBY),HEAP(INBZ),DEVI)
      QUPDAT=DEVI.GT.NBXTOL
      END IF
C
      IF (QUPDAT) THEN
C
      IF (QFSTNB) THEN
C
C call fast grid XPSEAR - PDA 11/4/94
      CALL FXPSEAR(N,NNPIG)
      ELSE
C
C call standard XPSEAR routine
      CALL XPSEAR(N,NNPIG)
      END IF
C
C set the appropriate atom pair list update flag to false
      IF (QENER(SSVDW).OR.QENER(SSELEC)) UPNBLS(N)=.FALSE.
      IF (QENER(SSPVDW).OR.QENER(SSPELE)) UPNBSS(N)=.FALSE.
      UPNBEX(N)=.FALSE.
      END IF
C
C check to see whether we have to compute intra-molecular interactions
      IF (QENER(SSVDW).OR.QENER(SSELEC)) THEN
C
C check to see whether we've to compute electrostatics and /or vdw
      IF (QENER(SSVDW)) THEN
      VSCALE=ONE
      ELSE
      VSCALE=ZERO
      END IF
      IF (QENER(SSELEC)) THEN
      ESCALE=ONE
      ELSE
      ESCALE=ZERO
      END IF
      VSCALE=VSCALE*PIGWGHT(N,SSVDW)
      ESCALE=ESCALE*PIGWGHT(N,SSELEC)
C
C set QIMAG to .FALSE. and pass X, Y, Z to XJ, YJ, ZJ, set
C E-14-scale to 1.0
      IF (NNNB(N).GT.0) THEN
C
C initialize reading from the intramolecular interaction list
      CALL RESNB(LIST(1,N))
C
      CALL ENBRAN(VDW,ELEC,VDWV,ELECV,NATOM,HEAP(IJNB(N)),
     &        LIST(1,N),CG,MAXCN,CNBA,CNBB,CNBVR,LOOKUP,
     &        ESCALE,VSCALE,PIGAVWT(N,SSELEC),PIGAVWT(N,SSVDW),
     &        .FALSE.,.FALSE.,X,Y,Z,1,ANALYS)
      IF (TIMER.GT.1) CALL DSPCPU(' ENBOND: after intra interactions')
      END IF
C
C now the 1-4 interactions
      IF (NNB14(N).GT.0) THEN
C
C initialize reading from the symmetry-related interaction list
      CALL RESNB(LIST14(1,N))
C
C take care of some weighting
      ESCALE=ESCALE*E14FAC
      EVWGHT=PIGAVWT(N,SSELEC)*E14FAC
      CALL ENBRAN(VDW,ELEC,VDWV,ELECV,NATOM,HEAP(IINB14(N)),
     &    LIST14(1,N),CG,MAXCN,CNBA14,CNBB14,CBVR14,LOOKUP,
     &    ESCALE,VSCALE,EVWGHT,PIGAVWT(N,SSVDW),
     &    .FALSE.,.FALSE.,X,Y,Z,1,ANALYS)
      IF (TIMER.GT.1) CALL DSPCPU(' ENBOND: after 1-4 interactions')
      END IF
      END IF
C
C check to see whether we've to compute symmetry interactions
      IF (QENER(SSPVDW).OR.QENER(SSPELE)) THEN
C
C check to see whether we've to compute electrostatics and /or vdw
      IF (QENER(SSPVDW)) THEN
      VSCALE=ONE
      ELSE
      VSCALE=ZERO
      END IF
      IF (QENER(SSPELE)) THEN
      ESCALE=ONE
      ELSE
      ESCALE=ZERO
      END IF
      VSCALE=VSCALE*PIGWGHT(N,SSPVDW)
      ESCALE=ESCALE*PIGWGHT(N,SSPELE)
C
C allocate space for temporary XF, YF, ZF arrays.
      IF (NNSB(N).GT.0) THEN
C
C initialize reading from the symmetry interaction list
      CALL RESNB(SLIST(1,N))
      XF=ALLHP(IREAL8(NATOM))
      YF=ALLHP(IREAL8(NATOM))
      ZF=ALLHP(IREAL8(NATOM))
      CALL XPENE2(PVDW,PELEC,PVDWV,PELECV,HEAP(IJSB(N)),NATOM,
     &            SLIST(1,N),HEAP(XF),HEAP(YF),HEAP(ZF),
     &            ESCALE,VSCALE,PIGAVWT(N,SSPELE),PIGAVWT(N,SSPVDW),
     &            ANALYS)
      IF (TIMER.GT.1) CALL DSPCPU(' ENBOND: after symmetry interact.')
      CALL FREHP(ZF,IREAL8(NATOM))
      CALL FREHP(YF,IREAL8(NATOM))
      CALL FREHP(XF,IREAL8(NATOM))
C
      END IF
      END IF
C
      RETURN
      END
C===============================================================
      SUBROUTINE XPENE2(EVDW,ELEC,EVDWV,ELECV,JNB,LEN,SLISTL,XF,YF,ZF,
     &                  ESCALE,VSCALE,EVWGHT,VVWGHT,ANALYS)
C
C subroutine computes intermolecular repulsion energy
C Authors: Axel T. Brunger and William I Weis
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'param.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'ncs.inc'
      DOUBLE PRECISION EVDW, ELEC
      DOUBLE PRECISION EVDWV, ELECV
      DOUBLE PRECISION EVWGHT,VVWGHT
      INTEGER LEN, JNB(LEN,*), SLISTL(4)
      DOUBLE PRECISION XF(*), YF(*), ZF(*), ESCALE, VSCALE
      CHARACTER*4 ANALYS
C local
      INTEGER I, J, K, L, JSYM, JJSYM
C parameters
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C QIMAG (in ENBRAN) tells us to compute symmetry related molecules;
C LNCSST will tell us
C whether they are xtal images or NCS related molecules (mutually exclusive
C options).
C Since we are passing transformed coordinates to ENBRAN, for NCS we
C can simply take distances, as we do for intra interations,
C i.e. QIMAG=.FALSE.
C
C List was set with either xtal symmetry images or NCS images; choose the
C appropriate one.
      IF (.NOT.LNCSST) THEN
      JJSYM=XRNSYM
      ELSE
      JJSYM=NNSYM
      END IF
C
C loop over symmetry operators (either xtal or NCS)
      DO JSYM=1,JJSYM
C
      IF (.NOT.LNCSST) THEN
C Define transformation operator for xtal symmetry forces (F * S * F-1)
C this operator is used in XPFORC which is called from ENBRAN
C
      DO I=1,3
      DO J=1,3
      XPOPER(I,J)=ZERO
      DO K=1,3
      DO L=1,3
      XPOPER(I,J)=XPOPER(I,J) +XRINTR(I,K)*XRSYMM(JSYM,K,L)*XRTR(L,J)
      END DO
      END DO
      END DO
      END DO
C
C compute symmetry related molecule (in fractional coordinates)
      CALL XPSYMO(NATOM,X,Y,Z,JSYM,XF,YF,ZF)
C
      CALL ENBRAN(EVDW,ELEC,EVDWV,ELECV,NATOM,
     &            JNB(1,JSYM),SLISTL,
     &      CG,MAXCN,CNBA,CNBB,CNBVR,LOOKUP,ESCALE,
     &      VSCALE,EVWGHT,VVWGHT,
     &      .TRUE.,.TRUE.,XF,YF,ZF,JSYM,ANALYS)
C
      ELSE
C The transformation operator for NCS case is simply the NCS operator
      DO I=1,3
      DO J=1,3
      XPOPER(I,J)=NCSOP(JSYM,I,J)
      END DO
      END DO
C
C compute non-crystallographic symmetry related molecule (orthogonal coords)
      CALL NCSXYZ(NATOM,X,Y,Z,JSYM,XF,YF,ZF)
      CALL ENBRAN(EVDW,ELEC,EVDWV,ELECV,NATOM,
     &      JNB(1,JSYM),SLISTL,
     &      CG,MAXCN,CNBA,CNBB,CNBVR,LOOKUP,ESCALE,
     &      VSCALE,EVWGHT,VVWGHT,
     &      .FALSE.,.TRUE.,XF,YF,ZF,JSYM,ANALYS)
      END IF
      END DO
      RETURN
      END
C=================================================================
      SUBROUTINE NBUPDA
C
C fills LOOKUP array for type-based parameters.  Note: the
C LOOKUP array is already set for atom-based parameters.
C
C Author: Axel T. Brunger
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'param.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
C local
      INTEGER I, II, IFLAG
C begin
C
C get the Lennard-Jones type-based parameters
      IFLAG=0
      DO I=1,NATOM
C
C ignore the atom-based slots
      IF (QLOOKU(I)) THEN
C
C try to match the types
      II=1
      DO WHILE (CNAC(II).NE.IAC(I).AND.II.LT.NCN)
      II=II+1
      END DO
C
      IF (CNAC(II).EQ.IAC(I)) THEN
      LOOKUP(I)=II
      ELSE
C
      WRITE(6,'(A)')
     &' %NBUPDA-ERR: missing nonbonded Lennard-Jones parameters %%%%%%%'
      WRITE(6,'(9A)')
     & '  ATOM: SEGId="',SEGID(I),'",  RESId="',RESID(I),
     & '",  NAME="',TYPE(I),'",  CHEMical="',IAC(I),'"'
      IFLAG=IFLAG+1
      END IF
C
      END IF
      END DO
      IF (IFLAG.GT.0) THEN
      CALL WRNDIE(-5,'NBUPDA','program will be aborted.')
      ELSE
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)')
     & ' NBUPDA: type-based van der Waals parameters retrieved. '
      END IF
      END IF
C
      RETURN
      END
C=================================================================
      SUBROUTINE XPDEVI(NATOM,X,Y,Z,XPX,XPY,XPZ,DEVI)
C
C routine determines maximum deviation between X,Y,Z and XPX, XPY, XPZ
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER NATOM
      DOUBLE PRECISION X(*), Y(*), Z(*), XPX(*), XPY(*), XPZ(*), DEVI
C local
      INTEGER I
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      DEVI=ZERO
      DO I=1,NATOM
      DEVI=MAX(DEVI,
     &   (X(I)-XPX(I))**2 + (Y(I)-XPY(I))**2 + (Z(I)-XPZ(I))**2)
      END DO
      DEVI=SQRT(DEVI)
      RETURN
      END
C==================================================================
      SUBROUTINE ENBRAN(VDW,ELEC,VDWV,ELECV,NATOM,JNB,LISTL,
     &  CG,MAXCN,CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,EVWGHT,VVWGHT,
     &  QIMAG,QFORC,XJ,YJ,ZJ,JSYM,ANALYS)
C
C this routine branches into the various nonbond-option
C routines
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      DOUBLE PRECISION VDW, ELEC, VDWV, ELECV
      DOUBLE PRECISION EVWGHT, VVWGHT
      INTEGER NATOM, JNB(*), LISTL(4)
      DOUBLE PRECISION CG(*)
      INTEGER MAXCN
      DOUBLE PRECISION CNBA(MAXCN,MAXCN),CNBB(MAXCN,MAXCN)
      DOUBLE PRECISION CNBVR(MAXCN,MAXCN)
      INTEGER LOOKUP(*)
      DOUBLE PRECISION ESCALE, VSCALE
      LOGICAL QIMAG, QFORC
      DOUBLE PRECISION XJ(*), YJ(*), ZJ(*)
      INTEGER JSYM
      CHARACTER*4 ANALYS
C local
      INTEGER DIM, I
      LOGICAL LBSF
C pointer
      INTEGER CGL, JND, CNBAL, CNBBL, XJL, YJL, ZJL, CNBIL
C parameter
      DOUBLE PRECISION ZERO,ONE,FOUR
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,FOUR=4.0D0)
C begin
C
C determine maximum number of interactions between one atom
C and any other atom
      DIM=JNB(1)
      DO I=2,NATOM
      DIM=MAX(DIM,JNB(I))
      END DO
      DIM=DIM+2
C
C allocate space for the temporary vectorization arrays
      JND=ALLHP(INTEG4(DIM))
      CGL=ALLHP(IREAL8(DIM))
      CNBAL=ALLHP(IREAL8(DIM))
      CNBBL=ALLHP(IREAL8(DIM))
      CNBIL=ALLHP(IREAL8(DIM))
      XJL=ALLHP(IREAL8(DIM))
      YJL=ALLHP(IREAL8(DIM))
      ZJL=ALLHP(IREAL8(DIM))
C
C branch to appropriate routine corresponding to current
C non-bonded options
      LBSF=(LCONS.AND.(LSHFT.OR.LSHFO).AND..NOT.LVSHFT)
      IF (NATOM.GT.0) THEN
      IF (ANALYS.EQ.'ANAL') THEN
C
C call nonbonded distance analysis routine
      CALL ENBANL(VDW,ELEC,VDWV,ELECV,NATOM,JNB,LISTL,CG,MAXCN,
     &         CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,EVWGHT,VVWGHT,
     &         QIMAG,QFORC,XJ,YJ,ZJ,
     &         HEAP(JND),HEAP(CGL),HEAP(CNBAL),HEAP(CNBBL),HEAP(CNBIL),
     &         HEAP(XJL),HEAP(YJL),HEAP(ZJL),JSYM)
C
C added by JJK 8/2/95
C
C Rattract is the distance at which the attractive part of the
C attractive-repulsive term comes on
C
      ELSE IF (RATT.GT.RSMALL) THEN
C
C call nonbonding routine with attractive-repulsive energy
C
      CALL ENBAR(VDW,ELEC,VDWV,ELECV,NATOM,JNB,LISTL,CG,MAXCN,
     &       cnbb,cnba,LOOKUP,ESCALE,VSCALE,EVWGHT,VVWGHT,QIMAG,
     &       QFORC,XJ,YJ,ZJ,
     &       HEAP(JND),HEAP(CGL),HEAP(CNBAL),HEAP(CNBBL),HEAP(CNBIL),
     &       HEAP(XJL),HEAP(YJL),HEAP(ZJL))
      ELSE IF (FREPEL.GT.RSMALL) THEN
C
C call vectorized nonbonding routine with repel potential and
C no electrostatics
      CALL ENBREP(VDW,ELEC,VDWV,ELECV,NATOM,JNB,LISTL,CG,MAXCN,
     &       CNBVR,LOOKUP,ESCALE,VSCALE,EVWGHT,VVWGHT,QIMAG,
     &       QFORC,XJ,YJ,ZJ,
     &       HEAP(JND),HEAP(CGL),HEAP(CNBAL),HEAP(CNBBL),HEAP(CNBIL),
     &       HEAP(XJL),HEAP(YJL),HEAP(ZJL))
C
      ELSE IF (LCONS.AND.QTRUNC) THEN
C
C call vectorized nonbonding routine just with truncation (c-diel.)
      CALL ENBTR(VDW,ELEC,VDWV,ELECV,NATOM,JNB,LISTL,CG,MAXCN,
     &         CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,EVWGHT,VVWGHT,
     &         QIMAG,QFORC,XJ,YJ,ZJ,
     &         HEAP(JND),HEAP(CGL),HEAP(CNBAL),HEAP(CNBBL),HEAP(CNBIL),
     &         HEAP(XJL),HEAP(YJL),HEAP(ZJL))
      ELSE IF (LBSF) THEN
      IF (LSHFO.AND.(BSHFO.EQ.FOUR)) THEN
C
C call vectorized nonbonding routine with VdW switching and
C electrostatics force shifting (c-diel.)
      CALL ENBFS4(VDW,ELEC,VDWV,ELECV,NATOM,JNB,LISTL,CG,MAXCN,
     &       CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,EVWGHT,VVWGHT,
     &       QIMAG,QFORC,XJ,YJ,ZJ,
     &       HEAP(JND),HEAP(CGL),HEAP(CNBAL),HEAP(CNBBL),HEAP(CNBIL),
     &       HEAP(XJL),HEAP(YJL),HEAP(ZJL))
      ELSE IF (LSHFO.AND.(BSHFO.EQ.ONE)) THEN
C
C call vectorized nonbonding routine with VdW switching and
C electrostatics force shifting (c-diel.)
      CALL ENBFS1(VDW,ELEC,VDWV,ELECV,NATOM,JNB,LISTL,CG,MAXCN,
     &       CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,EVWGHT,VVWGHT,
     &       QIMAG,QFORC,XJ,YJ,ZJ,
     &       HEAP(JND),HEAP(CGL),HEAP(CNBAL),HEAP(CNBBL),HEAP(CNBIL),
     &       HEAP(XJL),HEAP(YJL),HEAP(ZJL))
      ELSE IF (LSHFO) THEN
C
C call vectorized nonbonding routine with VdW switching and
C electrostatics generalized force shifting (c-diel.)
      CALL ENBFSB(VDW,ELEC,VDWV,ELECV,NATOM,JNB,LISTL,CG,MAXCN,
     &       CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,EVWGHT,VVWGHT,
     &       QIMAG,QFORC,XJ,YJ,ZJ,
     &       HEAP(JND),HEAP(CGL),HEAP(CNBAL),HEAP(CNBBL),HEAP(CNBIL),
     &       HEAP(XJL),HEAP(YJL),HEAP(ZJL))
      ELSE
C
C call vectorized nonbonding routine with VdW switching and
C electrostatics shifting (c-diel.)
      CALL ENBSF(VDW,ELEC,VDWV,ELECV,NATOM,JNB,LISTL,CG,MAXCN,
     &       CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,EVWGHT,VVWGHT,
     &       QIMAG,QFORC,XJ,YJ,ZJ,
     &       HEAP(JND),HEAP(CGL),HEAP(CNBAL),HEAP(CNBBL),HEAP(CNBIL),
     &       HEAP(XJL),HEAP(YJL),HEAP(ZJL))
      END IF
      ELSE IF (.NOT.LCONS.AND..NOT.LSHFT.AND..NOT.LSHFO.AND.
     &         .NOT.LVSHFT) THEN
C
C call vectorized nonbonding routine with 1/r -dependent dielectric
C and switching on both electrostatics and VdW term.
      CALL ENBRD(VDW,ELEC,VDWV,ELECV,NATOM,JNB,LISTL,CG,MAXCN,
     &       CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,EVWGHT,VVWGHT,
     &       QIMAG,QFORC,XJ,YJ,ZJ,
     &       HEAP(JND),HEAP(CGL),HEAP(CNBAL),HEAP(CNBBL),HEAP(CNBIL),
     &       HEAP(XJL),HEAP(YJL),HEAP(ZJL))
      ELSE
      WRITE(6,'(A)') ' %ENBOND-ERR: non-bonding option not supported'
      CALL DIE
      END IF
      END IF
C
C de-allocate space for the temporary vectorization arrays
      CALL FREHP(ZJL,IREAL8(DIM))
      CALL FREHP(YJL,IREAL8(DIM))
      CALL FREHP(XJL,IREAL8(DIM))
      CALL FREHP(CNBIL,IREAL8(DIM))
      CALL FREHP(CNBBL,IREAL8(DIM))
      CALL FREHP(CNBAL,IREAL8(DIM))
      CALL FREHP(CGL,IREAL8(DIM))
      CALL FREHP(JND,INTEG4(DIM))
C
      RETURN
      END
C====================================================================
      SUBROUTINE ENBTR(EVDW,ELEC,EVDWV,ELECV,NATOM,JNB,LISTL,CG,
     &          MAXCN,CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,
     &          EVWGHT,VVWGHT,QIMAG,QFORC,XJ,YJ,ZJ,JND,CGL,
     &          CNBAL,CNBBL,CNBIL,XJL,YJL,ZJL)
C
C Features:
C           Lennard-Jones and free electrostatic potential
C           pure cutoff without shifting or switching
C
C By Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'heap.inc'
      DOUBLE PRECISION EVDW, ELEC
      DOUBLE PRECISION EVDWV, ELECV
      DOUBLE PRECISION EVWGHT, VVWGHT
      INTEGER NATOM, JNB(*), LISTL(4)
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
C local
      INTEGER I, J, N
      DOUBLE PRECISION EVDWL, EVDWK, ELECL, ELECK, DEVDWK, DELECK
      DOUBLE PRECISION S, R2, SIG6, SIG12, CGIL, CGF
      DOUBLE PRECISION DELTX, DELTY, DELTZ
      DOUBLE PRECISION CTOFN2
C parameters
      DOUBLE PRECISION ZERO, ONE, SIX, TWO
      PARAMETER (ZERO=0.D0, ONE=1.D0, SIX=6.D0, TWO=2.0D0)
C begin
C
      CGF=CCELEC/EPS
      CTOFN2=CTOFNB**2
C
C main loop over atom I
      DO I=1,NATOM
C
      N=JNB(I)
      IF (N.GT.0) THEN
C
C define local index list JND
      CALL GETNB(N,JND,LISTL)
C
C gather Lennard-Jones A, B and charge CG
C gather "j" coordinates
      DO J=1,N
      CNBBL(J)=CNBB(LOOKUP(JND(J)),LOOKUP(I))
      CNBAL(J)=CNBA(LOOKUP(JND(J)),LOOKUP(I))
      CNBIL(J)=(INHIBT*CNBVR(LOOKUP(JND(J)),(LOOKUP(I))))**2
      END DO
      DO J=1,N
      CGL(J)=CG(JND(J))
      XJL(J)=XJ(JND(J))
      YJL(J)=YJ(JND(J))
      ZJL(J)=ZJ(JND(J))
      END DO
C
C appropriately scale self-interactions
      IF (JND(1).EQ.I) THEN
      CNBBL(1)=CNBBL(1)/TWO
      CNBAL(1)=CNBAL(1)/TWO
      CGL(1)=CGL(1)/TWO
      END IF
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
      ELECL=ZERO
      EVDWL=ZERO
      CGIL=CG(I)*CGF
C
C compute Lennard Jones and electrostatic expression
      DO J=1,N
C
      DELTX=XJL(J)
      DELTY=YJL(J)
      DELTZ=ZJL(J)
      S=MAX(CNBIL(J),DELTX**2+DELTY**2+DELTZ**2)
C
      R2=ONE/S
C
C evaluate Lennard-Jones
C
C               A      B
C     E   =  ( ---  - --- )
C      VdW       12     6
C               R      R
C
      SIG6=R2*R2*R2
      SIG12=CNBAL(J)*SIG6**2
      SIG6=CNBBL(J)*SIG6
      EVDWK=(SIG12-SIG6)
      DEVDWK=VSCALE*SIX*R2*(SIG6-SIG12-SIG12)
C
C evaluate electrostatic potential
C            C
C     Q Q  -----
C      i j EPS R
C
      ELECK=CGIL*CGL(J)*SQRT(R2)
      DELECK=-R2*ELECK*ESCALE
C
C forces (including PIG's weights)
      DELTX=DELTX*(DEVDWK+DELECK)
      DELTY=DELTY*(DEVDWK+DELECK)
      DELTZ=DELTZ*(DEVDWK+DELECK)
C
      DX(I)=DX(I)+DELTX
      DY(I)=DY(I)+DELTY
      DZ(I)=DZ(I)+DELTZ
C
      XJL(J)=DELTX
      YJL(J)=DELTY
      ZJL(J)=DELTZ
C
C local accumulation
      EVDWL=EVDWL+EVDWK
      ELECL=ELECL+ELECK
      END DO
C
C global accumulation
      EVDW=EVDW+EVDWL*VSCALE
      ELEC=ELEC+ELECL*ESCALE
      EVDWV=EVDWV+EVDWL*VVWGHT
      ELECV=ELECV+ELECL*EVWGHT
C
C transform "j" image forces
      IF (QFORC) CALL XPFORC(N,XJL,YJL,ZJL)
C accumulate "j" forces
C >>>>>>ignore vector dependencies<<<<<<<
C$DIR NO_RECURRENCE
CDIR$ IVDEP
C*$*ASSERT PERMUTATION (JND)
      DO J=1,N
      DX(JND(J))=DX(JND(J))-XJL(J)
      DY(JND(J))=DY(JND(J))-YJL(J)
      DZ(JND(J))=DZ(JND(J))-ZJL(J)
      END DO
C
      END IF
C
C end of main loop
      END DO
C
      RETURN
      END
C====================================================================
      SUBROUTINE ENBSF(EVDW,ELEC,EVDWV,ELECV,NATOM,JNB,LISTL,CG,
     &         MAXCN,CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,
     &         EVWGHT,VVWGHT,QIMAG,QFORC,XJ,YJ,ZJ,JND,CGL,
     &         CNBAL,CNBBL,CNBIL,XJL,YJL,ZJL)
C
C Features:
C    Lennard-Jones and free electrostatic potential
C    atom by atom switching Lennard-Jones and shifting electrostatics
C
C By Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'nbonds.inc'
      DOUBLE PRECISION EVDW, ELEC
      DOUBLE PRECISION EVDWV, ELECV
      DOUBLE PRECISION EVWGHT, VVWGHT
      INTEGER NATOM, JNB(*), LISTL(4)
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
C local
      INTEGER I, J, N
      DOUBLE PRECISION EN, EVDWL, EVDWK, ELECL, ELECK, DEVDWK, DELECK
      DOUBLE PRECISION R2, SIG6, SIG12, CGIL, CGF
      DOUBLE PRECISION G1, G3, C2OFNB, C2ROF2, CHROF2
      DOUBLE PRECISION C2ONNB, RUL3, RUL12, S, SS,  RIJL, RIJU, FUNCT
      DOUBLE PRECISION DELTX, DELTY, DELTZ, DFN
      DOUBLE PRECISION WSIX
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, SIX, TWELVE
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0)
      PARAMETER (SIX=6.D0, TWELVE=12.D0)
C begin
C
      WSIX=SIX*VSCALE
      CGF=CCELEC/EPS
      C2OFNB=CTOFNB*CTOFNB
      C2ONNB=CTONNB*CTONNB
      RUL3=ONE/(C2OFNB-C2ONNB)**3
      RUL12=RUL3*TWELVE*VSCALE
      C2ROF2=-ONE/C2OFNB
      CHROF2=-FOUR/C2OFNB
C main loop over atom I
      DO I=1,NATOM
C
      N=JNB(I)
      IF (N.GT.0) THEN
C
C define local index list JND
      CALL GETNB(N,JND,LISTL)
C
C gather Lennard-Jones A, B and charge CG
C gather "j" coordinates
      DO J=1,N
      CNBBL(J)=CNBB(LOOKUP(JND(J)),LOOKUP(I))
      CNBAL(J)=CNBA(LOOKUP(JND(J)),LOOKUP(I))
      CNBIL(J)=(INHIBT*CNBVR(LOOKUP(JND(J)),(LOOKUP(I))))**2
      END DO
      DO J=1,N
      CGL(J)=CG(JND(J))
      XJL(J)=XJ(JND(J))
      YJL(J)=YJ(JND(J))
      ZJL(J)=ZJ(JND(J))
      END DO
C
C appropriately scale self-interactions
      IF (JND(1).EQ.I) THEN
      CNBBL(1)=CNBBL(1)/TWO
      CNBAL(1)=CNBAL(1)/TWO
      CGL(1)=CGL(1)/TWO
      END IF
CC
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
      ELECL=ZERO
      EVDWL=ZERO
      CGIL=CG(I)*CGF
C
C compute Lennard Jones and electrostatic expression
      DO J=1,N
C
      DELTX=XJL(J)
      DELTY=YJL(J)
      DELTZ=ZJL(J)
      S=MAX(CNBIL(J),DELTX**2+DELTY**2+DELTZ**2)
      R2=ONE/S
C
C evaluate Lennard Jones terms multiplied by switching functions
C
C                 A       B
C      E    = (  ---  -  ---  )  *  Sw(R)
C       VdW        12      6
C                 R       R
C
C where the switching function Sw(R) is defined by
C
C
C            (                0                          R > R
C            (                                                off
C            (
C            (     2    2  2     2   2       2    2
C            (  ( R  - R  ) * ( R - R  - 3( R  - R  ) )
C            (     off          off          on
C     Sw(R)= (  -------------------------------------
C            (                2     2   3
C            (              (R   - R   )
C            (                off   on
C            (
C            (
C            (                1                          R < R
C                                                             on
C
      SS=MIN(MAX(S,C2ONNB),C2OFNB)
      RIJL=C2ONNB-SS
      RIJU=C2OFNB-SS
      FUNCT=RIJU**2*(RIJU-THREE*RIJL)*RUL3
      DFN=RIJL*RIJU*RUL12
C
      SIG6=R2*R2*R2
      SIG12=CNBAL(J)*SIG6*SIG6
      SIG6=CNBBL(J)*SIG6
      EN=SIG12-SIG6
      EVDWK=EN*FUNCT
      DEVDWK=EN*DFN-WSIX*R2*(EN+SIG12)*FUNCT
C
C evaluate shifted electrostatic potential
C                            2
C            C              R       2
C     Q Q  -----   (  1 -  ---    )
C      i j EPS R             2
C                           R
C                            off
C
C for R less than CTOFNB and zero otherwise
C
      G1=MAX(ZERO, ONE+S*C2ROF2)
      G3=CGIL*CGL(J)*G1*SQRT(R2)
      ELECK=G3*G1
      DELECK=(-ELECK*R2+CHROF2*G3)*ESCALE
C
C forces (including PIG's weights)
      DELTX=DELTX*(DEVDWK+DELECK)
      DELTY=DELTY*(DEVDWK+DELECK)
      DELTZ=DELTZ*(DEVDWK+DELECK)
C
      DX(I)=DX(I)+DELTX
      DY(I)=DY(I)+DELTY
      DZ(I)=DZ(I)+DELTZ
C
      XJL(J)=DELTX
      YJL(J)=DELTY
      ZJL(J)=DELTZ
C
C local accumulation
      EVDWL=EVDWL+EVDWK
      ELECL=ELECL+ELECK
      END DO
C
C global accumulation
      EVDW=EVDW+EVDWL*VSCALE
      ELEC=ELEC+ELECL*ESCALE
      EVDWV=EVDWV+EVDWL*VVWGHT
      ELECV=ELECV+ELECL*EVWGHT
C
C transform "j" image forces
      IF (QFORC) CALL XPFORC(N,XJL,YJL,ZJL)
C accumulate "j" forces
C >>>>>>ignore vector dependencies<<<<<<<
C$DIR NO_RECURRENCE
CDIR$ IVDEP
C*$*ASSERT PERMUTATION (JND)
      DO J=1,N
      DX(JND(J))=DX(JND(J))-XJL(J)
      DY(JND(J))=DY(JND(J))-YJL(J)
      DZ(JND(J))=DZ(JND(J))-ZJL(J)
      END DO
C
      END IF
C
C end of main loop
      END DO
C
      RETURN
      END
C====================================================================
      SUBROUTINE ENBFS1(EVDW,ELEC,EVDWV,ELECV,NATOM,JNB,LISTL,CG,
     &         MAXCN,CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,
     &         EVWGHT,VVWGHT,QIMAG,QFORC,XJ,YJ,ZJ,JND,CGL,
     &         CNBAL,CNBBL,CNBIL,XJL,YJL,ZJL)
C
C Features:
C    Lennard-Jones and free electrostatic potential
C    atom by atom switching Lennard-Jones and force shifting electrostatics
C    with B=1 in Eq 15 J Comp Chem 15 (1994) 667
C
C
C Author: Axel T. Brunger
C modified by Alexandre Bonvin, 11-SEP-95
C =======================================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'nbonds.inc'
      DOUBLE PRECISION EVDW, ELEC
      DOUBLE PRECISION EVDWV, ELECV
      DOUBLE PRECISION EVWGHT, VVWGHT
      INTEGER NATOM, JNB(*), LISTL(4)
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
C local
      INTEGER I, J, N
      DOUBLE PRECISION EN, EVDWL, EVDWK, ELECL, ELECK, DEVDWK, DELECK
      DOUBLE PRECISION R2, SIG6, SIG12, CGIL, CGF
      DOUBLE PRECISION G1, G3, C2OFNB, C2ROF2, CHROF2, CROFNB, CHROFN
      DOUBLE PRECISION C2ONNB, RUL3, RUL12, S, SS,  RIJL, RIJU, FUNCT
      DOUBLE PRECISION DELTX, DELTY, DELTZ, DFN, SQS, RSQS
      DOUBLE PRECISION WSIX
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, SIX, TWELVE
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0)
      PARAMETER (SIX=6.D0, TWELVE=12.D0)
C begin
C
      WSIX=SIX*VSCALE
      CGF=CCELEC/EPS
      C2OFNB=CTOFNB*CTOFNB
      CROFNB=-ONE/CTOFNB
      CHROFN=-TWO/CTOFNB
      C2ONNB=CTONNB*CTONNB
      RUL3=ONE/(C2OFNB-C2ONNB)**3
      RUL12=RUL3*TWELVE*VSCALE
      C2ROF2=-ONE/C2OFNB
      CHROF2=-FOUR/C2OFNB
C main loop over atom I
      DO I=1,NATOM
C
      N=JNB(I)
      IF (N.GT.0) THEN
C
C define local index list JND
      CALL GETNB(N,JND,LISTL)
C
C gather Lennard-Jones A, B and charge CG
C gather "j" coordinates
      DO J=1,N
      CNBBL(J)=CNBB(LOOKUP(JND(J)),LOOKUP(I))
      CNBAL(J)=CNBA(LOOKUP(JND(J)),LOOKUP(I))
      CNBIL(J)=(INHIBT*CNBVR(LOOKUP(JND(J)),LOOKUP(I)))**2
      END DO
      DO J=1,N
      CGL(J)=CG(JND(J))
      XJL(J)=XJ(JND(J))
      YJL(J)=YJ(JND(J))
      ZJL(J)=ZJ(JND(J))
      END DO
C
C appropriately scale self-interactions
      IF (JND(1).EQ.I) THEN
      CNBBL(1)=CNBBL(1)/TWO
      CNBAL(1)=CNBAL(1)/TWO
      CGL(1)=CGL(1)/TWO
      END IF
CC
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
      ELECL=ZERO
      EVDWL=ZERO
      CGIL=CG(I)*CGF
C
C compute Lennard Jones and electrostatic expression
      DO J=1,N
C
      DELTX=XJL(J)
      DELTY=YJL(J)
      DELTZ=ZJL(J)
      S=MAX(CNBIL(J),DELTX**2+DELTY**2+DELTZ**2)
      SQS=SQRT(S)
      RSQS=ONE/SQS
      R2=ONE/S
C
C evaluate Lennard Jones terms multiplied by switching functions
C
C                 A       B
C      E    = (  ---  -  ---  )  *  Sw(R)
C       VdW        12      6
C                 R       R
C
C where the switching function Sw(R) is defined by
C
C
C            (                0                          R > R
C            (                                                off
C            (
C            (     2    2  2     2   2       2    2
C            (  ( R  - R  ) * ( R - R  - 3( R  - R  ) )
C            (     off          off          on
C     Sw(R)= (  -------------------------------------
C            (                2     2   3
C            (              (R   - R   )
C            (                off   on
C            (
C            (
C            (                1                          R < R
C                                                             on
C
      SS=MIN(MAX(S,C2ONNB),C2OFNB)
      RIJL=C2ONNB-SS
      RIJU=C2OFNB-SS
      FUNCT=RIJU**2*(RIJU-THREE*RIJL)*RUL3
      DFN=RIJL*RIJU*RUL12
C
      SIG6=R2*R2*R2
      SIG12=CNBAL(J)*SIG6*SIG6
      SIG6=CNBBL(J)*SIG6
      EN=SIG12-SIG6
      EVDWK=EN*FUNCT
      DEVDWK=EN*DFN-WSIX*R2*(EN+SIG12)*FUNCT
C
C evaluate force shifted electrostatic potential
C
C            C              R      2
C     Q Q  -----   (  1 -  ---   )
C      i j EPS R            R
C                            off
C
C
C for R less than CTOFNB and zero otherwise
C
      G1=MAX(ZERO, ONE+SQS*CROFNB)
      G3=CGIL*CGL(J)*G1*RSQS
      ELECK=G3*G1
      DELECK=(-ELECK*R2+CHROFN*G3*RSQS)*ESCALE
C
C forces (including PIG's weights)
      DELTX=DELTX*(DEVDWK+DELECK)
      DELTY=DELTY*(DEVDWK+DELECK)
      DELTZ=DELTZ*(DEVDWK+DELECK)
C
      DX(I)=DX(I)+DELTX
      DY(I)=DY(I)+DELTY
      DZ(I)=DZ(I)+DELTZ
C
      XJL(J)=DELTX
      YJL(J)=DELTY
      ZJL(J)=DELTZ
C
C local accumulation
      EVDWL=EVDWL+EVDWK
      ELECL=ELECL+ELECK
      END DO
C
C global accumulation
      EVDW=EVDW+EVDWL*VSCALE
      ELEC=ELEC+ELECL*ESCALE
      EVDWV=EVDWV+EVDWL*VVWGHT
      ELECV=ELECV+ELECL*EVWGHT
C
C transform "j" image forces
      IF (QFORC) CALL XPFORC(N,XJL,YJL,ZJL)
C accumulate "j" forces
C >>>>>>ignore vector dependencies<<<<<<<
C$DIR NO_RECURRENCE
CDIR$ IVDEP
C*$*ASSERT PERMUTATION (JND)
      DO J=1,N
      DX(JND(J))=DX(JND(J))-XJL(J)
      DY(JND(J))=DY(JND(J))-YJL(J)
      DZ(JND(J))=DZ(JND(J))-ZJL(J)
      END DO
C
      END IF
C
C end of main loop
      END DO
C
      RETURN
      END
C====================================================================
      SUBROUTINE ENBFS4(EVDW,ELEC,EVDWV,ELECV,NATOM,JNB,LISTL,CG,
     &         MAXCN,CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,
     &         EVWGHT,VVWGHT,QIMAG,QFORC,XJ,YJ,ZJ,JND,CGL,
     &         CNBAL,CNBBL,CNBIL,XJL,YJL,ZJL)
C
C Features:
C    Lennard-Jones and free electrostatic potential
C    atom by atom switching Lennard-Jones and force shifting electrostatics
C    with B=4 in Eq 15 J Comp Chem 15 (1994) 667
C
C
C Author: Axel T. Brunger
C modified by Alexandre Bonvin Yale, 01-SEP-95
C =======================================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'nbonds.inc'
      DOUBLE PRECISION EVDW, ELEC
      DOUBLE PRECISION EVDWV, ELECV
      DOUBLE PRECISION EVWGHT, VVWGHT
      INTEGER NATOM, JNB(*), LISTL(4)
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
C local
      INTEGER I, J, N
      DOUBLE PRECISION EN, EVDWL, EVDWK, ELECL, ELECK, DEVDWK, DELECK
      DOUBLE PRECISION R2, SIG6, SIG12, CGIL, CGF
      DOUBLE PRECISION G1, G3, C2OFNB, C2ROF2, CHROF2, CROFNB, CHROFN
      DOUBLE PRECISION C2ONNB, RUL3, RUL12, S, SS,  RIJL, RIJU, FUNCT
      DOUBLE PRECISION DELTX, DELTY, DELTZ, DFN, SQS, RSQS
      DOUBLE PRECISION WSIX, CCSHFO, DDSHFO, BSHFOP, G2
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, SIX, TWELVE
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0)
      PARAMETER (SIX=6.D0, TWELVE=12.D0)
C begin
C
      WSIX=SIX*VSCALE
      CGF=CCELEC/EPS
      C2OFNB=CTOFNB*CTOFNB
      CROFNB=-ONE/CTOFNB
      BSHFOP=BSHFO+ONE
      CCSHFO=ONE/(BSHFO*CTOFNB**BSHFOP)
      DDSHFO=BSHFOP*CROFNB/BSHFO
      CHROFN=-TWO/CTOFNB
      C2ONNB=CTONNB*CTONNB
      RUL3=ONE/(C2OFNB-C2ONNB)**3
      RUL12=RUL3*TWELVE*VSCALE
      C2ROF2=-ONE/C2OFNB
      CHROF2=-FOUR/C2OFNB
C main loop over atom I
      DO I=1,NATOM
C
      N=JNB(I)
      IF (N.GT.0) THEN
C
C define local index list JND
      CALL GETNB(N,JND,LISTL)
C
C gather Lennard-Jones A, B and charge CG
C gather "j" coordinates
      DO J=1,N
      CNBBL(J)=CNBB(LOOKUP(JND(J)),LOOKUP(I))
      CNBAL(J)=CNBA(LOOKUP(JND(J)),LOOKUP(I))
      CNBIL(J)=(INHIBT*CNBVR(LOOKUP(JND(J)),LOOKUP(I)))**2
      END DO
      DO J=1,N
      CGL(J)=CG(JND(J))
      XJL(J)=XJ(JND(J))
      YJL(J)=YJ(JND(J))
      ZJL(J)=ZJ(JND(J))
      END DO
C
C appropriately scale self-interactions
      IF (JND(1).EQ.I) THEN
      CNBBL(1)=CNBBL(1)/TWO
      CNBAL(1)=CNBAL(1)/TWO
      CGL(1)=CGL(1)/TWO
      END IF
CC
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
      ELECL=ZERO
      EVDWL=ZERO
      CGIL=CG(I)*CGF
C
C compute Lennard Jones and electrostatic expression
      DO J=1,N
C
      DELTX=XJL(J)
      DELTY=YJL(J)
      DELTZ=ZJL(J)
      S=MAX(CNBIL(J),DELTX**2+DELTY**2+DELTZ**2)
      SQS=SQRT(S)
      RSQS=ONE/SQS
      R2=ONE/S
C
C evaluate Lennard Jones terms multiplied by switching functions
C
C                 A       B
C      E    = (  ---  -  ---  )  *  Sw(R)
C       VdW        12      6
C                 R       R
C
C where the switching function Sw(R) is defined by
C
C
C            (                0                          R > R
C            (                                                off
C            (
C            (     2    2  2     2   2       2    2
C            (  ( R  - R  ) * ( R - R  - 3( R  - R  ) )
C            (     off          off          on
C     Sw(R)= (  -------------------------------------
C            (                2     2   3
C            (              (R   - R   )
C            (                off   on
C            (
C            (
C            (                1                          R < R
C                                                             on
C
      SS=MIN(MAX(S,C2ONNB),C2OFNB)
      RIJL=C2ONNB-SS
      RIJU=C2OFNB-SS
      FUNCT=RIJU**2*(RIJU-THREE*RIJL)*RUL3
      DFN=RIJL*RIJU*RUL12
C
      SIG6=R2*R2*R2
      SIG12=CNBAL(J)*SIG6*SIG6
      SIG6=CNBBL(J)*SIG6
      EN=SIG12-SIG6
      EVDWK=EN*FUNCT
      DEVDWK=EN*DFN-WSIX*R2*(EN+SIG12)*FUNCT
C
C evaluate force shifted electrostatic potential
C                               b
C            C        1        R         b+1
C     Q Q  -----   ( ---- + -------- - --------  )
C      i j  EPS       R         b+1      b R
C                            b R            off
C                               off
C
C for R less than CTOFNB and zero otherwise
C
      G2=MAX(0,INT(TWO+SQS*CROFNB))
      G1=RSQS+CCSHFO*S*S+DDSHFO
      G3=CGIL*CGL(J)
      ELECK=G3*G1*G2
      DELECK=G2*G3*(FOUR*CCSHFO*S-RSQS*R2)*ESCALE
C
C forces (including PIG's weights)
      DELTX=DELTX*(DEVDWK+DELECK)
      DELTY=DELTY*(DEVDWK+DELECK)
      DELTZ=DELTZ*(DEVDWK+DELECK)
C
      DX(I)=DX(I)+DELTX
      DY(I)=DY(I)+DELTY
      DZ(I)=DZ(I)+DELTZ
C
      XJL(J)=DELTX
      YJL(J)=DELTY
      ZJL(J)=DELTZ
C
C local accumulation
      EVDWL=EVDWL+EVDWK
      ELECL=ELECL+ELECK
      END DO
C
C global accumulation
      EVDW=EVDW+EVDWL*VSCALE
      ELEC=ELEC+ELECL*ESCALE
      EVDWV=EVDWV+EVDWL*VVWGHT
      ELECV=ELECV+ELECL*EVWGHT
C
C transform "j" image forces
      IF (QFORC) CALL XPFORC(N,XJL,YJL,ZJL)
C accumulate "j" forces
C >>>>>>ignore vector dependencies<<<<<<<
C$DIR NO_RECURRENCE
CDIR$ IVDEP
C*$*ASSERT PERMUTATION (JND)
      DO J=1,N
      DX(JND(J))=DX(JND(J))-XJL(J)
      DY(JND(J))=DY(JND(J))-YJL(J)
      DZ(JND(J))=DZ(JND(J))-ZJL(J)
      END DO
C
      END IF
C
C end of main loop
      END DO
C
      RETURN
      END
C====================================================================
      SUBROUTINE ENBFSB(EVDW,ELEC,EVDWV,ELECV,NATOM,JNB,LISTL,CG,
     &         MAXCN,CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,
     &         EVWGHT,VVWGHT,QIMAG,QFORC,XJ,YJ,ZJ,JND,CGL,
     &         CNBAL,CNBBL,CNBIL,XJL,YJL,ZJL)
C
C Features:
C    Lennard-Jones and free electrostatic potential
C    atom by atom switching Lennard-Jones and force shifting electrostatics
C    as function of B in Eq 15 J Comp Chem 15 (1994) 667
C
C
C Author: Axel T. Brunger
C modified by Alexandre Bonvin Yale, 11-SEP-95
C =======================================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'nbonds.inc'
      DOUBLE PRECISION EVDW, ELEC
      DOUBLE PRECISION EVDWV, ELECV
      DOUBLE PRECISION EVWGHT, VVWGHT
      INTEGER NATOM, JNB(*), LISTL(4)
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
C local
      INTEGER I, J, N
      DOUBLE PRECISION EN, EVDWL, EVDWK, ELECL, ELECK, DEVDWK, DELECK
      DOUBLE PRECISION R2, SIG6, SIG12, CGIL, CGF
      DOUBLE PRECISION G1, G3, C2OFNB, C2ROF2, CHROF2, CROFNB, CHROFN
      DOUBLE PRECISION C2ONNB, RUL3, RUL12, S, SS,  RIJL, RIJU, FUNCT
      DOUBLE PRECISION DELTX, DELTY, DELTZ, DFN, SQS, RSQS
C      DOUBLE PRECISION BSHFO2, BSHFOM2
      DOUBLE PRECISION WSIX, CCSHFO, DDSHFO, BSHFOP, BSHFOM, G2
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, SIX, TWELVE
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0)
      PARAMETER (SIX=6.D0, TWELVE=12.D0)
C begin
C
      WSIX=SIX*VSCALE
      CGF=CCELEC/EPS
      C2OFNB=CTOFNB*CTOFNB
      CROFNB=-ONE/CTOFNB
      BSHFOP=BSHFO+ONE
      BSHFOM=BSHFO-TWO
C      BSHFO2=BSHFO/TWO
C      BSHFOM2=BSHFOM/TWO
      CCSHFO=ONE/(BSHFO*CTOFNB**BSHFOP)
      DDSHFO=BSHFOP*CROFNB/BSHFO
      CHROFN=-TWO/CTOFNB
      C2ONNB=CTONNB*CTONNB
      RUL3=ONE/(C2OFNB-C2ONNB)**3
      RUL12=RUL3*TWELVE*VSCALE
      C2ROF2=-ONE/C2OFNB
      CHROF2=-FOUR/C2OFNB
C main loop over atom I
      DO I=1,NATOM
C
      N=JNB(I)
      IF (N.GT.0) THEN
C
C define local index list JND
      CALL GETNB(N,JND,LISTL)
C
C gather Lennard-Jones A, B and charge CG
C gather "j" coordinates
      DO J=1,N
      CNBBL(J)=CNBB(LOOKUP(JND(J)),LOOKUP(I))
      CNBAL(J)=CNBA(LOOKUP(JND(J)),LOOKUP(I))
      CNBIL(J)=(INHIBT*CNBVR(LOOKUP(JND(J)),LOOKUP(I)))**2
      END DO
      DO J=1,N
      CGL(J)=CG(JND(J))
      XJL(J)=XJ(JND(J))
      YJL(J)=YJ(JND(J))
      ZJL(J)=ZJ(JND(J))
      END DO
C
C appropriately scale self-interactions
      IF (JND(1).EQ.I) THEN
      CNBBL(1)=CNBBL(1)/TWO
      CNBAL(1)=CNBAL(1)/TWO
      CGL(1)=CGL(1)/TWO
      END IF
CC
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
      ELECL=ZERO
      EVDWL=ZERO
      CGIL=CG(I)*CGF
C
C compute Lennard Jones and electrostatic expression
      DO J=1,N
C
      DELTX=XJL(J)
      DELTY=YJL(J)
      DELTZ=ZJL(J)
      S=MAX(CNBIL(J),DELTX**2+DELTY**2+DELTZ**2)
      SQS=SQRT(S)
      RSQS=ONE/SQS
      R2=ONE/S
C
C evaluate Lennard Jones terms multiplied by switching functions
C
C                 A       B
C      E    = (  ---  -  ---  )  *  Sw(R)
C       VdW        12      6
C                 R       R
C
C where the switching function Sw(R) is defined by
C
C
C            (                0                          R > R
C            (                                                off
C            (
C            (     2    2  2     2   2       2    2
C            (  ( R  - R  ) * ( R - R  - 3( R  - R  ) )
C            (     off          off          on
C     Sw(R)= (  -------------------------------------
C            (                2     2   3
C            (              (R   - R   )
C            (                off   on
C            (
C            (
C            (                1                          R < R
C                                                             on
C
      SS=MIN(MAX(S,C2ONNB),C2OFNB)
      RIJL=C2ONNB-SS
      RIJU=C2OFNB-SS
      FUNCT=RIJU**2*(RIJU-THREE*RIJL)*RUL3
      DFN=RIJL*RIJU*RUL12
C
      SIG6=R2*R2*R2
      SIG12=CNBAL(J)*SIG6*SIG6
      SIG6=CNBBL(J)*SIG6
      EN=SIG12-SIG6
      EVDWK=EN*FUNCT
      DEVDWK=EN*DFN-WSIX*R2*(EN+SIG12)*FUNCT
C
C evaluate force shifted electrostatic potential
C                               b
C            C        1        R         b+1
C     Q Q  -----   ( ---- + -------- - --------  )
C      i j  EPS       R         b+1      b R
C                            b R            off
C                               off
C
C for R less than CTOFNB and zero otherwise
C
      G2=MAX(0,INT(TWO+SQS*CROFNB))
C      G1=RSQS+CCSHFO*S**BSHFO2+DDSHFO
      G1=RSQS+CCSHFO*SQS**BSHFO+DDSHFO
      G3=CGIL*CGL(J)
      ELECK=G3*G1*G2
C      DELECK=G2*G3*(BSHFO*CCSHFO*S**BSHFOM2-RSQS*R2)*ESCALE
      DELECK=G2*G3*(BSHFO*CCSHFO*SQS**BSHFOM-RSQS*R2)*ESCALE
C
C forces (including PIG's weights)
      DELTX=DELTX*(DEVDWK+DELECK)
      DELTY=DELTY*(DEVDWK+DELECK)
      DELTZ=DELTZ*(DEVDWK+DELECK)
C
      DX(I)=DX(I)+DELTX
      DY(I)=DY(I)+DELTY
      DZ(I)=DZ(I)+DELTZ
C
      XJL(J)=DELTX
      YJL(J)=DELTY
      ZJL(J)=DELTZ
C
C local accumulation
      EVDWL=EVDWL+EVDWK
      ELECL=ELECL+ELECK
      END DO
C
C global accumulation
      EVDW=EVDW+EVDWL*VSCALE
      ELEC=ELEC+ELECL*ESCALE
      EVDWV=EVDWV+EVDWL*VVWGHT
      ELECV=ELECV+ELECL*EVWGHT
C
C transform "j" image forces
      IF (QFORC) CALL XPFORC(N,XJL,YJL,ZJL)
C accumulate "j" forces
C >>>>>>ignore vector dependencies<<<<<<<
C$DIR NO_RECURRENCE
CDIR$ IVDEP
C*$*ASSERT PERMUTATION (JND)
      DO J=1,N
      DX(JND(J))=DX(JND(J))-XJL(J)
      DY(JND(J))=DY(JND(J))-YJL(J)
      DZ(JND(J))=DZ(JND(J))-ZJL(J)
      END DO
C
      END IF
C
C end of main loop
      END DO
C
      RETURN
      END
C==================================================================
      SUBROUTINE ENBRD(EVDW,ELEC,EVDWV,ELECV,NATOM,JNB,LISTL,CG,
     &           MAXCN,CNBA,CNBB,CNBVR,LOOKUP,ESCALE,VSCALE,
     &           EVWGHT,VVWGHT,QIMAG,QFORC,XJ,YJ,ZJ,
     &           JND,CGL,CNBAL,CNBBL,CNBIL,XJL,YJL,ZJL)
C
C Features:
C    Lennard-Jones and 1/r-dependent dielectric electrostatic
C    potential. Atom by atom switching Lennard-Jones and
C    electrostatics.
C
C Author:  Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'nbonds.inc'
      DOUBLE PRECISION EVDW, ELEC
      DOUBLE PRECISION EVDWV, ELECV
      DOUBLE PRECISION EVWGHT, VVWGHT
      INTEGER NATOM, JNB(*), LISTL(4)
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
C local
      INTEGER J
      DOUBLE PRECISION EN
      DOUBLE PRECISION EVDWL, EVDWK, ELECL, ELECK, DEVDWK, DELECK
      DOUBLE PRECISION R2, SIG6, SIG12, CGIL, CGF
      DOUBLE PRECISION G3, C2OFNB
      DOUBLE PRECISION S, RIJL, RIJU, FUNCT
      DOUBLE PRECISION DELTX, DELTY, DELTZ, DFN
      DOUBLE PRECISION C2ONNB, RUL3, RUL12
      INTEGER I, N
      DOUBLE PRECISION WSIX
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, SIX, TWELVE
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0)
      PARAMETER (TWELVE=12.D0, SIX=6.D0)
C begin
C
      WSIX=SIX*VSCALE
      CGF=CCELEC/EPS
      C2OFNB=CTOFNB*CTOFNB
      C2ONNB=CTONNB*CTONNB
      RUL3=ONE/(C2OFNB-C2ONNB)**3
      RUL12=RUL3*TWELVE*VSCALE
C
C main loop over atom I
      DO I=1,NATOM
C
      N=JNB(I)
      IF (N.GT.0) THEN
C
C define local index list JND
      CALL GETNB(N,JND,LISTL)
C
C gather Lennard-Jones A, B and charge CG
C gather "j" coordinates
      DO J=1,N
      CNBBL(J)=CNBB(LOOKUP(JND(J)),LOOKUP(I))
      CNBAL(J)=CNBA(LOOKUP(JND(J)),LOOKUP(I))
      CNBIL(J)=(INHIBT*CNBVR(LOOKUP(JND(J)),(LOOKUP(I))))**2
      END DO
      DO J=1,N
      CGL(J)=CG(JND(J))
      XJL(J)=XJ(JND(J))
      YJL(J)=YJ(JND(J))
      ZJL(J)=ZJ(JND(J))
      END DO
C
C appropriately scale self-interactions
      IF (JND(1).EQ.I) THEN
      CNBBL(1)=CNBBL(1)/TWO
      CNBAL(1)=CNBAL(1)/TWO
      CGL(1)=CGL(1)/TWO
      END IF
CC
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
      ELECL=ZERO
      EVDWL=ZERO
      CGIL=CG(I)*CGF
C
C compute Lennard Jones and electrostatic expression
      DO J=1,N
C
      DELTX=XJL(J)
      DELTY=YJL(J)
      DELTZ=ZJL(J)
      S=MAX(CNBIL(J),DELTX**2+DELTY**2+DELTZ**2)
      R2=ONE/S
C
C evaluate Lennard Jones terms multiplied by switching functions
C
C                 A       B
C      E    = (  ---  -  ---  )  *  Sw(R)
C       VdW        12      6
C                 R       R
C
C where the switching function Sw(R) is defined by
C
C
C
C            (                0                          R > R
C            (                                                off
C            (
C            (     2    2  2     2   2       2    2
C            (  ( R  - R  ) * ( R - R  - 3( R  - R  ) )
C            (     off          off          on
C     Sw(R)= (  -------------------------------------
C            (                2     2   3
C            (              (R   - R   )
C            (                off   on
C            (
C            (
C            (                1                          R < R
C                                                             on
C
      FUNCT=MIN(MAX(S,C2ONNB),C2OFNB)
      RIJL=C2ONNB-FUNCT
      RIJU=C2OFNB-FUNCT
      FUNCT=RIJU**2*(RIJU-THREE*RIJL)*RUL3
      DFN=RIJL*RIJU*RUL12
C
      SIG6=R2*R2*R2
      SIG12=CNBAL(J)*SIG6*SIG6
      SIG6=CNBBL(J)*SIG6
      EN=SIG12-SIG6
      EVDWK=EN*FUNCT
      DEVDWK=EN*DFN-WSIX*R2*(EN+SIG12)*FUNCT
C
C evaluate switched , 1/r dielectric electrostatic potential
C
C            C
C
C     Q Q  -----     * SW(R)
C      i j       2
C          EPS R
C
C where SW(R) is the same as above.
C
C
      G3=CGIL*CGL(J)*R2
      ELECK=G3*FUNCT
      DELECK=(-TWO*R2*ELECK+DFN*G3)*ESCALE
C
C forces
      DELTX=DELTX*(DEVDWK+DELECK)
      DELTY=DELTY*(DEVDWK+DELECK)
      DELTZ=DELTZ*(DEVDWK+DELECK)
C
      DX(I)=DX(I)+DELTX
      DY(I)=DY(I)+DELTY
      DZ(I)=DZ(I)+DELTZ
C
      XJL(J)=DELTX
      YJL(J)=DELTY
      ZJL(J)=DELTZ
C
C local accumulation
      EVDWL=EVDWL+EVDWK
      ELECL=ELECL+ELECK
      END DO
C
C global accumulation
      EVDW=EVDW+EVDWL*VSCALE
      ELEC=ELEC+ELECL*ESCALE
      EVDWV=EVDWV+EVDWL*VVWGHT
      ELECV=ELECV+ELECL*EVWGHT
C
C transform "j" image forces
      IF (QFORC) CALL XPFORC(N,XJL,YJL,ZJL)
C accumulate "j" forces
C >>>>>>ignore vector dependencies<<<<<<<
C$DIR NO_RECURRENCE
CDIR$ IVDEP
C*$*ASSERT PERMUTATION (JND)
      DO J=1,N
      DX(JND(J))=DX(JND(J))-XJL(J)
      DY(JND(J))=DY(JND(J))-YJL(J)
      DZ(JND(J))=DZ(JND(J))-ZJL(J)
      END DO
C
      END IF
C
C end of main loop
      END DO
C
      RETURN
      END
C==================================================================
      SUBROUTINE ENBREP(EVDW,ELEC,EVDWV,ELECV,NATOM,JNB,LISTL,CG,
     &              MAXCN,CNBVR,LOOKUP,ESCALE,VSCALE,EVWGHT,VVWGHT,
     &              QIMAG,QFORC,XJ,YJ,ZJ,
     &              JND,CGL,CNBAL,CNBBL,CNBIL,XJL,YJL,ZJL)
C
C Features:
C           Harmonic vdW repel term, no electrostatic potential
C
C Authors: Michael Nilges and  Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'nbonds.inc'
      DOUBLE PRECISION EVDW, ELEC
      DOUBLE PRECISION EVDWV, ELECV
      DOUBLE PRECISION EVWGHT, VVWGHT
      INTEGER NATOM, JNB(*), LISTL(4)
      DOUBLE PRECISION CG(*)
      INTEGER MAXCN
      DOUBLE PRECISION CNBVR(MAXCN,MAXCN)
      INTEGER LOOKUP(*)
      DOUBLE PRECISION ESCALE, VSCALE
      LOGICAL QIMAG, QFORC
      DOUBLE PRECISION XJ(*), YJ(*), ZJ(*)
      INTEGER JND(*)
      DOUBLE PRECISION CGL(*), CNBAL(*), CNBBL(*), CNBIL(*), XJL(*)
      DOUBLE PRECISION YJL(*), ZJL(*)
C local
      INTEGER I, J, N
      DOUBLE PRECISION EVDWL, EVDWK, DEVDWK
      DOUBLE PRECISION S
      DOUBLE PRECISION DELTX, DELTY, DELTZ, RDIF, RDIF2
      DOUBLE PRECISION CREP2, CREP8, CREP4
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, SIX
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, FOUR=4.D0, SIX=6.D0)
C begin
      CREP2=TWO*CREP*VSCALE
      CREP4=TWO*CREP2
      CREP8=TWO*CREP4
C
C check things
      IF (EXIREP.NE.1.AND.EXIREP.NE.2) THEN
      WRITE(6,'(A)') ' %ENBOND-ERR: EXIREP can only be 1 or 2'
      CALL DIE
      END IF
      IF (EXPREP.NE.2.AND.EXPREP.NE.4) THEN
      WRITE(6,'(A)') ' %ENBOND-ERR: EXPREP can only be 2 or 4'
      CALL DIE
      END IF
C
C main loop over atom I
      DO I=1,NATOM
C
      N=JNB(I)
      IF (N.GT.0) THEN
C
C define local index list JND
      CALL GETNB(N,JND,LISTL)
C
C gather vdW radii and store in vdW - "A" parameter array
      IF (EXIREP.EQ.2) THEN
      DO J=1,N
      CNBAL(J)=(FREPEL*CNBVR(LOOKUP(JND(J)),(LOOKUP(I))))**2
      CNBIL(J)=(INHIBT*CNBVR(LOOKUP(JND(J)),(LOOKUP(I))))**2
      END DO
      ELSE
      DO J=1,N
      CNBAL(J)=FREPEL*CNBVR(LOOKUP(JND(J)),(LOOKUP(I)))
      CNBIL(J)=INHIBT*CNBVR(LOOKUP(JND(J)),(LOOKUP(I)))
      END DO
      END IF
C
C appropriately scale self-interactions
      IF (JND(1).EQ.I) THEN
      CNBIL(1)=CNBIL(1)/TWO
      CNBAL(1)=CNBAL(1)/TWO
      END IF
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
      EVDWL=ZERO
C
C
C compute repel term
C                      EXIREP     EXIREP  EXPREP
C     E   =  CREP*(RMIN      -   R       ) ,
C      VdW
C
C only exprep=2 and exprep=4 supported to speed things up
C only exirep=1 and exirep=2 supported to speed things up
C
      IF (EXPREP.EQ.2) THEN
      IF (EXIREP.EQ.2) THEN
C
      DO J=1,N
      DELTX=XJL(J)
      DELTY=YJL(J)
      DELTZ=ZJL(J)
C
      S=MAX(CNBIL(J),DELTX**2+DELTY**2+DELTZ**2)
      RDIF=MAX(ZERO,CNBAL(J)-S)
      DEVDWK=-CREP4*RDIF
      EVDWK=CREP*RDIF**2
C
C
C forces
      DELTX=DELTX*DEVDWK
      DELTY=DELTY*DEVDWK
      DELTZ=DELTZ*DEVDWK
C
      DX(I)=DX(I)+DELTX
      DY(I)=DY(I)+DELTY
      DZ(I)=DZ(I)+DELTZ
C
      XJL(J)=DELTX
      YJL(J)=DELTY
      ZJL(J)=DELTZ
C
C local accumulation
      EVDWL=EVDWL+EVDWK
      END DO
C
      ELSE IF (EXIREP.EQ.1) THEN
C
      DO J=1,N
      DELTX=XJL(J)
      DELTY=YJL(J)
      DELTZ=ZJL(J)
C
      S=MAX(CNBIL(J),SQRT(DELTX**2+DELTY**2+DELTZ**2))
      RDIF=MAX(ZERO,CNBAL(J)-S)
      DEVDWK=-CREP2*RDIF/S
      EVDWK=CREP*RDIF**2
C
C
C forces
      DELTX=DELTX*DEVDWK
      DELTY=DELTY*DEVDWK
      DELTZ=DELTZ*DEVDWK
C
      DX(I)=DX(I)+DELTX
      DY(I)=DY(I)+DELTY
      DZ(I)=DZ(I)+DELTZ
C
      XJL(J)=DELTX
      YJL(J)=DELTY
      ZJL(J)=DELTZ
C
C local accumulation
      EVDWL=EVDWL+EVDWK
      END DO
C
      END IF
      ELSE IF (EXPREP .EQ. 4) THEN
      IF (EXIREP .EQ. 2) THEN
C
      DO J=1,N
      DELTX=XJL(J)
      DELTY=YJL(J)
      DELTZ=ZJL(J)
C
      S=MAX(CNBIL(J),DELTX**2+DELTY**2+DELTZ**2)
      RDIF=MAX(ZERO,CNBAL(J)-S)
      RDIF2=RDIF*RDIF
      EVDWK=CREP*RDIF2*RDIF2
      DEVDWK=-CREP8*RDIF2*RDIF
C
C
C forces
      DELTX=DELTX*DEVDWK
      DELTY=DELTY*DEVDWK
      DELTZ=DELTZ*DEVDWK
C
      DX(I)=DX(I)+DELTX
      DY(I)=DY(I)+DELTY
      DZ(I)=DZ(I)+DELTZ
C
      XJL(J)=DELTX
      YJL(J)=DELTY
      ZJL(J)=DELTZ
C
C local accumulation
      EVDWL=EVDWL+EVDWK
      END DO
      ELSE IF (EXIREP .EQ. 1) THEN
C
      DO J=1,N
      DELTX=XJL(J)
      DELTY=YJL(J)
      DELTZ=ZJL(J)
C
      S=MAX(CNBIL(J),SQRT(DELTX**2+DELTY**2+DELTZ**2))
      RDIF=MAX(ZERO,CNBAL(J)-S)
      RDIF2=RDIF*RDIF
      EVDWK=CREP*RDIF2*RDIF2
      DEVDWK=-CREP4*RDIF2*RDIF/S
C
C
C forces
      DELTX=DELTX*DEVDWK
      DELTY=DELTY*DEVDWK
      DELTZ=DELTZ*DEVDWK
C
      DX(I)=DX(I)+DELTX
      DY(I)=DY(I)+DELTY
      DZ(I)=DZ(I)+DELTZ
C
      XJL(J)=DELTX
      YJL(J)=DELTY
      ZJL(J)=DELTZ
C
C local accumulation
      EVDWL=EVDWL+EVDWK
      END DO
      END IF
      END IF
C
C
C
C global accumulation
      EVDW=EVDW+EVDWL*VSCALE
      EVDWV=EVDWV+EVDWL*VVWGHT
C
C transform "j" image forces
      IF (QFORC) CALL XPFORC(N,XJL,YJL,ZJL)
C accumulate "j" forces
C >>>>>>ignore vector dependencies<<<<<<<
C$DIR NO_RECURRENCE
CDIR$ IVDEP
C*$*ASSERT PERMUTATION (JND)
      DO J=1,N
      DX(JND(J))=DX(JND(J))-XJL(J)
      DY(JND(J))=DY(JND(J))-YJL(J)
      DZ(JND(J))=DZ(JND(J))-ZJL(J)
      END DO
C
      END IF
C
C end of main loop
      END DO
C
      RETURN
      END
C==================================================================
      SUBROUTINE ENBAR(EVDW,ELEC,EVDWV,ELECV,NATOM,JNB,LISTL,CG,
     &              MAXCN,cnbb,cnba,LOOKUP,ESCALE,VSCALE,EVWGHT,
     &              VVWGHT,QIMAG,QFORC,XJ,YJ,ZJ,
     &              JND,CGL,CNBAL,CNBBL,CNBIL,XJL,YJL,ZJL)
C
C Non-bonding routine added by JJK
C
C Features:
C           Harmonic vdW attractive-repulsive term
C
C Authors: John Kuszewski, Michael Nilges, Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'nbonds.inc'
      DOUBLE PRECISION EVDW, ELEC
      DOUBLE PRECISION EVDWV, ELECV
      DOUBLE PRECISION EVWGHT, VVWGHT
      INTEGER NATOM, JNB(*), LISTL(4)
      DOUBLE PRECISION CG(*)
      INTEGER MAXCN
      DOUBLE PRECISION cnbb(maxcn,maxcn), cnba(maxcn,maxcn)
      INTEGER LOOKUP(*)
      DOUBLE PRECISION ESCALE, VSCALE
      LOGICAL QIMAG, QFORC
      DOUBLE PRECISION XJ(*), YJ(*), ZJ(*)
      INTEGER JND(*)
      DOUBLE PRECISION CGL(*), CNBAL(*), CNBBL(*), CNBIL(*), XJL(*)
      DOUBLE PRECISION YJL(*), ZJL(*)
C local
      INTEGER I, J, N
      DOUBLE PRECISION EVDWL, EVDWK
      DOUBLE PRECISION CREP2, CREP8, CREP4
C
      DOUBLE PRECISION R, KREPEL
      DOUBLE PRECISION O1, O2, O3, O4, O6, O8, O9, O10
      DOUBLE PRECISION DEVDWDR, DEVDWDDELTAX, DEVDWDDELTAY
      DOUBLE PRECISION DEVDWDDELTAZ, A, B, EPSILON, SIGMA
      DOUBLE PRECISION SIZE, RATTRACT, KREP
      DOUBLE PRECISION DELTAX, DELTAY, DELTAZ
      LOGICAL O5, O7
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, SIX
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, FOUR=4.D0, SIX=6.D0)
      DOUBLE PRECISION SIXTEEN
      PARAMETER (SIXTEEN=16.0D0)
C
C begin
C
      CREP2=TWO*CREP*VSCALE
      CREP4=TWO*CREP2
      CREP8=TWO*CREP4
C
C check things
      IF (EXIREP.NE.2) THEN
      WRITE(6,'(A)') ' %ENBOND-ERR: EXIREP can only be 2'
      CALL DIE
      END IF
      IF (EXPREP.NE.2) THEN
      WRITE(6,'(A)') ' %ENBOND-ERR: EXPREP can only be 2'
      CALL DIE
      END IF
C
C main loop over atom I
C
      DO I=1,NATOM
C
      N=JNB(I)
      IF (N.GT.0) THEN
C
C define local index list JND
C
      CALL GETNB(N,JND,LISTL)
C
C gather "j" coordinates
C
      DO J=1,N
      XJL(J)=XJ(JND(J))
      YJL(J)=YJ(JND(J))
      ZJL(J)=ZJ(JND(J))
      END DO
C
C compute differences to X(I), Y(I), Z(I)
C
      IF (QIMAG) THEN
      CALL XPIMAG(X(I),Y(I),Z(I),1,N,XJL,YJL,ZJL,
     &            XJL,YJL,ZJL,.FALSE.,0)
      ELSE
      DO J=1,N
      XJL(J)=X(I)-XJL(J)
      YJL(J)=Y(I)-YJL(J)
      ZJL(J)=Z(I)-ZJL(J)
      END DO
      END IF
C
C zero out total energy
C
      EVDWL=ZERO
      KREPEL = CREP * VSCALE
C
C compute attractive-repulsive term
C
      DO J=1,N
C
C get its distance from atom I
C
      DELTAX=XJL(J)
      DELTAY=YJL(J)
      DELTAZ=ZJL(J)
C
      R = SQRT(DELTAX**2 + DELTAY**2 + DELTAZ**2)
C
C take care of inhibition at small Rs
C
C Get the Lennard-Jones parameters from the
C L-J "a" and "b" arrays.
C
      A = CNBA(LOOKUP(I),LOOKUP(JND(J)))
      B = CNBB(LOOKUP(I),LOOKUP(JND(J)))
C
C appropriately scale self-interactions
      IF (JND(1).EQ.I) THEN
      A=A/TWO
      B=B/TWO
      END IF
      EPSILON = (B**2/SIXTEEN)/(A/FOUR)
      SIGMA = (B/(FOUR * EPSILON))**(ONE/SIX)
C
C get everything into order for the Mma-generated code
C
      SIZE = FREPEL
      RATTRACT = RATT
      KREP = KREPEL
C
C new attrep function from Mma
C
      O1 = (ONE/SIX)
      O2 = O1
      O3 = 2**O2
      O4 = O3*SIGMA*SIZE
      O5 = R.LT.O4
      O6 = O3*RATT*SIGMA*SIZE
      O7 = R.LT.O6
      O8 = O3*SIGMA
      O9 = -O4 + O8
      O10 = O6 - O8
C
C get energy
C
      IF (O5) THEN
      EVDWK =  KREP*(O4**2 - R**2)**2
      ELSE
      IF (O7) THEN
      IF (R.LT.O8) THEN
      EVDWK = -(DOW*EPSILON*(O9**2 - (-O8 + R)**2)**2/O9**4)
      ELSE
      EVDWK = -(DOW*EPSILON*(O10**2 - (-O8 +R)**2)**2/O10**4)
      END IF
      ELSE
      EVDWK =ZERO
      END IF
      END IF
C
C get derivative of energy wrt r
C
      IF (O5) THEN
      DEVDWDR = -4*KREP*R*(O4**2 - R**2)
      ELSE
      IF (O7) THEN
      IF (R.LT.O8) THEN
      DEVDWDR = FOUR*DOW*EPSILON*(-O8 + R)*(O9**2 -(-O8 + R)**2)/O9**4
      ELSE
      DEVDWDR = FOUR*DOW*EPSILON*(-O8 + R)*(O10**2 -(-O8 + R)**2)/O10**4
      END IF
      ELSE
      DEVDWDR =ZERO
      END IF
      END IF
C
C deal with self-interactions
C
      IF (JND(J).EQ.I) THEN
      EVDWK =ZERO
      DEVDWDR =ZERO
      END IF
C
C now get dE/dx, dE/dy, dE/dz
C from dE/dR by the chain rule
C
      DEVDWDDELTAX = DEVDWDR * TWO * DELTAX
      DEVDWDDELTAY = DEVDWDR * TWO * DELTAY
      DEVDWDDELTAZ = DEVDWDR * TWO * DELTAZ
C
C update the derivative arrays
C
      DX(I)=DX(I)+DEVDWDDELTAX
      DY(I)=DY(I)+DEVDWDDELTAY
      DZ(I)=DZ(I)+DEVDWDDELTAZ
C
C I have no idea why this is done
C
      XJL(J)=DEVDWDDELTAX
      YJL(J)=DEVDWDDELTAY
      ZJL(J)=DEVDWDDELTAZ
C
C local accumulation
C
      EVDWL=EVDWL+EVDWk
      END DO
C
C
C
C global accumulation
      EVDW=EVDW+EVDWL*VSCALE
      EVDWV=EVDWV+EVDWL*VVWGHT
C
C transform "j" image forces
      IF (QFORC) CALL XPFORC(N,XJL,YJL,ZJL)
C accumulate "j" forces
C >>>>>>ignore vector dependencies<<<<<<<
C$DIR NO_RECURRENCE
CDIR$ IVDEP
C*$*ASSERT PERMUTATION (JND)
      DO J=1,N
      DX(JND(J))=DX(JND(J))-XJL(J)
      DY(JND(J))=DY(JND(J))-YJL(J)
      DZ(JND(J))=DZ(JND(J))-ZJL(J)
      END DO
C
      END IF
C
C end of main loop
      END DO
C
      RETURN
      END
C================================================================
      SUBROUTINE XPFORC(N,DX,DY,DZ)
C
C routine transforms forces for "j" image atoms,
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'nbonds.inc'
      INTEGER N
      DOUBLE PRECISION DX(*), DY(*), DZ(*)
C local
      INTEGER I
      DOUBLE PRECISION XTEMP, YTEMP, ZTEMP
C begin
      DO I=1,N
      XTEMP=XPOPER(1,1)*DX(I)+XPOPER(2,1)*DY(I)+XPOPER(3,1)*DZ(I)
      YTEMP=XPOPER(1,2)*DX(I)+XPOPER(2,2)*DY(I)+XPOPER(3,2)*DZ(I)
      ZTEMP=XPOPER(1,3)*DX(I)+XPOPER(2,3)*DY(I)+XPOPER(3,3)*DZ(I)
      DX(I)=XTEMP
      DY(I)=YTEMP
      DZ(I)=ZTEMP
      END DO
      RETURN
      END
