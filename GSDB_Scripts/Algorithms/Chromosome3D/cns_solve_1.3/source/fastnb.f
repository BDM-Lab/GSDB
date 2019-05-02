      SUBROUTINE FNBSET(OPTION)
C
C Fast Nonbonded parsing routine - for grid based pair-list
C and PMTA interface
C
C Paul Adams 11/4/94
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'fastnb.inc'
      INCLUDE 'nbonds.inc'
      CHARACTER*(*) OPTION
C local
      INTEGER N
C
C init
      IF (OPTION.EQ.'INIT') THEN
      CALL FNBINI
C
C parsing
      ELSE IF (OPTION.EQ.'PARSE') THEN
      CALL PUSEND('FASTNB>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('FASTNB>')
      CALL MISCOM('FASTNB>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-fastnb')
C
      ELSE IF (WD(1:4).EQ.'GRID') THEN
      QFSTNB=.TRUE.
      ELSE IF (WD(1:4).EQ.'GROU') THEN
      QFSTNB=.FALSE.
      ELSE IF (WD(1:4).EQ.'COUL') THEN
      PMTAM='COUL'
      ELSE IF (WD(1:4).EQ.'MULT') THEN
      PMTAM='MULT'
      ELSE IF (WD(1:4).EQ.'TREE') THEN
      PMTAM='TREE'
      ELSE IF (WD(1:4).EQ.'THET') THEN
      CALL NEXTF('THETa=',PTHETA)
      ELSE IF (WD(1:4).EQ.'LEVE') THEN
      CALL NEXTI('LEVEls=',PLEVEL)
      ELSE IF (WD(1:4).EQ.'TERM') THEN
      CALL NEXTI('TERMs=',PTERMS)
      ELSE IF (WD(1:3).EQ.'FFT') THEN
      QPFFT=.TRUE.
      ELSE IF (WD(1:4).EQ.'DIRE') THEN
      QPFFT=.FALSE.
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      CALL FNBPRI
      ELSE
      CALL CHKEND('FASTNB>',DONE)
      END IF
C
      DO N=1,NPIGMAX
C update exclusion list, 1-4 interaction list, and initialize others
      UPNBEX(N)=.TRUE.
C update intra-molecular interaction list
      UPNBLS(N)=.TRUE.
C update symmetry-related interaction list
      UPNBSS(N)=.TRUE.
      ENDDO
C
      END IF
      END DO
      DONE=.FALSE.
C
      END IF
C
      RETURN
      END
C
C===================================================================
C
      SUBROUTINE FNBINI
C
C routine initializes information for fast nonbonded routines
C
C Paul Adams 11/4/94
C
      IMPLICIT NONE
C I/O
      INCLUDE 'fastnb.inc'
C begin
C pairlist options
      QFSTNB=.FALSE.
C PMTA options
      PMTAM='COUL'
      PTHETA=1.5D0
      PLEVEL=3
      PTERMS=8
      QPFFT=.TRUE.
C
      RETURN
      END
C
C===================================================================
C
      SUBROUTINE FNBPRI
C
C routine prints info about current fast nonbonded options
C
C Paul Adams 11/4/94
C
      IMPLICIT NONE
C I/O
      INCLUDE 'fastnb.inc'
C local
      CHARACTER*5 SGROUP
      CHARACTER*6 SPSUM
C begin
      IF (QFSTNB) THEN
      SGROUP='GRID '
      ELSE
      SGROUP='GROUp'
      END IF
      IF (QPFFT) THEN
      SPSUM='FFT   '
      ELSE
      SPSUM='DIRECT'
      END IF
C
      WRITE(6,'(2A)') ' -----fast-nonbonded-options-----------',
     & '--------------------'
      WRITE(6,'(3A)') ' | ',SGROUP,
     & '                                                  |'
      IF (PMTAM.NE.'COUL') THEN
      WRITE(6,'(2A)') ' -----PMTA-nonbonded-options-----------',
     & '--------------------'
      IF (PMTAM.EQ.'TREE') THEN
      WRITE(6,'(3A,F4.2,A)')
     & ' | ',PMTAM,' THETa=',PTHETA,
     & '                                        |'
      ELSE
      WRITE(6,'(3A)')
     & ' | ',PMTAM,
     & '                                                   |'
      END IF
      WRITE(6,'(A,I2,A,I2,3A)')
     & ' | LEVEls=',PLEVEL,' TERMs=',PTERMS,'   ',SPSUM,
     & '                            |'
      END IF
      WRITE(6,'(2A)') ' --------------------------------------',
     & '--------------------'
C
      RETURN
      END
C
C===================================================================
C
      SUBROUTINE FXPSEAR(N,NNPIG)
C
C Routine searches for nonbonded interactions.
C
C QSYMM vs QNCS mode (exclusive options):
C  QSYMM mode uses crystallographic symmetry operators to generate
C        images of the primary molecule stored in the structure file.
C        All symmetry related molecules are generated and the
C        distance vectors to the primary molecule are computed.
C        Then the minimum image convention is applied to all distance
C        vectors in the following way: The distance vector is
C        converted to fractional coordinates, then the standard
C        minimum image convention is applied, and then the distance
C        vector is converted back to orthogonal A coordinates.
C        Then the distance cutoff criterion is applied to the
C        distance vector.  The information about symmetry
C        operators and the unit-cell is communicated via xcrystal.inc.
C
C  QNCS mode uses non-crystallographic symmetry operators to generate
C       images of the primary molecule stored in the structure file.
C        All symmetry related molecules are generated and the
C        distance vectors to the primary molecule are computed.
C        Then the distance cutoff criterion is applied to the
C        distance vector.  The information about symmetry
C        operators is communicated via ncs.inc.
C
C The routine splits the interactions into intra-molecular interactions
C (JNB, LIST) and symmetry interactions (JSB, SLIST).  These lists
C correspond to the energy terms ELEC,VDW and PELE,PVDW, respectively.
C Depending on which energy terms are turned on, the appropriate lists
C are generated.  Only for intra-molecular interactions, the
C non-bonded exclusion list is checked.
C
C Intra-molecular interactions between primary molecule atoms
C are contained as "trivial" cases in both QSYMM and QNCS mode:
C the first symmetry operator is always the identity.
C The case that only intra-molecular interactions are required
C (PELE, PVDW false) is treated as a special QNCS case.
C
C GROUP vs. ATOM mode (exclusive options)
C  GROUP mode excludes or includes whole groups as specified in the
C  structure file.  Distances are checked between group centroids (after
C  applying the minimum image convention in QSYMM mode).  The standard
C  atom-atom interaction lists are generated for allowed group-group
C  interactions.  Then, the minimum image convention is applied to
C  individual atom-atom interactions, and they are put either in the
C  list of intra-molecular interactions (JNB, INBLO) or the list of symmetry
C  related interactions (JSB, ISBLO).  The minimum image convention
C  is used again during the actual energy calculation on an atom-by
C  atom basis.  This can cause problems for situations where one of the
C  unit-cell lenghts is less than twice the non-bonded cutoff.  There
C  is no work-around for this problem yet.
C
C ATOM mode excludes or includes individual atom-atom interactions.
C  This mode uses the groups to make the calculation more efficient
C  (rectangles around group-centroid criterion).   Note, that the
C  results are independent of the size of the groups!
C
C Modified to use grid based selection scheme at highest level.
C The grid method can be turned on and off with the GRID/NOGRid flag
C in nbonds. If the method is on but the system is not big enough to
C efficiently use the method it will automatically be turned off.
C The BOXES subroutine tries to partition the primary molecule into a series
C of cubes. The initial size of the cube is based on CUTNB and the maximum
C groupsize (max(XSIZ,YSIZ,ZSIZ). The final size of the cubes is dynamically
C determined so that the number of cubes does not exceed some preset maximum,
C and so that the ratio of groups to boxes is larger than some preset minimum.
C If the number of cubes is below a preset minimum the grid method is
C turned off.
C The coordinates of the centres of each box are stored. The symmetry
C operators and minimum image convention are applied to these coordinates.
C The program loops over groups i_{1} to i_{ngrp}, to maintain the sequential
C nature of the non-bonded pairlist. The distance between the box that contains
C group i and all other boxes j_{1} to j_{nboxes} is calculated. Only boxes
C which have centres within a certain cut-off value are selected. The groups
C within these selected boxes are then used in the usual group centre based
C search. If the grid option is turned off the normal method (groups j_{i} to
C j_{ngrp} is used - this vectorises better than a derivative of the
C grid based code.
C Paul Adams 16/4/93
C
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'fastnb.inc'
      INTEGER N, NNPIG
C local
      LOGICAL QSYMM, QINTRA, QNCS
      INTEGER NLEN, NSLEN
C pointers
      INTEGER QMIM
      INTEGER QTEST, QSTEST, XTEMP, YTEMP, ZTEMP, XF, YF, ZF
      INTEGER XJF, YJF, ZJF, INTERJ
      INTEGER XCENT, YCENT, ZCENT, XCENTF, YCENTF, ZCENTF
      INTEGER XSIZ, YSIZ, ZSIZ, SPOINT, INTERG, INTERC, ILOCAL
      INTEGER XBOX, YBOX, ZBOX, XBOXF, YBOXF, ZBOXF, BOXTEST, BOX
      INTEGER BOXSIZ, BOXMEM, BOXTEMP, XBTEMP, YBTEMP, ZBTEMP
C begin
      IF (NATOM.GT.0) THEN
C
      IF (TIMER.GT.1) CALL DSPCPU(' NBONDS: at entry')
C
C
C Determine whether we have to compute the intra or
C symmetry interactions (or
C non-crystallographic symmetry interactions).  If we are using strict
C non-crystallographic symmetry, and the packing terms are on, then
C crystallographic symmetry is automatically shut off, and
C the packing energies used for NCS related non-bonded interactions.
      IF (.NOT. LNCSST) THEN
      QSYMM=QENER(SSPVDW).OR.QENER(SSPELE)
      QNCS=.FALSE.
      ELSE
      QSYMM=.FALSE.
      QNCS=QENER(SSPVDW).OR.QENER(SSPELE)
      END IF
      QINTRA=QENER(SSVDW).OR.QENER(SSELEC)
C
C Determine lengths of some of the arrays
      IF (QINTRA) THEN
      NLEN=NATOM
      ELSE
      NLEN=0
      END IF
      IF (QSYMM .OR. QNCS) THEN
C Note that xtal packing and ncs images are mutually exclusive
      IF (QSYMM) THEN
      NSLEN=NATOM*XRNSYM
      ELSE
      NSLEN=NATOM*NNSYM
      END IF
      ELSE
      NSLEN=0
      END IF
C
C Initialization.
C ==============
      IF (UPNBEX(N)) THEN
      CALL RENBND(0,0,0,0,0,0,0,0,0,N)
C
      IF (QINTRA) THEN
C update exclusion list.
      WRITE(6,'(A,I2)')
     &' NBONDS: generating intra-molecular exclusion list with mode=',
     & NBXMOD
      CALL UPINB(N)
      IF (TIMER.GT.1) CALL DSPCPU(' NBONDS: after excl. ')
      END IF
C allocate space for various lists (the fragmentation of the
C linked lists is NATOM).  Note that
C for the last Pair of Interacting Groups (PIG), allocate the
C 'previous position' arrays.
      IF (N.EQ.NNPIG) THEN
      CALL RENBND(NLEN,NATOM,NSLEN,NATOM,NATOM,-1,-1,-1,-1,N)
      ELSE
      CALL RENBND(NLEN,NATOM,NSLEN,NATOM,-1,-1,-1,-1,-1,N)
      END IF
      END IF
C
C allocate temporary space for search routine
      QTEST=ALLHP(ILOGIC(NATOM))
      QSTEST=ALLHP(ILOGIC(NATOM))
      XTEMP=ALLHP(IREAL8(NATOM))
      YTEMP=ALLHP(IREAL8(NATOM))
      ZTEMP=ALLHP(IREAL8(NATOM))
      QMIM=ALLHP(INTEG4(NATOM))
      XF=ALLHP(IREAL8(NATOM))
      YF=ALLHP(IREAL8(NATOM))
      ZF=ALLHP(IREAL8(NATOM))
      XJF=ALLHP(IREAL8(NATOM))
      YJF=ALLHP(IREAL8(NATOM))
      ZJF=ALLHP(IREAL8(NATOM))
      INTERJ=ALLHP(INTEG4(NATOM))
      XCENT=ALLHP(IREAL8(NGRP))
      YCENT=ALLHP(IREAL8(NGRP))
      ZCENT=ALLHP(IREAL8(NGRP))
      XCENTF=ALLHP(IREAL8(NGRP))
      YCENTF=ALLHP(IREAL8(NGRP))
      ZCENTF=ALLHP(IREAL8(NGRP))
      XSIZ=ALLHP(IREAL8(NGRP))
      YSIZ=ALLHP(IREAL8(NGRP))
      ZSIZ=ALLHP(IREAL8(NGRP))
      SPOINT=ALLHP(INTEG4(NATOM))
      INTERG=ALLHP(INTEG4(NGRP))
      INTERC=ALLHP(INTEG4(NGRP))
      ILOCAL=ALLHP(INTEG4(NATOM))
      XBOX=ALLHP(IREAL8(MXBOX))
      YBOX=ALLHP(IREAL8(MXBOX))
      ZBOX=ALLHP(IREAL8(MXBOX))
      BOXTEST=ALLHP(ILOGIC(NGRP))
      XBOXF=ALLHP(IREAL8(MXBOX))
      YBOXF=ALLHP(IREAL8(MXBOX))
      ZBOXF=ALLHP(IREAL8(MXBOX))
      BOX=ALLHP(INTEG4(NGRP))
      BOXSIZ=ALLHP(INTEG4(MXBOX+1))
      BOXMEM=ALLHP(INTEG4(NGRP))
      BOXTEMP=ALLHP(INTEG4(NGRP))
      XBTEMP=ALLHP(IREAL8(NGRP))
      YBTEMP=ALLHP(IREAL8(NGRP))
      ZBTEMP=ALLHP(IREAL8(NGRP))
C
      CALL FXPSEA2(QSYMM,QINTRA,QNCS,HEAP(IJNB(N)),
     &            LIST(1,N),HEAP(IJSB(N)),NATOM,
     &            SLIST(1,N),HEAP(IIBLOE(N)),
     &            HEAP(IINBEX(N)),HEAP(QTEST),HEAP(QSTEST),
     &            HEAP(XTEMP),HEAP(YTEMP),HEAP(ZTEMP),HEAP(QMIM),
     &            HEAP(XF),HEAP(YF),HEAP(ZF),HEAP(XJF),HEAP(YJF),
     &            HEAP(ZJF),HEAP(INTERJ),
     &            HEAP(XCENT),HEAP(YCENT),HEAP(ZCENT),
     &            HEAP(XCENTF),HEAP(YCENTF),HEAP(ZCENTF),
     &            HEAP(XSIZ),HEAP(YSIZ),HEAP(ZSIZ),
     &            HEAP(XBOX),HEAP(YBOX),HEAP(ZBOX),HEAP(BOXTEST),
     &            HEAP(XBOXF),HEAP(YBOXF),HEAP(ZBOXF),HEAP(BOX),
     &            HEAP(BOXSIZ),HEAP(BOXMEM),
     &            HEAP(SPOINT),HEAP(BOXTEMP),
     &            HEAP(XBTEMP),HEAP(YBTEMP),HEAP(ZBTEMP),
     &            HEAP(INTERG),HEAP(INTERC),
     &            HEAP(IINTER(N)),NNNB(N),NNSB(N),
     &            CUTNB,WMIN,LGROUP,HEAP(ILOCAL),NBQSPC)
C
C free temporary space for search routine
      CALL FREHP(ZBTEMP,IREAL8(NGRP))
      CALL FREHP(YBTEMP,IREAL8(NGRP))
      CALL FREHP(XBTEMP,IREAL8(NGRP))
      CALL FREHP(BOXTEMP,INTEG4(NGRP))
      CALL FREHP(BOXMEM,INTEG4(NGRP))
      CALL FREHP(BOXSIZ,INTEG4(MXBOX+1))
      CALL FREHP(BOX,INTEG4(NGRP))
      CALL FREHP(ZBOXF,IREAL8(MXBOX))
      CALL FREHP(YBOXF,IREAL8(MXBOX))
      CALL FREHP(XBOXF,IREAL8(MXBOX))
      CALL FREHP(BOXTEST,ILOGIC(NGRP))
      CALL FREHP(ZBOX,IREAL8(MXBOX))
      CALL FREHP(YBOX,IREAL8(MXBOX))
      CALL FREHP(XBOX,IREAL8(MXBOX))
      CALL FREHP(ILOCAL,INTEG4(NATOM))
      CALL FREHP(INTERC,INTEG4(NGRP))
      CALL FREHP(INTERG,INTEG4(NGRP))
      CALL FREHP(SPOINT,INTEG4(NATOM))
      CALL FREHP(ZSIZ,IREAL8(NGRP))
      CALL FREHP(YSIZ,IREAL8(NGRP))
      CALL FREHP(XSIZ,IREAL8(NGRP))
      CALL FREHP(ZCENTF,IREAL8(NGRP))
      CALL FREHP(YCENTF,IREAL8(NGRP))
      CALL FREHP(XCENTF,IREAL8(NGRP))
      CALL FREHP(ZCENT,IREAL8(NGRP))
      CALL FREHP(YCENT,IREAL8(NGRP))
      CALL FREHP(XCENT,IREAL8(NGRP))
      CALL FREHP(INTERJ,INTEG4(NATOM))
      CALL FREHP(ZJF,IREAL8(NATOM))
      CALL FREHP(YJF,IREAL8(NATOM))
      CALL FREHP(XJF,IREAL8(NATOM))
      CALL FREHP(ZF,IREAL8(NATOM))
      CALL FREHP(YF,IREAL8(NATOM))
      CALL FREHP(XF,IREAL8(NATOM))
      CALL FREHP(QMIM,INTEG4(NATOM))
      CALL FREHP(ZTEMP,IREAL8(NATOM))
      CALL FREHP(YTEMP,IREAL8(NATOM))
      CALL FREHP(XTEMP,IREAL8(NATOM))
      CALL FREHP(QSTEST,ILOGIC(NATOM))
      CALL FREHP(QTEST,ILOGIC(NATOM))
C
      IF (TIMER.GT.1) CALL DSPCPU(' NBONDS: at exit')
C
C now update INBX, INBY, INBZ for next tolerance check
      IF (N.EQ.NNPIG)
     &  CALL COPYCO(NATOM,X,Y,Z,HEAP(INBX),HEAP(INBY),HEAP(INBZ))
C
      END IF
      RETURN
      END
C
      SUBROUTINE FXPSEA2(QSYMM,QINTRA,QNCS,JNB,LIST,
     &            JSB,LEN,SLIST,IBLOEX,INBEX,
     &            QTEST,QSTEST,XTEMP,YTEMP,ZTEMP,QMIM,
     &            XF,YF,ZF,XJF,YJF,ZJF,
     &            INTERJ,XCENT,YCENT,
     &            ZCENT,XCENTF,YCENTF,ZCENTF,XSIZ,YSIZ,ZSIZ,
     &            XBOX,YBOX,ZBOX,BOXTEST,
     &            XBOXF,YBOXF,ZBOXF,BOX,
     &            BOXSIZ,BOXMEM,
     &            SPOINT,BOXTEMP,
     &            XBTEMP,YBTEMP,ZBTEMP,
     &            INTERG,INTERC,
     &            INTERE,NNNB,NNSB,CUTNB,WMIN,LGROUP,
     &            ILOCAL,NBQSPC)
C
C see routine above
C
C Author: Axel T. Brunger
C
C modified to use a grid based selection scheme at the highest level
C Paul Adams 16/4/93
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      LOGICAL QSYMM, QINTRA, QNCS
      INTEGER JNB(*), LIST(4)
      INTEGER LEN, JSB(LEN, *), SLIST(4)
      INTEGER INBEX(*), IBLOEX(*)
      LOGICAL QTEST(*), QSTEST(*)
      DOUBLE PRECISION XTEMP(*), YTEMP(*), ZTEMP(*)
      INTEGER QMIM(*)
      DOUBLE PRECISION XF(*), YF(*), ZF(*), XJF(*), YJF(*), ZJF(*)
      INTEGER INTERJ(*)
      DOUBLE PRECISION XCENT(*), YCENT(*), ZCENT(*)
      DOUBLE PRECISION XCENTF(*), YCENTF(*), ZCENTF(*)
      DOUBLE PRECISION XSIZ(*), YSIZ(*), ZSIZ(*)
      DOUBLE PRECISION XBOX(*), YBOX(*), ZBOX(*)
      DOUBLE PRECISION XBOXF(*), YBOXF(*), ZBOXF(*)
      DOUBLE PRECISION XBTEMP(*), YBTEMP(*), ZBTEMP(*)
      INTEGER BOX(*), BOXSIZ(*), BOXMEM(*)
      LOGICAL BOXTEST(*)
      INTEGER SPOINT(*), BOXTEMP(*), INTERG(*), INTERC(*)
      INTEGER INTERE(*),NNNB,NNSB
      DOUBLE PRECISION CUTNB,WMIN
      LOGICAL LGROUP
      INTEGER ILOCAL(*)
      DOUBLE PRECISION NBQSPC
C local
      DOUBLE PRECISION CUTNB2, WMIN2, BOXCUT, BOXCUT2, GRDMMX
      DOUBLE PRECISION DBPREC
      DOUBLE COMPLEX DBCOMP
      LOGICAL CHCK, QBOXL
      INTEGER ISYM, I, IAT, J, JJ, NBSYM, NXI, NXIMAX
      INTEGER JSYM, NNIAT, III
      INTEGER NBOXES, IBOX, K, JGRP, NVIOL
C parameter
      DOUBLE PRECISION ZERO, EEPS, PADDNG, LARGE, ONE, TWO
      PARAMETER (ZERO=0.0D0, EEPS=0.1D0, PADDNG=3)
      PARAMETER (LARGE=9999.0D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
C begin
C
      NVIOL=0
      QBOXL=.TRUE.
      WMIN2=WMIN**2
      CUTNB2=CUTNB**2
C
C set initial number of interactions to zero
      NNSB=0
      NNNB=0
C
C reset the nonbonded linked list pointers
      CALL RESNB(LIST)
      CALL RESNB(SLIST)
C
C compute group centers, size and interaction criterium
C for all groups of asymmetric unit
      CALL XPGRUP(X,Y,Z,XCENT,YCENT,ZCENT,XSIZ,YSIZ,ZSIZ,INTERG,INTERC,
     &            INTERE)
C
C
C calculate box cutoff size
C get the maximum group size
      CALL GRPRAD(NGRP,XSIZ,YSIZ,ZSIZ,GRDMMX,INTERG)
C
C grid method on
      BOXCUT=(SQRT(3.0D0)/2.0D0)*(CUTNB+2.0D0*GRDMMX)
C
C first make the boxes - return number of boxes and centres
C then assign groups to boxes  - return array of length ngrp with a
C pointer to the box number for each group.
      CALL BOXES(BOXCUT,NGRP,XCENT,YCENT,ZCENT,
     &     NBOXES,XBOX,YBOX,ZBOX,BOX,BOXSIZ,BOXMEM,QBOXL,
     &     INTERG)
C
C get number of symmetry operators which we want to use.
      IF (.NOT.QSYMM.AND..NOT.QNCS) THEN
C
C no symmetry interactions required: just use the first symmetry operator
C which is the identity.
      JSYM=1
      ELSE IF (QNCS) THEN
C
C NCS symmetry
      JSYM=NNSYM
      ELSE
C
C crystallographic symmetry
      JSYM=XRNSYM
      END IF
C
      DO ISYM=1,JSYM
C
      IF (QBOXL) THEN
      IF (ISYM.EQ.1) THEN
C
C if this is the primary molecule then use the minimum box centre cutoff
      BOXCUT2=((4.0D0/SQRT(3.0D0))*BOXCUT)**2
      ELSE
C
C if we are looking at symmetry molecules use larger box centre cutoff
C - needed because of the possible geometric orientation of the boxes
C after symmetry has been applied -- Apr4 '95 (factor two)
      BOXCUT2=((2.0D0*BOXCUT)+(CUTNB+2.0D0*SQRT(3.0D0)*GRDMMX))**2
      ENDIF
      ENDIF
C
      IF (QSYMM) THEN
C
C crystallographic-symmetry-mode:
C we have to apply
C the symmetry operators to the j-set atoms and groups.  --> return in
C XF, YF, ZF, XCENTF, YCENTF, ZCENTF (fractional coordinates!).
      CALL XPSYMO(NATOM,X,Y,Z,ISYM,XF,YF,ZF)
      CALL XPSYMO(NGRP,XCENT,YCENT,ZCENT,ISYM,XCENTF,YCENTF,ZCENTF)
C
      IF (QBOXL) THEN
C
C symmetry applied to box centres
      CALL XPSYMO(NBOXES,XBOX,YBOX,ZBOX,ISYM,XBOXF,YBOXF,ZBOXF)
      ENDIF
C
      ELSE
C
C NCS-symmetry mode (includes pure-intramolecular mode):
C we have to apply the NCS symmetry operator
C to the j-set atoms and groups.  --> return in XF, YF, ZF, XCENTF,
C YCENTF, ZCENTF (orthogonal coordinates!).
      CALL NCSXYZ(NATOM,X,Y,Z,ISYM,XF,YF,ZF)
      CALL NCSXYZ(NGRP,XCENT,YCENT,ZCENT,ISYM,XCENTF,YCENTF,ZCENTF)
C
      IF (QBOXL) THEN
C
C symmetry applied to box centres
      CALL NCSXYZ(NBOXES,XBOX,YBOX,ZBOX,ISYM,XBOXF,YBOXF,ZBOXF)
      ENDIF
C
      END IF
C
C always must do calculations for the very first box
      IBOX=-9999
C
C loop over i groups (we only have to search groups j=i,...,ngrp
C because of symmetry, i.e., if G[i] and S G{j} are interacting
C then by application of another symmetry operator we must have
C S' G[i] and G[j] interacting.
C This is still the main loop (rather than over boxes i_{1} to i_{nboxes})
C so that the pairlist is still sequential with respect to atom number
      DO I=1,NGRP
C
C check whether group contains atoms that interact at all
      IF (ABS(INTERG(I)).LT.+9999) THEN
C
C only do this code if we are using the box option and there are enough boxes
      IF (QBOXL) THEN
C
C if we are looking at a new box for group i then calculate interacting
C boxes (to give possible groups j). Otherwise we just use the old group
C list - as it is the same for each member of a box
      IF (BOX(I).NE.IBOX) THEN
C
C get box that group i is in
      IBOX=BOX(I)
C
C get distances between box i and j box centres
      IF (QSYMM) THEN
C
C crystallographic symmetry mode:
C compute minimum image distance between box centres i and all boxes j
         CALL XPIMAG(XBOX(IBOX),YBOX(IBOX),ZBOX(IBOX),1,NBOXES,
     &        XBOXF,YBOXF,ZBOXF,XTEMP,YTEMP,ZTEMP,.FALSE.,QMIM)
C
      ELSE
C
C NCS symmetry mode (includes pure intra-molecular mode):
C compute direct distance between box centres i and all boxes j
         DO J=1,NBOXES
            XTEMP(J)=XBOX(IBOX)-XBOXF(J)
            YTEMP(J)=YBOX(IBOX)-YBOXF(J)
            ZTEMP(J)=ZBOX(IBOX)-ZBOXF(J)
         END DO
      END IF
C
C squared distance between box centres
      DO J=1,NBOXES
         XBTEMP(J)=XTEMP(J)**2+YTEMP(J)**2+ZTEMP(J)**2
      END DO
C
C this should fill the boxtest flag array with trues for those groups j
C which are close enough to interact with group i
C note: all j are included (even j<i), we will just ignore these later
      DO J = 1,NGRP
         BOXTEST(J) = .FALSE.
      END DO
C
      DO J = 1,NBOXES
         IF (XBTEMP(J).LT.BOXCUT2) THEN
            DO K = BOXSIZ(J),(BOXSIZ(J+1)-1)
               BOXTEST(BOXMEM(K)) = .TRUE.
            END DO
         END IF
      END DO
C
      ENDIF
C
C now we pack the j groups for group i into smaller arrays
C boxtest i must always be true to comply with the common code later on
      BOXTEST(I)=.TRUE.
      JGRP=0
      DO J=I,NGRP
         IF (BOXTEST(J)) THEN
            JGRP=JGRP+1
            BOXTEMP(JGRP)=J
         END IF
      END DO
      DO J=1,JGRP
         XBTEMP(J)=XCENTF(BOXTEMP(J))
         YBTEMP(J)=YCENTF(BOXTEMP(J))
         ZBTEMP(J)=ZCENTF(BOXTEMP(J))
      END DO
C
C now check the distance between group i and all the groups j found in the
C box based search
      IF (QSYMM) THEN
C
C crystallographic symmetry mode:
C compute minimum image distance between groups i,j
         CALL XPIMAG(XCENT(I),YCENT(I),ZCENT(I),1,JGRP,
     &        XBTEMP,YBTEMP,ZBTEMP,XTEMP,YTEMP,ZTEMP,.FALSE.,QMIM)
C
      ELSE
C
C NCS symmetry mode (includes pure intra-molecular mode):
C compute direct distance between groups i, j.
         DO J=1,JGRP
            XTEMP(J)=XCENT(I)-XBTEMP(J)
            YTEMP(J)=YCENT(I)-YBTEMP(J)
            ZTEMP(J)=ZCENT(I)-ZBTEMP(J)
         END DO
      END IF
C
      IF (LGROUP) THEN
C
C in group-by-group-cutoff mode compute distances between group centers
      DO J=1,JGRP
      XJF(J)=XTEMP(J)**2+YTEMP(J)**2+ZTEMP(J)**2
      END DO
      ELSE
C
C in atom-by-atom-cutoff mode
C compute minimum distances between solid rectangles
      DO J=1,JGRP
      XJF(J)= MAX(ZERO,ABS(XTEMP(J))-XSIZ(I)-XSIZ(BOXTEMP(J)))**2
     &     +  MAX(ZERO,ABS(YTEMP(J))-YSIZ(I)-YSIZ(BOXTEMP(J)))**2
     &     +  MAX(ZERO,ABS(ZTEMP(J))-ZSIZ(I)-ZSIZ(BOXTEMP(J)))**2
      END DO
      END IF
C
C check cutoff criterion and interaction criterion between groups i,j
      DO J=1,JGRP
      QSTEST(J)=XJF(J).LT.CUTNB2 .AND.
     &         ABS(INTERG(I)+INTERG(BOXTEMP(J))) .LT.
     &         INTERC(I)+INTERC(BOXTEMP(J))
      END DO
C
C set QSTEST(I) always to true (this is needed below)
      QSTEST(1)=.TRUE.
C
C now that we know all groups j which are interacting with group i
C we generate the index list SPOINT of all atoms of groups j that
C are (possibly) interacting with group i.
      NBSYM=0
      DO J=1,JGRP
      IF (QSTEST(J)) THEN
      DO JJ=IGPBS(BOXTEMP(J))+1,IGPBS(BOXTEMP(J)+1)
      NBSYM=NBSYM+1
      SPOINT(NBSYM)=JJ
      END DO
      END IF
      END DO
C
C not using boxes then check all groups j_{i} to j_{ngrp} - vectorises better
      ELSE
C
      IF (QSYMM) THEN
C
C crystallographic symmetry mode:
C compute minimum image distance between groups i,j
      CALL XPIMAG(XCENT(I),YCENT(I),ZCENT(I),I,NGRP,
     &            XCENTF,YCENTF,ZCENTF,XTEMP,YTEMP,ZTEMP,.FALSE.,QMIM)
C
      ELSE
C
C NCS symmetry mode (includes pure intra-molecular mode):
C compute direct distance between groups i, j.
      DO J=I,NGRP
      XTEMP(J)=XCENT(I)-XCENTF(J)
      YTEMP(J)=YCENT(I)-YCENTF(J)
      ZTEMP(J)=ZCENT(I)-ZCENTF(J)
      END DO
      END IF
C
      IF (LGROUP) THEN
C
C in group-by-group-cutoff mode compute distances between group centers
      DO J=I,NGRP
      XJF(J)=XTEMP(J)**2+YTEMP(J)**2+ZTEMP(J)**2
      END DO
      ELSE
C
C in atom-by-atom-cutoff mode
C compute minimum distances between solid rectangles
      DO J=I,NGRP
      XJF(J)= MAX(ZERO,ABS(XTEMP(J))-XSIZ(I)-XSIZ(J))**2
     &     +  MAX(ZERO,ABS(YTEMP(J))-YSIZ(I)-YSIZ(J))**2
     &     +  MAX(ZERO,ABS(ZTEMP(J))-ZSIZ(I)-ZSIZ(J))**2
      END DO
      END IF
C
C check cutoff criterion and interaction criterion between groups i,j
      DO J=I,NGRP
      QSTEST(J)=XJF(J).LT.CUTNB2 .AND.
     &         ABS(INTERG(I)+INTERG(J)) .LT. INTERC(I)+INTERC(J)
      END DO
C
C set QSTEST(I) always to true (this is needed below)
      QSTEST(I)=.TRUE.
C
C now that we know all groups j which are interacting with group i
C we generate the index list SPOINT of all atoms of groups j that
C are (possibly) interacting with group i.
      NBSYM=0
      DO J=I,NGRP
      IF (QSTEST(J)) THEN
      DO JJ=IGPBS(J)+1,IGPBS(J+1)
      NBSYM=NBSYM+1
      SPOINT(NBSYM)=JJ
      END DO
      END IF
      END DO
C
      ENDIF
C
C == The rest of the code is independent of the initial search method ==
C to make the necessary atom-atom search vectorizable we
C gather the "j" interacting atom coordinates and
C interaction array
      DO J=1,NBSYM
      XJF(J)=XF(SPOINT(J))
      YJF(J)=YF(SPOINT(J))
      ZJF(J)=ZF(SPOINT(J))
      INTERJ(J)=INTERE(SPOINT(J))
      END DO
C
C now we have to search the individual atom i to atom j interactions
C
C loop over atoms of group i
      III=0
      DO IAT=IGPBS(I)+1,IGPBS(I+1)
C
C keep track of the number of atoms in this group
C we use this to exclude self-interactions within the
C same group (NOTE: QSTEST(I) is always true)
      III=III+1
C
C If ISYM is 1, the symmetry operation is the identity.  Therefore
C interactions between atoms i, j have to be check whether they are
C subject to the exclusions.
C We update the exclusion list pointers for atom iat.
CCC      IF (ISYM.EQ.1) THEN  ! ATB 2/18/93
      IF (QINTRA.AND.ISYM.EQ.1) THEN
      IF (IAT.GT.1) THEN
      NXI=IBLOEX(IAT-1)+1
      ELSE
      NXI=1
      END IF
      NXIMAX=IBLOEX(IAT)
      END IF
C
      IF (QSYMM) THEN
C
C crystallographic-symmetry mode:
C determine minimum image distance vectors to any "j" atom
      CHCK=ISYM.EQ.1
      CALL XPIMAG(X(IAT),Y(IAT),Z(IAT),III,NBSYM,
     &            XJF,YJF,ZJF,XTEMP,YTEMP,ZTEMP,CHCK,QMIM)
      ELSE
C
C NCS-symmetry mode (includes pure intramolecular mode):
C direct distance vectors to any "j" atom
      DO J=III,NBSYM
      XTEMP(J)=X(IAT)-XJF(J)
      YTEMP(J)=Y(IAT)-YJF(J)
      ZTEMP(J)=Z(IAT)-ZJF(J)
      END DO
C
      END IF
C
C compute distance vector lenghts
      DO J=III,NBSYM
      XTEMP(J)=XTEMP(J)**2+YTEMP(J)**2+ZTEMP(J)**2
      END DO
C
      IF (LGROUP) THEN
C
C group-by-group-cutoff mode:
C   just check interaction criterion
      DO J=III,NBSYM
      QSTEST(J)=ABS(INTERE(IAT)+INTERJ(J)).LE.+1
      END DO
      ELSE
C
C atom-by-atom-cutoff mode:
C    check cutoff and interaction criterion
      DO J=III,NBSYM
      QSTEST(J)=XTEMP(J).LT.CUTNB2
     &           .AND. ABS(INTERE(IAT)+INTERJ(J)).LE.+1
      END DO
      END IF
C
C for the first symmetry operator we have to separate
C intra-molecular and inter-molecular interactions
      IF (ISYM.EQ.1) THEN
C
C crystallographic-symmetry mode:
      IF (QSYMM) THEN
      DO J=III,NBSYM
      QTEST(J)=QSTEST(J).AND.QMIM(J).EQ.0
      QSTEST(J)=QSTEST(J).AND..NOT.QTEST(J)
      END DO
      ELSE
C
C NCS-symmetry mode (includes pure intramolecular mode):
      DO J=III,NBSYM
      QTEST(J)=QSTEST(J)
      QSTEST(J)=.FALSE.
      END DO
      END IF
C
C if required, compute intra-molecular interactions
      IF (QINTRA) THEN
C
C set number of interactions with atom IAT to zero
      NNIAT=0
C
C exclude self-interactions  (III+1!!)
      DO J=III+1,NBSYM
      IF (QTEST(J)) THEN
C
C Check for exclusions
35    IF (NXI.GT.NXIMAX) GOTO 40
      IF (SPOINT(J).EQ.INBEX(NXI)) GOTO 55
      IF (SPOINT(J).LT.INBEX(NXI)) GOTO 40
      NXI=NXI+1
      GOTO 35
C
40    CONTINUE
C found an atom-atom interaction.
      NNNB=NNNB+1
      NNIAT=NNIAT+1
C
C store this interaction in the temporary array
      ILOCAL(NNIAT)=SPOINT(J)
C
C give warning message if appropriate
      IF (XTEMP(J).LT.WMIN2) THEN
      NVIOL=NVIOL+1
      WRITE(PUNIT,'(17A,F5.2,A)')
     &  ' %atoms "', SEGID(IAT),'-',RESID(IAT),'-',RES(IAT),
     &  '-',TYPE(IAT),
     & '" and "',SEGID(SPOINT(J)),'-',RESID(SPOINT(J)),'-',
     & RES(SPOINT(J)),'-',TYPE(SPOINT(J)),
     & '" only ',SQRT(XTEMP(J)),' A apart'
      END IF
C     Skip interaction label
55    CONTINUE
      END IF
      END DO
C
C copy number of interactions with atom IAT into JNB
      JNB(IAT)=NNIAT
C copy list of interactions with atom IAT into LIST
      IF (NNIAT.GT.0) THEN
      CALL PUTNB(1,NNIAT,ILOCAL,LIST)
      END IF
      END IF
      END IF
C
C if required, compute symmetry related interactions
      IF (QSYMM.OR.QNCS) THEN
C
C test for special crystallographic positions
      IF (NBQSPC.GT.ZERO) THEN
      IF (ISYM.NE.1) THEN
      IF (ABS(XTEMP(III)).LT.NBQSPC) THEN
      QSTEST(III)=.FALSE.
      END IF
      END IF
      END IF
C
C set number of interactions with atom IAT to zero
      NNIAT=0
C
C include self-interactions (III !!)
      DO J=III,NBSYM
      IF (QSTEST(J)) THEN
C
C found a symmetry-related atom-atom interaction
      NNSB=NNSB+1
      NNIAT=NNIAT+1
C
C store this interaction in the temporary array
      ILOCAL(NNIAT)=SPOINT(J)
C
C give warning message if appropriate
      IF (XTEMP(J).LT.WMIN2) THEN
      NVIOL=NVIOL+1
C
      IF (QSYMM) WRITE(PUNIT,'(17A,I3,A,F5.2,A)')
     &  ' %atoms "', SEGID(IAT),'-',RESID(IAT),
     &  '-',RES(IAT),'-',TYPE(IAT),
     & '" and "',SEGID(SPOINT(J)),'-',RESID(SPOINT(J)),
     & '-',RES(SPOINT(J)),'-',TYPE(SPOINT(J)),
     & '"(XSYM#',ISYM,') only ',SQRT(XTEMP(J)),' A apart'
C
      IF (QNCS) WRITE(PUNIT,'(17A,I3,A,F5.2,A)')
     &  ' %atoms "', SEGID(IAT),'-',RESID(IAT),
     &  '-',RES(IAT),'-',TYPE(IAT),
     & '" and "',SEGID(SPOINT(J)),'-',RESID(SPOINT(J)),
     & '-',RES(SPOINT(J)),'-',TYPE(SPOINT(J)),
     & '"(NCS#',ISYM,') only ',SQRT(XTEMP(J)),' A apart'
      END IF
      END IF
      END DO
C
C copy number of interactions with atom IAT into JSB
      JSB(IAT,ISYM)=NNIAT
C copy list of interactions with atom IAT into SLIST
      IF (NNIAT.GT.0) THEN
      CALL PUTNB(1,NNIAT,ILOCAL,SLIST)
      END IF
      END IF
      END DO
C
      ELSE
C
C this group has no interactions whatsoever, but we still have
C to fill the JNB, JSB arrays for the atoms of this group
      IF (QSYMM.OR.QNCS) THEN
      DO IAT=IGPBS(I)+1,IGPBS(I+1)
      JSB(IAT,ISYM)=0
      END DO
      END IF
      IF (QINTRA.AND.ISYM.EQ.1) THEN
      DO IAT=IGPBS(I)+1,IGPBS(I+1)
      JNB(IAT)=0
      END DO
      END IF
      END IF
C
      END DO
      END DO
C
      IF (QINTRA) THEN
      WRITE(PUNIT,'(A,I8,A,I8,A)')
     & ' NBONDS: found ',NNNB,' intra-atom interactions'
      END IF
      IF (QSYMM) THEN
      WRITE(PUNIT,'(A,I8,A,I8,A)')
     & ' NBONDS: found ',NNSB,' symmetry-atom interactions'
      END IF
      IF (QNCS) THEN
      WRITE(PUNIT,'(A,I8,A,I8,A)')
     & ' NBONDS: found ',NNSB,' NCS-atom interactions'
      END IF
C
C declare number of violations
      IF (QINTRA.OR.QSYMM.OR.QNCS) THEN
      IF (NVIOL.GT.0) THEN
      WRITE(PUNIT,'(A,I8,A,I8,A)')
     & ' NBONDS: found ',NVIOL,' nonbonded violations'
      END IF
      DBPREC=NVIOL
      CALL DECLAR( 'VIOLATIONS', 'DP', ' ', DBCOMP, DBPREC )
      END IF
C
      RETURN
      END
C
C=================================================================
C
      SUBROUTINE GRPRAD(NGRP,XSIZ,YSIZ,ZSIZ,GRDMMX,INTERG)
C
C Paul Adams 22/4/93
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      DOUBLE PRECISION GRDMMX
      INTEGER NGRP, INTERG(*)
      DOUBLE PRECISION XSIZ(*), YSIZ(*), ZSIZ(*)
C local
      INTEGER I
      DOUBLE PRECISION MAXX, MAXY, MAXZ
C
      MAXX=RSMALL
      MAXY=RSMALL
      MAXZ=RSMALL
C
C get maximum sizes
      DO I=1,NGRP
      IF (ABS(INTERG(I)).LT.+9999) THEN
         IF(XSIZ(I).GT.MAXX) MAXX=XSIZ(I)
         IF(YSIZ(I).GT.MAXY) MAXY=YSIZ(I)
         IF(ZSIZ(I).GT.MAXZ) MAXZ=ZSIZ(I)
      END IF
      END DO
C
      GRDMMX=MAX(MAXX,MAXY,MAXZ)
C
C If point groups then set grdmmx to something sensible
      IF (GRDMMX.LT.0.5D0) THEN
         GRDMMX=0.5D0
      END IF
C
      RETURN
      END
C=================================================================
C
      SUBROUTINE BOXES(BOXCUT,NGRP,XCENT,YCENT,ZCENT,
     &                 NBOXES,XBOX,YBOX,ZBOX,BOX,BOXSIZ,BOXMEM,QBOXL,
     &                 INTERG)
C
C to calculate the boxes for the grid based part of xpsea2
C nboxes = number of boxes created
C xbox,ybox,zbox = coordinates of centres of boxes
C box = list giving box number for groups {1....ngrp}
C boxsiz = list giving the pointer to the first group in a box
C boxmem = continuous list of members of box i_{1}...i_{nboxes} referenced
C          by pointer list boxsiz
C
C Paul Adams 16/4/93
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'fastnb.inc'
C
      DOUBLE PRECISION BOXCUT
      INTEGER NGRP
      DOUBLE PRECISION XCENT(*), YCENT(*), ZCENT(*)
      INTEGER NBOXES
      DOUBLE PRECISION XBOX(*), YBOX(*), ZBOX(*)
      INTEGER BOX(*), BOXSIZ(*), BOXMEM(*)
      INTEGER INTERG(*)
      LOGICAL QBOXL
C local
      DOUBLE PRECISION MINX, MINY, MINZ, MAXX, MAXY, MAXZ
      DOUBLE PRECISION ROOT2, ROOT3
      DOUBLE PRECISION CGX, CGY, CGZ, RATIO
      DOUBLE PRECISION STEP, XSTART, YSTART, ZSTART
      DOUBLE PRECISION VECX, VECY, VECZ
      INTEGER I, J, K, XPOS, YPOS, ZPOS, COUNTER
      INTEGER SIZEX, SIZEY, SIZEZ, BOXTOT
      LOGICAL BNDFLG, CMPLTD, RATDON, BOXLIM
C parameters
      DOUBLE PRECISION LARGE
      PARAMETER (LARGE=9999.0D0)
      INTEGER ZERO
      PARAMETER (ZERO=0)
C
      ROOT2=SQRT(2.0D0)
      ROOT3=SQRT(3.0D0)
C
      CMPLTD=.FALSE.
      BNDFLG=.FALSE.
      RATDON=.FALSE.
      BOXLIM=.FALSE.
C
      MINX=R4BIG
      MINY=R4BIG
      MINZ=R4BIG
      MAXX=RSMALL
      MAXY=RSMALL
      MAXZ=RSMALL
C
C first get limits of group centre coordinates
      DO I=1,NGRP
         IF (ABS(INTERG(I)).LT.+9999) THEN
            IF(XCENT(I).LT.MINX) MINX=XCENT(I)
            IF(YCENT(I).LT.MINY) MINY=YCENT(I)
            IF(ZCENT(I).LT.MINZ) MINZ=ZCENT(I)
            IF(XCENT(I).GT.MAXX) MAXX=XCENT(I)
            IF(YCENT(I).GT.MAXY) MAXY=YCENT(I)
            IF(ZCENT(I).GT.MAXZ) MAXZ=ZCENT(I)
         END IF
      END DO
C
      CGX=(MAXX+MINX)/2.0
      CGY=(MAXY+MINY)/2.0
      CGZ=(MAXZ+MINZ)/2.0
      VECX=(MAXX-MINX)
      VECY=(MAXY-MINY)
      VECZ=(MAXZ-MINZ)
C
C repeat this until the min and max number of boxes is satisfied
      DO WHILE (.NOT.BOXLIM)
C
C step size in orthogonal cartesian space - calculated from length of cube
C diagonal
         STEP=2.0*(BOXCUT/ROOT3)
         SIZEX=INT(VECX/STEP)+1
         SIZEY=INT(VECY/STEP)+1
         SIZEZ=INT(VECZ/STEP)+1
         BOXTOT=SIZEX*SIZEY*SIZEZ
         IF ((BOXTOT.GT.MNBOX).AND.(BOXTOT.LE.MXBOX)) THEN
            BOXLIM=.TRUE.
C
C loop until the minimum group/box ratio is satisfied
            DO WHILE (.NOT.RATDON)
               STEP=2.0*(BOXCUT/ROOT3)
C
C loop until all groups lie within boxes
               DO WHILE (.NOT.CMPLTD)
                  SIZEX=INT(VECX/STEP)+1
                  SIZEY=INT(VECY/STEP)+1
                  SIZEZ=INT(VECZ/STEP)+1
                  VECX=SIZEX*STEP
                  VECY=SIZEY*STEP
                  VECZ=SIZEZ*STEP
C
                  XSTART=CGX-(VECX/2.0)-(STEP/2.0)
                  YSTART=CGY-(VECY/2.0)-(STEP/2.0)
                  ZSTART=CGZ-(VECZ/2.0)-(STEP/2.0)
C
C first we generate a box of cubes around the molecule
                  NBOXES=0
                  DO I=1,SIZEZ
                     DO J=1,SIZEY
                        DO K=1,SIZEX
                           NBOXES=NBOXES+1
                           XBOX(NBOXES)=XSTART+(K*STEP)
                           YBOX(NBOXES)=YSTART+(J*STEP)
                           ZBOX(NBOXES)=ZSTART+(I*STEP)
                           BOXSIZ(NBOXES)=0
                        END DO
                     END DO
                  END DO
C
C this determines which box a group belongs to
                  CMPLTD=.TRUE.
                  BNDFLG=.FALSE.
                  DO I=1,NGRP
                  IF (ABS(INTERG(I)).LT.+9999) THEN
                  XPOS=1+INT((XCENT(I)-(XSTART+(STEP/2.0)))/STEP)
                  CALL BNDCHK(XPOS,1,SIZEX,BNDFLG)
                  YPOS=1+INT((YCENT(I)-(YSTART+(STEP/2.0)))/STEP)
                  CALL BNDCHK(YPOS,1,SIZEY,BNDFLG)
                  ZPOS=1+INT((ZCENT(I)-(ZSTART+(STEP/2.0)))/STEP)
                  CALL BNDCHK(ZPOS,1,SIZEZ,BNDFLG)
                  BOX(I)=
     &              1+(XPOS-1)+(YPOS-1)*SIZEX+(ZPOS-1)*(SIZEX*SIZEY)
                  BOXSIZ(BOX(I))=BOXSIZ(BOX(I))+1
                  ELSE
                  BOX(I)=ZERO
                  END IF
                  END DO
                  IF(BNDFLG) THEN
                     IF(TIMER.GT.1) THEN
                     WRITE(6,'(2A)')'BOX-BND-SIZ% ',
     &           'There are groups outside the macro-box - RESIZING'
                     CMPLTD=.FALSE.
                     END IF
                  END IF
               END DO
C
C this should have got us the filled box array
C now we have to take any empty boxes out of the lists
               COUNTER=0
               DO I=1,NBOXES
                  IF(BOXSIZ(I).NE.0) THEN
                     COUNTER=COUNTER+1
                     XBOX(COUNTER)=XBOX(I)
                     YBOX(COUNTER)=YBOX(I)
                     ZBOX(COUNTER)=ZBOX(I)
                     DO J=1,NGRP
                        IF(BOX(J).EQ.I) BOX(J)=COUNTER
                     END DO
                  END IF
               END DO
               NBOXES=COUNTER
C
               RATIO=NGRP/NBOXES
               IF (RATIO.LT.RATMIN) THEN
                  RATDON=.FALSE.
                  CMPLTD=.FALSE.
                  BOXCUT=BOXCUT*ROOT2
               ELSE
                  RATDON=.TRUE.
               ENDIF
            END DO
         ELSE
            IF (BOXTOT.LT.MNBOX) THEN
C
C If there are too few boxes we turn the box method off
               QBOXL=.FALSE.
               BOXLIM=.TRUE.
            ELSE
               BOXCUT=BOXCUT*ROOT2
               BOXLIM=.FALSE.
            ENDIF
         ENDIF
      ENDDO
C
C now fill up the box member list and number of members per box
      IF (QBOXL) THEN
         BOXSIZ(1)=1
         COUNTER=0
         DO I=1,NBOXES
            DO J=1,NGRP
               IF(BOX(J).EQ.I) THEN
                  COUNTER=COUNTER+1
                  BOXMEM(COUNTER)=J
               END IF
            END DO
            BOXSIZ(I+1)=COUNTER+1
         END DO
      ENDIF
C
C only print this if the box method is actually being used
C and timer is set appropriately
      IF(TIMER.GT.1) THEN
      IF(QBOXL.AND.(NBOXES.GT.1)) THEN
      WRITE(6,'(A,3I3)')' BOX: Macro box size ',SIZEX,SIZEY,SIZEZ
      WRITE(6,'(A,I3)')' BOX: Number of boxes filled ',NBOXES
      ENDIF
      ENDIF
C
      RETURN
      END
C=================================================================
      SUBROUTINE BNDCHK(X,LOWER,UPPER,BNDFLG)
C
C routine to check that x is within the limits specified by
C lower and upper - arguments are assumed integers!!
C
C Paul Adams 16/4/93
C
      INTEGER X, LOWER, UPPER
      LOGICAL BNDFLG
C
      IF (X.LT.LOWER) THEN
         BNDFLG=.TRUE.
      ELSE IF (X.GT.UPPER) THEN
         BNDFLG=.TRUE.
      END IF
C
      RETURN
      END
C=================================================================
