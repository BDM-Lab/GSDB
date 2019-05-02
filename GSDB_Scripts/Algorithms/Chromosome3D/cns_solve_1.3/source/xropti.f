      SUBROUTINE XROPTI
C
C Restrained, isotropic B-factor optimizer.
C
C All free (as opposed to "fixed") atoms with SCATter entries
C are refined.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'xtarget.inc'
C pointer
      INTEGER FLAGS
      INTEGER IBGROUP, IAGROUP
      INTEGER INDEX
C begin
      FLAGS=ALLHP(INTEG4(NATOM))
      IBSIGM=ALLHP(IREAL8(NATOM))
      IASIGM=ALLHP(IREAL8(NATOM))
      IBGROUP=ALLHP(INTEG4(NATOM))
      IAGROUP=ALLHP(INTEG4(NATOM))
      INDEX=ALLHP(INTEG4(NATOM))
      ILIST=ALLHP(INTEG4(NATOM))
      IPOINT=ALLHP(INTEG4(NATOM+1))
      IOMAIN=ALLHP(IREAL8(NATOM))
C
      CALL XROPT2(HEAP(FLAGS),HEAP(IBSIGM),HEAP(IASIGM),
     &            HEAP(IBGROUP),HEAP(IAGROUP),
     &            HEAP(INDEX),HEAP(ILIST),HEAP(IPOINT),
     &            HEAP(HPATOF),HEAP(IOMAIN))
C  
      CALL FREHP(IOMAIN,IREAL8(NATOM))   
      CALL FREHP(IPOINT,INTEG4(NATOM+1))
      CALL FREHP(ILIST,INTEG4(NATOM))
      CALL FREHP(INDEX,INTEG4(NATOM))
      CALL FREHP(IAGROUP,INTEG4(NATOM))
      CALL FREHP(IBGROUP,INTEG4(NATOM))
      CALL FREHP(IASIGM,IREAL8(NATOM))
      CALL FREHP(IBSIGM,IREAL8(NATOM))
      CALL FREHP(FLAGS,INTEG4(NATOM))
      RETURN
      END
C=====================================================================
      SUBROUTINE XROPT2(FLAGS,BSIGMA,ASIGMA,BGROUP,AGROUP,
     &                  INDEX,LIST,POINT,XRATOF,OMAIN)
C
C See XROPTI above
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'timer.inc'
      INTEGER FLAGS(*)
      DOUBLE PRECISION BSIGMA(*), ASIGMA(*)
      INTEGER BGROUP(*), AGROUP(*)
      INTEGER INDEX(*), LIST(*), POINT(*), XRATOF(*)
      DOUBLE PRECISION OMAIN(*)
C local
      INTEGER CGSTEP, IAT, NFLAGS, NCGRP
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION STEP, TOLER, TEMP
      LOGICAL DEBUG
      INTEGER NBGRP, NAGRP
      INTEGER NLIST, I, MAXFEV
      CHARACTER*6 MINTYPE
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO, R15
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, R15=1.5D0)
C begin
C
C default values
      CGSTEP=0
C
C maximum number of function evaluation for lbgfs minimizer
      MAXFEV=10
C Initial step size for conjugate gradient minimizer
      STEP=10.0D0
C Tolerance on gradient
      TOLER=0.000001D0
C B-factor restraint weight
      RWGHT=-ONE
      QRESTR=.TRUE.
C bond restraint sigma
      DO IAT=1,NATOM
      BSIGMA(IAT)=R15
      END DO
C angle restraint sigma
      DO IAT=1,NATOM
      ASIGMA(IAT)=TWO
      END DO
C minimum allowed B-factor
      BFMIN=2.0D0
C maximum allowed B-factor
      BFMAX=100.D0
C
C number of groups for bond restraints
      NBGRP=0
C atomic selections of groups for bond restraints
      DO IAT=1,NATOM
      BGROUP(IAT)=0
      END DO
C number of groups for angle restraints
      NAGRP=0
C atomic selections of groups for angle restraints
      DO IAT=1,NATOM
      AGROUP(IAT)=0
      END DO
C initialize group list
      NCGRP=0
      NLIST=0
      POINT(1)=0
C
C debugging flag
      DEBUG=.FALSE.
C
C minimization method
      MINTYPE='POWELL'
C
C parsing
      CALL PUSEND('XROPTI>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('XROPTI>')
      CALL MISCOM('XROPTI>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-optimize-bfactor')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NSTE') THEN
      CALL NEXTI('Number of steps= ',CGSTEP)
C======================================================================
      ELSE IF (WD(1:4).EQ.'STEP'.OR.WD(1:4).EQ.'DROP') THEN
      CALL NEXTF('DROP=',STEP)
C======================================================================
      ELSE IF (WD(1:4).EQ.'TOLE') THEN
      CALL NEXTF('TOLErance=',TOLER)
C======================================================================
      ELSE IF (WD(1:4).EQ.'MAXF') THEN
      CALL NEXTI('MAXFev=',MAXFEV)
C======================================================================
      ELSE IF (WD(1:4).EQ.'METH') THEN
      CALL NEXTST('METHod=',MINTYPE)
C======================================================================
      ELSE IF (WD(1:4).EQ.'BOND') THEN
      TEMP=ZERO
      CALL NEXTF('Bond B-factor Sigma=',TEMP)
      DO IAT=1,NATOM
      BSIGMA(IAT)=TEMP
      END DO
C======================================================================
      ELSE IF (WD(1:4).EQ.'BSIG') THEN
      CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
      CALL MAKIND(FLAGS,NATOM,NFLAGS)
      CALL NEXTF('BSIGma=',TEMP)
      DO IAT=1,NFLAGS
      BSIGMA(FLAGS(IAT))=TEMP
      END DO
      IF (NFLAGS.GT.0) THEN
      NBGRP=NBGRP+1
      DO IAT=1,NFLAGS
      BGROUP(FLAGS(IAT))=NBGRP
      END DO
      ELSE
      WRITE(6,'(A)') 'XROPT-WRN: zero atom is selected!'
      END IF
C======================================================================
      ELSE IF (WD(1:4).EQ.'ANGL') THEN
      TEMP=ZERO
      CALL NEXTF('Angle B-factor Sigma=',TEMP)
      DO IAT=1,NATOM
      ASIGMA(IAT)=TEMP
      END DO
C======================================================================
      ELSE IF (WD(1:4).EQ.'ASIG') THEN
      CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
      CALL MAKIND(FLAGS,NATOM,NFLAGS)
      CALL NEXTF('ASIGma=',TEMP)
      DO IAT=1,NFLAGS
      ASIGMA(FLAGS(IAT))=TEMP
      END DO
      IF (NFLAGS.GT.0) THEN
      NAGRP=NAGRP+1
      DO IAT=1,NFLAGS
      AGROUP(FLAGS(IAT))=NAGRP
      END DO
      ELSE
      WRITE(6,'(A)') 'XROPT-WRN: zero atom is selected!'
      END IF
C======================================================================
      ELSE IF (WD(1:4).EQ.'BFMI'.OR.WD(1:4).EQ.'BMIN') THEN
      CALL NEXTF('BMIN=',BFMIN)
C======================================================================
      ELSE IF (WD(1:4).EQ.'BFMA'.OR.WD(1:4).EQ.'BMAX') THEN
      CALL NEXTF('BMAX=',BFMAX)
C======================================================================
      ELSE IF (WD(1:4).EQ.'RWEI') THEN
      CALL NEXTF('B-factor Restraint Weight=',RWGHT)
      QRESTR=ABS(RWGHT).GT.RSMALL
C
C begin modification, 12/23/08, ATB
C======================================================================
      ELSE IF (WD(1:4).EQ.'GROU') THEN
      CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
      NFLAGS=0
      DO IAT=1,XRNATF
      IF (FLAGS(XRATOF(IAT)).EQ.1.AND.IMOVE(XRATOF(IAT)).EQ.0) THEN
      NFLAGS=NFLAGS+1
      INDEX(NFLAGS)=XRATOF(IAT)
      END IF
      END DO
C
      IF (NFLAGS.GT.0) THEN
C make a new group
      IF (NCGRP.GE.NATOM) THEN
      WRITE(6,'(A)') ' %XROPTI-ERR: too many groups specified'
      GOTO 9999
      ELSE
      NCGRP=NCGRP+1
      DO IAT=1,NFLAGS
      IF (NLIST.GE.NATOM) THEN
      WRITE(6,'(A)')
     &' %XROPTI-ERR: too many atoms selected--check for group overlaps'
      GOTO 9999
      ELSE
      NLIST=NLIST+1
      LIST(NLIST)=INDEX(IAT)
      END IF
      END DO
      POINT(NCGRP+1)=NLIST
      END IF
      END IF
C======================================================================
C end modification
C
      ELSE IF (WD(1:4).EQ.'DEBU') THEN
      CALL NEXTLO('DEBUg=',DEBUG)
C======================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' ---------------------------',
     & 'B-factor-refinement-parameters----------------------'
      WRITE(6,'(A,I6,A,F12.5,A,F12.5,A,F12.5,A,F12.5)')
     &' | NSTEp=',CGSTEP,' DROP=',STEP,' TOLErance=',TOLER,
     &' BMIN=',BFMIN, ' BMAX=',BFMAX
      WRITE(6,'(A,I5)') ' | MAXFevaluations=',MAXFEV
      WRITE(6,'(2A)') ' | METHod=',MINTYPE
      IF (QRESTR) THEN
      WRITE(6,'(A)')
     &' | restrained B-factor refinement'
C
      IF (RWGHT.LT.ZERO) THEN
      WRITE(6,'(2A)') ' | Weights for B-factor restraints will be set',
     & ' first time by comparing gradients'
      ELSE
      WRITE(6,'(A,F12.5)')
     & ' | Weight for B-factor restraints: ',RWGHT
      END IF
      ELSE
      WRITE(6,'(A)') ' | B-factor restraints turned off'
      END IF
      IF (NBGRP.GT.0) WRITE(6,'(A,I5,A)')
     & ' | There are ',NBGRP,' groups of bond-restraints'
      IF (NAGRP.GT.0) WRITE(6,'(A,I5,A)')
     & ' | There are ',NAGRP,' groups of angle-restraints'
      WRITE(6,'(2A)') ' ---------------------------',
     & '----------------------------------------------------'
C=====================================================================
      ELSE
      CALL CHKEND('XROPTI>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
C - begin modification, 12/23/08, ATB
C check for overlapping group definitions
      ERROR=.FALSE.
      DO IAT=1,NATOM
      FLAGS(IAT)=0
      END DO
      DO I=1,NCGRP
      DO IAT=POINT(I)+1,POINT(I+1)
      FLAGS(LIST(IAT))=FLAGS(LIST(IAT))+1
      END DO
      END DO
      DO IAT=1,NATOM
      IF (FLAGS(IAT).GT.1) THEN
      WRITE(6,'(9A)')
     & ' %XGROUP-ERR: B-GROUP definition overlap for ( ',
     & SEGID(IAT),' ',RES(IAT),' ',RESID(IAT),' ',TYPE(IAT),' )'
      ERROR=.TRUE.
      END IF
      END DO
C
C now add individual atoms to the group list (as groups with one atom)
      I=POINT(NCGRP+1)
      DO IAT=1,XRNATF
      IF (IMOVE(XRATOF(IAT)).EQ.0) THEN
      IF (FLAGS(XRATOF(IAT)).EQ.0) THEN
      NCGRP=NCGRP+1
      I=I+1
      POINT(NCGRP+1)=I
      LIST(I)=XRATOF(IAT)
      END IF
      END IF
      END DO
C
C debug
      IF (WRNLEV.GE.15) THEN
      DO I=1,NCGRP
      DO IAT=POINT(I)+1,POINT(I+1)
      WRITE(6,'(A,I6,9A)')
     & ' group ',I,' ( ',
     & SEGID(LIST(IAT)),' ',
     & RES(LIST(IAT)),' ',RESID(LIST(IAT)),
     & ' ',TYPE(LIST(IAT)),' )'      
      END DO
      END DO
      END IF 
C end debug 
C
      IF (.NOT.ERROR) THEN
C end modification
C
C
C carry out minimization or just do debugging of derivatives
      FWMAIN=ALLHP(IREAL8(NATOM))
      CALL XBREF(CGSTEP,STEP,TOLER,DEBUG,BGROUP,AGROUP,NBGRP,NAGRP,
     &           HEAP(FWMAIN),HEAP(HPATOF),OMAIN,POINT,NCGRP,MINTYPE,
     &           MAXFEV)
      CALL FREHP(FWMAIN,IREAL8(NATOM))
C
C put the restraint weight into a variable
      CALL DECLAR('B_RWEIGHT','DP',' ',DBCOMP,RWGHT)
      END IF
C
9999  CONTINUE
C
      RETURN
      END
C======================================================================
      SUBROUTINE XBREF(CGSTEP,STEP,TOLER,DEBUG,
     &                    BGROUP,AGROUP,NBGRP,NAGRP,MWMAIN,XRATOF,OMAIN,
     &                    POINT,NCGRP,MINTYPE,MAXFEV)
C
C routine minimizes target function or checks consistency of derivatives
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'numbers.inc'
      INTEGER CGSTEP
      DOUBLE PRECISION STEP
      DOUBLE PRECISION TOLER
      LOGICAL DEBUG
      INTEGER BGROUP(*), AGROUP(*)
      INTEGER NBGRP, NAGRP
      DOUBLE PRECISION MWMAIN(*)
      INTEGER XRATOF(*)
      DOUBLE PRECISION OMAIN(*)
      INTEGER POINT(*), NCGRP
      CHARACTER*(*) MINTYPE
      INTEGER MAXFEV
C local
      DOUBLE PRECISION TARGET
      INTEGER I, II, III
      INTEGER IER, IAT
      INTEGER M, NPRINT, IW
C pointer
      INTEGER BRMSD, ARMSD, NBRG, NTRG
      INTEGER PARAM, GRDENT, WORK, DIAG, BEST
C external
      EXTERNAL XBFUNC, XBFUNC2
C
C determine number of B-factors to be refined
      II=0
      III=0
      DO I=1,NCGRP
      IF (POINT(I+1)-POINT(I)-1.GT.1) THEN
      II=II+1
      END IF
      DO IAT=POINT(I)+1,POINT(I+1)
      III=III+1
      END DO
      END DO
      IF (II.EQ.0) THEN
      WRITE (6,'(A,I8,A)')
     & ' XBREF: ',III,' individual B-factors will be refined.'
      ELSE
      WRITE (6,'(A,I8,A,/,A,I8,A)')
     & ' XBREF: ',II,' constrained B-factor groups will be refined.',    
     & ' XBREF: ',NCGRP-II,' individual B-factors will be refined,'
      END IF
C
      IF (NCGRP.GT.0) THEN
C
C allocate space
      PARAM=ALLHP(IREAL8(NCGRP))
      GRDENT=ALLHP(IREAL8(NCGRP))
      IIBR=ALLHP(INTEG4(NBOND))
      IJBR=ALLHP(INTEG4(NBOND))
      IITR=ALLHP(INTEG4(NTHETA))
      IKTR=ALLHP(INTEG4(NTHETA))
      IFLAG=ALLHP(INTEG4(NATOM))
      ITEMP=ALLHP(IREAL8(XRNATF))
      BRMSD=ALLHP(IREAL8(NBGRP+2))
      ARMSD=ALLHP(IREAL8(NAGRP+2))
      NBRG=ALLHP(INTEG4(NBGRP+2))
      NTRG=ALLHP(INTEG4(NAGRP+2))
C
C skip if NSTEP<0 (CGSTEP<0) - JSJ, 10/19/95
      IF (CGSTEP.GE.0) THEN
C
C if restraints are requested, set up the lists
      IF (QRESTR) THEN
      CALL XBSET(HEAP(IIBR),HEAP(IJBR),HEAP(IITR),HEAP(IKTR),
     &     HEAP(IFLAG),HEAP(HPATOF),HEAP(ITEMP),HEAP(IBSIGM),
     &     HEAP(IASIGM))
C
C calculate initial RMS deviations based on bond and angle list
      CALL XBRMS(HEAP(IIBR),HEAP(IJBR),HEAP(IITR),HEAP(IKTR),
     &     HEAP(IFLAG),HEAP(HPATOF),HEAP(ITEMP),
     &     BGROUP,AGROUP,NBGRP,NAGRP,HEAP(BRMSD),HEAP(ARMSD),
     &     HEAP(NBRG),HEAP(NTRG))
      END IF
C
C initialize the group b factor shifts to zero
      CALL FILLR8(HEAP(PARAM),NCGRP,ZERO)
C
C note: the inidividual b factor will be the sum of the b
C factor shifts and the previously assigned b individual b factor
C which is stored in OMAIN
      DO I=1,NATOM
      OMAIN(I)=WMAIN(I)
      END DO
C
      IF (DEBUG) THEN
      CALL XBDEBG(HEAP(PARAM),HEAP(GRDENT),STEP,NCGRP)
      ELSE
C
C do minimization
      ICYCLE=0
C
C Set tolerance to the sum of squares
      TOLER=(TOLER**2)*NCGRP
C
C
C initialize the lowest target value
      MTARG=R4BIG
      MCYCLE=0
C
C create array that will contain the B-factors that produce 
C the lowest target value
      I=0
      DO IAT=1,XRNATF
      IF (IMOVE(XRATOF(IAT)).EQ.0) THEN
      I=I+1
      MWMAIN(XRATOF(IAT))=WMAIN(XRATOF(IAT))
      END IF
      END DO
C
      IF (MINTYPE.EQ.'POWELL') THEN
C
C POWELL-method minimizer
      WORK=ALLHP(IREAL8(6*NCGRP))
      CALL ZXCGR(XBFUNC, NCGRP, TOLER, CGSTEP, STEP, HEAP(PARAM),
     &           HEAP(GRDENT), TARGET, HEAP(WORK), IER)
      CALL FREHP(WORK,IREAL8(6*NCGRP))
C
C print termination conditions
      IF (IER.EQ.0) THEN
      WRITE(6,'(A)')' XBREF: gradient converged'
      ELSE IF (IER.EQ.129) THEN
      WRITE(6,'(A)')' XBREF: Line search terminated'
      ELSE IF (IER.EQ.130) THEN
      WRITE(6,'(A)')' %XBREF-ERR: Search Direction Uphill'
      ELSE IF (IER.EQ.131) THEN
      WRITE(6,'(A)')' XBREF: step limit reached'
      ELSE IF (IER.EQ.132) THEN
      WRITE(6,'(A)')' %XBREF-ERR: Failure to reduce E'
      END IF
C
      ELSE IF (MINTYPE.EQ.'LBFGS') THEN
C
C LBFGS minimizer
      M=6
      NPRINT=1
      IW=NCGRP*(2*M+1)+2*M
      DIAG=ALLHP(IREAL8(NCGRP))
      BEST=ALLHP(IREAL8(NCGRP))
      WORK=ALLHP(IREAL8(IW))
      CALL ZXLBFGS(XBFUNC2,NCGRP,M,HEAP(PARAM),HEAP(GRDENT),TARGET,
     &             .FALSE.,HEAP(DIAG),NPRINT,TOLER,CGSTEP,HEAP(WORK),
     &             HEAP(BEST),MAXFEV)
      CALL FREHP(WORK,IREAL8(IW))
      CALL FREHP(BEST,IREAL8(NCGRP))
      CALL FREHP(DIAG,IREAL8(NCGRP))
      ELSE
      CALL WRNDIE(-5,'XBREF','unknown minimization method.')
      END IF
C
C set parameters to the values that produces the lowest target value
      I=0
      DO IAT=1,XRNATF
      IF (IMOVE(XRATOF(IAT)).EQ.0) THEN
      I=I+1
      WMAIN(XRATOF(IAT))=MWMAIN(XRATOF(IAT))
      END IF
      END DO
      WRITE(6,'(A,I6,A)')
     & ' XBREF: using the parameters (at cycle ', MCYCLE,
     & ') which produce the minimum target. '
C
C
      END IF
C
      END IF
C
C calculate final RMS deviations based on bond and angle list
      IF (QRESTR) THEN
      CALL XBRMS(HEAP(IIBR),HEAP(IJBR),HEAP(IITR),HEAP(IKTR),
     &     HEAP(IFLAG),HEAP(HPATOF),HEAP(ITEMP),
     &     BGROUP,AGROUP,NBGRP,NAGRP,HEAP(BRMSD),HEAP(ARMSD),
     &     HEAP(NBRG),HEAP(NTRG))
      END IF
C
C release dynamic HEAP space
      CALL FREHP(NTRG,INTEG4(NAGRP+2))
      CALL FREHP(NBRG,INTEG4(NBGRP+2))
      CALL FREHP(ARMSD,IREAL8(NAGRP+2))
      CALL FREHP(BRMSD,IREAL8(NBGRP+2))
      CALL FREHP(ITEMP,IREAL8(XRNATF))
      CALL FREHP(IFLAG,INTEG4(NATOM))
      CALL FREHP(IKTR,INTEG4(NTHETA))
      CALL FREHP(IITR,INTEG4(NTHETA))
      CALL FREHP(IJBR,INTEG4(NBOND))
      CALL FREHP(IIBR,INTEG4(NBOND))
      CALL FREHP(GRDENT,IREAL8(NCGRP))
      CALL FREHP(PARAM,IREAL8(NCGRP))
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XBDEBG(PARAM,GRDENT,STEP,NCGRP)
C
C check consistency of B-factor derivatives
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'xropti.inc'
      DOUBLE PRECISION PARAM(*), GRDENT(*), STEP
      INTEGER NCGRP
C local
      INTEGER IAT, J
      DOUBLE PRECISION TARGET, NUMGRD
C parameters
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
C begin
      ICYCLE=-1
C
      DO J=1,NCGRP
      PARAM(J)=PARAM(J)+STEP
      CALL XBFUNC(NCGRP,PARAM,TARGET,GRDENT)
      NUMGRD=TARGET
      PARAM(J)=PARAM(J)-TWO*STEP
      CALL XBFUNC(NCGRP,PARAM,TARGET,GRDENT)
      NUMGRD=NUMGRD-TARGET
      NUMGRD=NUMGRD/(TWO*STEP)
      PARAM(J)=PARAM(J)+STEP
      CALL XBFUNC(NCGRP,PARAM,TARGET,GRDENT)
      WRITE(6,'(A,I6,A,F12.5,A,F12.5,A,F12.5)')
     &' XBDEBUG: index # ',J,' analyt.=',GRDENT(J),
     &' finite=',NUMGRD
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XBFUNC(NCGRP,PARAM,TARGET,GRDENT)
C
C ZXCGR (Powell minimizer) target function routine for restrained, isotropic B-factor
C minimization.  Front-end routine which passes the pointers
C of some dynamic arrays to the next routine.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'mtf.inc'
      INTEGER NCGRP
      DOUBLE PRECISION PARAM(*), TARGET, GRDENT(*)
C begin
      CALL XBFUN2(NCGRP,PARAM,TARGET,GRDENT,HEAP(IIBR),HEAP(IJBR),
     &            HEAP(IITR),HEAP(IKTR),HEAP(ITEMP),HEAP(HPATOF),
     &            HEAP(HPDT),HEAP(IFLAG),HEAP(IBSIGM),HEAP(IASIGM),
     &            HEAP(FWMAIN),HEAP(IOMAIN),
     &            HEAP(ILIST),HEAP(IPOINT))
      RETURN
      END
C======================================================================
      SUBROUTINE XBFUNC2(NCGRP,PARAM,TARGET,GRDENT,ITER,QPRINT)
C
C ZXLBFGS (LBGFGS minimizer) target function routine for restrained, isotropic B-factor
C minimization.  Front-end routine which passes the pointers
C of some dynamic arrays to the next routine.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'mtf.inc'
      INTEGER NCGRP
      DOUBLE PRECISION PARAM(*), TARGET, GRDENT(*)
      INTEGER ITER
      LOGICAL QPRINT
C begin
      CALL XBFUN2(NCGRP,PARAM,TARGET,GRDENT,HEAP(IIBR),HEAP(IJBR),
     &            HEAP(IITR),HEAP(IKTR),HEAP(ITEMP),HEAP(HPATOF),
     &            HEAP(HPDT),HEAP(IFLAG),HEAP(IBSIGM),HEAP(IASIGM),
     &            HEAP(FWMAIN),HEAP(IOMAIN),
     &            HEAP(ILIST),HEAP(IPOINT))
      RETURN
      END
C=====================================================================
      SUBROUTINE XBFUN2(NCGRP,PARAM,TARGET,GRDENT,IBR,JBR,ITR,KTR,TEMP,
     &                  XRATOF,XRDT,FLAG,BSIGMA,ASIGMA,MWMAIN,OMAIN,
     &                  LIST,POINT)
C
C see routine above
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'heap.inc'
      INTEGER NCGRP
      DOUBLE PRECISION PARAM(*), TARGET, GRDENT(*)
      INTEGER IBR(*), JBR(*), ITR(*), KTR(*)
      DOUBLE PRECISION TEMP(*)
      INTEGER XRATOF(*)
      DOUBLE PRECISION XRDT(*)
      INTEGER FLAG(*)
      DOUBLE PRECISION BSIGMA(*), ASIGMA(*), MWMAIN(*)
      DOUBLE PRECISION OMAIN(*)
      INTEGER LIST(*), POINT(*)
C local
      DOUBLE PRECISION BENER, RENER, DB, DNORM1, DNORM2, GRAD
      INTEGER I, II, IAT, K, NNCSX
C parameters
      DOUBLE PRECISION ZERO, TWO, FOUR, EIGHT
      PARAMETER (ZERO=0.0D0, TWO=2.0D0, FOUR=4.0D0, EIGHT=8.0D0)
C begin
C
C Initialize the gradient vector
      DO I=1,NCGRP
      GRDENT(I)=ZERO
      END DO
C
C add parameter shifts to OMAIN, apply B-factor truncations
C store current parameters
      DO I=1,NCGRP
      DO IAT=POINT(I)+1,POINT(I+1)
      WMAIN(LIST(IAT))=OMAIN(LIST(IAT))+PARAM(I)
      WMAIN(LIST(IAT))=MIN(MAX(WMAIN(LIST(IAT)),
     &                         BFMIN),BFMAX)
      END DO
      END DO    
C
C compute the diffraction energy term and gradient
C ================================================
      CALL XCALCS(.FALSE.,.TRUE.,XRE,
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
      BENER=XRE
C
C
C if required then compute B-factor restraints
C ============================================
      RENER=ZERO
      IF (QRESTR) THEN
C
C initialize restraint derivatives
      DO IAT=1,XRNATF
      TEMP(IAT)=ZERO
      END DO
C
C compute bond restraint energy and derivatives w.r.t. b-factors
      DO I=1,NBR
      RENER=RENER
     &      +(WMAIN(XRATOF(IBR(I)))-WMAIN(XRATOF(JBR(I))))**2 *FOUR/
     &        (BSIGMA(XRATOF(IBR(I)))+BSIGMA(XRATOF(JBR(I))))**2
      DB=EIGHT*(WMAIN(XRATOF(IBR(I)))-WMAIN(XRATOF(JBR(I)))) /
     &        (BSIGMA(XRATOF(IBR(I)))+BSIGMA(XRATOF(JBR(I))))**2
      TEMP(IBR(I))=TEMP(IBR(I))+DB
      TEMP(JBR(I))=TEMP(JBR(I))-DB
      END DO
C
C compute angle restraint energy and derivatives w.r.t. b-factors
      DO I=1,NTR
      RENER=RENER
     &      +(WMAIN(XRATOF(ITR(I)))-WMAIN(XRATOF(KTR(I))))**2 *FOUR/
     &        (ASIGMA(XRATOF(ITR(I)))+ASIGMA(XRATOF(KTR(I))))**2
      DB=EIGHT*(WMAIN(XRATOF(ITR(I)))-WMAIN(XRATOF(KTR(I))))/
     &        (ASIGMA(XRATOF(ITR(I)))+ASIGMA(XRATOF(KTR(I))))**2
      TEMP(ITR(I))=TEMP(ITR(I))+DB
      TEMP(KTR(I))=TEMP(KTR(I))-DB
      END DO
C
      IF (LNCSRE) THEN
C
C compute NCS B-factor restraints
      DO K=1,NGROUP
      CALL BNCS(RENER,TEMP,NATOM,WMAIN,IMOVE,FLAG,XRATOF,XRNATF,
     &          NATNCS(K),NUMEQV(K),HEAP(HPSPNT(K)),SIGBNC(K),NNCSX)
      END DO
      END IF
C
C if required, compare the gradient of restraints and diffraction
C terms
      IF (RWGHT.LT.ZERO) THEN
      DNORM1=ZERO
      DNORM2=ZERO
      DO IAT=1,XRNATF
      IF (IMOVE(XRATOF(IAT)).EQ.0) THEN
      DNORM1=DNORM1+XRDT(XRATOF(IAT))**2
      DNORM2=DNORM2+TEMP(IAT)**2
      END IF
      END DO
      DNORM1=SQRT(DNORM1)
      DNORM2=SQRT(DNORM2)
C
C  Set the scale factor: we want (RWGHT*DNORM2)/DNORM1 = 1.0
      IF (DNORM2.GT.RSMALL) THEN
      RWGHT=DNORM1/DNORM2
      ELSE
      RWGHT=ZERO
      END IF
      WRITE(6,'(A,F15.5)')
     &' XBFUNC: restraint weight has been set to ',RWGHT
      END IF
C
C accumulate derivatives and energy
      DO IAT=1,XRNATF
      XRDT(XRATOF(IAT))=XRDT(XRATOF(IAT))+RWGHT*TEMP(IAT)
      END DO
      ELSE
      RENER=ZERO
      END IF
      TARGET=BENER+RWGHT*RENER
C
C store parameters for the cycle with the lowest target value
      IF (TARGET.LT.MTARG) THEN
      I=0
      DO IAT=1,XRNATF
      IF (IMOVE(XRATOF(IAT)).EQ.0) THEN
      I=I+1
      MWMAIN(XRATOF(IAT))=WMAIN(XRATOF(IAT))
      END IF
      END DO
      END IF
C
C Fill the gradient vector, apply B-factor truncation correction
      DO I=1,NCGRP
      GRDENT(I)=ZERO
      END DO
      DO I=1,NCGRP
      DO IAT=POINT(I)+1,POINT(I+1)
      IF (WMAIN(LIST(IAT)).GT.BFMIN.AND.
     &    WMAIN(LIST(IAT)).LT.BFMAX) THEN
      GRDENT(I)=GRDENT(I)+XRDT(LIST(IAT))
      END IF
      END DO
      END DO
C
      GRAD=ZERO
      DO I=1,NCGRP
      GRAD=GRAD+GRDENT(I)**2
      END DO
      GRAD=SQRT(GRAD/NCGRP)
C
C print some stuff
      IF (ICYCLE.GE.0) THEN
      ICYCLE=ICYCLE+1
      WRITE(6,'(A,I6,A,F10.4,A)')
     & ' --------------- cycle=',ICYCLE,
     & ' --------------------------------------------------'
      WRITE(6,'(A,E10.3,A,E10.3,A,E10.3,A,E10.3,A)')
     & ' | Target=',TARGET,
     & ' E(XREF)=',BENER,'  E(BRES)=',RENER,'  E(grad)=',GRAD,
     & '|'
      WRITE(6,'(2A)') ' --------------------------------------------',
     & '-----------------------------------'
      END IF
C
C store the cycle number and value for the lowest energy cycle
      IF (TARGET.LT.MTARG) THEN
      MCYCLE=ICYCLE
      MTARG=TARGET
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XBSET(IBR,JBR,ITR,KTR,FLAG,XRATOF,TEMP,BSIGMA,ASIGMA)
C
C routine sets up bond and angle B-factor restraint list
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'heap.inc'
      INTEGER IBR(*), JBR(*), ITR(*), KTR(*), FLAG(*)
      INTEGER XRATOF(*)
      DOUBLE PRECISION TEMP(*), BSIGMA(*), ASIGMA(*)
C local
      INTEGER I, K, NNCSX
      DOUBLE PRECISION DELTA, RENER
C parameters
      DOUBLE PRECISION ZERO, FOUR
      PARAMETER (ZERO=0.0D0, FOUR=4.0D0)
C begin
C
C fill atom flags
      DO I=1,NATOM
      FLAG(I)=0
      END DO
      DO I=1,XRNATF
      IF (IMOVE(XRATOF(I)).EQ.0) FLAG(XRATOF(I))=I
      END DO
C
C fill the bond list
      NBR=0
      DO I=1,NBOND
      IF (FLAG(IB(I)).NE.0.AND.FLAG(JB(I)).NE.0) THEN
      NBR=NBR+1
      IBR(NBR)=FLAG(IB(I))
      JBR(NBR)=FLAG(JB(I))
      END IF
      END DO
C
C fill the angle list
      NTR=0
      DO I=1,NTHETA
      IF (FLAG(IT(I)).NE.0.AND.FLAG(KT(I)).NE.0) THEN
      NTR=NTR+1
      ITR(NTR)=FLAG(IT(I))
      KTR(NTR)=FLAG(KT(I))
      END IF
      END DO
C
C initial statistics
      DELTA=ZERO
      DO I=1,NBR
      DELTA=DELTA+(WMAIN(XRATOF(IBR(I)))-WMAIN(XRATOF(JBR(I))))**2
     &       *FOUR/(BSIGMA(XRATOF(IBR(I)))+BSIGMA(XRATOF(JBR(I))))**2
      END DO
      IF (NBR.GT.0) DELTA=SQRT(DELTA/NBR)
      WRITE(6,'(A,I5,A,/,A,F10.5)')
     & ' XBSET: There are ',NBR,' bond restraints on B-factors',
     & '        Initial sum (bi-bj)/bsigma^2: ',DELTA
C
      DELTA=ZERO
      DO I=1,NTR
      DELTA=DELTA+(WMAIN(XRATOF(ITR(I)))-WMAIN(XRATOF(KTR(I))))**2
     &       *FOUR/(ASIGMA(XRATOF(ITR(I)))+ASIGMA(XRATOF(KTR(I))))**2
      END DO
      IF (NTR.GT.0) DELTA=SQRT(DELTA/NTR)
      WRITE(6,'(A,I5,A,/,A,F10.5)')
     & ' XBSET: There are ',NTR,' angle restraints on B-factors',
     & '        Initial sum (bi-bj)/asigma^2: ',DELTA
C
      IF (LNCSRE) THEN
C
C compute initial NCS B-factor restraints (setup a fake call to BNCS)
      DO K=1,NGROUP
      RENER=ZERO
      DO I=1,XRNATF
      TEMP(I)=ZERO
      END DO
      CALL BNCS(RENER,TEMP,NATOM,WMAIN,IMOVE,FLAG,XRATOF,XRNATF,
     &          NATNCS(K),NUMEQV(K),HEAP(HPSPNT(K)),SIGBNC(K),NNCSX)
      DELTA=SQRT(RENER/NNCSX)*SIGBNC(K)
      WRITE(6,'(A,I5,A,A,I3,/,A,F10.5,A,F10.5)')
     & ' XBSET: There are ',NNCSX,' NCS restraints on B-factors',
     & ' for group ',K,
     & '        Initial rms deviation of B-factors is: ',DELTA,
     & ' target=',SIGBNC(K)
      END DO
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XBRMS(IBR,JBR,ITR,KTR,FLAG,XRATOF,TEMP,
     &          BGROUP,AGROUP,NBGRP,NAGRP,BRMSD,ARMSD,NBRG,NTRG)
C
C Evaluate RMS deviations of B values for bond and angle
C restraints.
C
C Symbols are declared as follows:
C
C  $BRMS_BOND     overall rmsd for bond B restraints
C  $BN_BOND       total number of bond B restraints
C  $NGROUP_BOND   number of groups of bond restraints
C  $BRMS_BOND_$k  rmsd for bond B restraints for group k
C  $BN_BOND_$k    number of bond B restraints for group k
C  $NGROUP_ANGL   number of groups of angle restraints
C  $BRMS_ANGL     overall rmsd for angle B restraints
C  $BN_ANGL       total number of angle restraints
C  $BRMS_ANGL_$k  rmsd for angle B restraints for group k
C  $BN_ANGL_$k    number of angle B restraints for group k
C  $NGROUP_NCS    number of groups of NCS restraints
C  $BRMS_NCS_$k   rmsd for NCS B restraints for group k
C  $BN_NCS_$k     number of NCS B restraints for group k
C
C Authors: J.-S. Jiang and Axel T. Brunger
C =======================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'heap.inc'
      INTEGER IBR(*), JBR(*), ITR(*), KTR(*), FLAG(*)
      INTEGER XRATOF(*)
      DOUBLE PRECISION TEMP(*)
      INTEGER BGROUP(*), AGROUP(*)
      INTEGER NBGRP, NAGRP
      DOUBLE PRECISION BRMSD(*), ARMSD(*)
      INTEGER NBRG(*), NTRG(*)
C local
      INTEGER I, K, K1, NNCSX
      DOUBLE PRECISION DELTA, RENER
      DOUBLE COMPLEX DBCOMP
      CHARACTER*32 WD, CTEMP
      INTEGER WDMAX, WDLEN
      PARAMETER (WDMAX=20)
C parameters
      DOUBLE PRECISION ZERO, FOUR
      PARAMETER (ZERO=0.0D0, FOUR=4.0D0)
C begin
C
C fill atom flags
      DO I=1,NATOM
      FLAG(I)=0
      END DO
      DO I=1,XRNATF
      IF (IMOVE(XRATOF(I)).EQ.0) FLAG(XRATOF(I))=I
      END DO
C
C fill the bond list
      NBR=0
      DO I=1,NBOND
      IF (FLAG(IB(I)).NE.0.AND.FLAG(JB(I)).NE.0) THEN
      NBR=NBR+1
      IBR(NBR)=FLAG(IB(I))
      JBR(NBR)=FLAG(JB(I))
      END IF
      END DO
C
C fill the angle list
      NTR=0
      DO I=1,NTHETA
      IF (FLAG(IT(I)).NE.0.AND.FLAG(KT(I)).NE.0) THEN
      NTR=NTR+1
      ITR(NTR)=FLAG(IT(I))
      KTR(NTR)=FLAG(KT(I))
      END IF
      END DO
C
C compute rms deviations on the bond list
      DELTA=NBR
      CALL DECLAR('BN_BOND','DP',' ',DBCOMP,DELTA)
      WRITE(6,'(A,I8,A)')
     & ' XBRMS: total number of bond B restraints = ',NBR,
     & ' ( $BN_BOND   )'
      DELTA=ZERO
      DO I=1,NBR
      DELTA=DELTA+(WMAIN(XRATOF(IBR(I)))-WMAIN(XRATOF(JBR(I))))**2
      END DO
      IF (NBR.GT.0) DELTA=SQRT(DELTA/NBR)
      CALL DECLAR('BRMS_BOND','DP',' ',DBCOMP,DELTA)
      WRITE(6,'(A,F8.3,A)')
     & ' XBRMS: overall rmsd for bond B restraints = ',DELTA,
     & ' ( $BRMS_BOND )'
C
C bond restraints for groups
      DELTA=NBGRP
      CALL DECLAR('NGROUP_BOND','DP',' ',DBCOMP,DELTA)
      IF (NBGRP.GT.0) THEN
      DO K=1,NBGRP
      BRMSD(K)=ZERO
      NBRG(K)=0
      END DO
      DO I=1,NBR
      K=BGROUP(XRATOF(IBR(I)))
      K1=BGROUP(XRATOF(JBR(I)))
      K=MAX(K,K1)
      IF (K.GT.0) THEN
      NBRG(K)=NBRG(K)+1
      BRMSD(K)=BRMSD(K)+(WMAIN(XRATOF(IBR(I)))
     &     -WMAIN(XRATOF(JBR(I))))**2
      END IF
      END DO
      DO K=1,NBGRP
      IF (NBRG(K).GT.0) THEN
      BRMSD(K)=SQRT(BRMSD(K)/NBRG(K))
      ELSE
      BRMSD(K)=ZERO
      END IF
      CALL ENCODI(K,WD,WDMAX,WDLEN)
      CTEMP='BN_BOND_'//WD(1:WDLEN)
      DELTA=NBRG(K)
      CALL DECLAR(CTEMP(1:WDLEN+8),'DP',' ',DBCOMP,DELTA)
      WRITE(6,'(A,I3,A,I8,A,A,A)')
     & ' XBRMS: number of bond B restraints for group ',K,
     & ' = ',NBRG(K),
     & ' ( $',CTEMP(1:WDLEN+8),'   )'
      CTEMP='BRMS_BOND_'//WD(1:WDLEN)
      CALL DECLAR(CTEMP(1:WDLEN+10),'DP',' ',DBCOMP,BRMSD(K))
      WRITE(6,'(A,I3,A,F8.3,A,A,A)')
     & ' XBRMS: rmsd for bond B restraints for group ',K,
     & ' = ',BRMSD(K),
     & ' ( $',CTEMP(1:WDLEN+10),' )'
      END DO
      END IF
C
C compute rms deviations on the angle list
      DELTA=NTR
      CALL DECLAR('BN_ANGL','DP',' ',DBCOMP,DELTA)
      WRITE(6,'(A,I8,A)')
     & ' XBRMS: total number of angle B restraints = ',NTR,
     & ' ( $BN_ANGL   )'
      DELTA=ZERO
      DO I=1,NTR
      DELTA=DELTA+(WMAIN(XRATOF(ITR(I)))-WMAIN(XRATOF(KTR(I))))**2
      END DO
      IF (NTR.GT.0) DELTA=SQRT(DELTA/NTR)
      CALL DECLAR('BRMS_ANGL','DP',' ',DBCOMP,DELTA)
      WRITE(6,'(A,F8.3,A)')
     & ' XBRMS: overall rmsd for angle B restraints = ',DELTA,
     & ' ( $BRMS_ANGL )'
C
C angle restraints for groups
      DELTA=NAGRP
      CALL DECLAR('NGROUP_ANGL','DP',' ',DBCOMP,DELTA)
      IF (NAGRP.GT.0) THEN
      DO K=1,NAGRP
      ARMSD(K)=ZERO
      NTRG(K)=0
      END DO
      DO I=1,NTR
      K=AGROUP(XRATOF(ITR(I)))
      K1=AGROUP(XRATOF(KTR(I)))
      K=MAX(K,K1)
      IF (K.GT.0) THEN
      NTRG(K)=NTRG(K)+1
      ARMSD(K)=ARMSD(K)+
     &    (WMAIN(XRATOF(ITR(I)))-WMAIN(XRATOF(KTR(I))))**2
      END IF
      END DO
      DO K=1,NAGRP
      IF (NTRG(K).GT.0) THEN
      ARMSD(K)=SQRT(ARMSD(K)/NTRG(K))
      ELSE
      ARMSD(K)=ZERO
      END IF
      CALL ENCODI(K,WD,WDMAX,WDLEN)
      CTEMP='BN_ANGL_'//WD(1:WDLEN)
      DELTA=NTRG(K)
      CALL DECLAR(CTEMP(1:WDLEN+8),'DP',' ',DBCOMP,DELTA)
      WRITE(6,'(A,I3,A,I8,A,A,A)')
     & ' XBRMS: number of angle B restraints for group ',K,
     &  ' = ',NTRG(K),
     & ' ( $',CTEMP(1:WDLEN+8),'   )'
      CTEMP='BRMS_ANGL_'//WD(1:WDLEN)
      CALL DECLAR(CTEMP(1:WDLEN+10),'DP',' ',DBCOMP,ARMSD(K))
      WRITE(6,'(A,I3,A,F8.3,A,A,A)')
     & ' XBRMS: rmsd for angle B restraints for group ',K,
     &  ' = ',ARMSD(K),' ( $',CTEMP(1:WDLEN+10),' )'
      END DO
      END IF
C
C compute rms deviations on the NCS list
      DELTA=ZERO
      CALL DECLAR('NGROUP_NCS','DP',' ',DBCOMP,DELTA)
C
      IF (LNCSRE) THEN
C
      DELTA=NGROUP
      CALL DECLAR('NGROUP_NCS','DP',' ',DBCOMP,DELTA)
      DO K=1,NGROUP
      RENER=ZERO
      DO I=1,XRNATF
      TEMP(I)=ZERO
      END DO
      CALL BNCS(RENER,TEMP,NATOM,WMAIN,IMOVE,FLAG,XRATOF,XRNATF,
     &          NATNCS(K),NUMEQV(K),HEAP(HPSPNT(K)),SIGBNC(K),NNCSX)
      CALL ENCODI(K,WD,WDMAX,WDLEN)
      CTEMP='BN_NCS_'//WD(1:WDLEN)
      DELTA=NNCSX
      CALL DECLAR(CTEMP(1:WDLEN+7),'DP',' ',DBCOMP,DELTA)
      WRITE(6,'(A,I3,A,I8,A,A,A)')
     & ' XBRMS: number of NCS B restraints for group ',K,' = ',NNCSX,
     & ' ( $',CTEMP(1:WDLEN+7),'   )'
      DELTA=SQRT(RENER/NNCSX)*SIGBNC(K)
      CTEMP='BRMS_NCS_'//WD(1:WDLEN)
      CALL DECLAR(CTEMP(1:WDLEN+9),'DP',' ',DBCOMP,DELTA)
      WRITE(6,'(A,I3,A,F8.3,A,A,A)')
     & ' XBRMS: rmsd for NCS B restraints for group ',K,' = ',DELTA,
     & ' ( $',CTEMP(1:WDLEN+9),' )'
      END DO
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XOVERB(XRH,XRK,XRL,TSEL)
C
C This routine refines an k exp(-B s^2/4) scale factor for
C the reciprocal space object with specified name (default: fcalc)
C in order to minimize the specified TARGET.  Uses
C conjugate gradient minimization.  dtarget(<object_name>) must be defined.
C
C Overall scale factor (QKSCALE=TRUE) and/or B values (QBSCALE=TRUE)
C are refined.
C
C B-refinement can be isotropic, anisotropic (including isotropic), or
C anisotropic while contraining the overall isotropic B-value to the
C initial value.
C
C Scale factors (k, Bs) are returned as symbols.
C
C Author: A. T. Brunger
C =====================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'xropti.inc'
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER TSEL(*)
C local
      INTEGER DIM
      PARAMETER (DIM=7)
      DOUBLE PRECISION PARAM(DIM), GRDENT(DIM)
      DOUBLE PRECISION DROP, TOLER
      DOUBLE PRECISION STEP, RAD
      INTEGER NSTEP
      LOGICAL QDEBUG
      DOUBLE COMPLEX DBCOMP
      INTEGER IFCALC, I
      CHARACTER*(WDMAX) OVERNAME
      INTEGER FNDIOBJ
      EXTERNAL FNDIOBJ
      CHARACTER*6 MINTYPE
      INTEGER MAXFEV
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, FORTY, DM5, SMAL
      DOUBLE PRECISION THREE
      PARAMETER (ZERO=0.0D0, TWO=2.0D0, FOUR=4.0D0, FORTY=0.4D0)
      PARAMETER (SMAL=0.001D0, DM5=1.0D-5, ONE=1.0D0)
      PARAMETER (RAD=PI/180.0D0,THREE=3.0D0)
C begin
C
C set default values
      XSCK=ONE
      XSCB=ZERO
      XSCB11=ZERO
      XSCB22=ZERO
      XSCB33=ZERO
      XSCB12=ZERO
      XSCB13=ZERO
      XSCB23=ZERO
      RESTRC='ALL '
      NSTEP=20
      MAXFEV=10
      DROP=FORTY
      STEP=SMAL
      TOLER=DM5
      QDEBUG=.FALSE.
      OVERNAME='FCALC'
C
C minimization method
      MINTYPE='POWELL'
C
      QKSCALE=.FALSE.
      QBSCALE=.TRUE.
      QANISO=.FALSE.
      QISO=.TRUE.
C
      CALL PUSEND ('OVERALL_B>')
C
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('OVERALL_B>')
      CALL MISCOM('OVERALL_B>',USED)
      IF (.NOT. USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-optimize-overall')
C
      ELSE IF (WD(1:4).EQ.'NSTE') THEN
      CALL NEXTI('Number of steps= ',NSTEP)
      ELSE IF (WD(1:4).EQ.'MAXF') THEN
      CALL NEXTI('MAXFev=',MAXFEV)
      ELSE IF (WD(1:4).EQ.'DROP') THEN
      CALL NEXTF('DROP=',DROP)
      ELSE IF (WD(1:4).EQ.'TOLE') THEN
      CALL NEXTF('TOLErance=',TOLER)
      ELSE IF (WD(1:4).EQ.'METH') THEN
      CALL NEXTST('METHod=',MINTYPE)
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' ---------------------------',
     & 'overall-B-refinement-parameters---------------------'
      WRITE(6,'(A,I6,A,F12.5,A,F12.5,A,F8.3)')
     &' | NSTEp=',NSTEP,' DROP=',DROP,' TOLErance=',TOLER
      WRITE(6,'(A,I5)') ' | MAXFevaluations=',MAXFEV
      WRITE(6,'(2A)') ' | METHod=',MINTYPE
      IF (QANISO) THEN
      WRITE(6,'(6(A,F5.2),A,F8.4,/,A)')
     &  ' | B11=',XSCB11,' B22=',XSCB22,
     &  ' B33=',XSCB33,' B12=',XSCB12,' B13=',XSCB13,
     &  ' B23=',XSCB23,'  K=',XSCK,
     &  ' | RESTriction=',RESTRC
       ELSE
      WRITE(6,'(A,F10.5,A,F12.5)') ' | K=',XSCK,'  B=',XSCB
      END IF
      WRITE(6,'(2A)') ' ---------------------------',
     & '----------------------------------------------------'
C
      ELSE IF (WD(1:WDLEN).EQ.'NAME') THEN
      CALL NEXTST('NAME=',OVERNAME)
      ELSE IF (WD(1:4).EQ.'KSCA') THEN
      CALL NEXTLO('KSCAle=',QKSCALE)
      ELSE IF (WD(1:4).EQ.'BSCA') THEN
      CALL NEXTLO('BSCALe=',QBSCALE)
      ELSE IF (WD(1:4).EQ.'K   ') THEN
      CALL NEXTF('K=',XSCK)
      ELSE IF (WD(1:4).EQ.'B   ') THEN
      CALL NEXTF('B=',XSCB)
      ELSE IF (WD(1:4).EQ.'B11 ') THEN
      CALL NEXTF('B11=',XSCB11)
      ELSE IF (WD(1:4).EQ.'B22 ') THEN
      CALL NEXTF('B22=',XSCB22)
      ELSE IF (WD(1:4).EQ.'B33 ') THEN
      CALL NEXTF('B33=',XSCB33)
      ELSE IF (WD(1:4).EQ.'B12 ') THEN
      CALL NEXTF('B12=',XSCB12)
      ELSE IF (WD(1:4).EQ.'B13 ') THEN
      CALL NEXTF('B13=',XSCB13)
      ELSE IF (WD(1:4).EQ.'B23 ') THEN
      CALL NEXTF('B23=',XSCB23)
      ELSE IF (WD(1:4).EQ.'STEP') THEN
      CALL NEXTF('STEP=',STEP)
      ELSE IF (WD(1:4).EQ.'DEBU') THEN
      CALL NEXTLO('DEBUg=',QDEBUG)
      ELSE IF (WD(1:4).EQ.'ANIS') THEN
      CALL NEXTLO('ANISotropic=',QANISO)
      ELSE IF (WD(1:4).EQ.'ISOT') THEN
      CALL NEXTLO('ISOTropic=',QISO)
      ELSE IF (WD(1:4).EQ.'REST') THEN
      CALL NEXTA4('RESTriction=',RESTRC)
      ELSE
      CALL CHKEND('OVERALL_B>',DONE)
      END IF
      END IF
      END DO
      DONE = .FALSE.
C
      IF (QBSCALE.AND..NOT.QANISO.AND..NOT.QISO) THEN
      WRITE(6,'(A)')
     & ' %XOVERB-ERR: Inconsistent options: QBSCAle=TRUE, QISO=FALSE',
     & '              and QANISO=FALSE.  QBSCAle set to FALSE.'
      QBSCALE=.FALSE.
      END IF
C
      IF (QKSCALE.OR.QBSCALE) THEN
      IF (XRNREF.GT.0) THEN
C
C get pointer for reciprocal space array
      IFCALC = FNDIOBJ(XSFNUM, XSFNAM, OVERNAME)
      IF (IFCALC .LT. 1) THEN
      WRITE(6, '(3A)')
     &  ' %XOFUNC-ERR: reciprocal space object "',
     &  'FCALC', '" not declared.'
      CALL WRNDIE(0,'XOFUNC', 'object undeclared.')
      IFCALC = -1
      END IF
      CALL CHKSFTYPE(IFCALC, XSFNAM, XSFTYPE, 'COMP', 'XOFUNC')
      IF (IFCALC.GT.0) THEN
      IF (HPSF(IFCALC).EQ.0) THEN
      WRITE(6, '(3A)')
     &  ' %XOFUNC-ERR: reciprocal space object "',OVERNAME,
     &  '" has not been defined.'
      CALL WRNDIE(0,'XOFUNC', 'object undefined.')
      ELSE
      OVERPFCALC=HPSF(IFCALC)
C
C look up which derivative (dtarget) definition corresponds
C to FCALC
C check if ASSOciate object is defined
      OVERII=0
      DO I=1,IHPFLAG
      IF (XASSOC(I).EQ.OVERNAME) THEN
      OVERII=I
      END IF
      END DO
C
      IF (OVERII.EQ.0) THEN
      WRITE(6,'(3A)') ' %XOFUNC-ERR: Object ',OVERNAME,
     &                ' has no dtarget expression associated.'
      CALL WRNDIE(-5,'XRAY',
     & 'object not associated with dtarget expression.')
      END IF
      END IF
      END IF
C
      IF (IFCALC.GT.0.AND.OVERPFCALC.GT.0.AND.OVERII.GT.0) THEN
C call the optimization routine
      CALL XOISOT(NSTEP,PARAM,GRDENT,DROP,STEP,TOLER,QDEBUG,
     &            ERROR,XRH,XRK,XRL,MINTYPE,MAXFEV)
C
      IF (.NOT.QDEBUG.AND..NOT.ERROR) THEN
C
C declare symbols
      CALL DECLAR('KSCALE','DP',' ',DBCOMP,XSCK)
      WRITE(6,'(A,F10.5,A)')
     & ' XOVERB: overall scale K ( $KSCALE=',XSCK,' )'
      IF (QBSCALE.AND.QANISO) THEN
      CALL DECLAR('B_11','DP',' ',DBCOMP,XSCB11)
      CALL DECLAR('B_22','DP',' ',DBCOMP,XSCB22)
      CALL DECLAR('B_33','DP',' ',DBCOMP,XSCB33)
      CALL DECLAR('B_12','DP',' ',DBCOMP,XSCB12)
      CALL DECLAR('B_13','DP',' ',DBCOMP,XSCB13)
      CALL DECLAR('B_23','DP',' ',DBCOMP,XSCB23)
      WRITE(6,'(A,6(/,A,F10.5,A))')
     & ' XOVERB: overall anisotropic B values: ',
     & '         ( $B_11=',XSCB11,' )',
     & '         ( $B_22=',XSCB22,' )',
     & '         ( $B_33=',XSCB33,' )',
     & '         ( $B_12=',XSCB12,' )',
     & '         ( $B_13=',XSCB13,' )',
     & '         ( $B_23=',XSCB23,' )'
      WRITE(6,'(A)')
     & ' XOVERB: To apply use',
     & '  $KSCALE *exp(-( $B_11*h*h*$astar*$astar',
     & '                 +$B_22*k*k*$bstar*$bstar',
     & '                 +$B_33*l*l*$cstar*$cstar',
     & '                 +2.0*$B_12*h*k*$astar*$bstar',
     & '                 +2.0*$B_13*h*l*$astar*$cstar',
     & '                 +2.0*$B_23*k*l*$bstar*$cstar',
     & '                 )/4.0)) * <structure-factor-array>'
      ELSE
      CALL DECLAR('BSCALE','DP',' ',DBCOMP,XSCB)
      WRITE(6,'(A,F10.5,A)')
     & ' XOVERB: overall B-factor ( $BSCALE=',XSCB,' )'
      WRITE(6,'(A)')
     & ' XOVERB: To apply use',
     & '  $KSCALE*exp(- $BSCALE s^2 )/4.0))*<structure-factor-array>'
      END IF
C
      END IF
      END IF
      END IF
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE XOISOT(NSTEP,PARAM,GRDENT,DROP,STEP,TOLER,DEBUG,
     &            ERROR,XRH,XRK,XRL,MINTYPE,MAXFEV)
C
C Refines an overall isotropic B-factor b or anisotropic tensor
C (b11, b22, b33, b12, b13, b23) such that the target defined
C by the general target routine XTARGETS is minimal.
C
C Author: A. T. Brunger
C =====================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'consta.inc'
      INTEGER NSTEP
      DOUBLE PRECISION PARAM(*), GRDENT(*)
      DOUBLE PRECISION DROP, STEP, TOLER
      LOGICAL DEBUG, ERROR
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER DIM
      DOUBLE PRECISION TARGET, NUMGRD, TOL2
      INTEGER IER, J
      CHARACTER*(*) MINTYPE
      INTEGER MAXFEV
C local
      INTEGER M, NPRINT, IW
C pointer
      INTEGER WORK, DIAG, BEST
C parameter
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
C external
      EXTERNAL XOFUNC, XOFUNC2
C begin
C
C
C initialize the lowest target value
      MTARG=R4BIG
      MCYCLE=0
C
C
C copy those parameters that are being refined
      DIM=0
      IF (QKSCALE) THEN
      DIM=DIM+1
      PARAM(DIM)=XSCK
      END IF
      IF (QBSCALE) THEN
      IF (QANISO) THEN
C
      DIM=DIM+1
      PARAM(DIM)=XSCB11
      DIM=DIM+1
      PARAM(DIM)=XSCB22
      DIM=DIM+1
      PARAM(DIM)=XSCB33
      DIM=DIM+1
      PARAM(DIM)=XSCB12
      DIM=DIM+1
      PARAM(DIM)=XSCB13
      DIM=DIM+1
      PARAM(DIM)=XSCB23
C
      ELSE
      DIM=DIM+1
      PARAM(DIM)=XSCB
      END IF
      END IF
C
      IF (DEBUG) THEN
C
C check consistency of first derivatives
      ICYCLE=-1
      DO J=1,DIM
      PARAM(J)=PARAM(J)+STEP
      CALL XOFUNC(DIM,PARAM,TARGET,GRDENT)
      NUMGRD=TARGET
      PARAM(J)=PARAM(J)-TWO*STEP
      CALL XOFUNC(DIM,PARAM,TARGET,GRDENT)
      NUMGRD=NUMGRD-TARGET
      NUMGRD=NUMGRD/(TWO*STEP)
      PARAM(J)=PARAM(J)+STEP
      CALL XOFUNC(DIM,PARAM,TARGET,GRDENT)
      WRITE(6,'(A,I6,A,E12.5,A,E12.5,A,E12.5)')
     &' XBDEBUG: parameter # ',J,' analyt.=',GRDENT(J),
     &' finite=',NUMGRD
      END DO
      ELSE
C
C carry out minimzation.
C
C set cycle counter to zero
      ICYCLE=0
C
C store current values of parameters
      XSCKO=XSCK
      XSCBO=XSCB
      XSCB11O=XSCB11
      XSCB22O=XSCB22
      XSCB33O=XSCB33
      XSCB12O=XSCB12
      XSCB13O=XSCB13
      XSCB23O=XSCB23
C
C Set tolerance to the sum of squares
      TOL2=(TOLER**2)*DIM
C
      ERROR=.FALSE.
      IF (MINTYPE.EQ.'POWELL') THEN
C
C POWELL-method minimizer
      WORK=ALLHP(IREAL8(DIM*6))
      CALL ZXCGR(XOFUNC, DIM, TOL2, NSTEP, DROP, PARAM,
     &           GRDENT, TARGET, HEAP(WORK), IER)
      CALL FREHP(WORK,IREAL8(DIM*6))
C check termniation conditions
      IF (IER.EQ.0) THEN
      WRITE(6,'(A)')' XBREF: gradient converged'
      ELSE IF (IER.EQ.129) THEN
      WRITE(6,'(A)')' XBREF: Line search terminated'
      ELSE IF (IER.EQ.130) THEN
      WRITE(6,'(A)')' %XBREF-ERR: Search Direction Uphill'
      ERROR=.TRUE.
      ELSE IF (IER.EQ.131) THEN
      WRITE(6,'(A)')' XBREF: step limit reached'
      ELSE IF (IER.EQ.132) THEN
      WRITE(6,'(A)')' %XBREF-ERR: Failure to reduce E'
      ERROR=.TRUE.
      END IF
      ELSE IF (MINTYPE.EQ.'LBFGS') THEN
C
C LBFGS-method minimizer
      M=6
      NPRINT=1
      IW=DIM*(2*M+1)+2*M
      DIAG=ALLHP(IREAL8(DIM))
      BEST=ALLHP(IREAL8(DIM))
      WORK=ALLHP(IREAL8(IW))
      CALL ZXLBFGS(XOFUNC2,DIM,M,PARAM,GRDENT,TARGET,
     &             .FALSE.,HEAP(DIAG),NPRINT,TOL2,NSTEP,HEAP(WORK),
     &              HEAP(BEST),MAXFEV)
      CALL FREHP(WORK,IREAL8(IW))
      CALL FREHP(BEST,IREAL8(DIM))
      CALL FREHP(DIAG,IREAL8(DIM))
      ELSE
      CALL WRNDIE(-5,'XOISOT','unknown minimization method.')
      END IF
C
C set parameters to the values that produces the lowest target value
      XSCK=XSCKO
      XSCB=XSCBO
      XSCB11=XSCB11O
      XSCB22=XSCB22O
      XSCB33=XSCB33O
      XSCB12=XSCB12O
      XSCB13=XSCB13O
      XSCB23=XSCB23O
      WRITE(6,'(A,I6,A)')
     & ' XOVERB: using the parameters (at cycle ', MCYCLE,
     & ') which produce the minimum target. '
C
      END IF
C
      RETURN
      END
C
C======================================================================
      SUBROUTINE XOFUNC(DIM,PARAM,TARGET,GRDENT)
C
C ZXCGR target function routine for overall B-factor (isotropic or
C anisotropic) minimization.  Front-end routine which passes the pointers
C of some dynamic arrays to the next routine.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INTEGER DIM
      DOUBLE PRECISION PARAM(*), TARGET, GRDENT(*)
C pointers
      INTEGER FCALCL, HPXDERIV
C begin
      FCALCL=ALLHP(ICPLX8(XRNREF))
      HPXDERIV=ALLHP(ICPLX8(XRNREF))
C
      CALL XOFUN2(DIM,PARAM,TARGET,GRDENT,HEAP(HPH),HEAP(HPK),
     &            HEAP(HPL),
     &            HEAP(OVERPFCALC),HEAP(FCALCL),
     &            HEAP(HPXDERIV),HEAP(HPTSEL),OVERII)
      CALL FREHP(HPXDERIV,ICPLX8(XRNREF))
      CALL FREHP(FCALCL,ICPLX8(XRNREF))
      RETURN
      END
C======================================================================
      SUBROUTINE XOFUNC2(DIM,PARAM,TARGET,GRDENT,ITER,QPRINT)
C
C ZXLBFGS target function routine for overall B-factor (isotropic or
C anisotropic) minimization.  Front-end routine which passes the pointers
C of some dynamic arrays to the next routine.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INTEGER DIM
      DOUBLE PRECISION PARAM(*), TARGET, GRDENT(*)
      INTEGER ITER
      LOGICAL QPRINT
C pointers
      INTEGER FCALCL, HPXDERIV
C begin
      FCALCL=ALLHP(ICPLX8(XRNREF))
      HPXDERIV=ALLHP(ICPLX8(XRNREF))
C
      CALL XOFUN2(DIM,PARAM,TARGET,GRDENT,HEAP(HPH),HEAP(HPK),
     &            HEAP(HPL),
     &            HEAP(OVERPFCALC),HEAP(FCALCL),
     &            HEAP(HPXDERIV),HEAP(HPTSEL),OVERII)
      CALL FREHP(HPXDERIV,ICPLX8(XRNREF))
      CALL FREHP(FCALCL,ICPLX8(XRNREF))
      RETURN
      END
C======================================================================
      SUBROUTINE XOFUN2(DIM,PARAM,TARGET,GRDENT,XRH,XRK,XRL,
     &            FCALC,FCALCL,
     &            XDERIV,TSEL,II)
C
C Routine calculates the first derivatives of the target function
C with respect to the overall B-factor(s) (isotropic or anisotropic).
C
C   Target(fobs, G*fcalc+fpart )
C
C G is the Debye-Waller factor
C
C   G = K * exp ( - B * s**2/4.0)
C
C where s = 1/d.
C
C
C Authors: Axel T. Brunger, J.-S. Jiang, J. Kuriyan
C =================================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'xropti.inc'
      INTEGER DIM
      DOUBLE PRECISION PARAM(*), TARGET, GRDENT(*)
      INTEGER XRH(*), XRK(*), XRL(*)
      DOUBLE COMPLEX FCALC(*)
      DOUBLE COMPLEX FCALCL(*)
      DOUBLE COMPLEX XDERIV(*)
      INTEGER TSEL(*), II
C local
      DOUBLE PRECISION SSQ, G, GRAD, S1, S2, S3
      INTEGER REFLCT, HH, KK, LL, I
      DOUBLE PRECISION CABG1,CABG2,CABG3,SABG1,SABG2,SABG3,RAD
      DOUBLE PRECISION CABGS1,CABGS2,CABGS3,ABCS1,ABCS2,ABCS3
      DOUBLE PRECISION DB11, DB22, DB33, DB12, DB13, DB23, TEMP
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, CSMALL, THREE
      PARAMETER(ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, FOUR=4.0D0)
      PARAMETER(RAD=PI/180.0D0, CSMALL=0.001D0,THREE=3.0D0)
C begin
C
C safe current FCALC array
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).NE.0) THEN
      FCALCL(REFLCT)=FCALC(REFLCT)
      END IF
      END DO
C
C copy those parameters that are being refined
      DIM=0
      IF (QKSCALE) THEN
      DIM=DIM+1
      XSCK=PARAM(DIM)
      END IF
      IF (QBSCALE) THEN
      IF (QANISO) THEN
C
      DIM=DIM+1
      XSCB11=PARAM(DIM)
      DIM=DIM+1
      XSCB22=PARAM(DIM)
      DIM=DIM+1
      XSCB33=PARAM(DIM)
      DIM=DIM+1
      XSCB12=PARAM(DIM)
      DIM=DIM+1
      XSCB13=PARAM(DIM)
      DIM=DIM+1
      XSCB23=PARAM(DIM)
C
      ELSE
      DIM=DIM+1
      XSCB=PARAM(DIM)
      END IF
      END IF
C
C compute G*FCALC and store in temporary array
      IF (QBSCALE.AND.QANISO) THEN
C anisotropic b-factor
C compute reciprocal unit cell dimensions
      CABG1=COS(XRCELL(4)*RAD)
      CABG2=COS(XRCELL(5)*RAD)
      CABG3=COS(XRCELL(6)*RAD)
      SABG1=SIN(XRCELL(4)*RAD)
      SABG2=SIN(XRCELL(5)*RAD)
      SABG3=SIN(XRCELL(6)*RAD)
      CABGS1=(CABG2*CABG3-CABG1)/(SABG2*SABG3)
      CABGS2=(CABG3*CABG1-CABG2)/(SABG3*SABG1)
      CABGS3=(CABG1*CABG2-CABG3)/(SABG1*SABG2)
      XRVOL=XRCELL(1)*XRCELL(2)*XRCELL(3)*
     &                   SQRT(ONE+TWO*CABG1*CABG2*CABG3
     &                  -CABG1**2-CABG2**2-CABG3**2)
      ABCS1=XRCELL(2)*XRCELL(3)*SABG1/XRVOL
      ABCS2=XRCELL(1)*XRCELL(3)*SABG2/XRVOL
      ABCS3=XRCELL(1)*XRCELL(2)*SABG3/XRVOL
C restrictions due to symmetry according to the crystal system
      IF (RESTRC.EQ.'ALL ') THEN
      IF (ABCS1.EQ.ABCS2.AND.ABCS2.EQ.ABCS3) THEN
      XSCB11=(XSCB11+XSCB22+XSCB33)/THREE
      XSCB22=XSCB11
      XSCB33=XSCB11
      ELSE IF (ABCS1.EQ.ABCS2) THEN
      XSCB11=(XSCB11+XSCB22)/TWO
      XSCB22=XSCB11
      ELSE IF (ABCS2.EQ.ABCS3) THEN
      XSCB22=(XSCB22+XSCB33)/TWO
      XSCB33=XSCB22
      ELSE IF (ABCS3.EQ.ABCS1) THEN
      XSCB33=(XSCB33+XSCB11)/TWO
      XSCB11=XSCB33
      END IF
      END IF
      IF (RESTRC.EQ.'ALL '.OR.RESTRC.EQ.'OFFD') THEN
      IF (ABS(CABGS1).LT.CSMALL) XSCB23=ZERO
      IF (ABS(CABGS2).LT.CSMALL) XSCB13=ZERO
      IF (ABS(CABGS3).LT.CSMALL) XSCB12=ZERO
C identify Rhombohedral for restriction
      IF (ABCS1.EQ.ABCS2.AND.ABCS2.EQ.ABCS3.AND.
     &    CABGS1.EQ.CABGS2.AND.CABGS2.EQ.CABGS3.AND.
     &    ABS(CABGS1).GT.CSMALL) THEN
      XSCB12=(XSCB12+XSCB13+XSCB23)/THREE
      XSCB13=XSCB12
      XSCB23=XSCB12
      END IF
      END IF
      END IF
C
C zero trace restriction
      IF (.NOT.QISO) THEN
      TEMP=(XSCB11+XSCB22+XSCB33)/THREE
      XSCB11=XSCB11-TEMP
      XSCB22=XSCB22-TEMP
      XSCB33=XSCB33-TEMP
      END IF
C
      IF (QBSCALE.AND.QANISO) THEN
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).NE.0) THEN
      HH=XRH(REFLCT)
      KK=XRK(REFLCT)
      LL=XRL(REFLCT)
      S1= ABCS1*HH
      S2= ABCS2*KK
      S3= ABCS3*LL
      FCALC(REFLCT)=XSCK*EXP(-(
     &  XSCB11*S1**2+XSCB22*S2**2+XSCB33*S3**2+
     &  TWO*(XSCB12*S1*S2+XSCB13*S1*S3+XSCB23*S2*S3))
     &  /FOUR)*FCALCL(REFLCT)
      END IF
      END DO
      ELSE
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).NE.0) THEN
      HH=XRH(REFLCT)
      KK=XRK(REFLCT)
      LL=XRL(REFLCT)
      SSQ=  (XRTR(1,1)*HH + XRTR(2,1)*KK + XRTR(3,1)*LL)**2
     &    + (XRTR(1,2)*HH + XRTR(2,2)*KK + XRTR(3,2)*LL)**2
     &    + (XRTR(1,3)*HH + XRTR(2,3)*KK + XRTR(3,3)*LL)**2
      FCALC(REFLCT)=XSCK*EXP(-XSCB*SSQ/FOUR)
     &      *FCALCL(REFLCT)
      END IF
      END DO
      END IF
C
C compute target function and first derivatives. No printout.
C working set only.
      CALL XTARGETS(.TRUE.,.FALSE.,1,XRE,XDERIV,
     &           HPTSEL,XRNREF,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           DRPNMX,DRPNX,DRPN(1,1,II),DRPNL(1,1,II),
     &           DRPNN(II),DRPNDB(1,1,II),
     &           DRPNMLT(1,II),DRPNLEV(1,II),DRPNTYP(1,II),
     &           DRPNDOM(1,II),DDEPTH(II),
     &           MRPNMX,MRPNX,MRPN,MRPNL,MRPNN,MRPNDB,
     &           MRPNMLT,MRPNLEV,MRPNTYP,MRPNDOM,MDEPTH,
     &           HPH,HPK,HPL,XSFNUM,XSFNAM,
     &           XSFTYPE,HPSF,HPMULT,HPTYPE,XRMREF,
     &           QHERM,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,1,
     &           XRCELL,XRVOL)
C
      TARGET=XRE
C
C store parameters for the cycle with the lowest target value
      IF (TARGET.LT.MTARG) THEN
      XSCKO=XSCK
      XSCBO=XSCB
      XSCB11O=XSCB11
      XSCB22O=XSCB22
      XSCB33O=XSCB33
      XSCB12O=XSCB12
      XSCB13O=XSCB13
      XSCB23O=XSCB23
      END IF
C
C now compute derivatives
      DO I=1,DIM
      GRDENT(I)=ZERO
      END DO
C
      DIM=0
      IF (QKSCALE) THEN
      DIM=DIM+1
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).NE.0) THEN
      IF (QBSCALE.AND.QANISO) THEN
      HH=XRH(REFLCT)
      KK=XRK(REFLCT)
      LL=XRL(REFLCT)
      S1= ABCS1*HH
      S2= ABCS2*KK
      S3= ABCS3*LL
      G=XSCK*EXP(-(
     &  XSCB11*S1**2+XSCB22*S2**2+XSCB33*S3**2+
     &  TWO*(XSCB12*S1*S2+XSCB13*S1*S3+XSCB23*S2*S3))
     &  /FOUR)
      ELSE
      G=XSCK*EXP(-XSCB*SSQ/FOUR)
      END IF
      GRDENT(DIM)=GRDENT(DIM)+G/XSCK*
     &(DBLE(XDERIV(REFLCT))*(DBLE(FCALCL(REFLCT)))
     &+DIMAG(XDERIV(REFLCT))*(DIMAG(FCALCL(REFLCT))))
      END IF
      END DO
      END IF
C
      IF (QBSCALE) THEN
      IF (QANISO) THEN
C anisotropic b-factor
      DB11=ZERO
      DB22=ZERO
      DB33=ZERO
      DB12=ZERO
      DB13=ZERO
      DB23=ZERO
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).NE.0) THEN
      HH=XRH(REFLCT)
      KK=XRK(REFLCT)
      LL=XRL(REFLCT)
      S1= ABCS1*HH
      S2= ABCS2*KK
      S3= ABCS3*LL
      G=XSCK*EXP(-(
     &  XSCB11*S1**2+XSCB22*S2**2+XSCB33*S3**2+
     &  TWO*(XSCB12*S1*S2+XSCB13*S1*S3+XSCB23*S2*S3))
     &  /FOUR)
      DB11=DB11-S1**2*G/FOUR*
     &(DBLE(XDERIV(REFLCT))*(DBLE(FCALCL(REFLCT)))
     &+DIMAG(XDERIV(REFLCT))*(DIMAG(FCALCL(REFLCT))))
      DB22=DB22-S2**2*G/FOUR*
     &(DBLE(XDERIV(REFLCT))*(DBLE(FCALCL(REFLCT)))
     &+DIMAG(XDERIV(REFLCT))*(DIMAG(FCALCL(REFLCT))))
      DB33=DB33-S3**2*G/FOUR*
     &(DBLE(XDERIV(REFLCT))*(DBLE(FCALCL(REFLCT)))
     &+DIMAG(XDERIV(REFLCT))*(DIMAG(FCALCL(REFLCT))))
      DB12=DB12-S1*S2*G/TWO*
     &(DBLE(XDERIV(REFLCT))*(DBLE(FCALCL(REFLCT)))
     &+DIMAG(XDERIV(REFLCT))*(DIMAG(FCALCL(REFLCT))))
      DB13=DB13-S1*S3*G/TWO*
     &(DBLE(XDERIV(REFLCT))*(DBLE(FCALCL(REFLCT)))
     &+DIMAG(XDERIV(REFLCT))*(DIMAG(FCALCL(REFLCT))))
      DB23=DB23-S2*S3*G/TWO*
     &(DBLE(XDERIV(REFLCT))*(DBLE(FCALCL(REFLCT)))
     &+DIMAG(XDERIV(REFLCT))*(DIMAG(FCALCL(REFLCT))))
      END IF
      END DO
C restrictions due to symmetry according to the crystal system
      IF (RESTRC.EQ.'ALL ') THEN
      IF (ABCS1.EQ.ABCS2.AND.ABCS2.EQ.ABCS3) THEN
      DB11=(DB11+DB22+DB33)/THREE
      DB22=DB11
      DB33=DB11
      ELSE IF (ABCS1.EQ.ABCS2) THEN
      DB11=(DB11+DB22)/TWO
      DB22=DB11
      ELSE IF (ABCS2.EQ.ABCS3) THEN
      DB22=(DB22+DB33)/TWO
      DB33=DB22
      ELSE IF (ABCS3.EQ.ABCS1) THEN
      DB33=(DB33+DB11)/TWO
      DB11=DB33
      END IF
      END IF
      IF ((RESTRC.EQ.'ALL '.OR.RESTRC.EQ.'OFFD')) THEN
      IF (ABS(CABGS1).LT.CSMALL) DB23=ZERO
      IF (ABS(CABGS2).LT.CSMALL) DB13=ZERO
      IF (ABS(CABGS3).LT.CSMALL) DB12=ZERO
C identify Rhombohedral for restriction
      IF (ABCS1.EQ.ABCS2.AND.ABCS2.EQ.ABCS3.AND.
     &    CABGS1.EQ.CABGS2.AND.CABGS2.EQ.CABGS3.AND.
     &    ABS(CABGS1).GT.CSMALL) THEN
      DB12=(DB12+DB13+DB23)/THREE
      DB13=DB12
      DB23=DB12
      END IF
      END IF
C
C zero trace restriction
      IF (.NOT.QISO) THEN
      TEMP=(DB11+DB22+DB33)/THREE
      DB11=DB11-TEMP
      DB22=DB22-TEMP
      DB33=DB33-TEMP
      END IF
C
      DIM=DIM+1
      GRDENT(DIM)=DB11
      DIM=DIM+1
      GRDENT(DIM)=DB22
      DIM=DIM+1
      GRDENT(DIM)=DB33
      DIM=DIM+1
      GRDENT(DIM)=DB12
      DIM=DIM+1
      GRDENT(DIM)=DB13
      DIM=DIM+1
      GRDENT(DIM)=DB23
      ELSE
C isotropic b-factor
      DIM=DIM+1
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).NE.0) THEN
      HH=XRH(REFLCT)
      KK=XRK(REFLCT)
      LL=XRL(REFLCT)
      SSQ=  (XRTR(1,1)*HH + XRTR(2,1)*KK + XRTR(3,1)*LL)**2
     &    + (XRTR(1,2)*HH + XRTR(2,2)*KK + XRTR(3,2)*LL)**2
     &    + (XRTR(1,3)*HH + XRTR(2,3)*KK + XRTR(3,3)*LL)**2
      G=XSCK*EXP(-XSCB*SSQ/FOUR)
      GRDENT(DIM)=GRDENT(DIM)-(SSQ/FOUR)*G*
     &(DBLE(XDERIV(REFLCT))*(DBLE(FCALCL(REFLCT)))
     &+DIMAG(XDERIV(REFLCT))*(DIMAG(FCALCL(REFLCT))))
      END IF
      END DO
      END IF
      END IF
C
C print some stuff
      IF (ICYCLE.GE.0) THEN
      ICYCLE=ICYCLE+1
      GRAD=ZERO
      DO I=1,DIM
      GRAD=GRAD+GRDENT(I)**2
      END DO
      GRAD=SQRT(MAX(ZERO,GRAD/DIM))
C
      SSQ=ONE/(FOUR*TWO*PI**2)
      WRITE(6,'(A,I6,A)')
     & ' --------------- cycle=',ICYCLE,
     & ' --------------------------------------------------'
      IF (QBSCALE.AND.QANISO) THEN
      WRITE(6,'(6(A,F8.4),A,/,6(A,F8.4),A,/,2(A,E10.3),A,F10.5,A)')
     & ' |B11=',XSCB11,' B22=',XSCB22,
     &  ' B33=',XSCB33,' B12=',XSCB12,' B13=',XSCB13,
     &  ' B23=',XSCB23,'|',
     & ' |U11=',XSCB11*SSQ,' U22=',XSCB22*SSQ,
     &  ' U33=',XSCB33*SSQ,' U12=',XSCB12*SSQ,
     &  ' U13=',XSCB13*SSQ,' U23=',XSCB23*SSQ,'|',
     & ' | E(XREF)=',TARGET,'  grad(E)=',GRAD,
     & '  K=',XSCK,'                        |'
      ELSE
      WRITE(6,'(A,F10.5,A,F10.5,A,E10.3,A,E10.3,A)')
     & ' | K=',XSCK,'  B=',XSCB,
     & ' E(XREF)=',TARGET,'  grad(E)=',GRAD,
     & '           |'
      END IF
      WRITE(6,'(2A)') ' --------------------------------------------',
     & '-----------------------------------'
      END IF
C
C store the cycle number and value for the lowest energy cycle
      IF (TARGET.LT.MTARG) THEN
      MCYCLE=ICYCLE
      MTARG=TARGET
      END IF
C
C
C re-establish FCALC array
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).NE.0) THEN
      FCALC(REFLCT)=FCALCL(REFLCT)
      END IF
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XGROUP
C
C Front end to parsing routine for group b-factor,
C occupancy, f', and f'' optimization.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
C local
      INTEGER FLAGS, INDEX, MPOINT, MLIST
C begin
      MPOINT=4*(NATOM+2)
      MLIST=4*(NATOM+2)
      FLAGS=ALLHP(INTEG4(NATOM))
      INDEX=ALLHP(INTEG4(NATOM))
      ILIST=ALLHP(INTEG4(MLIST))
      IPOINT=ALLHP(INTEG4(MPOINT))
      ITYPE=ALLHP(INTEG4(MPOINT))
      CALL XGRUP2(HEAP(FLAGS),HEAP(INDEX),
     &             MLIST,HEAP(ILIST),MPOINT,HEAP(IPOINT),
     &             HEAP(ITYPE),
     &             HEAP(HPATOF),HEAP(HPINDF))
      CALL FREHP(ITYPE,INTEG4(MPOINT))
      CALL FREHP(IPOINT,INTEG4(MPOINT))
      CALL FREHP(ILIST,INTEG4(MLIST))
      CALL FREHP(INDEX,INTEG4(NATOM))
      CALL FREHP(FLAGS,INTEG4(NATOM))
      RETURN
      END
C======================================================================
      SUBROUTINE XGRUP2(FLAGS,INDEX,MLIST,LIST,MPOINT,POINT,GTYPE,
     &                   XRATOF,XRINDF)
C
C Parsing routine for group b-factor and occupancy optimization.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INTEGER FLAGS(*), INDEX(*), MLIST, LIST(*), MPOINT, POINT(*)
      INTEGER GTYPE(*), XRATOF(*), XRINDF(*)
C local
      LOGICAL DEBUG
      DOUBLE PRECISION DROP, TOLER
      INTEGER NSTEP, NGROUP, NLIST, GROUP, GRADNT, NFLAGS, IAT
      INTEGER I, SCATTX
      CHARACTER*6 MINTYPE
C parameters
C begin
C
C
C defaults
      NSTEP=0
C Initial step size for conjugate gradient minimizer
      DROP=10.0D0
C Tolerance on gradient
      TOLER=0.000001D0
C minimum allowed B-factor
      BFMIN=2.0D0
C maximum allowed B-factor
      BFMAX=100.D0
C minimum allowed Q
      QFMIN=0.01D0
C maximum allowed Q
      QFMAX=10.D0
C minimum allowed f'
      FPMIN=-100.D0
C maximum allowed f'
      FPMAX=100.D0
C minimum allowed f''
      FDPMIN=-100.D0
C maximum allowed f''
      FDPMAX=100.D0
C minimization method
      MINTYPE='POWELL'
C
      DEBUG=.FALSE.
C
C pointer list
      NGROUP=0
      NLIST=0
      POINT(1)=0
C
      CALL PUSEND('GROUP>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('GROUP>')
      CALL MISCOM('GROUP>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-optimize-group')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NSTE') THEN
      CALL NEXTI('Number of steps= ',NSTEP)
      ELSE IF (WD(1:4).EQ.'DROP') THEN
      CALL NEXTF('DROP=',DROP)
      ELSE IF (WD(1:4).EQ.'TOLE') THEN
      CALL NEXTF('TOLErance=',TOLER)
      ELSE IF (WD(1:4).EQ.'DEBU') THEN
      CALL NEXTLO('DEBUg=',DEBUG)
C======================================================================
      ELSE IF (WD(1:4).EQ.'METH') THEN
      CALL NEXTST('METHod=',MINTYPE)
C======================================================================
      ELSE IF (WD(1:4).EQ.'BFMI'.OR.WD(1:WDLEN).EQ.'BMIN') THEN
      CALL NEXTF('BMIN=',BFMIN)
C======================================================================
      ELSE IF (WD(1:4).EQ.'BFMA'.OR.WD(1:WDLEN).EQ.'BMAX') THEN
      CALL NEXTF('BMAX=',BFMAX)
C======================================================================
      ELSE IF (WD(1:WDLEN).EQ.'QMIN') THEN
      CALL NEXTF('QMIN=',QFMIN)
C======================================================================
      ELSE IF (WD(1:WDLEN).EQ.'QMAX') THEN
      CALL NEXTF('QMAX=',QFMAX)
C======================================================================
      ELSE IF (WD(1:WDLEN).EQ.'FPMIN') THEN
      CALL NEXTF('FPMIN=',FPMIN)
C======================================================================
      ELSE IF (WD(1:WDLEN).EQ.'FPMAX') THEN
      CALL NEXTF('FPMAX=',FPMAX)
C======================================================================
      ELSE IF (WD(1:WDLEN).EQ.'FDPMIN') THEN
      CALL NEXTF('FDPMIN=',FDPMIN)
C======================================================================
      ELSE IF (WD(1:WDLEN).EQ.'FDPMAX') THEN
      CALL NEXTF('FDPMAX=',FDPMAX)
C======================================================================
      ELSE IF (WD(1:WDLEN).EQ.'B') THEN
      CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
      NFLAGS=0
      DO IAT=1,XRNATF
      IF (FLAGS(XRATOF(IAT)).EQ.1) THEN
      NFLAGS=NFLAGS+1
      INDEX(NFLAGS)=IAT
      END IF
      END DO
C
      IF (NFLAGS.GT.0) THEN
C make a new group
      IF (NGROUP+1.GE.MPOINT) THEN
      WRITE(6,'(A)') ' %GROUP-ERR: too many groups specified'
      GOTO 9999
      ELSE
      NGROUP=NGROUP+1
      DO IAT=1,NFLAGS
      IF (NLIST.GE.MLIST) THEN
      WRITE(6,'(A)')
     &' %GROUP-ERR: too many atoms selected--check for group overlaps'
      ELSE
      NLIST=NLIST+1
      LIST(NLIST)=INDEX(IAT)
      END IF
      END DO
      POINT(NGROUP+1)=NLIST
      GTYPE(NGROUP)=1
      END IF
      END IF
C======================================================================
      ELSE IF (WD(1:WDLEN).EQ.'FP') THEN
      CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
      NFLAGS=0
      DO IAT=1,XRNATF
      IF (FLAGS(XRATOF(IAT)).EQ.1) THEN
      NFLAGS=NFLAGS+1
      INDEX(NFLAGS)=IAT
      END IF
      END DO
C
      IF (NFLAGS.GT.0) THEN
C make a new group
      IF (NGROUP+1.GE.MPOINT) THEN
      WRITE(6,'(A)') ' %GROUP-ERR: too many groups specified'
      GOTO 9999
      ELSE
      NGROUP=NGROUP+1
      DO IAT=1,NFLAGS
      IF (NLIST.GE.MLIST) THEN
      WRITE(6,'(A)')
     &' %GROUP-ERR: too many atoms selected--check for group overlaps'
      ELSE
      NLIST=NLIST+1
      LIST(NLIST)=INDEX(IAT)
      END IF
      END DO
      POINT(NGROUP+1)=NLIST
      GTYPE(NGROUP)=3
      END IF
      END IF
C======================================================================
      ELSE IF (WD(1:WDLEN).EQ.'FDP') THEN
      CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
      NFLAGS=0
      DO IAT=1,XRNATF
      IF (FLAGS(XRATOF(IAT)).EQ.1) THEN
      NFLAGS=NFLAGS+1
      INDEX(NFLAGS)=IAT
      END IF
      END DO
C
      IF (NFLAGS.GT.0) THEN
C make a new group
      IF (NGROUP+1.GE.MPOINT) THEN
      WRITE(6,'(A)') ' %GROUP-ERR: too many groups specified'
      GOTO 9999
      ELSE
      NGROUP=NGROUP+1
      DO IAT=1,NFLAGS
      IF (NLIST.GE.MLIST) THEN
      WRITE(6,'(A)')
     &' %GROUP-ERR: too many atoms selected--check for group overlaps'
      ELSE
      NLIST=NLIST+1
      LIST(NLIST)=INDEX(IAT)
      END IF
      END DO
      POINT(NGROUP+1)=NLIST
      GTYPE(NGROUP)=4
      END IF
      END IF
C======================================================================
      ELSE IF (WD(1:WDLEN).EQ.'Q') THEN
      CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
      NFLAGS=0
      DO IAT=1,XRNATF
      IF (FLAGS(XRATOF(IAT)).EQ.1) THEN
      NFLAGS=NFLAGS+1
      INDEX(NFLAGS)=IAT
      END IF
      END DO
C
      IF (NFLAGS.GT.0) THEN
C make a new group
      IF (NGROUP+1.GE.MPOINT) THEN
      WRITE(6,'(A)') ' %GROUP-ERR: too many groups specified'
      GOTO 9999
      ELSE
      NGROUP=NGROUP+1
      DO IAT=1,NFLAGS
      IF (NLIST.GE.MLIST) THEN
      WRITE(6,'(A)')
     &' %GROUP-ERR: too many atoms selected--check for group overlaps'
      ELSE
      NLIST=NLIST+1
      LIST(NLIST)=INDEX(IAT)
      END IF
      END DO
      POINT(NGROUP+1)=NLIST
      GTYPE(NGROUP)=2
      END IF
      END IF
C
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' ---------------------------',
     & 'GROUP-B-factor-refinement-parameters----------------'
      WRITE(6,
     & '(A,I6,A,F12.5,A,F12.5,/,4(A,F12.5),/,4(A,F12.5),/,A,I6)')
     &' | NSTEp=',NSTEP,' DROP=',DROP,' TOLErance=',TOLER,
     &' | BMIN=',BFMIN,' BMAX=',BFMAX,
     &' QMIN=',QFMIN,' QMAX=',QFMAX,
     &' | FPMIN=',FPMIN,' FPMAX=',FDPMAX,
     &' FDPMIN=',FDPMIN,' FDPMAX=',FDPMAX,
     &' | number of B-factor and occupancy groups=',NGROUP
      WRITE(6,'(2A)') ' | METHod=',MINTYPE
      WRITE(6,'(2A)') ' ---------------------------',
     & '----------------------------------------------------'
      ELSE
      CALL CHKEND('GROUP>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
C check for overlapping definitions
      ERROR=.FALSE.
      DO IAT=1,NATOM
      FLAGS(IAT)=0
      END DO
      DO I=1,NGROUP
      IF (GTYPE(I).EQ.1) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      FLAGS(LIST(IAT))=FLAGS(LIST(IAT))+1
      END DO
      END IF
      END DO
      DO IAT=1,NATOM
      IF (FLAGS(IAT).GT.1) THEN
      WRITE(6,'(9A)')
     & ' %XGROUP-ERR: B-GROUP definition overlap for ( ',
     & SEGID(IAT),' ',RES(IAT),' ',RESID(IAT),' ',TYPE(IAT),' )'
      ERROR=.TRUE.
      END IF
      END DO
C
      DO IAT=1,NATOM
      FLAGS(IAT)=0
      END DO
      DO I=1,NGROUP
      IF (GTYPE(I).EQ.2) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      FLAGS(XRATOF(LIST(IAT)))=FLAGS(XRATOF(LIST(IAT)))+1
      END DO
      END IF
      END DO
      DO IAT=1,NATOM
      IF (FLAGS(IAT).GT.1) THEN
      WRITE(6,'(9A)')
     & ' %XGROUP-ERR: Q-GROUP definition overlap for ( ',
     & SEGID(IAT),' ',RES(IAT),' ',RESID(IAT),' ',TYPE(IAT),' )'
      ERROR=.TRUE.
      END IF
      END DO
C
      DO IAT=1,NATOM
      FLAGS(IAT)=0
      END DO
      DO I=1,NGROUP
      IF (GTYPE(I).EQ.3) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      FLAGS(XRATOF(LIST(IAT)))=FLAGS(XRATOF(LIST(IAT)))+1
      END DO
      END IF
      END DO
      DO IAT=1,NATOM
      IF (FLAGS(IAT).GT.1) THEN
      WRITE(6,'(9A)')
     & ' %XGROUP-ERR: FP definition overlap for ( ',
     & SEGID(IAT),' ',RES(IAT),' ',RESID(IAT),' ',TYPE(IAT),' )'
      ERROR=.TRUE.
      END IF
      END DO
C
      DO IAT=1,NATOM
      FLAGS(IAT)=0
      END DO
      DO I=1,NGROUP
      IF (GTYPE(I).EQ.4) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      FLAGS(XRATOF(LIST(IAT)))=FLAGS(XRATOF(LIST(IAT)))+1
      END DO
      END IF
      END DO
      DO IAT=1,NATOM
      IF (FLAGS(IAT).GT.1) THEN
      WRITE(6,'(9A)')
     & ' %XGROUP-ERR: FDP-GROUP definition overlap for ( ',
     & SEGID(IAT),' ',RES(IAT),' ',RESID(IAT),' ',TYPE(IAT),' )'
      ERROR=.TRUE.
      END IF
      END DO
C
C for FP refinement we have to check that the
C definition of groups is compatible with the definition
C of atomic scatterers
      DO IAT=1,NATOM
      FLAGS(IAT)=0
      END DO
      DO I=1,NGROUP
      IF (GTYPE(I).EQ.3) THEN
      SCATTX=XRINDF(LIST(POINT(I)+1))
      DO IAT=POINT(I)+1,POINT(I+1)
      IF (FLAGS(XRATOF(LIST(IAT))).NE.-9999) THEN
      FLAGS(XRATOF(LIST(IAT)))=SCATTX
      END IF
      IF (XRINDF(LIST(IAT)).NE.SCATTX) THEN
      FLAGS(XRATOF(LIST(IAT)))=-9999
      END IF
      END DO
      DO IAT=1,XRNATF
      IF (XRINDF(IAT).EQ.SCATTX.AND.
     &         FLAGS(XRATOF(IAT)).NE.SCATTX) THEN
      FLAGS(XRATOF(IAT))=-9999
      END IF
      END DO
      END IF
      END DO
      DO IAT=1,NATOM
      IF (FLAGS(IAT).LT.0) THEN
      WRITE(6,'(9A)')
     & ' %XGROUP-ERR: FP-GROUP incompatible with SCATter for ( ',
     & SEGID(IAT),' ',RES(IAT),' ',RESID(IAT),' ',TYPE(IAT),' )'
      ERROR=.TRUE.
      END IF
      END DO
C
C
C for FDP refinement we have to check that the
C definition of groups is compatible with the definition
C of atomic scatterers
      DO IAT=1,NATOM
      FLAGS(IAT)=0
      END DO
      DO I=1,NGROUP
      IF (GTYPE(I).EQ.4) THEN
      SCATTX=XRINDF(LIST(POINT(I)+1))
      DO IAT=POINT(I)+1,POINT(I+1)
      IF (FLAGS(XRATOF(LIST(IAT))).NE.-9999) THEN
      FLAGS(XRATOF(LIST(IAT)))=SCATTX
      END IF
      IF (XRINDF(LIST(IAT)).NE.SCATTX) THEN
      FLAGS(XRATOF(LIST(IAT)))=-9999
      END IF
      END DO
      DO IAT=1,XRNATF
      IF (XRINDF(IAT).EQ.SCATTX.AND.
     &         FLAGS(XRATOF(IAT)).NE.SCATTX) THEN
      FLAGS(XRATOF(IAT))=-9999
      END IF
      END DO
      END IF
      END DO
      DO IAT=1,NATOM
      IF (FLAGS(IAT).LT.0) THEN
      WRITE(6,'(9A)')
     & ' %XGROUP-ERR: FDP-GROUP incompatible with SCATter for ( ',
     & SEGID(IAT),' ',RES(IAT),' ',RESID(IAT),' ',TYPE(IAT),' )'
      ERROR=.TRUE.
      END IF
      END DO
C
      IF (.NOT.ERROR.AND.NGROUP.GT.0.AND.NATOM.GT.0) THEN
C
C initialize the lowest target value
      MTARG=R4BIG
      MCYCLE=0
C
C
C there are two sets of temporary parameter arrays.
C The OMAIN, QOLD, XPOLD, and XDPOLD sets are necessary
C to store the starting values of the various parameters.
C These arrays are necessary since the minimization routine
C (ZXCGR) only produces shifts relative to the starting
C values and the main arrays need to be updated at
C each minimization step in order to compute the target value.
C
C The MWMAIN, MQMAIN, MXRFP, MXRFDP arrays are needed to
C store the values that correspond to the cycle with the
C lowest target value. This is necessary since the ZXCGR routine
C doesn't necessarily exit with the lowest value.
C
C allocate space
      GROUP=ALLHP(IREAL8(NGROUP))
      GRADNT=ALLHP(IREAL8(NGROUP))
      IOMAIN=ALLHP(IREAL8(NATOM))
      IQOLD=ALLHP(IREAL8(NATOM))
      FPOLD=ALLHP(IREAL8(XRSN))
      FDPOLD=ALLHP(IREAL8(XRSN))
      FWMAIN=ALLHP(IREAL8(NATOM))
      FQMAIN=ALLHP(IREAL8(NATOM))
      FXRFP=ALLHP(IREAL8(XRSN))
      FXRFDP=ALLHP(IREAL8(XRSN))
      IF (DEBUG) THEN
      CALL XGDEBG(HEAP(GROUP),HEAP(GRADNT),HEAP(IOMAIN),
     &            HEAP(IQOLD),HEAP(FPOLD),HEAP(FDPOLD),DROP,NGROUP,
     &            XRFP,XRFDP,XRINDF,XRNATF)
      ELSE
      CALL XGRUP3(NSTEP,DROP,TOLER,NGROUP,POINT,LIST,
     &             HEAP(GROUP),HEAP(GRADNT),HEAP(IOMAIN),
     &             HEAP(IQOLD),HEAP(FPOLD),HEAP(FDPOLD),XRFP,XRFDP,
     &             XRINDF,XRNATF,
     &             HEAP(FWMAIN),HEAP(FQMAIN),HEAP(FXRFP),
     &             HEAP(FXRFDP),GTYPE,HEAP(HPATOF),MINTYPE)
      END IF
      CALL FREHP(FXRFDP,IREAL8(XRSN))
      CALL FREHP(FXRFP,IREAL8(XRSN))
      CALL FREHP(FQMAIN,IREAL8(NATOM))
      CALL FREHP(FWMAIN,IREAL8(NATOM))
      CALL FREHP(FPOLD,IREAL8(XRSN))
      CALL FREHP(FDPOLD,IREAL8(XRSN))
      CALL FREHP(IQOLD,IREAL8(NATOM))
      CALL FREHP(IOMAIN,IREAL8(NATOM))
      CALL FREHP(GRADNT,IREAL8(NGROUP))
      CALL FREHP(GROUP,IREAL8(NGROUP))
      END IF
C
C error label
9999  CONTINUE
      RETURN
      END
C======================================================================
      SUBROUTINE XGDEBG(GROUP,GRADNT,OMAIN,QOLD,
     &                  XPOLD,XDPOLD,STEP,NGROUP,
     &                  XRFP,XRFDP,XRINDF,XRNATF)
C
C check consistency of B-factor derivatives
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      DOUBLE PRECISION GROUP(*), GRADNT(*), OMAIN(*), QOLD(*)
      DOUBLE PRECISION XPOLD(*), XDPOLD(*), STEP
      INTEGER NGROUP
      DOUBLE PRECISION XRFP(*), XRFDP(*)
      INTEGER XRINDF(*), XRNATF
C local
      INTEGER J, I
      DOUBLE PRECISION TARGET, NUMGRD
C parameters
      DOUBLE PRECISION ZERO, TWO
      PARAMETER (ZERO=0.0D0, TWO=2.0D0)
C begin
      ICYCLE=-1
      DO J=1,NATOM
      OMAIN(J)=WMAIN(J)
      QOLD(J)=QMAIN(J)
      END DO
C
      DO I=1,XRNATF
      XPOLD(XRINDF(I))=XRFP(XRINDF(I))
      XDPOLD(XRINDF(I))=XRFDP(XRINDF(I))
      END DO
C
      DO J=1,NGROUP
      GROUP(J)=ZERO
      END DO
C
C loop over all groups
      DO J=1,NGROUP
      GROUP(J)=STEP
      CALL XGFUNC(NGROUP,GROUP,TARGET,GRADNT)
      NUMGRD=TARGET
      GROUP(J)=-STEP
      CALL XGFUNC(NGROUP,GROUP,TARGET,GRADNT)
      NUMGRD=NUMGRD-TARGET
      NUMGRD=NUMGRD/(TWO*STEP)
      GROUP(J)=ZERO
      CALL XGFUNC(NGROUP,GROUP,TARGET,GRADNT)
      WRITE(6,'(A,I6,A,F12.5,A,F12.5,A,F12.5)')
     &' XGDEBUG: group # ',J,' analyt.=',GRADNT(J),
     &' finite=',NUMGRD
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XGRUP3(NSTEP,DROP,TOLER,NGROUP,POINT,LIST,
     &                   GROUP,GRADNT,OMAIN,QOLD,XPOLD,XDPOLD,
     &                   XRFP,XRFDP,XRINDF,XRNATF,
     &                   MWMAIN,MQMAIN,MXRFP,MXRFDP,GTYPE,XRATOF,
     &                   MINTYPE)
C
C Routine for group b-factor optimization.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER NSTEP
      DOUBLE PRECISION DROP, TOLER
      INTEGER NGROUP, POINT(*), LIST(*)
      DOUBLE PRECISION GROUP(*), GRADNT(*), OMAIN(*), QOLD(*)
      DOUBLE PRECISION XPOLD(*), XDPOLD(*)
      DOUBLE PRECISION XRFP(*), XRFDP(*)
      INTEGER XRINDF(*), XRNATF
      DOUBLE PRECISION MWMAIN(*), MQMAIN(*)
      DOUBLE PRECISION MXRFP(*), MXRFDP(*)
      INTEGER GTYPE(*)
      INTEGER XRATOF(*)
      CHARACTER*(*) MINTYPE
C local
      INTEGER I, IER, IAT
      DOUBLE PRECISION TARGET
      INTEGER M, NPRINT, IW
C pointer
      INTEGER WORK, DIAG, BEST
C external
      EXTERNAL XGFUNC, XGFUNC2
C parameters
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      TOLER=(TOLER**2)*6
C
C initialize the group b factor shifts
C note: the inidividual b factor will be the sum of the group b
C factor shifts and the previously assigned b individual b factor
C
      DO I=1,NGROUP
      GROUP(I)=ZERO
      END DO
      DO I=1,NATOM
      OMAIN(I)=WMAIN(I)
      QOLD(I)=QMAIN(I)
      END DO
C
      DO I=1,XRNATF
      XPOLD(XRINDF(I))=XRFP(XRINDF(I))
      XDPOLD(XRINDF(I))=XRFDP(XRINDF(I))
      END DO
C
      ICYCLE=0
C
C
C store current parameters
      DO I=1,NGROUP
      IF (GTYPE(I).EQ.1) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      MWMAIN(XRATOF(LIST(IAT)))=WMAIN(XRATOF(LIST(IAT)))
      END DO
      ELSEIF (GTYPE(I).EQ.2) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      MQMAIN(XRATOF(LIST(IAT)))=QMAIN(XRATOF(LIST(IAT)))
      END DO
      ELSEIF (GTYPE(I).EQ.3) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      MXRFP(XRINDF(LIST(IAT)))=XRFP(XRINDF(LIST(IAT)))
      END DO
      ELSEIF (GTYPE(I).EQ.4) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      MXRFDP(XRINDF(LIST(IAT)))=XRFDP(XRINDF(LIST(IAT)))
      END DO
      END IF
      END DO
C
      IF (MINTYPE.EQ.'POWELL') THEN
C
C Powell-method minimizer
      WORK=ALLHP(IREAL8(NGROUP*6))
      CALL ZXCGR(XGFUNC,NGROUP,TOLER,NSTEP,DROP,GROUP,GRADNT,
     &           TARGET,HEAP(WORK),IER)
      CALL FREHP(WORK,IREAL8(NGROUP*6))
      IF (IER.EQ.0) THEN
      WRITE(6,'(A)') ' ZXCGR: gradient converged'
      ELSE IF (IER.EQ.129) THEN
      WRITE(6,'(A)') ' ZXCGR: Line search terminated'
      ELSE IF (IER.EQ.130) THEN
      WRITE(6,'(A)') ' %ZXCGR-ERR: Search direction uphill'
      ELSE IF (IER.EQ.131) THEN
      WRITE(6,'(A)') ' ZXCGR: NSTEP limit reached'
      ELSE IF (IER.EQ.132) THEN
      WRITE(6,'(A)') ' %ZXCGR-ERR: failure to reduce E'
      END IF
      ELSE IF (MINTYPE.EQ.'LBFGS') THEN
C
C LBFGS-method minimizer
      M=6
      NPRINT=1
      IW=NGROUP*(2*M+1)+2*M
      DIAG=ALLHP(IREAL8(NGROUP))
      BEST=ALLHP(IREAL8(NGROUP))
      WORK=ALLHP(IREAL8(IW))
      CALL ZXLBFGS(XGFUNC2,NGROUP,M,GROUP,GRADNT,TARGET,
     &             .FALSE.,HEAP(DIAG),NPRINT,TOLER,NSTEP,HEAP(WORK),
     &             HEAP(BEST),10)
      CALL FREHP(WORK,IREAL8(IW))
      CALL FREHP(BEST,IREAL8(NGROUP))
      CALL FREHP(DIAG,IREAL8(NGROUP))
      ELSE
      CALL WRNDIE(-5,'XGROUP','unknown minimization method.')
      END IF
C
C set parameters to the values that produces the lowest target value
      DO I=1,NGROUP
      IF (GTYPE(I).EQ.1) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      WMAIN(XRATOF(LIST(IAT)))=MWMAIN(XRATOF(LIST(IAT)))
      END DO
      ELSEIF (GTYPE(I).EQ.2) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      QMAIN(XRATOF(LIST(IAT)))=MQMAIN(XRATOF(LIST(IAT)))
      END DO
      ELSEIF (GTYPE(I).EQ.3) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      XRFP(XRINDF(LIST(IAT)))=MXRFP(XRINDF(LIST(IAT)))
      END DO
      ELSEIF (GTYPE(I).EQ.4) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      XRFDP(XRINDF(LIST(IAT)))=MXRFDP(XRINDF(LIST(IAT)))
      END DO
      END IF
      END DO
C
      WRITE(6,'(A,I6,A)')
     & ' XGROUP: using the parameters (at cycle ', MCYCLE,
     & ') which produce the minimum target. '
C
      RETURN
      END
C======================================================================
      SUBROUTINE XGFUNC(NGROUP,GROUP,TARGET,GRADNT)
C
C Target function for powell method minimizer for group
C b factor optimization.  Front-end routine which passes the
C pionters of some arrays to the next routine
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INTEGER NGROUP
      DOUBLE PRECISION GROUP(*), TARGET, GRADNT(*)
C begin
      CALL XGFUN2(NGROUP,GROUP,TARGET,GRADNT,HEAP(ILIST),HEAP(IPOINT),
     &            HEAP(ITYPE),HEAP(IOMAIN),HEAP(IQOLD),HEAP(FPOLD),
     &            HEAP(FDPOLD),
     &            HEAP(HPATOF),HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),
     &            HEAP(HPDT),HEAP(HPDQ),HEAP(HPINDF),
     &            HEAP(FWMAIN),HEAP(FQMAIN),HEAP(FXRFP),
     &            HEAP(FXRFDP))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XGFUNC2(NGROUP,GROUP,TARGET,GRADNT,ITER,QPRINT)
C
C Target function for lbfgs method minimizer for group
C b factor optimization.  Front-end routine which passes the
C pionters of some arrays to the next routine
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'xropti.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INTEGER NGROUP
      DOUBLE PRECISION GROUP(*), TARGET, GRADNT(*)
      INTEGER ITER
      LOGICAL QPRINT
C begin
      CALL XGFUN2(NGROUP,GROUP,TARGET,GRADNT,HEAP(ILIST),HEAP(IPOINT),
     &            HEAP(ITYPE),HEAP(IOMAIN),HEAP(IQOLD),HEAP(FPOLD),
     &            HEAP(FDPOLD),
     &            HEAP(HPATOF),HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),
     &            HEAP(HPDT),HEAP(HPDQ),HEAP(HPINDF),
     &            HEAP(FWMAIN),HEAP(FQMAIN),HEAP(FXRFP),
     &            HEAP(FXRFDP))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XGFUN2(NGROUP,GROUP,TARGET,GRADNT,LIST,POINT,GTYPE,
     &           OMAIN,QOLD,XPOLD,XDPOLD,
     &           XRATOF,XRDX,XRDY,XRDZ,XRDT,XRDQ,
     &           XRINDF,
     &           MWMAIN,MQMAIN,MXRFP,MXRFDP)
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'xropti.inc'
      INTEGER NGROUP
      DOUBLE PRECISION GROUP(*), TARGET, GRADNT(*)
      INTEGER LIST(*), POINT(*), GTYPE(*)
      DOUBLE PRECISION OMAIN(*), QOLD(*), XPOLD(*), XDPOLD(*)
      INTEGER XRATOF(*)
      DOUBLE PRECISION XRDX(*), XRDY(*), XRDZ(*), XRDT(*), XRDQ(*)
      INTEGER XRINDF(*)
      DOUBLE PRECISION MWMAIN(*), MQMAIN(*)
      DOUBLE PRECISION MXRFP(*), MXRFDP(*)
C local
      INTEGER I, IAT
      DOUBLE PRECISION GRAD
C parameters
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
C Initialize the gradient vector
      DO I=1,NGROUP
      GRADNT(I)=ZERO
      END DO
C
C truncate B-factors to BFMIN...BFMAX and copy into WMAIN
      DO I=1,NGROUP
      IF (GTYPE(I).EQ.1) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      WMAIN(XRATOF(LIST(IAT)))=
     &   MIN(MAX(OMAIN(XRATOF(LIST(IAT)))+GROUP(I),BFMIN),BFMAX)
      END DO
      ELSEIF (GTYPE(I).EQ.2) THEN
C truncate Qs to QMIN...QMAX and copy into QMAIN
      DO IAT=POINT(I)+1,POINT(I+1)
      QMAIN(XRATOF(LIST(IAT)))=
     &    MIN(MAX(QOLD(XRATOF(LIST(IAT)))+GROUP(I),QFMIN),QFMAX)
      END DO
      ELSEIF (GTYPE(I).EQ.3) THEN
C truncate FPs to FPMIN...FPMAX and copy into XRFP
      DO IAT=POINT(I)+1,POINT(I+1)
      XRFP(XRINDF(LIST(IAT)))=
     &    MIN(MAX(XPOLD(XRINDF(LIST(IAT)))+GROUP(I),FPMIN),
     &    FPMAX)
      END DO
      ELSEIF (GTYPE(I).EQ.4) THEN
C truncate FDPs to FDPMIN...FDPMAX and copy into XRFDP
      DO IAT=POINT(I)+1,POINT(I+1)
      XRFDP(XRINDF(LIST(IAT)))=
     &    MIN(MAX(XDPOLD(XRINDF(LIST(IAT)))+GROUP(I),FDPMIN),
     &    FDPMAX)
      END DO
C
      END IF
      END DO
C
C compute the diffraction energy term and gradient
C ================================================
      CALL XCALCS(.FALSE.,.TRUE.,XRE,
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
C store coordinates, b-factors, and occupancies for
C lowest target value
      IF (XRE.LT.MTARG) THEN
C
C truncate B-factors to BFMIN...BFMAX and copy into WMAIN
      DO I=1,NGROUP
      IF (GTYPE(I).EQ.1) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      MWMAIN(XRATOF(LIST(IAT)))=WMAIN(XRATOF(LIST(IAT)))
      END DO
      ELSEIF (GTYPE(I).EQ.2) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      MQMAIN(XRATOF(LIST(IAT)))=QMAIN(XRATOF(LIST(IAT)))
      END DO
      ELSEIF (GTYPE(I).EQ.3) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      MXRFP(XRINDF(LIST(IAT)))=XRFP(XRINDF(LIST(IAT)))
      END DO
      ELSEIF (GTYPE(I).EQ.4) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      MXRFDP(XRINDF(LIST(IAT)))=XRFDP(XRINDF(LIST(IAT)))
      END DO
      END IF
      END DO
      END IF
C
      TARGET=XRE
C
C fill the gradient vector
      GRAD=ZERO
      DO I=1,NGROUP
      GRADNT(I)=ZERO
C
      IF (GTYPE(I).EQ.1) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      IF (WMAIN(XRATOF(LIST(IAT))).GT.BFMIN.AND.
     &    WMAIN(XRATOF(LIST(IAT))).LT.BFMAX) THEN
      GRADNT(I)=GRADNT(I)+XRDT(XRATOF(LIST(IAT)))
C modification, ATB, 12/26/08
CCC      ELSE IF (WMAIN(XRATOF(LIST(IAT))).LE.BFMIN
CCC     &         .AND.XRDT(XRATOF(LIST(IAT))).LT.ZERO) THEN
CCC      GRADNT(I)=GRADNT(I)+XRDT(XRATOF(LIST(IAT)))
CCC      ELSE IF (WMAIN(XRATOF(LIST(IAT))).GE.BFMAX
CCC     &         .AND.XRDT(XRATOF(LIST(IAT))).GT.ZERO) THEN
CCC      GRADNT(I)=GRADNT(I)+XRDT(XRATOF(LIST(IAT)))
      END IF
      END DO
C
      ELSEIF (GTYPE(I).EQ.2) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      IF (QMAIN(XRATOF(LIST(IAT))).GT.QFMIN.AND.
     &    QMAIN(XRATOF(LIST(IAT))).LT.QFMAX) THEN
      GRADNT(I)=GRADNT(I)+XRDQ(XRATOF(LIST(IAT)))
C modification, ATB, 12/26/08
CCC      ELSE IF (QMAIN(XRATOF(LIST(IAT))).LE.QFMIN
CCC   &         .AND.XRDQ(XRATOF(LIST(IAT))).LT.ZERO) THEN
CCC      GRADNT(I)=GRADNT(I)+XRDQ(XRATOF(LIST(IAT)))
CCC      ELSE IF (QMAIN(XRATOF(LIST(IAT))).GE.QFMAX
CCC     &         .AND.XRDQ(XRATOF(LIST(IAT))).GT.ZERO) THEN
CCC      GRADNT(I)=GRADNT(I)+XRDQ(XRATOF(LIST(IAT)))
      END IF
      END DO
C
      ELSEIF (GTYPE(I).EQ.3) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      IF (XRFP(XRINDF(LIST(IAT))).GT.FPMIN.AND.
     &    XRFP(XRINDF(LIST(IAT))).LT.FPMAX) THEN
      GRADNT(I)=GRADNT(I)+XRDX(XRATOF(LIST(IAT)))
C modification, ATB, 12/26/08
CCC      ELSE IF (XRFP(XRINDF(LIST(IAT))).LE.FPMIN
CCC     &         .AND.XRDX(XRATOF(LIST(IAT))).LT.ZERO) THEN
CCC      GRADNT(I)=GRADNT(I)+XRDX(XRATOF(LIST(IAT)))
CCC      ELSE IF (XRFP(XRINDF(LIST(IAT))).GE.FPMAX
CCC     &         .AND.XRDX(XRATOF(LIST(IAT))).GT.ZERO) THEN
CCC      GRADNT(I)=GRADNT(I)+XRDX(XRATOF(LIST(IAT)))
      END IF
      END DO
C
      ELSEIF (GTYPE(I).EQ.4) THEN
      DO IAT=POINT(I)+1,POINT(I+1)
      IF (XRFDP(XRINDF(LIST(IAT))).GT.FDPMIN.AND.
     &    XRFDP(XRINDF(LIST(IAT))).LT.FDPMAX) THEN
      GRADNT(I)=GRADNT(I)+XRDY(XRATOF(LIST(IAT)))
C modification, ATB, 12/26/08
CCC      ELSE IF (XRFDP(XRINDF(LIST(IAT))).LE.FDPMIN
CCC     &         .AND.XRDY(XRATOF(LIST(IAT))).LT.ZERO) THEN
CCC      GRADNT(I)=GRADNT(I)+XRDY(XRATOF(LIST(IAT)))
CCC      ELSE IF (XRFDP(XRINDF(LIST(IAT))).GE.FDPMAX
CCC     &         .AND.XRDY(XRATOF(LIST(IAT))).GT.ZERO) THEN
CCC      GRADNT(I)=GRADNT(I)+XRDY(XRATOF(LIST(IAT)))
      END IF
      END DO
C
      END IF
      GRAD=GRAD+GRADNT(I)**2
      END DO
C
C print some stuff
      IF (ICYCLE.GE.0) THEN
      GRAD=SQRT(GRAD/NGROUP)
      ICYCLE=ICYCLE+1
      WRITE(6,'(A,I6,A,F10.4,A)')
     & ' --------------- cycle=',ICYCLE,
     & ' --------------------------------------------------'
      WRITE(6,'(A,E10.3,A,E10.3,A)')
     & ' | Target=E(XREF)=',TARGET,'  grad(E)=',GRAD,
     & '                               |'
      WRITE(6,'(2A)') ' --------------------------------------------',
     & '-----------------------------------'
      END IF
C
      IF (XRE.LT.MTARG) THEN
      MTARG=XRE
      MCYCLE=ICYCLE
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE BNCS(RENER,TEMP,NATOM,WMAIN,IMOVE,FLAG,XRATOF,
     &          XRNATF,NATNCS,NUMEQV,POINTR,SIGBNC,NNCSX)
C
C routine computes NCS B-factor restraints for a group
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      DOUBLE PRECISION RENER, TEMP(*)
      INTEGER NATOM
      DOUBLE PRECISION WMAIN(*)
      INTEGER IMOVE(*), FLAG(*), XRATOF(*), XRNATF
      INTEGER NATNCS, NUMEQV, POINTR(*)
      DOUBLE PRECISION SIGBNC
      INTEGER NNCSX
C local
      INTEGER BBAR
C begin
      BBAR=ALLHP(IREAL8(NATOM))
      CALL BNCS2(RENER,TEMP,NATOM,WMAIN,IMOVE,FLAG,XRATOF,
     &     XRNATF,NATNCS,NUMEQV,POINTR,SIGBNC,NNCSX,HEAP(BBAR))
      CALL FREHP(BBAR,IREAL8(NATOM))
      RETURN
      END
C=====================================================================
      SUBROUTINE BNCS2(RENER,TEMP,NATOM,WMAIN,IMOVE,FLAG,XRATOF,
     &          XRNATF,NATNCS,NUMEQV,POINTR,SIGBNC,NNCSX,BBAR)
C
C routine computes NCS B-factor restraints for a group
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION RENER, TEMP(*)
      INTEGER NATOM
      DOUBLE PRECISION WMAIN(*)
      INTEGER IMOVE(*), FLAG(*), XRATOF(*), XRNATF
      INTEGER NATNCS, NUMEQV, POINTR(*)
      DOUBLE PRECISION SIGBNC
      INTEGER NNCSX
      DOUBLE PRECISION BBAR(*)
C local
      INTEGER I, J, EQPOS, IPT
      DOUBLE PRECISION DELB, DB
C parameter
      DOUBLE PRECISION ZERO, TWO
      PARAMETER (ZERO=0.0D0, TWO=2.0D0)
C
C define pointer array
C fill atom flags
      DO I=1,NATOM
      FLAG(I)=0
      END DO
      DO I=1,XRNATF
      IF (IMOVE(XRATOF(I)).EQ.0) FLAG(XRATOF(I))=I
      END DO
C
      DO I=1,NATNCS
      BBAR(I)=ZERO
      END DO
C
      DO J=1,NUMEQV
      EQPOS=NATNCS*(J-1)
      DO I=1,NATNCS
      BBAR(I)=BBAR(I) + WMAIN(POINTR(I+EQPOS))
      END DO
      END DO
C
      DO I=1,NATNCS
      BBAR(I)=BBAR(I)/NUMEQV
      END DO
C
C Compute energy and derivatives.
      NNCSX=0
      DO J=1,NUMEQV
      EQPOS=NATNCS*(J-1)
      DO I=1,NATNCS
      IPT=FLAG(POINTR(I+EQPOS))
      IF (IPT.NE.0) THEN
      NNCSX=NNCSX+1
      DELB=WMAIN(POINTR(I+EQPOS))-BBAR(I)
      RENER=RENER + DELB**2/SIGBNC**2
      DB=TWO*DELB/SIGBNC**2
      TEMP(IPT)=TEMP(IPT) + DB
      END IF
      END DO
      END DO
      RETURN
      END
