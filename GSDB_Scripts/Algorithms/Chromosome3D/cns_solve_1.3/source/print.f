C======================================================================
      SUBROUTINE PRINT
C
C parse PRINT qualifiers
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'noe.inc'
C local
      CHARACTER*4 SOBJ
      DOUBLE PRECISION THRESH
C begin
      CALL NEXTA4('PRINT-object=',SOBJ)
      IF (SOBJ.EQ.'THRE') THEN
      CALL NEXTF('THREshold=',THRESH)
      CALL NEXTA4('PRINT-object=',SOBJ)
      ELSE
      THRESH=-0.0001D0
      END IF
C
      IF (SOBJ.EQ.'HELP') THEN
C
      CALL CNSHELP('cns-print')
C
      ELSE IF (SOBJ.EQ.'BOND') THEN
      CALL PRPICK('BOND',THRESH)
      ELSE IF (SOBJ.EQ.'ANGL') THEN
      CALL PRPICK('ANGL',THRESH)
      ELSE IF (SOBJ.EQ.'DIHE') THEN
      CALL PRPICK('DIHE',THRESH)
      ELSE IF (SOBJ.EQ.'IMPR') THEN
      CALL PRPICK('IMPR',THRESH)
      ELSE IF (SOBJ.EQ.'NOE ') THEN
      CALL NOEPRI(THRESH,0)
      IF (NOEICV.GT.0) CALL NOEPRI(THRESH,1)
      ELSE IF (SOBJ.EQ.'CDIH') THEN
      CALL CDIPRI(PUNIT,CNPHI,CNSCA,HEAP(ICIP),HEAP(ICJP),
     &        HEAP(ICKP),HEAP(ICLP),HEAP(ICICP),HEAP(ICCPD),
     &        HEAP(ICCPE),HEAP(ICCPB),HEAP(ICCPC),HEAP(ICCPO),
     &        THRESH,HEAP(HPDCV),0,DIHICV)
      IF (DIHICV.GT.0) THEN
      CALL CDIPRI(PUNIT,CNPHI,CNSCA,HEAP(ICIP),HEAP(ICJP),
     &            HEAP(ICKP),HEAP(ICLP),HEAP(ICICP),HEAP(ICCPD),
     &            HEAP(ICCPE),HEAP(ICCPB),HEAP(ICCPC),HEAP(ICCPO),
     &            THRESH,HEAP(HPDCV),1,DIHICV)
      END IF
      ELSE
      CALL DSPERR('PRINT','unknown value')
      END IF
      RETURN
      END
C======================================================================
