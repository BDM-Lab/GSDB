C data structure MMDG
C =====================
C
C Common Block distance geometry module
C
C heap pointers for bounds matrix
      INTEGER DUBPTR
C lengths for bounds matrix
      INTEGER DUBLEN
C number of selected atoms for dg energy term
      INTEGER NDGSELE
C pointer for selection array for dg energy term
      INTEGER DGSELE
C length for selection array for dg energy term
      INTEGER DGSLEN
C scale factor for dg energy term restraints
      DOUBLE PRECISION DGSCALE
C exponent for dg energy term
      INTEGER DGEXPO
C flag indicating that vdw bounds have been added to bounds matrix
      INTEGER DGQVDW
C
C
C     integer block
C
      COMMON/MMDGC/ DUBPTR, DUBLEN, NDGSELE, DGSELE, DGSLEN, DGEXPO,
     &              DGQVDW
C
C double precision block
C
      COMMON/MMDGD/ DGSCALE
      SAVE /MMDGC/
      SAVE /MMDGD/
