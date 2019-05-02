C
C xmask.inc
C
      INTEGER MXMAVG
      PARAMETER (MXMAVG=200)
C
C constants for packing and unpacking masks
C storage is: LTX LTY LTZ SYMMOP
C
      INTEGER PMX,PMY,PMZ,MXTRAN
      PARAMETER (PMX=100000,PMY=10000)
      PARAMETER (PMZ=1000,MXTRAN=4)
C
