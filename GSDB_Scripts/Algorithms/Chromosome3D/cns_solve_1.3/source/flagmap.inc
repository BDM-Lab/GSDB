C Storage for local symmetry
      INTEGER   MSYMMX
      PARAMETER(MSYMMX = 192)
      INTEGER  SYMMX(MSYMMX, 3, 4), ITSYMMX(MSYMMX, 3, 3), NSYMMX
      COMMON /FMSYMM/ SYMMX, ITSYMMX, NSYMMX
C
C Storage for additional generators
      INTEGER   MADDLGEN
      PARAMETER(MADDLGEN = 6)
      INTEGER  ADDLGEN(MADDLGEN, 3, 4), ITADDLGEN(MADDLGEN, 3, 3)
      INTEGER  N1ADDLGEN, N2ADDLGEN
      COMMON /FMADDL/ ADDLGEN, ITADDLGEN, N1ADDLGEN, N2ADDLGEN
C
C Structure-seminvariant vectors & moduli
      INTEGER  N_SSVM, SSV(3, 3), SSM(3)
C Dimensions of current FlagMap
      INTEGER  RNFM(3)
C Number of independent points in FlagMap
      INTEGER  NINDFM
C Size of FlagMap
      INTEGER  NFLAGMAP
C Pointer for FlagMap
      INTEGER  PTFLAGMAP
C Dimensions of corresponding asymmetric unit
      INTEGER  NASYFM(2, 3)
C Parameters and internal flag ValidFlagMap
      LOGICAL  USESYM, USE_SS, USEADDL, VALIDFLAGMAP
C
      COMMON /FMINT/ N_SSVM, SSV, SSM,
     &               RNFM,
     &               NINDFM, NFLAGMAP,
     &               NASYFM
      COMMON /FMPTR/ PTFLAGMAP
      COMMON /FMLOG/ USESYM, USE_SS, USEADDL, VALIDFLAGMAP
C
      SAVE /FMSYMM/
      SAVE /FMADDL/
      SAVE /FMINT/
      SAVE /FMPTR/
      SAVE /FMLOG/
