C seed.inc
C ========
C
C random number generator seed and parameters for random number
C generator.  The parameters always have to be computed as double
C precision variables.
      D O U B L E  P R E C I S I O N SEED, D2P31M, D2P31
C D2P31M=(2**31) - 1
C D2P31 =(2**31)(OR AN ADJUSTED VALUE)
      PARAMETER (D2P31M=2147483647.0D0)
      PARAMETER (D2P31=2147483711.0D0)
C
      COMMON /CCSEED/ SEED
      SAVE /CCSEED/
