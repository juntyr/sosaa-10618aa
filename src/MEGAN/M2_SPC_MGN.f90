!=======================================================================
!  SPC_MGN.EXT
!  This include file contains MEGAN species
!
!
!
!  MEGAN v2.0
!
!  Created by Tan 12/02/06
!=======================================================================

MODULE M2_SPC_MGN

      INTEGER, parameter   ::     N_MGN_SPC=22
      CHARACTER*16 ::  MGN_SPC(22)
      REAL        ::   MGN_MWT(22)

      DATA     MGN_SPC(  1), MGN_MWT(  1) / 'ISOP            ', 1.0    /

      DATA     MGN_SPC(  2), MGN_MWT(  2) / 'MYRC            ', 1.0    /
      DATA     MGN_SPC(  3), MGN_MWT(  3) / 'SABI            ', 1.0    /
      DATA     MGN_SPC(  4), MGN_MWT(  4) / 'LIMO            ', 1.0    /
      DATA     MGN_SPC(  5), MGN_MWT(  5) / '3CAR            ', 1.0    /
      DATA     MGN_SPC(  6), MGN_MWT(  6) / 'OCIM            ', 1.0    /
      DATA     MGN_SPC(  7), MGN_MWT(  7) / 'BPIN            ', 1.0    /
      DATA     MGN_SPC(  8), MGN_MWT(  8) / 'APIN            ', 1.0    /

      DATA     MGN_SPC(  9), MGN_MWT(  9) / 'OMTP            ', 1.0    /

      DATA     MGN_SPC( 10), MGN_MWT( 10) / 'FARN            ', 1.0    /
      DATA     MGN_SPC( 11), MGN_MWT( 11) / 'BCAR            ', 1.0    /

      DATA     MGN_SPC( 12), MGN_MWT( 12) / 'OSQT            ', 1.0    /

      DATA     MGN_SPC( 13), MGN_MWT( 13) / 'MBO             ', 1.0    /
      DATA     MGN_SPC( 14), MGN_MWT( 14) / 'MEOH            ', 1.0    /
      DATA     MGN_SPC( 15), MGN_MWT( 15) / 'ACTO            ', 1.0    /
      DATA     MGN_SPC( 16), MGN_MWT( 16) / 'CH4             ', 1.0    /
      DATA     MGN_SPC( 17), MGN_MWT( 17) / 'NO              ', 1.0    /

      DATA     MGN_SPC( 18), MGN_MWT( 18) / 'ACTA            ', 1.0    /

      DATA     MGN_SPC( 19), MGN_MWT( 19) / 'FORM            ', 1.0    /

      DATA     MGN_SPC( 20), MGN_MWT( 20) / 'CO              ', 1.0    /

      DATA     MGN_SPC( 21), MGN_MWT( 21) / 'CINE            ', 1.0    /

      DATA     MGN_SPC( 22), MGN_MWT( 22) / 'LINA            ', 1.0    /
      
END MODULE M2_SPC_MGN
