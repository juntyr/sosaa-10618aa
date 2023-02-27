!=======================================================================
!  LD_FCT.EXT
!  This include file contains "light dependent" factors.
!
!  This is an input file for MEGAN v2.1 as according to Table 4 in Guenther et
!  al., Geosci. Model Dev., 5, 1471-1492, 2012
!=======================================================================

Module M2_LD_FCT

      INTEGER, parameter ::   N_LDF_SPC=22
      CHARACTER*16       ::   LDF_SPC(22)
      REAL               ::   LDF_FCT(22)
      INTEGER            ::   LDF_MAP(22)

      DATA     LDF_SPC(  1)      , LDF_FCT(  1), LDF_MAP(  1) &
            / 'ISOP            ', 0.9999      , 1            /

      DATA     LDF_SPC(  2)      , LDF_FCT(  2), LDF_MAP(  2) &
            / 'MYRC            ', 0.1        , 2            /

      DATA     LDF_SPC(  3)      , LDF_FCT(  3), LDF_MAP(  3) &
            / 'SABI            ', 0.8         , 3            /

      DATA     LDF_SPC(  4)      , LDF_FCT(  4), LDF_MAP(  4) &
            / 'LIMO            ', 0.1        , 4            /

      DATA     LDF_SPC(  5)      , LDF_FCT(  5), LDF_MAP(  5) &
            / '3CAR            ', 0.1        , 5            /

      DATA     LDF_SPC(  6)      , LDF_FCT(  6), LDF_MAP(  6) &
            / 'OCIM            ', 0.8         , 6            /

      DATA     LDF_SPC(  7)      , LDF_FCT(  7), LDF_MAP(  7) &
            / 'BPIN            ', 0.1        , 7            /

      DATA     LDF_SPC(  8)      , LDF_FCT(  8), LDF_MAP(  8) &
            / 'APIN            ', 0.1        , 8            /

      DATA     LDF_SPC(  9)      , LDF_FCT(  9), LDF_MAP(  9) &
            / 'OMTP            ', 0.1        , 9            /

      DATA     LDF_SPC( 10)      , LDF_FCT( 10), LDF_MAP( 10) &
            / 'FARN            ', 0.1         , 10           /

      DATA     LDF_SPC( 11)      , LDF_FCT( 11), LDF_MAP( 11) &
            / 'BCAR            ', 0.1         , 11           /

      DATA     LDF_SPC( 12)      , LDF_FCT( 12), LDF_MAP( 12) &
            / 'OSQT            ', 0.1         , 12           /

      DATA     LDF_SPC( 13)      , LDF_FCT( 13), LDF_MAP( 13) &
            / 'MBO             ', 0.9999      , 13           /

      DATA     LDF_SPC( 14)      , LDF_FCT( 14), LDF_MAP( 14) &
            / 'MEOH            ', 0.75        , 14           /

      DATA     LDF_SPC( 15)      , LDF_FCT( 15), LDF_MAP( 15) &
            / 'ACTO            ', 0.25        , 15           /

      DATA     LDF_SPC( 16)      , LDF_FCT( 16), LDF_MAP( 16) &
            / 'CH4             ', 0.75        , 16           /

      DATA     LDF_SPC( 17)      , LDF_FCT( 17), LDF_MAP( 17) &
            / 'NO              ', 0.0         , 17           /

      DATA     LDF_SPC( 18)      , LDF_FCT( 18), LDF_MAP( 18) &
            / 'ACTA            ', 0.5         , 18           /

      DATA     LDF_SPC( 19)      , LDF_FCT( 19), LDF_MAP( 19) &
            / 'FORM            ', 0.5         , 19           /

      DATA     LDF_SPC( 20)      , LDF_FCT( 20), LDF_MAP( 20) &
             / 'CO              ', 1.0         , 20           /

      DATA     LDF_SPC( 21)      , LDF_FCT( 21), LDF_MAP( 21) &
             / 'CINE            ', 0.1        , 21           /

      DATA     LDF_SPC( 22)      , LDF_FCT( 22), LDF_MAP( 22) &
             / 'LINA            ', 0.1         , 22           /

end module M2_LD_FCT
