!=======================================================================
!  REL_EM_ACT.EXT
!
!  This file contains activity factors that describe the relative emission 
!  rates for new, growing, mature, and senescing leaves, respectively for
!  individual compounds.
!
!  MEGAN v2.1
!=======================================================================

MODULE M2_REL_EM_ACT

      INTEGER, parameter ::      N_CAT=5
      REAL    ::       Anew(5)
      REAL    ::       Agro(5)
      REAL    ::       Amat(5)
      REAL    ::       Aold(5)

      DATA    Anew(  1),  Agro(  1),  Amat(  1),  Aold(  1) &
          /  1.0      ,  1.0      ,  1.0      ,  1.0       /

      DATA    Anew(  2),  Agro(  2),  Amat(  2),  Aold(  2) &
          /  2.0      ,  1.8      ,  1.0      ,  1.05       /

      DATA    Anew(  3),  Agro(  3),  Amat(  3),  Aold(  3) &
          /  0.4      ,  0.6      ,  1.0      ,  0.95      /

      DATA    Anew(  4),  Agro(  4),  Amat(  4),  Aold(  4) &
          /  3.5      ,  3.0      ,  1.0      ,  1.2       /

      DATA    Anew(  5),  Agro(  5),  Amat(  5),  Aold(  5) &
          /  0.05     ,  0.6      ,  1.0      ,  0.9       /
     
END MODULE M2_REL_EM_ACT
