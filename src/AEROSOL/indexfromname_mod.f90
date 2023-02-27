MODULE INDEXFROMNAME_MOD
IMPLICIT NONE

contains
    !==============================================================================
    ! Function finds out the internal index using the proper name. E.g. IndexFromName('APINENE')
    ! will return 150 in current version.
    ! Input:
    ! name: the name (trimmed or untrimmed) of the variable. This name must match the name in NAMES.DAT exactly.
    ! Output:
    ! integer
    !..............................................................................
    PURE INTEGER FUNCTION IndexFromName(NAME, list_of_names)
      IMPLICIT NONE
      character(*), INTENT(IN) :: NAME
      character(*), optional, INTENT(IN) :: list_of_names(:)
      integer :: i,m
      if (PRESENT(list_of_names)) then
        m = size(list_of_names)
        DO i=1, m
          if (TRIM(NAME) == TRIM(list_of_names(I))) EXIT
        END DO
      END IF

      if (i == m+1) i=0
      IndexFromName = i

    END FUNCTION IndexFromName


END MODULE INDEXFROMNAME_MOD
