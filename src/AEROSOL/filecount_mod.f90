MODULE FILECOUNT_MOD
IMPLICIT NONE

contains
!==============================================================================
! Function to count number of rows in a file. Takes in the UNIT (file pointer)
! and returns number of rows in a file (integer). The function can be used at
! any point of file reading, the UNIT will be moved back where it was in the
! beginning.
! NOTE:
! The function does not check if the file has equal amounts of colums in each row
! Empty rows are also counted
!..............................................................................
INTEGER FUNCTION ROWCOUNT(file_id,hdr)
  IMPLICIT NONE
  INTEGER :: file_id, ioi
  INTEGER(8) :: OFFSET
  character(1), optional :: hdr
  character(20) :: buffer
  OFFSET = FTELL(file_id)
  REWIND(file_id)
  ioi = 0
  ROWCOUNT = 0
  DO WHILE (ioi==0)
    READ(file_id, *, iostat=ioi) buffer
      if (present(hdr)) THEN
        if (trim(buffer(1:1)) .ne. hdr) ROWCOUNT = ROWCOUNT + 1
      ELSE
        ROWCOUNT = ROWCOUNT + 1
      END IF
  END DO
  ROWCOUNT = ROWCOUNT-1
  REWIND(file_id)
  CALL FSEEK(file_id, OFFSET, 0)
END FUNCTION ROWCOUNT


!==============================================================================
! Function to count number of colums in a file. Takes in the
! UNIT (file pointer) and returns number of columns in the FIRST row in file
! (integer). The function can be used at any point of file reading, the UNIT
! will be moved back where it was in the beginning.
! NOTE:
! Optionally separator can be sent in, e.g. for tab: COLCOUNT(file_id, achar(9))
! This function will still count all consecutive tabs as 1!
! The function does not check if the file has equal amounts of colums in each row
! Empty row or row with only spaces will return 0
!..............................................................................
INTEGER FUNCTION COLCOUNT(file_id, separator)
  IMPLICIT NONE
  INTEGER             :: file_id,ioi,i
  INTEGER(8)          :: OFFSET
  CHARACTER(60000)    :: buffer
  CHARACTER,OPTIONAL  :: separator
  CHARACTER           :: sep
  logical             :: invalue
  invalue = .false.
  if(present(separator)) then
    sep = separator
  else
    sep = " "
  end if
  ioi = 0
  COLCOUNT = 0
  OFFSET = FTELL(file_id)
  REWIND(file_id)
  do i=1,2
    READ(file_id, '(a)', iostat=ioi) buffer
  end do
  if (ioi /= 0) THEN
    COLCOUNT = -9999
  else
    DO i=1,LEN_TRIM(buffer)
      if (((buffer(i:i) /= sep).and.(buffer(i:i) /= achar(9))) .and. (invalue .neqv. .true.)) THEN
        if (COLCOUNT == 0) COLCOUNT = 1
        invalue = .true.
      elseif (((buffer(i:i) == sep).or.(buffer(i:i) == achar(9))).and. (invalue) ) THEN
        invalue = .false.
        COLCOUNT = COLCOUNT +1
      end if
    END DO
  end if

  REWIND(file_id)
  CALL FSEEK(file_id, OFFSET, 0)
END FUNCTION COLCOUNT



END MODULE FILECOUNT_MOD
