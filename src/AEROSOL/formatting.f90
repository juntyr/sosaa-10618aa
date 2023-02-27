MODULE FORMATTING
IMPLICIT NONE
CHARACTER(56) :: FMT_TIME    = '("+.",t18,82("."),t3,"Time: ",a,t100, "+")'
CHARACTER(56) :: FMT_CVU     = '("| ",t3, a,es10.3,t100, "|")'
CHARACTER(56) :: FMT10_CVU   = '("| ",t3, a,t13, es10.3,a,t100, "|")'
CHARACTER(56) :: FMT30_CVU   = '("| ",t3, a,t33, es10.3,a,t100, "|")'
CHARACTER(86) :: FMT10_2CVU  = '("| ",t3, a,t13, es10.3,a,t35,a,es10.3,a, t100, "|")'
CHARACTER(86) :: FMT10_3CVU  = '("| ",t3, a,t13, es10.3,a,t35,a,es10.3,a,t65,a,es10.3,a, t100, "|")'
CHARACTER(86) :: FMT_INTRMT  = '(" +",t4, 96("."), t4, a, t100, "+")'
CHARACTER(86) :: FMT_LEND    = '("+",t2, 98("."), t100, "+")'
CHARACTER(86) :: FMT_LOOPEND = '("+",t2, 49("--"), t100, "+")'
CHARACTER(86) :: FMT_SUB     = '("| ",t5,":",8("."), a,t100, "|")'
CHARACTER(86) :: FMT_MSG     = '("| ",t3,a,t100, "|")'
CHARACTER(86) :: FMT_HDR     = '("+",t2, 98(".") t3,a,t100, "+")'
CHARACTER(86) :: FMT_WARN0   = '("| WARNING:",t12,88("~"),"+", t12, a)'
CHARACTER(86) :: FMT_WARN1   = '("| WARNING:",t12,88("~"),"+", t12, a,f0.4)'
CHARACTER(86) :: FMT_NOTE0   = '("| NOTE:",t9,91("~"),"+", t9, a)'
CHARACTER(86) :: FMT_NOTE1   = '("| NOTE:",t9,91("~"),"+", t9, a,f0.4)'
CHARACTER(86) :: FMT_FAT0    = '("| FATAL ERROR:",t16,84("~"),"+", t16, a)'

contains
    
subroutine handle_file_io(ioi, file, halt)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: ioi
    CHARACTER(*), INTENT(in) :: file
    CHARACTER(100) :: fmt
    CHARACTER(*), OPTIONAL, INTENT(in) :: halt
    if (PRESENT(halt)) THEN
        fmt = FMT_FAT0
    else
        fmt = FMT_WARN0
    end if
    if (ioi /= 0) THEN
        print fmt, 'Could not open '//TRIM(file)//', does it exist?'
        if (PRESENT(halt)) THEN
            print FMT_SUB, halt
            print FMT_LEND
            stop
        END IF
    END if
end subroutine handle_file_io


END MODULE FORMATTING
