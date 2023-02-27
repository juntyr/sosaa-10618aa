MODULE multimodal_mod
IMPLICIT NONE

contains

    !====================================================================================
    ! Calculates multimodal particle size distribution, used for intialization. Modevector
    ! is of length <number of modes>*3 and contains the ount median diameter, (CMD),
    ! geometric standard deviation (GSD) and relative size in total particle count
    ! (relative sizes must add up to 1). diameters is the diameter vector of the model
    ! (the bins), psd is the output vector and N is the total particle count in the output
    ! vector
    !....................................................................................
    subroutine Multimodal(modevector, diameters, psd, N)
        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)
        real(dp), PARAMETER :: pi = ACOS(-1d0)

        real(dp), INTENT(in)    :: modevector(:)
        real(dp), INTENT(in)    :: diameters(:)
        real(dp), INTENT(inout) :: psd(:)
        real(dp), allocatable   :: x(:), sumv(:)
        real(dp)                :: N ,mu,sig           ! size factor and total count
        INTEGER                 :: ii,nModes, uprInd, nb

        nb = size(psd)

        IF (MODULO(size(modevector),3) .ne. 0) THEN
            print*, 'The vector for modes is not correct'
            STOP
        ELSE
            nModes = size(modevector)/3
            print*, 'Building PSD from ',nModes,' modes'
        END IF

        ALLOCATE(x(nb))
        ALLOCATE(sumv(nb))
        sumv = 0d0

        x = LOG10(diameters)

        DO ii = 1,nModes
            mu   = LOG10(modevector((ii-1)*3+1))
            sig  = modevector((ii-1)*3+2)
            sumv = sumv + modevector((ii-1)*3+3) * gauss(x, mu, sig)
        END DO

        psd = N * sumv / (sum( sumv ))
        uprInd = nb+1 - min(nb-1, 8)
        WHERE (psd(uprInd:) > 1d-12) psd(uprInd:) = 0d0

    end subroutine Multimodal


    pure function gauss(x,mu,sig) result(zz)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)
        real(dp), PARAMETER :: pi = ACOS(-1d0)

        real(dp), INTENT(IN) :: x(:)
        real(dp), INTENT(IN) :: mu,sig
        real(dp), allocatable :: zz(:)
        allocate(zz(size(x)))
        zz = exp(-(x-mu)**2/(2*sig**2))/sqrt(2*pi*sig**2)
    end function gauss



END MODULE multimodal_mod
