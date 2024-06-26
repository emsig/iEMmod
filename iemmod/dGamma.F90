SUBROUTINE dGAMMArad (gam, dgam, gamA, nk, cor, etaV, etaH, zetaH, zeta)
! This function computes Gamma, the vertical wavenumber, for one layer.
!
! Calling arguments:
! gam:  Gamma, the vertical wavenumber, complex array for all horizontal wavenumbers nk
! dgam: Derivative of Gamma, complex array for all horizontal wavenumbers nk
! gamA: Small gamma squared (zetaH*etaV), complex scalar
! nk:   Amount of horizontal wavenumbers, integer scalar
! cor:  Coordinates for the horizontal wavenumber, real array of size nk
! etaV: Material parameter eta for the vertical direction, complex scalar
! etaH: Material parameter eta for the horizontal direction, complex scalar
! zetaH:Material parameter zeta for the horizontal direction, complex scalar
! zeta: Fixed material parameter needed for derivatives, DO NOT CHANGE
!
! Note: In order to compute Gamma Bar, the same function is used, but instead of etaV, etaH, and gamA
!       the arguments zetaV, zetaH and gamB are used.

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk
COMPLEX,INTENT(INOUT) :: gam(nk), dgam(nk), gamA
COMPLEX,INTENT(IN)    :: etaV, etaH, zetaH, zeta
REAL   ,INTENT(IN)    :: cor(nk)
REAL                  :: kx, kx2
INTEGER               :: ikx, iom
COMPLEX               :: fact

gamA=(zetaH*etaV)

fact=(etaH/etaV)

do ikx=1,nk
    kx = cor(ikx);
    kx2=kx*kx;
    gam(ikx) = sqrt(fact*(kx2+gamA));
    dgam(ikx) = (zeta/(2.0*gam(ikx)));
end do

END SUBROUTINE dGAMMArad
