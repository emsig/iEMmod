SUBROUTINE dgin33 (Gin, zrcv, zsrc, etaH, etaV, zetaH, zetaV, gamA, xpos, nlogx);
! Calculate the gradient of the direct field in the space domain of component 33
!
! Calling arguments:
! Gin:      Direct field in the space domain, complex vector of size nlogx
! zrcv:     Depth of receivers, real scalar
! zsrc:     Depth of source, real scalar
! etaH:     Material parameter eta for the horizontal direction, complex scalar
! etaV:     Material parameter eta for the vertical direction, complex scalar
! zetaH:    Material parameter zeta for the horizontal direction, complex scalar
! zetaV:    Material parameter zeta for the vertical direction, complex scalar
! gamA:     Small gamma squared (zetaH*etaV), complex scalar
! xpos:     Radial horizontal distance from source, real vector of size nlogx
! nlogx:    Amount of samples in distance vector, integer scalar

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nlogx
COMPLEX,INTENT(INOUT) :: Gin(nlogx)
COMPLEX, INTENT(IN)   :: etaH, etaV, zetaH, zetaV, gamA
REAL   ,INTENT(IN)    :: xpos(nlogx), zrcv, zsrc
COMPLEX               :: Rbig, G, term1, term2, dGdg, dgds, dGAds
COMPLEX               :: dt1dg, dt1ds, dt2ds
REAL                  :: pi
INTEGER               :: ix

pi = 2*acos(0.0);

do ix=1,nlogx
    Rbig = sqrt(xpos(ix)**2.0+(zrcv-zsrc)**2.0);
    if (abs(Rbig) <= 10**(-6.0)) then
        Rbig = 10**(-6.0);
    endif
    G = exp(-1.0*sqrt(gamA)*Rbig)/(4.0*pi*Rbig);
    dGdg = (-1.0*exp(-1.0*sqrt(gamA)*Rbig))/(4.0*pi);
    dgds = 0.5*sqrt(zetaH/etaH);
    dGAds = (dGdg*dgds);

    term1 = ((zrcv-zsrc)**2.0*(sqrt(gamA)/Rbig+2.0/(Rbig**2.0)))*(sqrt(gamA)/Rbig+Rbig**(-2.0));
    dt1dg = ((zrcv-zsrc)**2.0)*((2.0*sqrt(gamA))/(Rbig**2.0)+2.0/(Rbig**2.0));
    dt1ds = dt1dg*dgds;

    term2 = ((zrcv-zsrc)/(Rbig**2.0))**2.0-etaH*zetaH;
    dt2ds = -1.0*zetaH;

    Gin(ix) = (G/etaV)*(dt1ds+dt2ds) + dGAds*(term1 + term2)/etaV - (term1 + term2)*G/(etaV**2.0);
end do

END SUBROUTINE dgin33
