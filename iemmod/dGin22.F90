SUBROUTINE dgin22 (Gin, zrcv, zsrc, etaH, etaV, zetaH, zetaV, gamA, gamB, xpos, y, nlogx)
! Calculate the gradient of the direct field in the space domain of component 22
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
! gamB:     Small gamma bar squared (zetaV*etaH), complex scalar
! xpos:     Radial horizontal distance from source, real vector of size nlogx
! y:        Crossline horizontal distance from source, real vector of size nlogx
! nlogx:    Amount of samples in distance vector, integer scalar

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nlogx
COMPLEX,INTENT(INOUT) :: Gin(nlogx)
COMPLEX, INTENT(IN)   :: etaH, etaV, zetaH, zetaV, gamA, gamB
REAL   ,INTENT(IN)    :: xpos(nlogx), y(nlogx), zrcv, zsrc
COMPLEX               :: Rbig, G, dGdg, dgds, dGAds, dt1det, dt1dg, dt1dGA
COMPLEX               :: dt1ds, dt2dGA, dt2ds
REAL                  :: pi, rad
INTEGER               :: ix

pi = 2*acos(0.0);

do ix=1,nlogx
    if (abs(xpos(ix))<= 10**(-6.0)) then
        rad = 10**(-6.0);
    else
        rad = xpos(ix);
    endif
    Rbig = sqrt(rad**2.0+((zrcv-zsrc)**2.0));
    G = exp(-1.0*sqrt(gamA)*Rbig)/(4.0*pi*Rbig);
    dGdg = (-1.0*exp(-1.0*sqrt(gamA)*Rbig))/(4.0*pi);
    dgds = 0.5*sqrt(zetaH/etaH);
    dGAds = (dGdg*dgds);
    dt1det = ((3.0*y(ix)*y(ix)/(Rbig**2.0)-1.0)*(-1.0/((etaH*Rbig)**2.0)&
        &-sqrt(gamA)/((etaH**2.0)*Rbig)))*G;
    dt1dg = ((3.0*y(ix)*y(ix)/(Rbig**2.0)-1.0)*(+1.0/(etaH*Rbig)))*G;
    dt1dGA = ((3.0*y(ix)*y(ix)/(Rbig**2.0)-1.0)*(1.0/(etaV*Rbig**2.0)&
        &+sqrt(gamA)/(etaV*Rbig))+zetaH*y(ix)*y(ix)/(Rbig**2.0));
    dt1ds = dt1det + (dt1dg*dgds) + (dt1dGA*dGAds);
    dt2dGA = -1.0*zetaH;
    dt2ds = (dt2dGA*dGAds);
    Gin(ix) = dt1ds + dt2ds;
end do

END SUBROUTINE dgin22
