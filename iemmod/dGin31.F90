SUBROUTINE dgin31 (Gin, zrcv, zsrc, etaH, etaV, zetaH, zetaV, gamA, xpos, x, nlogx);
! Calculate the gradient of the direct field in the space domain of component 31
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
! x:        Inline horizontal distance from source, real vector of size nlogx
! nlogx:    Amount of samples in distance vector, integer scalar

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nlogx
COMPLEX,INTENT(INOUT) :: Gin(nlogx)
COMPLEX, INTENT(IN)   :: etaH, etaV, zetaH, zetaV, gamA
REAL   ,INTENT(IN)    :: xpos(nlogx), x(nlogx), zrcv, zsrc
COMPLEX               :: Rbig, G, dGdg, dgds, dGAds, dTdGA, dTdg, dTdet
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

    dTdGA = x(ix)*(zrcv-zsrc)*(gamA/(Rbig**2.0)+3.0*sqrt(gamA)/(Rbig**3.0)+3.0/(Rbig**4.0))/etaV;
    dTdg = x(ix)*(zrcv-zsrc)*(2.0*sqrt(gamA)/(Rbig**2.0)+3.0/(Rbig**3.0))*G/etaV;
    dTdet = -1.0*x(ix)*(zrcv-zsrc)*(gamA/(Rbig**2.0)+3.0*sqrt(gamA)/(Rbig**3.0)+3.0/(Rbig**4.0))*G/(etaV**2.0);

    Gin(ix) = (dTdGA*dGAds) + (dTdg*dgds) + dTdet;
end do

END SUBROUTINE dgin31
