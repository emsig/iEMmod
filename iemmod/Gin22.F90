SUBROUTINE gin22 (Gin, zrcv, zsrc, etaH, etaV, zetaH, zetaV, gamA, gamB, xpos, y, nlogx)
! Calculate the direct field in the space domain of component 22
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
COMPLEX               :: RbigA, RbigB, GA, GB, term1, term2, term3, term4
REAL                  :: pi, rad
INTEGER               :: ix

pi = 2*acos(0.0);

do ix=1,nlogx
    if (abs(xpos(ix))<= 10**(-6.0)) then
        rad = 10**(-6.0);
    else
        rad = xpos(ix);
    endif
    RbigA = sqrt(rad**2.0+(etaH*(zrcv-zsrc)**2.0/etaV));
    !if (abs(RbigA) <= 10**(-6.0)) then
    !    RbigA = 10**(-6.0);
    !endif
    RbigB = sqrt(rad**2.0+(zetaH*(zrcv-zsrc)**2.0/zetaV));
    !if (abs(RbigB) <= 10**(-6.0)) then
    !    RbigB = 10**(-6.0);
    !endif

    GA = exp(-1.0*sqrt(gamA)*RbigA)/(4.0*pi*RbigA*sqrt(etaH/etaV));
    GB = exp(-1.0*sqrt(gamB)*RbigB)/(4.0*pi*RbigB*sqrt(zetaH/zetaV));
    term1 = ((3.0*y(ix)*y(ix)/(RbigA**2.0)-1.0)*(1.0/(etaV*RbigA**2.0)&
        &+sqrt(gamA)/(etaV*RbigA))+zetaH*y(ix)*y(ix)/(RbigA**2.0))*GA;
    term2 = -1.0*zetaH*GB;
    term3 = -1.0*zetaH*y(ix)*y(ix)/(rad*rad)*(GA-GB);
    term4 = -1.0*sqrt(zetaH)*(2.0*y(ix)*y(ix)/(rad*rad)-1.0)*((exp(-1.0*sqrt(gamA)*RbigA)&
        &-exp(-1.0*sqrt(gamB)*RbigB))/(4.0*pi*sqrt(etaH)*rad*rad));
    Gin(ix) = term1 + term2 + term3 + term4;
end do

END SUBROUTINE gin22
