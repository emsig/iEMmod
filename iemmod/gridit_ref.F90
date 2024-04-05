SUBROUTINE gridit_ref_lin (Ptot,xtotrad1,dx,nxh,dy,nyh,xpos,nlogx)
! Compute the data of one quadrant of the grid based on the radial data using linear interpolation
! for component 77 and 88 (TM-mode and TE-mode reflection response, respectively).
! 
! Calling arguments:
! Ptot:     After completion of this function, this array will contain one quadrant of the reflection response in the space domain, 
!           complex array of size nxh times nyh
! xtotrad1: 1st Hankel transformed field-term, complex array of size nlogx
! dx:       Sampling in inline direction (x-direction), real scalar
! nxh:      Number of samples of one quadrant in inline direction (x-direction), integer scalar
! dy:       Sampling in crossline direction (y-direction), real scalar
! nyh:      Number of samples of one quadrant in crossline direction (y-direction), integer scalar
! xpos:     Logarithmic coordinate vector in the space domain, real array of size nlogx
! nlogx:    Amount of samples in logarithmic coordinate vector, integer scalar

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: Ptot(nxh,nyh)
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx)
REAL   ,INTENT(IN)    :: xpos(nlogx), dx, dy
INTEGER,INTENT(IN)    :: nxh, nyh, nlogx
REAL                  :: x, y, rad, radfilo, radfihi, weightlo, weighthi, fact
REAL                  :: thisdx
INTEGER               :: ix, iy, lopos, hipos, minel
INTEGER, DIMENSION(1) :: temp

do iy=1,nyh
    y = (iy-1)*dy;
    do ix=1,nxh
        x = (ix-1)*dx;
        rad = sqrt(x*x+y*y);
        temp = minloc(abs(xpos-rad));
        minel = temp(1);
        if (xpos(minel)<=rad) then
            lopos = minel;
            hipos = minel + 1;
        else
            hipos = minel;
            lopos = minel - 1;
        endif
        if (lopos == hipos) then
            hipos = hipos + 1;
        endif
        thisdx = xpos(hipos) - xpos(lopos);
        weightlo = (xpos(hipos)-rad)/thisdx;
        weighthi = (rad-xpos(lopos))/thisdx;
        Ptot(ix,iy) = weightlo*xtotrad1(lopos)+weighthi*xtotrad1(hipos);
    end do
end do
END SUBROUTINE gridit_ref_lin
