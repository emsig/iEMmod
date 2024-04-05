SUBROUTINE gridit_xx_lin (Ptot,xtotrad1,xtotrad2,dx,nxh,dy,nyh,xpos,nlogx,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,gamB,above,xdirect)
! Compute the data of one quadrant of the grid based on the radial data using linear interpolation
! for component 11 (inline oriented electric receiver, inline oriented electric source).
! 
! Calling arguments:
! Ptot:     After completion of this function, this array will contain one quadrant of the EM-field in the space domain, 
!           complex array of size nxh times nyh
! xtotrad1: 1st Hankel transformed field-term, complex array of size nlogx
! xtotrad2: 2nd Hankel transformed field-term, complex array of size nlogx
! dx:       Sampling in inline direction (x-direction), real scalar
! nxh:      Number of samples of one quadrant in inline direction (x-direction), integer scalar
! dy:       Sampling in crossline direction (y-direction), real scalar
! nyh:      Number of samples of one quadrant in crossline direction (y-direction), integer scalar
! xpos:     Logarithmic coordinate vector in the space domain, real array of size nlogx
! nlogx:    Amount of samples in logarithmic coordinate vector, integer scalar
! zrcv:     Depth of receivers, real scalar
! zsrc:     Depth of source, real scalar
! etaH:     Material parameter eta for the horizontal direction, complex scalar
! etaV:     Material parameter eta for the vertical direction, complex scalar
! zetaH:    Material parameter zeta for the horizontal direction, complex scalar
! zetaV:    Material parameter zeta for the vertical direction, complex scalar
! gamA:     Small gamma squared (zetaH*etaV), complex scalar
! gamB:     Small gamma bar squared (zetaV*etaH), complex scalar
! above:    Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!           in the same layer (above=0), integer scalar
! xdirect:  Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: Ptot(nxh,nyh)
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx), xtotrad2(nlogx)
COMPLEX,INTENT(IN)    :: etaH, etaV, zetaH, zetaV, gamA, gamB
REAL   ,INTENT(IN)    :: xpos(nlogx), dx, dy, zrcv, zsrc
INTEGER,INTENT(IN)    :: nxh, nyh, nlogx, above, xdirect
COMPLEX               :: gin
REAL                  :: x, y, rad, radfilo, radfihi, weightlo, weighthi, fact
REAL                  :: thisdx
INTEGER               :: ix, iy, ipos, lopos, hipos

do iy=1,nyh
    ! This is the current y-position for which we want to compute the field
    y = (iy-1)*dy;
    lopos = 1;
    do ix=1,nxh
        ! This is the currend x-position for which we want to compute the field
        x = (ix-1)*dx;
        ! Using the x- and y-position we get the radial distance which can be projected on the line of data in the space domain.
        rad = sqrt(x*x+y*y);
        ! Search for the elements on the Hankel-transformed line which flank the point on which we want to know the field.
        do ipos=lopos,nlogx
            ! In the first iteration this loop starts at lopos=1 as lopos was initialized like that.
            ! The next iteration will look for a point with a greater radial distance. There is therefore no need
            ! to start looking again from lopos=1. It has to be further away than the previous point. 
            if (xpos(ipos) > rad .OR. ipos == nlogx) then
                hipos = ipos;
                lopos = ipos -1;
                exit
            endif
        end do
        ! Prepare the linear interpolation of the two points to get the field at the desired point.
        thisdx = xpos(hipos) - xpos(lopos);
        weightlo = (xpos(hipos)-rad)/thisdx;
        weighthi = (rad-xpos(lopos))/thisdx;
        ! Fact is a symmetry-parameter that needs to be taken into account.
        if (ix==1) then
            fact = -1.0;
        else
            fact = cos(2.0*atan(y/x));
        endif
        ! Linear interpolation and symmetry factor are combined
        Ptot(ix,iy) = (weightlo*xtotrad1(lopos)+weighthi*xtotrad1(hipos));
        Ptot(ix,iy) = Ptot(ix,iy) - fact*(weightlo*xtotrad2(lopos)+weighthi*xtotrad2(hipos));
        ! If xdirect=1 the direct field has not been computed in the wavenumber domain and needs to be added now. 
        if (above==0 .AND. xdirect==1) then
            ! The function gin11 computes the direct field in the space domain. 
            CALL gin11(gin,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,gamB,rad,x,1);
            ! Add the direct field to the Hankel-transformed and regridded space-domain solution.
            Ptot(ix,iy) = Ptot(ix,iy) + gin;
        endif
    end do
end do
END SUBROUTINE gridit_xx_lin
      
SUBROUTINE gridit_xy_lin (Ptot,xtotrad,dx,nxh,dy,nyh,xpos,nlogx,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,gamB,above,xdirect)
! Compute the data of one quadrant of the grid based on the radial data using linear interpolation
! for component 12 (inline oriented electric receiver, crossline oriented electric source).
! 
! Calling arguments:
! Ptot:    After completion of this function, this array will contain one quadrant of the EM-field in the space domain, 
!          complex array of size nxh times nyh
! xtotrad: 1st Hankel transformed field-term, complex array of size nlogx
! dx:      Sampling in inline direction (x-direction), real scalar
! nxh:     Number of samples of one quadrant in inline direction (x-direction), integer scalar
! dy:      Sampling in crossline direction (y-direction), real scalar
! nyh:     Number of samples of one quadrant in crossline direction (y-direction), integer scalar
! xpos:    Logarithmic coordinate vector in the space domain, real array of size nlogx
! nlogx:   Amount of samples in logarithmic coordinate vector, integer scalar
! zrcv:    Depth of receivers, real scalar
! zsrc:    Depth of source, real scalar
! etaH:    Material parameter eta for the horizontal direction, complex scalar
! etaV:    Material parameter eta for the vertical direction, complex scalar
! zetaH:   Material parameter zeta for the horizontal direction, complex scalar
! zetaV:   Material parameter zeta for the vertical direction, complex scalar
! gamA:    Small gamma squared (zetaH*etaV), complex scalar
! gamB:    Small gamma bar squared (zetaV*etaH), complex scalar
! above:   Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!          in the same layer (above=0), integer scalar
! xdirect: Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: Ptot(nxh,nyh)
COMPLEX,INTENT(IN)    :: xtotrad(nlogx)
COMPLEX,INTENT(IN)    :: etaH, etaV, zetaH, zetaV, gamA, gamB
REAL   ,INTENT(IN)    :: xpos(nlogx), dx, dy, zrcv, zsrc
INTEGER,INTENT(IN)    :: nxh, nyh, nlogx, above, xdirect
COMPLEX               :: gin
REAL                  :: x, y, rad, radfilo, radfihi, weightlo, weighthi, fact
REAL                  :: thisdx
INTEGER               :: ix, iy, ipos, lopos, hipos

do iy=1,nyh
    y = (iy-1)*dy;
    lopos = 1;
    do ix=1,nxh
        x = (ix-1)*dx;
        rad = sqrt(x*x+y*y);
        do ipos=lopos,nlogx
            if (xpos(ipos) > rad .OR. ipos == nlogx) then
                hipos = ipos;
                lopos = ipos -1;
                exit
            endif
        end do
        thisdx = xpos(hipos) - xpos(lopos);
        weightlo = (xpos(hipos)-rad)/thisdx;
        weighthi = (rad-xpos(lopos))/thisdx;
        if (ix==1) then
            fact = 0.0;
        else
            fact = sin(2.0*atan(y/x));
        endif
        Ptot(ix,iy) = (-1.0)*fact*(weightlo*xtotrad(lopos)+weighthi*xtotrad(hipos));
        if (above==0 .AND. xdirect==1) then
            CALL gin12(gin,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,gamB,rad,x,y,1);
            Ptot(ix,iy) = Ptot(ix,iy) + gin;
        endif
    end do
end do
END SUBROUTINE gridit_xy_lin

SUBROUTINE gridit_xz_lin (Ptot,xtotrad1,dx,nxh,dy,nyh,xpos,nlogx,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,above,xdirect)
! Compute the data of one quadrant of the grid based on the radial data using linear interpolation
! for component 13 (inline oriented electric receiver, vertical electric source).
! 
! Calling arguments:
! Ptot:     After completion of this function, this array will contain one quadrant of the EM-field in the space domain, 
!           complex array of size nxh times nyh
! xtotrad1: 1st Hankel transformed field-term, complex array of size nlogx
! dx:       Sampling in inline direction (x-direction), real scalar
! nxh:      Number of samples of one quadrant in inline direction (x-direction), integer scalar
! dy:       Sampling in crossline direction (y-direction), real scalar
! nyh:      Number of samples of one quadrant in crossline direction (y-direction), integer scalar
! xpos:     Logarithmic coordinate vector in the space domain, real array of size nlogx
! nlogx:    Amount of samples in logarithmic coordinate vector, integer scalar
! zrcv:     Depth of receivers, real scalar
! zsrc:     Depth of source, real scalar
! etaH:     Material parameter eta for the horizontal direction, complex scalar
! etaV:     Material parameter eta for the vertical direction, complex scalar
! zetaH:    Material parameter zeta for the horizontal direction, complex scalar
! zetaV:    Material parameter zeta for the vertical direction, complex scalar
! gamA:     Small gamma squared (zetaH*etaV), complex scalar
! above:    Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!           in the same layer (above=0), integer scalar
! xdirect:  Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: Ptot(nxh,nyh)
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx), etaH, etaV, zetaH, zetaV, gamA
REAL   ,INTENT(IN)    :: xpos(nlogx), dx, dy, zrcv, zsrc
INTEGER,INTENT(IN)    :: nxh, nyh, nlogx, above, xdirect
COMPLEX               :: gin
REAL                  :: x, y, rad, radfilo, radfihi, weightlo, weighthi, fact
REAL                  :: thisdx
INTEGER               :: ix, iy, ipos, lopos, hipos

do iy=1,nyh
    y = (iy-1)*dy;
    lopos = 1;
    do ix=1,nxh
        x = (ix-1)*dx;
        rad = sqrt(x*x+y*y);
        do ipos=lopos,nlogx
            if (xpos(ipos) > rad .OR. ipos == nlogx) then
                hipos = ipos;
                lopos = ipos -1;
                exit
            endif
        end do
        thisdx = xpos(hipos) - xpos(lopos);
        weightlo = (xpos(hipos)-rad)/thisdx;
        weighthi = (rad-xpos(lopos))/thisdx;
        if (ix==1) then
            fact = 0.0;
        else
            fact = cos(atan(y/x));
        endif
        Ptot(ix,iy) = fact*(weightlo*xtotrad1(lopos)+weighthi*xtotrad1(hipos));
        if (above==0 .AND. xdirect==1) then
            CALL gin13(gin,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,rad,x,1);
            Ptot(ix,iy) = Ptot(ix,iy) - gin;
        endif
    end do
end do
END SUBROUTINE gridit_xz_lin
