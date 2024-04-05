SUBROUTINE dgridit_zx_lin (Ptot,xtotrad1,dx,nxh,dy,nyh,xpos,nlogx,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,above,xdirect,df)
! Compute the data of one quadrant of the grid based on the radial data using linear interpolation
! for component 31 (vertical electric receiver, inline oriented electric source).
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
COMPLEX               :: dgin
REAL                  :: x, y, rad, radfilo, radfihi, weightlo, weighthi, fact
REAL                  :: thisdx
INTEGER               :: ix, iy, ipos, lopos, hipos, df

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
            fact = -1.0*cos(atan(y/x));
        endif
        Ptot(ix,iy) = fact*(weightlo*xtotrad1(lopos)+weighthi*xtotrad1(hipos));
        if (df == 1) then
            if (above==0 .AND. xdirect==1) then
                CALL dgin31(dgin,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,rad,x,1);
                Ptot(ix,iy) = Ptot(ix,iy) - dgin;
            endif
        else
            
        endif
    end do
end do
END SUBROUTINE dgridit_zx_lin

SUBROUTINE dgridit_zy_lin (Ptot,xtotrad1,dx,nxh,dy,nyh,xpos,nlogx,zrcv,zsrc,&
    etaH,etaV,zetaH,zetaV,gamA,above,xdirect,df)
! Compute the data of one quadrant of the grid based on the radial data using linear interpolation
! for component 32 (vertical electric receiver, crossline oriented electric source).
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
INTEGER,INTENT(IN)    :: nxh, nyh, nlogx, above, xdirect, df
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
            fact = -1.0;
        else
            fact = -1.0*sin(atan(y/x));
        endif
        Ptot(ix,iy) = fact*(weightlo*xtotrad1(lopos)+weighthi*xtotrad1(hipos));
        if (df == 1) then
            if (above==0 .AND. xdirect==1) then
                CALL dgin32(gin,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,rad,y,1);
                Ptot(ix,iy) = Ptot(ix,iy) - gin;
            endif
        endif
    end do
end do
END SUBROUTINE dgridit_zy_lin

SUBROUTINE dgridit_zz_lin (Ptot,xtotrad1,dx,nxh,dy,nyh,xpos,nlogx,zrcv,zsrc,&
    etaH,etaV,zetaH,zetaV,gamA,above,xdirect,df)
! Compute the data of one quadrant of the grid based on the radial data using linear interpolation
! for component 33 (vertical electric receiver, vertical electric source).
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
INTEGER,INTENT(IN)    :: nxh, nyh, nlogx, above, xdirect, df
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
        Ptot(ix,iy) = weightlo*xtotrad1(lopos)+weighthi*xtotrad1(hipos);
        if (df == 1) then
            if (above==0 .AND. xdirect==1) then
                CALL dgin33(gin,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,rad,1);
                Ptot(ix,iy) = Ptot(ix,iy) + gin;
            endif
        endif
    end do
end do
END SUBROUTINE dgridit_zz_lin
