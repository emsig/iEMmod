SUBROUTINE Rplusderiv (Rp, dRp, nk, cor, etaH, Gam, dGam, nz, z, zrcvlay, zsrclay, nlayers, above, zetaH)
! Calculate global reflection coefficient of layers below the receivers
!
! Calling arguments:
! Rp:      Upgoing global reflection coefficient, complex array of size nk times nlayers
! nk:      Amount of wavenumbers, integer scalar
! cor:     Coordinates in the wavenumber domain, real array of size nk
! etaH:    Material parameter eta for the horizontal direction, complex array of size nz
! Gam:     Gamma, the vertical wavenumber, complex array of size nk times nz
! dGam:    Derivative of Gamma, complex array of size nk times nz
! nz:      Amount of layers, integer scalar
! z:       Depth of interfaces, real array of size nz
! zrcvlay: The number of the layer where the receivers are located, integer scalar
! zsrclay: The number of the layer where the source is located, integer scalar
! nlayers: Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! nR:      Amount of layers for dRp needed for the calculation of P, changes according to above, integer scalar
! above:   Indicates if the receivers are above the source (above=1), below the source (above=-1) or
!          in the same layer (above=0), integer scalar
! zetaH:   Always zetaH, do not replace, checks whether Rp or Rpbar is computed
!
! Note: To compute Rpbar replace etaH with zetaH and Gamma with Gamma bar.

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zrcvlay, zsrclay, nlayers, above
COMPLEX,INTENT(INOUT) :: Rp(nk,nz), dRp(nk,nz,nz)
COMPLEX,INTENT(IN)    :: etaH(nz), Gam(nk,nz), dGam(nk,nz), zetaH(nz)
REAL   ,INTENT(IN)    :: cor(nk), z(nz)
REAL                  :: kx
INTEGER               :: iz, izz, izl, izzl, izout, ikx, zrcvlayfort, zsrclayfort, bar, startlayer, finallayer
COMPLEX               :: term1, term2, term3, term4, term5, term6, term7
COMPLEX, DIMENSION(:,:), ALLOCATABLE :: rloc
COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: drloc

if (etaH(1) == zetaH(1)) then
    bar=1;
else
    bar=0;
endif

zrcvlayfort = zrcvlay + 1;
zsrclayfort = zsrclay + 1;

ALLOCATE(rloc (nk,nz-1))
ALLOCATE(drloc (nk,nz,2))

do ikx=1,nk
    do iz=1,nz
        do izz=1,nz
            Rp(ikx,izz) = cmplx(0.0,0.0);
            dRp(ikx,iz,izz) = cmplx(0.0,0.0);
        end do
    end do
end do

do iz=nz-1,1,-1
    do ikx=1,nk
        kx = cor(ikx);
        rloc(ikx,iz) = (etaH(iz+1)*Gam(ikx,iz)-etaH(iz)*Gam(ikx,iz+1))/(etaH(iz+1)*Gam(ikx,iz)+etaH(iz)*Gam(ikx,iz+1));
        if (iz == nz-1) then
            term1 = cmplx(0.0,0.0);
            term2 = cmplx(0.0,0.0);
        else 
            term1 = Rp(ikx,iz+1)*exp(-2.0*Gam(ikx,iz+1)*(z(iz+2)-z(iz+1)));
            term2 = rloc(ikx,iz)*Rp(ikx,iz+1)*exp(-2.0*Gam(ikx,iz+1)*(z(iz+2)-z(iz+1)));
        endif
        Rp(ikx,iz) = (rloc(ikx,iz)+term1)/(1.0+term2);
    end do
end do
do ikx=1,nk
    if (bar == 0) then
        drloc(ikx,nz,2) = ((-etaH(nz-1)*dGam(ikx,nz))*(1.0+rloc(ikx,nz-1))+(Gam(ikx,nz-1)*(1.0-rloc(ikx,nz-1))))&
            /((etaH(nz)*Gam(ikx,nz-1))+(etaH(nz-1)*Gam(ikx,nz)));
    elseif (bar == 1) then
        drloc(ikx,nz,2) = ((-etaH(nz-1)*dGam(ikx,nz))*(1.0+rloc(ikx,nz-1)))&
            /((etaH(nz)*Gam(ikx,nz-1))+(etaH(nz-1)*Gam(ikx,nz)));
    endif
end do
do iz=nz-1,1,-1
    do ikx=1,nk
        if (bar == 0) then
            drloc(ikx,iz,1) = ((etaH(iz+1)*dGam(ikx,iz))*(1.0-rloc(ikx,iz))-(Gam(ikx,iz+1)*(1.0+rloc(ikx,iz))))&
            /((etaH(iz+1)*Gam(ikx,iz))+(etaH(iz)*Gam(ikx,iz+1)));
        elseif (bar == 1) then
            drloc(ikx,iz,1) = ((etaH(iz+1)*dGam(ikx,iz))*(1.0-rloc(ikx,iz)))&
            /((etaH(iz+1)*Gam(ikx,iz))+(etaH(iz)*Gam(ikx,iz+1)));
        endif
        if (iz == 1) then
            drloc(ikx,iz,2)=cmplx(0.0,0.0);
        elseif (iz > 1 .AND. bar == 0) then
            drloc(ikx,iz,2) = ((-etaH(iz-1)*dGam(ikx,iz))*(1.0+rloc(ikx,iz-1))+(Gam(ikx,iz-1)*(1.0-rloc(ikx,iz-1))))&
            /((etaH(iz)*Gam(ikx,iz-1))+(etaH(iz-1)*Gam(ikx,iz)));
        elseif (iz > 1 .AND. bar == 1) then
            drloc(ikx,iz,2) = ((-etaH(iz-1)*dGam(ikx,iz))*(1.0+rloc(ikx,iz-1)))&
            /((etaH(iz)*Gam(ikx,iz-1))+(etaH(iz-1)*Gam(ikx,iz)));
        endif
    end do
end do
do iz=nz,1,-1 !sigma layer number
    do izz=nz-1,1,-1 !Rp layer number
        do ikx=1,nk !wavenumber
            if (izz > iz) then !Rp is deeper than sigma
                dRp(ikx,iz,izz) = cmplx(0.0,0.0);
            elseif (iz == nz-1 .AND. izz == nz-1) then !Rp and sigma of deepest layer with a bottom border
                dRp(ikx,iz,izz) = drloc(ikx,izz,1);
            elseif (izz == iz .AND. iz < nz-1) then !sigma and Rp are in the same layer
                term3 = rloc(ikx,izz)*Rp(ikx,izz+1)*exp(-2.0*Gam(ikx,izz+1)*(z(izz+2)-z(izz+1)));
                term4 = Rp(ikx,izz)*Rp(ikx,izz+1)*exp(-2.0*Gam(ikx,izz+1)*(z(izz+2)-z(izz+1)));
                dRp(ikx,iz,izz)=((1.0-term4)/(1.0+term3))*drloc(ikx,izz,1);
            elseif (iz-izz == 1) then !Rp is one layer above sigma
                if (iz == nz) then !Derivative wrt to lowest layer
                    dRp(ikx,iz,izz) = drloc(ikx,iz,2);
                else
                    term3 = rloc(ikx,izz)*Rp(ikx,izz+1)*exp(-2.0*Gam(ikx,izz+1)*(z(izz+2)-z(izz+1)));
                    term4 = Rp(ikx,izz)*Rp(ikx,izz+1)*exp(-2.0*Gam(ikx,izz+1)*(z(izz+2)-z(izz+1)));
                    term5 = (1.0-(rloc(ikx,izz)*Rp(ikx,izz)))*exp(-2.0*Gam(ikx,izz+1)*(z(izz+2)-z(izz+1)));
                    term6 = Rp(ikx,izz+1)*(1.0-(rloc(ikx,izz)*Rp(ikx,izz)));
                    term7 = (-2.0*(z(izz+2)-z(izz+1)))*exp(-2.0*Gam(ikx,izz+1)*(z(izz+2)-z(izz+1)))*dGam(ikx,izz+1);
                    dRp(ikx,iz,izz) = (((1.0-term4)/(1.0+term3))*drloc(ikx,iz,2))+(((term5)/(1.0+term3))*dRp(ikx,iz,izz+1))&
                        +(((term6)/(1.0+term3))*term7);
                endif
            elseif (iz-izz > 1) then !Rp is at least two layers above sigma
                term3 = rloc(ikx,izz)*Rp(ikx,izz+1)*exp(-2.0*Gam(ikx,izz+1)*(z(izz+2)-z(izz+1)));
                term5 = (1.0-(rloc(ikx,izz)*Rp(ikx,izz)))*exp(-2.0*Gam(ikx,izz+1)*(z(izz+2)-z(izz+1)));
                dRp(ikx,iz,izz) = (((term5)/(1.0+term3))*dRp(ikx,iz,izz+1));
            endif
        end do
    end do
end do
DEALLOCATE(rloc)
DEALLOCATE(drloc)
END SUBROUTINE Rplusderiv
