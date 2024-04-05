SUBROUTINE Rminderiv (Rm, dRm, nk, cor, etaH, Gam, dGam, nz, z, zrcvlay, zsrclay, nlayers, above, zetaH)
! Calculate global reflection coefficient of layers below the receivers
!
! Calling arguments:
! Rm:      Downgoing global reflection coefficient, complex array of size nk times nlayers
! dRm:     Derivative of downgoing global reflection coefficient, complex array of size nk times nlayers times nlayers
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
! above:   Indicates if the receivers are above the source (above=1), below the source (above=-1) or
!          in the same layer (above=0), integer scalar
! zetaH:   Always zetaH, do not replace, checks whether Rm or Rmbar is computed
!
! Note: To compute Rmbar replace etaH with zetaH and Gamma with Gamma bar.

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zrcvlay, zsrclay, nlayers, above
COMPLEX,INTENT(INOUT) :: Rm(nk,nz), dRm(nk,nz,nz)
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

zrcvlayfort = zrcvlay + 1
zsrclayfort = zsrclay + 1

ALLOCATE(rloc (nk,nz))
ALLOCATE(drloc (nk,nz,2))

Rm(:,:)    = cmplx(0.0,0.0)
dRm(:,:,:) = cmplx(0.0,0.0)

do iz=2,nz,1
    do ikx=1,nk
        kx = cor(ikx);
        rloc(ikx,iz) = (etaH(iz-1)*Gam(ikx,iz)-etaH(iz)*Gam(ikx,iz-1))/(etaH(iz-1)*Gam(ikx,iz)+etaH(iz)*Gam(ikx,iz-1));
        if (iz == 2) then
            term1 = cmplx(0.0,0.0);
            term2 = cmplx(0.0,0.0);
        else 
            term1 = Rm(ikx,iz-1)*exp(-2.0*Gam(ikx,iz-1)*(z(iz)-z(iz-1)));
            term2 = rloc(ikx,iz)*Rm(ikx,iz-1)*exp(-2.0*Gam(ikx,iz-1)*(z(iz)-z(iz-1)));
        endif
        Rm(ikx,iz) = (rloc(ikx,iz)+term1)/(1.0+term2);
    end do
end do
do iz=1,nz,1
    do ikx=1,nk
        if (iz == 1) then
            if (iz < nz .AND. bar == 0) then
                drloc(ikx,iz,2) = (((-etaH(iz+1)*dGam(ikx,iz)*(1.0+rloc(ikx,iz+1))))+(Gam(ikx,iz+1)*(1.0-rloc(ikx,iz+1))))&
                    /((etaH(iz)*Gam(ikx,iz+1))+(etaH(iz+1)*Gam(ikx,iz)));
            elseif (iz < nz .AND. bar == 1) then
                drloc(ikx,iz,2) = ((-etaH(iz+1)*dGam(ikx,iz)*(1.0+rloc(ikx,iz+1))))&
                    /((etaH(iz)*Gam(ikx,iz+1))+(etaH(iz+1)*Gam(ikx,iz)));
            endif
        else
            if (bar == 0) then
                drloc(ikx,iz,1) = (((etaH(iz-1)*dGam(ikx,iz))*(1.0-rloc(ikx,iz)))-(Gam(ikx,iz-1)*(1.0+rloc(ikx,iz))))&
                    /((etaH(iz-1)*Gam(ikx,iz))+(etaH(iz)*Gam(ikx,iz-1)));
            elseif (bar == 1) then
                drloc(ikx,iz,1) = ((etaH(iz-1)*dGam(ikx,iz)*(1.0-rloc(ikx,iz))))&
                    /((etaH(iz-1)*Gam(ikx,iz))+(etaH(iz)*Gam(ikx,iz-1)));
            endif
            if (iz == nz) then
                drloc(ikx,iz,2)=cmplx(0.0,0.0);
            elseif (iz < nz .AND. bar == 0) then
                drloc(ikx,iz,2) = (((-etaH(iz+1)*dGam(ikx,iz)*(1.0+rloc(ikx,iz+1))))+(Gam(ikx,iz+1)*(1.0-rloc(ikx,iz+1))))&
                    /((etaH(iz)*Gam(ikx,iz+1))+(etaH(iz+1)*Gam(ikx,iz)));
            elseif (iz < nz .AND. bar == 1) then
                drloc(ikx,iz,2) = ((-etaH(iz+1)*dGam(ikx,iz)*(1.0+rloc(ikx,iz+1))))&
                    /((etaH(iz)*Gam(ikx,iz+1))+(etaH(iz+1)*Gam(ikx,iz)));
            endif
        endif
    end do
end do
do iz=1,nz,1 !sigma layer number
    do izz=2,nz,1 !Rm layer number
        do ikx=1,nk !wavenumber
            if (izz < iz) then !Rm is less deep than sigma
                dRm(ikx,iz,izz) = cmplx(0.0,0.0);
            elseif (iz == 2 .AND. izz == 2) then !Rm and sigma of the most shallow layer
                dRm(ikx,iz,izz) = drloc(ikx,izz,1);
            elseif (izz == iz .AND. iz > 2) then !sigma and Rm are in the same layer
                term3 = rloc(ikx,izz)*Rm(ikx,izz-1)*exp(-2.0*Gam(ikx,izz-1)*(z(izz)-z(izz-1)));
                term4 = Rm(ikx,izz)*Rm(ikx,izz-1)*exp(-2.0*Gam(ikx,izz-1)*(z(izz)-z(izz-1)));
                dRm(ikx,iz,izz)=((1.0-term4)/(1.0+term3))*drloc(ikx,izz,1);
            elseif (izz-iz == 1) then !Rm is one layer below sigma
                if (iz == 1) then !Derivative wrt to top layer
                    dRm(ikx,iz,izz) = drloc(ikx,iz,2);
                else
                    term3 = rloc(ikx,izz)*Rm(ikx,izz-1)*exp(-2.0*Gam(ikx,izz-1)*(z(izz)-z(izz-1)));
                    term4 = Rm(ikx,izz)*Rm(ikx,izz-1)*exp(-2.0*Gam(ikx,izz-1)*(z(izz)-z(izz-1)));
                    term5 = (1.0-(rloc(ikx,izz)*Rm(ikx,izz)))*exp(-2.0*Gam(ikx,izz-1)*(z(izz)-z(izz-1)));
                    term6 = Rm(ikx,izz-1)*(1.0-(rloc(ikx,izz)*Rm(ikx,izz)));
                    term7 = (-2.0*(z(izz)-z(izz-1)))*exp(-2.0*Gam(ikx,izz-1)*(z(izz)-z(izz-1)))*dGam(ikx,izz-1);
                    dRm(ikx,iz,izz) = (((1.0-term4)/(1.0+term3))*drloc(ikx,iz,2))+(((term5)/(1.0+term3))*dRm(ikx,iz,izz-1))&
                        +(((term6)/(1.0+term3))*term7);
                endif
            elseif (izz-iz > 1) then !Rm is at least two layers below sigma
                term3 = rloc(ikx,izz)*Rm(ikx,izz-1)*exp(-2.0*Gam(ikx,izz-1)*(z(izz)-z(izz-1)));
                term5 = (1.0-(rloc(ikx,izz)*Rm(ikx,izz)))*exp(-2.0*Gam(ikx,izz-1)*(z(izz)-z(izz-1)));
                dRm(ikx,iz,izz) = (((term5)/(1.0+term3))*dRm(ikx,iz,izz-1));
            endif
        end do
    end do
end do
DEALLOCATE(rloc)
DEALLOCATE(drloc)
END SUBROUTINE Rminderiv
