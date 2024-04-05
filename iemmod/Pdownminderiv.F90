SUBROUTINE Pdownminderiv (Pdown, dPdown, Rp, Rm, dRp, dRm, Gam, dGam, nk, cor, nz, z, zsrc, zrcvlay, zsrclay, nlayers, above);
! Calculate Pdown-
!
! Calling arguments:
! Pdown:   Downgoing field, complex array of size nk
! dPdown:  Derivative of downgoing field, complex array of size nk times (nlayers+1)
! Rp:      Upgoing global reflection coefficient, complex array of size nk times (nlayers+1)
! Rm:      Downgoing global reflection coefficient, complex array of size nk times (nlayers+1)
!          NOTE: both Rm and Rp have an extra reflection coefficient added for a calculation.
!          In Rm this is the first entry and in Rp the last therefore Rp(ikx,1) and Rm(ikx,2) are the same layer
! dRp:     Derivative of upgoing reflection coefficient, complex array of size nk times (nlayers+1) times (nlayers+1)
! dRm:     Derivative of downgoing reflection coefficient, complex array of size nk times (nlayers+1) times (nlayers+1)
!          NOTE: Just like Rm and Rp there is an extra coefficient added that is needed for calculation
!          Therefore dRp(ikx,1,1) is in the same layer as dRm(ikx,2,2)
! Gam:     Gamma, the vertical wavenumber, complex array of size nk times nz
! dGam:    Derivative of the vertical wavenumber, complex array of size nk times nz
! nk:      Amount of wavenumbers, integer scalar
! cor:     Coordinates in the wavenumber domain, real array of size nk
! nz:      Amount of layers, integer scalar
! z:       Depth of interfaces, real array of size nz
! zsrc:    Depth of source, real scalar
! zrcvlay: The number of the layer where the receivers are located, integer scalar
! zsrclay: The number of the layer where the source is located, integer scalar
! nlayers: Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! above:   Indicates if the receivers are above the source (above=1), below the source (above=-1) or
!          in the same layer (above=0), integer scalar

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zrcvlay, zsrclay, nlayers, above
COMPLEX,INTENT(INOUT) :: Pdown(nk), dPdown(nk,nz), Rp(nk,nz), Rm(nk,nz), dRp(nk,nz,nz), dRm(nk,nz,nz)
COMPLEX,INTENT(IN)    :: Gam(nk,nz), dGam(nk,nz)
REAL   ,INTENT(IN)    :: cor(nk), z(nz)
REAL                  :: zsrc, ds, dp, dm
INTEGER               :: iz, ikx, zsrclayfort, zrcvlayfort
COMPLEX               :: Ms, term1, term2, term3, term4, term5, term6, term7, term8, term9, term10
COMPLEX               :: term11, term12, term13


zsrclayfort = zsrclay + 1;
zrcvlayfort = zrcvlay + 1;
ds = z(zsrclayfort+1)-z(zsrclayfort);
dp = z(zsrclayfort+1)-zsrc;
dm = zsrc-z(zsrclayfort);

if (above == 0) then ! receivers in the source layer
    do ikx=1,nk
        if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
            Pdown(ikx) = 0.0;
        elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
            Ms = 1.0;
            Pdown(ikx) = Rm(ikx,zsrclayfort)/Ms*(exp(-1.0*Gam(ikx,zsrclayfort)*dm));
        else ! If the source and the receivers are in any layer in between
            Ms = 1.0-Rm(ikx,zsrclayfort)*Rp(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
            Pdown(ikx) = Rm(ikx,zsrclayfort)/Ms*(exp(-1.0*Gam(ikx,zsrclayfort)*dm)&
                -Rp(ikx,zsrclayfort)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dp)));
        endif
        do iz=1,nz
            if (iz == zsrclayfort) then
                if (zsrclayfort == 1) then
                    dPdown(ikx,iz) = 0.0;
                elseif (zsrclayfort == nz) then
                    Ms = 1.0;
                    term1  = (1.0/Ms)*(exp(-1.0*Gam(ikx,zsrclayfort)*dm));
                    term2  = ((-1.0*Rm(ikx,zsrclayfort))/(Ms*Ms))*(exp(-1.0*Gam(ikx,zsrclayfort)*dm));
                    term3  = ((-1.0*Rm(ikx,zsrclayfort))/Ms)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dp));
                    term4  = Rm(ikx,iz)/Ms;
                    term7  = -Rm(ikx,iz)*exp(-2.0*Gam(ikx,iz)*ds);
                    term9  = -2.0*ds*exp(-2.0*Gam(ikx,iz)*ds)*dGam(ikx,iz);
                    term10 = -1.0*dm*exp(-1.0*Gam(ikx,iz)*dm)*dGam(ikx,iz);
                    term11 = -1.0*(ds+dp)*exp(-1.0*Gam(ikx,iz)*(ds+dp))*dGam(ikx,iz);
                    dPdown(ikx,iz) = (term1*dRm(ikx,iz,iz))+(term4*term10);
                else
                    Ms = 1.0-Rm(ikx,zsrclayfort)*Rp(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
                    term1  = (1.0/Ms)*(exp(-1.0*Gam(ikx,zsrclayfort)*dm)&
                        -Rp(ikx,zsrclayfort)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dp)));
                    term2  = ((-1.0*Rm(ikx,zsrclayfort))/(Ms*Ms))*(exp(-1.0*Gam(ikx,zsrclayfort)*dm)&
                        -Rp(ikx,zsrclayfort)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dp)));
                    term3  = ((-1.0*Rm(ikx,zsrclayfort))/Ms)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dp));
                    term4  = Rm(ikx,iz)/Ms;
                    term5  = -1.0*((Rp(ikx,iz)*Rm(ikx,iz))/Ms);
                    term6  = -Rp(ikx,iz)*exp(-2.0*Gam(ikx,iz)*ds);
                    term7  = -Rm(ikx,iz)*exp(-2.0*Gam(ikx,iz)*ds);
                    term8  = -1.0*Rm(ikx,iz)*Rp(ikx,iz);
                    term9  = -2.0*ds*exp(-2.0*Gam(ikx,iz)*ds)*dGam(ikx,iz);
                    term10 = -1.0*dm*exp(-1.0*Gam(ikx,iz)*dm)*dGam(ikx,iz);
                    term11 = -1.0*(ds+dp)*exp(-1.0*Gam(ikx,iz)*(ds+dp))*dGam(ikx,iz);
                    term12 = (term7*dRp(ikx,iz,iz))&
                        +(term6*dRm(ikx,iz,iz))+(term8*term9);
                    dPdown(ikx,iz) = (term1*dRm(ikx,iz,iz))+&
                        (term3*dRp(ikx,iz,iz))+(term2*term12)&
                        +(term4*term10)+(term5*term11);
                endif
            else
                if (iz < zsrclayfort) then ! No contribution by dRp
                    if (zsrclayfort == 1) then
                        dPdown(ikx,iz) = 0.0;
                    elseif (zsrclayfort == nz) then
                        Ms = 1.0;
                        term1  = (1.0/Ms)*(exp(-1.0*Gam(ikx,zsrclayfort)*dm));
                        dPdown(ikx,iz) = (term1*dRm(ikx,iz,zsrclayfort));
                    else
                        Ms = 1.0-Rm(ikx,zsrclayfort)*Rp(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
                        term1  = (1.0/Ms)*(exp(-1.0*Gam(ikx,zsrclayfort)*dm)&
                            -Rp(ikx,zsrclayfort)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dp)));
                        term2  = ((-1.0*Rm(ikx,zsrclayfort))/(Ms*Ms))*(exp(-1.0*Gam(ikx,zsrclayfort)*dm)&
                            -Rp(ikx,zsrclayfort)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dp)));
                        term6  = -Rp(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
                        term13 = (term6*dRm(ikx,iz,zsrclayfort));
                        dPdown(ikx,iz) = (term1*dRm(ikx,iz,zsrclayfort))+(term2*term13);
                    endif
                elseif (iz > zsrclayfort) then !No contribution by dRm
                    if (zsrclayfort == 1) then
                        dPdown(ikx,iz) = 0.0;
                    elseif (zsrclayfort == nz) then
                        dPdown(ikx,iz) = 0.0;
                    else
                        Ms = 1.0-Rm(ikx,zsrclayfort)*Rp(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
                        term1  = (1.0/Ms)*(exp(-1.0*Gam(ikx,zsrclayfort)*dm)&
                            -Rp(ikx,zsrclayfort)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dp)));
                        term2  = ((-1.0*Rm(ikx,zsrclayfort))/(Ms*Ms))*(exp(-1.0*Gam(ikx,zsrclayfort)*dm)&
                            -Rp(ikx,zsrclayfort)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dp)));
                        term7  = -Rm(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
                        term13 = (term7*dRp(ikx,iz,zsrclayfort));
                        dPdown(ikx,iz) = (term3*dRp(ikx,iz,zsrclayfort))+(term2*term13);
                    endif
                endif
            endif
        end do
    end do
elseif (above == 1) then ! receivers above the source layer: Pdownplus is not used
    do ikx=1,nk
        Pdown(ikx) = 0.0;
        do iz = 1,nlayers+1,1
            dPdown(ikx,iz) = 0.0
        end do
    end do
elseif (above == -1) then ! receivers below the source layer
    ! First compute P_{s+1}
    do ikx=1,nk
        if (zsrclayfort == 1) then ! If the source is in the top most layer
            Ms = 1.0;
            Pdown(ikx) = (1+Rp(ikx,1))*(exp(-1.0*Gam(ikx,zsrclayfort)*dp));
        else
            Ms = 1.0-Rm(ikx,1)*Rp(ikx,1)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
            Pdown(ikx) = (1+Rp(ikx,1))*(exp(-1.0*Gam(ikx,zsrclayfort)*dp)&
                -Rm(ikx,1)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm)));
        endif
        if (zsrclayfort+1 < nz) then ! if the source is in the last but one layer, Rp(:,:,2) is zero.
            Pdown(ikx) = Pdown(ikx)/(Ms*(1+Rp(ikx,2)&
                *exp(-2.0*Gam(ikx,zsrclayfort+1)*(z(zsrclayfort+2)-z(zsrclayfort+1)))));
        else
            Pdown(ikx) = Pdown(ikx)/Ms;
        end if
    end do
    if (nlayers > 2) then
        ! Second compute P for all other layers
        do iz = 3,nlayers
            do ikx=1,nk
                Pdown(ikx) = Pdown(ikx)*(1+Rp(ikx,iz-1))&
                    *exp(-1.0*Gam(ikx,zsrclayfort+iz-2)*(z(zsrclayfort+iz-1)-z(zsrclayfort+iz-2)));
                !if (zrcvlayfort+iz <= nz) then ! If the receiver is NOT in the last layer
                !if (zsrclayfort+iz <= nz) then ! If the receiver is NOT in the last layer
                if (zsrclayfort+iz-1 /= nz) then ! If the receiver is NOT in the last layer
                    Pdown(ikx) = Pdown(ikx)/(1+Rp(ikx,iz)&
                        *exp(-2.0*Gam(ikx,zsrclayfort+iz-1)*(z(zsrclayfort+iz)-z(zsrclayfort+iz-1))));
                endif
            end do
        end do
    endif
endif

END SUBROUTINE Pdownminderiv
