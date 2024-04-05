SUBROUTINE Pupplusderiv (Pup, dPup, Rp, Rm, dRp, dRm, Gam, dGam, nk, cor, nz, z, zsrc, zrcvlay, zsrclay, nlayers, above);
! Calculate Pup+
!
! Calling arguments:
! Pup:     Upgoing field, complex array of size nk
! dPup:    Derivative of upgoing field, complex array of size nk times (nlayers+1)
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
COMPLEX,INTENT(INOUT) :: Pup(nk), Rp(nk,nz), Rm(nk,nz), dPup(nk,nz), dRp(nk,nz,nz)
COMPLEX,INTENT(INOUT) :: dRm(nk,nz,nz)
COMPLEX,INTENT(IN)    :: Gam(nk,nz), dGam(nk,nz)
REAL   ,INTENT(IN)    :: cor(nk), z(nz)
REAL                  :: zsrc, ds, dp, dm
INTEGER               :: iz, ikx, zsrclayfort, zrcvlayfort
COMPLEX               :: Ms, term1, term2, term3, term4, term5, term6, term7, term8, term9, term10
COMPLEX               :: term11, term12, term13

zsrclayfort = zsrclay + 1; ! First array index in C is 0, but in Fortran it is 1!
zrcvlayfort = zrcvlay + 1;
ds = z(zsrclayfort+1)-z(zsrclayfort);
dp = z(zsrclayfort+1)-zsrc;
dm = zsrc-z(zsrclayfort);
if (above == 0) then ! receivers in the source layer
    do ikx=1,nk
        if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
            Ms = 1.0;
            Pup(ikx) = Rp(ikx,zsrclayfort)/Ms*(exp(-1.0*Gam(ikx,zsrclayfort)*dp));
        elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
            Pup(ikx) = 0.0;
        else ! If the source and the receivers are in any layer in between
            Ms = 1.0-Rm(ikx,zsrclayfort)*Rp(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
            Pup(ikx) = Rp(ikx,zsrclayfort)/Ms*(exp(-1.0*Gam(ikx,zsrclayfort)*dp)&
                +Rm(ikx,zsrclayfort)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm)));
        endif
        do iz=1,nz
            if (iz == zsrclayfort) then
                if (zsrclayfort == nz) then
                    dPup(ikx,iz) = 0.0;
                elseif (zsrclayfort == 1) then
                    Ms = 1.0;
                    term1  = (1.0/Ms)*(exp(-1.0*Gam(ikx,zsrclayfort)*dp));
                    term2  = ((-1.0*Rp(ikx,zsrclayfort))/(Ms*Ms))*(exp(-1.0*Gam(ikx,zsrclayfort)*dp));
                    term3  = ((1.0*Rp(ikx,zsrclayfort))/Ms)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm));
                    term4  = Rp(ikx,zsrclayfort)/Ms;
                    term6  = -Rp(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
                    term9  = -2.0*ds*exp(-2.0*Gam(ikx,zsrclayfort)*ds)*dGam(ikx,zsrclayfort);
                    term10 = -1.0*dp*exp(-1.0*Gam(ikx,zsrclayfort)*dp)*dGam(ikx,zsrclayfort);
                    term11 = -1.0*(ds+dm)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm))*dGam(ikx,zsrclayfort);
                    dPup(ikx,iz) = (term1*dRp(ikx,iz,iz))+(term4*term10);
                else ! If the source and the receivers are in any layer in between or the top layer
                    Ms = 1.0-Rm(ikx,zsrclayfort)*Rp(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
                    term1  = (1.0/Ms)*(exp(-1.0*Gam(ikx,zsrclayfort)*dp)&
                        +Rm(ikx,zsrclayfort)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm)));
                    term2  = ((-1.0*Rp(ikx,zsrclayfort))/(Ms*Ms))*(exp(-1.0*Gam(ikx,zsrclayfort)*dp)&
                        +Rm(ikx,zsrclayfort)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm)));
                    term3  = ((1.0*Rp(ikx,zsrclayfort))/Ms)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm));
                    term4  = Rp(ikx,zsrclayfort)/Ms;
                    term5  = 1.0*((Rp(ikx,zsrclayfort)*Rm(ikx,zsrclayfort))/Ms);
                    term6  = -Rp(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
                    term7  = -Rm(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
                    term8  = -1.0*Rm(ikx,zsrclayfort)*Rp(ikx,zsrclayfort);
                    term9  = -2.0*ds*exp(-2.0*Gam(ikx,zsrclayfort)*ds)*dGam(ikx,zsrclayfort);
                    term10 = -1.0*dp*exp(-1.0*Gam(ikx,zsrclayfort)*dp)*dGam(ikx,zsrclayfort);
                    term11 = -1.0*(ds+dm)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm))*dGam(ikx,zsrclayfort);
                    term12 = (term6*dRm(ikx,zsrclayfort,zsrclayfort))&
                        +(term7*dRp(ikx,zsrclayfort,zsrclayfort))+(term8*term9);
                    dPup(ikx,iz) = (term1*dRp(ikx,iz,iz))+(term3*dRm(ikx,iz,iz))+(term2*term12)&
                        +(term4*term10)+(term5*term11);
                endif
            else
                if (iz < zsrclayfort) then ! No contribution by dRp
                    if (zsrclayfort == 1) then
                        dPup(ikx,iz) = 0.0;
                    elseif (zsrclayfort == nz) then
                        dPup(ikx,iz) = 0.0;
                    else
                        Ms = 1.0-Rm(ikx,zsrclayfort)*Rp(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
                        term2  = ((-1.0*Rp(ikx,zsrclayfort))/(Ms*Ms))*(exp(-1.0*Gam(ikx,zsrclayfort)*dp)&
                            +Rm(ikx,zsrclayfort)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm)));
                        term3  = ((1.0*Rp(ikx,zsrclayfort))/Ms)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm));
                        term6  = -Rp(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
                        term13 = (term6*dRm(ikx,iz,zsrclayfort));
                        dPup(ikx,iz) = (term3*dRm(ikx,iz,zsrclayfort))+(term2*term13);
                    endif
                elseif (iz > zsrclayfort) then !No contribution by dRm
                    if (zsrclayfort == 1) then
                        Ms = 1.0;
                        term1  = (1.0/Ms)*(exp(-1.0*Gam(ikx,zsrclayfort)*dp));
                        dPup(ikx,iz) = (term1*dRp(ikx,iz,zsrclayfort));
                    elseif (zsrclayfort == nz) then
                        dPup(ikx,iz) = 0.0;
                    else
                        Ms = 1.0-Rm(ikx,zsrclayfort)*Rp(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
                        term1  = (1.0/Ms)*(exp(-1.0*Gam(ikx,zsrclayfort)*dp)&
                            +Rm(ikx,zsrclayfort)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm)));
                        term2  = ((-1.0*Rp(ikx,zsrclayfort))/(Ms*Ms))*(exp(-1.0*Gam(ikx,zsrclayfort)*dp)&
                            +Rm(ikx,zsrclayfort)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm)));
                        term7  = -Rm(ikx,zsrclayfort)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
                        term13 = (term7*dRp(ikx,iz,zsrclayfort));
                        dPup(ikx,iz) = (term1*dRp(ikx,iz,zsrclayfort))+(term2*term13);
                    endif
                endif
            endif
        end do
        ! Note that if the source and the receivers are in the top most layer, Rm is zero. 
        ! The if statement is introduced to avoid the exponential function to produce values
        ! like infinity. Even if Rm is zero, a product with infinity causes problems. 
    end do
elseif (above == 1) then ! receivers above the source layer
    ! First compute P_{s-1}
    do ikx=1,nk
        if (zsrclayfort == nz) then ! If the source is in the lower most layer
            Ms = 1.0;
            Pup(ikx) = (1+Rm(ikx,nlayers+1))*(exp(-1.0*Gam(ikx,zsrclayfort)*dm));
        else
            Ms = 1.0-Rm(ikx,nlayers+1)*Rp(ikx,nlayers)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
            Pup(ikx) = (1+Rm(ikx,nlayers+1))*(exp(-1.0*Gam(ikx,zsrclayfort)*dm)&
                +Rp(ikx,nlayers)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dp)));
        endif
        Pup(ikx) = Pup(ikx)/(Ms*(1+Rm(ikx,nlayers)&
            *exp(-2.0*Gam(ikx,zsrclayfort-1)*(z(zsrclayfort)-z(zsrclayfort-1)))));
    end do
    if (nlayers > 2) then
        ! Second compute P for all other layers
        do iz = 1,nlayers-2
            do ikx=1,nk
                Pup(ikx) = Pup(ikx)*(1+Rm(ikx,iz+2))&
                    *exp(-1.0*Gam(ikx,zrcvlayfort+iz)*(z(zrcvlayfort+iz+1)-z(zrcvlayfort+iz)));
                if (zrcvlayfort+iz-1 /= 1) then ! If the receivers are NOT in the first layer
                    Pup(ikx) = Pup(ikx)/(1+Rm(ikx,iz+1)&
                        *exp(-2.0*Gam(ikx,zrcvlayfort+iz-1)*(z(zrcvlayfort+iz)-z(zrcvlayfort+iz-1)))); 
                    ! Note that if the receivers are in the first layer, Rm is zero and this line is just 
                    ! a division by 1. The if statement is introduced to avoid the exponential function
                    ! to cause a problem, because z(zrcvlayfort+iz-1) has only a dummy value.
                endif
            end do
        end do
    endif
elseif (above == -1) then ! receivers below the source layer: Pupplus is not used
    do ikx=1,nk
        Pup(ikx) = 0.0;
        do iz = 1,nlayers+1,1
            dPup(ikx,iz) = 0.0;
        end do
    end do
endif

END SUBROUTINE Pupplusderiv
