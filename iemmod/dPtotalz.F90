SUBROUTINE dPtotalzz (Ptot, dPtot, Pdown, dPdown, Pup, dPup, Wup, dWup, Wdown, dWdown, Gam, dGam, &
    Rp, Rm, etaH, etaV, nk, cor, nz, z, zrcv, zsrc, zsrclay, zrcvlay, gamA, nlayers, above, xdirect);
! Calculate total field for the 33 component (vertical electric receiver, vertical electric source)
!
! Calling arguments:
! Ptot:    Total field, complex array of size nk
! Pdown:   TM-mode part of the downgoing field, complex array of size nk
! Pup:     TM-mode part of the upgoing field, complex array of size nk
! Wup:     TM-mode propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Wdown:   TM-mode propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nk
! Gam:     Vertical wavenumber Gamma, complex array of size nk
! Rp:      TM-mode upgoing global reflection coefficient, complex array of size nk times nlayers
! Rm:      TM-mode downgoing global reflection coefficient, complex array of size nk times nlayers
! etaH:    Material parameter eta for the horizontal direction, complex array of size nz
! etaV:    Material parameter eta for the vertical direction, complex array of size nz
! nk:      Amount of wavenumbers, integer scalar
! cor:     Coordinates in the wavenumber domain, real array of size nk
! nz:      Amount of layers, integer scalar
! z:       Depth of interfaces, real array of size nz
! zrcv:    Depth of receivers, real scalar
! zsrc:    Depth of source, real scalar
! zsrclay: The number of the layer where the source is located, integer scalar
! zrcvlay: The number of the layer where the receivers are located, integer scala
! gamA:    Small gamma squared (zetaH*etaV), complex array of size nz
! nlayers: Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! above:   Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!          in the same layer (above=0), integer scalar
! xdirect: Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zsrclay, zrcvlay, nlayers, above, xdirect
COMPLEX,INTENT(INOUT) :: Ptot(nk), dPtot(nk*nz)
COMPLEX,INTENT(IN)    :: Gam(nk,nz), dGam(nk,nz), Pup(nk), dPup(nk,nz)
COMPLEX,INTENT(IN)    :: Pdown(nk), dPdown(nk,nz), Wup(nk), dWup(nk)
COMPLEX,INTENT(IN)    :: Rp(nk,nz), Rm(nk,nz), etaH(nz), etaV(nz), gamA(nz), Wdown(nk), dWdown(nk)
REAL   ,INTENT(IN)    :: cor(nk), z(nz), zsrc, zrcv
REAL                  :: d
INTEGER               :: ikx, iz, zsrclayfort, zrcvlayfort, signd
COMPLEX               :: factor, dfactor, temp1, temp2, Gin

zsrclayfort = zsrclay + 1; ! First array index in C is 0, but in Fortran it is 1!
zrcvlayfort = zrcvlay + 1;
d = zsrc-zrcv;
if (d < 0) d = -1.0*d;

if (above == 0) then ! receivers in the source layer
    do ikx=1,nk
        factor = 1.0/(2.0*Gam(ikx,zsrclayfort)*etaV(zsrclayfort));
        dfactor = (-1.0/(2.0*(Gam(ikx,zsrclayfort)**2.0)*etaV(zsrclayfort)))*dGam(ikx,zsrclayfort)-&
            1.0/(2.0*Gam(ikx,zsrclayfort)*(etaV(zsrclayfort)**2.0));
        if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
            Ptot(ikx) = factor*Pup(ikx)*Wup(ikx);
        elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
            Ptot(ikx) = factor*Pdown(ikx)*Wdown(ikx);
        else ! If the source and the receivers are in any layer in beteen 
            Ptot(ikx) = factor*(Pup(ikx)*Wup(ikx)+Pdown(ikx)*Wdown(ikx));
        endif
        if (xdirect /= 1) then ! The direct field is NOT computed in the space domain
            Gin = factor*exp(-1.0*Gam(ikx,zsrclayfort)*d);
            if (zrcv <= zsrc+0.001 .AND. zrcv >= zsrc-0.001) then !I guess zrcv and zsrc are never exactly the same
                Gin = Gin-1.0/etaV(zsrclayfort);
            endif
            Ptot(ikx) = Ptot(ikx)+Gin;
        endif
        do iz=1,nz
            if (xdirect ==1) then ! The direct field is computed in the space domain
                if (iz == zsrclayfort) then
                    factor = 1.0/(2.0*Gam(ikx,zsrclayfort)*etaV(zsrclayfort));
                    dfactor = (-1.0/(2.0*(Gam(ikx,zsrclayfort)**2.0)*etaV(zsrclayfort)))*dGam(ikx,zsrclayfort)-&
                        1.0/(2.0*Gam(ikx,zsrclayfort)*(etaV(zsrclayfort)**2.0));
                    if (zsrclayfort == 1) then !Conductivity wrt to top layer
                        temp1 = Pup(ikx)*Wup(ikx);
                        temp2 = (Wup(ikx)*dPup(ikx,iz)+Pup(ikx)*dWup(ikx));
                        dPtot(((iz-1)*nk)+ikx) = (temp1*dfactor)+(factor*temp2);
                    elseif (zsrclayfort == nz) then !Conductivity wrt to bottom layer
                        temp1 = (-1.0*Pdown(ikx)*Wdown(ikx));
                        temp2 = (-1.0*Wdown(ikx)*dPdown(ikx,iz)-Pdown(ikx)*dWdown(ikx));
                        dPtot(((iz-1)*nk)+ikx) = (temp1*dfactor)+(factor*temp2);
                    else !Conductivity wrt to other layers
                        temp1 = (Pup(ikx)*Wup(ikx)-Pdown(ikx)*Wdown(ikx));
                        temp2 = (Wup(ikx)*dPup(ikx,iz)+Pup(ikx)*dWup(ikx)&
                            -Wdown(ikx)*dPdown(ikx,iz)-Pdown(ikx)*dWdown(ikx));
                        dPtot(((iz-1)*nk)+ikx) = (temp1*dfactor)+(factor*temp2);
                    endif
                else
                    if (zsrclayfort == 1) then
                        temp1 = Pup(ikx)*Wup(ikx);
                        temp2 = (Wup(ikx)*dPup(ikx,iz));
                    elseif (zsrclayfort == nz) then
                        temp1 = (-1.0*Pdown(ikx)*Wdown(ikx));
                        temp2 = (-1.0*Wdown(ikx)*dPdown(ikx,iz));
                    else
                        temp1 = (Pup(ikx)*Wup(ikx)-Pdown(ikx)*Wdown(ikx));
                        temp2 = (Wup(ikx)*dPup(ikx,iz)-Wdown(ikx)*dPdown(ikx,iz));
                    endif
                    factor = 1.0/(2.0*Gam(ikx,zsrclayfort)*etaV(zsrclayfort));
                    dPtot(((iz-1)*nk)+ikx) = (factor*temp2);
                endif
            else
                
            endif
        end do
    end do
elseif (above == 1) then ! receivers above the source layer: Pdownplus is not used
    do ikx=1,nk
        factor = etaH(zsrclayfort)&
            /(2.0*etaV(zsrclayfort)*etaV(zrcvlayfort)*Gam(ikx,zsrclayfort));
        if (zrcvlayfort == 1) then ! If the receivers are in the top most layer
            Ptot(ikx) = factor*Pup(ikx)*Wup(ikx);
        else
            Ptot(ikx) = factor*Pup(ikx)*(Wup(ikx)+Rm(ikx,1)&
                *exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wdown(ikx));
        endif
    end do
elseif (above == -1) then ! receivers below the source layer
    do ikx=1,nk
        factor = etaH(zsrclayfort)&
            /(2.0*etaV(zsrclayfort)*etaV(zrcvlayfort)*Gam(ikx,zsrclayfort));
        if (zrcvlayfort == nz) then ! If the receivers are in the lower most layer
            Ptot(ikx) = factor*Pdown(ikx)*Wdown(ikx);
        else
            Ptot(ikx) = factor*Pdown(ikx)*(Wdown(ikx)+Rp(ikx,nlayers)&
                *exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wup(ikx));
        endif
    end do
endif

END SUBROUTINE dPtotalzz

SUBROUTINE dPtotalzx (Ptot, dPtot, Pdown, dPdown, Pup, dPup, Wup, dWup, Wdown, dWdown, Gam, &
    Rp, Rm, etaH, etaV, nk, cor, nz, z, zrcv, zsrc, zsrclay, zrcvlay, gamA, nlayers, above, xdirect);
! Calculate total field for the 31 component (vertical electric receiver, inline oriented electric source)
!
! Calling arguments:
! Ptot:    Total field, complex array of size nk
! Pdown:   TM-mode part of the downgoing field, complex array of size nk
! Pup:     TM-mode part of the upgoing field, complex array of size nk
! Wup:     TM-mode propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Wdown:   TM-mode propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nk
! Gam:     Vertical wavenumber Gamma, complex array of size nk
! Rp:      TM-mode upgoing global reflection coefficient, complex array of size nk times nlayers
! Rm:      TM-mode downgoing global reflection coefficient, complex array of size nk times nlayers
! etaH:    Material parameter eta for the horizontal direction, complex array of size nz
! etaV:    Material parameter eta for the vertical direction, complex array of size nz
! nk:      Amount of wavenumbers, integer scalar
! cor:     Coordinates in the wavenumber domain, real array of size nk
! nz:      Amount of layers, integer scalar
! z:       Depth of interfaces, real array of size nz
! zrcv:    Depth of receivers, real scalar
! zsrc:    Depth of source, real scalar
! zsrclay: The number of the layer where the source is located, integer scalar
! zrcvlay: The number of the layer where the receivers are located, integer scala
! gamA:    Small gamma squared (zetaH*etaV), complex array of size nz
! nlayers: Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! above:   Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!          in the same layer (above=0), integer scalar
! xdirect: Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zsrclay, zrcvlay, nlayers, above, xdirect
COMPLEX,INTENT(INOUT) :: Ptot(nk), dPtot(nk*nz)
COMPLEX,INTENT(IN)    :: Gam(nk,nz), Pup(nk), dPup(nk,nz), Pdown(nk), dPdown(nk,nz), Wup(nk), dWup(nk)
COMPLEX,INTENT(IN)    :: Rp(nk,nz), Rm(nk,nz), etaH(nz), etaV(nz), gamA(nz), Wdown(nk), dWdown(nk)
REAL   ,INTENT(IN)    :: cor(nk), z(nz), zsrc, zrcv
REAL                  :: d
INTEGER               :: ikx, iz, zsrclayfort, zrcvlayfort, signd
COMPLEX               :: factor, dfactor, temp1, temp2

zsrclayfort = zsrclay + 1; ! First array index in C is 0, but in Fortran it is 1!
zrcvlayfort = zrcvlay + 1;
d = zrcv-zsrc;
if (d < 0) then
    signd = -1.0;
else 
    signd = 1.0;
endif
if (d < 0) d = -1.0*d; ! replaces abs(d)

if (above == 0) then ! receivers in the source layer
    do ikx=1,nk
        factor = 1.0/(2.0*etaV(zsrclayfort));
        dfactor = -1.0/(2.0*(etaV(zsrclayfort)**2.0));
        if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
            if (xdirect == 1) then ! The direct field is computed in the space domain
                temp1 = Pup(ikx)*Wup(ikx);
                Ptot(ikx) = factor*temp1;
            else
                Ptot(ikx) = factor*(Pup(ikx)*Wup(ikx)+signd*exp(-1.0*Gam(ikx,zsrclayfort)*d));
            endif
        elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
            if (xdirect == 1) then ! The direct field is computed in the space domain
                temp1 = (-1.0*Pdown(ikx)*Wdown(ikx));
                Ptot(ikx) = factor*temp1;
            else
                Ptot(ikx) = factor*(-1.0*Pdown(ikx)*Wdown(ikx)+signd*exp(-1.0*Gam(ikx,zsrclayfort)*d));
            endif
        else ! If the source and the receivers are in any layer in between. 
            if (xdirect == 1) then ! The direct field is computed in the space domain
                temp1 = (Pup(ikx)*Wup(ikx)-Pdown(ikx)*Wdown(ikx));
                Ptot(ikx) = factor*temp1;
            else
                Ptot(ikx) = factor*(Pup(ikx)*Wup(ikx)-Pdown(ikx)*Wdown(ikx)&
                    &+signd*exp(-1.0*Gam(ikx,zsrclayfort)*d));
            endif
        endif
        do iz=1,nz
            if (xdirect ==1) then ! The direct field is computed in the space domain
                if (iz == zsrclayfort) then
                    factor = 1.0/(2.0*etaV(zsrclayfort));
                    if (zsrclayfort == 1) then !Conductivity wrt to top layer
                        temp1 = Pup(ikx)*Wup(ikx);
                        temp2 = (Wup(ikx)*dPup(ikx,iz)+Pup(ikx)*dWup(ikx));
                        dPtot(((iz-1)*nk)+ikx) = (temp1*dfactor)+(factor*temp2);
                    elseif (zsrclayfort == nz) then !Conductivity wrt to bottom layer
                        temp1 = (-1.0*Pdown(ikx)*Wdown(ikx));
                        temp2 = (-1.0*Wdown(ikx)*dPdown(ikx,iz)-Pdown(ikx)*dWdown(ikx));
                        dPtot(((iz-1)*nk)+ikx) = (temp1*dfactor)+(factor*temp2);
                    else !Conductivity wrt to other layers
                        temp1 = (Pup(ikx)*Wup(ikx)-Pdown(ikx)*Wdown(ikx));
                        temp2 = (Wup(ikx)*dPup(ikx,iz)+Pup(ikx)*dWup(ikx)&
                            -Wdown(ikx)*dPdown(ikx,iz)-Pdown(ikx)*dWdown(ikx));
                        dPtot(((iz-1)*nk)+ikx) = (temp1*dfactor)+(factor*temp2);
                    endif
                else
                    if (zsrclayfort == 1) then
                        temp1 = Pup(ikx)*Wup(ikx);
                        temp2 = (Wup(ikx)*dPup(ikx,iz));
                    elseif (zsrclayfort == nz) then
                        temp1 = (-1.0*Pdown(ikx)*Wdown(ikx));
                        temp2 = (-1.0*Wdown(ikx)*dPdown(ikx,iz));
                    else
                        temp1 = (Pup(ikx)*Wup(ikx)-Pdown(ikx)*Wdown(ikx));
                        temp2 = (Wup(ikx)*dPup(ikx,iz)-Wdown(ikx)*dPdown(ikx,iz));
                    endif
                    factor = 1.0/(2.0*etaV(zsrclayfort));
                    dPtot(((iz-1)*nk)+ikx) = (factor*temp2);
                endif
            else
                
            endif
        end do
    end do
elseif (above == 1) then ! receivers above the source layer: Pdownplus is not used
    do ikx=1,nk
        factor = 1.0/(2.0*etaV(zrcvlayfort));
        if (zrcvlayfort == 1) then ! If the receivers are in the top most layer
            Ptot(ikx) = factor*Pup(ikx)*Wup(ikx);
        else
            Ptot(ikx) = factor*Pup(ikx)*(Wup(ikx)+Rm(ikx,2)&
                *exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wdown(ikx));
        endif
    end do
elseif (above == -1) then ! receivers below the source layer
    do ikx=1,nk
        factor = 1.0/(2.0*etaV(zrcvlayfort));
        if (zrcvlayfort == nz) then ! If the receivers are in the lower most layer
            Ptot(ikx) = factor*Pdown(ikx)*Wdown(ikx);
        else
            Ptot(ikx) = factor*Pdown(ikx)*(Wdown(ikx)+Rp(ikx,nlayers)&
                *exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wup(ikx));
        endif
    end do
endif

END SUBROUTINE dPtotalzx

SUBROUTINE dPtotalzy (Ptot, dPtot, Pdown, dPdown, Pup, dPup, Wup, dWup, Wdown, dWdown, Gam, &
    Rp, Rm, etaH, etaV, nk, cor, nz, z, zrcv, zsrc, zsrclay, zrcvlay, gamA, nlayers, above, xdirect);
! Calculate total field for the 32 component (vertical electric receiver, crossline oriented electric source)
!
! Calling arguments:
! Ptot:    Total field, complex array of size nk
! Pdown:   TM-mode part of the downgoing field, complex array of size nk
! Pup:     TM-mode part of the upgoing field, complex array of size nk
! Wup:     TM-mode propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Wdown:   TM-mode propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nk
! Gam:     Vertical wavenumber Gamma, complex array of size nk
! Rp:      TM-mode upgoing global reflection coefficient, complex array of size nk times nlayers
! Rm:      TM-mode downgoing global reflection coefficient, complex array of size nk times nlayers
! etaH:    Material parameter eta for the horizontal direction, complex array of size nz
! etaV:    Material parameter eta for the vertical direction, complex array of size nz
! nk:      Amount of wavenumbers, integer scalar
! cor:     Coordinates in the wavenumber domain, real array of size nk
! nz:      Amount of layers, integer scalar
! z:       Depth of interfaces, real array of size nz
! zrcv:    Depth of receivers, real scalar
! zsrc:    Depth of source, real scalar
! zsrclay: The number of the layer where the source is located, integer scalar
! zrcvlay: The number of the layer where the receivers are located, integer scala
! gamA:    Small gamma squared (zetaH*etaV), complex array of size nz
! nlayers: Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! above:   Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!          in the same layer (above=0), integer scalar
! xdirect: Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zsrclay, zrcvlay, nlayers, above, xdirect
COMPLEX,INTENT(INOUT) :: Ptot(nk), dPtot(nk*nz)
COMPLEX,INTENT(IN)    :: Gam(nk,nz), Pup(nk), dPup(nk,nz), Pdown(nk), dPdown(nk,nz), Wup(nk), dWup(nk)
COMPLEX,INTENT(IN)    :: Rp(nk,nz), Rm(nk,nz), etaH(nz), etaV(nz), gamA(nz), Wdown(nk), dWdown(nk)
REAL   ,INTENT(IN)    :: cor(nk), z(nz), zsrc, zrcv
REAL                  :: d
INTEGER               :: ikx, iz, zsrclayfort, zrcvlayfort, signd
COMPLEX               :: factor, dfactor, temp1, temp2

zsrclayfort = zsrclay + 1; ! First array index in C is 0, but in Fortran it is 1!
zrcvlayfort = zrcvlay + 1;
d = zrcv-zsrc;
if (d < 0) then
    signd = -1.0;
else 
    signd = 1.0;
endif
if (d < 0) d = -1.0*d; ! replaces abs(d)

if (above == 0) then ! receivers in the source layer
    do ikx=1,nk
        factor = 1.0/(2.0*etaV(zsrclayfort));
        dfactor = -1.0/(2.0*(etaV(zsrclayfort)**2.0));
        if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
            if (xdirect == 1) then ! The direct field is computed in the space domain
                Ptot(ikx) = factor*(Pup(ikx)*Wup(ikx));
            else
                Ptot(ikx) = factor*(Pup(ikx)*Wup(ikx)+signd*exp(-1.0*Gam(ikx,zsrclayfort)*d));
            endif
        elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
            if (xdirect == 1) then ! The direct field is computed in the space domain
                Ptot(ikx) = factor*(-1.0*Pdown(ikx)*Wdown(ikx));
            else
                Ptot(ikx) = factor*(-1.0*Pdown(ikx)*Wdown(ikx)+signd*exp(-1.0*Gam(ikx,zsrclayfort)*d));
            endif
        else ! If the source and the receivers are in any layer in between. 
            if (xdirect == 1) then ! The direct field is computed in the space domain
                Ptot(ikx) = factor*(Pup(ikx)*Wup(ikx)-Pdown(ikx)*Wdown(ikx));
            else
                Ptot(ikx) = factor*(Pup(ikx)*Wup(ikx)-Pdown(ikx)*Wdown(ikx)&
                    &+signd*exp(-1.0*Gam(ikx,zsrclayfort)*d));
            endif
        endif
        do iz=1,nz
            if (xdirect ==1) then ! The direct field is computed in the space domain
                if (iz == zsrclayfort) then
                    factor = 1.0/(2.0*etaV(zsrclayfort));
                    if (zsrclayfort == 1) then !Conductivity wrt to top layer
                        temp1 = Pup(ikx)*Wup(ikx);
                        temp2 = (Wup(ikx)*dPup(ikx,iz)+Pup(ikx)*dWup(ikx));
                        dPtot(((iz-1)*nk)+ikx) = (temp1*dfactor)+(factor*temp2);
                    elseif (zsrclayfort == nz) then !Conductivity wrt to bottom layer
                        temp1 = (-1.0*Pdown(ikx)*Wdown(ikx));
                        temp2 = (-1.0*Wdown(ikx)*dPdown(ikx,iz)-Pdown(ikx)*dWdown(ikx));
                        dPtot(((iz-1)*nk)+ikx) = (temp1*dfactor)+(factor*temp2);
                    else !Conductivity wrt to other layers
                        temp1 = (Pup(ikx)*Wup(ikx)-Pdown(ikx)*Wdown(ikx));
                        temp2 = (Wup(ikx)*dPup(ikx,iz)+Pup(ikx)*dWup(ikx)&
                            -Wdown(ikx)*dPdown(ikx,iz)-Pdown(ikx)*dWdown(ikx));
                        dPtot(((iz-1)*nk)+ikx) = (temp1*dfactor)+(factor*temp2);
                    endif
                else
                    if (zsrclayfort == 1) then
                        temp1 = Pup(ikx)*Wup(ikx);
                        temp2 = (Wup(ikx)*dPup(ikx,iz));
                    elseif (zsrclayfort == nz) then
                        temp1 = (-1.0*Pdown(ikx)*Wdown(ikx));
                        temp2 = (-1.0*Wdown(ikx)*dPdown(ikx,iz));
                    else
                        temp1 = (Pup(ikx)*Wup(ikx)-Pdown(ikx)*Wdown(ikx));
                        temp2 = (Wup(ikx)*dPup(ikx,iz)-Wdown(ikx)*dPdown(ikx,iz));
                    endif
                    factor = 1.0/(2.0*etaV(zsrclayfort));
                    dPtot(((iz-1)*nk)+ikx) = (factor*temp2);
                endif
            else
                
            endif
        end do
    end do
elseif (above == 1) then ! receivers above the source layer: Pdownplus is not used
    do ikx=1,nk
        factor = 1.0/(2.0*etaV(zrcvlayfort));
        if (zrcvlayfort == 1) then ! If the receivers are in the top most layer
            Ptot(ikx) = factor*Pup(ikx)*Wup(ikx);
        else
            Ptot(ikx) = factor*Pup(ikx)*(Wup(ikx)+Rm(ikx,1)&
                *exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wdown(ikx));
        endif
    end do
elseif (above == -1) then ! receivers below the source layer
    do ikx=1,nk
        factor = 1.0/(2.0*etaV(zrcvlayfort));
        if (zrcvlayfort == nz) then ! If the receivers are in the lower most layer
            Ptot(ikx) = factor*Pdown(ikx)*Wdown(ikx);
        else
            Ptot(ikx) = factor*Pdown(ikx)*(Wdown(ikx)+Rp(ikx,nlayers)&
                *exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wup(ikx));
        endif
    end do
endif

END SUBROUTINE dPtotalzy
