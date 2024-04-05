SUBROUTINE dWprop (Wup, Wdown, dWup, dWdown, Gam, dGam, nkx, cor, nz, z, zrcv, zrcvlay)
! Calculate the propagation of the EM-field from the layer boundaries to the receiver position
!
! Calling arguments:
! Wup:     Propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nkx
! Wdown:   Propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nkx
! dWup:    Derivative of Wup, complex array of size nkx
! dWdown:  Derivative of Wdown, complex array of size nkx
! Gam:     Vertical wavenumber Gamma, complex array of size nkx*nz
! dGam:    Derivative of Gam, complex array of size nkx*nz
! nkx:     Amount of wavenumbers, integer scalar
! cor:     Coordinates in the wavenumber domain, real array of size nkx
! nz:      Amount of layers, integer scalar
! z:       Depth of interfaces, real array of size nz
! zrcv:    Depth of receiver, real scalar
! zrcvlay: The number of the layer where the receivers are located, integer scalar

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nkx, nz, zrcvlay
COMPLEX,INTENT(INOUT) :: Wup(nkx), dWup(nkx), Wdown(nkx), dWdown(nkx)
COMPLEX,INTENT(IN)    :: Gam(nkx,nz), dGam(nkx,nz)
REAL   ,INTENT(IN)    :: cor(nkx), z(nz), zrcv
REAL                  :: dp, dm
INTEGER               :: ikx, zrcvlayfort

zrcvlayfort = zrcvlay + 1;
dp = z(zrcvlayfort+1)-zrcv;
dm = zrcv-z(zrcvlayfort);

do ikx=1,nkx
    Wup(ikx) = exp(-1.0*Gam(ikx,zrcvlayfort)*dp);
    Wdown(ikx) = exp(-1.0*Gam(ikx,zrcvlayfort)*dm);
    dWup(ikx) = -1.0*dp*exp(-1.0*Gam(ikx,zrcvlayfort)*dp)*dGam(ikx,zrcvlayfort);
    dWdown(ikx) = -1.0*dm*exp(-1.0*Gam(ikx,zrcvlayfort)*dm)*dGam(ikx,zrcvlayfort);
end do

END SUBROUTINE dWprop
