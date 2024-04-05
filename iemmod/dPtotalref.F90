SUBROUTINE dPtotalref (Ptot, dPtot, Wup, dWup, Rp, dRp, nk, nz, zrcvlay);
! calculate the reflection response at the receiver level
!
! Calling arguments:
! Ptot: Total reflection response, complex array of size nk
! Wup:  Propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Rp:   Upgoing global reflection coefficient, complex array of size nk
! nk:   Amount of wavenumbers, integer scalar

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zrcvlay
COMPLEX,INTENT(INOUT) :: Ptot(nk), dPtot(nk*nz)
COMPLEX,INTENT(IN)    :: Wup(nk), dWup(nk)
COMPLEX,INTENT(IN)    :: Rp(nk*nz), dRp(nk,nz,nz)
INTEGER               :: ikx, iz, zrcvlayfort

zrcvlayfort = zrcvlay + 1;
do ikx=1,nk
    Ptot(ikx) = Rp((zrcvlayfort-1)*nk+ikx)*Wup(ikx);
    do iz=1,nz
        if (iz == zrcvlayfort) then
            dPtot((iz-1)*nk+ikx) = dRp(ikx,iz,zrcvlayfort)*Wup(ikx) + Rp((zrcvlayfort-1)*nk+ikx)*dWup(ikx);
        else
            dPtot((iz-1)*nk+ikx) = dRp(ikx,iz,zrcvlayfort)*Wup(ikx);
        endif
    end do
end do

END SUBROUTINE dPtotalref
