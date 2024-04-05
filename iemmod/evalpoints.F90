SUBROUTINE evalpoints_lin (xposnew,xtotrad1new,xtotrad2new,markernew,nlogxnew,nlogxdata,newit,xpos,xtotrad1,xtotrad2,nlogx,&
    c1,c2)
! Evaluate which points are necessary to achieve the precision defined with the parameters c1 and c2 
! for linear interpolation and Hankel transformations with two terms.
!
! Calling arguments:
! xposnew:     After completion of the function, this array contains the coordinates in the space domain 
!              for all datapoints for which xtotrad1 and xtotrad2 have been computed 
!              as well as all the possible new datapoints, real array of size nlogxnew 
! xtotrad1new: After completion of this function, this array contains the values of xtotrad1 (Hankel transformed EM-field) 
!              and zeros on the location where the Hankel transformation still needs to be computed, complex array of size nlogxnew
! xtotrad2new: After completion of this function, this array contains the values of xtotrad2 (Hankel transformed EM-field) 
!              and zeros on the location where the Hankel transformation still needs to be computed, complex array of size nlogxnew
! markernew:   After completion of this function, this array hold the indication for which locations in space 
!              the Hankel transformation still needs to be computed, integer array of size nlogxnew
! nlogxnew:    Maximum amount of values in the space domain that are possible. which is  
!              nlogx points of data plus nlogx-1 points that could be added, integer scalar
!              Note: In the worst case, between every point a new datapoint is required.
! nlogxdata:   After completion of this function, this parameter contains the amount of actual datapoints, which is
!              nlogx points of data plus the amount of points that were actually added, integer scalar
!              Note: This value can not be larger than nlogxnew.
! newit:       After completion of this function, this parameter indicates if a new iteration is required, i.e., 
!              if new datapoints have to be added, integer scalar
! xpos:        Logarithmic coordinate vector in the space domain, real array of size nlogx
! xtotrad1:    1st Hankel transformed field-term, complex array of size nlogx
! xtotrad2:    2nd Hankel transformed field-term, complex array of size nlogx
! nlogx:       Amount of samples in logarithmic coordinate vector, integer scalar
! c1:          1st precision parameter, real scalar
! c2:          2nd precision parameter, real scalar

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: xtotrad1new(nlogxnew), xtotrad2new(nlogxnew)
REAL,INTENT(INOUT)    :: xposnew(nlogxnew)
INTEGER,INTENT(INOUT) :: markernew(nlogxnew), newit, nlogxdata
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx), xtotrad2(nlogx)
REAL   ,INTENT(IN)    :: xpos(nlogx), c1, c2
INTEGER,INTENT(IN)    :: nlogx, nlogxnew
COMPLEX               :: test1, test2, lobound, hibound
REAL                  :: thisdx, weightlo, weighthi, lim
INTEGER, DIMENSION(:), ALLOCATABLE :: marker
INTEGER               :: nlogxfix, nlogxint, ix, ixs, lopos, hipos, ierr, next
LOGICAL               :: domark

! mark elements of xpos with 1 that are kept and  with 0 that are tested
domark = .true.;
allocate(marker(1:nlogx))
do ix=1,nlogx
    if (domark .eqv. .true.) then
        marker(ix) = 1; 
    else
        marker(ix) = 0;
    endif
    domark = .not.domark; 
end do
! the last element of xpos should always be kept
marker(nlogx) = 1;
! count the elements that are kept
nlogxfix = sum(marker);
! compute the elements that need to be tested
nlogxint = nlogx-nlogxfix;
! interpolate and test the points in question
do ix=1,nlogx
    ! If the marker for a specific point is 0, it needs to be tested
    if (marker(ix) == 0) then
        ! The value of the field at this point is computed
        ! based on a linear interpolation between the two adjacents points.
        ! The interpolated value is stored in the variable test1
        lopos = ix-1;
        hipos = ix+1;
        thisdx = xpos(hipos) - xpos(lopos);
        weightlo = (xpos(hipos)-xpos(ix))/thisdx;
        weighthi = (xpos(ix)-xpos(lopos))/thisdx;
        test1 = weightlo*xtotrad1(lopos)+weighthi*xtotrad1(hipos);
        ! Based on the correct value at that position and the precision defined with parameters c1 and c2
        ! the confindence area is computed, which has a lowerbound lobound and an upper bound hibound.
        lim = c1*log10(abs(real(xtotrad1(ix))))+c2;
        lobound = log10(abs(real(xtotrad1(ix))))-lim*abs(log10(abs(real(xtotrad1(ix)))));
        hibound = log10(abs(real(xtotrad1(ix))))+lim*abs(log10(abs(real(xtotrad1(ix)))));
        ! Test if the interpolated value falls within the confidence area
        if (real(log10(abs(real(test1))))>=real(lobound) .AND. real(log10(abs(real(test1))))<=real(hibound)) then
            ! Yes, there precision is sufficient at this point.
            marker(ix) = 1;
        else 
            ! No, more points need to be added at this point.
            marker(ix) = 2;
        endif
        ! The test above was only for the real part. If the real part was considered precise enough, 
        ! the same test is also carried out for the imaginary part. If the real part was not precise enough,
        ! and it was decided to add more points, there is no need to test the imaginary part as well.
        if (marker(ix) == 1) then
            lim = c1*log10(abs(aimag(xtotrad1(ix))))+c2;
            lobound = log10(abs(aimag(xtotrad1(ix))))-lim*abs(log10(abs(aimag(xtotrad1(ix)))));
            hibound = log10(abs(aimag(xtotrad1(ix))))+lim*abs(log10(abs(aimag(xtotrad1(ix)))));
            if (real(log10(abs(aimag(test1))))>=real(lobound) .AND. real(log10(abs(aimag(test1))))<=real(hibound)) then
                marker(ix) = 1;
            else 
                marker(ix) = 2;
            endif
        endif
        ! The field consists of two different parts, xtotrad1 and xtotrad2. If xtotrad1 was precise enough and no new point
        ! was added, it is tested if more points are necessary for xtotrad2. This is done first for the real part and secondly
        ! for the imaginary part.
        if (marker(ix) == 1) then
            test2 = weightlo*xtotrad2(lopos)+weighthi*xtotrad2(hipos);
            lim = c1*log10(abs(real(xtotrad2(ix))))+c2;
            lobound = log10(abs(real(xtotrad2(ix))))-lim*abs(log10(abs(real(xtotrad2(ix)))));
            hibound = log10(abs(real(xtotrad2(ix))))+lim*abs(log10(abs(real(xtotrad2(ix)))));
            if (real(log10(abs(real(test2))))>=real(lobound) .AND. real(log10(abs(real(test2))))<=real(hibound)) then
                marker(ix) = 1;
            else 
                marker(ix) = 2;
            endif
        endif
        if (marker(ix) == 1) then
            lim = c1*log10(abs(aimag(xtotrad2(ix))))+c2;
            lobound = log10(abs(aimag(xtotrad2(ix))))-lim*abs(log10(abs(aimag(xtotrad2(ix)))));
            hibound = log10(abs(aimag(xtotrad2(ix))))+lim*abs(log10(abs(aimag(xtotrad2(ix)))));
            if (real(log10(abs(aimag(test2))))>=real(lobound) .AND. real(log10(abs(aimag(test2))))<=real(hibound)) then
                marker(ix) = 1;
            else 
                marker(ix) = 2;
            endif
        endif
    endif
end do
! Set up the new coordinate and data vector
ixs = 1;
newit = 0;
do ix = 1,nlogx
    if (marker(ix) == 1) then
        ! The precision at this point was high enough. 
        ! Thus, the point is simply copied to the new coordinate vector.
        markernew(ixs) = 1;
        xposnew(ixs) = xpos(ix);
        xtotrad1new(ixs) = xtotrad1(ix);
        xtotrad2new(ixs) = xtotrad2(ix);
        ixs = ixs + 1;
    else if (marker(ix) == 2) then
        ! The precision at this point was not high enough.
        ! A new point is introduced before this point.
        newit = 1;
        markernew(ixs) = 2;
        xposnew(ixs) = (xpos(ix)-xpos(ix-1))/2+xpos(ix-1);
        xtotrad1new(ixs) = 0.0;
        xtotrad2new(ixs) = 0.0;
        ixs = ixs + 1;
        ! The point itself is also copied to the new coordinate vector.
        markernew(ixs) = 1;
        xposnew(ixs) = xpos(ix);
        xtotrad1new(ixs) = xtotrad1(ix);
        xtotrad2new(ixs) = xtotrad2(ix);
        ixs = ixs + 1;
        ! A new point is introduced after this point. 
        markernew(ixs) = 2;
        xposnew(ixs) = (xpos(ix+1)-xpos(ix))/2+xpos(ix);
        xtotrad1new(ixs) = 0.0;
        xtotrad2new(ixs) = 0.0;
        ixs = ixs + 1;
    endif
end do
nlogxdata = ixs-1;
deallocate(marker)
END SUBROUTINE evalpoints_lin

SUBROUTINE evalpoints_lin_mono (xposnew,xtotrad1new,markernew,nlogxnew,nlogxdata,newit,xpos,xtotrad1,nlogx,&
    c1,c2)
! Evaluate which points are necessary to achieve the precision defined with the parameters c1 and c2 
! for linear interpolation and Hankel transformations with one term. 
!
! For detailed comments see function evalpoints_lin (further above in this file).
!
! Calling arguments:
! xposnew:     After completion of the function, this array contains the coordinates in the space domain 
!              for all datapoints for which xtotrad1 and xtotrad2 have been computed 
!              as well as all the possible new datapoints, real array of size nlogxnew 
! xtotrad1new: After completion of this function, this array contains the values of xtotrad1 (Hankel transformed EM-field) 
!              and zeros on the location where the Hankel transformation still needs to be computed, complex array of size nlogxnew
! markernew:   After completion of this function, this array hold the indication for which locations in space 
!              the Hankel transformation still needs to be computed, integer array of size nlogxnew
! nlogxnew:    Maximum amount of values in the space domain that are possible. which is  
!              nlogx points of data plus nlogx-1 points that could be added, integer scalar
!              Note: In the worst case, between every point a new datapoint is required.
! nlogxdata:   After completion of this function, this parameter contains the amount of actual datapoints, which is
!              nlogx points of data plus the amount of points that were actually added, integer scalar
!              Note: This value can not be larger than nlogxnew.
! newit:       After completion of this function, this parameter indicates if a new iteration is required, i.e., 
!              if new datapoints have to be added, integer scalar
! xpos:        Logarithmic coordinate vector in the space domain, real array of size nlogx
! xtotrad1:    1st Hankel transformed field-term, complex array of size nlogx
! nlogx:       Amount of samples in logarithmic coordinate vector, integer scalar
! c1:          1st precision parameter, real scalar
! c2:          2nd precision parameter, real scalar

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: xtotrad1new(nlogxnew)
REAL,INTENT(INOUT)    :: xposnew(nlogxnew)
INTEGER,INTENT(INOUT) :: markernew(nlogxnew), newit, nlogxdata
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx)
REAL   ,INTENT(IN)    :: xpos(nlogx), c1, c2
INTEGER,INTENT(IN)    :: nlogx, nlogxnew
COMPLEX               :: test1, test2, lobound, hibound
REAL                  :: thisdx, weightlo, weighthi, lim
INTEGER, DIMENSION(:), ALLOCATABLE :: marker
INTEGER               :: nlogxfix, nlogxint, ix, ixs, lopos, hipos, ierr, next
LOGICAL               :: domark

! mark elements of xpos with 1 that are kept and  with 0 that are tested
domark = .true.;
allocate(marker(1:nlogx))
do ix=1,nlogx
    if (domark .eqv. .true.) then
        marker(ix) = 1; 
    else
        marker(ix) = 0;
    endif
    domark = .not.domark; 
end do
! the last element of xpos should always be kept
marker(nlogx) = 1;
! count the elements that are kept
nlogxfix = sum(marker);
! compute the elements that need to be tested
nlogxint = nlogx-nlogxfix;
! interpolate and test the points in question
do ix=1,nlogx
    if (marker(ix) == 0) then
        lopos = ix-1;
        hipos = ix+1;
        thisdx = xpos(hipos) - xpos(lopos);
        weightlo = (xpos(hipos)-xpos(ix))/thisdx;
        weighthi = (xpos(ix)-xpos(lopos))/thisdx;
        test1 = weightlo*xtotrad1(lopos)+weighthi*xtotrad1(hipos);
        lim = c1*log10(abs(real(xtotrad1(ix))))+c2;
        lobound = log10(abs(real(xtotrad1(ix))))-lim*abs(log10(abs(real(xtotrad1(ix)))));
        hibound = log10(abs(real(xtotrad1(ix))))+lim*abs(log10(abs(real(xtotrad1(ix)))));
        if (real(log10(abs(real(test1))))>=real(lobound) .AND. real(log10(abs(real(test1))))<=real(hibound)) then
            marker(ix) = 1;
        else
            marker(ix) = 2;
        endif
        ! Testing of the imaginary part is only required if no datapoints need to be added for the real part.
        if (marker(ix) == 1) then
            lim = c1*log10(abs(aimag(xtotrad1(ix))))+c2;
            lobound = log10(abs(aimag(xtotrad1(ix))))-lim*abs(log10(abs(aimag(xtotrad1(ix)))));
            hibound = log10(abs(aimag(xtotrad1(ix))))+lim*abs(log10(abs(aimag(xtotrad1(ix)))));
            if (real(log10(abs(aimag(test1))))>=real(lobound) .AND. real(log10(abs(aimag(test1))))<=real(hibound)) then
                marker(ix) = 1;
            else 
                marker(ix) = 2;
            endif
        endif
    endif
end do
! Set up the new coordinate and data vector
ixs = 1;
newit = 0;
do ix = 1,nlogx
    if (marker(ix) == 1) then
        markernew(ixs) = 1;
        xposnew(ixs) = xpos(ix);
        xtotrad1new(ixs) = xtotrad1(ix);
        ixs = ixs + 1;
    else if (marker(ix) == 2) then
        newit = 1;
        markernew(ixs) = 2;
        xposnew(ixs) = (xpos(ix)-xpos(ix-1))/2+xpos(ix-1);
        xtotrad1new(ixs) = 0.0;
        ixs = ixs + 1;
        markernew(ixs) = 1;
        xposnew(ixs) = xpos(ix);
        xtotrad1new(ixs) = xtotrad1(ix);
        ixs = ixs + 1;
        markernew(ixs) = 2;
        xposnew(ixs) = (xpos(ix+1)-xpos(ix))/2+xpos(ix);
        xtotrad1new(ixs) = 0.0;
        ixs = ixs + 1;
    endif
end do
nlogxdata = ixs-1;
deallocate(marker)
END SUBROUTINE evalpoints_lin_mono
