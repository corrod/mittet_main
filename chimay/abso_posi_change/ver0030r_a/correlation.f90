!>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>> MAIN PROGRAM >>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>

program output_correlation
 implicit none
 integer, parameter :: n=3
 integer i, j
 real(8) x(n), y(n)
 real(8) xm, ym, sigx, sigy, cor
 real(8) a, b

 x(1) = 50000.0d0
 x(2) = 60000.0d0
 x(3) = 30000.0d0

 y(1) = 120.0d0
 y(2) = 130.0d0
 y(3) =  90.0d0

 call correlation(x, y, n, xm, ym, sigx, sigy, cor)
 write(*,*) 'Correlation = ', cor

end program output_correlation

!>>>>>>>>>>>>>>>>>>>>>>
!>>>>> SUBROUTINE >>>>>
!>>>>>>>>>>>>>>>>>>>>>>

subroutine correlation(x, y, n, xm, ym, sigx, sigy, cor)
 implicit none
 integer n, i
 real(8) x(n), y(n)
 real(8) xm, ym, sigx, sigy, cor
 real(8) sx, sy, sxy, sx2, sy2

 sx=0.0d0; sy=0.0d0; sxy=0.0d0; sx2=0.0d0; sy2=0.0d0

 if( n > 2 ) then
  do i=1, n
   sx = sx +x(i)
   sy = sy +y(i)
  end do
  xm = sx/n; ym = sy/n
  do i=1, n
   sx2 = sx2 + (x(i)-xm)*(x(i)-xm)
   sy2 = sy2 + (y(i)-ym)*(y(i)-ym)
   sxy = sxy + (x(i)-xm)*(y(i)-ym)
  end do

  sigx = sqrt(sx2/n); sigy = sqrt(sy2/n)
  if( sigx > 0 .AND. sigy > 0 ) then
  cor = (sxy/n)/(sigx*sigy)
  end if
 end if
end subroutine correlation



subroutine rearrange_Rdata

end subroutine rearrange_Rdata


subroutine rearrange_modeldata

end subroutine rearrange_modeldata

! !>>>>>>>>>>>>>>>>>>>>>>>>
! !>>>>> MAIN PROGRAM >>>>>
! !>>>>>>>>>>>>>>>>>>>>>>>>

! program output_correlation
!  implicit none
!  integer, parameter :: n=3
!  integer i, j
!  real(8) x(n), y(n)
!  real(8) xm, ym, sigx, sigy, cor
!  real(8) a, b

!  x(1) = 50000.0d0
!  x(2) = 60000.0d0
!  x(3) = 30000.0d0

!  y(1) = 120.0d0
!  y(2) = 130.0d0
!  y(3) =  90.0d0

!  call correlation(x, y, n, xm, ym, sigx, sigy, cor)
!  write(*,*) 'Correlation = ', cor

! end program output_correlation

! !>>>>>>>>>>>>>>>>>>>>>>
! !>>>>> SUBROUTINE >>>>>
! !>>>>>>>>>>>>>>>>>>>>>>

! subroutine correlation(x, y, n, xm, ym, sigx, sigy, cor)
!  implicit none
!  integer n, i
!  real(8) x(n), y(n)
!  real(8) xm, ym, sigx, sigy, cor
!  real(8) sx, sy, sxy, sx2, sy2

!  sx=0.0d0; sy=0.0d0; sxy=0.0d0; sx2=0.0d0; sy2=0.0d0

!  if( n > 2 ) then
!   do i=1, n
!    sx = sx +x(i)
!    sy = sy +y(i)
!   end do
!   xm = sx/n; ym = sy/n
!   do i=1, n
!    sx2 = sx2 + (x(i)-xm)*(x(i)-xm)
!    sy2 = sy2 + (y(i)-ym)*(y(i)-ym)
!    sxy = sxy + (x(i)-xm)*(y(i)-ym)
!   end do

!   sigx = sqrt(sx2/n); sigy = sqrt(sy2/n)
!   if( sigx > 0 .AND. sigy > 0 ) then
!   cor = (sxy/n)/(sigx*sigy)
!   end if
!  end if
! end subroutine correlation



! subroutine rearrange_Rdata

! end subroutine rearrange_Rdata


! subroutine rearrange_modeldata

! end subroutine rearrange_modeldata
