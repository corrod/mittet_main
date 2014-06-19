module fft_lib
!============================================================================
   !******** physical constant ******** 
   complex(kind(0d0)),parameter :: complex_j=(0.0d0,1.0d0) ! imaginary number
   real(8),parameter :: pi=3.1415927d0          ! circular constant
   !******** time data ********
   real(8),dimension(:),allocatable :: h        ! input time domain data
   !******** frequency spectrum data ********
   complex(kind(0d0)),dimension(:),allocatable :: x,a,b        ! output spectrum
   !******** fft parameters ********
   integer :: s       ! total sampling numbers
   integer :: gamma   ! s=2**gamma -> gamma=log10(s)/log10(2) 
   real(8) :: t       ! delta time (t0/s) [sec]
   real(8) :: t0      ! assumed cycle [sec]
   !******** input file ********
   integer :: code          ! 
   integer :: count=0       !
   integer :: idm     ! dummy data in integer type
   real(8) :: rdm,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15     ! dummy data in real type
   real(8) :: col1    ! data under consideration
   character(90) file1      !
!============================================================================
end module fft_lib  

program fft_ver3
!***********************************************************************     
!   Fast Fourier Transform based on the 2**gamma  version 1.00   
!   coded by Yusuke Kusama   <kusama@dt.takuma-ct.ac.jp>
!
!   fisrt date    : 2003/10/06 !
!   last modified : 2005/08/31 !
!***********************************************************************
   use fft_lib
   implicit none
   integer :: n,k !  
   integer :: ts,te,t_rate,t_max,diff
  call system_clock(ts)

   call initial
   call read_data
   call CDFT_ftod
   call output

  call system_clock(te,t_rate,t_max)
  if (te<ts) then
     diff=t_max-ts+te
  else
     diff=te-ts
  end if
  print "(A,F10.3)","time it took was:", diff/dble(t_rate)
end program fft_ver3

!*************************************************************************
! initial data
!*************************************************************************
subroutine initial 
   use fft_lib
   implicit none
   integer :: n !   

   !******** read input data information ********
   file1= 'func_gauss1f.dat'
   !==== sampling numbers ====
   print *, "enter total sampling numbers N"
   s=16384
!   s=65536  !2**16
!   s=131072  !2**17
!   s=262144
!   s=10000   ! 2**13
   !==== sampling interval T ====
   print *, "enter sampling width T"
   t=1.14315354400754132d-3
   gamma=nint(log10(real(s))/log10(real(2)))
   t0=s*t

   !******** memory allocation ********
   allocate(h(0:s-1),x(0:s-1),a(0:s-1),b(0:s-1),stat=code)
   if(code/=0) stop "*** not enough memory <fft_v1> ***"
   !******** array initialize ********
   do n=0,s-1
      h(n)=0.0d0
      x(n)=(0.0d0,0.0d0) 
      a(n)=(0.0d0,0.0d0)
      b(n)=(0.0d0,0.0d0)
   end do

   !******** initial data confirmation ********
   write(50,*)
   write(50,'(1x, "sampling number     :",i6)') s   
   write(50,'(1x, "gamma               :",i6)') gamma
   write(50,'(1x, "sampling interval   :",d18.10)') t
   write(50,'(1x, "assumed period t0   :",d18.10)') t0
   write(50,'(1x, "delta freq          :",d18.10)') 1.0d0/t0

end subroutine initial

!*************************************************************************
! read data
!*************************************************************************
subroutine read_data 
   use fft_lib
   implicit none
   integer :: n !

   !******** read data from file1 ********
   open(10,file=file1)
     do
        read(10,*,iostat=code) idm,rdm,r1
        if(code<0)exit
        h(count)=r1
        count=count+1 
     end do
   close(10)

   !******** input data to x0(n) ********
   do n=0,s-1
      a(n)=h(n)+complex_j*0.0d0
   end do

   !******** read data re-confirm!! ******** 
   open(unit=22,file="readdata_fic_dif_j")
      do n=0,s-1
         write(22,*) n,real(a(n)),aimag(a(n))
         a(n)=h(n)+complex_j*0.0d0
      end do   
   close(unit=22)


   write(6,*) 'data input succeeded...'


end subroutine read_data
  
!*****************************************************************************
! the fictitious domain to the diffusive domain 
!****************************************************************************
   
subroutine CDFT_ftod
  use fft_lib
  implicit none
  integer :: i,j,k
  do j=0,s-1
     do k=0,s-1

!        x(j) = x(j)+sqrt((-2.0d0*s*t)/(complex_j*j))*exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))&
!             *exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a(k)
        x(j) = x(j)+sqrt((-2.0d0*s*t)/(complex_j*j))*&
                    exp(-2.0d0*pi*sqrt(j*t/s)*k)*&
                    exp(2.0d0*complex_j*pi*sqrt(j*t/s)*k)*a(k)


     end do
!     print *, j,x(j)
  end do
  x(0)=1.0d4+2.0d3*complex_j
 end subroutine   
     

!************************************************************************
! output spectrum
!************************************************************************
subroutine output
   use fft_lib
   implicit none
   integer :: n 

   !******** output ********
   open(unit=11,file="dft_gauss1f.dat")
      do n=0, s-1
          write(11,*) n,n/(s*t),cdabs(x(n)),real(x(n)),aimag(x(n))
      end do

!   open(unit=12,file="dft4_hagen_iJ.dat")
!      do n=1, s-1
!         b(n)=1/x(n)
!      end do
      
!      do n=0,s-1
!      write(12,*) n,n/(s*t),cdabs(b(n)),real(b(n)),aimag(b(n))
!      end do


  close(unit=11)
  close(unit=12)

end subroutine output


