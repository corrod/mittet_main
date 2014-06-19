module fft_lib
!============================================================================
   !******** physical constant ******** 
   complex,parameter :: complex_j=(0.0d0,1.0d0) ! imaginary number
   real(8),parameter :: pi=3.1415927d0          ! circular constant
   !******** time data ********
   real(8),dimension(:),allocatable :: h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15        ! input time domain data
   !******** frequency spectrum data ********
   complex(kind(0d0)),dimension(:),allocatable :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15        
   complex(kind(0d0)),dimension(:),allocatable :: x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15        ! output spectrum

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
   use fft_lib
   implicit none
   integer :: n,k !  
   integer :: ts,te,t_rate,t_max,diff
  call system_clock(ts)

   call initial
   call read_data
   call CDFT
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
   file1= 'receiver_0108_m6hEy.dat'
   !==== sampling numbers ====
   print *, "enter total sampling numbers N"
   s=4096
!   s=65536  !2**16
!   s=131072  !2**17
!   s=262144
!   s=10000 ! 2**13
   !==== sampling interval T ====
   print *, "enter sampling width T"
   t=1.14315354400754132d-3
!   t=3.61496891435735771E-003 
   gamma=nint(log10(real(s))/log10(real(2)))
   t0=s*t

   !******** memory allocation ********
   allocate(h1(0:s-1),h2(0:s-1),h3(0:s-1),h4(0:s-1),h5(0:s-1),h6(0:s-1),h7(0:s-1),h8(0:s-1),h9(0:s-1),h10(0:s-1),h11(0:s-1)&
,h12(0:s-1),h13(0:s-1),h14(0:s-1),h15(0:s-1),stat=code)
   allocate(x1(0:s-1),x2(0:s-1),x3(0:s-1),x4(0:s-1),x5(0:s-1),x6(0:s-1),x7(0:s-1),x8(0:s-1),x9(0:s-1),x10(0:s-1),x11(0:s-1)&
,x12(0:s-1),x13(0:s-1),x14(0:s-1),x15(0:s-1),stat=code)
   allocate(a1(0:s-1),a2(0:s-1),a3(0:s-1),a4(0:s-1),a5(0:s-1),a6(0:s-1),a7(0:s-1),a8(0:s-1),a9(0:s-1),a10(0:s-1),a11(0:s-1)&
,a12(0:s-1),a13(0:s-1),a14(0:s-1),a15(0:s-1),stat=code)

   if(code/=0) stop "*** not enough memory <fft_v1> ***"
   !******** array initialize ********
   do n=0,s-1
      h1(n)=0.0d0
   end do

   do n=0,s-1
      x1(n)=(0.0d0,0.0d0) 
   end do

   do n=0,s-1
     a1(n)=(0.0d0,0.0d0)
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
        read(10,*,iostat=code) idm,rdm,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15
        if(code<0)exit
        h1(count)=r11
        count=count+1 
     end do
   close(10)

   !******** input data to x0(n) ********
   do n=0,s-1
      a1(n)=h1(n)+complex_j*0.0d0
   end do

   !******** read data re-confirm!! ******** 
   open(unit=22,file="fic4_t6_11.dat")
!    open(unit=22,file="readdata")
      do n=0,s-1
         write(22,*) n,n*t,real(a1(n)),aimag(a1(n)),cdabs(a1(n))
         a1(n)=h1(n)+complex_j*0.0d0
      end do   
   close(unit=22)


   write(6,*) 'data input succeeded...'


end subroutine read_data
  
!*****************************************************************************
! DFT
!*****************************************************************************
   
subroutine CDFT
  use fft_lib
  implicit none
  integer i,j,k
  do j=0,s-1
     do k=0,s-1
        x1(j) = x1(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a1(k)
     end do
     print *, j
  end do
  end subroutine   
     

!************************************************************************
! output spectrum
!************************************************************************
subroutine output
   use fft_lib
   implicit none
   integer :: n 

   !******** output ********
   open(unit=11,file="dif4_m6_11")
  do n=0, s-1
          write(11,*) n,n/(s*t),cdabs(x1(n)),real(x1(n)),aimag(x1(n))
      end do
  close(11)
  open(unit=12,file="dif4_m6_11a")
  do n=0, s-1
          write(12,*) x1(n)
      end do
  close(12)

end subroutine output


