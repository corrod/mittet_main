module fft_lib
!============================================================================
   !******** physical constant ******** 
   complex(kind(0d0)),parameter :: complex_j=(0.0d0,1.0d0) ! imaginary number
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
   file1= 'receiver_0108_m7hEy.dat'
   !==== sampling numbers ====
   print *, "enter total sampling numbers N"
   s=16384
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
      h2(n)=0.0d0
      h3(n)=0.0d0
      h4(n)=0.0d0
      h5(n)=0.0d0
      h6(n)=0.0d0
      h7(n)=0.0d0
      h8(n)=0.0d0
      h9(n)=0.0d0
      h10(n)=0.0d0
      h11(n)=0.0d0
      h12(n)=0.0d0
      h13(n)=0.0d0
      h14(n)=0.0d0
      h15(n)=0.0d0
   end do

   do n=0,s-1
      x1(n)=(0.0d0,0.0d0) 
      x2(n)=(0.0d0,0.0d0) 
      x3(n)=(0.0d0,0.0d0) 
      x4(n)=(0.0d0,0.0d0) 
      x5(n)=(0.0d0,0.0d0) 
      x6(n)=(0.0d0,0.0d0) 
      x7(n)=(0.0d0,0.0d0) 
      x8(n)=(0.0d0,0.0d0) 
      x9(n)=(0.0d0,0.0d0) 
      x10(n)=(0.0d0,0.0d0) 
      x11(n)=(0.0d0,0.0d0) 
      x12(n)=(0.0d0,0.0d0) 
      x13(n)=(0.0d0,0.0d0) 
      x14(n)=(0.0d0,0.0d0) 
      x15(n)=(0.0d0,0.0d0) 
   end do

   do n=0,s-1
     a1(n)=(0.0d0,0.0d0)
     a2(n)=(0.0d0,0.0d0) 
     a3(n)=(0.0d0,0.0d0)
     a4(n)=(0.0d0,0.0d0) 
     a5(n)=(0.0d0,0.0d0)
     a6(n)=(0.0d0,0.0d0) 
     a7(n)=(0.0d0,0.0d0)
     a8(n)=(0.0d0,0.0d0) 
     a9(n)=(0.0d0,0.0d0)
     a10(n)=(0.0d0,0.0d0) 
     a11(n)=(0.0d0,0.0d0)
     a12(n)=(0.0d0,0.0d0) 
     a13(n)=(0.0d0,0.0d0)
     a14(n)=(0.0d0,0.0d0) 
     a15(n)=(0.0d0,0.0d0)
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
        h1(count)=r1
        h2(count)=r2
        h3(count)=r3
        h4(count)=r4
        h5(count)=r5
        h6(count)=r6
        h7(count)=r7
        h8(count)=r8
        h9(count)=r9
        h10(count)=r10
        h11(count)=r11
        h12(count)=r12
        h13(count)=r13
        h14(count)=r14
        h15(count)=r15
       count=count+1 
     end do
   close(10)

   !******** input data to x0(n) ********
   do n=0,s-1
      a1(n)=h1(n)+complex_j*0.0d0
      a2(n)=h2(n)+complex_j*0.0d0
      a3(n)=h3(n)+complex_j*0.0d0
      a4(n)=h4(n)+complex_j*0.0d0
      a5(n)=h5(n)+complex_j*0.0d0
      a6(n)=h6(n)+complex_j*0.0d0
      a7(n)=h7(n)+complex_j*0.0d0
      a8(n)=h8(n)+complex_j*0.0d0
      a9(n)=h9(n)+complex_j*0.0d0
      a10(n)=h10(n)+complex_j*0.0d0
      a11(n)=h11(n)+complex_j*0.0d0
      a12(n)=h12(n)+complex_j*0.0d0
      a13(n)=h13(n)+complex_j*0.0d0
      a14(n)=h14(n)+complex_j*0.0d0
      a15(n)=h15(n)+complex_j*0.0d0
   end do

   !******** read data re-confirm!! ******** 
   open(unit=22,file="readdata_dft")
      do n=0,s-1
         write(22,*) n,real(a1(n)),aimag(a1(n))
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
        x2(j) = x2(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a2(k)
        x3(j) = x3(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a3(k)
        x4(j) = x4(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a4(k)
        x5(j) = x5(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a5(k)
        x6(j) = x6(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a6(k)
        x7(j) = x7(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a7(k)
        x8(j) = x8(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a8(k)
        x9(j) = x9(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a9(k)
        x10(j) = x10(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a10(k)
        x11(j) = x11(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a11(k)
        x12(j) = x12(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a12(k)
        x13(j) = x13(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a13(k)
        x14(j) = x14(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a14(k)
        x15(j) = x15(j)+exp(-(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt(real(j)*s*t)*k/s))*a15(k)


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
   open(unit=11,file="dif_m7_1 ")
   open(unit=12,file="dif_m7_2 ")
   open(unit=13,file="dif_m7_3 ")
   open(unit=14,file="dif_m7_4 ")
   open(unit=15,file="dif_m7_5 ")
   open(unit=16,file="dif_m7_6 ")
   open(unit=17,file="dif_m7_7 ")
   open(unit=18,file="dif_m7_8 ")
   open(unit=19,file="dif_m7_9 ")
   open(unit=20,file="dif_m7_10 ")
   open(unit=21,file="dif_m7_11 ")
   open(unit=22,file="dif_m7_12 ")
   open(unit=23,file="dif_m7_13")
   open(unit=24,file="dif_m7_14 ")
   open(unit=25,file="dif_m7_15 ")
      do n=0, s-1
          write(11,*) n,n/(s*t),cdabs(x1(n)),real(x1(n)),aimag(x1(n))
          write(12,*) n,n/(s*t),cdabs(x2(n)),real(x2(n)),aimag(x2(n))
          write(13,*) n,n/(s*t),cdabs(x3(n)),real(x3(n)),aimag(x3(n))
          write(14,*) n,n/(s*t),cdabs(x4(n)),real(x4(n)),aimag(x4(n))
          write(15,*) n,n/(s*t),cdabs(x5(n)),real(x5(n)),aimag(x5(n))
          write(16,*) n,n/(s*t),cdabs(x6(n)),real(x6(n)),aimag(x6(n))
          write(17,*) n,n/(s*t),cdabs(x7(n)),real(x7(n)),aimag(x7(n))
          write(18,*) n,n/(s*t),cdabs(x8(n)),real(x8(n)),aimag(x8(n))
          write(19,*) n,n/(s*t),cdabs(x9(n)),real(x9(n)),aimag(x9(n))
          write(20,*) n,n/(s*t),cdabs(x10(n)),real(x10(n)),aimag(x10(n))
          write(21,*) n,n/(s*t),cdabs(x11(n)),real(x11(n)),aimag(x11(n))
          write(22,*) n,n/(s*t),cdabs(x12(n)),real(x12(n)),aimag(x12(n))
          write(23,*) n,n/(s*t),cdabs(x13(n)),real(x13(n)),aimag(x13(n))
          write(24,*) n,n/(s*t),cdabs(x14(n)),real(x14(n)),aimag(x14(n))
          write(25,*) n,n/(s*t),cdabs(x15(n)),real(x15(n)),aimag(x15(n))
      end do
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
  close(17)
  close(18)
  close(19)
  close(20)
  close(21)
  close(22)
  close(23)
  close(24)
  close(25)

end subroutine output


