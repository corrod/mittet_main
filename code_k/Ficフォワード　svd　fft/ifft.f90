module fft_lib
!============================================================================
   !******** physical constant ******** 
   complex,parameter :: complex_j=(0.0d0,1.0d0) ! imaginary number
   real(8),parameter :: pi=3.1415927d0          ! circular constant
   !******** time data ********
   real(8),dimension(:),allocatable :: h        ! input time domain data
   real(8),dimension(:),allocatable :: c        ! input time domain data
   !******** frequency spectrum data ********
   complex,dimension(:),allocatable :: x        ! output spectrum
   !******** fft parameters ********
   integer :: s       ! total sampling numbers
   integer :: gamma   ! s=2**gamma -> gamma=log10(s)/log10(2) 
   real(8) :: t       ! delta time (t0/s) [sec]
   real(8) :: t0      ! assumed cycle [sec]
   integer :: n2      ! binominal point when row is l  eq.(10-21)
   integer :: nu1     ! bit shift factor
   integer :: m       ! k which bit has shifted by a factor of bits
   complex :: wp      ! phase rotation 
   complex :: t1      !  // with x
   integer :: p       ! phase rotation factor of wp
   integer :: skip    ! skip counter
   complex :: ctmp    ! temporary complex variable
   integer :: rtmp    ! temporary real variable 
   !******** input file ********
   integer :: code          ! 
   integer :: count=0       !
   integer :: idm     ! dummy data in integer type
   real(8) :: rdm     ! dummy data in real type
   real(8) :: col1,d,col2    ! data under consideration
   character(90) file1      !
!============================================================================
end module fft_lib  

program fft_ver3
!***********************************************************************     
!   Fast Fourier Transform 
!***********************************************************************
   use fft_lib
   implicit none
   integer :: n,k 
   integer :: ts,te,t_rate,t_max,diff
  call system_clock(ts)

   call initial
   call read_data
   call fft_main
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
   file1= 'fft4_Eobst_m7_11.dat'
   open(unit=11,file="dif4_t7_cut.dat")

   !==== sampling numbers ====
   s=4096  !2**17
   !==== sampling interval T ====
   print *, "enter sampling width T"
 
 !t=6.67398930899843271E-003 
 !  t=2.11050081998109818E-003 
 !  t=1.68840065598487854E-002 
 !  t=5.33919148147106171E-002
   t=0.21356765788794985 
   !t=0.85427063703536987   
   !==== gamma ====
   gamma=nint(log10(real(s))/log10(real(2)))
   ! assumed cycle 
   t0=s*t
   print *, 'dt=', 1.0d0/t0

   !******** memory allocation ********
   allocate(h(0:s-1),x(0:s-1),c(0:s-1),stat=code)
   if(code/=0) stop "*** not enough memory <fft_v1> ***"
   !******** array initialize ********
   do n=0,s-1
      h(n)=0.0d0
      c(n)=0.0d0
      x(n)=(0.0d0,0.0d0) 
   end do
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
        read(10,*,iostat=code) idm,rdm,d,col1,col2
        if(code<0)exit
        h(count)=col1
        c(count)=col2
        count=count+1 
     end do
   close(10)

   !******** input data to x0(n) ********
   do n=0,s-1
      if (n .le. 85) then
      x(n)=conjg(h(n)+complex_j*c(n))
   else
      x(n)=(0.0d0,0.0d0)
      end if
   end do

   !******** read data re-confirm!! ******** 
   open(unit=22,file="E4_m7_11cut.dat")
      do n=0,s-1
         write(22,*) n,n*t,cabs(x(n)),real(x(n)),aimag(x(n))
         if (n.le.85) then
         x(n)=h(n)+complex_j*c(n)
         else
         x(n)=(0.0d0,0.0d0)
         end if
      end do   
   close(unit=22)
end subroutine read_data

!*************************************************************************
! fft main
!*************************************************************************
subroutine fft_main
   use fft_lib
   implicit none
   integer :: ibr     ! bit order reversal function
   integer :: l,k,i,j  !

   !******** column calculation from l=1 to gamma ********
   do l=1,gamma

      write(6,*) l ,gamma
      n2=s/2**l !2,12
      nu1=gamma-l
   
      !******** row calculation from n=1 to s-1 ********
      k=0
      skip=1 !4

      k_loop: do
         if(k>s-1)exit

         m=int(k/2**nu1) !5
         p=ibr(m,gamma)

         !******** dual node pair ********
         wp=exp(-complex_j*(2.0d0*pi*p/s)) !6
         t1=wp*x(k+n2)
         x(k+n2)=x(k)-t1
         x(k)=x(k)+t1
         k=k+1 !7
         !******** skip condition ********
         if(skip==n2)then !8
            k=k+n2
            skip=1
         else
            skip=skip+1 !10           
         end if
      end do k_loop
   end do ! l loop

   !******** re-arrangement of the spectrum order ********
   do i=0,s-1
      j=ibr(i,gamma) !
      if(j>=i)then !14
         ctmp=x(i)
         x(i)=x(j)
         x(j)=ctmp
      end if
   end do

end subroutine fft_main

!************************************************************************
! output spectrum
!************************************************************************
subroutine output
   use fft_lib
   implicit none
   integer :: n !

   !******** output ********
!   open(unit=11,file="time_yz_m7_2_7.dat")
      do n=0,s-1
!          write(11,*) n,n/(s*t),real(conjg(x(n)))/s,aimag(conjg(x(n)))/s,cabs(conjg(x(n)))/s 
          write(11,*) n,n/(s*t),real(conjg(x(n)))/s,aimag(conjg(x(n)))/s,cabs(conjg(x(n)))/s 
   end do
   close(unit=11)

end subroutine output

!************************************************************************
! bit order reversal function 
!************************************************************************
integer function ibr(m,gamma)
   implicit none
   integer :: i,m,im !
   integer :: gamma !
   integer :: code !
   integer,dimension(:),allocatable :: a ! 

   !******** memory allocation ********
   allocate(a(1:gamma),stat=code)
   if(code/=0) stop "*** not enough memory <bit> ***"
   !******** array initialize ********
   do i=1,gamma 
      a(i)=0
   end do

   !******** 10 to 2 transform ********
   do i=1,gamma
      a(i)=int(m/2**(i-1))-2*(int(m/2**i))
   end do
   !******** 2 to 10 transform with reversal ********
   im=0
   do i=1,gamma 
      im=im+a(i)*2**(gamma-i)       
   end do   

   ibr=im

end function ibr



