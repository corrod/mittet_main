module fft_lib
!============================================================================
   !******** physical constant ********
   complex(kind(0d0)),parameter :: complex_j=(0.0d0,1.0d0) ! imaginary number
   real(8),parameter :: pi=3.14159265358979323846d0          ! circular constant
   !******** time data ********
   real(8),dimension(:),allocatable :: h,rw,iw        ! input time domain data
   !******** frequency spectrum data ********
   complex(kind(0d0)),dimension(:),allocatable :: x,w        ! output spectrum
   !******** fft parameters ********
   integer :: s,u,v       ! total sampling numbers
   integer :: gamma   ! s=2**gamma -> gamma=log10(s)/log10(2)
   real(8) :: t,d       ! delta time (t0/s) [sec]
   real(8) :: t0      ! assumed cycle [sec]
   integer :: n2      ! binominal point when row is l  eq.(10-21)
   integer :: nu1     ! bit shift factor
   integer :: m       ! k which bit has shifted by a factor of bits
   complex(kind(0d0)) :: wp      ! phase rotation
   complex(kind(0d0)) :: t1      !  // with x
   integer :: p       ! phase rotation factor of wp
   integer :: skip    ! skip counter
   complex(kind(0d0)) :: ctmp    ! temporary complex variable
   integer :: rtmp    ! temporary real variable
   !******** input file ********
   integer :: code          !
   integer :: count=0       !
   integer :: idm     ! dummy data in integer type
   real(8) :: rdm,r1,r2,r3     ! dummy data in real type
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
!   print *, "enter input data file name"
!   read(5,*) file1
   file1='func_gauss1d.dat'
   !==== sampling numbers ====
   print *, "enter total sampling numbers N"

!   s=131072
   s=16384  !2**14
   u=16384
   v=200
   !==== sampling interval T ====
   print *, "enter sampling width T"
   t=1.14315354400754132d-3
!   t=3.61496891435735771E-003

   !==== gamma ====
   gamma=nint(log10(real(s))/log10(real(2)))
   ! assumed cycle
   t0=s*t

   !******** memory allocation ********
   allocate(h(0:s-1),x(0:s-1),rw(0:s-1),iw(0:s-1),w(0:s-1),stat=code)
   if(code/=0) stop "*** not enough memory <fft_v1> ***"
   !******** array initialize ********
   do n=0,s-1
      h(n)=0.0d0
      rw(n)=0.0d0
      iw(n)=0.0d0
      x(n)=(0.0d0,0.0d0)
      w(n)=(0.0d0,0.0d0)
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
      x(n)=h(n)+complex_j*0.0d0
   end do

   do n=0,s-1
!      rw(n)=0.5d0-0.5d0*dcos(2.0d0*pi*(1.0d0*n/s))

      if(n.le.v) then
         rw(n)=0.5d0-0.5d0*dcos(2.0d0*pi*(1.0d0*n/(2.0d0*v)))
      else if (n.ge.s-1-v) then
         rw(n)=0.5d0-0.5d0*dcos(2.0d0*pi*(1.0d0*(s-1-n)/(2.0d0*v)))
      else
         rw(n)=1.0d0
      end if
   end do

   do n=0,s-1
     ! iw(n)=0.5d0-0.5d0*dcos(2.0d0*pi*(1.0d0*n/s))
      if(n.le.v) then
         iw(n)=0.5d0-0.5d0*dcos(2.0d0*pi*(1.0d0*n/(2.0d0*v)))
      else if (n.ge.s-1-v) then
         iw(n)=0.5d0-0.5d0*dcos(2.0d0*pi*(1.0d0*(s-1-n)/(2.0d0*v)))
      else
         iw(n)=1.0d0
      end if

   end do

   open(30,file='hanning.dat')
   do n=0,s-1
      w(n)=rw(n)+complex_j*iw(n)
      write(30,*) n,rw(n),iw(n),real(w(n))
   end do
   close(30)
!-------------------------------------------------------------------------------
!hanning window
!-------------------------------------------------------------------------------
!   do n=0,s-1
!      x(n)=x(n)*w(n)
!   end do

   !******** read data re-confirm!! ********
!   open(unit=22,file="Eobst4_m7_11.dat")
!     open(unit=22,file="dif_t7_str2.dat")
   open(unit=22,file="readdata")
      do n=0,s-1
         write(22,*) n,n*t,real(x(n)),aimag(x(n))
!         x(n)=h(n)+complex_j*0.0d0
      end do
   close(unit=22)


   write(6,*) 'data input succeeded...'


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
         wp=exp(complex_j*(2.0d0*pi*p/s)) !6
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

   write(6,*) 'main spectrum calculation has succeeded....'


   !******** re-arrangement of the spectrum order ********
   do i=0,s-1
!      write(6,*) i, s-1
      j=ibr(i,gamma) !13

      if(j>=i)then !14
         ctmp=x(i)
         x(i)=x(j)
         x(j)=ctmp
      end if

   end do

   write(6,*) 're-arrangement of the spectrum has succeeded....'
   write(6,*) 'following is your information....'
   write(6,'(1x, "sampling number     :",i6)') s
   write(6,'(1x, "gamma               :",i6)') gamma
   write(6,'(1x, "sampling interval   :",d18.10)') t
   write(6,'(1x, "assumed period t0   :",d18.10)') t0
   write(6,'(1x, "delta freq          :",d18.10)') 1.0d0/t0


end subroutine fft_main

!************************************************************************
! output spectrum
!************************************************************************
subroutine output
   use fft_lib
   implicit none
   integer :: n !

   !******** output ********
   open(unit=11,file="fft_gauss1d.dat")
      do n=0, s-1
          write(11,*) n,n/(s*t),cdabs(x(n)),real(x(n)),aimag(x(n))
!         write(11,*) n,n/(s*t),2.0d0*cdabs(x(n)),2.0d0*real(x(n)),2.0d0*aimag(x(n))
    print *, n
      end do
   close(unit=11)

end subroutine output


!************************************************************************
! bit order reversal function
!
! data number :: 2**r
! a  :: binary digit (bit)
! m  :: decimal digit
! im :: inversed decimal digit
!
! m = a(r) a(r-1) a(r-2).....  a(3)   a(2)   a(1)
! im= a(1) a(2)   a(3)  ...... a(r-2) a(r-1) a(r)
!
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

