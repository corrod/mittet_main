module input
  implicit none
  integer,parameter :: s=4096
  integer :: idm
  real(8) :: g1(s),g2(s),g3(s),g4(s),g5(s),g6(s),g7(s),g8(s),g9(s),g10(s),g11(s),g12(s),g13(s),g14(s),g15(s)
  real(8) :: h1(s),h2(s),h3(s),h4(s),h5(s),h6(s),h7(s),h8(s),h9(s),h10(s),h11(s),h12(s),h13(s),h14(s),h15(s)
  real(8) :: i1(s),i2(s),i3(s),i4(s),i5(s),i6(s),i7(s),i8(s),i9(s),i10(s),i11(s),i12(s),i13(s),i14(s),i15(s)
  real(8) :: rdm,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,d ! dummy data in real type
  real(8) :: t
end module input


program gap
  use input
  implicit none
  integer :: n

  t=1.14315354400754132d-3
 !t=0.21356765788794985 

   open(10,file="fic4_t7_11t14.dat")


   do n=0,s-1

   read(10,*) idm,rdm,r1,d


      i1(n)=r1


   end do

   close(10)
 
   open(30,file="fic4_t6_11t14.dat")

   do n=0,s-1

   read(30,*) idm,rdm,r1,d


      h1(n)=r1
 

   end do

   close(30)
 
   do n=0,s-1
      g1(n)=i1(n)-h1(n)

   end do


   open(unit=22,file="gap4_11t14.dat")
   do n=0,s-1
      write(22,*) n,n*t,g1(n)
   end do
   close(unit=22)

!   open(23,file="gapreal_amp.dat")
!      write(23,*) maxval(real(g1))-minval(real(g1))
!      write(23,*) maxval(real(g2))-minval(real(g2))
!      write(23,*) maxval(real(g3))-minval(real(g3))
!      write(23,*) maxval(real(g4))-minval(real(g4))
!      write(23,*) maxval(real(g5))-minval(real(g5))
!      write(23,*) maxval(real(g6))-minval(real(g6))
!      write(23,*) maxval(real(g7))-minval(real(g7))
!      write(23,*) maxval(real(g8))-minval(real(g8))
!      write(23,*) maxval(real(g9))-minval(real(g9))
!      write(23,*) maxval(real(g10))-minval(real(g10))
!      write(23,*) maxval(real(g11))-minval(real(g11))
!     write(23,*) maxval(real(g12))-minval(real(g12))
!     write(23,*) maxval(real(g13))-minval(real(g13))
!    write(23,*) maxval(real(g14))-minval(real(g14))
!    write(23,*) maxval(real(g15))-minval(real(g15))
! close(23)
   


end program

