module input
  implicit none
  integer,parameter :: s=16384
  integer :: idm
  real(8) :: re1(s),re2(s)
  real(8) :: im1(s),im2(s)
  complex :: e1(s),e2(s),e3(s),e4(s),e5(s),e6(s),e7(s),e8(s),e9(s),e10(s),e11(s),e12(s),e13(s),e14(s),e15(s)
  complex :: j(s),k(s)
  complex :: G1(s),G2(s),G3(s),G4(s),G5(s),G6(s),G7(s),G8(s),G9(s),G10(s),G11(s),G12(s),G13(s),G14(s),G15(s)
  complex,parameter :: complex_j=(0.0d0,1.0d0) ! imaginary number
  real(8) :: rdm,r1,r2,d,i1,i2 !
  real(8) :: t
end module input


program gap
  use input
  implicit none
  integer :: n

  !t=6.67398930899843271E-003
  t=5.33919144719874617E-002
 !t=2.11050072960483441E-003 
 ! t=3.61496891435735771E-003 
! t=1.68840058368386753E-002  

   open(10,file="dif_m6_1_1")
   open(11,file="dif_m6_2_1")
   open(12,file="dif_m6_3_1")
   open(13,file="dif_m6_4_1")
   open(14,file="dif_m6_5_1")
   open(15,file="dif_m6_6_1")
   open(16,file="dif_m6_7_1")
   open(17,file="dif_m6_8_1")
   open(18,file="dif_m6_9_1")
   open(19,file="dif_m6_10_1")
   open(20,file="dif_m6_11_1")
   open(21,file="dif_m6_12_1")
   open(22,file="dif_m6_13_1")
   open(23,file="dif_m6_14_1")
   open(24,file="dif_m6_15_1")

   do n=0,s-1
   read(10,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e1(n)=re1(n)+complex_j*im1(n)   
   end do
   close(10)

   do n=0,s-1
   read(11,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e2(n)=re1(n)+complex_j*im1(n)   
   end do
   close(11)

   do n=0,s-1
   read(12,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e3(n)=re1(n)+complex_j*im1(n)   
   end do
   close(12)

   do n=0,s-1
   read(13,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e4(n)=re1(n)+complex_j*im1(n)   
   end do
   close(13)

   do n=0,s-1
   read(14,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e5(n)=re1(n)+complex_j*im1(n)   
   end do
   close(14)

   do n=0,s-1
   read(15,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e6(n)=re1(n)+complex_j*im1(n)   
   end do
   close(15)

   do n=0,s-1
   read(16,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e7(n)=re1(n)+complex_j*im1(n)   
   end do
   close(16)

   do n=0,s-1
   read(17,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e8(n)=re1(n)+complex_j*im1(n)   
   end do
   close(17)

   do n=0,s-1
   read(18,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e9(n)=re1(n)+complex_j*im1(n)   
   end do
   close(18)

   do n=0,s-1
   read(19,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e10(n)=re1(n)+complex_j*im1(n)   
   end do
   close(19)

   do n=0,s-1
   read(20,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e11(n)=re1(n)+complex_j*im1(n)   
   end do
   close(20)

   do n=0,s-1
   read(21,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e12(n)=re1(n)+complex_j*im1(n)   
   end do
   close(21)

   do n=0,s-1
   read(22,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e13(n)=re1(n)+complex_j*im1(n)   
   end do
   close(22)

   do n=0,s-1
   read(23,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e14(n)=re1(n)+complex_j*im1(n)   
   end do
   close(23)

   do n=0,s-1
   read(24,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e15(n)=re1(n)+complex_j*im1(n)   
   end do
   close(24)





   open(50,file="dft_J.dat")

   do n=0,s-1
   read(50,*) idm,rdm,d,r2,i2
      re2(n)=r2
      im2(n)=i2
      j(n)=re2(n)+complex_j*im2(n)    
   end do
   close(50)

   open(50,file="fft_gauss1.dat")

   do n=0,s-1
   read(50,*) idm,rdm,d,r2,i2
      re2(n)=r2
      im2(n)=i2
      k(n)=re2(n)+complex_j*im2(n)    
   end do
   close(50)
 
   do n=0,s-1
      G1(n)=e1(n)/j(n)
      G2(n)=e2(n)/j(n)
      G3(n)=e3(n)/j(n)
      G4(n)=e4(n)/j(n)
      G5(n)=e5(n)/j(n)
      G6(n)=e6(n)/j(n)
      G7(n)=e7(n)/j(n)
      G8(n)=e8(n)/j(n)
      G9(n)=e9(n)/j(n)
      G10(n)=e10(n)/j(n)
      G11(n)=e11(n)/j(n)
      G12(n)=e12(n)/j(n)
      G13(n)=e13(n)/j(n)
      G14(n)=e14(n)/j(n)
      G15(n)=e15(n)/j(n)
   end do

   open(30,file="Eobsf_m6hez_1.dat") 
   open(31,file="Eobsf_m6hez_2.dat") 
   open(32,file="Eobsf_m6hez_3.dat") 
   open(33,file="Eobsf_m6hez_4.dat") 
   open(34,file="Eobsf_m6hez_5.dat") 
   open(35,file="Eobsf_m6hez_6.dat") 
   open(36,file="Eobsf_m6hez_7.dat") 
   open(37,file="Eobsf_m6hez_8.dat") 
   open(38,file="Eobsf_m6hez_9.dat") 
   open(39,file="Eobsf_m6hez_10.dat") 
   open(40,file="Eobsf_m6hez_11.dat") 
   open(41,file="Eobsf_m6hez_12.dat") 
   open(42,file="Eobsf_m6hez_13.dat") 
   open(43,file="Eobsf_m6hez_14.dat") 
   open(44,file="Eobsf_m6hez_15.dat")
  

   do n=0,s-1
      e1(n)=G1(n)*k(n)
      e2(n)=G2(n)*k(n)
      e3(n)=G3(n)*k(n)
      e4(n)=G4(n)*k(n)
      e5(n)=G5(n)*k(n)
      e6(n)=G6(n)*k(n)
      e7(n)=G7(n)*k(n)
      e8(n)=G8(n)*k(n)
      e9(n)=G9(n)*k(n)
      e10(n)=G10(n)*k(n)
      e11(n)=G11(n)*k(n)
      e12(n)=G12(n)*k(n)
      e13(n)=G13(n)*k(n)
      e14(n)=G14(n)*k(n)
      e15(n)=G15(n)*k(n)


      write(30,*) n,n*t,cabs(e1(n)),real(e1(n)),aimag(e1(n))
      write(31,*) n,n*t,cabs(e2(n)),real(e2(n)),aimag(e2(n))
      write(32,*) n,n*t,cabs(e3(n)),real(e3(n)),aimag(e3(n))
      write(33,*) n,n*t,cabs(e4(n)),real(e4(n)),aimag(e4(n))
      write(34,*) n,n*t,cabs(e5(n)),real(e5(n)),aimag(e5(n))
      write(35,*) n,n*t,cabs(e6(n)),real(e6(n)),aimag(e6(n))
      write(36,*) n,n*t,cabs(e7(n)),real(e7(n)),aimag(e7(n))
      write(37,*) n,n*t,cabs(e8(n)),real(e8(n)),aimag(e8(n))
      write(38,*) n,n*t,cabs(e9(n)),real(e9(n)),aimag(e9(n))
      write(39,*) n,n*t,cabs(e10(n)),real(e10(n)),aimag(e10(n))
      write(40,*) n,n*t,cabs(e11(n)),real(e11(n)),aimag(e11(n))
      write(41,*) n,n*t,cabs(e12(n)),real(e12(n)),aimag(e12(n))
      write(42,*) n,n*t,cabs(e13(n)),real(e13(n)),aimag(e13(n))
      write(43,*) n,n*t,cabs(e14(n)),real(e14(n)),aimag(e14(n))
      write(44,*) n,n*t,cabs(e15(n)),real(e15(n)),aimag(e15(n))


   end do

   close(30)
   close(31)
   close(32)
   close(33)
   close(34)
   close(35)
   close(36)
   close(37)
   close(38)
   close(39)
   close(40)
   close(41)
   close(42)
   close(43)
   close(44)


end program

