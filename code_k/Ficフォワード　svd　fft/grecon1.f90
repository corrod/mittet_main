module input
  implicit none
  integer,parameter :: s=16384
  integer :: idm
  real(8) :: re1(s),re2(s)
  real(8) :: im1(s),im2(s)
  complex(kind(0d0)) :: e1(s),e2(s),e3(s),e4(s),e5(s),e6(s),e7(s),e8(s),e9(s),e10(s),e11(s),e12(s),e13(s),e14(s),e15(s)
  complex(kind(0d0)) :: j(s),k(s)
  complex(kind(0d0)) :: G1(s),G2(s),G3(s),G4(s),G5(s),G6(s),G7(s),G8(s),G9(s),G10(s),G11(s),G12(s),G13(s),G14(s),G15(s)
  complex(kind(0d0)),parameter :: complex_j=(0.0d0,1.0d0) ! imaginary number
  real(8) :: rdm,r1,r2,d,i1,i2 !
  real(8) :: t
end module input


program gap
  use input
  implicit none
  integer :: n

  !t=6.67398930899843271E-003
  t=5.33919144719874617E-002
 ! t=0.21356765788794985 
 !t=2.11050072960483441E-003 


   open(10,file="fft_Eobst_m7_11.dat")

   do n=0,s-1
   read(10,*) idm,rdm,d,r1,i1
      re1(n)=r1
      im1(n)=i1
      e1(n)=re1(n)+complex_j*im1(n)   
   end do
   close(10)





   open(50,file="fft_gauss1.dat")

   do n=0,s-1
   read(50,*) idm,rdm,d,r2,i2
      re2(n)=r2
      im2(n)=i2
      j(n)=re2(n)+complex_j*im2(n)    
   end do
   close(50)

   open(50,file="dft_J.dat")

   do n=0,s-1
   read(50,*) idm,rdm,d,r2,i2
      re2(n)=r2
      im2(n)=i2
      k(n)=re2(n)+complex_j*im2(n)    
   end do
   close(50)
 
   do n=0,s-1
      G1(n)=e1(n)/j(n)
   end do

   open(30,file="Ef_obsgc_11.dat") 
   do n=0,s-1
      e1(n)=G1(n)*k(n)
      write(30,*) n,n*t,cdabs(e1(n)),real(e1(n)),aimag(e1(n))
   end do
   close(30)
 
   open(31,file="Ef_obsgc_11a.dat") 
    do n=0,s-1
      write(31,*) e1(n)
   end do
  close(31)


end program

