!******************************************************************************
!          
!******************************************************************************
subroutine output2
  use fdtd
  implicit none

  integer :: i,j,k

   open (40,file="receiver_0108_m7hEy.dat")
   i=nx/2
   j=ny/2
   k=nz/2+15+1
   write(40,*) n,n*dt,ey(i,j-35,k),ey(i,j-30,k),ey(i,j-25,k),ey(i,j-20,k),ey(i,j-15,k),ey(i,j-10,k),ey(i,j-5,k),ey(i,j,k)&
,ey(i,j+5,k),ey(i,j+10,k),ey(i,j+15,k),ey(i,j+20,k),ey(i,j+25,k),ey(i,j+30,k),ey(i,j+35,k)

   open (41,file="receiver_0108_m7hEz.dat")
   i=nx/2
   j=ny/2
   k=nz/2+15+1
   write(41,*) n,n*dt,ez(i,j-35,k),ez(i,j-30,k),ez(i,j-25,k),ez(i,j-20,k),ez(i,j-15,k),ez(i,j-10,k),ez(i,j-5,k),ez(i,j,k)&
,ez(i,j+5,k),ez(i,j+10,k),ez(i,j+15,k),ez(i,j+20,k),ez(i,j+25,k),ez(i,j+30,k),ez(i,j+35,k)

   open (42,file="receiver_0108_m7hEx.dat")
   i=nx/2
   j=ny/2
   k=nz/2+15+1
   write(42,*) n,n*dt,ex(i,j-35,k),ex(i,j-30,k),ex(i,j-25,k),ex(i,j-20,k),ex(i,j-15,k),ex(i,j-10,k),ex(i,j-5,k),ex(i,j,k)&
,ex(i,j+5,k),ex(i,j+10,k),ex(i,j+15,k),ex(i,j+20,k),ex(i,j+25,k),ex(i,j+30,k),ex(i,j+35,k)

 open (43,file="receiver_0108_m7hHx.dat")
   i=nx/2
   j=ny/2
   k=nz/2+15+1
   write(43,*) n,n*dt,hx(i,j-35,k),hx(i,j-30,k),hx(i,j-25,k),hx(i,j-20,k),hx(i,j-15,k),hx(i,j-10,k),hx(i,j-5,k),hx(i,j,k)&
,hx(i,j+5,k),hx(i,j+10,k),hx(i,j+15,k),hx(i,j+20,k),hx(i,j+25,k),hx(i,j+30,k),hx(i,j+35,k)

 open (44,file="receiver_0108_m7hHy.dat")
   i=nx/2
   j=ny/2
   k=nz/2+15+1
   write(44,*) n,n*dt,hy(i,j-35,k),hy(i,j-30,k),hy(i,j-25,k),hy(i,j-20,k),hy(i,j-15,k),hy(i,j-10,k),hy(i,j-5,k),hy(i,j,k)&
,hy(i,j+5,k),hy(i,j+10,k),hy(i,j+15,k),hy(i,j+20,k),hy(i,j+25,k),hy(i,j+30,k),hy(i,j+35,k)

 open (45,file="receiver_0108_m7hHz.dat")
   i=nx/2
   j=ny/2
   k=nz/2+15+1
   write(45,*) n,n*dt,hz(i,j-35,k),hz(i,j-30,k),hz(i,j-25,k),hz(i,j-20,k),hz(i,j-15,k),hz(i,j-10,k),hz(i,j-5,k),hz(i,j,k)&
,hz(i,j+5,k),hz(i,j+10,k),hz(i,j+15,k),hz(i,j+20,k),hz(i,j+25,k),hz(i,j+30,k),hz(i,j+35,k)

   return
end subroutine
