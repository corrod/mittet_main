!////////////////////////////////////////////////////////////////////////////////////////////////

!////////////////////////////////////////////////////////////////////////////////////////////////

subroutine e_pml(ex,ey,ez,hx,hy,hz)
  use const_para
  implicit none

  integer,parameter :: npml = 5 !PMML層数
  integer,parameter :: MM = 3 !2~4
  real(8),parameter :: RRcoef = 0.01d0!-120.0d0 !反射係数 |R(0)| !R should be [10^-2, 10^-12]
  integer :: l
  integer :: i0,i1,j0,j1,k0,k1,l1,l2,l3,lpmlst(6)
  real(8) :: lpmlii(6,2),lpmljj(6,2),lpmlkk(6,2)
  real(8) :: exy(nx,ny,nz),exz(nx,ny,nz),eyx(nx,ny,nz),eyz(nx,ny,nz),ezx(nx,ny,nz),ezy(nx,ny,nz)
  real(8) :: hxy(nx,ny,nz),hxz(nx,ny,nz),hyx(nx,ny,nz),hyz(nx,ny,nz),hzx(nx,ny,nz),hzy(nx,ny,nz)
  real(8) :: cxe(nx),cye(ny),cze(nz)
  real(8) :: cxh(nx),cyh(ny),czh(nz)
  real(8) :: cxel(nx),cyel(ny),czel(nz)
  real(8) :: cxhl(nx),cyhl(ny),czhl(nz)
  real(8) :: shx,shy,shz,sex,sey,sez
  real(8) :: sig_max,sig(npml)
  complex(kind(0d0)),intent(inout) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
  complex(kind(0d0)),intent(in) :: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)

    sig_max = - (MM+1)*epsi0*cmax/2.0d0/dble(npml)/dx *log(1.0d0/RRcoef)

    do i=1,npml
    sig(i) = sig_max * (dble(i)-1.0d0/2.0d0)/dble(npml)**m
    enddo


    do k=1,npml
      do j=1,npml
        do i=1,npml
      cxe(i) = (1.0d0 - sig(i)*dt/2.0d0/epsi0) / (1.0d0 + sig(i)*dt/2.0d0/epsi0)
      cye(j) = (1.0d0 - sig(j)*dt/2.0d0/epsi0) / (1.0d0 + sig(j)*dt/2.0d0/epsi0)
      cze(k) = (1.0d0 - sig(k)*dt/2.0d0/epsi0) / (1.0d0 + sig(k)*dt/2.0d0/epsi0)

      cxel(i) = (dt/epsi0) / (1.0d0+sig(i)*dt/2.0d0/epsi0)/dx
      cyel(j) = (dt/epsi0) / (1.0d0+sig(j)*dt/2.0d0/epsi0)/dy
      czel(k) = (dt/epsi0) / (1.0d0+sig(k)*dt/2.0d0/epsi0)/dz
        enddo
      enddo
    enddo

   lpmlii(1,1)=1
   lpmljj(1,1)=1
   lpmlkk(1,1)=1
   lpmlii(1,2)=npml+1
   lpmljj(1,2)=ny
   lpmlkk(1,2)=nz
   lpmlii(2,1)=nx-npml
   lpmljj(2,1)=1
   lpmlkk(2,1)=1
   lpmlii(2,2)=nx
   lpmljj(2,2)=ny
   lpmlkk(2,2)=nz

   lpmlii(3,1)=1
   lpmljj(3,1)=1
   lpmlkk(3,1)=1
   lpmlii(3,2)=nx
   lpmljj(3,2)=npml+1
   lpmlkk(3,2)=nz
   lpmlii(4,1)=1
   lpmljj(4,1)=ny-npml
   lpmlkk(4,1)=1
   lpmlii(4,2)=nx
   lpmljj(4,2)=ny
   lpmlkk(4,2)=nz

   lpmlii(5,1)=1
   lpmljj(5,1)=1
   lpmlkk(5,1)=1
   lpmlii(5,2)=nx
   lpmljj(5,2)=ny
   lpmlkk(5,2)=npml+1
   lpmlii(6,1)=1
   lpmljj(6,1)=1
   lpmlkk(6,1)=nz-npml
   lpmlii(6,2)=nx
   lpmljj(6,2)=ny
   lpmlkk(6,2)=nz

  do l=1,6

     i0=lpmlii(l,1)
     i1=lpmlii(l,2)
     j0=lpmljj(l,1)
     j1=lpmljj(l,2)
     k0=lpmlkk(l,1)
     k1=lpmlkk(l,2)

!     l1=lpmlst(l)
     do i=i0,i1-1
        do j=j0+1,j1-1
           do k=k0+1,k1-1
              exy(i,j,k)=cye(j)*exy(i,j,k)+cyel(j)*(hz(i,j,k)-hz(i,j-1,k))
              exz(i,j,k)=cze(k)*exz(i,j,k)+czel(k)*(hy(i,j,k-1)-hy(i,j,k))
              ex(i,j,k)=exy(i,j,k)+exz(i,j,k)
!              l1=l1+1
           enddo
        enddo
     enddo

!     l2=lpmlst(l)
     do i=i0+1,i1-1
        do j=j0,j1-1
           do k=k0+1,k1-1
              eyx(i,j,k)=cxe(i)*eyx(i,j,k)+cxel(i)*(hz(i-1,j,k)-hz(i,j,k))
              eyz(i,j,k)=cze(k)*eyz(i,j,k)+cyel(k)*(hx(i,j,k)-hx(i,j,k-1))
              ey(i,j,k)=eyx(i,j,k)+eyz(i,j,k)
!              l2=l2+1
           enddo
        enddo
     enddo

!      l3=lpmlst(l)
     do i=i0+1,i1-1
        do j=j0+1,j1-1
           do k=k0,k1-1
              ezx(i,j,k)=cxe(i)*ezx(i,j,k)+cxel(i)*(hy(i,j,k)-hy(i-1,j,k))
              ezy(i,j,k)=cye(j)*ezy(i,j,k)+cyel(j)*(hx(i,j-1,k)-hx(i,j,k))
              ez(i,j,k)=ezx(i,j,k)+ezy(i,j,k)
!              l3=l3+1
           enddo
        enddo
     enddo

enddo

!return !!!***

end subroutine e_pml



!////////////////////////////////////////////////////////////////////////////////////////////////

!/////////////////////////////////////////////////////////////////////////////////////////////////

subroutine h_pml(ex,ey,ez,hx,hy,hz)
  use const_para
  implicit none

  integer,parameter :: npml = 5 !PMML層数 4~16
  integer,parameter :: MM = 3 !2~4
  real(8),parameter :: RRcoef = 0.01d0!-120.0d0 !反射係数 |R(0)| !R should be [10^-2, 10^-12]
  integer :: l
  integer :: i0,i1,j0,j1,k0,k1,l1,l2,l3,lpmlst(6)
  real(8) :: lpmlii(6,2),lpmljj(6,2),lpmlkk(6,2)
  real(8) :: exy(nx,ny,nz),exz(nx,ny,nz),eyx(nx,ny,nz),eyz(nx,ny,nz),ezx(nx,ny,nz),ezy(nx,ny,nz)
  real(8) :: hxy(nx,ny,nz),hxz(nx,ny,nz),hyx(nx,ny,nz),hyz(nx,ny,nz),hzx(nx,ny,nz),hzy(nx,ny,nz)
  real(8) :: cxe(nx),cye(ny),cze(nz)
  real(8) :: cxh(nx),cyh(ny),czh(nz)
  real(8) :: cxel(nx),cyel(ny),czel(nz)
  real(8) :: cxhl(nx),cyhl(ny),czhl(nz)
  real(8) :: shx,shy,shz,sex,sey,sez
  real(8) :: sig_max,sig(npml)
  complex(kind(0d0)),intent(in) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
  complex(kind(0d0)),intent(inout) :: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)

    sig_max = - (MM+1)*epsi0*cmax/2.0d0/dble(npml)/dx *log(1.0d0/RRcoef)

    do i=1,npml
    sig(i) = sig_max * (dble(i)-1.0d0/2.0d0)/dble(npml)**MM
    enddo


    do i=1,npml
      do j=1,npml
        do k=1,npml
      cxh(i) = (1.0d0 - sig(i)*dt/2.0d0/epsi0) /(1.0d0 + sig(i)*dt/2.0d0/epsi0)
      cyh(j) = (1.0d0 - sig(j)*dt/2.0d0/epsi0) /(1.0d0 + sig(j)*dt/2.0d0/epsi0)
      czh(k) = (1.0d0 - sig(k)*dt/2.0d0/epsi0) /(1.0d0 + sig(k)*dt/2.0d0/epsi0)

      cxhl(i) = - (dt/MU0) / (1.0d0 + sig(i)*dt/2.0d0/epsi0)/dx
      cyhl(j) = - (dt/MU0) / (1.0d0 + sig(j)*dt/2.0d0/epsi0)/dy
      czhl(k) = - (dt/MU0) / (1.0d0 + sig(k)*dt/2.0d0/epsi0)/dz
        enddo
      enddo
    enddo

    !i=1
   lpmlii(1,1)=1
   lpmljj(1,1)=1
   lpmlkk(1,1)=1
   lpmlii(1,2)=npml+1
   lpmljj(1,2)=ny
   lpmlkk(1,2)=nz
   lpmlii(2,1)=nx-npml
   lpmljj(2,1)=1
   lpmlkk(2,1)=1
   lpmlii(2,2)=nx
   lpmljj(2,2)=ny
   lpmlkk(2,2)=nz

   lpmlii(3,1)=1
   lpmljj(3,1)=1
   lpmlkk(3,1)=1
   lpmlii(3,2)=nx
   lpmljj(3,2)=npml+1
   lpmlkk(3,2)=nz
   lpmlii(4,1)=1
   lpmljj(4,1)=ny-npml
   lpmlkk(4,1)=1
   lpmlii(4,2)=nx
   lpmljj(4,2)=ny
   lpmlkk(4,2)=nz

   lpmlii(5,1)=1
   lpmljj(5,1)=1
   lpmlkk(5,1)=1
   lpmlii(5,2)=nx
   lpmljj(5,2)=ny
   lpmlkk(5,2)=npml+1
   lpmlii(6,1)=1
   lpmljj(6,1)=1
   lpmlkk(6,1)=nz-npml
   lpmlii(6,2)=nx
   lpmljj(6,2)=ny
   lpmlkk(6,2)=nz


  do l=1,6

     i0=lpmlii(l,1) !l=1 : i-1, l=2 : i=nxに接する層
     i1=lpmlii(l,2)
     j0=lpmljj(l,1)
     j1=lpmljj(l,2)
     k0=lpmlkk(l,1)
     k1=lpmlkk(l,2)

!     l1=lpmlst(l)
     do i=i0,i1-1
        do j=j0+1,j1-1
           do k=k0+1,k1-1
              hxy(i,j,k)=cyh(j)*hxy(i,j,k)+cyhl(j)*(ez(i,j,k)-ez(i,j-1,k))
              hxz(i,j,k)=czh(k)*hxz(i,j,k)+czhl(k)*(ey(i,j,k-1)-ey(i,j,k))
              hx(i,j,k)=hxy(i,j,k)+hxz(i,j,k)
!              l1=l1+1
           enddo
        enddo
     enddo

!     l2=lpmlst(l)
     do i=i0+1,i1-1
        do j=j0,j1-1
           do k=k0+1,k1-1
              hyx(i,j,k)=cxh(i)*hyx(i,j,k)+cxhl(i)*(ez(i-1,j,k)-ez(i,j,k))
              hyz(i,j,k)=czh(k)*hyz(i,j,k)+cyhl(k)*(ex(i,j,k)-ex(i,j,k-1))
              hy(i,j,k)=hyx(i,j,k)+hyz(i,j,k)
!              l2=l2+1
           enddo
        enddo
     enddo

!      l3=lpmlst(l)
     do i=i0+1,i1-1
        do j=j0+1,j1-1
           do k=k0,k1-1
              hzx(i,j,k)=cxh(i)*hzx(i,j,k)+cxhl(i)*(ey(i,j,k)-ey(i-1,j,k))
              hzy(i,j,k)=cyh(j)*hzy(i,j,k)+cyhl(j)*(ex(i,j-1,k)-ex(i,j,k))
              hz(i,j,k)=hzx(i,j,k)+hzy(i,j,k)
!              l3=l3+1
           enddo
        enddo
     enddo

enddo  !!!***

!return  !!!***

end subroutine h_pml
