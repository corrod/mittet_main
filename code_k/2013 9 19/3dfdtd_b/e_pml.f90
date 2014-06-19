subroutine e_pml
  use fdtd
  implicit none
  integer :: l,i,j,k
  
  do l=1,6

     i0=lpmlii(l,1)
     i1=lpmlii(l,2)
     j0=lpmljj(l,1)
     j1=lpmljj(l,2)
     k0=lpmlkk(l,1)
     k1=lpmlkk(l,2)

     do i=i0,i1-1
        do j=j0+1,j1-1
           do k=k0+1,k1-1
              exy(i,j,k)=cye(j)*exy(i,j,k)+cyel(j)*(hz(i,j,k)-hz(i,j-1,k))
              exz(i,j,k)=cze(k)*exz(i,j,k)+czel(k)*(hy(i,j,k-1)-hy(i,j,k))
              ex(i,j,k)=exy(i,j,k)+exz(i,j,k)
           end do
        end do
     end do

     do i=i0+1,i1-1
        do j=j0,j1-1
           do k=k0+1,k1-1
              eyx(i,j,k)=cxe(i)*eyx(i,j,k)+cxel(i)*(hz(i-1,j,k)-hz(i,j,k))
              eyz(i,j,k)=cze(k)*eyz(i,j,k)+czel(k)*(hx(i,j,k)-hx(i,j,k-1))
              ey(i,j,k)=eyx(i,j,k)+eyz(i,j,k)
           end do
        end do
     end do

     do i=i0+1,i1-1
        do j=j0+1,j1-1
           do k=k0,k1-1
              ezx(i,j,k)=cxe(i)*ezx(i,j,k)+cxel(i)*(hy(i,j,k)-hy(i-1,j,k))
              ezy(i,j,k)=cye(j)*ezy(i,j,k)+cyel(j)*(hx(i,j-1,k)-hx(i,j,k))
              ez(i,j,k)=ezx(i,j,k)+ezy(i,j,k)
           end do
        end do
     end do

  open(1,file='pml_model')
  write(1,*) l,i0,i1,j0,j1,k0,k1

  end do
  close(1)
end subroutine
