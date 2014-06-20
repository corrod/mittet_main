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

!     l1=lpmlst(l)
     do i=i0,i1-1
        do j=j0+1,j1-1
           do k=k0+1,k1-1
              exy(i,j,k)=cye(j)*exy(i,j,k)+cyel(j)*(hz(i,j,k)-hz(i,j-1,k))
              exz(i,j,k)=cze(k)*exz(i,j,k)+czel(k)*(hy(i,j,k-1)-hy(i,j,k))
              ex(i,j,k)=exy(i,j,k)+exz(i,j,k)
!              l1=l1+1
           end do
        end do
     end do

!     l2=lpmlst(l)
     do i=i0+1,i1-1
        do j=j0,j1-1
           do k=k0+1,k1-1
              eyx(i,j,k)=cxe(i)*eyx(i,j,k)+cxel(i)*(hz(i-1,j,k)-hz(i,j,k))
              eyz(i,j,k)=cze(k)*eyz(i,j,k)+cyel(k)*(hx(i,j,k)-hx(i,j,k-1))
              ey(i,j,k)=eyx(i,j,k)+eyz(i,j,k)
!              l2=l2+1
           end do
        end do
     end do

!      l3=lpmlst(l)
     do i=i0+1,i1-1
        do j=j0+1,j1-1
           do k=k0,k1-1
              ezx(i,j,k)=cxe(i)*ezx(i,j,k)+cxel(i)*(hy(i,j,k)-hy(i-1,j,k))
              ezy(i,j,k)=cye(j)*ezy(i,j,k)+cyel(j)*(hx(i,j-1,k)-hx(i,j,k))
              ez(i,j,k)=ezx(i,j,k)+ezy(i,j,k)
!              l3=l3+1
           end do
        end do
     end do
  end do

return

end subroutine
