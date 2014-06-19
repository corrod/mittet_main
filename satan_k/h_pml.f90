subroutine h_pml
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
              hxy(i,j,k)=cyh(j)*hxy(i,j,k)+cyhl(j)*(ez(i,j,k)-ez(i,j-1,k))
              hxz(i,j,k)=czh(k)*hxz(i,j,k)+czhl(k)*(ey(i,j,k-1)-ey(i,j,k))
              hx(i,j,k)=hxy(i,j,k)+hxz(i,j,k)
!              l1=l1+1
           end do
        end do
     end do

!     l2=lpmlst(l)
     do i=i0+1,i1-1
        do j=j0,j1-1
           do k=k0+1,k1-1
              hyx(i,j,k)=cxh(i)*hyx(i,j,k)+cxhl(i)*(ez(i-1,j,k)-ez(i,j,k))
              hyz(i,j,k)=czh(k)*hyz(i,j,k)+cyhl(k)*(ex(i,j,k)-ex(i,j,k-1))
              hy(i,j,k)=hyx(i,j,k)+hyz(i,j,k)
!              l2=l2+1
           end do
        end do
     end do

!      l3=lpmlst(l)
     do i=i0+1,i1-1
        do j=j0+1,j1-1
           do k=k0,k1-1
              hzx(i,j,k)=cxh(i)*hzx(i,j,k)+cxhl(i)*(ey(i,j,k)-ey(i-1,j,k))
              hzy(i,j,k)=cyh(j)*hzy(i,j,k)+cyhl(j)*(ex(i,j-1,k)-ex(i,j,k))
              hz(i,j,k)=hzx(i,j,k)+hzy(i,j,k)
!              l3=l3+1
           end do
        end do
     end do
  end do

end subroutine
