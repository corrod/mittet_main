
subroutine e_pml_3d
!****************************************************************************
! electric field calculation by fd method
!****************************************************************************
   use fdtd
   implicit none

   integer :: i,j,k,id


   ! rectangular pml domain 1
   pmlzs=1
   pmlze=lpml
   pmlys=1
   pmlye=lpml
   pmlxs=1
   pmlxe=lpml
   call ex_pml_calc
   call ey_pml_calc
   call ez_pml_calc

   ! rectangular pml domain 2
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=1
   pmlye=lpml
   pmlxs=1
   pmlxe=lpml
   call ey_pml_calc
   call ex_pml_calc
   call ez_pml_calc

   ! rectangular pml domain 3
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=1
   pmlxe=lpml
   call ex_pml_calc
   call ey_pml_calc
   call ez_pml_calc

   ! rectangular pml domain 4
   pmlzs=1
   pmlze=lpml
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=1
   pmlxe=lpml
   call ex_pml_calc
   call ey_pml_calc
   call ez_pml_calc

   ! rectangular pml domain 5
   pmlzs=1
   pmlze=lpml
   pmlys=1
   pmlye=lpml
   pmlxs=nx-lpml
   pmlxe=nx
   call ex_pml_calc
   call ey_pml_calc
   call ez_pml_calc

   ! rectangular pml domain 6
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=1
   pmlye=lpml
   pmlxs=nx-lpml
   pmlxe=nx
   call ex_pml_calc
   call ey_pml_calc
   call ez_pml_calc

   ! rectangular pml domain 7
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=nx-lpml
   pmlxe=nx
   call ez_pml_calc
   call ey_pml_calc
   call ex_pml_calc

   ! rectangular pml domain 8
   pmlzs=1
   pmlze=lpml
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=nx-lpml
   pmlxe=nx
   call ex_pml_calc
   call ez_pml_calc
   call ey_pml_calc

   ! rectangular pml domain 9
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=1
   pmlye=lpml
   pmlxs=1
   pmlxe=lpml
   call ex_pml_calc
    call ez_pml_calc
   call ey_pml_calc

   ! rectangular pml domain 10
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=1
   pmlxe=lpml
   call ex_pml_calc
   call ez_pml_calc
    call ey_pml_calc
 
   ! rectangular pml domain 11
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=1
   pmlxe=lpml
   call ex_pml_calc
    call ez_pml_calc
   call ey_pml_calc
 
   ! rectangular pml domain 12
   pmlzs=1
   pmlze=lpml
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=1
   pmlxe=lpml
   call ex_pml_calc
   call ez_pml_calc
    call ey_pml_calc

   ! rectangular pml domain 13
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=1
   pmlye=lpml
   pmlxs=nx-lpml
   pmlxe=nx
   call ex_pml_calc
    call ez_pml_calc
   call ey_pml_calc
 
   ! rectangular pml domain 14
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=nx-lpml
   pmlxe=nx
   call ex_pml_calc
   call ez_pml_calc
    call ey_pml_calc

   ! rectangular pml domain 15
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=nx-lpml
   pmlxe=nx
   call ex_pml_calc
    call ez_pml_calc
   call ey_pml_calc

   ! rectangular pml domain 16
   pmlzs=1
   pmlze=lpml
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=nx-lpml
   pmlxe=nx
   call ex_pml_calc
   call ez_pml_calc
    call ey_pml_calc

   ! rectangular pml domain 17
   pmlzs=1
   pmlze=lpml
   pmlys=1
   pmlye=lpml
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call ex_pml_calc
   call ez_pml_calc
   call ey_pml_calc

   ! rectangular pml domain 18
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=1
   pmlye=lpml
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call ex_pml_calc
   call ez_pml_calc
   call ey_pml_calc

   ! rectangular pml domain 19
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call ex_pml_calc
   call ez_pml_calc
   call ey_pml_calc

   ! rectangular pml domain 20
   pmlzs=1
   pmlze=lpml
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call ex_pml_calc
   call ez_pml_calc
   call ey_pml_calc

   ! rectangular pml domain 21
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=1
   pmlye=lpml
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call ex_pml_calc
    call ez_pml_calc
   call ey_pml_calc

   ! rectangular pml domain 22
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call ex_pml_calc
   call ez_pml_calc
    call ey_pml_calc

   ! rectangular pml domain 23
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call ex_pml_calc
    call ez_pml_calc
   call ey_pml_calc

   ! rectangular pml domain 24
   pmlzs=1
   pmlze=lpml
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call ex_pml_calc
   call ez_pml_calc
    call ey_pml_calc

   ! rectangular pml domain 25
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=1
   pmlxe=lpml
   call ex_pml_calc
    call ez_pml_calc
    call ey_pml_calc

   ! rectangular pml domain 26
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=nx-lpml
   pmlxe=nx
   call ex_pml_calc
    call ez_pml_calc
    call ey_pml_calc  

end subroutine e_pml_3d

subroutine h_pml_3d
!****************************************************************************
! electric field calculation by fd method
!****************************************************************************
   use fdtd
   implicit none

   integer :: i,j,k,id


   ! rectangular pml domain 1
   pmlzs=1
   pmlze=lpml
   pmlys=1
   pmlye=lpml
   pmlxs=1
   pmlxe=lpml
   call hx_pml_calc
   call hy_pml_calc
   call hz_pml_calc

   ! rectangular pml domain 2
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=1
   pmlye=lpml
   pmlxs=1
   pmlxe=lpml
   call hy_pml_calc
   call hx_pml_calc
   call hz_pml_calc

   ! rectangular pml domain 3
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=1
   pmlxe=lpml
   call hx_pml_calc
   call hy_pml_calc
   call hz_pml_calc

   ! rectangular pml domain 4
   pmlzs=1
   pmlze=lpml
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=1
   pmlxe=lpml
   call hx_pml_calc
   call hy_pml_calc
   call hz_pml_calc

   ! rectangular pml domain 5
   pmlzs=1
   pmlze=lpml
   pmlys=1
   pmlye=lpml
   pmlxs=nx-lpml
   pmlxe=nx
   call hx_pml_calc
   call hy_pml_calc
   call hz_pml_calc

   ! rectangular pml domain 6
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=1
   pmlye=lpml
   pmlxs=nx-lpml
   pmlxe=nx
   call hx_pml_calc
   call hy_pml_calc
   call hz_pml_calc

   ! rectangular pml domain 7
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=nx-lpml
   pmlxe=nx
   call hz_pml_calc
   call hy_pml_calc
   call hx_pml_calc

   ! rectangular pml domain 8
   pmlzs=1
   pmlze=lpml
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=nx-lpml
   pmlxe=nx
   call hx_pml_calc
   call hz_pml_calc
   call hy_pml_calc

   ! rectangular pml domain 9
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=1
   pmlye=lpml
   pmlxs=1
   pmlxe=lpml
   call hx_pml_calc
    call hz_pml_calc
   call hy_pml_calc

   ! rectangular pml domain 10
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=1
   pmlxe=lpml
   call hx_pml_calc
   call hz_pml_calc
    call hy_pml_calc

   ! rectangular pml domain 11
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=1
   pmlxe=lpml
   call hx_pml_calc
    call hz_pml_calc
   call hy_pml_calc
 
   ! rectangular pml domain 12
   pmlzs=1
   pmlze=lpml
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=1
   pmlxe=lpml
   call hx_pml_calc
   call hz_pml_calc
    call hy_pml_calc

   ! rectangular pml domain 13
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=1
   pmlye=lpml
   pmlxs=nx-lpml
   pmlxe=nx
   call hx_pml_calc
    call hz_pml_calc
   call hy_pml_calc
 
   ! rectangular pml domain 14
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=nx-lpml
   pmlxe=nx
   call hx_pml_calc
   call hz_pml_calc
    call hy_pml_calc

   ! rectangular pml domain 15
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=nx-lpml
   pmlxe=nx
   call hx_pml_calc
    call hz_pml_calc
   call hy_pml_calc

   ! rectangular pml domain 16
   pmlzs=1
   pmlze=lpml
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=nx-lpml
   pmlxe=nx
   call hx_pml_calc
   call hz_pml_calc
    call hy_pml_calc

   ! rectangular pml domain 17
   pmlzs=1
   pmlze=lpml
   pmlys=1
   pmlye=lpml
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call hx_pml_calc
   call hz_pml_calc
   call hy_pml_calc

   ! rectangular pml domain 18
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=1
   pmlye=lpml
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call hx_pml_calc
   call hz_pml_calc
   call hy_pml_calc

   ! rectangular pml domain 19
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call hx_pml_calc
   call hz_pml_calc
   call hy_pml_calc

   ! rectangular pml domain 20
   pmlzs=1
   pmlze=lpml
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call hx_pml_calc
   call hz_pml_calc
   call hy_pml_calc

   ! rectangular pml domain 21
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=1
   pmlye=lpml
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call hx_pml_calc
    call hz_pml_calc
   call hy_pml_calc

   ! rectangular pml domain 22
   pmlzs=nz-lpml
   pmlze=nz
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call hx_pml_calc
   call hz_pml_calc
    call hy_pml_calc

   ! rectangular pml domain 23
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=ny-lpml
   pmlye=ny
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call hx_pml_calc
    call hz_pml_calc
   call hy_pml_calc

   ! rectangular pml domain 24
   pmlzs=1
   pmlze=lpml
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=lpml+1
   pmlxe=nx-lpml-1
    call hx_pml_calc
   call hz_pml_calc
    call hy_pml_calc

   ! rectangular pml domain 25
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=1
   pmlxe=lpml
   call hx_pml_calc
    call hz_pml_calc
    call hy_pml_calc

   ! rectangular pml domain 26
   pmlzs=lpml+1
   pmlze=nz-lpml-1
   pmlys=lpml+1
   pmlye=ny-lpml-1
   pmlxs=nx-lpml
   pmlxe=nx
   call hx_pml_calc
    call hz_pml_calc
    call hy_pml_calc  

end subroutine h_pml_3d




subroutine ex_pml_calc
!****************************************************************************
! electric field calculation in pml by fd method
!****************************************************************************
   use fdtd
   implicit none

   integer :: i,j,k,id

   !******** ex ********
   do k=pmlzs+1,pmlze-1
      do j=pmlys+1,pmlye-1
         do i=pmlxs,pmlxe-1

            id=media_id(i,j,k)

            if(id.eq.0) then
            ! 0: vacume
               exy(i,j,k)=cye(j)*exy(i,j,k)+cyel(j)*(hz(i,j,k)-hz(i,j-1,k))
               exz(i,j,k)=cze(k)*exz(i,j,k)+czel(k)*(hy(i,j,k-1)-hy(i,j,k))
            else
            ! 1: ideal conductor
               exy(i,j,k)=0.0e0
               exz(i,j,k)=0.0e0
            end if

               ex(i,j,k)=exy(i,j,k)+exz(i,j,k)
            end do
      end do
   end do

end subroutine ex_pml_calc


subroutine ey_pml_calc
!****************************************************************************
! electric field calculation in pml by fd method
!****************************************************************************
   use fdtd
   implicit none

   integer :: i,j,k,id

   !******** ey ********
   do k=pmlzs+1,pmlze-1
      do j=pmlys,pmlye-1
         do i=pmlxs+1,pmlxe-1

            id=media_id(i,j,k)

            if(id.eq.0) then
            ! 0: vacume
               eyx(i,j,k)=cxe(i)*eyx(i,j,k)+cxel(i)*(hz(i-1,j,k)-hz(i,j,k))
               eyz(i,j,k)=cze(k)*eyz(i,j,k)+czel(k)*(hx(i,j,k)-hx(i,j,k-1))
            else
            ! 1: ideal conductor
               eyx(i,j,k)=0.0e0
               eyz(i,j,k)=0.0e0
            end if

               ey(i,j,k)=eyx(i,j,k)+eyz(i,j,k)
            end do
      end do
   end do

end subroutine ey_pml_calc


subroutine ez_pml_calc
!****************************************************************************
! electric field calculation in pml by fd method
!****************************************************************************
   use fdtd
   implicit none

   integer :: i,j,k,id

   !******** ez ********
   do k=pmlzs,pmlze-1
      do j=pmlys+1,pmlye-1
         do i=pmlxs+1,pmlxe-1

            id=media_id(i,j,k)

            if(id.eq.0) then
            ! 0: vacume
               ezx(i,j,k)=cxe(i)*ezx(i,j,k)+cxel(i)*(hy(i,j,k)-hy(i-1,j,k))
               ezy(i,j,k)=cye(j)*ezy(i,j,k)+cyel(j)*(hx(i,j-1,k)-hx(i,j,k))
            else
            ! 1: ideal conductor
               ezy(i,j,k)=0.0e0
               ezx(i,j,k)=0.0e0
               
            end if

               ez(i,j,k)=ezy(i,j,k)+ezx(i,j,k)
           end do
      end do
   end do

end subroutine ez_pml_calc


subroutine hx_pml_calc
!****************************************************************************
! electric field calculation in pml by fd method
!****************************************************************************
   use fdtd
   implicit none

   integer :: i,j,k,id

   !******** hx ********
   do k=pmlzs,pmlze-1
      do j=pmlys,pmlye-1
         do i=pmlxs,pmlxe-1

            id=media_id(i,j,k)

            if(id.eq.0) then
            ! 0: vacume
               hxy(i,j,k)=cyh(j)*hxy(i,j,k)-cyhl(j)*(ez(i,j+1,k)-ez(i,j,k))
               hxz(i,j,k)=czh(k)*hxz(i,j,k)+czhl(k)*(ey(i,j,k+1)-ey(i,j,k))
            else
            ! 1: magnetic ideal conductor
               hxy(i,j,k)=0.0e0
               hxz(i,j,k)=0.0e0
            end if

               hx(i,j,k)=hxy(i,j,k)+hxz(i,j,k)
            end do
      end do
   end do

end subroutine hx_pml_calc


subroutine hy_pml_calc
!****************************************************************************
! electric field calculation in pml by fd method
!****************************************************************************
   use fdtd
   implicit none

   integer :: i,j,k,id

   !******** hy ********
   do k=pmlzs,pmlze-1
      do j=pmlys,pmlye-1
         do i=pmlxs,pmlxe-1

            id=media_id(i,j,k)

            if(id.eq.0) then
            ! 0: vacume
               hyx(i,j,k)=cxh(i)*hyx(i,j,k)-cxhl(i)*(ez(i+1,j,k)-ez(i,j,k))
               hyz(i,j,k)=czh(k)*hyz(i,j,k)-czhl(k)*(ex(i,j,k)-ex(i,j,k+1))
            else
            ! 1: magnetic ideal conductor
               hyx(i,j,k)=0.0e0
               hyz(i,j,k)=0.0e0
            end if

               hy(i,j,k)=hyx(i,j,k)+hyz(i,j,k)
           end do
      end do
   end do

end subroutine hy_pml_calc

subroutine hz_pml_calc
!****************************************************************************
! electric field calculation in pml by fd method
!****************************************************************************
   use fdtd
   implicit none

   integer :: i,j,k,id

   !******** hz ********
   do k=pmlzs,pmlze-1
      do j=pmlys,pmlye-1
         do i=pmlxs,pmlxe-1

            id=media_id(i,j,k)

            if(id.eq.0) then
            ! 0: vacume
               hzx(i,j,k)=cxh(i)*hzx(i,j,k)-cxhl(i)*(ey(i+1,j,k)-ez(i,j,k))
               hzy(i,j,k)=cyh(j)*hzy(i,j,k)-cyhl(j)*(ex(i,j,k)-ey(i,j+1,k))
            else
            ! 1: magnetic ideal conductor
               hzx(i,j,k)=0.0e0
               hzy(i,j,k)=0.0e0
            end if

               hz(i,j,k)=hzx(i,j,k)+hzy(i,j,k)
            end do
      end do
   end do

end subroutine hz_pml_calc