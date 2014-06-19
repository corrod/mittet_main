subroutine init_eh_pml_wg_2dtm
!**************************************************************************
!   subroutine initial for setting initial condition
!   path: main/init_eh_pml_wg
!**************************************************************************
   use fdtd_lib_2dtm
   implicit none

   integer i,j,k

   !******** number of lattice ********
   my=iy
   mz=iz
   ! reference point
   nypoint=4
   nzpoint=4
   !******** memory allocation ********
   ! for field 
   allocate(eyz(my,mz),ezy(my,mz), &
            hxy(my,mz),hxz(my,mz), &
            stat=code)
   if(code /= 0) stop "*** memory <pml_field> ***"

   ! for conductivity
   allocate(esigy(my),esigz(mz), &
            msigy(my),msigz(mz), &
            stat=code)
   if(code /= 0) stop "*** memory <pml_sigma> ***"

   ! for field coefficient 
   allocate(ceyz(mz), &
            cezy(my), &
            chxy(my),chxz(mz), &
            ceyzrz(mz), &
            cezyry(my), &
            chxyry(my),chxzrz(mz), &
            stat=code)
   if(code /= 0) stop "*** memory <pml_field_coeff> ***"

   ! for reference point
   allocate(zip(nzpoint),yip(nypoint), &
            stat=code)
   if(code /= 0) stop "*** memory <pml_number> ***"

   !******** set all field to zero ********
   do j=1,my
      do k=1,mz
         eyz(j,k)=0.0d0
         ezy(j,k)=0.0d0
         hxy(j,k)=0.0d0
         hxz(j,k)=0.0d0
      end do
   end do

   !******** set all conductivity to zero ********
   do j=1,my
      esigy(j)=0.0e0
      msigy(j)=0.0e0
   end do
   do k=1,mz
      esigz(k)=0.0e0
      msigz(k)=0.0e0
   end do

   !******** initial value for pml ********
   dlp=dl
   dyp=dl
   dzp=dl

   r0=1.0d-5  !(2.90)
   order=2
!   esig_max=-(real(order+1)*eps0*vc*dlog(r0))/(2.0d0*ncpml*dl)  !(2.91)
   esig_max=10.0d0
   msig_max=(mu0/eps0)*esig_max

   zip(1)=zi(1)
   zip(2)=zi(2)
   zip(3)=zi(7)
   zip(4)=zi(8)

   yip(1)=yi(1)
   yip(2)=yi(2)
   yip(3)=yi(7)
   yip(4)=yi(8)


   write(1,*) "pml" 
   write(1,'(1x, "esig_max:",d18.10)') esig_max
   write(1,'(1x, "msig_max:",d18.10)') msig_max
   write(1,'(1x, "zip(1):",i4)') zip(1)    
   write(1,'(1x, "zip(2):",i4)') zip(2)    
   write(1,'(1x, "zip(3):",i4)') zip(3)    
   write(1,'(1x, "zip(4):",i4)') zip(4) 
   write(1,'(1x, "zip(1):",i4)') yip(1)    
   write(1,'(1x, "zip(2):",i4)') yip(2)    
   write(1,'(1x, "zip(3):",i4)') yip(3)    
   write(1,'(1x, "zip(4):",i4)') yip(4) 
   write(1,*)

   !******** conductivity for z propagation ********
   ! for sigmaz (electric conductivity)
   do k=zip(1),zip(2)
      esigz(k)=esig_max*((real(ncpml)-real(k)+1.0d0)/real(ncpml))**order
   end do
   do k=zip(3),zip(4)
      esigz(k)=esig_max*real((real(k)-(real(zip(4))-real(ncpml))) & 
               /real(ncpml))**order
   end do
   ! for sigma*z (magnetic conductivity)
   do k=zip(1),zip(2)-1
      msigz(k)=msig_max*((real(ncpml)-real(k)+0.5d0)/real(ncpml))**order
   end do
   do k=zip(3),zip(4)-1
      msigz(k)=msig_max*((real(k)-(real(zip(4))-real(ncpml)-0.5d0)) & 
               /real(ncpml))**order
   end do

   open(unit=2,file='pml_conductivity')
   do k=1,iz
      write(2,*) k,esigz(k),msigz(k)
   end do
   close(2)


   !******** conductivity for y propagation ********
   ! for sigmay (electric conductivity)
   do j=yip(1),yip(2)
      esigy(j)=esig_max*((real(ncpml)-real(j)+1.0d0)/real(ncpml))**order
   end do
   do j=yip(3),yip(4)
      esigy(j)=esig_max*real((real(j)-(real(yip(4))-real(ncpml))) &
               /real(ncpml))**order
   end do
  ! for sigma*y (magnetic conductivity)
   do j=yip(1),yip(2)-1
      msigy(j)=msig_max*((real(ncpml)-real(j)+0.5d0)/real(ncpml))**order
   end do
   do j=yip(3),yip(4)-1
      msigy(j)=msig_max*((real(j)-(real(yip(4))-real(ncpml)-0.5d0)) &
               /real(ncpml))**order
   end do

   open(unit=3,file='pml_conductivity_y')
   do j=1,iy
      write(3,*) j,esigy(j),msigy(j)
   end do
   close(3)

   !******** e_pml field coefficient ********
   do j=1,my
      do k=1,mz   
         ceyz(k)=(1.0d0-((esigz(k)*dt)/(2.0d0*eps0))) &
                /(1.0d0+((esigz(k)*dt)/(2.0d0*eps0)))
         cezy(j)=(1.0d0-((esigy(j)*dt)/(2.0d0*eps0))) &
                /(1.0d0+((esigy(j)*dt)/(2.0d0*eps0)))
 
         ceyzrz(k)=(dt/eps0)/(1.0d0+((esigz(k)*dt)/(2.0d0*eps0)))/dzp
         cezyry(j)=(dt/eps0)/(1.0d0+((esigy(j)*dt)/(2.0d0*eps0)))/dyp
      end do
   end do

   !******** h_pml field coefficient ********
   do j=1,my
      do k=1,mz 
         chxy(j)=(1.0d0-((msigy(j)*dt)/(2.0d0*mu0))) &
                /(1.0d0+((msigy(j)*dt)/(2.0d0*mu0)))
         chxz(k)=(1.0d0-((msigz(k)*dt)/(2.0d0*mu0))) &
                /(1.0d0+((msigz(k)*dt)/(2.0d0*mu0)))

         chxyry(j)=(dt/mu0)/(1.0d0+((msigy(j)*dt)/(2.0d0*mu0)))/dyp
         chxzrz(k)=(dt/mu0)/(1.0d0+((msigz(k)*dt)/(2.0d0*mu0)))/dzp
      end do
   end do

   open(unit=30,file='ce_z')
   open(unit=31,file='ch_z')
   do k=1,iz
      write(30,'(1x,i4,4d18.10)') k,ceyz(k),ceyzrz(k)
      write(31,'(1x,i4,4d18.10)') k,chxz(k),chxzrz(k)
   end do
   close(30)
   close(31)

   open(unit=34,file='ce_y')
   open(unit=35,file='ch_y')
   do j=1,iy
      write(34,'(1x,i4,4d18.10)') j,cezy(j),cezyry(j)
      write(35,'(1x,i4,4d18.10)') j,chxy(j),chxyry(j)
   end do
   close(34)
   close(35)

end subroutine init_eh_pml_wg_2dtm


subroutine e_pml_wg_2dtm
!****************************************************************************
! electric field calculation by fd method
!****************************************************************************
   use fdtd_lib_2dtm
   implicit none

   integer :: i,j,k,id

!******** PML domain  ******** 
!
!   +---+---------------+---+
!   | 1 |       8       | 6 |
!   +---+---------------+---+
!   |   |               |   |
!   |   |               |   |
!   |   |               |   |
!   | 2 |               | 5 |
!   |   |               |   |
!   |   |               |   |
!   |   |               |   |
!   +---+---------------+---+
!   | 3 |       7       | 4 |
!   +---+---------------+---+


   ! rectangular pml domain 1
   pmlzs=zip(1)
   pmlze=zip(2)
   pmlys=yip(1)
   pmlye=yip(2)
   call ey_pml_calc
   call ez_pml_calc

   ! rectangular pml domain 2
   pmlzs=zip(1)
   pmlze=zip(2)
   pmlys=yip(2)
   pmlye=yip(3)
   call ey_pml_calc

   ! rectangular pml domain 3
   pmlzs=zip(1)
   pmlze=zip(2)
   pmlys=yip(3)
   pmlye=yip(4)
   call ey_pml_calc
   call ez_pml_calc


   ! rectangular pml domain 4
   pmlzs=zip(3)
   pmlze=zip(4)
   pmlys=yip(1)
   pmlye=yip(2)
   call ey_pml_calc
   call ez_pml_calc

   ! rectangular pml domain 5
   pmlzs=zip(3)
   pmlze=zip(4)
   pmlys=yip(2)
   pmlye=yip(3)
   call ey_pml_calc

   ! rectangular pml domain 6
   pmlzs=zip(3)
   pmlze=zip(4)
   pmlys=yip(3)
   pmlye=yip(4)
   call ey_pml_calc
   call ez_pml_calc

   ! rectangular pml domain 7
   pmlzs=zip(2)
   pmlze=zip(3)
   pmlys=yip(1)
   pmlye=yip(2)
   call ez_pml_calc
!   call ey_pml_calc

   ! rectangular pml domain 8
   pmlzs=zip(2)
   pmlze=zip(3)
   pmlys=yip(3)
   pmlye=yip(4)
   call ez_pml_calc
!   call ey_pml_calc

end subroutine e_pml_wg_2dtm


subroutine h_pml_wg_2dtm
!****************************************************************************
! magnetic field calculation by fd method
!****************************************************************************
   use fdtd_lib_2dtm
   implicit none

   integer :: i,j,k,id

   ! rectangular pml domain 1
   pmlzs=zip(1)
   pmlze=zip(2)
   pmlys=yip(1)
   pmlye=yip(4)
   call h_pml_calc

   ! rectangular pml domain 2
   pmlzs=zip(3)
   pmlze=zip(4)
   pmlys=yip(1)
   pmlye=yip(4)
   call h_pml_calc

   ! rectangular pml domain 3
   pmlzs=zip(2)
   pmlze=zip(3)
   pmlys=yip(1)
   pmlye=yip(2)
   call h_pml_calc

   ! rectangular pml domain 4
   pmlzs=zip(2)
   pmlze=zip(3)
   pmlys=yip(3)
   pmlye=yip(4)
   call h_pml_calc


end subroutine h_pml_wg_2dtm



subroutine ey_pml_calc
!****************************************************************************
! electric field calculation in pml by fd method
!****************************************************************************
   use fdtd_lib_2dtm
   implicit none

   integer :: i,j,k,id

   !******** ey ********
   do k=pmlzs+1,pmlze-1
      do j=pmlys,pmlye-1

            id=id_ey(j,k)

            if(id.eq.0) then
            ! 0: vacume
               eyz(j,k)=ceyz(k)*eyz(j,k) &
                       +ceyzrz(k)*(hx(j,k)-hx(j,k-1))
            else
            ! 1: ideal conductor
               eyz(j,k)=0.0e0
            end if

               ey(j,k)=eyz(j,k)

      end do
   end do

end subroutine ey_pml_calc


subroutine ez_pml_calc
!****************************************************************************
! electric field calculation in pml by fd method
!****************************************************************************
   use fdtd_lib_2dtm
   implicit none

   integer :: i,j,k,id

   !******** ez ********
   do k=pmlzs,pmlze-1
      do j=pmlys+1,pmlye-1

            id=id_ez(j,k)

            if(id.eq.0) then
            ! 0: vacume
               ezy(j,k)=cezy(j)*ezy(j,k) &
                       -cezyry(j)*(hx(j,k)-hx(j-1,k))
            else
            ! 1: ideal conductor
               ezy(j,k)=0.0e0
            end if

               ez(j,k)=ezy(j,k)

      end do
   end do

end subroutine ez_pml_calc


subroutine h_pml_calc
!****************************************************************************
! electric field calculation in pml by fd method
!****************************************************************************
   use fdtd_lib_2dtm
   implicit none

   integer :: i,j,k,id

   !******** hx ********
   do k=pmlzs,pmlze-1
      do j=pmlys,pmlye-1

            id=id_hx(j,k)

            if(id.eq.0) then
            ! 0: vacume
               hxy(j,k)=chxy(j)*hxy(j,k) &
                       -chxyry(j)*(ez(j+1,k)-ez(j,k))
               hxz(j,k)=chxz(k)*hxz(j,k) &
                       +chxzrz(k)*(ey(j,k+1)-ey(j,k))
            else
            ! 1: magnetic ideal conductor
               hxy(j,k)=0.0e0
               hxz(j,k)=0.0e0
            end if

               hx(j,k)=hxy(j,k)+hxz(j,k)

      end do
   end do

end subroutine h_pml_calc


