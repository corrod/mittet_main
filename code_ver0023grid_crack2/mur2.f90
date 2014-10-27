!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Murの吸収境界
!
!eyx2 eyx1を引数にいれないといけないかも
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!(x軸方向に伝播i=1,nx) 486 ok  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mur_yz(ex,ey,ez)
	use const_para
	implicit none


! 	real(8) :: cxd,cxu,cxx
! 	real(8) :: cxfyd,cxfzd
! 	complex(kind(0d0)) :: eyx1(nx,ny,nz),eyx2(nx,ny,nz),ezx1(nx,ny,nz),ezx2(nx,ny,nz)
! 	complex(kind(0d0)) :: exy1(nx,ny,nz),exy2(nx,ny,nz),ezy1(nx,ny,nz),ezy2(nx,ny,nz)
! 	complex(kind(0d0)) :: exz1(nx,ny,nz),exz2(nx,ny,nz),eyz1(nx,ny,nz),eyz2(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)

	cxd = (cmax*dt-dx)/(cmax*dt+dx)
	cxu = (cmax*dt-dx)/(cmax*dt+dx)
	cxx = (2.0d0*dx)/(cmax*dt+dx)
	cxfyd = (dx*(cmax*dt)**2.0d0) / (2.0d0*(dy*dy)*(cmax*dt+dx))
	cxfzd = (dx*(cmax+dt)**2.0d0) / (2.0d0*(dz*dz)*(cmax*dt+dx))


!-----------eyに対して未-----------------------------------------------------------------
	!1次吸収境界条件
	do k=2,nz-1
		j=1
		ey(1,j,k) = eyx1(2,j,k)+cxd*(ey(2,j,k)-eyx1(1,j,k))
		ey(nx,j,k)= eyx1(3,j,k)+cxu*(ey(nx-1,j,k)-eyx1(4,j,k))

		j=ny-1
		ey(1,j,k) = eyx1(2,j,k)+cxd*(ey(2,j,k)-eyx1(1,j,k))
		ey(nx,j,k)= eyx1(3,j,k)+cxu*(ey(nx-1,j,k)-eyx1(4,j,k))
	enddo

	do j=2,ny-2
		k=2
		ey(1,j,k) = eyx1(2,j,k)+cxd*(ey(2,j,k)-eyx1(1,j,k))
		ey(nx,j,k)= eyx1(3,j,k)+cxu*(ey(nx-1,j,k)-eyx1(4,j,k))

		k=nz-1
		ey(1,j,k) = eyx1(2,j,k)+cxd*(ey(2,j,k)-eyx1(1,j,k))
		ey(nx,j,k)= eyx1(3,j,k)+cxu*(ey(nx-1,j,k)-eyx1(4,j,k))
	enddo

	!2次吸収境界条件
	do k=3,nz-2
		do j=2,ny-2
			ey(1,j,k) = -eyx2(2,j,k)+cxd*(ey(2,j,k)+eyx2(1,j,k))&
						+cxx*(eyx1(1,j,k)+eyx1(2,j,k))&
						+cxfyd*(eyx1(1,j+1,k)-2.0d0*eyx1(1,j,k)+eyx1(1,j-1,k)&
							   +eyx1(2,j+1,k)-2.0d0*eyx1(2,j,k)+eyx1(2,j-1,k))&
						+cxfzd*(eyx1(1,j,k+1)-2.0d0*eyx1(1,j,k)+eyx1(1,j,k-1)&
							   +eyx1(2,j,k+1)-2.0d0*eyx1(2,j,k)+eyx1(2,j,k-1))
			ey(nx,j,k) = -eyx2(3,j,k)+cxd*(ey(nx-1,j,k)+eyx2(4,j,k))&
						+cxx*(eyx1(4,j,k)+eyx1(3,j,k))&
						+cxfyd*(eyx1(4,j+1,k)-2.0d0*eyx1(4,j,k)+eyx1(4,j-1,k)&
						       +eyx1(3,j+1,k)-2.0d0*eyx1(3,j,k)+eyx1(3,j-1,k))&
						+cxfzd*(eyx1(4,j,k+1)-2.0d0*eyx1(4,j,k)+eyx1(4,j,k-1)&
							   +eyx1(3,j,k+1)-2.0d0*eyx1(3,j,k)+eyx1(3,j,k-1))
		enddo
	enddo

	!過去の値の更新
	do k=2,nz-1
		do j=1,ny-1
			eyx2(1,j,k)=eyx1(1,j,k)
			eyx2(2,j,k)=eyx1(2,j,k)
			eyx2(3,j,k)=eyx1(3,j,k)
			eyx2(4,j,k)=eyx1(4,j,k)

			eyx1(1,j,k)=ey(1,j,k)
			eyx1(2,j,k)=ey(2,j,k)
			eyx1(3,j,k)=ey(nx-1,j,k)
			eyx1(4,j,k)=ey(nx,j,k)
		enddo
	enddo

!---------ezに対して-----------------------------------------------------------------
	!1次吸収境界条件
	do k=1,nz-1
		j=2
		ez(1,j,k) = ezx1(2,j,k)+cxd*(ez(2,j,k)-ezx1(1,j,k))
		ez(nx,j,k)= ezx1(3,j,k)+cxu*(ez(nx-1,j,k)-ezx1(4,j,k))

		j=ny-1
		ez(1,j,k) = ezx1(2,j,k)+cxd*(ez(2,j,k)-ezx1(1,j,k))
		ez(nx,j,k)= ezx1(3,j,k)+cxu*(ez(nx-1,j,k)-ezx1(4,j,k))
	enddo

	do j=3,ny-2
		k=1
		ez(1,j,k) = ezx1(2,j,k)+cxd*(ez(2,j,k)-ezx1(1,j,k))
		ez(nx,j,k)= ezx1(3,j,k)+cxu*(ez(nx-1,j,k)-ezx1(4,j,k))

		k=nz-1
		ez(1,j,k) = ezx1(2,j,k)+cxd*(ez(2,j,k)-ezx1(1,j,k))
		ez(nx,j,k)= ezx1(3,j,k)+cxu*(ez(nx-1,j,k)-ezx1(4,j,k))
	enddo

	!2次吸収境界条件
	do k=2,nz-2
		do j=3,ny-2
			ez(1,j,k)  = -ezx2(2,j,k)+cxd*(ez(2,j,k)+ezx2(1,j,k))&
						+cxx*(ezx1(1,j,k)+ezx1(2,j,k))&
						+cxfyd*(ezx1(1,j+1,k)-2.0d0*ezx1(1,j,k)+ezx1(1,j-1,k)&
					    	   +ezx1(2,j+1,k)-2.0d0*ezx1(2,j,k)+ezx1(2,j-1,k))&
						+cxfzd*(ezx1(1,j,k+1)-2.0d0*ezx1(1,j,k)+ezx1(1,j,k-1)&
							   +ezx1(2,j,k+1)-2.0d0*ezx1(2,j,k)+ezx1(2,j,k-1))
			ez(nx,j,k) = -ezx2(3,j,k)+cxd*(ez(nx-1,j,k)+ezx2(4,j,k))&
						+cxx*(ezx1(4,j,k)+ezx1(3,j,k))&
						+cxfyd*(ezx1(4,j+1,k)-2.0d0*ezx1(4,j,k)+ezx1(4,j-1,k)&
							   +ezx1(3,j+1,k)-2.0d0*ezx1(3,j,k)+ezx1(3,j-1,k))&
						+cxfzd*(ezx1(4,j,k+1)-2.0d0*ezx1(4,j,k)+ezx1(4,j,k-1)&
							   +ezx1(3,j,k+1)-2.0d0*ezx1(3,j,k)+ezx1(3,j,k-1))
		enddo
	enddo

	!過去の値の更新
	do k=1,nz-1
		do j=2,ny-1
			ezx2(1,j,k)=ezx1(1,j,k)
			ezx2(2,j,k)=ezx1(2,j,k)
			ezx2(3,j,k)=ezx1(3,j,k)
			ezx2(4,j,k)=ezx1(4,j,k)

			ezx1(1,j,k)=ez(1,j,k)
			ezx1(2,j,k)=ez(2,j,k)
			ezx1(3,j,k)=ez(nx-1,j,k)
			ezx1(4,j,k)=ez(nx,j,k)
		enddo
	enddo

	end subroutine mur_yz






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!(y軸方向に伝播) 未  引数を合わせていない do範囲!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine mur_zx(ex,ey,ez)
		use const_para
	implicit none
! 	real(8), intent(in) :: cmax

! 	real(8) :: cyd,cyu,cyy
! 	real(8) :: cyfxd,cyfzd
! 	complex(kind(0d0)) :: eyx1(nx,ny,nz),eyx2(nx,ny,nz),ezx1(nx,ny,nz),ezx2(nx,ny,nz)
! 	complex(kind(0d0)) :: exy1(nx,ny,nz),exy2(nx,ny,nz),ezy1(nx,ny,nz),ezy2(nx,ny,nz)
! 	complex(kind(0d0)) :: exz1(nx,ny,nz),exz2(nx,ny,nz),eyz1(nx,ny,nz),eyz2(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)

	cyd = (cmax*dt-dy)/(cmax*dt+dy)
	cyu = (cmax*dt-dy)/(cmax*dt+dy)
	cyy = (2.0d0*dy)/(cmax*dt+dy)
	cyfxd = (dy*(cmax*dt)**2.0d0) / (2.0d0*(dx*dx)*(cmax*dt+dy))
	cyfzd = (dy*(cmax+dt)**2.0d0) / (2.0d0*(dz*dz)*(cmax*dt+dy))


!--------exに対して-----------------------------------------------------------------
	!1次吸収境界条件
	do k=2,nz-1
		i=1
		ex(i,1,k) = exy1(i,2,k)+cyd*(ex(i,2,k)-exy1(i,1,k))
		ex(i,ny,k)= exy1(i,3,k)+cyu*(ex(i,ny-1,k)-exy1(i,4,k))

		i=nx-1
		ex(i,1,k) = exy1(i,2,k)+cyd*(ex(i,2,k)-exy1(i,1,k))
		ex(i,ny,k)= exy1(i,3,k)+cyu*(ex(i,ny-1,k)-exy1(i,4,k))
	enddo

	do i=2,nx-2
		k=2
		ex(i,1,k) = exy1(i,2,k)+cyd*(ex(i,2,k)-exy1(i,1,k))
		ex(i,ny,k) = exy1(i,3,k)+cyu*(ex(i,ny-1,k)-exy1(i,4,k))

		k=nz-1
		ex(i,1,k) = exy1(i,2,k)+cyd*(ex(i,2,k)-exy1(i,1,k))
		ex(i,ny,k)= exy1(i,3,k)+cyu*(ex(i,ny-1,k)-exy1(i,4,k))
	enddo

	!2次吸収境界条件
	do k=3,nz-2
		do i=2,nx-2
			ex(i,1,k)  = -exy2(i,2,k)+cyd*(ex(i,2,k)+exy2(i,1,k))&
						+cyy*(exy1(i,1,k)+exy1(i,2,k))&
						+cyfxd*(exy1(i+1,1,k)-2.0d0*exy1(i,1,k)+exy1(i-1,1,k)&
							   +exy1(i+1,2,k)-2.0d0*exy1(i,2,k)+exy1(i-1,2,k))&
						+cyfzd*(exy1(i,1,k+1)-2.0d0*exy1(i,1,k)+exy1(i,1,k-1)&
							   +exy1(i,2,k+1)-2.0d0*exy1(i,2,k)+exy1(i,2,k-1))
			ex(i,ny,k) = -exy2(i,3,k)+cyd*(ex(i,ny-1,k)+exy2(i,4,k))&
						+cyy*(exy1(i,4,k)+exy1(i,3,k))&
						+cyfxd*(exy1(i+1,4,k)-2.0d0*exy1(i,4,k)+exy1(i-1,4,k)&
							   +exy1(i+1,3,k)-2.0d0*exy1(i,3,k)+exy1(i-1,3,k))&
						+cyfzd*(exy1(i,4,k+1)-2.0d0*exy1(i,4,k)+exy1(i,4,k-1)&
						       +exy1(i,3,k+1)-2.0d0*exy1(i,3,k)+exy1(i,3,k-1))
		enddo
	enddo

	!過去の値の更新
	do k=2,nz-1
		do i=1,nx-1
			exy2(i,1,k)=exy1(i,1,k)
			exy2(i,2,k)=exy1(i,2,k)
			exy2(i,3,k)=exy1(i,3,k)
			exy2(i,4,k)=exy1(i,4,k)

			exy1(i,1,k)=ex(i,1,k)
			exy1(i,2,k)=ex(i,2,k)
			exy1(i,3,k)=ex(i,ny-1,k)
			exy1(i,4,k)=ex(i,ny,k)
		enddo
	enddo

!---------ezに対して-----------------------------------------------------------------
	!1次吸収境界条件

	do k=1,nz-1
		i=2
		ez(i,1,k) = ezy1(i,2,k)+cyd*(ez(i,2,k)-ezy1(i,1,k))
		ez(i,ny,k)= ezy1(i,3,k)+cyu*(ez(i,ny-1,k)-ezy1(i,4,k))

		i=nx-1
		ez(i,1,k) = ezy1(i,2,k)+cyd*(ez(i,2,k)-ezy1(i,1,k))
		ez(i,ny,k)= ezy1(i,3,k)+cyu*(ez(i,ny-1,k)-ezy1(i,4,k))
	enddo

	do i=3,nx-2
		k=1
		ez(i,1,k) = ezy1(i,2,k)+cyd*(ez(i,2,k)-ezy1(i,1,k))
		ez(i,ny,k)= ezy1(i,3,k)+cyu*(ez(i,ny-1,k)-ezy1(i,4,k))

		k=nz-1
		ez(i,1,k) = ezy1(i,2,k)+cyd*(ez(i,2,k)-ezy1(i,1,k))
		ez(i,ny,k)= ezy1(i,3,k)+cyu*(ez(i,ny-1,k)-ezy1(i,4,k))
	enddo

	!2次吸収境界条件
	do k=2,nz-2
		do i=3,nx-2
			ez(i,1,k)  = -ezy2(i,2,k)+cyd*(ez(i,2,k)+ezy2(i,1,k))&
						+cyy*(ezy1(i,1,k)+ezy1(i,2,k))&
						+cyfxd*(ezy1(i+1,1,k)-2.0d0*ezy1(i,1,k)+ezy1(i-1,1,k)&
							   +ezy1(i+1,2,k)-2.0d0*ezy1(i,2,k)+ezy1(i-1,2,k))&
						+cyfzd*(ezy1(i,1,k+1)-2.0d0*ezy1(i,1,k)+ezy1(i,1,k-1)&
							   +ezy1(i,2,k+1)-2.0d0*ezy1(i,2,k)+ezy1(i,2,k-1))
			ez(i,ny,k) = -ezy2(i,3,k)+cyd*(ez(i,ny-1,k)+ezy2(i,4,k))&
						+cyy*(ezy1(i,4,k)+ezy1(i,3,k))&
						+cyfxd*(ezy1(i+1,4,k)-2.0d0*ezy1(i,4,k)+ezy1(i-1,4,k)&
							   +ezy1(i+1,3,k)-2.0d0*ezy1(i,3,k)+ezy1(i-1,3,k))&
						+cyfzd*(ezy1(i,4,k+1)-2.0d0*ezy1(i,4,k)+ezy1(i,4,k-1)&
						       +ezy1(i,3,k+1)-2.0d0*ezy1(i,3,k)+ezy1(i,3,k-1))
		enddo
	enddo

	!過去の値の更新
	do k=1,nz-1
		do i=2,nx-1
			ezy2(i,1,k)=ezy1(i,1,k)
			ezy2(i,2,k)=ezy1(i,2,k)
			ezy2(i,3,k)=ezy1(i,3,k)
			ezy2(i,4,k)=ezy1(i,4,k)

			ezy1(i,1,k)=ez(i,1,k)
			ezy1(i,2,k)=ez(i,2,k)
			ezy1(i,3,k)=ez(i,ny-1,k)
			ezy1(i,4,k)=ez(i,ny,k)
		enddo
	enddo

	end subroutine mur_zx







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!(z軸方向に伝播) 486 未 引数をあわせていない  do範囲 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mur_xy(ex,ey,ez)
	use const_para
	implicit none

! 	real(8) :: czd,czu,czz
! 	real(8) :: czfxd,czfyd
! 	complex(kind(0d0)) :: eyx1(nx,ny,nz),eyx2(nx,ny,nz),ezx1(nx,ny,nz),ezx2(nx,ny,nz)
! 	complex(kind(0d0)) :: exy1(nx,ny,nz),exy2(nx,ny,nz),ezy1(nx,ny,nz),ezy2(nx,ny,nz)
! 	complex(kind(0d0)) :: exz1(nx,ny,nz),exz2(nx,ny,nz),eyz1(nx,ny,nz),eyz2(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)

	czd = (cmax*dt-dz)/(cmax*dt+dz)
	czu = (cmax*dt-dz)/(cmax*dt+dz)
	czz = (2.0d0*dz)/(cmax*dt+dz)
	czfxd = (dz*(cmax*dt)**2.0d0) / (2.0d0*(dx*dx)*(cmax*dt+dz))
	czfyd = (dz*(cmax+dt)**2.0d0) / (2.0d0*(dy*dy)*(cmax*dt+dz))


!-----------exに対して-----------------------------------------------------------------
	!1次吸収境界条件
	do j=2,ny-1
		i=1
		ex(i,j,1) = exz1(i,j,2)+czd*(ex(i,j,2)-exz1(i,j,1))
		ex(i,j,nz) =exz1(i,j,3)+czu*(ex(i,j,nz-1)-exz1(i,j,4))

		i=nx-1
		ex(i,j,1) = exz1(i,j,2)+czd*(ex(i,j,2)-exz1(i,j,1))
		ex(i,j,nz)= exz1(i,j,3)+czu*(ex(i,j,ny-1)-exz1(i,j,4))
	enddo

	do i=2,nx-2
		j=2
		ex(i,j,1) = exz1(i,j,2)+czd*(ex(i,j,2)-exz1(i,j,1))
		ex(i,j,nz)= exz1(i,j,3)+czu*(ex(i,j,nz-1)-exz1(i,j,4))

		j=ny-1
		ex(i,j,1) = exz1(i,j,2)+czd*(ex(i,j,2)-exz1(i,j,1))
		ex(i,j,nz)= exz1(i,j,3)+czu*(ex(i,j,nz-1)-exz1(i,j,4))
	enddo

	!2次吸収境界条件
	do j=3,ny-2
		do i=2,nx-2
			ex(i,j,1)  = -exz2(i,j,2)+czd*(ex(i,j,2)+exz2(i,j,1))&
						+czz*(exz1(i,j,1)+exz1(i,j,2))&
						+czfxd*(exz1(i+1,j,1)-2.0d0*exz1(i,j,1)+exz1(i-1,j,1)&
						       +exz1(i+1,j,2)-2.0d0*exz1(i,j,2)+exz1(i-1,j,2))&
						+czfyd*(exz1(i,j+1,1)-2.0d0*exz1(i,j,1)+exz1(i,j-1,1)&
						       +exz1(i,j+1,2)-2.0d0*exz1(i,j,2)+exz1(i,j-1,2))
			ex(i,j,nz) = -exz2(i,j,3)+czd*(ex(i,j,ny-1)+exz2(i,j,4))&
						+czz*(exz1(i,j,4)+exz1(i,j,3))&
						+czfxd*(exz1(i+1,j,4)-2.0d0*exz1(i,j,4)+exz1(i-1,j,4)&
						       +exz1(i+1,j,3)-2.0d0*exz1(i,j,3)+exz1(i-1,j,3))&
						+czfyd*(exz1(i,j+1,4)-2.0d0*exz1(i,j,4)+exz1(i,j-1,4)&
						       +exz1(i,j+1,3)-2.0d0*exz1(i,j,3)+exz1(i,j-1,3))
		enddo
	enddo

!過去の値の更
 do j=2,ny-1
	do i=1,nx-1
			exz2(i,j,1)=exz1(i,j,1)
			exz2(i,j,2)=exz1(i,j,2)
			exz2(i,j,3)=exz1(i,j,3)
			exz2(i,j,4)=exz1(i,j,4)

			exz1(i,j,1)=ex(j,j,1)
			exz1(i,j,2)=ex(i,j,2)
			exz1(i,j,3)=ex(i,j,nz-1)
			exz1(i,j,4)=ex(i,j,nz)
		enddo
	enddo

!---------eyに対して-----------------------------------------------------------------
	!1次吸収境界条件

	do j=1,ny-1
		i=2
		ey(i,j,1) = eyz1(i,j,2)+czd*(ey(i,j,2)-eyz1(i,j,1))
		ey(i,j,nz)= eyz1(i,j,3)+czu*(ey(i,j,nz-1)-eyz1(i,j,4))

		i=nx-1
		ey(i,j,1) = eyz1(i,j,2)+czd*(ey(i,j,2)-eyz1(i,j,1))
		ey(i,j,nz)= eyz1(i,j,3)+czu*(ey(i,j,nz-1)-eyz1(i,j,4))
	enddo

	do i=3,nx-2
		j=1
		ey(i,j,1) = eyz1(i,j,2)+czd*(ey(i,j,2)-eyz1(i,j,1))
		ey(i,j,nz)= eyz1(i,j,3)+czu*(ey(i,j,nz-1)-eyz1(i,j,4))

		j=ny-1
		ey(i,j,1) = eyz1(i,j,2)+czd*(ey(i,j,2)-eyz1(i,j,1))
		ey(i,j,nz)= eyz1(i,j,3)+czu*(ey(i,j,nz-1)-eyz1(i,j,4))
	enddo

	!2次吸収境界条件
	do j=2,ny-2
		do i=3,nx-2
			ey(i,j,1)  = -eyz2(i,j,2)+czd*(ey(i,j,2)+eyz2(i,j,1))&
						+czz*(eyz1(i,j,1)+eyz1(i,j,2))&
						+czfxd*(eyz1(i+1,j,1)-2.0d0*eyz1(i,j,1)+eyz1(i-1,j,1)&
						       +eyz1(i+1,j,2)-2.0d0*eyz1(i,j,2)+eyz1(i-1,j,2))&
						+czfyd*(eyz1(i,j+1,1)-2.0d0*eyz1(i,j,1)+eyz1(i,j-1,i)&
						       +eyz1(i,j+1,2)-2.0d0*eyz1(i,j,2)+eyz1(i,j-1,2))
			ey(i,j,nz) = -eyz2(i,j,3)+czd*(ey(i,j,ny-1)+eyz2(i,j,4))&
						+czz*(eyz1(i,j,4)+eyz1(i,j,3))&
						+czfxd*(eyz1(i+1,j,4)-2.0d0*eyz1(i,j,4)+eyz1(i-1,j,4)&
						       +eyz1(i+1,j,3)-2.0d0*eyz1(i,j,3)+eyz1(i-1,j,3))&
						+czfyd*(eyz1(i,j+1,4)-2.0d0*eyz1(i,j,4)+eyz1(i,j-1,4)&
						       +eyz1(i,j+1,3)-2.0d0*eyz1(i,j,3)+eyz1(i,j-1,3))
		enddo
	enddo

	!過去の値の更新
	do j=1,ny-1
		do i=2,nx-1
			eyz2(i,j,1)=eyz1(i,j,1)
			eyz2(i,j,2)=eyz1(i,j,2)
			eyz2(i,j,3)=eyz1(i,j,3)
			eyz2(i,j,4)=eyz1(i,j,4)

			eyz1(i,j,1)=ey(i,j,1)
			eyz1(i,j,2)=ey(j,j,2)
			eyz1(i,j,3)=ey(i,j,nz-1)
			eyz1(i,j,4)=ey(i,j,nz)
		enddo
	enddo

	   	end subroutine mur_xy