!(x軸方向に伝播) 486 ok  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mur_yz
	implicit none

	real(8) :: cxd = (cmax*dt-dx)/(cmax*dt+dx)
	real(8) :: cxu = (cmax*dt-dx)/(cmax*dt+dx)
	real(8) :: cxx = (2.0d0*dx)/(cmax*dt+dx)
	real(8) :: cxfyd = (dx*(cmax*dt)**2.0d0) / (2.0d0*(dy*dy)*(cmax*dt+dx))
	real(8) :: cxfzd = (dx*(cmax+dt)**2.0d0) / (2.0d0*(dz*dz)*(cmax*dt+dx))


!eyに対して-----------------------------------------------------------------
	!1次吸収境界条件
	do k=3,nz-2
		j=1
		ey(1,j,k) = eyx1(2,j,k)+cxd*(ey(2,j,k)-eyx1(1,j,k))
		ey(nx,j,k) =  eyx1(3,j,k)+cxu*(ey(nx-1,j,k)-eyx1(4,j,k))
		j=ny-1
		ey(1,j,k) = eyx1(2,j,k)+cxd*(ey(2,j,k)-eyx1(1,j,k))
		ey(nx,j,k)= eyx1(3,j,k)+cxu*(ey(nx-1,j,k)-eyx1(4,j,k))
	enddo

	do j=1,ny-1
		k=2
		ey(1,j,k) = eyx1(2,j,k)+cxd*(ey(2,j,k)-eyx1(1,j,k))
		ey(nx,j,k) = eyx1(3,j,k)+cxu*(ey(nx-1,j,k)-eyx1(4,j,k))
		k=nz-1
		ey(1,j,k) = eyx1(2,j,k)+cxd*(ey(2,j,k)-eyx1(1,j,k))
		ey(nx,j,k)= eyx1(3,j,k)+cxu*(ey(nx-1,j,k)-eyx1(4,j,k))
	enddo
	!2次吸収境界条件
	do k=3,nz-2
		do j=2,ny-2
			ey(1,j,k) = -eyx2(2,j,k)+cxd*(ey(2,j,k)+eyx2(1,j,k))&
						+cxx*(eyx1(1,j,k)+eyx1(2,j,k))&
						+cxfyd*(eyx1(1,j,k+1)-2.0d0*eyx1(1,j,k)+eyx1(1,j,k-1)&
						+eyx1(2,j,k+1)-2.0d0*eyx1(2,j,k)+eyx1(2,j,k-1))&
						+cxfzd*(eyx1(1,j+1,k)-2.0d0*eyx1(1,j,k)+eyx1(1,j-1,k)&
						+eyx1(2,j+1,k)-2.0d0*eyx1(2,j,k)+eyx1(2,j-1,k))
			ey(nx,j,k) = -eyx2(3,j,k)+cxd*(ey(nx-1,j,k)+eyx2(4,j,k))&
						+cxx*(eyx1(4,j,k)+eyx1(3,j,k))&
						+cxfyd*(eyx1(4,j,k+1)-2.0d0*eyx1(4,j,k)+eyx1(4,j,k-1)&
						+eyx1(3,j,k+1)-2.0d0*eyx1(3,j,k)+eyx1(3,j,k-1))&
						+cxfzd*(eyx1(4,j+1,k)-2.0d0*eyx1(4,j,k)+eyx1(4,j-1,k)&
						+eyx1(3,j+1,k)-2.0d0*eyx1(3,j,k)+eyx1(3,j-1,k))
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

!ezに対して-----------------------------------------------------------------
	!1次吸収境界条件
	do k=1,nz-1
		j=2
		ez(1,j,k) = ezx1(2,j,k)+cxd*(ez(2,j,k)-ezx1(1,j,k))
		ez(nx,j,k) =  ezx1(3,j,k)+cxu*(ez(nx-1,j,k)-ezx1(4,j,k))
		j=ny-1
		ez(1,j,k) = ezx1(2,j,k)+cxd*(ez(2,j,k)-ezx1(1,j,k))
		ez(nx,j,k)= ezx1(3,j,k)+cxu*(ez(nx-1,j,k)-ezx1(4,j,k))
	enddo

	do j=3,ny-2
		k=1
		ez(1,j,k) = ezx1(2,j,k)+cxd*(ez(2,j,k)-ezx1(1,j,k))
		ez(nx,j,k) = ezx1(3,j,k)+cxu*(ez(nx-1,j,k)-ezx1(4,j,k))
		k=nz-1
		ez(1,j,k) = ezx1(2,j,k)+cxd*(ez(2,j,k)-ezx1(1,j,k))
		ez(nx,j,k)= ezx1(3,j,k)+cxu*(ez(nx-1,j,k)-ezx1(4,j,k))
	enddo

	!2次吸収境界条件
	do k=2,nz-2
		do j=3,ny-2
			ez(1,j,k) = -ezx2(2,j,k)+cxd*(ez(2,j,k)+ezx2(1,j,k))&
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







!(y軸方向に伝播) 未  引数を合わせていない do範囲!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine mur_zx
	implicit none

	real(8) :: cyd = (cmax*dt-dy)/(cmax*dt+dy)
	real(8) :: cyu = (cmax*dt-dy)/(cmax*dt+dy)
	real(8) :: cyy = (2.0d0*dy)/(cmax*dt+dy)
	real(8) :: cyfxd = (dy*(cmax*dt)**2.0d0) / (2.0d0*(dx*dx)*(cmax*dt+dy))
	real(8) :: cyfzd = (dy*(cmax+dt)**2.0d0) / (2.0d0*(dz*dz)*(cmax*dt+dy))


!exに対して-----------------------------------------------------------------
	!1次吸収境界条件
	do k=3,nz-2
		j=1
		ex(1,j,k) = exy1(2,j,k)+cyd*(ex(2,j,k)-exy1(1,j,k))
		ex(nx,j,k) =  exy1(3,j,k)+cyu*(ex(nx-1,j,k)-exy1(4,j,k))
		j=ny-1
		ex(1,j,k) = exy1(2,j,k)+cyd*(ex(2,j,k)-exy1(1,j,k))
		ex(nx,j,k)= exy1(3,j,k)+cyu*(ex(nx-1,j,k)-exy1(4,j,k))
	enddo

	do i=1,ny-1
		k=2
		ex(1,j,k) = exy1(2,j,k)+cyd*(ex(2,j,k)-exy1(1,j,k))
		ex(nx,j,k) = exy1(3,j,k)+cyu*(ex(nx-1,j,k)-exy1(4,j,k))
		k=nz-1
		ex(1,j,k) = exy1(2,j,k)+cyd*(ex(2,j,k)-exy1(1,j,k))
		ex(nx,j,k)= exy1(3,j,k)+cyu*(ex(nx-1,j,k)-exy1(4,j,k))
	enddo
	!2次吸収境界条件
	do k=3,nz-2
		do j=2,ny-2
			ex(1,j,k) = -exy2(2,j,k)+cyd*(ex(2,j,k)+exy2(1,j,k))&
						+cyy*(exy1(1,j,k)+exy1(2,j,k))&
						+cyfxd*(exy1(1,j,k+1)-2.0d0*exy1(1,j,k)+exy1(1,j,k-1)&
						+exy1(2,j,k+1)-2.0d0*exy1(2,j,k)+exy1(2,j,k-1))&
						+cyfzd*(exy1(1,j+1,k)-2.0d0*exy1(1,j,k)+exy1(1,j-1,k)&
						+exy1(2,j+1,k)-2.0d0*exy1(2,j,k)+exy1(2,j-1,k))
			ex(nx,j,k) = -exy2(3,j,k)+cyd*(ex(nx-1,j,k)+exy2(4,j,k))&
						+cyy*(exy1(4,j,k)+exy1(3,j,k))&
						+cyfxd*(exy1(4,j,k+1)-2.0d0*exy1(4,j,k)+exy1(4,j,k-1)&
						+exy1(3,j,k+1)-2.0d0*exy1(3,j,k)+exy1(3,j,k-1))&
						+cyfzd*(exy1(4,j+1,k)-2.0d0*exy1(4,j,k)+exy1(4,j-1,k)&
						+exy1(3,j+1,k)-2.0d0*exy1(3,j,k)+exy1(3,j-1,k))
		enddo
	enddo
	!過去の値の更新
	do k=2,nz-1
		do j=1,ny-1
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

!ezに対して-----------------------------------------------------------------
	!1次吸収境界条件
	do k=1,nz-1
		j=2
		ez(1,j,k) = ezy1(2,j,k)+cyd*(ez(2,j,k)-ezy1(1,j,k))
		ez(nx,j,k) =  ezy1(3,j,k)+cyu*(ez(nx-1,j,k)-ezy1(4,j,k))
		j=ny-1
		ez(1,j,k) = ezy1(2,j,k)+cyd*(ez(2,j,k)-ezy1(1,j,k))
		ez(nx,j,k)= ezy1(3,j,k)+cyu*(ez(nx-1,j,k)-ezy1(4,j,k))
	enddo

	do j=3,ny-2
		k=1
		ez(1,j,k) = ezy1(2,j,k)+cyd*(ez(2,j,k)-ezy1(1,j,k))
		ez(nx,j,k) = ezy1(3,j,k)+cyu*(ez(nx-1,j,k)-ezy1(4,j,k))
		k=nz-1
		ez(1,j,k) = ezy1(2,j,k)+cyd*(ez(2,j,k)-ezy1(1,j,k))
		ez(nx,j,k)= ezy1(3,j,k)+cyu*(ez(nx-1,j,k)-ezy1(4,j,k))
	enddo

	!2次吸収境界条件
	do k=2,nz-2
		do j=3,ny-2
			ez(1,j,k) = -ezy2(2,j,k)+cyd*(ez(2,j,k)+ezy2(1,j,k))&
						+cyy*(ezy1(1,j,k)+ezy1(2,j,k))&
						+cyfxd*(ezy1(1,j+1,k)-2.0d0*ezy1(1,j,k)+ezy1(1,j-1,k)&
						+ezy1(2,j+1,k)-2.0d0*ezy1(2,j,k)+ezy1(2,j-1,k))&
						+cyfzd*(ezy1(1,j,k+1)-2.0d0*ezy1(1,j,k)+ezy1(1,j,k-1)&
						+ezy1(2,j,k+1)-2.0d0*ezy1(2,j,k)+ezy1(2,j,k-1))
			ez(nx,j,k) = -ezy2(3,j,k)+cyd*(ez(nx-1,j,k)+ezy2(4,j,k))&
						+cyy*(ezy1(4,j,k)+ezy1(3,j,k))&
						+cyfxd*(ezy1(4,j+1,k)-2.0d0*ezy1(4,j,k)+ezy1(4,j-1,k)&
						+ezy1(3,j+1,k)-2.0d0*ezy1(3,j,k)+ezy1(3,j-1,k))&
						+cyfzd*(ezy1(4,j,k+1)-2.0d0*ezy1(4,j,k)+ezy1(4,j,k-1)&
						+ezy1(3,j,k+1)-2.0d0*ezy1(3,j,k)+ezy1(3,j,k-1))
		enddo
	enddo
	!過去の値の更新
	do k=1,nz-1
		do j=2,ny-1
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








!(z軸方向に伝播) 486 未 引数をあわせていない  do範囲 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mur_xy
	implicit none

	real(8) :: czd = (cmax*dt-dz)/(cmax*dt+dz)
	real(8) :: czu = (cmax*dt-dz)/(cmax*dt+dz)
	real(8) :: czz = (2.0d0*dz)/(cmax*dt+dz)
	real(8) :: czfxd = (dz*(cmax*dt)**2.0d0) / (2.0d0*(dx*dx)*(cmax*dt+dz))
	real(8) :: czfyd = (dz*(cmax+dt)**2.0d0) / (2.0d0*(dy*dy)*(cmax*dt+dz))


!exに対して-----------------------------------------------------------------
	!1次吸収境界条件
	do k=3,nz-2
		j=1
		ex(1,j,k) = exz1(2,j,k)+czd*(ex(2,j,k)-exz1(1,j,k))
		ex(nx,j,k) =  exz1(3,j,k)+czu*(ex(nx-1,j,k)-exz1(4,j,k))
		j=ny-1
		ex(1,j,k) = exz1(2,j,k)+czd*(ex(2,j,k)-exz1(1,j,k))
		ex(nx,j,k)= exz1(3,j,k)+czu*(ex(nx-1,j,k)-exz1(4,j,k))
	enddo

	do j=1,ny-1
		k=2
		ex(1,j,k) = exz1(2,j,k)+czd*(ex(2,j,k)-exz1(1,j,k))
		ex(nx,j,k) = exz1(3,j,k)+czu*(ex(nx-1,j,k)-exz1(4,j,k))
		k=nz-1
		ex(1,j,k) = exz1(2,j,k)+czd*(ex(2,j,k)-exz1(1,j,k))
		ex(nx,j,k)= exz1(3,j,k)+czu*(ex(nx-1,j,k)-exz1(4,j,k))
	enddo
	!2次吸収境界条件
	do k=3,nz-2
		do j=2,ny-2
			ex(1,j,k) = -exz2(2,j,k)+czd*(ex(2,j,k)+exz2(1,j,k))&
						+czz*(exz1(1,j,k)+exz1(2,j,k))&
						+czfxd*(exz1(1,j,k+1)-2.0d0*exz1(1,j,k)+exz1(1,j,k-1)&
						+exz1(2,j,k+1)-2.0d0*exz1(2,j,k)+exz1(2,j,k-1))&
						+czfyd*(exz1(1,j+1,k)-2.0d0*exz1(1,j,k)+exz1(1,j-1,k)&
						+exz1(2,j+1,k)-2.0d0*exz1(2,j,k)+exz1(2,j-1,k))
			ex(nx,j,k) = -exz2(3,j,k)+czd*(ex(nx-1,j,k)+exz2(4,j,k))&
						+czz*(exz1(4,j,k)+exz1(3,j,k))&
						+czfxd*(exz1(4,j,k+1)-2.0d0*exz1(4,j,k)+exz1(4,j,k-1)&
						+exz1(3,j,k+1)-2.0d0*exz1(3,j,k)+exz1(3,j,k-1))&
						+czfyd*(exz1(4,j+1,k)-2.0d0*exz1(4,j,k)+exz1(4,j-1,k)&
						+exz1(3,j+1,k)-2.0d0*exz1(3,j,k)+exz1(3,j-1,k))
		enddo
	enddo
!過去の値の更i
do k=2,nz-i
	do j=1,ny-i
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

!eyに対して-----------------------------------------------------------------
	!1次吸収境界条件
	do k=1,nz-1
		j=2
		ey(1,j,k) = eyz1(2,j,k)+czd*(ey(2,j,k)-eyz1(1,j,k))
		ey(nx,j,k) =  eyz1(3,j,k)+czu*(ey(nx-1,j,k)-eyz1(4,j,k))
		j=ny-1
		ey(1,j,k) = eyz1(2,j,k)+czd*(ey(2,j,k)-eyz1(1,j,k))
		ey(nx,j,k)= eyz1(3,j,k)+czu*(ey(nx-1,j,k)-eyz1(4,j,k))
	enddo

	do j=3,ny-2
		k=1
		ey(1,j,k) = eyz1(2,j,k)+czd*(ey(2,j,k)-eyz1(1,j,k))
		ey(nx,j,k) = eyz1(3,j,k)+czu*(ey(nx-1,j,k)-eyz1(4,j,k))
		k=nz-1
		ey(1,j,k) = eyz1(2,j,k)+czd*(ey(2,j,k)-eyz1(1,j,k))
		ey(nx,j,k)= eyz1(3,j,k)+czu*(ey(nx-1,j,k)-eyz1(4,j,k))
	enddo

	!2次吸収境界条件
	do k=2,nz-2
		do j=3,ny-2
			ey(1,j,k) = -eyz2(2,j,k)+czd*(ey(2,j,k)+eyz2(1,j,k))&
						+czz*(eyz1(1,j,k)+eyz1(2,j,k))&
						+czfxd*(eyz1(1,j+1,k)-2.0d0*eyz1(1,j,k)+eyz1(1,j-1,k)&
						+eyz1(2,j+1,k)-2.0d0*eyz1(2,j,k)+eyz1(2,j-1,k))&
						+czfyd*(eyz1(1,j,k+1)-2.0d0*eyz1(1,j,k)+eyz1(1,j,k-1)&
						+eyz1(2,j,k+1)-2.0d0*eyz1(2,j,k)+eyz1(2,j,k-1))
			ey(nx,j,k) = -eyz2(3,j,k)+czd*(ey(nx-1,j,k)+eyz2(4,j,k))&
						+czz*(eyz1(4,j,k)+eyz1(3,j,k))&
						+czfxd*(eyz1(4,j+1,k)-2.0d0*eyz1(4,j,k)+eyz1(4,j-1,k)&
						+eyz1(3,j+1,k)-2.0d0*eyz1(3,j,k)+eyz1(3,j-1,k))&
						+czfyd*(eyz1(4,j,k+1)-2.0d0*eyz1(4,j,k)+eyz1(4,j,k-1)&
						+eyz1(3,j,k+1)-2.0d0*eyz1(3,j,k)+eyz1(3,j,k-1))
		enddo
	enddo
	!過去の値の更新
	do k=1,nz-1
		do j=2,ny-1
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