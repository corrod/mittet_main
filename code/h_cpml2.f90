!Convolutional PML_H !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sigmamax amax kappamax の求め方 sigma*の求め方
! 導入の仕方
! メイン部分の計算は通常の電磁波伝播と同じでプサイのぶぶんだけCPML？
! subrouitne 分ける必要ないのかも
!kappa >=1 real
!sigmai>0 real
!ai=alphai>0 real
!psi loop あと３枚必要？
!割り算のとき倍精度d0に注意!!
!db_x,db_y,db_zを仮想領域にあわせる
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CPML_H(ex,ey,ez,hx,hy,hz,sigma,myu)
    use const_para
    implicit none

    integer, parameter :: m = 4, ma = 4
    integer, parameter :: nxpml1 = 10, nypml1 = 10, nzpml1 = 10 !pmlの厚さ
    real(8), parameter :: kappa_max = 7.0d0 !!!
    real(8), parameter :: a_max = 0.2d0     !!!
    real(8), intent(in) :: myu(nx,ny,nz)
    real(8), intent(in) :: sigma(nx,ny,nz)
    real(8) :: sigma_opt
    real(8) :: sigma_max!!!
    real(8) :: msigma_x(nx),msigma_y(ny),msigma_z(nz)
!     real(8), parameter :: lnR0 = -100.0d0  !ln|R(0)|
    real(8) :: epsi(nx,ny,nz)!1.0d0
    real(8) :: epsir = 1.0d0
    real(8) :: am_x(nx),am_y(ny),am_z(nz)
    real(8) :: mkappa_x(nx),mkappa_y(ny),mkappa_z(nz)
    real(8) :: khdx(nx),khdy(ny),khdz(nz)
    real(8) :: da_x(nx,ny,nz),da_y(nx,ny,nz),da_z(nx,ny,nz)
    real(8) :: db_x(nx,ny,nz),db_y(nx,ny,nz),db_z(nx,ny,nz)
    real(8) :: bh_x(nxpml1),bh_y(nypml1),bh_z(nzpml1)
    real(8) :: ch_x(nxpml1),ch_y(nypml1),ch_z(nzpml1)
    complex(kind(0d0)) :: psi_Hzx1(nx,ny,nz),psi_Hyx1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Hxy1(nx,ny,nz),psi_Hzy1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Hyz1(nx,ny,nz),psi_Hxz1(nx,ny,nz)
    complex(kind(0d0)), intent(in)    :: ex(nx,ny,nz),ey(nx,ny,nz),ez(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: hx(nx,ny,nz),hy(nx,ny,nz),hz(nx,ny,nz)

    epsi(1:nx,1:ny,1:nz)=sigma(1:nx,1:ny,1:nz)/(2.0d0*omega0)

!!!    sigma_max = -(m+1)*lnR0 / (2.0d0*(sqrt(myu/epsi))*nxpml1*dx)  !ln(R(0));反射係数!!!
    sigma_opt = (dble(m)+1.0d0) / (150.0d0*pai*sqrt(epsir)*dx)
    sigma_max = 0.7d0*sigma_opt
    

!係数の設定

    do i = 1,nx
if (i<=nxpml1) then
    msigma_x(i) = sigma_max* ((dble(nxpml1)-dble(i)-0.5d0)/(dble(nxpml1)-1.0d0))**dble(m)  !!!-i-1/2の取り扱い
    mkappa_x(i) = 1.0d0 + (kappa_max-1.0d0)*((dble(nxpml1)-dble(i)-0.5d0)/(dble(nxpml1)-1.0d0))**dble(m)  !!!-i-1/2の取扱い
    am_x(i)     = a_max* ((dble(i)-0.5d0)/(dble(nxpml1)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

    bh_x(i)     = exp(-(msigma_x(i)/mkappa_x(i)+am_x(i)) *dt) !/epsi0)
    ch_x(i)     = msigma_x(i)*(bh_x(i)-1.0d0) / (msigma_x(i) + mkappa_x(i)*am_x(i)) / mkappa_x(i)
    khdx(i)     = mkappa_x(i)*dx !!!(i-1/2)dxの取り扱い


else if(i>=nx-nxpml1+1) then
    msigma_x(i) = sigma_max* ((dble(i)-dble(nx)+0.5d0+dble(nxpml1))/(dble(nxpml1)-1.0d0))**dble(m)  !!!-i-1/2の取り扱い
    mkappa_x(i) = 1.0d0 + (kappa_max-1.0d0)*((dble(i)-dble(nx)+0.5d0+dble(nxpml1))/(dble(nxpml1)-1.0d0))**dble(m)  !!!-i-1/2の取扱い
    am_x(i)     = a_max* ((dble(-i)+nx-0.5d0)/(dble(nxpml1)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

    bh_x(i)     = exp(-(msigma_x(i)/mkappa_x(i)+am_x(i)) *dt) !/epsi0)
    ch_x(i)     = msigma_x(i)*(bh_x(i)-1.0d0) / (msigma_x(i) + mkappa_x(i)*am_x(i)) / mkappa_x(i) 
    khdx(i)     = mkappa_x(i)*dx !!!(i-1/2)dxの取り扱い
else
    msigma_x(i) = 0.0d0
    mkappa_x(i) = 1.0d0
    am_x(i) = 0.0d0   
    bh_x(i) = 0.0d0   
    ch_x(i) = 0.0d0   
    khdx(i) = mkappa_x(i)*dx
    endif
        enddo




do j = 1,ny
if (j<=nypml1) then
    msigma_y(j) = sigma_max* ((dble(nypml1)-dble(j)-0.5d0)/(dble(nypml1)-1.0d0))**dble(m)  !!!-i-1/2の取り扱い
    mkappa_y(j) = 1.0d0 + (kappa_max-1.0d0)*((dble(nypml1)-dble(j)-0.5d0)/(dble(nypml1)-1.0d0))**dble(m)  !!!-i-1/2の取扱い
    am_y(j)     = a_max* ((dble(j)-0.5d0)/(dble(nypml1)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

    bh_y(j)     = exp(-(msigma_y(j)/mkappa_y(j)+am_y(j)) *dt) !/epsi0)
    ch_y(j)     = msigma_y(j)*(bh_y(j)-1.0d0) / (msigma_y(j) + mkappa_y(i)*am_y(j)) / mkappa_y(j)
    khdy(j)     = mkappa_y(j)*dx !!!(i-1/2)dxの取り扱い

else if(j>=ny-nypml1+1) then
    msigma_y(j) = sigma_max* ((dble(j)-dble(ny)+0.5d0+dble(nypml1))/(dble(nypml1)-1.0d0))**dble(m)  !!!-i-1/2の取り扱い
    mkappa_y(j) = 1.0d0 + (kappa_max-1.0d0)*((dble(j)-dble(ny)+0.5d0+dble(nypml1))/(dble(nypml1)-1.0d0))**dble(m)  !!!-i-1/2の取扱い
    am_y(j)     = a_max* ((dble(-j)+ny-0.5d0)/(dble(nypml1)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

    bh_y(j)     = exp(-(msigma_y(j)/mkappa_y(j)+am_y(j)) *dt) !/epsi0)
    ch_y(j)     = msigma_y(j)*(bh_y(j)-1.0d0) / (msigma_y(j) + mkappa_y(j)*am_y(j)) / mkappa_y(j) 
    khdy(j)     = mkappa_y(j)*dy !!!(i-1/2)dxの取り扱い
else
    msigma_y(j) = 0.0d0
    mkappa_y(j) = 1.0d0
    am_y(j) = 0.0d0   
    bh_y(j) = 0.0d0   
    ch_y(j) = 0.0d0   
    khdy(j) = mkappa_y(j)*dy
    endif
        enddo



do k = 1,nz
if (k<=nzpml1) then
    msigma_z(k) = sigma_max* ((dble(nzpml1)-dble(k)-0.5d0)/(dble(nzpml1)-1.0d0))**dble(m)  !!!-i-1/2の取り扱い
    mkappa_z(k) = 1.0d0 + (kappa_max-1.0d0)*((dble(nzpml1)-dble(k)-0.5d0)/(dble(nzpml1)-1.0d0))**dble(m)  !!!-i-1/2の取扱い
    am_z(k)     = a_max* ((dble(k)-0.5d0)/(dble(nzpml1)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

    bh_z(k)     = exp(-(msigma_z(k)/mkappa_z(k)+am_z(k)) *dt) !/epsi0)
    ch_z(k)     = msigma_z(k)*(bh_z(k)-1.0d0) / (msigma_z(k) + mkappa_z(k)*am_z(k)) / mkappa_z(k)
    khdz(k)     = mkappa_z(k)*dz !!!(i-1/2)dxの取り扱い

else if(k>=nz-nzpml1+1) then
    msigma_z(k) = sigma_max* ((dble(k)-dble(nz)+0.5d0+dble(nzpml1))/(dble(nzpml1)-1.0d0))**dble(m)  !!!-i-1/2の取り扱い
    mkappa_z(k) = 1.0d0 + (kappa_max-1.0d0)*((dble(k)-dble(nz)+0.5d0+dble(nzpml1))/(dble(nzpml1)-1.0d0))**dble(m)  !!!-i-1/2の取扱い
    am_z(k)     = a_max* ((dble(-k)+nz-0.5d0)/(dble(nzpml1)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

    bh_z(k)     = exp(-(msigma_z(k)/mkappa_z(k)+am_z(k)) *dt) !/epsi0)
    ch_z(k)     = msigma_z(k)*(bh_z(k)-1.0d0) / (msigma_z(k) + mkappa_z(k)*am_z(k)) / mkappa_z(k) 
    khdz(k)     = mkappa_z(k)*dz !!!(i-1/2)dxの取り扱い
else
    msigma_z(k) = 0.0d0
    mkappa_z(k) = 1.0d0
    am_z(k) = 0.0d0   
    bh_z(k) = 0.0d0   
    ch_z(k) = 0.0d0   
    khdz(k) = mkappa_z(k)*dz
    endif
        enddo
























!     (+)の係数
!     do i=1,nxpml1
!         msigma_x(nx-(i-1))= msigma_x(i)
!         mkappa_x(nx-(i-1))=mkappa_x(i)
!         am_x(nx-(i-1)) = am_x(i)
!         khdx(nx-(i-1)) = khdx(i)
!         bh_x(nx-(i-1)) = bh_x(i)
!         ch_x(nx-(i-1)) = ch_x(i)
!     enddo

!     do j=1,nypml1
!         msigma_y(ny-(j-1))= msigma_y(j)
!         mkappa_y(ny-(j-1))=mkappa_y(j)
!         am_y(ny-(j-1)) = am_y(j)
!         khdy(ny-(j-1)) = khdy(j)
!         bh_y(ny-(j-1)) = bh_y(j)
!         ch_y(ny-(j-1)) = ch_y(j)
!     enddo

!     do k=1,nzpml1
!         msigma_z(nz-(k-1))= msigma_z(k)
!         mkappa_z(nz-(k-1))=mkappa_z(k)
!         am_z(nz-(k-1)) = am_z(k)
!         khdz(nz-(k-1)) = khdz(k)
!         bh_z(nz-(k-1)) = bh_z(k)
!         ch_z(nz-(k-1)) = ch_z(k)
!     enddo








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!psi update!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!x-PML loop(-)
    do k = 1,nz-1
        do j = 1,ny-1
            do i = 1,nxpml1-1
                psi_Hzx1(i,j,k) = bh_x(i)*psi_Hzx1(i,j,k)&
                                 +ch_x(i)*(ey(i+1,j,k)-ey(i,j,k)) / dx
                psi_Hyx1(i,j,k) = bh_x(i)*psi_Hyx1(i,j,k)&
                                 +ch_x(i)*(ez(i+1,j,k)-ez(i,j,k)) / dx
                hz(i,j,k) = hz(i,j,k) - db_z(i,j,k)*psi_Hzx1(i,j,k)
                hy(i,j,k) = hy(i,j,k) + db_y(i,j,k)*psi_Hyx1(i,j,k)
                   enddo
                        enddo
                            enddo
!x-PML loop(+)
    do k = 1,nz-1
        do j = 1,ny-1
            do i = nx-nxpml1+1,nx-1
                psi_Hzx1(i,j,k) = bh_x(i)*psi_Hzx1(i,j,k)&
                                 +ch_x(i)*(ey(i+1,j,k)-ey(i,j,k)) / dx
                psi_Hyx1(i,j,k) = bh_x(i)*psi_Hyx1(i,j,k)&
                                 +ch_x(i)*(ez(i+1,j,k)-ez(i,j,k)) / dx
                hz(i,j,k) = hz(i,j,k) - db_z(i,j,k)*psi_Hzx1(i,j,k)
                hy(i,j,k) = hy(i,j,k) + db_y(i,j,k)*psi_Hyx1(i,j,k)
                   enddo
                        enddo
                            enddo





!y-PML loop(-)
    do k = 1,nz-1
        do j = 1,nypml1-1
            do i = 1,nx-1
                psi_Hxy1(i,j,k) = bh_y(j)*psi_Hxy1(i,j,k)&
                                 +ch_y(j)*(ez(i,j+1,k)-ez(i,j,k)) / dy
                psi_Hzy1(i,j,k) = bh_y(j)*psi_Hzy1(i,j,k)&
                                 +ch_y(j)*(ex(i,j+1,k)-ex(i,j,k)) / dy
                hx(i,j,k) = hx(i,j,k) - db_x(i,j,k)*psi_Hxy1(i,j,k)
                hz(i,j,k) = hz(i,j,k) + db_z(i,j,k)*psi_Hzy1(i,j,k)
                  enddo
                       enddo
                           enddo
!!y-PML loop(+)
    do k = 1,nz-1
        do j = ny-nypml1+1,ny-1
            do i = 1,nx-1
                psi_Hxy1(i,j,k) = bh_y(j)*psi_Hxy1(i,j,k)&
                                 +ch_y(j)*(ez(i,j+1,k)-ez(i,j,k)) / dy
                psi_Hzy1(i,j,k) = bh_y(j)*psi_Hzy1(i,j,k)&
                                 +ch_y(j)*(ex(i,j+1,k)-ex(i,j,k)) / dy
                hx(i,j,k) = hx(i,j,k) - db_x(i,j,k)*psi_Hxy1(i,j,k)
                hz(i,j,k) = hz(i,j,k) + db_z(i,j,k)*psi_Hzy1(i,j,k)
                  enddo
                       enddo
                           enddo




!z-PML loop(-)
    do k = 1,nzpml1-1
        do j = 1,ny-1
            do i = 1,nx-1
                psi_Hyz1(i,j,k) = bh_z(k)*psi_Hyz1(i,j,k)&
                                 +ch_z(k)*(ex(i,j,k+1)-ex(i,j,k)) / dz
                psi_Hxz1(i,j,k) = bh_y(k)*psi_Hxz1(i,j,k)&
                                 +ch_z(k)*(ey(i,j,k+1)-ey(i,j,k)) / dz
                hy(i,j,k) = hy(i,j,k) - db_y(i,j,k)*psi_Hyz1(i,j,k)
                hx(i,j,k) = hx(i,j,k) + db_x(i,j,k)*psi_Hxz1(i,j,k)
                   enddo
                        enddo
                            enddo
!z-PML loop(+)
    do k = nz-nzpml1+1,nz-1
        do j = 1,ny-1
            do i = 1,nx-1
                psi_Hyz1(i,j,k) = bh_z(k)*psi_Hyz1(i,j,k)&
                                 +ch_z(k)*(ex(i,j,k+1)-ex(i,j,k)) / dz
                psi_Hxz1(i,j,k) = bh_y(k)*psi_Hxz1(i,j,k)&
                                 +ch_z(k)*(ey(i,j,k+1)-ey(i,j,k)) / dz
                hy(i,j,k) = hy(i,j,k) - db_y(i,j,k)*psi_Hyz1(i,j,k)
                hx(i,j,k) = hx(i,j,k) + db_x(i,j,k)*psi_Hxz1(i,j,k)
                   enddo
                        enddo
                            enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!field update loop!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do k=1,nz
!  do j=1,ny
!    do i=1,nx
!    da_x(i,j,k) = (1.0d0-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) !sigma=σ*
!    da_y(i,j,k) = (1.0d0-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))!導磁率σ
!    da_z(i,j,k) = (1.0d0-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))

!    db_x(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
!    db_y(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
!    db_z(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
!enddo
!enddo
!enddo
!  !Hx
!  do k=1,nz-1
!        do j=1,ny-1
!            do i=2,nx-1
!                hx(i,j,k) = da_x(i,j,k)*hx(i,j,k)&
!                           -db_x(i,j,k)*((ez(i,j+1,k)-ez(i,j,k))/khdy(j)&
!                                        -(ey(i,j,k+1)-ey(i,j,k))/khdz(k))
!                           enddo
!                                enddo
!                                    enddo
!  !Hy
!  do k=1,nz-1
!        do j=2,ny-1
!            do i=1,nx-1
!                hy(i,j,k) = da_y(i,j,k)*hy(i,j,k)&
!                           -db_y(i,j,k)*((ex(i,j,k+1)-ex(i,j,k))/khdz(k)&
!                                        -(ez(i+1,j,k)-ez(i,j,k))/khdx(i))
!                            enddo
!                                enddo
!                                    enddo
!  !Hz
! do k=2,nz-1
!        do j=1,ny-1
!            do i=1,nx-1
!                hz(i,j,k) = da_z(i,j,k)*hz(i,j,k)&
!                           -db_z(i,j,k)*((ey(i+1,j,k)-ey(i,j,k))/khdx(i)&
!                                        -(ex(i,j+1,k)-ex(i,j,k))/khdy(j))
!                            enddo
!                                enddo
!                                    enddo
            end subroutine CPML_H
