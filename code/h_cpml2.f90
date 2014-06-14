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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CPML_H(ex,ey,ez,hx,hy,hz,sigma,myu)
    use const_para
    implicit none

    integer, parameter :: m = 4, ma = 4
    integer, parameter :: nxpml1 = 10, nypml1 = 10, nzpml1 = 10 !pmlの厚さ
    real(8), parameter :: kappa_max = 7.0d0 !!!
    real(8), parameter :: a_max = 0.05d0     !!!
    real(8), intent(in) :: myu(nx,ny,nz)
    real(8), intent(in) :: sigma(nx,ny,nz)
    real(8) :: sigma_opt
    real(8) :: sigma_max!!!
    real(8) :: sigma_x(nxpml1),sigma_y(nypml1),sigma_z(nzpml1)
!     real(8), parameter :: lnR0 = -100.0d0  !ln|R(0)|
    real(8) :: epsi(nx,ny,nz)!1.0d0
    real(8) :: epsir = 1.0d0
    real(8) :: a_x(nxpml1),a_y(nypml1),a_z(nzpml1)
    real(8) :: kappa_x(nxpml1),kappa_y(nypml1),kappa_z(nzpml1)
    real(8) :: khdx(nxpml1),khdy(nypml1),khdz(nzpml1)
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
    sigma_opt = (m+1.0d0) / (150.0d0*pai*sqrt(epsir)*dx)
    sigma_max = 0.7d0*sigma_opt
    

!係数の設定
do k = 1,nzpml1
  do j = 1,nypml1
    do i = 1,nxpml1

    sigma_x(i) = sigma_max* ((nxpml1-i)/(nxpml1-1.0d0))**m  !!!-i-1/2の取り扱い
    sigma_y(j) = sigma_max* ((nypml1-j)/(nypml1-1.0d0))**m
    sigma_z(k) = sigma_max* ((nzpml1-k)/(nzpml1-1.0d0))**m    
!!!    kappa_max =    !!!導出要確認
    kappa_x(i) = 1.0d0 + (kappa_max-1.0d0)*((nxpml1-i)/(nxpml1-1.0d0))**m  !!!-i-1/2の取扱い
    kappa_y(j) = 1.0d0 + (kappa_max-1.0d0)*((nypml1-j)/(nypml1-1.0d0))**m
    kappa_z(k) = 1.0d0 + (kappa_max-1.0d0)*((nzpml1-k)/(nzpml1-1.0d0))**m

!!!    a_max =         !!導出要確認
    a_x(i) = a_max* (i/(nxpml1-1.0d0))**ma !!!-i-1/2の取り扱い
    a_y(j) = a_max* (j/(nypml1-1.0d0))**ma
    a_z(k) = a_max* (k/(nzpml1-1.0d0))**ma

    khdx(i) = kappa_x(i)*dx !!!(i-1/2)dxの取り扱い
    khdy(j) = kappa_y(j)*dy
    khdz(k) = kappa_z(k)*dz


    bh_x(i) = exp(-(sigma_x(i)/kappa_x(i)+a_x(i)) *dt/epsi0)
    bh_y(j) = exp(-(sigma_y(j)/kappa_y(j)+a_y(j)) *dt/epsi0)
    bh_z(k) = exp(-(sigma_z(k)/kappa_z(k)+a_z(k)) *dt/epsi0)

    ch_x(i) = sigma_x(j)*(bh_x(i)-1.0d0) / (sigma_x(i) + kappa_x(i)*a_x(i)) / kappa_x(i)
    ch_y(j) = sigma_y(j)*(bh_y(j)-1.0d0) / (sigma_y(j) + kappa_y(j)*a_y(j)) / kappa_y(j)
    ch_z(j) = sigma_z(k)*(bh_z(k)-1.0d0) / (sigma_z(k) + kappa_z(k)*a_z(k)) / kappa_z(k)

    enddo
      enddo
        enddo

do k=1,nz
  do j=1,ny
    do i=1,nx
    da_x(i,j,k) = (1.0d0-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) !sigma=σ*
    da_y(i,j,k) = (1.0d0-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))!導磁率σ
    da_z(i,j,k) = (1.0d0-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))

    db_x(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
    db_y(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
    db_z(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
enddo
enddo
enddo

! write(*,*) khdx(1), bh_x(1), ch_x(1)
! write(*,*) khdx(5), bh_x(5), ch_x(5)
! write(*,*) khdx(10), bh_x(10), ch_x(10)





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!field update loop!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


!  write(*,*) "hz333",hz(3,3,3)
!  write(*,*) "hy333",hy(3,3,3)
!  write(*,*) "hx333",hx(3,3,3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!psi update!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!x-PML loop
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

!y-PML loop
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
!z-PML loop
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
            end subroutine CPML_H
