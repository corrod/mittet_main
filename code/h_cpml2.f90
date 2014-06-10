!Convolutional PML_H ******************************************************************
! sigmamax amax kappamax の求め方 sigma*の求め方
! 導入の仕方
! メイン部分の計算は通常の電磁波伝播と同じでプサイのぶぶんだけCPML？
! subrouitne 分ける必要ないのかも
!**************************************************************************************
subroutine CPML_H(ex,ey,ez,hx,hy,hz)
    use const_para
    implicit none

    integer,parameter :: m=4, ma=4
    integer :: nxpml1=5, nypml1=5, nzpml1=5 !pmlの厚さ
    real(8),parameter :: lnR0= -100d0  !ln|R(0)|
    real(8)           :: sigma_opt
    real(8)           :: sigma_max!!!
    real(8),parameter :: kappa_max=1.0d0 !!!
    real(8),parameter :: a_max=0.1d0     !!!
    real(8) :: epsir = 1.0d0
    real(8) :: myu(nx,ny,nz),epsi(nx,ny,nz)
    real(8) :: da_x(nx,ny,nz),da_y(nx,ny,nz),da_z(nx,ny,nz)
    real(8) :: db_x(nx,ny,nz),db_y(nx,ny,nz),db_z(nx,ny,nz)
    real(8) :: khdx(nx),khdy(ny),khdz(nz)
    real(8) :: bh_x(nx),bh_y(ny),bh_z(dz)
    real(8) :: ch_x(nx),ch_y(ny),ch_z(nz)
    real(8) :: sigma(nx,ny,nz)
    real(8) :: sigma_x(nx),sigma_y(ny),sigma_z(nz)
    real(8) :: kappa_x(nx),kappa_y(ny),kappa_z(nz)
    real(8) :: a_x(nx),a_y(ny),a_z(nz)
    complex(kind(0d0)) :: psi_hzx1(nx,ny,nz),psi_hyx1(nx,ny,nz)
    complex(kind(0d0)) :: psi_hxy1(nx,ny,nz),psi_hzy1(nx,ny,nz)
    complex(kind(0d0)) :: psi_hyz1(nx,ny,nz),psi_hxz1(nx,ny,nz)
    complex(kind(0d0)), intent(in) :: ex(nx,ny,nz),ey(nx,ny,nz),ez(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: hx(nx,ny,nz),hy(nx,ny,nz),hz(nx,ny,nz)

!!!    sigma_max = -(m+1)*lnR0 / (2.0d0*(sqrt(myu/epsi))*nxpml1*dx)  !ln(R(0));反射係数!!!
    sigma_opt = (m+1)/(150d0*pai*sqrt(epsir)*dx)
    sigma_max = 0.7d0*sigma_opt
    
    do k=1,nz
    do j=1,nz
    do i=1,nz

    sigma_x(i) = sigma_max* ((nxpml1-i)/(nxpml1-1)) **m  !!!-i-1/2の取り扱い
    sigma_y(j) = sigma_max* ((nypml1-j)/(nypml1-1)) **m
    sigma_z(k) = sigma_max* ((nzpml1-k)/(nzpml1-1)) **m    
!!!    kappa_max =    !!!導出要確認
    kappa_x(i) = 1 + (kappa_max-1)*((nxpml1-i)/(nxpml1-1)) **m  !!!-i-1/2の取扱い
    kappa_y(j) = 1 + (kappa_max-1)*((nypml1-j)/(nypml1-1)) **m
    kappa_z(k) = 1 + (kappa_max-1)*((nzpml1-k)/(nzpml1-1)) **m

!!!    a_max =         !!導出要確認
    a_x(i) = a_max* (i/(nxpml1-1))**ma !!!-i-1/2の取り扱い
    a_y(j) = a_max* (j/(nypml1-1))**ma
    a_z(k) = a_max* (k/(nzpml1-1))**ma

    khdx(i) = kappa_x(i)*dx !!!(i-1/2)dxの取り扱い
    khdy(j) = kappa_y(j)*dy
    khdz(k) = kappa_z(k)*dz

!係数の設定
    da_x(i,j,k) = (1.0d0-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) !sigma=σ*
    da_y(i,j,k) = (1.0d0-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))!導磁率σ
    da_z(i,j,k) = (1.0d0-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))

    db_x(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
    db_y(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
    db_z(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))

    bh_x(i) = exp(-(sigma_x(i)/kappa_x(i)+a_x(i)) *dt/epsi0)
    bh_y(j) = exp(-(sigma_y(j)/kappa_y(j)+a_y(j)) *dt/epsi0)
    bh_z(k) = exp(-(sigma_z(k)/kappa_z(k)+a_z(k)) *dt/epsi0)

    ch_x(i) = sigma_x(j)*(bh_x(i)-1) / (sigma_x(i) + kappa_x(i)*a_x(i)) / kappa_x(i)
    ch_y(j) = sigma_y(j)*(bh_y(j)-1) / (sigma_y(j) + kappa_y(j)*a_y(j)) / kappa_y(j)
    ch_z(j) = sigma_z(k)*(bh_z(k)-1) / (sigma_z(k) + kappa_z(k)*a_z(k)) / kappa_z(k)
    enddo
    enddo
    enddo




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!field update loop!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !Hx
! do k=1,nz-1
 !       do j=1,ny-1
  !          do i=2,nx-1
   !             hx(i,j,k) = da_x(i,j,k)*hx(i,j,k)&
    !                       -db_x(i,j,k)*((ez(i,j+1,k)-ez(i,j,k))/khdy(j)&
     !                                   -(ey(i,j,k+1)-ey(i,j,k))/khdz(k))
      !                     enddo
       !                         enddo
        !                            enddo
 !Hy
 !do k=1,nz-1
  !      do j=2,ny-1
   !         do i=1,nx-1
    !            hy(i,j,k) = da_y(i,j,k)*hy(i,j,k)&
     !                      -db_y(i,j,k)*((ex(i,j,k+1)-ex(i,j,k))/khdz(k)&
      !                                  -(ez(i+1,j,k)-ez(i,j,k))/khdx(i))
       !                     enddo
        !                        enddo
         !                           enddo
 !Hz
! do k=2,nz-1
 !       do j=1,ny-1
  !          do i=1,nx-1
   !             hz(i,j,k) = da_z(i,j,k)*hz(i,j,k)&
    !                       -db_z(i,j,k)*((ey(i+1,j,k)-ey(i,j,k))/khdx(i)&
     !                                   -(ex(i,j+1,k)-ex(i,j,k))/khdy(j))
      !                      enddo
       !                         enddo
        !                            enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!psi update!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!x-PML loop
    do k=1,nz-1
        do j=1,ny-1
            do i=1,nxpml1-1
                psi_hzx1(i,j,k) = bh_x(i)*psi_hzx1(i,j,k)&
                                 +ch_x(i)*(ey(i+1,j,k)-ey(i,j,k)) / dx
                psi_hyx1(i,j,k) = bh_x(i)*psi_hyx1(i,j,k)&
                                 +ch_x(i)*(ez(i+1,j,k)-ez(i,j,k)) / dx
                hz(i,j,k) = hz(i,j,k)-db_z(i,j,k)*psi_hzx1(i,j,k)
                hy(i,j,k) = hy(i,j,k)+db_y(i,j,k)*psi_hyx1(i,j,k)
                   enddo
                        enddo
                            enddo

!y-PML loop
    do k=1,nz-1
        do j=1,nypml1-1
            do i=1,nx-1
                psi_hxy1(i,j,k) = bh_y(j)*psi_hxy1(i,j,k)&
                                 +ch_y(j)*(ez(i,j+1,k)-ez(i,j,k)) / dy
                psi_hzy1(i,j,k) = bh_y(j)*psi_hzy1(i,j,k)&
                                 +ch_y(j)*(ex(i,j+1,k)-ex(i,j,k)) / dy
                hx(i,j,k) = hx(i,j,k)-db_x(i,j,k)*psi_hxy1(i,j,k)
                hz(i,j,k) = hz(i,j,k)+db_z(i,j,k)*psi_hzy1(i,j,k)
                  enddo
                       enddo
                           enddo
!z-PML loop
    do k=1,nzpml1-1
        do j=1,ny-1
            do i=1,nx-1
                psi_hyz1(i,j,k) = bh_z(k)*psi_hyz1(i,j,k)&
                                 +ch_z(k)*(ex(i,j,k+1)-ex(i,j,k)) / dz
                psi_hxz1(i,j,k) = bh_y(k)*psi_hxz1(i,j,k)&
                                 +ch_z(k)*(ey(i,j,k+1)-ey(i,j,k)) / dz
                hy(i,j,k) = hy(i,j,k)-db_y(i,j,k)*psi_hyz1(i,j,k)
                hx(i,j,k) = hx(i,j,k)+db_x(i,j,k)*psi_hxz1(i,j,k)
                   enddo
                        enddo
                            enddo
            end subroutine CPML_H
