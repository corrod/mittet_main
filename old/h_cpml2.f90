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
!epsi0=1 from imamu
!sigma* の求め方
!db_z[ijk] = dt/MU0 /(1.f+(msigz[k]*dt)/(2.f*eps2));
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CPML_H(Ex,Ey,Ez,Hx,Hy,Hz,sigma,myu,cmax)
    use const_para
    implicit none

  !  integer, parameter :: nxpml1 = 10, nypml1 = 10, nzpml1 = 10 !pmlの厚さ
    integer, parameter  :: m         = 4, ma = 1 !mの代わりにnn
    real(8), parameter  :: kappa_max = 1.0d0 !!!
    real(8), parameter  :: a_max     = 0.2d0     !!!
    real(8), parameter  :: nn        = 3.0d0 !nn should be [2,6]
    real(8), parameter  :: order     = 0.0d0 !order should be (0,3]
    real(8), parameter  :: optToMax  = 10.0d0
    real(8), parameter  :: Rcoef     = 0.01d0 !R should be [10^-2, 10^-12]
    real(8), parameter  :: c1        = 1.125d0, c2 = -0.04167d0 !pml4の係数 from taylor expansion
    real(8),parameter   :: epsir     = 1.0d0
    real(8)             :: delta     = nxpml1*dx
    real(8), intent(in) :: cmax
    real(8), intent(in) :: myu(nx,ny,nz)
    real(8), intent(in) :: sigma(nx,ny,nz)
    real(8)             :: sigma_opt
    real(8)             :: sigma_max!!!
    real(8)             :: msigma_max
    real(8)             :: msigma_x(nx),msigma_y(ny),msigma_z(nz)
    !     real(8), parameter :: lnR0 = -100.0d0  !ln|R(0)|
    real(8)             :: epsi(nx,ny,nz)!1.0d0
    real(8)             :: am_x(nx),am_y(ny),am_z(nz)
    real(8)             :: mkappa_x(nx),mkappa_y(ny),mkappa_z(nz)
    real(8)             :: khdx(nx),khdy(ny),khdz(nz)
    real(8)             :: bh_x(nx),bh_y(ny),bh_z(nz)
    real(8)             :: ch_x(nx),ch_y(ny),ch_z(nz)
    real(8)             :: da_x(nx,ny,nz),da_y(nx,ny,nz),da_z(nx,ny,nz)
    real(8)             :: db_x(nx,ny,nz),db_y(nx,ny,nz),db_z(nx,ny,nz)
    complex(kind(0d0))  :: psi_Hzx1(nx,ny,nz),psi_Hyx1(nx,ny,nz)
    complex(kind(0d0))  :: psi_Hxy1(nx,ny,nz),psi_Hzy1(nx,ny,nz)
    complex(kind(0d0))  :: psi_Hyz1(nx,ny,nz),psi_Hxz1(nx,ny,nz)
    complex(kind(0d0)), intent(in)    :: Ex(1:nx,1:ny,1:nz),Ey(1:nx,1:ny,1:nz),Ez(1:nx,1:ny,1:nz)
    complex(kind(0d0)), intent(inout) :: Hx(1:nx,1:ny,1:nz),Hy(1:nx,1:ny,1:nz),Hz(1:nx,1:ny,1:nz)

   !Holberg optimization scheme
  !  alpha(1,1)   = 1.00235d0
    !alpha(2,1:2) = (/1.14443d0,-0.04886d0/)
    !alpha(3,1:3) = (/1.20282d0,-0.08276d0,0.00950d0/)
    !alpha(4,1:4) = (/1.23041d0,-0.10313d0,0.02005d0,-0.00331d0/)

    !Taylor expansion
!     alpha(1,1)   = 1.0d0
!     alpha(2,1:2) = (/1.12500d0,-0.04167d0/)
!     alpha(3,1:3) = (/1.17188d0,-0.06510d0,0.00469d0/)
!     alpha(4,1:4) = (/1.19629d0,-0.07975d0,-0.00070d0/)

    epsi(1:nx,1:ny,1:nz)=sigma(1:nx,1:ny,1:nz)/2.0d0/omega0

!!!    sigma_max = -(m+1)*lnR0 / (2.0d0*(sqrt(myu/epsi))*nxpml1*dx)  !ln(R(0));反射係数!!!
    sigma_max = (nn+order+1.0d0)*cmax*log(1.0d0/Rcoef) / (2.0d0*delta) * optToMax  !!x方向だけ？
    !msigma_max = sigma_max*MU0/epsi2

   ! sigma_opt = (dble(m)+1.0d0) / (150.0d0*pai*sqrt(epsir)*dx)
  !  sigma_max = 0.7d0*sigma_opt
    

!係数の設定xh

do i = 1,nx
if (i<=nxpml1) then
    msigma_x(i) = sigma_max* ((dble(nxpml1)-dble(i)-0.5d0)/(dble(nxpml1)-1.0d0))**dble(nn)  !!!-i-1/2の取り扱い
    mkappa_x(i) = 1.0d0 + (kappa_max-1.0d0)*((dble(nxpml1)-dble(i)-0.5d0)/(dble(nxpml1)-1.0d0))**dble(nn)  !!!-i-1/2の取扱い
    am_x(i)     = a_max* ((dble(i)-0.5d0)/(dble(nxpml1)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

    bh_x(i)     = exp(-(msigma_x(i)/mkappa_x(i)+am_x(i)) *dt) !/epsi0)
    ch_x(i)     = msigma_x(i)*(bh_x(i)-1.0d0) / (msigma_x(i) + mkappa_x(i)*am_x(i)) / mkappa_x(i)
    khdx(i)     = mkappa_x(i)*dx !!!(i-1/2)dxの取り扱い


else if(i>=nx-nxpml1+1) then
    msigma_x(i) = sigma_max* ((dble(i)-dble(nx)+0.5d0+dble(nxpml1))/(dble(nxpml1)-1.0d0))**dble(nn)  !!!-i-1/2の取り扱い
    mkappa_x(i) = 1.0d0 + (kappa_max-1.0d0)*((dble(i)-dble(nx)+0.5d0+dble(nxpml1))/(dble(nxpml1)-1.0d0))**dble(nn)  !!!-i-1/2の取扱い
    am_x(i)     = a_max* ((dble(-i)+nx+0.5d0)/(dble(nxpml1)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

    bh_x(i)     = exp(-(msigma_x(i)/mkappa_x(i)+am_x(i)) *dt) !/epsi0)
    ch_x(i)     = msigma_x(i)*(bh_x(i)-1.0d0) / (msigma_x(i) + mkappa_x(i)*am_x(i)) / mkappa_x(i) 
    khdx(i)     = mkappa_x(i)*dx !!!(i-1/2)dxの取り扱い

else
    msigma_x(i) = 0.0d0
    mkappa_x(i) = 1.0d0
    am_x(i)     = 0.0d0   
    bh_x(i)     = 0.0d0   
    ch_x(i)     = 0.0d0   
    khdx(i)     = mkappa_x(i)*dx
    endif
        enddo


!係数の設定yh

do j = 1,ny
if (j<=nypml1) then
    msigma_y(j) = sigma_max* ((dble(nypml1)-dble(j)-0.5d0)/(dble(nypml1)-1.0d0))**dble(nn)  !!!-i-1/2の取り扱い
    mkappa_y(j) = 1.0d0 + (kappa_max-1.0d0)*((dble(nypml1)-dble(j)-0.5d0)/(dble(nypml1)-1.0d0))**dble(nn)  !!!-i-1/2の取扱い
    am_y(j)     = a_max* ((dble(j)-0.5d0)/(dble(nypml1)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

    bh_y(j)     = exp(-(msigma_y(j)/mkappa_y(j)+am_y(j)) *dt) !/epsi0)
    ch_y(j)     = msigma_y(j)*(bh_y(j)-1.0d0) / (msigma_y(j) + mkappa_y(j)*am_y(j)) / mkappa_y(j)
    khdy(j)     = mkappa_y(j)*dx !!!(i-1/2)dxの取り扱い

else if(j>=ny-nypml1+1) then
    msigma_y(j) = sigma_max* ((dble(j)-dble(ny)+0.5d0+dble(nypml1))/(dble(nypml1)-1.0d0))**dble(nn)  !!!-i-1/2の取り扱い
    mkappa_y(j) = 1.0d0 + (kappa_max-1.0d0)*((dble(j)-dble(ny)+0.5d0+dble(nypml1))/(dble(nypml1)-1.0d0))**dble(nn)  !!!-i-1/2の取扱い
    am_y(j)     = a_max* ((dble(-j)+ny+0.5d0)/(dble(nypml1)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

    bh_y(j)     = exp(-(msigma_y(j)/mkappa_y(j)+am_y(j)) *dt) !/epsi0)
    ch_y(j)     = msigma_y(j)*(bh_y(j)-1.0d0) / (msigma_y(j) + mkappa_y(j)*am_y(j)) / mkappa_y(j) 
    khdy(j)     = mkappa_y(j)*dy !!!(i-1/2)dxの取り扱い

else
    msigma_y(j) = 0.0d0
    mkappa_y(j) = 1.0d0
    am_y(j)     = 0.0d0   
    bh_y(j)     = 0.0d0   
    ch_y(j)     = 0.0d0   
    khdy(j)     = mkappa_y(j)*dy
    endif
        enddo


!係数の設定zh

do k = 1,nz
if (k<=nzpml1) then
    msigma_z(k) = sigma_max* ((dble(nzpml1)-dble(k)-0.5d0)/(dble(nzpml1)-1.0d0))**dble(nn)  !!!-i-1/2の取り扱い
    mkappa_z(k) = 1.0d0 + (kappa_max-1.0d0)*((dble(nzpml1)-dble(k)-0.5d0)/(dble(nzpml1)-1.0d0))**dble(nn)  !!!-i-1/2の取扱い
    am_z(k)     = a_max* ((dble(k)-0.5d0)/(dble(nzpml1)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

    bh_z(k)     = exp(-(msigma_z(k)/mkappa_z(k)+am_z(k)) *dt) !/epsi0)
    ch_z(k)     = msigma_z(k)*(bh_z(k)-1.0d0) / (msigma_z(k) + mkappa_z(k)*am_z(k)) / mkappa_z(k)
    khdz(k)     = mkappa_z(k)*dz !!!(i-1/2)dxの取り扱い

else if(k>=nz-nzpml1+1) then
    msigma_z(k) = sigma_max* ((dble(k)-dble(nz)+0.5d0+dble(nzpml1))/(dble(nzpml1)-1.0d0))**dble(nn)  !!!-i-1/2の取り扱い
    mkappa_z(k) = 1.0d0 + (kappa_max-1.0d0)*((dble(k)-dble(nz)+0.5d0+dble(nzpml1))/(dble(nzpml1)-1.0d0))**dble(nn)  !!!-i-1/2の取扱い
    am_z(k)     = a_max* ((dble(-k)+nz+0.5d0)/(dble(nzpml1)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

    bh_z(k)     = exp(-(msigma_z(k)/mkappa_z(k)+am_z(k)) *dt) !/epsi0)
    ch_z(k)     = msigma_z(k)*(bh_z(k)-1.0d0) / (msigma_z(k) + mkappa_z(k)*am_z(k)) / mkappa_z(k) 
    khdz(k)     = mkappa_z(k)*dz !!!(i-1/2)dxの取り扱い

else
    msigma_z(k) = 0.0d0
    mkappa_z(k) = 1.0d0
    am_z(k)     = 0.0d0   
    bh_z(k)     = 0.0d0   
    ch_z(k)     = 0.0d0   
    khdz(k)     = mkappa_z(k)*dz
    endif
        enddo

!!!sigma=simga*

do k=1,nz
    do j=1,ny
        do i=1,nx
        !imamu system
        da_x(i,j,k) = (1.0d0 - ((msigma_x(i)*dt)/(2.0d0*epsi(i,j,k)))) / (1.0d0 + ((msigma_x(i)*dt)/(2.0d0*epsi(i,j,k))))
        da_y(i,j,k) = (1.0d0 - ((msigma_y(j)*dt)/(2.0d0*epsi(i,j,k)))) / (1.0d0 + ((msigma_y(j)*dt)/(2.0d0*epsi(i,j,k))))
        da_z(i,j,k) = (1.0d0 - ((msigma_z(k)*dt)/(2.0d0*epsi(i,j,k)))) / (1.0d0 + ((msigma_z(k)*dt)/(2.0d0*epsi(i,j,k))))
        !!!****imamu systemではmyu→MU0 myu→eps2
        db_x(i,j,k) = dt/MU0/(1.0d0+(msigma_x(i)*dt)/(2.0d0*epsi(i,j,k)))
        db_y(i,j,k) = dt/MU0/(1.0d0+(msigma_y(j)*dt)/(2.0d0*epsi(i,j,k)))
        db_z(i,j,k) = dt/MU0/(1.0d0+(msigma_z(k)*dt)/(2.0d0*epsi(i,j,k)))

         !saito system   
    !   da_x(i,j,k) = (1.0d0-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) !sigma=σ*
    !   da_y(i,j,k) = (1.0d0-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))!導磁率σ
    !   da_z(i,j,k) = (1.0d0-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
    !   saito system
    !   db_x(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
    !   db_y(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
    !   db_z(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
!!!sigma_x(i)?sigma(i,j,k)?
        enddo
    enddo
enddo


! eps2 = sig2[ijk] /2.f /omega_0;
!            // CPML Coefficient
!            // coef1
!            ca_x[ijk]   = (1.0f - ((esigx[i]*dt)/(2.f*eps2))) \
!                        / (1.0f + ((esigx[i]*dt)/(2.f*eps2)));
!            ca_y[ijk]   = (1.0f - ((esigy[j]*dt)/(2.f*eps2))) \
!                        / (1.0f + ((esigy[j]*dt)/(2.f*eps2)));
!            ca_z[ijk]   = (1.0f - ((esigz[k]*dt)/(2.f*eps2))) \
!                        / (1.0f + ((esigz[k]*dt)/(2.f*eps2)));
!            // coef2
!            da_x[ijk]   = (1.0f - ((msigx[i]*dt)/(2.f*eps2))) \
!                        / (1.0f + ((msigx[i]*dt)/(2.f*eps2)));
!            da_y[ijk]   = (1.0f - ((msigy[j]*dt)/(2.f*eps2))) \
!                        / (1.0f + ((msigy[j]*dt)/(2.f*eps2)));
!            da_z[ijk]   = (1.0f - ((msigz[k]*dt)/(2.f*eps2))) \
!                        / (1.0f + ((msigz[k]*dt)/(2.f*eps2)));
!            // coef3
!            cb_x[ijk] = dt/eps2 /(1.f+(esigx[i]*dt)/(2.f*eps2));
!            cb_y[ijk] = dt/eps2 /(1.f+(esigy[j]*dt)/(2.f*eps2));
!            cb_z[ijk] = dt/eps2 /(1.f+(esigz[k]*dt)/(2.f*eps2));
!            // coef4
!            db_x[ijk] = dt/MU0 /(1.f+(msigx[i]*dt)/(2.f*eps2));
!            db_y[ijk] = dt/MU0 /(1.f+(msigy[j]*dt)/(2.f*eps2));
!            db_z[ijk] = dt/MU0 /(1.f+(msigz[k]*dt)/(2.f*eps2));


!
!msig_max=esig_max*MU0/epsi2
!  scaler = 0.01f * sigmax2/gradmax;
!  //scaler = 0.01f * sigmax2;
!sig2[ijk] = sig[ijk] - scaler * grad[ijk];


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





!44444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
!xh-PML4 loop(-)
    do k = 1,nz-1
        do j = 1,ny-1
            do i = 2,nxpml1-1
                psi_Hzx1(i,j,k) = bh_x(i)*psi_Hzx1(i,j,k)&
                                 +ch_x(i)*( c1*Ey(i+1,j,k)-c1*Ey(i,j,k) + c2*Ey(i+2,j,k)-c2*Ey(i-1,j,k) ) / dx
                psi_Hyx1(i,j,k) = bh_x(i)*psi_Hyx1(i,j,k)&
                                 +ch_x(i)*( c1*Ez(i+1,j,k)-c2*Ez(i,j,k) + c2*Ez(i+2,j,k)-c2*Ez(i-1,j,k) ) / dx
                Hz(i,j,k) = Hz(i,j,k) - db_z(i,j,k)*psi_Hzx1(i,j,k)
                Hy(i,j,k) = Hy(i,j,k) + db_y(i,j,k)*psi_Hyx1(i,j,k)
                   enddo
                        enddo
                            enddo
!xh-PML4 loop(+)
    do k = 1,nz-1
        do j = 1,ny-1
            do i = nx-nxpml1+1,nx-2
                psi_Hzx1(i,j,k) = bh_x(i)*psi_Hzx1(i,j,k)&
                                 +ch_x(i)*( c1*Ey(i+1,j,k)-c1*Ey(i,j,k) + c2*Ey(i+2,j,k)-c2*Ey(i-1,j,k) ) / dx
                psi_Hyx1(i,j,k) = bh_x(i)*psi_Hyx1(i,j,k)&
                                 +ch_x(i)*( c1*Ez(i+1,j,k)-c2*Ez(i,j,k) + c2*Ez(i+2,j,k)-c2*Ez(i-1,j,k) ) / dx
                Hz(i,j,k) = Hz(i,j,k) - db_z(i,j,k)*psi_Hzx1(i,j,k)
                Hy(i,j,k) = Hy(i,j,k) + db_y(i,j,k)*psi_Hyx1(i,j,k)
                   enddo
                        enddo
                            enddo

!yh-PML4 loop(-)
    do k = 1,nz-1
        do j = 2,nypml1-1
            do i = 1,nx-1
                psi_Hxy1(i,j,k) = bh_y(j)*psi_Hxy1(i,j,k)&
                                 +ch_y(j)*( c1*Ez(i,j+1,k)-c1*Ez(i,j,k) + c2*Ez(i,j+2,k)-c2*Ez(i,j-1,k) ) / dy
                psi_Hzy1(i,j,k) = bh_y(j)*psi_Hzy1(i,j,k)&
                                 +ch_y(j)*( c1*Ex(i,j+1,k)-c1*Ex(i,j,k) + c2*Ex(i,j+2,k)-c2*Ex(i,j-1,k) ) / dy
                Hx(i,j,k) = Hx(i,j,k) - db_x(i,j,k)*psi_Hxy1(i,j,k)
                Hz(i,j,k) = Hz(i,j,k) + db_z(i,j,k)*psi_Hzy1(i,j,k)
                  enddo
                       enddo
                           enddo
!!yh-PML4 loop(+)
    do k = 1,nz-1
        do j = ny-nypml1+1,ny-2
            do i = 1,nx-1
                psi_Hxy1(i,j,k) = bh_y(j)*psi_Hxy1(i,j,k)&
                                 +ch_y(j)*( c1*Ez(i,j+1,k)-c1*Ez(i,j,k) + c2*Ez(i,j+2,k)-c2*Ez(i,j-1,k) ) / dy
                psi_Hzy1(i,j,k) = bh_y(j)*psi_Hzy1(i,j,k)&
                                 +ch_y(j)*( c1*Ex(i,j+1,k)-c1*Ex(i,j,k) + c2*Ex(i,j+2,k)-c2*Ex(i,j-1,k) ) / dy
                Hx(i,j,k) = Hx(i,j,k) - db_x(i,j,k)*psi_Hxy1(i,j,k)
                Hz(i,j,k) = Hz(i,j,k) + db_z(i,j,k)*psi_Hzy1(i,j,k)
                  enddo
                       enddo
                           enddo

!zh-PML4 loop(-)
    do k = 2,nzpml1-1
        do j = 1,ny-1
            do i = 1,nx-1
                psi_Hyz1(i,j,k) = bh_z(k)*psi_Hyz1(i,j,k)&
                                 +ch_z(k)*( c1*Ex(i,j,k+1)-c1*Ex(i,j,k) + c2*Ex(i,j,k+2)-c2*Ex(i,j,k-1) ) / dz
                psi_Hxz1(i,j,k) = bh_y(k)*psi_Hxz1(i,j,k)&
                                 +ch_z(k)*( c1*Ey(i,j,k+1)-c1*Ey(i,j,k) + c2*Ey(i,j,k+2)-c2*Ey(i,j,k-1) ) / dz
                Hy(i,j,k) = Hy(i,j,k) - db_y(i,j,k)*psi_Hyz1(i,j,k)
                Hx(i,j,k) = Hx(i,j,k) + db_x(i,j,k)*psi_Hxz1(i,j,k)
                   enddo
                        enddo
                            enddo
!zh-PML4 loop(+)
    do k = nz-nzpml1+1,nz-2
        do j = 1,ny-1
            do i = 1,nx-1
                psi_Hyz1(i,j,k) = bh_z(k)*psi_Hyz1(i,j,k)&
                                 +ch_z(k)*( c1*Ex(i,j,k+1)-c1*Ex(i,j,k) + c2*Ex(i,j,k+2)-c2*Ex(i,j,k-1) ) / dz
                psi_Hxz1(i,j,k) = bh_y(k)*psi_Hxz1(i,j,k)&
                                 +ch_z(k)*( c1*Ey(i,j,k+1)-c1*Ey(i,j,k) + c2*Ey(i,j,k+2)-c2*Ey(i,j,k-1) ) / dz
                Hy(i,j,k) = Hy(i,j,k) - db_y(i,j,k)*psi_Hyz1(i,j,k)
                Hx(i,j,k) = Hx(i,j,k) + db_x(i,j,k)*psi_Hxz1(i,j,k)
                   enddo
                        enddo
                            enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!psi update!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!xh-PML loop(-)
!    do k = 1,nz-1
!        do j = 1,ny-1
!            do i = 1,nxpml1-1
!                psi_Hzx1(i,j,k) = bh_x(i)*psi_Hzx1(i,j,k)&
!                                 +ch_x(i)*(Ey(i+1,j,k)-Ey(i,j,k)) / dx
!                psi_Hyx1(i,j,k) = bh_x(i)*psi_Hyx1(i,j,k)&
!                                 +ch_x(i)*(Ez(i+1,j,k)-Ez(i,j,k)) / dx
!                Hz(i,j,k) = Hz(i,j,k) - db_z(i,j,k)*psi_Hzx1(i,j,k)
!                Hy(i,j,k) = Hy(i,j,k) + db_y(i,j,k)*psi_Hyx1(i,j,k)
!                   enddo
!                        enddo
!                            enddo
!!xh-PML loop(+)
!    do k = 1,nz-1
!        do j = 1,ny-1
!            do i = nx-nxpml1+1,nx-1
!                psi_Hzx1(i,j,k) = bh_x(i)*psi_Hzx1(i,j,k)&
!                                 +ch_x(i)*(Ey(i+1,j,k)-Ey(i,j,k)) / dx
!                psi_Hyx1(i,j,k) = bh_x(i)*psi_Hyx1(i,j,k)&
!                                 +ch_x(i)*(Ez(i+1,j,k)-Ez(i,j,k)) / dx
!                Hz(i,j,k) = Hz(i,j,k) - db_z(i,j,k)*psi_Hzx1(i,j,k)
!                Hy(i,j,k) = Hy(i,j,k) + db_y(i,j,k)*psi_Hyx1(i,j,k)
!                   enddo
!                        enddo
!                            enddo
!
!
!
!
!
!!yh-PML loop(-)
!    do k = 1,nz-1
!        do j = 1,nypml1-1
!            do i = 1,nx-1
!                psi_Hxy1(i,j,k) = bh_y(j)*psi_Hxy1(i,j,k)&
!                                 +ch_y(j)*(Ez(i,j+1,k)-Ez(i,j,k)) / dy
!                psi_Hzy1(i,j,k) = bh_y(j)*psi_Hzy1(i,j,k)&
!                                 +ch_y(j)*(Ex(i,j+1,k)-Ex(i,j,k)) / dy
!                Hx(i,j,k) = Hx(i,j,k) - db_x(i,j,k)*psi_Hxy1(i,j,k)
!                Hz(i,j,k) = Hz(i,j,k) + db_z(i,j,k)*psi_Hzy1(i,j,k)
!                  enddo
!                       enddo
!                           enddo
!!!yh-PML loop(+)
!    do k = 1,nz-1
!        do j = ny-nypml1+1,ny-1
!            do i = 1,nx-1
!                psi_Hxy1(i,j,k) = bh_y(j)*psi_Hxy1(i,j,k)&
!                                 +ch_y(j)*(Ez(i,j+1,k)-Ez(i,j,k)) / dy
!                psi_Hzy1(i,j,k) = bh_y(j)*psi_Hzy1(i,j,k)&
!                                 +ch_y(j)*(Ex(i,j+1,k)-Ex(i,j,k)) / dy
!                Hx(i,j,k) = Hx(i,j,k) - db_x(i,j,k)*psi_Hxy1(i,j,k)
!                Hz(i,j,k) = Hz(i,j,k) + db_z(i,j,k)*psi_Hzy1(i,j,k)
!                  enddo
!                       enddo
!                           enddo
!
!
!
!
!!zh-PML loop(-)
!    do k = 1,nzpml1-1
!        do j = 1,ny-1
!            do i = 1,nx-1
!                psi_Hyz1(i,j,k) = bh_z(k)*psi_Hyz1(i,j,k)&
!                                 +ch_z(k)*(Ex(i,j,k+1)-Ex(i,j,k)) / dz
!                psi_Hxz1(i,j,k) = bh_y(k)*psi_Hxz1(i,j,k)&
!                                 +ch_z(k)*(Ey(i,j,k+1)-Ey(i,j,k)) / dz
!                Hy(i,j,k) = Hy(i,j,k) - db_y(i,j,k)*psi_Hyz1(i,j,k)
!                Hx(i,j,k) = Hx(i,j,k) + db_x(i,j,k)*psi_Hxz1(i,j,k)
!                   enddo
!                        enddo
!                            enddo
!!zh-PML loop(+)
!    do k = nz-nzpml1+1,nz-1
!        do j = 1,ny-1
!            do i = 1,nx-1
!                psi_Hyz1(i,j,k) = bh_z(k)*psi_Hyz1(i,j,k)&
!                                 +ch_z(k)*(Ex(i,j,k+1)-Ex(i,j,k)) / dz
!                psi_Hxz1(i,j,k) = bh_y(k)*psi_Hxz1(i,j,k)&
!                                 +ch_z(k)*(Ey(i,j,k+1)-Ey(i,j,k)) / dz
!                Hy(i,j,k) = Hy(i,j,k) - db_y(i,j,k)*psi_Hyz1(i,j,k)
!                Hx(i,j,k) = Hx(i,j,k) + db_x(i,j,k)*psi_Hxz1(i,j,k)
!                   enddo
!                        enddo
!                            enddo
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!field update loop!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  !Hx
!  do k=1,nz-1
!        do j=1,ny-1
!            do i=2,nx-1
!                Hx(i,j,k) = da_x(i,j,k)*Hx(i,j,k)&
!                           -db_x(i,j,k)*((Ez(i,j+1,k)-Ez(i,j,k))/khdy(j)&
!                                        -(Ey(i,j,k+1)-Ey(i,j,k))/khdz(k))
!                           enddo
!                                enddo
!                                    enddo
!  !Hy
!  do k=1,nz-1
!        do j=2,ny-1
!            do i=1,nx-1
!                Hy(i,j,k) = da_y(i,j,k)*Hy(i,j,k)&
!                           -db_y(i,j,k)*((Ex(i,j,k+1)-Ex(i,j,k))/khdz(k)&
!                                        -(Ez(i+1,j,k)-Ez(i,j,k))/khdx(i))
!                            enddo
!                                enddo
!                                    enddo
!  !Hz
! do k=2,nz-1
!        do j=1,ny-1
!            do i=1,nx-1
!                Hz(i,j,k) = da_z(i,j,k)*Hz(i,j,k)&
!                           -db_z(i,j,k)*((Ey(i+1,j,k)-Ey(i,j,k))/khdx(i)&
!                                        -(Ex(i,j+1,k)-Ex(i,j,k))/khdy(j))
!                            enddo
!                                enddo
!                                    enddo
            end subroutine CPML_H
