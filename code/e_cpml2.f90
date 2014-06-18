!!!Convolutional PML_E !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sigmax amax kappamax の求め方
! メイン部分の計算は通常の電磁波伝播と同じでプサイのぶぶんだけCPML？
! subrouitne 分ける必要ないのかも
!psi部分だけPMLバージョン。
!kappa >=1 real
!sigi>0 real
!ai=alphai>0 real
!psi loop あと３枚必要？
!演算中の倍精度の表現に注意.d0,
!be_x(5)が0になってしまう問題.１になるはず,epsi0=8.854d-12
!epsi0が小さすぎる!!:q
!cb_x,cb_y,cb_zを仮想領域にあわせる
!係数が1-nxpml1, nx-(nx-nxpml)でことなるので調整
!epsi0=1 from imamu
!sig_xの値をどう取るか
!sig(i,j,k) or sig_x(i) ??
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CPML_E(Ex,Ey,Ez,Hx,Hy,Hz,sig,cmax)
    use const_para
    implicit none

   ! integer, parameter :: nxpml1 = 10, nypml1 = 10, nzpml1 = 10 !PML層数
    integer, parameter  :: m         = 4, ma = 1
    real(8), parameter  :: kappa_max = 1.0d0 !!!kappa should be [0,10]
    real(8), parameter  :: a_max     = 0.2d0     !!!if a_max=0, CPML changes to UPML
    real(8), parameter  :: nn        = 3.0d0 !nn should be [2,6]
    real(8), parameter  :: order     = 0.0d0 !order should be (0,3]
    real(8), parameter  :: optToMax  = 10.0d0
    real(8), parameter  :: Rcoef     = 0.01d0 !R should be [10^-2, 10^-12]
    real(8), parameter  :: c1        = 1.125d0, c2 = -0.04167d0 !pml4の係数 from taylor Expansion
    real(8),parameter   :: epsir     = 1.0d0
    real(8)             :: delta = nxpml1*dx
    real(8), intent(in) :: cmax
    real(8), intent(in) :: sig(nx,ny,nz)
    real(8)             :: sig_opt !!!
    real(8)             :: sig_max !!! 導出法確認
    real(8)             :: sig_x(nx),sig_y(ny),sig_z(nz)
    real(8)             :: a_x(nx),a_y(ny),a_z(nz)
    real(8)             :: kappa_x(nx),kappa_y(ny),kappa_z(nz)
    real(8)             :: kedx(nx),kedy(ny),kedz(nz) 
    real(8)             :: epsi(nx,ny,nz)!1.0d0
    real(8)             :: ca_x(nx,ny,nz),ca_y(nx,ny,nz),ca_z(nx,ny,nz)
    real(8)             :: cb_x(nx,ny,nz),cb_y(nx,ny,nz),cb_z(nx,ny,nz)
    real(8)             :: be_x(nx),be_y(ny),be_z(nz)
    real(8)             :: ce_x(nx),ce_y(ny),ce_z(nz)
    complex(kind(0d0))  :: psi_Ezx1(nx,ny,nz),psi_Eyx1(nx,ny,nz)
    complex(kind(0d0))  :: psi_Exy1(nx,ny,nz),psi_Ezy1(nx,ny,nz)
    complex(kind(0d0))  :: psi_Eyz1(nx,ny,nz),psi_Exz1(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ex(1:nx,1:ny,1:nz),Ey(1:nx,1:ny,1:nz),Ez(1:nx,1:ny,1:nz)
    complex(kind(0d0)), intent(in)    :: Hx(1:nx,1:ny,1:nz),Hy(1:nx,1:ny,1:nz),Hz(1:nx,1:ny,1:nz)

    epsi(1:nx,1:ny,1:nz)=sig(1:nx,1:ny,1:nz)/(2.0d0*omega0)
    sig_max = (nn+order+1.0d0)*cmax*log(1.0d0/Rcoef) / (2.0d0*delta) * optToMax  !!x方向だけ？
!!!    sig_max = -(m+1)*lnR0 / (2.0d0*(sqrt(myu/epsi))*nxpml1*dx)  !ln(R(0));反射係数!!!
!    sig_opt = (dble(m)+1.0d0) / (150.0d0*pai*sqrt(epsir)*dx)
!    sig_max = 0.7d0*sig_opt


!    Holberg optimization scheme
!    alpha(1,1)   = 1.00235d0
!     alpha(2,1:2) = (/1.14443d0,-0.04886d0/)
!     alpha(3,1:3) = (/1.20282d0,-0.08276d0,0.00950d0/)
!     alpha(4,1:4) = (/1.23041d0,-0.10313d0,0.02005d0,-0.00331d0/)

!     Taylor expansion
!     alpha(1,1)   = 1.0d0
!     alpha(2,1:2) = (/1.12500d0,-0.04167d0/)
!     alpha(3,1:3) = (/1.17188d0,-0.06510d0,0.00469d0/)
!     alpha(4,1:4) = (/1.19629d0,-0.07975d0,-0.00070d0/)







!係数の設定xe
 do i = 1,nx
    if(i<=nxpml1) then
      sig_x(i) = sig_max * ((dble(nxpml1)-dble(i))/(dble(nxpml1)-1.0d0))**dble(nn+order)
      kappa_x(i) = 1.0d0 + (kappa_max-1.0d0)*( (dble(nxpml1)-dble(i))/(dble(nxpml1)-1.0d0) )**dble(nn)
      a_x(i)     = a_max * ((dble(i)-1.0d0)/(dble(nxpml1)-1.0d0))**dble(ma)

      be_x(i)    = exp(-(sig_x(i)/kappa_x(i)+a_x(i))*dt) !/epsi0)
      ce_x(i)    = sig_x(i)*(be_x(i)-1.0d0) / (sig_x(i)+kappa_x(i)*a_x(i)) / kappa_x(i)
      kedx(i)    = kappa_x(i)*dx

    else if(i>=nx-nxpml1+1) then
      sig_x(i) = sig_max * ((dble(i)-dble(nx)-1.0d0+dble(nxpml1))/(dble(nxpml1)-1.0d0))**dble(nn+order)
      kappa_x(i) = 1.0d0 + (kappa_max-1.0d0)*( (dble(i)-dble(nx)-1.0d0+dble(nxpml1))/(dble(nxpml1)-1.0d0) )**dble(nn)
      a_x(i)     = a_max * ((dble(-i)+dble(nx)  )/(dble(nxpml1)-1.0d0))**dble(ma)

      be_x(i)    = exp(-(sig_x(i)/kappa_x(i)+a_x(i))*dt) !/epsi0)
      ce_x(i)    = sig_x(i)*(be_x(i)-1.0d0) / (sig_x(i)+kappa_x(i)*a_x(i)) / kappa_x(i)
      kedx(i)    = kappa_x(i)*dx

    else
      sig_x(i) = 0.0d0
      kappa_x(i) = 1.0d0
      a_x(i)     = 0.0d0
      be_x(i)    = 0.0d0
      ce_x(i)    = 0.0d0
      kedx(i)    = kappa_x(i)*dx
    endif
enddo


!係数の設定ye

do j = 1,ny
  if(j<=nypml1) then
    sig_y(j) = sig_max * ((dble(nypml1)-dble(j))/(dble(nypml1)-1.0d0))**dble(nn+order)
    kappa_y(j) = 1.0d0 + (kappa_max-1.0d0)*( (dble(nypml1)-dble(j))/(dble(nypml1)-1.0d0) )**dble(nn)
    a_y(j)     = a_max * ((dble(j)-1.0d0)/(dble(nypml1)-1.0d0))**dble(ma)

    be_y(j)    = exp(-(sig_y(j)/kappa_y(j)+a_y(j))*dt) !/epsi0)
    ce_y(j)    = sig_y(j)*(be_y(j)-1.0d0) / (sig_y(j)+kappa_y(j)*a_y(j)) / kappa_y(j)
    kedy(j)    = kappa_y(j)*dy  !!!

  else if(j>=ny-nypml1+1) then
    sig_y(j) = sig_max * ((dble(j)-dble(ny)-1.0d0+dble(nypml1))/(dble(nypml1)-1.0d0))**dble(nn+order)
    kappa_y(j) = 1.0d0 + (kappa_max-1.0d0)*( (dble(j)-dble(ny)-1.0d0+dble(nypml1))/(dble(nypml1)-1.0d0) )**dble(nn)
    a_y(j)     = a_max * ((dble(-j)+dble(ny)  )/(dble(nypml1)-1.0d0))**dble(ma)

    be_y(j)    = exp(-(sig_y(j)/kappa_y(j)+a_y(j))*dt) !/epsi0)
    ce_y(j)    = sig_y(j)*(be_y(j)-1.0d0) / (sig_y(j)+kappa_y(j)*a_y(j)) / kappa_y(j)
    kedy(j)    = kappa_y(j)*dy  !!!

  else
    sig_y(j) = 0.0d0
    kappa_y(j) = 1.0d0
    a_y(j)     = 0.0d0
    be_y(j)    = 0.0d0
    ce_y(j)    = 0.0d0
    kedy(j)    = kappa_y(j)*dy
      endif
  enddo


!係数の設定ze

do k = 1,nz
  if(k<=nzpml1) then
    sig_z(k) = sig_max * ((dble(nzpml1)-dble(k))/(dble(nzpml1)-1.0d0))**dble(nn+order)
    kappa_z(k) = 1.0d0 + (kappa_max-1.0d0)*( (dble(nzpml1)-dble(k))/(dble(nzpml1)-1.0d0) )**dble(nn)
    a_z(k)     = a_max*((dble(k)-1.0d0)/(dble(nzpml1)-1.0d0))**dble(ma)

    be_z(k)    = exp(-(sig_z(k)/kappa_z(k)+a_z(k))*dt) !/epsi0)
    ce_z(k)    = sig_z(k)*(be_z(k)-1.0d0) / (sig_z(k)+kappa_z(k)*a_z(k)) / kappa_z(k)
    kedz(k)    = kappa_z(k)*dz  !!!

  else if(k>=nz-nzpml1+1) then
    sig_z(k) = sig_max * ((dble(k)-dble(nz)-1.0d0+dble(nzpml1))/(dble(nzpml1)-1.0d0))**dble(nn+order)
    kappa_z(k) = 1.0d0 + (kappa_max-1.0d0)*( (dble(k)-dble(nz)-1.0d0+dble(nzpml1))/(dble(nzpml1)-1.0d0) )**dble(nn)
    a_z(k)     = a_max * ((dble(-k)+dble(nz)  )/(dble(nzpml1)-1.0d0))**dble(ma)

    be_z(k)    = exp(-(sig_z(k)/kappa_z(k)+a_z(k))*dt) !/epsi0)
    ce_z(k)    = sig_z(k)*(be_z(k)-1.0d0) / (sig_z(k)+kappa_z(k)*a_z(k)) / kappa_z(k)
    kedz(k)    = kappa_z(k)*dz  !!!

  else
    sig_z(k) = 0.0d0
    kappa_z(k) = 1.0d0
    a_z(k)     = 0.0d0
    be_z(k)    = 0.0d0
    ce_z(k)    = 0.0d0
    kedz(k)    = kappa_z(k)*dz  
  endif
enddo


do k=1,nz
  do j=1,ny
    do i=1,nx
    !imamu system
    ca_x(i,j,k) = (1.0d0-sig_x(i)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sig_x(i)*dt/(2.0d0*epsi(i,j,k)))
    ca_y(i,j,k) = (1.0d0-sig_y(j)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sig_y(j)*dt/(2.0d0*epsi(i,j,k)))
    ca_z(i,j,k) = (1.0d0-sig_z(k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sig_z(k)*dt/(2.0d0*epsi(i,j,k)))
    cb_x(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+sig_x(i)*dt/(2.0d0*epsi(i,j,k)))
    cb_y(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+sig_y(j)*dt/(2.0d0*epsi(i,j,k)))
    cb_z(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+sig_z(k)*dt/(2.0d0*epsi(i,j,k)))

    !saito system
   ! ca_x(i,j,k) = (1.0d0-sig(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sig(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
   ! ca_y(i,j,k) = (1.0d0-sig(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sig(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
   ! ca_z(i,j,k) = (1.0d0-sig(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sig(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
   ! cb_x(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+sig(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
   ! cb_y(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+sig(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
   ! cb_z(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+sig(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
    enddo
  enddo
enddo

!for(i=0;i<ix;i++){
!            ijk  = k*iy*ix + j*ix + i;
!            eps2 = sig2[ijk] /2.f /omega_0;
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


! int i,j,k,ijk;
!  double sigx2,gradmax;
!  FILE *ophi;
!  sigx2 = sig[0];
!  gradmax = grad[0];
!  for(k=0;k<iz;k++){
!    for(j=0;j<iy;j++){
!      for(i=0;i<ix;i++){
!        ijk = k*ix*iy + j*ix + i;
!        if(k<=iseabed){
!          if(sig[ijk]  > sigx2)  sigx2 = sig[ijk];
!          if(grad[ijk] > gradmax)  gradmax = grad[ijk];
!        }
!      }
!    }
!  }

!  scaler = 0.01f * sigx2/gradmax;
!  //scaler = 0.01f * sigx2;
!sig2[ijk] = sig[ijk] - scaler * grad[ijk];



  !(+)の係数
!   do i=1,nxpml1
!     sig_x(nx-(i-1))=sig_x(i)
!     kappa_x(nx-(i-1))=kappa_x(i)
!     a_x(nx-(i-1))=a_x(i)
!     kedx(nx-(i-1))=kedx(i)
!     be_x(nx-(i-1))=be_x(i)
!     ce_x(nx-(i-1))=ce_x(i)
!   enddo

!     do j=1,nypml1
!     sig_y(ny-(j-1))=sig_y(j)
!     kappa_y(ny-(j-1))=kappa_y(j)
!     a_y(ny-(j-1))=a_y(j)
!     kedy(ny-(j-1))=kedy(j)
!     be_y(ny-(j-1))=be_y(j)
!     ce_y(ny-(j-1))=ce_y(j)
!   enddo

!       do k=1,nzpml1
!     sig_z(nz-(k-1))=sig_z(k)
!     kappa_z(nz-(k-1))=kappa_z(k)
!     a_z(nz-(k-1))=a_z(k)
!     kedz(nz-(k-1))=kedz(k)
!     be_z(nz-(k-1))=be_z(k)
!     ce_z(nz-(k-1))=ce_z(k)
!   enddo






!psi-update
!444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
    !xe-PML4 loop(-)
    do k = 1,nz-1
        do j = 1,ny-1
            do i = 3,nxpml1
                psi_Ezx1(i,j,k) = be_x(i)*psi_Ezx1(i,j,k)&
                                 +ce_x(i)*(c1*Hy(i,j,k)-c1*Hy(i-1,j,k) + c2*Hy(i+1,j,k)-c2*Hy(i-2,j,k)) / dx
                psi_Eyx1(i,j,k) = be_x(i)*psi_Eyx1(i,j,k)&
                                 +ce_x(i)*(c1*Hz(i,j,k)-c1*Hz(i-1,j,k) + c2*Hz(i+1,j,k)-c2*Hz(i-2,j,k)) / dx
                Ez(i,j,k) = Ez(i,j,k) + cb_z(i,j,k)*psi_Ezx1(i,j,k)
                Ey(i,j,k) = Ey(i,j,k) - cb_y(i,j,k)*psi_Eyx1(i,j,k)
            enddo
        enddo
    enddo
     !xe-PML4 loop(+)
    do k = 1,nz-1
        do j = 1,ny-1
            do i = nx-nxpml1+1,nx-1
                psi_Ezx1(i,j,k) = be_x(i)*psi_Ezx1(i,j,k)&
                                 +ce_x(i)*(c1*Hy(i,j,k)-c1*Hy(i-1,j,k) + c2*Hy(i+1,j,k)-c2*Hy(i-2,j,k)) / dx
                psi_Eyx1(i,j,k) = be_x(i)*psi_Eyx1(i,j,k)&
                                 +ce_x(i)*(c1*Hz(i,j,k)-c1*Hz(i-1,j,k) + c2*Hz(i+1,j,k)-c2*Hz(i-2,j,k)) / dx
                Ez(i,j,k) = Ez(i,j,k) + cb_z(i,j,k)*psi_Ezx1(i,j,k)
                Ey(i,j,k) = Ey(i,j,k) - cb_y(i,j,k)*psi_Eyx1(i,j,k)
            enddo
        enddo
    enddo

    !ye-PML4 loop(-)
    do k = 1,nz-1
        do j = 3,nypml1
            do i = 1,nx-1
                psi_Exy1(i,j,k) = be_y(j)*psi_Exy1(i,j,k)&
                                 +ce_y(j)*(c1*Hz(i,j,k)-c1*Hz(i,j-1,k) + c2*Hz(i,j+1,k)-c2*Hz(i,j-2,k)) / dy
                psi_Ezy1(i,j,k) = be_y(j)*psi_Ezy1(i,j,k)&
                                 +ce_y(j)*(c1*Hx(i,j,k)-c1*Hx(i,j-1,k) + c2*Hx(i,j+1,k)-c2*Hx(i,j-2,k)) / dy
                Ex(i,j,k) = Ex(i,j,k) + cb_x(i,j,k)*psi_Exy1(i,j,k)
                Ez(i,j,k) = Ez(i,j,k) - cb_z(i,j,k)*psi_Ezy1(i,j,k)
            enddo
        enddo
    enddo
   !ye-PML4 loop(+)
    do k = 1,nz-1
        do j = ny-nypml1+1,ny-1
            do i = 1,nx-1
                psi_Exy1(i,j,k) = be_y(j)*psi_Exy1(i,j,k)&
                                 +ce_y(j)*(c1*Hz(i,j,k)-c1*Hz(i,j-1,k) + c2*Hz(i,j+1,k)-c2*Hz(i,j-2,k)) / dy
                psi_Ezy1(i,j,k) = be_y(j)*psi_Ezy1(i,j,k)&
                                 +ce_y(j)*(c1*Hx(i,j,k)-c1*Hx(i,j-1,k) + c2*Hx(i,j+1,k)-c2*Hx(i,j-2,k)) / dy
                Ex(i,j,k) = Ex(i,j,k) + cb_x(i,j,k)*psi_Exy1(i,j,k)
                Ez(i,j,k) = Ez(i,j,k) - cb_z(i,j,k)*psi_Ezy1(i,j,k)
            enddo
        enddo
    enddo

    !ze-PML4 loop(-)
    do k = 3,nzpml1
        do j = 1,ny-1
            do i = 1,nx-1
                psi_Eyz1(i,j,k) = be_z(k)*psi_Eyz1(i,j,k)&
                                 +ce_z(k)*(c1*Hx(i,j,k)-c1*Hx(i,j,k-1) + c2*Hx(i,j,k+1)-c2*Hx(i,j,k-2)) / dz
                psi_Exz1(i,j,k) = be_z(k)*psi_Exz1(i,j,k)&
                                 +ce_z(k)*(c1*Hy(i,j,k)-c1*Hy(i,j,k-1) + c2*Hy(i,j,k+1)-c2*Hy(i,j,k-2)) / dz
                Ey(i,j,k) = Ey(i,j,k) + cb_y(i,j,k)*psi_Eyz1(i,j,k)
                Ex(i,j,k) = Ex(i,j,k) - cb_x(i,j,k)*psi_Exz1(i,j,k)
            enddo
        enddo
    enddo
    !ze-PML4 loop(+)
    do k = nz-nzpml1+1,nz-1
        do j = 1,ny-1
            do i = 1,nx-1
                psi_Eyz1(i,j,k) = be_z(k)*psi_Eyz1(i,j,k)&
                                 +ce_z(k)*(c1*Hx(i,j,k)-c1*Hx(i,j,k-1) + c2*Hx(i,j,k+1)-c2*Hx(i,j,k-2)) / dz
                psi_Exz1(i,j,k) = be_z(k)*psi_Exz1(i,j,k)&
                                 +ce_z(k)*(c1*Hy(i,j,k)-c1*Hy(i,j,k-1) + c2*Hy(i,j,k+1)-c2*Hy(i,j,k-2)) / dz
                Ey(i,j,k) = Ey(i,j,k) + cb_y(i,j,k)*psi_Eyz1(i,j,k)
                Ex(i,j,k) = Ex(i,j,k) - cb_x(i,j,k)*psi_Exz1(i,j,k)
            enddo
        enddo
    enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!psi-update!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !xe-PML loop(-)
!   do k = 1,nz-1
!       do j = 1,ny-1
!           do i = 2,nxpml1
!               psi_Ezx1(i,j,k) = be_x(i)*psi_Ezx1(i,j,k)&
!                                +ce_x(i)*(Hy(i,j,k)-Hy(i-1,j,k)) / dx
!               psi_Eyx1(i,j,k) = be_x(i)*psi_Eyx1(i,j,k)&
!                                +ce_x(i)*(Hz(i,j,k)-Hz(i-1,j,k)) / dx
!               Ez(i,j,k) = Ez(i,j,k) + cb_z(i,j,k)*psi_Ezx1(i,j,k)
!               Ey(i,j,k) = Ey(i,j,k) - cb_y(i,j,k)*psi_Eyx1(i,j,k)
!           enddo
!       enddo
!   enddo
!   !xe-PML loop(+)
!   do k = 1,nz-1
!       do j = 1,ny-1
!           do i = nx-nxpml1+1,nx-1
!               psi_Ezx1(i,j,k) = be_x(i)*psi_Ezx1(i,j,k)&
!                                +ce_x(i)*(Hy(i,j,k)-Hy(i-1,j,k)) / dx
!               psi_Eyx1(i,j,k) = be_x(i)*psi_Eyx1(i,j,k)&
!                                +ce_x(i)*(Hz(i,j,k)-Hz(i-1,j,k)) / dx
!               Ez(i,j,k) = Ez(i,j,k) + cb_z(i,j,k)*psi_Ezx1(i,j,k)
!               Ey(i,j,k) = Ey(i,j,k) - cb_y(i,j,k)*psi_Eyx1(i,j,k)
!           enddo
!       enddo
!   enddo

!   !ye-PML loop(-)
!   do k = 1,nz-1
!       do j = 2,nypml1
!           do i = 1,nx-1
!               psi_Exy1(i,j,k) = be_y(j)*psi_Exy1(i,j,k)&
!                                +ce_y(j)*(Hz(i,j,k)-Hz(i,j-1,k)) / dy
!               psi_Ezy1(i,j,k) = be_y(j)*psi_Ezy1(i,j,k)&
!                                +ce_y(j)*(Hx(i,j,k)-Hx(i,j-1,k)) / dy
!               Ex(i,j,k) = Ex(i,j,k) + cb_x(i,j,k)*psi_Exy1(i,j,k)
!               Ez(i,j,k) = Ez(i,j,k) - cb_z(i,j,k)*psi_Ezy1(i,j,k)
!           enddo
!       enddo
!   enddo
!  !ye-PML loop(+)
!   do k = 1,nz-1
!       do j = ny-nypml1+1,ny-1
!           do i = 1,nx-1
!               psi_Exy1(i,j,k) = be_y(j)*psi_Exy1(i,j,k)&
!                                +ce_y(j)*(Hz(i,j,k)-Hz(i,j-1,k)) / dy
!               psi_Ezy1(i,j,k) = be_y(j)*psi_Ezy1(i,j,k)&
!                                +ce_y(j)*(Hx(i,j,k)-Hx(i,j-1,k)) / dy
!               Ex(i,j,k) = Ex(i,j,k) + cb_x(i,j,k)*psi_Exy1(i,j,k)
!               Ez(i,j,k) = Ez(i,j,k) - cb_z(i,j,k)*psi_Ezy1(i,j,k)
!           enddo
!       enddo
!   enddo

!   !ze-PML loop(-)
!   do k = 2,nzpml1
!       do j = 1,ny-1
!           do i = 1,nx-1
!               psi_Eyz1(i,j,k) = be_z(k)*psi_Eyz1(i,j,k)&
!                                +ce_z(k)*(Hx(i,j,k)-Hx(i,j,k-1)) / dz
!               psi_Exz1(i,j,k) = be_z(k)*psi_Exz1(i,j,k)&
!                                +ce_z(k)*(Hy(i,j,k)-Hy(i,j,k-1)) / dz
!               Ey(i,j,k) = Ey(i,j,k) + cb_y(i,j,k)*psi_Eyz1(i,j,k)
!               Ex(i,j,k) = Ex(i,j,k) - cb_x(i,j,k)*psi_Exz1(i,j,k)
!           enddo
!       enddo
!   enddo
!   !ze-PML loop(+)
!   do k = nz-nzpml1+1,nz-1
!       do j = 1,ny-1
!           do i = 1,nx-1
!               psi_Eyz1(i,j,k) = be_z(k)*psi_Eyz1(i,j,k)&
!                                +ce_z(k)*(Hx(i,j,k)-Hx(i,j,k-1)) / dz
!               psi_Exz1(i,j,k) = be_z(k)*psi_Exz1(i,j,k)&
!                                +ce_z(k)*(Hy(i,j,k)-Hy(i,j,k-1)) / dz
!               Ey(i,j,k) = Ey(i,j,k) + cb_y(i,j,k)*psi_Eyz1(i,j,k)
!               Ex(i,j,k) = Ex(i,j,k) - cb_x(i,j,k)*psi_Exz1(i,j,k)
!           enddo
!       enddo
!   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!field-update loop!!!!!!!!!!!!!!!!!!!!!!!!!
!    
     !x-update
!    do  k=2,nz-1
!        do  j=2,ny-1
!           do i=1,nx-1
!               Ex(i,j,k) = ca_x(i,j,k)*Ex(i,j,k)&
!                           +cb_x(i,j,k)*((Hz(i,j,k)-Hz(i,j-1,k))/kedy(j) &
!                                       - (Hy(i,j,k)-Hy(i,j,k-1))/kedz(k))
!            enddo
!        enddo
!    enddo

!     !y-update
!    do  k=2,nz-1
!        do  j=1,ny-1
!            do i=2,nx-1
!                Ey(i,j,k) = ca_y(i,j,k)*Ey(i,j,k)&
!                           +cb_y(i,j,k)*((Hx(i,j,k)-Hx(i,j,k-1))/kedz(k) &  
!                                       - (Hz(i,j,k)-Hz(i-1,j,k))/kedx(i))
!            enddo
!        enddo
!    enddo
   
!     !z-update
!     do k=1,nz-1
!       do  j=2,ny-1
!            do  i=2,nx-1
!                Ez(i,j,k) = ca_z(i,j,k)*Ez(i,j,k)&
!                           +cb_z(i,j,k)*((Hy(i,j,k)-Hy(i-1,j,k))/kedx(i) &  
!                                       - (Hx(i,j,k)-Hx(i,j-1,k))/kedy(j)) 
!            enddo
!        enddo
!     enddo
        end subroutine CPML_E
