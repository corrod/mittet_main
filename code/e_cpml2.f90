!!!Convolutional PML_E !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sigmamax amax kappamax の求め方
! メイン部分の計算は通常の電磁波伝播と同じでプサイのぶぶんだけCPML？
! subrouitne 分ける必要ないのかも
!psi部分だけPMLバージョン。
!kappa >=1 real
!sigmai>0 real
!ai=alphai>0 real
!psi loop あと３枚必要？
!演算中の倍精度の表現に注意.d0,
!be_x(5)が0になってしまう問題.１になるはず,epsi0=8.854d-12
!epsi0が小さすぎる!!:q
!cb_x,cb_y,cb_zを仮想領域にあわせる
!係数が1-nxpml1, nx-(nx-nxpml)でことなるので調整
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CPML_E(ex,ey,ez,hx,hy,hz,sigma,cmax)
    use const_para
    implicit none

    integer, parameter :: m = 4, ma = 1.0d0
    integer, parameter :: nxpml1 = 10, nypml1 = 10, nzpml1 = 10 !PML層数
    real(8) :: delta=nxpml1*dx
    real(8), parameter :: kappa_max = 1.0d0 !!!kappa should be [0,10]
    real(8), parameter :: a_max = 0.2d0     !!!if a_max=0, CPML changes to UPML
    real(8), parameter :: nn = 3.0d0 !nn should be [2,6]
    real(8), parameter :: order = 0.0d0 !order should be (0,3]
    real(8), parameter :: optToMax = 10.0d0
    real(8), parameter :: Rcoef=0.01d0 !R should be [10^-2, 10^-12]
    real(8), intent(in) :: cmax
    real(8), intent(in) :: sigma(nx,ny,nz)
    real(8) :: sigma_opt !!!
    real(8) :: sigma_max !!! 導出法確認
    real(8) :: sigma_x(nx),sigma_y(ny),sigma_z(nz)
    real(8) :: a_x(nx),a_y(ny),a_z(nz)
    real(8) :: kappa_x(nx),kappa_y(ny),kappa_z(nz)
    real(8) :: kedx(nx),kedy(ny),kedz(nz) 
    real(8) :: epsi(nx,ny,nz)!1.0d0
    real(8) :: epsir = 1.0d0
    real(8) :: ca_x(nx,ny,nz),ca_y(nx,ny,nz),ca_z(nx,ny,nz)
    real(8) :: cb_x(nx,ny,nz),cb_y(nx,ny,nz),cb_z(nx,ny,nz)
    real(8) :: be_x(nx),be_y(ny),be_z(nz)
    real(8) :: ce_x(nx),ce_y(ny),ce_z(nz)
    complex(kind(0d0)) :: psi_Ezx1(nx,ny,nz),psi_Eyx1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Exy1(nx,ny,nz),psi_Ezy1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Eyz1(nx,ny,nz),psi_Exz1(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: ex(nx,ny,nz),ey(nx,ny,nz),ez(nx,ny,nz)
    complex(kind(0d0)), intent(in)    :: hx(nx,ny,nz),hy(nx,ny,nz),hz(nx,ny,nz)

    epsi(1:nx,1:ny,1:nz)=sigma(1:nx,1:ny,1:nz)/(2.0d0*omega0)


!!!    sigma_max = -(m+1)*lnR0 / (2.0d0*(sqrt(myu/epsi))*nxpml1*dx)  !ln(R(0));反射係数!!!
    sigma_max = (nn+order+1.0d0)*cmax+log(1.0d0/Rcoef) / (2.0d0*delta) * optToMax  !!x方向だけ？

!    sigma_opt = (dble(m)+1.0d0) / (150.0d0*pai*sqrt(epsir)*dx)
!    sigma_max = 0.7d0*sigma_opt


!係数の設定xe
 do i = 1,nx
    if(i<=nxpml1) then
      sigma_x(i) = sigma_max * ((dble(nxpml1)-dble(i))/(dble(nxpml1)-1.0d0))**dble(nn+order)
      kappa_x(i) = 1.0d0 + (kappa_max-1.0d0)*((dble(nxpml1)-dble(i))/(dble(nxpml1)-1.0d0))**dble(nn)
      a_x(i)     = a_max * ((dble(i)-1.0d0)/(dble(nxpml1)-1.0d0))**dble(ma)

      be_x(i)    = exp(-(sigma_x(i)/kappa_x(i)+a_x(i))*dt) !/epsi0)
      ce_x(i)    = sigma_x(i)*(be_x(i)-1.0d0) / (sigma_x(i)+kappa_x(i)*a_x(i)) / kappa_x(i)
      kedx(i)    = kappa_x(i)*dx

    else if(i>=nx-nxpml1+1) then
      sigma_x(i) = sigma_max * ((dble(i)-dble(nx)+1.0d0+dble(nxpml1))/(dble(nxpml1)-1.0d0))**dble(nn+order)
      kappa_x(i) = 1.0d0 + (kappa_max-1.0d0)*((dble(i)-dble(nx)+1.0d0+dble(nxpml1))/(dble(nxpml1)-1.0d0))**dble(nn)
      a_x(i)     = a_max * ((dble(-i)+dble(nx)-1.0d0)/(dble(nxpml1)-1.0d0))**dble(ma)

      be_x(i)    = exp(-(sigma_x(i)/kappa_x(i)+a_x(i))*dt) !/epsi0)
      ce_x(i)    = sigma_x(i)*(be_x(i)-1.0d0) / (sigma_x(i)+kappa_x(i)*a_x(i)) / kappa_x(i)
      kedx(i)    = kappa_x(i)*dx

    else
      sigma_x(i) = 0.0d0
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
    sigma_y(j) = sigma_max * ((dble(nypml1)-dble(j))/(dble(nypml1)-1.0d0))**dble(nn+order)
    kappa_y(j) = 1.0d0 + (kappa_max-1.0d0)*((dble(nypml1)-dble(j))/(dble(nypml1)-1.0d0))**dble(nn)
    a_y(j)     = a_max * ((dble(j)-1.0d0)/(dble(nypml1)-1.0d0))**dble(ma)

    be_y(j)    = exp(-(sigma_y(j)/kappa_y(j)+a_y(j))*dt) !/epsi0)
    ce_y(j)    = sigma_y(j)*(be_y(j)-1.0d0) / (sigma_y(j)+kappa_y(j)*a_y(j)) / kappa_y(j)
    kedy(j)    = kappa_y(j)*dy  !!!

  else if(j>=ny-nypml1+1) then
    sigma_y(j) = sigma_max * ((dble(j)-dble(ny)+1.0d0+dble(nypml1))/(dble(nypml1)-1.0d0))**dble(nn+order)
    kappa_y(j) = 1.0d0 + (kappa_max-1.0d0)*((dble(j)-dble(ny)+1.0d0+dble(nypml1))/(dble(nypml1)-1.0d0))**dble(nn)
    a_y(j)     = a_max * ((dble(-j)+dble(ny)-1.0d0)/(dble(nypml1)-1.0d0))**dble(ma)

    be_y(j)    = exp(-(sigma_y(j)/kappa_y(j)+a_y(j))*dt) !/epsi0)
    ce_y(j)    = sigma_y(j)*(be_y(j)-1.0d0) / (sigma_y(j)+kappa_y(j)*a_y(j)) / kappa_y(j)
    kedy(j)    = kappa_y(j)*dy  !!!

  else
    sigma_y(j) = 0.0d0
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
    sigma_z(k) = sigma_max * ((dble(nzpml1)-dble(k))/(dble(nzpml1)-1.0d0))**dble(nn+order)
    kappa_z(k) = 1.0d0 + (kappa_max-1.0d0)*((dble(nzpml1)-dble(k))/(dble(nzpml1)-1.0d0))**dble(nn)
    a_z(k)     = a_max*((dble(k)-1.0d0)/(dble(nzpml1)-1.0d0))**dble(ma)

    be_z(k)    = exp(-(sigma_z(k)/kappa_z(k)+a_z(k))*dt) !/epsi0)
    ce_z(k)    = sigma_z(k)*(be_z(k)-1.0d0) / (sigma_z(k)+kappa_z(k)*a_z(k)) / kappa_z(k)
    kedz(k)    = kappa_z(k)*dz  !!!

  else if(k>=nz-nzpml1+1) then
    sigma_z(k) = sigma_max * ((dble(k)-dble(nz)+1.0d0+dble(nzpml1))/(dble(nzpml1)-1.0d0))**dble(nn+order)
    kappa_z(k) = 1.0d0 + (kappa_max-1.0d0)*((dble(k)-dble(nz)+1.0d0+dble(nzpml1))/(dble(nzpml1)-1.0d0))**dble(nn)
    a_z(k)     = a_max * ((dble(-k)+dble(nz)-1.0d0)/(dble(nzpml1)-1.0d0))**dble(ma)

    be_z(k)    = exp(-(sigma_z(k)/kappa_z(k)+a_z(k))*dt) !/epsi0)
    ce_z(k)    = sigma_z(k)*(be_z(k)-1.0d0) / (sigma_z(k)+kappa_z(k)*a_z(k)) / kappa_z(k)
    kedz(k)    = kappa_z(k)*dz  !!!

  else
    sigma_z(k) = 0.0d0
    kappa_z(k) = 1.0d0
    a_z(k)     = 0.0d0
    be_z(k)    = 0.0d0
    ce_z(k)    = 0.0d0
    kedz(k)    = kappa_z(k)*dz  
  endif
enddo






    

  !(+)の係数
!   do i=1,nxpml1
!     sigma_x(nx-(i-1))=sigma_x(i)
!     kappa_x(nx-(i-1))=kappa_x(i)
!     a_x(nx-(i-1))=a_x(i)
!     kedx(nx-(i-1))=kedx(i)
!     be_x(nx-(i-1))=be_x(i)
!     ce_x(nx-(i-1))=ce_x(i)
!   enddo

!     do j=1,nypml1
!     sigma_y(ny-(j-1))=sigma_y(j)
!     kappa_y(ny-(j-1))=kappa_y(j)
!     a_y(ny-(j-1))=a_y(j)
!     kedy(ny-(j-1))=kedy(j)
!     be_y(ny-(j-1))=be_y(j)
!     ce_y(ny-(j-1))=ce_y(j)
!   enddo

!       do k=1,nzpml1
!     sigma_z(nz-(k-1))=sigma_z(k)
!     kappa_z(nz-(k-1))=kappa_z(k)
!     a_z(nz-(k-1))=a_z(k)
!     kedz(nz-(k-1))=kedz(k)
!     be_z(nz-(k-1))=be_z(k)
!     ce_z(nz-(k-1))=ce_z(k)
!   enddo





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!psi-update!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !xe-PML loop(-)
    do k = 1,nz-1
        do j = 1,ny-1
            do i = 2,nxpml1
                psi_Ezx1(i,j,k) = be_x(i)*psi_Ezx1(i,j,k)&
                                 +ce_x(i)*(hy(i,j,k)-hy(i-1,j,k)) / dx
                psi_Eyx1(i,j,k) = be_x(i)*psi_Eyx1(i,j,k)&
                                 +ce_x(i)*(hz(i,j,k)-hz(i-1,j,k)) / dx
                ez(i,j,k) = ez(i,j,k) + cb_z(i,j,k)*psi_Ezx1(i,j,k)
                ey(i,j,k) = ey(i,j,k) - cb_y(i,j,k)*psi_Eyx1(i,j,k)
            enddo
        enddo
    enddo

     !xe-PML loop(+)
    do k = 1,nz-1
        do j = 1,ny-1
            do i = nx-nxpml1+1,nx-1
                psi_Ezx1(i,j,k) = be_x(i)*psi_Ezx1(i,j,k)&
                                 +ce_x(i)*(hy(i,j,k)-hy(i-1,j,k)) / dx
                psi_Eyx1(i,j,k) = be_x(i)*psi_Eyx1(i,j,k)&
                                 +ce_x(i)*(hz(i,j,k)-hz(i-1,j,k)) / dx
                ez(i,j,k) = ez(i,j,k) + cb_z(i,j,k)*psi_Ezx1(i,j,k)
                ey(i,j,k) = ey(i,j,k) - cb_y(i,j,k)*psi_Eyx1(i,j,k)
            enddo
        enddo
    enddo




    !ye-PML loop(-)
    do k = 1,nz-1
        do j = 2,nypml1
            do i = 1,nx-1
                psi_Exy1(i,j,k) = be_y(j)*psi_Exy1(i,j,k)&
                                 +ce_y(j)*(hz(i,j,k)-hz(i,j-1,k)) / dy
                psi_Ezy1(i,j,k) = be_y(j)*psi_Ezy1(i,j,k)&
                                 +ce_y(j)*(hx(i,j,k)-hx(i,j-1,k)) / dy
                ex(i,j,k) = ex(i,j,k) + cb_x(i,j,k)*psi_Exy1(i,j,k)
                ez(i,j,k) = ez(i,j,k) - cb_z(i,j,k)*psi_Ezy1(i,j,k)
            enddo
        enddo
    enddo

   !ye-PML loop(+)
    do k = 1,nz-1
        do j = ny-nypml1+1,ny-1
            do i = 1,nx-1
                psi_Exy1(i,j,k) = be_y(j)*psi_Exy1(i,j,k)&
                                 +ce_y(j)*(hz(i,j,k)-hz(i,j-1,k)) / dy
                psi_Ezy1(i,j,k) = be_y(j)*psi_Ezy1(i,j,k)&
                                 +ce_y(j)*(hx(i,j,k)-hx(i,j-1,k)) / dy
                ex(i,j,k) = ex(i,j,k) + cb_x(i,j,k)*psi_Exy1(i,j,k)
                ez(i,j,k) = ez(i,j,k) - cb_z(i,j,k)*psi_Ezy1(i,j,k)
            enddo
        enddo
    enddo




    !ze-PML loop(-)
    do k = 2,nzpml1
        do j = 1,ny-1
            do i = 1,nx-1
                psi_Eyz1(i,j,k) = be_z(k)*psi_Eyz1(i,j,k)&
                                 +ce_z(k)*(hx(i,j,k)-hx(i,j,k-1)) / dz
                psi_Exz1(i,j,k) = be_z(k)*psi_Exz1(i,j,k)&
                                 +ce_z(k)*(hy(i,j,k)-hy(i,j,k-1)) / dz
                ey(i,j,k) = ey(i,j,k) + cb_y(i,j,k)*psi_Eyz1(i,j,k)
                ex(i,j,k) = ex(i,j,k) - cb_x(i,j,k)*psi_Exz1(i,j,k)
            enddo
        enddo
    enddo

    !ze-PML loop(+)
    do k = nz-nzpml1+1,nz-1
        do j = 1,ny-1
            do i = 1,nx-1
                psi_Eyz1(i,j,k) = be_z(k)*psi_Eyz1(i,j,k)&
                                 +ce_z(k)*(hx(i,j,k)-hx(i,j,k-1)) / dz
                psi_Exz1(i,j,k) = be_z(k)*psi_Exz1(i,j,k)&
                                 +ce_z(k)*(hy(i,j,k)-hy(i,j,k-1)) / dz
                ey(i,j,k) = ey(i,j,k) + cb_y(i,j,k)*psi_Eyz1(i,j,k)
                ex(i,j,k) = ex(i,j,k) - cb_x(i,j,k)*psi_Exz1(i,j,k)
            enddo
        enddo
    enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!field-update loop!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k=1,nz
!    do j=1,ny
!    do i=1,nx
!    ca_x(i,j,k) = (1.0d0-sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
!    ca_y(i,j,k) = (1.0d0-sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
!    ca_z(i,j,k) = (1.0d0-sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
!
!    cb_x(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
!    cb_y(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
!    cb_z(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
!    enddo
!    enddo
!    enddo
!    
     !x-update
!    do  k=2,nz-1
!        do  j=2,ny-1
!           do i=1,nx-1
!               ex(i,j,k) = ca_x(i,j,k)*ex(i,j,k)&
!                           +cb_x(i,j,k)*((hz(i,j,k)-hz(i,j-1,k))/kedy(j) &
!                                       - (hy(i,j,k)-hy(i,j,k-1))/kedz(k))
!            enddo
!        enddo
!    enddo

!     !y-update
!    do  k=2,nz-1
!        do  j=1,ny-1
!            do i=2,nx-1
!                ey(i,j,k) = ca_y(i,j,k)*ey(i,j,k)&
!                           +cb_y(i,j,k)*((hx(i,j,k)-hx(i,j,k-1))/kedz(k) &  
!                                       - (hz(i,j,k)-hz(i-1,j,k))/kedx(i))
!            enddo
!        enddo
!    enddo
   
!     !z-update
!     do k=1,nz-1
!       do  j=2,ny-1
!            do  i=2,nx-1
!                ez(i,j,k) = ca_z(i,j,k)*ez(i,j,k)&
!                           +cb_z(i,j,k)*((hy(i,j,k)-hy(i-1,j,k))/kedx(i) &  
!                                       - (hx(i,j,k)-hx(i,j-1,k))/kedy(j)) 
!            enddo
!        enddo
!     enddo
        end subroutine CPML_E
