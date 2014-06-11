!!!Convolutional PML_E *************************************************************
! sigmamax amax kappamax の求め方
! メイン部分の計算は通常の電磁波伝播と同じでプサイのぶぶんだけCPML？l
! subrouitne 分ける必要ないのかも
!psi部分だけPMLバージョン
!***********************************************************************************
subroutine CPML_E(ex,ey,ez,hx,hy,hz)
    use const_para
    implicit none

    integer,parameter :: m = 4,ma = 4
    integer,parameter :: nxpml1 = 5,nypml1 = 5,nzpml1 = 5 !PML層数
    real(8),parameter :: kappa_max = 1.0d0!!! 整数？
    real(8),parameter :: a_max = 0.1d0     !!!
    real(8) :: sigma_opt !!!
    real(8) :: sigma_max !!! 導出法確認
    real(8) :: sigma(nx,ny,nz)
    real(8) :: sigma_x(nx),sigma_y(ny),sigma_z(nz)
    real(8) :: ca_x(nx,ny,nz),ca_y(nx,ny,nz),ca_z(nx,ny,nz)
    real(8) :: cb_x(nx,ny,nz),cb_y(nx,ny,nz),cb_z(nx,ny,nz)
    real(8) :: kappa_x(nx),kappa_y(ny),kappa_z(nz)
    real(8) :: kedx(nx),kedy(ny),kedz(nz) !kedx(nxpml1),kedy(nypml1),kedz(nzpml1)
    real(8) :: epsi(nx,ny,nz)
    real(8) :: epsir = 1.0d0
    real(8) :: a_x(nx),a_y(ny),a_z(nz)
    real(8) :: be_x(nx),be_y(ny),be_z(nz)
    real(8) :: ce_x(nx),ce_y(ny),ce_z(nz)
    complex(kind(0d0)) :: psi_ezx1(nx,ny,nz),psi_eyx1(nx,ny,nz)
    complex(kind(0d0)) :: psi_exy1(nx,ny,nz),psi_ezy1(nx,ny,nz)
    complex(kind(0d0)) :: psi_eyz1(nx,ny,nz),psi_exz1(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: ex(nx,ny,nz),ey(nx,ny,nz),ez(nx,ny,nz)
    complex(kind(0d0)), intent(in) :: hx(nx,ny,nz),hy(nx,ny,nz),hz(nx,ny,nz)

    sigma_opt = (m+1)/(150d0*pai*sqrt(epsir)*dx)
    sigma_max = 0.7d0*sigma_opt


   do k = 1,nz
   do j = 1,ny
   do i = 1,nx
    sigma_x(i) = sigma_max*((nxpml1-i)/(nxpml1-1))**m
    sigma_y(j) = sigma_max*((nypml1-j)/(nypml1-1))**m
    sigma_z(k) = sigma_max*((nzpml1-k)/(nzpml1-1))**m

    kappa_x(i) = 1.0d0 + (kappa_max-1)*((nxpml1-i)/(nxpml1-1))**m
    kappa_y(j) = 1.0d0 + (kappa_max-1)*((nypml1-j)/(nypml1-1))**m
    kappa_z(k) = 1.0d0 + (kappa_max-1)*((nzpml1-k)/(nzpml1-1))**m

    a_x(i) = a_max*((i-1)/(nxpml1-1))**ma
    a_y(j) = a_max*((j-1)/(nypml1-1))**ma
    a_z(k) = a_max*((k-1)/(nzpml1-1))**ma

    kedx(i) = kappa_x(i)*dx  !kedy=kappa(jdy)dy
    kedy(j) = kappa_y(j)*dy  !!!
    kedz(k) = kappa_z(k)*dz  !!!

!係数の設定
    ca_x(i,j,k) = (1.0d0-sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
    ca_y(i,j,k) = (1.0d0-sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
    ca_z(i,j,k) = (1.0d0-sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))

    cb_x(i,j,k) = (dt/epsi(i,j,k)) / (1+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
    cb_y(i,j,k) = (dt/epsi(i,j,k)) / (1+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
    cb_z(i,j,k) = (dt/epsi(i,j,k)) / (1+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))

    be_x(i) = exp(-(sigma_x(i)/kappa_x(i)+a_x(i))*dt/epsi0)
    be_y(j) = exp(-(sigma_y(j)/kappa_y(j)+a_y(j))*dt/epsi0)
    be_z(k) = exp(-(sigma_z(k)/kappa_z(k)+a_z(k))*dt/epsi0)

    ce_x(i) = sigma_x(i)*(be_x(i)-1) / (sigma_x(i)+kappa_x(i)*a_x(i)) / kappa_x(i)
    ce_y(j) = sigma_y(j)*(be_y(j)-1) / (sigma_y(j)+kappa_y(j)*a_y(j)) / kappa_y(j)
    ce_z(k) = sigma_z(k)*(be_z(k)-1) / (sigma_z(k)+kappa_z(k)*a_z(k)) / kappa_z(k)
   enddo
    enddo
    enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!field-update loop!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !x-update
!    do  k=2,nz-1
 !       do  j=2,ny-1
  !             ex(i,j,k) = ca_x(i,j,k)*ex(i,j,k)&
   !                        +cb_x(i,j,k)*((hz(i,j,k)-hz(i,j-1,k))/kedy(j) &
    !                                   - (hy(i,j,k)-hy(i,j,k-1))/kedz(k))
     !       enddo
      !  enddo
 !   enddo
    !y-update
 !   do  k=2,nz-1
 !       do  j=1,ny-1
 !           do i=2,nx-1
 !               ey(i,j,k) = ca_y(i,j,k)*ey(i,j,k)&
  !                         +cb_y(i,j,k)*((hx(i,j,k)-hx(i,j,k-1))/kedz(k) &  
  !                                     - (hz(i,j,k)-hz(i-1,j,k))/kedx(i))
  !          enddo
   !     enddo
   ! enddo
    !z-update
  !     do  j=2,ny-1
   !         do  i=2,nx-1
    !            ez(i,j,k) = ca_z(i,j,k)*ez(i,j,k)&
     !                      +cb_z(i,j,k)*((hy(i,j,k)-hy(i-1,j,k))/kedx(i) &  
      !                                 - (hx(i,j,k)-hx(i,j-1,k))/kedy(j)) 
 !           enddo
!        enddo
  !   enddo



  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!psi-update!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !x-PML loop
    do k = 1,nz-1
        do j = 1,ny-1
            do i = 2,nxpml1
                psi_ezx1(i,j,k) = be_x(i)*psi_ezx1(i,j,k)&
                                 +ce_x(i)*(hy(i,j,k)-hy(i-1,j,k)) / dx
                psi_eyx1(i,j,k) = be_x(i)*psi_eyx1(i,j,k)&
                                 +ce_x(i)*(hz(i,j,k)-hz(i-1,j,k)) / dx
                ez(i,j,k) = ez(i,j,k)+cb_z(i,j,k)*psi_ezx1(i,j,k)
                ey(i,j,k) = ey(i,j,k)-cb_y(i,j,k)*psi_eyx1(i,j,k)
            enddo
        enddo
    enddo

    !y-PML loop
    do k = 1,nz-1
        do j = 2,nypml1
            do i = 1,nx-1
                psi_exy1(i,j,k) = be_y(j)*psi_exy1(i,j,k)&
                                 +ce_y(j)*(hz(i,j,k)-hz(i,j-1,k)) / dy
                psi_ezy1(i,j,k) = be_y(j)*psi_ezy1(i,j,k)&
                                 +ce_y(j)*(hx(i,j,k)-hx(i,j-1,k)) / dy
                ex(i,j,k) = ex(i,j,k)+cb_x(i,j,k)*psi_exy1(i,j,k)
                ez(i,j,k) = ez(i,j,k)-cb_z(i,j,k)*psi_ezy1(i,j,k)
            enddo
        enddo
    enddo

    !z-PML loop
    do k = 2,nzpml1
        do j = 1,ny-1
            do i = 1,nx-1
                psi_eyz1(i,j,k) = be_z(k)*psi_eyz1(i,j,k)&
                                 +ce_z(k)*(hx(i,j,k)-hx(i,j,k-1)) / dz
                psi_exz1(i,j,k) = be_z(k)*psi_exz1(i,j,k)&
                                 +ce_z(k)*(hy(i,j,k)-hy(i,j,k-1)) / dz
                ey(i,j,k) = ey(i,j,k)+cb_y(i,j,k)*psi_eyz1(i,j,k)
                ex(i,j,k) = ex(i,j,k)-cb_x(i,j,k)*psi_exz1(i,j,k)
            enddo
        enddo
    enddo
        end subroutine CPML_E
