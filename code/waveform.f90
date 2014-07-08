!//////////////////////////////////////////////////////////////////////////
! 送信源の設定
!/////////////////////////////////////////////////////////////////////////
subroutine read_source_3d(istep,t,Hz,Je,Jh)
! subroutine read_source_3d(istep,t,sig,myu,Hz,Je,Jh)

    use const_para
    implicit none

    integer, intent(in)  :: istep
    real(8), intent(in)  :: t
    complex(kind(0d0))   :: signal(nstep)!1st derivatice gaussian
    complex(kind(0d0)), intent(out) :: Je(nstep)!電場ソース
    complex(kind(0d0)), intent(out) :: Jh(nstep)!磁場ソース
!     complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz)!電場ソース
    complex(kind(0d0)), intent(inout) :: Hz(nx,ny,nz)!磁場ソース
!     real(8), intent(in)  :: sig(nx,ny,nz)
!     real(8), intent(in)  :: myu(nx,ny,nz)
    real(8)              :: etaxx(x0,y0,z0)
    real(8), parameter   :: t0 = pi/fmax
    real(8), parameter   :: beta = pi*(fmax**2)

    !1st_derivative gaussian
    signal(istep) = -(2.0d0*beta*(istep*dt-t0)*sqrt(beta/pi))*exp(-beta*(istep*dt-t0)**2.0d0)

    !電場ソースの設定
    etaxx(x0,y0,z0) = (2.0d0*omega0) / sig(x0,y0,z0)
    Je(istep) = dt*etaxx(x0,y0,z0)*signal(istep) /dx/dy/dz

!     Ex(x0,y0,z0) = Ex(x0,y0,z0) &
!                     - dt*etaxx(x0,y0,z0)*signal(istep) /dx/dy/dz

    !磁場ソースの設定
    Jh(istep) = signal(istep)*dt / myu(x0,y0,z0) /dx/dy/dz
    Hz(x0,y0,z0) = Hz(x0,y0,z0) &
                    - signal(istep)*dt / myu(x0,y0,z0) /dx/dy/dz
        end subroutine read_source_3d




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!送信源the first derivative of Gaussian!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Jh(istep) = signal(istep)*dt / myu(x0,y0,z0) /dx/dy/dz ←割る必要あるのか？ ある！！
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine firstderiv_gauss(istep,t,Je,Jh,sig,myu)
!     use const_para
!     implicit none

!     integer, intent(in)  :: istep
!     real(8), intent(in)  :: t
!     complex(kind(0d0))              :: signal(nstep)!1st derivatice gaussian
!     complex(kind(0d0)), intent(out) :: Je(nstep)!電場ソース
!     complex(kind(0d0)), intent(out) :: Jh(nstep)!磁場ソース
!     real(8), intent(in)  :: sig(nx,ny,nz)
!     real(8), intent(in)  :: myu(nx,ny,nz)
!     real(8)              :: etaxx(x0,y0,z0)
!     real(8), parameter   :: t0 = pi/fmax
!     real(8), parameter   :: beta = pi*(fmax**2)

!     !1st_derivative gaussian
! !     Jn(istep) = -(2.0d0*beta*(t-t0)*sqrt(beta/pi))*exp(-beta*(t-t0)**2.0d0)
!     signal(istep) = -(2.0d0*beta*(istep*dt-t0)*sqrt(beta/pi))*exp(-beta*(istep*dt-t0)**2.0d0)

!     !電場ソースの設定
!     etaxx(x0,y0,z0) = (2.0d0*omega0) / sig(x0,y0,z0)

!     Je(istep) = dt*etaxx(x0,y0,z0)*signal(istep) /dx/dy/dz

!     !磁場ソースの設定
!     Jh(istep) = signal(istep)*dt / myu(x0,y0,z0) /dx/dy/dz
!         end subroutine firstderiv_gauss








! !1st derivative gauss////////////////////////////
! subroutine firstderiv_gauss(istep,t,signal)
!     use const_para
!     implicit none
!     integer, intent(in) :: istep
!     real(8), intent(in) :: t
!     real(8), parameter   :: t0 = pi/fmax
!     real(8), parameter   :: beta = pi*(fmax**2)
!     complex(kind(0d0)),intent(out) :: signal(nstep)

!     signal(istep) = -(2.0d0*beta*(istep*dt-t0)*sqrt(beta/pi))*exp(-beta*(istep*dt-t0)**2.0d0)

! end subroutine firstderiv_gauss


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine sinw(istep,t)
!     use const_para
!     implicit none
!     real(8) :: signal(nstep)
!     real(8),parameter ::freq=100.0d0
!     signal(istep) = sin(2.0d0*pi*freq*istep*dt)
! end subroutine sinw

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine sinc(istep,t)
!     use const_para
!     implicit none
!     real(8) :: signal(nstep)
!     real(8),parameter :: freq=1000.0d0
!     signal(istep) = sin(2.0d0*pi*freq*(istep-nstep/2.0d0)*dt)/(2.0d0*pi*freq*istep*dt)
! end subroutine sinc

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine kukei(istep,t)
!     use const_para
!     implicit none
!     real(8) :: signal(nstep)
!     integer, parameter :: fp=1500
!     integer :: fp_r

!     fp_r=mod(istep,fp*4)
!     if(fp_r<fp) then
!         signal(istep)=0.0d0
!         else if(fp*1<=fp_r .and. fp_r<fp*2) then
!             signal(istep) = -1.0d0
!             else if(fp*2<=fp_r .and. fp_r<fp*3) then
!             signal(istep) = 0.0d0
!              else if(fp*3<=fp_r .and. fp_r<fp*4) then
!                 signal(istep) = 1.0d0
!             endif
! end subroutine kukei

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine ricker(istep,t)
!     use const_para
!     implicit none
!         real(8) :: alpha
!         real(8) :: beta
!         real(8) :: signal(nstep)
!         real(8),parameter :: fp =1.0d1
!         real(8),parameter :: Md =1.0d0
!         real(8),parameter :: phase = Md/fp

!         alpha = (pi*fp)**2.0d0*(dble(istep-200)*dt-phase)**2.0d0
!         beta=2.0d0*alpha - 1.0d0
!         signal(istep) = beta *exp(-alpha)
!     end subroutine ricker

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine sin2(istep,t)
!     use const_para
!     implicit none
!         real(8),parameter :: freq = 1.06d0
!         real(8),parameter :: om = 2.0d0*pi*freq
!         real(8),parameter :: Tw = 4.0d0*pi/om
!         real(8) :: signal(nstep)
!         integer           :: iTw
!         iTw = int(Tw/dt)

!         if(istep<iTw) then
!             signal(istep) = 0.5d0*(1.0d0-cos(pi*istep*dt/Tw)*sin(om*istep*dt))
!             elseif(iTw<=istep<nstep) then
!                 signal(istep) = sin(om*istep*dt)
!             endif
!         end subroutine sin2

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine gauss_sin(istep,t)
!     use const_para
!     implicit none
!     real(8),parameter :: om = 10.0d0
!     real(8),parameter :: t0 = pi/om
!     real(8),parameter :: alpha = (2.0d0/t0)**2.0d0
!     real(8),parameter :: wc = 1.0d0*pi/2.0d0/t0
!     real(8) :: signal(nstep)
!     integer,parameter :: it0 = int(2.0d0*t0/dt)

!     if(istep<it0) then
!         signal(istep) = exp(-alpha*(istep*dt-t0)**2) * cos(wc*(istep*dt-t0))
!     elseif(2.0d0*dt<ntstep) then
!         signal(istep) =0.0d0
!     endif
! end subroutine gauss_sin

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine fourier1(istep,t)
!     use const_para
!     implicit none
!     real(8) :: ReF,ImF,om
!     real(8) :: signal(nstep)
!     integer :: n

!     om = 1.0d0/nstep/dt
!     do n=1,nstep
!         do k=1,nstep
!     if(n<nstep) then
!         ReF=0.0d0
!         Imf=0.0d0
!         elseif(k<nstep) then
!             ReF = ReF + signal(k) * cos(2.0d0*pi*k*n/nstep)
!             ImF = ImF + signal(k) * sin(2.0d0*pi*k*n/nstep)
!         endif
!             enddo
!                 enddo
!         end subroutine fourier1

