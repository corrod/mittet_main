!//////////////////////////////////////////////////////////////////////////
! 送信源の設定
!/////////////////////////////////////////////////////////////////////////
program read_source_3d

!     use const_para
    implicit none

    integer  :: istep
    real(8)  :: t
    integer,parameter :: nstep = 500
    real(8),parameter :: fmax = 1.0d3
    real(8),parameter :: pi =3.14159265358979d0
    real(8),parameter :: dt=3.05d-5
    complex(kind(0d0))   :: signal(nstep)!1st derivatice gaussian
    complex(kind(0d0)) :: Je(nstep)!電場ソース
    complex(kind(0d0)) :: Jh(nstep)!磁場ソース
!     complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz)!電場ソース
!     complex(kind(0d0)) :: Hz(nx,ny,nz)!磁場ソース
!     real(8)              :: etaxx(x0,y0,z0)
    real(8), parameter   :: t0 = pi/fmax
    real(8), parameter   :: beta = pi*(fmax**2)

    !1st_derivative gaussian
    do istep =1,nstep
    signal(istep) = -(2.0d0*beta*(istep*dt-t0)*sqrt(beta/pi))*exp(-beta*(istep*dt-t0)**2.0d0)
    write(*,*) signal(istep)
    enddo

    !電場ソースの設定
!     etaxx(x0,y0,z0) = (2.0d0*omega0) / sig(x0,y0,z0)

!     Je(istep) = dt*etaxx(x0,y0,z0)*signal(istep) /dx/dy/dz

!     Ex(x0,y0,z0) = Ex(x0,y0,z0) &
!                     - dt*etaxx(x0,y0,z0)*signal(istep) /dx/dy/dz

    !磁場ソースの設定
!     Jh(istep) = signal(istep)*dt / myu(x0,y0,z0) !/dx/dy/dz　　　
!     Hz(x0,y0,z0) = Hz(x0,y0,z0) &
!                     - signal(istep)*dt / myu(x0,y0,z0) !/dx/dy/dz　　　
        end program read_source_3d
