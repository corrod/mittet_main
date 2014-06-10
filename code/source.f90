!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!送信源the first derivative of Gaussian!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gaussian(istep,t,Je,Jh,sigma,myu)
    use const_para
    implicit none

    integer, intent(in)  :: istep
    real(8), intent(in)  :: t
    real(8)              :: Jn(nstep)!1st derivatice gaussian
    real(8), intent(out) :: Je(nstep)
    real(8), intent(out) :: Jh(nstep)
    real(8), intent(in)  :: sigma(nx,ny,nz)
    real(8), intent(in)  :: myu(nx,ny,nz)
    real(8)              :: etaxx(x0,y0,z0)
    real(8), parameter   :: t0 = pai/fmax
    real(8), parameter   :: beta = pai*(fmax**2)

    Jn(istep) = -(2.0d0*beta*(t-t0)*sqrt(beta/pai))*exp(-beta*(t-t0)**2)

    !電場ソースの設定
    etaxx(x0,y0,z0) = (2.0d0*omega0)/sigma(x0,y0,z0)

    Je(istep) = dt*etaxx(x0,y0,z0)*Jn(istep)/dx/dy/dz

    !磁場ソースの設定
    Jh(istep) = Jn(istep)*dt/myu(x0,y0,z0)/dx/dy/dz

        end subroutine gaussian
