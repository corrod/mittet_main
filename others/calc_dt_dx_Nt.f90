program main
    implicit none
    integer :: i,j,k
    integer :: Nx
    real(8) :: f0,omega0
    real(8) :: dt,dx,fmax,Nt
    real(8) :: sigwa,sigfe,cmin,cmax
    real(8) :: myuwa_r,myufe_r,myuwa,myufe
!     real(8),parameter :: cmin = 1.156d0
!     real(8),parameter :: cmax = 1772.0d0
    real(8),parameter :: myu0 = 1.2566370614d-6
    real(8),parameter :: Glim = 10.4d0
    real(8),parameter :: pi = 3.141592d0

    write(*,*) 'f0=?'
    read(*,*) f0
    omega0 = 2.0d0*pi*f0

    write(*,*) 'sigwa,sigfe=?'
    read(*,*) sigwa,sigfe
    write(*,*) 'sigwa,sigfe',sigwa,sigfe

    write(*,*) 'myuwa_r,myufe_r=?'
    read(*,*) myuwa_r,myufe_r
    myuwa = myu0 *myuwa_r
    myufe = myu0 *myufe_r

    cmax = sqrt(2.0d0*omega0/myuwa/sigwa)
    cmin = sqrt(2.0d0*omega0/myufe/sigfe)

write(*,*) 'cmax,cmin',cmax,cmin
    write(*,*) 'fmax=?'
    read(*,*)  fmax

    dx = 0.999d0*cmin/fmax/Glim
    write(*,*) 'dx<',dx

    dt = (2.0d0*dx)/((3.0d0**0.5d0)*pi*cmax)
    write(*,*) 'dt<',dt

    write(*,*) 'Nx=?'
    read(*,*) Nx
    Nt = 0.6*(3.0d0**0.5d0)*Nx * cmax/cmin!sqrt(sigfe/sigwa)!cmax/cmin
    write(*,*) 'Nt>',Nt

end program main

