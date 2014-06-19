!磁場計算********************************************************************
!Hx-field
subroutine  HXFIELD(istep,t,Jh,Hx,Ey,Ez,myu)!myu追加
    use const_para
    implicit none

    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
    real(8), intent(in) :: Jh(nstep)
    complex(kind(0d0)), intent(inout) :: Hx(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ey(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ez(nx,ny,nz)
    real(8), intent(in) :: myu(nx,ny,nz)
    real(8) :: CHXLY(nx,ny,nz)
    real(8) :: CHXLZ(nx,ny,nz)

    !係数の設定
    do k=1,nz
       do j=1,ny
           do i=1,nx
              CHXLY(i,j,k) = - dt / myu(i,j,k) / dy
              CHXLZ(i,j,k) = dt / myu(i,j,k) / dz
          enddo
      enddo
    enddo

    !磁場ソースの設定
   ! Jh(istep) = Jn(istep) * dt /myu(x0,y0,z0) /dx/dy/dz


    !波動伝播計
    do k=1,nz-1
        do j=2,ny-1
            do i=2,nx-1
                  Hx(i,j,k) = Hx(i,j,k)&
                        + CHXLY(i,j,k) * (Ez(i,j+1,k) - Ez(i,j,k))&
                        + CHXLZ(i,j,k) * (Ey(i,j,k+1) - Ey(i,j,k))
            enddo
        enddo
    enddo

    !ソース項
!   Hx(x0,y0,z0) = Hx(x0,y0,z0) - Jh(istep)
            endsubroutine HXFIELD




!Hy-field-------------------------------------------------------------
subroutine HYFIELD(istep,t,Jh,Hy,Ex,Ez,myu) !myu追加
    use const_para
    implicit none

    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
    real(8), intent(in) :: myu(nx,ny,nz)
    real(8), intent(in) :: Jh(nstep)
    complex(kind(0d0)), intent(inout) :: Hy(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ez(nx,ny,nz)
    real(8) :: CHYLZ(nx,ny,nz)
    real(8) :: CHYLX(nx,ny,nz)

    !係数の設定
    do k=1,nz
        do j=1,ny
            do i=1,nx
             CHYLZ(i,j,k) = - dt / myu(i,j,k) / dz
             CHYLX(i,j,k) = dt / myu(i,j,k) / dx
          enddo
      enddo
    enddo

    !磁場ソースの設定
  !  Jh(istep) = Jn(istep) * dt /myu(x0,y0,z0) /dx/dy/dz

    !波動伝播計算
    do k=1,nz-1
        do j=2,ny-1
            do i=2,nx-1
                  Hy(i,j,k) = Hy(i,j,k)&
                        + CHYLZ(i,j,k) * (Ex(i,j,k+1) - Ex(i,j,k))&
                        + CHYLX(i,j,k) * (Ez(i+1,j,k) - Ez(i,j,k))
            enddo
        enddo
    enddo

    !ソース項
!   Hy(x0,y0,z0) = Hy(x0,y0,z0) - Jh(istep)
            endsubroutine HYFIELD




!Hz-field--------------------------------------------------------------------------
subroutine HZFIELD(istep,t,Jh,Hz,Ex,Ey,myu) !myu追加
    use const_para
    implicit none

    integer :: l !スクリプト用
    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
    real(8), intent(in) :: myu(nx,ny,nz)
    real(8), intent(in) :: Jh(nstep)
    complex(kind(0d0)), intent(inout) :: Hz(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ey(nx,ny,nz)
    real(8) :: CHZLX(nx,ny,nz)
    real(8) :: CHZLY(nx,ny,nz)
    character(8) :: name

    !係数の設定
    do k=1,nz
       do j=1,ny
         do i=1,nx
            CHZLX(i,j,k) = - dt / myu(i,j,k) / dz
            CHZLY(i,j,k) = dt / myu(i,j,k) / dx
        enddo
      enddo
    enddo

    !磁場ソースの設定
   ! Jh(istep) = Jn(istep) * dt /myu(x0,y0,z0) /dx/dy/dz

    !波動伝播計算
    do k=1,nz-1
        do j=2,ny-1
            do i=2,nx-1
                 Hz(i,j,k) = Hz(i,j,k)&
                        + CHZLX(i,j,k) * (Ey(i+1,j,k) - Ey(i,j,k))&
                        + CHZLY(i,j,k) * (Ex(i,j+1,k) - Ex(i,j,k))
            enddo
        enddo
    enddo

    !ソース項
    Hz(x0,y0,z0) = Hz(x0,y0,z0) - Jh(istep)


    !ソース位置、ソースから離れた位置でのhzの磁場の時間分布
   write(8,*) 't,jh,hz(x0,y0,z0+10)',t,Jh(istep),real(hz(x0,y0,z0+10))

!--------シェル用出力------
!   if (mod(istep,5)==0) then
!    l=10000+istep/5
!    write(name,"(I5)") l
!    open(9,file=name//".d")
!    do j=1,ny
!        do i=1,nx
!            write(2,*) i,j,real(hz(i,j,z0))
!        enddo
!    enddo
!    close(9)
!   endif
            endsubroutine HZFIELD
