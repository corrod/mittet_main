!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ファイル出力用サブルーチン
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_EH_J(istep,t,Je,Jh,Ex,Ey,Ez,Hx,Hy,Hz)

    use const_para
    implicit none

    integer :: l
    integer, intent(in) :: istep
    real(8), intent(in) :: t
    complex(kind(0d0)), intent(in) :: Jh(nstep)
    complex(kind(0d0)), intent(in) :: Je(nstep)
    complex(kind(0d0)), intent(in) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)), intent(in) :: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)
    character(5) :: name

    write(13,*) t, real(hz(x0,y0,nz-1)),  aimag(hz(x0,y0,nz-1))    !hz1050.d
    write(14,*) t, real(hz(x0,y0,z0+10)), aimag(hz(x0,y0,z0+10)) !hz1010.d
    write(15,*) t, real(hz(x0,y0,z0+20)), aimag(hz(x0,y0,z0+20))   !hz1030.d
    write(16,*) t, real(hz(x0,y0,z0+30)), aimag(hz(x0,y0,z0+30))
    write(17,*) t, real(Je(istep))      , aimag(Je(istep))       !je_fic.d
    write(18,*) t, real(Jh(istep))      , aimag(Jh(istep))       !jh_fic.d

    write(20,*) t, real(ex(x0,y0,z0)),    aimag(ex(x0,y0,z0))    !ex1000.d
    write(21,*) t, real(ex(x0+10,y0,z0)), aimag(ex(x0+10,y0,z0)) !ex1010.d
    write(22,*) t, real(ex(x0+20,y0,z0)), aimag(ex(x0+20,y0,z0)) !ex1020.d
    write(23,*) t, real(ex(x0+30,y0,z0)), aimag(ex(x0+30,y0,z0)) !ex1030.d



!-----------------シェル用出力-------------------------------
    if (mod(istep,50)==0) then
   l=10000+istep/50
    write(name,"(I5)") l
    open(7,file="hz"//name//".d")
    open(8,file="ex"//name//".d")

        do j=1,ny
            do i=1,nx
                 write(7,*) t,i,j,real(hz(i,j,z0))
                 write(8,*) t,i,j,real(ex(i,j,z0))
            enddo
        enddo
    close(7)
    close(8)
    endif

!!!!境界の監視
    if (mod(istep,50)==0) then
   l=10000+istep/50
    write(name,"(I5)") l
    open(9,file="bd"//name//".d")
            do i=1,nx
                 write(9,*) t,i,real(hz(i,y0,z0))
            enddo
    close(9)
    endif
            end subroutine output_EH_J


 !   open(1,file='Exfield.d', position='append')
 !   open(2,file='Eyfield.d', position='append')
 !   open(3,file='Ezfield.d', position='append')
  !  open(4,file='Hxfield.d', position='append')
  !  open(5,file='Hyfield.d', position='append')
  !  open(6,file='Hzfield.d', position='append')

!    do k=1,nz
!         do j=1,ny
!             do i=1,nx
!               write(1,'(i5,3i4,e12.4)') istep,i,j,k,real(Ex(i,j,k))
!               write(2,'(i5,3i4,e12.4)') istep,i,j,k,real(Ey(i,j,k))
!               write(3,'(i5,3i4,e12.4)') istep,i,j,k,real(Ez(i,j,k))
!               write(4,'(i5,3i4,e12.4)') istep,i,j,k,real(Hx(i,j,k))
!               write(5,'(i5,3i4,e12.4)') istep,i,j,k,real(Hx(i,j,k))
!               write(6,'(i5,3i4,e12.4)') istep,i,j,k,real(Hx(i,j,k))
!            enddo
!          enddo
!    enddo

 !   close(1)
  !  close(2)
   ! close(3)
!    close(4)
!    close(5)
!    close(6)

