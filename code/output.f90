!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!ファイル出力用サブルーチン
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_EH(istep,t,Jh,Ex,Ey,Ez,Hx,Hy,Hz)
    use const_para
    implicit none

    integer :: l
    integer, intent(in) :: istep
    real(8), intent(in) :: t
    real(8), intent(in) :: Jh(nstep)
    complex(kind(0d0)), intent(in) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)), intent(in) :: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)
    character(5) :: name


    write(13,*) t, real(hz(x0,y0,z0))
    write(14,*) t, real(hz(x0,y0,z0+10)) 
    write(15,*) t, real(hz(x0,y0,z0+20)) 
    write(16,*) t, real(hz(x0,y0,z0+30)) 
    write(17,*) Jh(istep)



!-----------------シェル用出力-------------------------------
    if (mod(istep,20)==0) then
   l=10000+istep/20
    write(name,"(I5)") l
    open(7,file=name//".d")
        do j=1,ny
            do i=1,nx
                 write(7,*) t,i,j,real(Ex(i,j,z0))
            enddo
        enddo
    close(7)
    endif

    

            end subroutine output_EH

            
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
