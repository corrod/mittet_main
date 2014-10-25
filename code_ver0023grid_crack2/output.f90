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

    !hzレシーバー
    write(31,*) t, real(hz(x1,y0,z_source)),   aimag(hz(x1,y0,z_source))    !hz1000.d
    write(32,*) t, real(hz(x1+5,y0,z_source)), aimag(hz(x1+5,y0,z_source))  !hz1010.d 　
    write(33,*) t, real(hz(x1+10,y0,z_source)), aimag(hz(x1+10,y0,z_source))  !hz1020.d 　
    write(34,*) t, real(hz(x1+15,y0,z_source)), aimag(hz(x1+15,y0,z_source))  !hz1030.d 　
    write(35,*) t, real(hz(x1+20,y0,z_source)), aimag(hz(x1+20,y0,z_source))  !hz1040.d
    write(36,*) t, real(hz(x1+20,y0,z_source)), aimag(hz(x1+20,y0,z_source))  !hz1050.d
    !exレシーバー
    write(20,*) t, real(ex(x1,y0,z_source)),   aimag(ex(x1,y0,z_source))     !ex1000.d
    write(21,*) t, real(ex(x1+5,y0,z_source)), aimag(ex(x1+5,y0,z_source)) !ex1010.d 　
    write(22,*) t, real(ex(x1+10,y0,z_source)), aimag(ex(x1+10,y0,z_source)) !ex1020.d 　
    write(23,*) t, real(ex(x1+15,y0,z_source)), aimag(ex(x1+15,y0,z_source)) !ex1030.d　

    !対称性確認
    write(24,*) t, real(hz(x1-2,y0,z_source)),  aimag(hz(x1+2,y0,z_source))  !hzleft1.d
    write(25,*) t, real(hz(nx-11,y0,z_source)), aimag(hz(nx-11,y0,z_source))  !hzright1.d
    write(26,*) t, real(hz(x1,y0,z_source)), aimag(hz(x1-5,y0,z_source))    !hzleft2.d
    write(27,*) t, real(hz(x1,y0,z_source)), aimag(hz(x1+5,y0,z_source))    !hzright2.d
    write(28,*) t, real(hz(x1,y0,z_source)), aimag(hz(x1-10,y0,z_source))   !hzleft3.d
    write(29,*) t, real(hz(x1,y0,z_source)), aimag(hz(x1+10,y0,z_source))   !hzright3.d


    !自己比較、相互比較用のアウトプット パターン１
    !①
    !② ③
    write(51,*) t, real(hz(x1,y1,z1)), aimag(hz(x1,y1,z1)) !① patern1_1.d
    write(52,*) t, real(hz(x2,y2,z2)), aimag(hz(x2,y2,z2)) !② patern1_2.d
    write(53,*) t, real(hz(x3,y3,z3)), aimag(hz(x3,y3,z3)) !③ patern1_3.d

    !自己比較、相互比較用のアウトプット パターン２
    !①
    !② ③
    write(55,*) t, real(hz(xx1,yy1,zz1)), aimag(hz(xx1,yy1,zz1)) !① patern2_1.d
    write(56,*) t, real(hz(xx2,yy2,zz2)), aimag(hz(xx2,yy2,zz2)) !② patern2_2.d
    write(57,*) t, real(hz(xx3,yy3,zz3)), aimag(hz(xx3,yy3,zz3)) !③ patern2_3.d



!-----------------シェル用出力 MOVIE-------------------------------
    if (mod(istep,50)==0) then
   l=10000+istep/50
    write(name,"(I5)") l
    open(7,file="./out1/hz"//name//".d")
    open(8,file="./out2/ex"//name//".d")
        do k =1,nz
            do i =1,nx
                 write(7,*) t,i,k,real(hz(i,y0,k)) !x-z
            enddo
        enddo

        do j=1,ny
            do i=1,nx
                 write(8,*) t,i,j,real(ex(i,j,plate-1)) !x-y
            enddo
        enddo
    close(7)
    close(8)
    endif

!!!!境界の監視
!     if (mod(istep,50)==0) then
!    l=10000+istep/50
!     write(name,"(I5)") l
!     open(9,file="bd"//name//".d")
!             do i=1,nx
!                  write(9,*) t,i,real(hz(i,y0,z0))
!             enddo
!     close(9)
!     endif
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

