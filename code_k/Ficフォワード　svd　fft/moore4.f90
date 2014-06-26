
!**********************************************************************
! 擬似逆行列の概念を使って行列方程式 A x=b を解く
!
! 言語: FORTRAN90
! コンパイル環境: Compaq Visual Fortran Version 6.5
!                                                   By Takuichi Hirano
!**********************************************************************

module mod_matrix
   integer,parameter :: mm=4096
   integer,parameter :: nn=4096
   real(8),parameter :: t=1.14315354400754132d-3
   integer :: m,n
   complex(kind(0d0)) :: a(mm,nn),b(mm)
end module

!**********************************************************************
! メインプログラム
!**********************************************************************
program main
   use mod_matrix
   implicit none
  integer :: ts,te,t_rate,t_max,diff
  call system_clock(ts)



   call input
   call matrix_pseudo_inversion

   call system_clock(te,t_rate,t_max)
   if (te<ts) then
      diff=t_max-ts+te
  else
     diff=te-ts
  end if
  print "(A,F10.3)","time it took was:", diff/dble(t_rate)

   stop
end program

!----------------------------------------------------------------------
! 行列をファイルから読み込む
!
! [入力ファイル例]
! (matrixa.dat)
! 4 4
! 1 0 0 -1
! 0 1 -1 0
! 0 0 1 1
! 1 1 0 0
!
! (matrixb.dat)
! 4
! 1
! 0
! 0
! 0
!
!----------------------------------------------------------------------
subroutine input
   use mod_matrix
   implicit none

   integer :: i,j
   m=4096
   n=4096

   ! a を読む
   open(10,file='fic_difmat2.dat')
!   read(10,*) m,n
   do i=1,m
      do j=1,m
         read(10,*) a(i,j)
      end do
   end do
   close(10)

   ! b を読む
   open(10,file='dif4_m7_11a')
!   read(10,*) m
   do i=1,m
      read(10,*) b(i)
   end do
   close(10)

   do i=1,m
      if (i.le.85) then
      b(i)=b(i)
      else
      b(i)=(0.0d0,0.0d0)
      end if
      print *, b(i)
   end do

   ! TEST
   open(10,file='test.txt')
   do i=1,m
      do j=1,m
         write(10,*) a(i,j)
      end do
   end do
   write(10,"(/)")
   do i=1,m
      write(10,*) b(i)
   end do
   close(10)
   print *, 'Input succeed.'
   return
end subroutine

!----------------------------------------------------------------------
! 擬似逆行列（一般逆行列, ムーア・ペンローズ逆行列）
! Pseudo Inverse (Generalized Inverse, Moore-Penrose Inverse)
! を使って一般の行列方程式
! A x=b
! の解を計算する。
! det(A)=0
! でも最小２乗誤差 min ||A x'- b|| を満たす解 x' が求まる。
!
! 一般逆行列を計算するルーチンはLapackで見つからないから
! 特異値分解 (Singular Value Decomposition, SVD)
! ルーチンを利用して擬似逆行列による逆演算を行う。
!
! [参考文献]
! G. ストラング著, 線形代数とその応用, 産業図書, 平成4年
!----------------------------------------------------------------------
subroutine   matrix_pseudo_inversion
   use mod_matrix
   implicit none

   integer :: i,j
   integer :: info
   real(8) :: dsum,cond,sing_val(mm),dwork(5*mm)
   complex(kind(0d0)) :: u(mm,nn),vt(mm,nn),zwork(3*mm)
   complex(kind(0d0)) :: a_wk(mm,nn),b_wk1(mm),b_wk2(mm),zsum

   ! a は zgesvd で破壊されてしまうから、コピーして使う
   do i=1,m
      do j=1,n
         a_wk(i,j)=a(i,j)
      end do
   end do

   ! -------- 特異値分解 --------
   ! Singular Value Decomposition (SVD)
   !
   ! Netlib (http://www.netlib.org/)
   ! LAPACK (Linear Algebra PACKage) (http://www.netlib.org/lapack/)
   ! lapack/complex16 (http://www.netlib.org/lapack/complex16/)
   ! "zgesvd.f plus dependencies"
   !
   call zgesvd('A', 'A', m, n, a_wk, mm,&
               sing_val, u, mm, vt, mm, zwork, 3*mm, dwork, info)

   if(info.ne.0) then
      write(*,*) "Sub(pseudo_inverse): argument is illegal for ''zgesvd''."
      write(*,*) "error code",info
      stop
   end if

   ! -------- 特異値 (singular value) 表示 --------
!   open(12,file='sing_val.dat')
!   do i=1,min(m,n)
!   open(12,file='sing_val_test.dat')
!   do i=1,m
!      write(12,*) i,sing_val(i)
!   end do
!   close(12)
!   stop

   ! -------- 条件数 (condition number) 表示 --------
   ! 条件数=(最大特異値)/(最小特異値)

   cond=sing_val(1)
   do i=1,n
      if(sing_val(i) > 1.0d-13) then
         cond=sing_val(i)
      end if
      print *,i
   end do
   cond=sing_val(1)/cond

 !  open(11,file="cond_t7_11t13.dat")
     open(11,file="cond_test.dat")
   write(11,*) "CONDITION NUMBER=",cond
   close(11)

   ! SVD: A    =U    *S    *VT
   !      [m,n] [m,m] [m,n] [n,n]
   ! pseudo_inverse(A)=VT^H*S^(-1)*U^H
   ! M^H means Hermitian conjugate of matrix M

   ! A x = b
   ! x = pseudo_inverse(A) b
   !   = VT^H*S^(-1)*U^H b

   ! b_wk1:= b
   ! [m,1]   [m,1]
   do i=1,m
      b_wk1(i)=b(i)
   end do

   ! b_wk := U^H   b
   ! [m,1]   [m,m] [m,1]
   do i=1,m
      zsum=dcmplx(0.0d0,0.0d0)
      do j=1,m
         zsum=zsum+dconjg(u(j,i))*b_wk1(j)
      end do
      b_wk2(i)=zsum
   end do

   ! b_wk1:= b_wk2
   ! [m,1]   [m,1]
   do i=1,m
      b_wk1(i)=b_wk2(i)
   end do

   ! b_wk := S^(-1)*U^H b
   ! [n,1]   [n,m]  [m,1]
   do i=1,n
      if(sing_val(i) > 1.0d-13) then
         b_wk2(i)=b_wk1(i)/sing_val(i)
      end if
   end do

   ! b_wk1:= b_wk2
   ! [n,1]   [n,1]
   do i=1,n
      b_wk1(i)=b_wk2(i)
   end do

   ! b_wk := VT^H *S^(-1)*U^H b
   ! [n,1]   [n,n] [n,1]
   do i=1,m
      zsum=dcmplx(0.0d0,0.0d0)
      do j=1,m
         zsum=zsum+dconjg(vt(j,i))*b_wk1(j)
      end do
      b_wk2(i)=zsum
   end do

   ! b_wk1:= b_wk2
   ! [n,1]   [n,1]
   do i=1,n
      b_wk1(i)=b_wk2(i)
      b(i)=b_wk2(i)        ! 解を右辺ベクトルに上書き
   end do

   ! -------- ２乗誤差の計算 --------
   ! A x' の計算
   !do i=1,m
   !   zsum=dcmplx(0.0d0)
   !   do j=1,n
   !      zsum=zsum+a(i,j)*b_wk1(j)
   !   end do
   !   b_wk2(i)=zsum
   !end do
   ! || A x'-b || の計算
   !dsum=0.0d0
   !do i=1,n
   !   dsum=dsum+cdabs(b_wk2(i)-b(i))**2
   !end do
   !write(*,*) dsum
   !stop

   ! 結果の表示

   open(15,file='fic4_t7_11t13_85.dat')
   do i=1,n
      write(15,*) i,(i-1)*t,real(b_wk1(i)),aimag(b_wk1(i))
   end do
   close(15)

   return
end subroutine

!
! End of File
!

