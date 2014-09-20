!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!モデル設定
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine model
    use const_para
    implicit none

! ! model1 海水一様モデル______________________
! ! !海水
! sig(1:nx, 1:ny, 1:nz) = sigwa
! myu(1:nx, 1:ny, 1:nz) = myuwa



! model2 欠陥なしモデル_____________________
! !海水
! sig(1:nx, 1:ny, 1:nz) = sigwa
! myu(1:nx, 1:ny, 1:nz) = myuwa

! !鉄板
! sig(1:nx, 1:ny, z0+10:nz) = sigfe
! myu(1:nx, 1:ny, z0+10:nz) = myufe



! !model3 立方体形状欠陥モデル________________
! 海水
sig(1:nx, 1:ny, 1:nz) = sigwa
myu(1:nx, 1:ny, 1:nz) = myuwa

!鉄板
sig(1:nx, 1:ny, z0+10:nz) = sigfe
myu(1:nx, 1:ny, z0+10:nz) = myufe

!欠陥
sig(x0:x0+3, y0-3:y0+3, z0+10:z0+15) = sigwa
myu(x0:x0+3, y0-3:y0+3, z0+10:z0+15) = myuwa



! ! model4 コーティング層挟む____________________
! !海水
! sig(1:nx, 1:ny, 1:nz) = sigwa
! myu(1:nx, 1:ny, 1:nz) = myuwa

! !コーティング
! sig(1:nx, 1:ny, x0+7:x0+9) = sig!　
! myu(1:nx, 1:ny, x0+7:x0+9) = myu!　

! !鉄板
! sig(1:nx, 1:ny, x0+10:nz) = sigfe
! myu(1:nx, 1:ny, x0+10:nz) = myufe

! !欠陥
! sig(x0-3:x0+3, y0-3:y0+3, x0+10:x0+15) = sigwa
! myu(x0-3:x0+3, y0-3:y0+3, x0+10:x0+15) = myuwa



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!モデル出力
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !x-z 鉛直断面
        open(100,file='model_sig_xz.dat')
            do k = 1, nz
                 j = (ny+1)/2
                    do i = 1, nx
                        write(100,*) i,j,k,sig(i,j,k)
                    enddo
            enddo
        close(100)

        open(101,file='model_myu_xz.dat')
            do k = 1, nz
                 j = (ny+1)/2
                    do i = 1, nx
                        write(101,*) i,j,k,myu(i,j,k)
                    enddo
            enddo
        close(101)

        !x-y 水平断面
        open(102,file='model_sig_xy.dat')
            k = nz-10
            do j = 1, ny
                do i = 1, nx
                    write(102,*) i,j,k,sig(i,j,k)
                enddo
            enddo
        close(102)

        open(103,file='model_myu_xy.dat')
            k = nz-10
            do j = 1, ny
                do i = 1, nx
                    write(103,*) i,j,k,myu(i,j,k)
                enddo
            enddo
        close(103)

end subroutine model
